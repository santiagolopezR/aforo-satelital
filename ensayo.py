"""
Ensayo: serie larga de una rotación para ver cómo sube y baja el forraje con los pastoreos.

    python ensayo.py "LA HOLANDA" 150

Junta todas las imágenes posibles de los últimos N días:
  - Sentinel-2 (10 m, cada 5 días)
  - Landsat 8 y 9 (30 m, cada 8 días) para llenar huecos de nubes
  - Sentinel-1 radar (ve a través de las nubes; VH sube con la biomasa)
y detecta los pastoreos (bajones fuertes) y el tiempo entre uno y otro (días del ciclo).

Salida: docs/ensayo/<rotación>.json, que muestra docs/ensayo/index.html. No toca el mapa principal.
"""
import datetime as dt
import json
import re
import statistics
import sys

import numpy as np
import planetary_computer
import pystac_client
import rasterio
from rasterio.features import geometry_mask
from shapely.geometry import box, mapping
from shapely.ops import unary_union

import aforo

CARPETA = aforo.RAIZ / "docs" / "ensayo"
aforo.guardar_pixeles = lambda *a, **k: None  # el ensayo no cambia las imágenes del mapa de detalle
log = aforo.log


def procesar_landsat(item, potreros, config, exclusiones):
    """NDVI por potrero con Landsat Collection 2 nivel 2 (reflectancia = DN · 0,0000275 − 0,2)."""
    epsg = aforo.epsg_de(item, "red")
    geoms = aforo.preparar_geoms(potreros, epsg, config["recorte_borde_m"])
    bounds = unary_union(list(geoms.values())).bounds
    try:
        rojo, transf, caja = aforo.leer_banda(item.assets["red"].href, bounds)
    except aforo.WindowError:
        return []
    nir = aforo.leer_banda(item.assets["nir08"].href, caja)[0]
    qa = aforo.leer_banda(item.assets["qa_pixel"].href, caja)[0]
    if nir.shape != rojo.shape or qa.shape != rojo.shape:
        return []
    r = rojo.astype("float32") * 2.75e-5 - 0.2
    n = nir.astype("float32") * 2.75e-5 - 0.2
    # Bits de calidad: 1 nube dilatada, 3 nube, 4 sombra de nube; 6 = despejado.
    malo = (qa & ((1 << 1) | (1 << 3) | (1 << 4))) != 0
    validos = (rojo > 0) & (nir > 0) & ~malo & ((qa & (1 << 6)) != 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ndvi = (n - r) / (n + r)
    validos &= np.isfinite(ndvi)
    excluido = aforo.mascara_excluida(rojo.shape, transf, epsg, exclusiones, config)
    filas = []
    for p in potreros:
        g = geoms[p["nombre"]]
        if not g.intersects(box(*caja)):
            continue
        dentro = geometry_mask([mapping(g)], out_shape=rojo.shape, transform=transf, invert=True)
        pasto = dentro & ~excluido
        total = int(pasto.sum())
        if total < 3:
            continue
        ok = pasto & validos
        if ok.sum() / total < config["min_pct_valido"]:
            continue
        filas.append({"fecha": item.datetime.date().isoformat(), "potrero": p["nombre"], "especie": p["especie"],
                      "sensor": "L", "escena": item.id, "pct_valido": round(float(ok.sum() / total), 3),
                      "ndvi": round(float(np.nanmean(ndvi[ok])), 4), "hist": aforo.histograma(ndvi[ok])})
    return filas


def armonizar_landsat(mejor):
    """
    Landsat lee distinto que Sentinel-2. Con los pares del mismo potrero y ±1 día (ambos despejados)
    se ajusta NDVI_S2 = a + b · NDVI_Landsat y se aplica a Landsat.
    """
    s2 = {(r["potrero"], r["fecha"]): r for r in mejor.values() if r["sensor"] == "S2" and float(r["pct_valido"]) >= 0.8}
    xs, ys = [], []
    for r in mejor.values():
        if r["sensor"] != "L" or float(r["pct_valido"]) < 0.8:
            continue
        f = dt.date.fromisoformat(r["fecha"])
        for d in (0, -1, 1):
            o = s2.get((r["potrero"], (f + dt.timedelta(days=d)).isoformat()))
            if o:
                xs.append(float(r["ndvi"])); ys.append(float(o["ndvi"]))
                break
    if len(xs) < 10:
        return {"a": 0.0, "b": 1.0, "pares": len(xs)}
    b, a = np.polyfit(xs, ys, 1)
    return {"a": round(float(a), 4), "b": round(float(b), 4), "pares": len(xs)}


def hist_armonizado(hist, a, b):
    """Pasa cada clase del histograma de NDVI por la recta de armonización."""
    nuevo = {}
    for par in hist.split():
        i, c = (int(v) for v in par.split(":"))
        centro = a + b * (aforo.HIST_MIN + (i + 0.5) * aforo.HIST_PASO)
        j = int(min(max((centro - aforo.HIST_MIN) / aforo.HIST_PASO, 0), aforo.HIST_N - 1))
        nuevo[j] = nuevo.get(j, 0) + c
    return " ".join(f"{j}:{c}" for j, c in sorted(nuevo.items()))


def pastoreos_vistos(serie, caida):
    """
    Bajones fuertes entre observaciones despejadas seguidas = pastoreo entre esas fechas.
    Se mide contra la mediana de la rotación ese día ("rel"): si toda la rotación baja junta es la imagen
    (bruma, atmósfera), no el ganado, porque solo unos pocos potreros se pastorean a la vez.
    """
    eventos = []
    for a, b in zip(serie, serie[1:]):
        if b["rel"] - a["rel"] <= -caida:
            fa, fb = dt.date.fromisoformat(a["f"]), dt.date.fromisoformat(b["f"])
            medio = fa + (fb - fa) / 2
            if eventos and (medio - dt.date.fromisoformat(eventos[-1]["medio"])).days < 12:
                continue  # el mismo pastoreo visto en dos bajones seguidos
            eventos.append({"desde": a["f"], "hasta": b["f"], "medio": medio.isoformat(),
                            "antes": a["kg"], "despues": b["kg"]})
    for e1, e2 in zip(eventos, eventos[1:]):
        e2["dias_desde_anterior"] = (dt.date.fromisoformat(e2["medio"]) - dt.date.fromisoformat(e1["medio"])).days
    return eventos


def main(rotacion, dias):
    config = json.loads((aforo.RAIZ / "config.json").read_text(encoding="utf-8"))
    potreros = [p for p in aforo.cargar_potreros(config) if p["rotacion"] == rotacion.upper()]
    if not potreros:
        sys.exit(f"No hay potreros de la rotación {rotacion}")
    log(f"Ensayo {rotacion}: {len(potreros)} potreros, últimos {dias} días")
    exclusiones, huella = aforo.cargar_exclusiones(), aforo.huella_mascara(config)
    hoy = dt.date.today()
    desde = hoy - dt.timedelta(days=dias)
    bbox = list(unary_union([p["geom"] for p in potreros]).bounds)
    cat = pystac_client.Client.open(aforo.STAC_URL, modifier=planetary_computer.sign_inplace)
    filas = []

    def recorrer(items, nombre, fn):
        log(f"{nombre}: {len(items)} escenas")
        for it in sorted(items, key=lambda i: i.datetime):
            try:
                nuevas = fn(it)
            except Exception as e:
                log(f"  ! {it.id}: {e}")
                continue
            if nuevas:
                log(f"  {it.datetime.date()} {nombre}: {len(nuevas)} potreros")
            filas.extend(nuevas)

    recorrer(aforo.buscar(cat, "sentinel-2-l2a", bbox, desde, hoy, {"eo:cloud_cover": {"lt": config["max_nubes_escena"]}}),
             "Sentinel-2", lambda it: aforo.procesar_s2(it, potreros, config, exclusiones, huella))
    recorrer(aforo.buscar(cat, "landsat-c2-l2", bbox, desde, hoy,
                          {"eo:cloud_cover": {"lt": config["max_nubes_escena"]}, "platform": {"in": ["landsat-8", "landsat-9"]}}),
             "Landsat", lambda it: procesar_landsat(it, potreros, config, exclusiones))
    recorrer(aforo.buscar(cat, "sentinel-1-rtc", bbox, desde, hoy),
             "Sentinel-1", lambda it: aforo.procesar_s1(it, potreros, config, exclusiones, huella))

    # Un dato por día, sensor y potrero (el de más píxeles despejados).
    mejor = {}
    for r in filas:
        k = (r["fecha"], r["potrero"], r["sensor"])
        if k not in mejor or float(r["pct_valido"]) > float(mejor[k]["pct_valido"]):
            mejor[k] = r
    aforos = aforo.leer_aforos(potreros)
    arm = armonizar_landsat(mejor)
    log(f"Armonización Landsat -> Sentinel-2: NDVI_S2 = {arm['a']} + {arm['b']} · NDVI_L ({arm['pares']} pares)")
    for r in mejor.values():
        if r["sensor"] == "L" and r.get("hist"):
            r["hist"] = hist_armonizado(r["hist"], arm["a"], arm["b"])
            r["ndvi"] = round(arm["a"] + arm["b"] * float(r["ndvi"]), 4)
        if r["sensor"] in ("S2", "L") and r.get("hist"):
            r["kg"] = aforo.estimar_por_pixel(r["hist"], r["especie"], config)[1]

    # Mediana de la rotación por fecha y sensor (solo fechas con suficientes potreros despejados).
    por_fecha = {}
    for r in mejor.values():
        if "kg" in r and float(r["pct_valido"]) >= 0.8:
            por_fecha.setdefault((r["fecha"], r["sensor"]), []).append(r["kg"])
    minimo = max(5, len(potreros) // 3)
    mediana = {k: statistics.median(v) for k, v in por_fecha.items() if len(v) >= minimo}

    salida = {"rotacion": rotacion.upper(), "desde": desde.isoformat(), "hasta": hoy.isoformat(),
              "generado": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"), "potreros": [],
              "armonizacion": arm,
              "mediana_rotacion": [{"f": f, "s": sen, "kg": round(v, 3)} for (f, sen), v in sorted(mediana.items())]}
    ciclos = []
    for p in potreros:
        mias = sorted((r for r in mejor.values() if r["potrero"] == p["nombre"]), key=lambda r: r["fecha"])
        optico = []
        for r in mias:
            if "kg" in r:
                med = mediana.get((r["fecha"], r["sensor"]))
                optico.append({"f": r["fecha"], "s": r["sensor"], "kg": r["kg"], "ndvi": float(r["ndvi"]),
                               "v": float(r["pct_valido"]),
                               "rel": round(r["kg"] - med, 3) if med is not None else None})
        radar = [{"f": r["fecha"], "vh": float(r["vh_db"]), "vv": float(r["vv_db"])}
                 for r in mias if r["sensor"] == "S1" and r.get("vh_db") not in ("", None)]
        limpias = [o for o in optico if o["v"] >= 0.8 and o["rel"] is not None]
        eventos = pastoreos_vistos(limpias, config["pastoreo"]["caida_minima_kg"])
        ciclos += [e["dias_desde_anterior"] for e in eventos if "dias_desde_anterior" in e]
        salida["potreros"].append({
            "nombre": p["nombre"], "especie": p["especie"], "optico": optico, "radar": radar,
            "aforos": [{"f": a["fecha"], "kg": a["kg_fv_m2"]} for a in aforos if a["potrero"] == p["nombre"]],
            "pastoreos": eventos,
        })
    salida["ciclos"] = ciclos
    salida["ciclo_mediana"] = statistics.median(ciclos) if ciclos else None

    CARPETA.mkdir(parents=True, exist_ok=True)
    nombre = re.sub(r"[^a-z0-9]+", "_", aforo.sin_tildes(rotacion).lower()).strip("_")
    (CARPETA / f"{nombre}.json").write_text(json.dumps(salida, ensure_ascii=False), encoding="utf-8")
    indice = CARPETA / "indice.json"
    lista = json.loads(indice.read_text(encoding="utf-8")) if indice.exists() else {}
    lista[nombre] = {"rotacion": rotacion.upper(), "dias": dias, "generado": salida["generado"]}
    indice.write_text(json.dumps(lista, ensure_ascii=False, indent=1), encoding="utf-8")
    n_opt = sum(len(x["optico"]) for x in salida["potreros"])
    log(f"Listo: {n_opt} datos ópticos, {sum(len(x['radar']) for x in salida['potreros'])} de radar, "
        f"{sum(len(x['pastoreos']) for x in salida['potreros'])} pastoreos vistos, ciclo mediano {salida['ciclo_mediana']} días "
        f"({len(ciclos)} intervalos)")


if __name__ == "__main__":
    with rasterio.Env(**aforo.OPCIONES_GDAL):
        main(sys.argv[1] if len(sys.argv) > 1 else "LA HOLANDA", int(sys.argv[2]) if len(sys.argv) > 2 else 150)
