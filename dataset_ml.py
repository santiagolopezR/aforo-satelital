"""
Arma la tabla para el modelo de machine learning: una fila por aforo de campo con todo lo que se pudo
medir antes de ese aforo.

    python dataset_ml.py altos_de_medina [minutos]

Fuentes (todas gratuitas):
  - Sentinel-2: 10 bandas (azul, verde, rojo, 3 de borde rojo, infrarrojo cercano x2, onda corta x2) e índices.
  - Sentinel-1 radar: VV, VH, VH−VV, RVI (ve a través de las nubes).
  - Landsat 8/9: rojo, infrarrojo, onda corta y temperatura de superficie (térmica).
  - Clima (Open-Meteo, reanálisis ERA5): lluvia, temperatura y evapotranspiración de los días previos.
  - Manejo (del registro de aforos): días de descanso, mes.

Cada fuente se guarda por escena en data/<finca>/ml_<fuente>.csv (incremental, se detiene a los N minutos).
La tabla final queda en data/<finca>/dataset_ml.csv.
"""
import csv
import datetime as dt
import json
import sys
import time

import numpy as np
import planetary_computer
import pystac_client
import rasterio
import requests
from rasterio.features import geometry_mask
from shapely.geometry import box, mapping
from shapely.ops import unary_union

import aforo

log = aforo.log
BANDAS_S2 = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
VENTANA = {"s2": 5, "s1": 12, "l": 8}   # días antes del aforo en que se busca cada sensor


def _indices(b):
    """Índices por píxel a partir de reflectancias (diccionario banda -> arreglo)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        n, r, g, bl = b["B08"], b["B04"], b["B03"], b["B02"]
        return {
            "ndvi": (n - r) / (n + r),
            "evi": 2.5 * (n - r) / (n + 6 * r - 7.5 * bl + 1),
            "savi": 1.5 * (n - r) / (n + r + 0.5),
            "msavi": (2 * n + 1 - np.sqrt(np.maximum((2 * n + 1) ** 2 - 8 * (n - r), 0))) / 2,
            "gndvi": (n - g) / (n + g),
            "ndre1": (b["B8A"] - b["B05"]) / (b["B8A"] + b["B05"]),
            "ndre2": (b["B8A"] - b["B06"]) / (b["B8A"] + b["B06"]),
            "cire": b["B07"] / b["B05"] - 1,
            "ndmi": (b["B8A"] - b["B11"]) / (b["B8A"] + b["B11"]),   # humedad / materia seca
            "nbr": (b["B8A"] - b["B12"]) / (b["B8A"] + b["B12"]),
            "ndwi": (g - n) / (g + n),
            "stvi": b["B11"] / b["B12"],                              # residuo seco vs verde
        }


def features_s2(item, potreros, config, exclusiones):
    epsg = aforo.epsg_de(item, "B04")
    geoms = aforo.preparar_geoms(potreros, epsg, config["recorte_borde_m"])
    bounds = unary_union(list(geoms.values())).bounds
    try:
        rojo, transf, caja = aforo.leer_banda(item.assets["B04"].href, bounds)
    except aforo.WindowError:
        return []
    forma = rojo.shape
    base = float(item.properties.get("s2:processing_baseline", "0") or 0)
    desfase = 1000 if base >= 4.0 else 0
    crudo = {"B04": rojo}
    for bnd in BANDAS_S2:
        if bnd != "B04":
            crudo[bnd] = aforo.leer_banda(item.assets[bnd].href, caja, forma)[0]
    scl = aforo.leer_banda(item.assets["SCL"].href, caja, forma)[0]
    refl = {k: (v.astype("float32") - desfase) / 10000 for k, v in crudo.items()}
    validos = np.isin(scl, aforo.SCL_VALIDAS) & (rojo > 0) & (crudo["B08"] > 0)
    idx = _indices(refl)
    validos &= np.isfinite(idx["ndvi"])
    excluido = aforo.mascara_excluida(forma, transf, epsg, exclusiones, config)
    filas = []
    for p in potreros:
        g = geoms[p["nombre"]]
        if not g.intersects(box(*caja)):
            continue
        pasto = geometry_mask([mapping(g)], out_shape=forma, transform=transf, invert=True) & ~excluido
        total = int(pasto.sum())
        if total < config["min_pixeles"]:
            continue
        ok = pasto & validos
        if ok.sum() / total < config.get("min_despejado_calibrar", 0.8):
            continue
        fila = {"fecha": item.datetime.date().isoformat(), "potrero": p["nombre"], "escena": item.id,
                "s2_pixeles": int(ok.sum())}
        for k, v in refl.items():
            fila[f"s2_{k}"] = round(float(np.nanmean(v[ok])), 5)
        for k, v in idx.items():
            vals = v[ok]
            vals = vals[np.isfinite(vals)]
            if len(vals):
                fila[f"s2_{k}"] = round(float(np.mean(vals)), 5)
                if k in ("ndvi", "ndmi", "ndre1"):
                    for q in (10, 50, 90):
                        fila[f"s2_{k}_p{q}"] = round(float(np.percentile(vals, q)), 5)
                    fila[f"s2_{k}_std"] = round(float(np.std(vals)), 5)
        filas.append(fila)
    return filas


def features_landsat(item, potreros, config, exclusiones):
    epsg = aforo.epsg_de(item, "red")
    geoms = aforo.preparar_geoms(potreros, epsg, config["recorte_borde_m"])
    bounds = unary_union(list(geoms.values())).bounds
    try:
        rojo, transf, caja = aforo.leer_banda(item.assets["red"].href, bounds)
    except aforo.WindowError:
        return []
    forma = rojo.shape
    lee = lambda a: aforo.leer_banda(item.assets[a].href, caja, forma)[0]
    nir, sw1, sw2, qa = lee("nir08"), lee("swir16"), lee("swir22"), lee("qa_pixel")
    termica = lee("lwir11") if "lwir11" in item.assets else None
    sr = lambda x: x.astype("float32") * 2.75e-5 - 0.2
    r, n, s1, s2 = sr(rojo), sr(nir), sr(sw1), sr(sw2)
    malo = (qa & ((1 << 1) | (1 << 3) | (1 << 4))) != 0
    validos = (rojo > 0) & (nir > 0) & ~malo & ((qa & (1 << 6)) != 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ind = {"ndvi": (n - r) / (n + r), "ndmi": (n - s1) / (n + s1), "nbr": (n - s2) / (n + s2)}
    excluido = aforo.mascara_excluida(forma, transf, epsg, exclusiones, config)
    filas = []
    for p in potreros:
        g = geoms[p["nombre"]]
        if not g.intersects(box(*caja)):
            continue
        pasto = geometry_mask([mapping(g)], out_shape=forma, transform=transf, invert=True, all_touched=True) & ~excluido
        total = int(pasto.sum())
        if total < 1:
            continue
        ok = pasto & validos
        if ok.sum() / total < 0.8:
            continue
        fila = {"fecha": item.datetime.date().isoformat(), "potrero": p["nombre"], "escena": item.id}
        for k, v in ind.items():
            fila[f"l_{k}"] = round(float(np.nanmean(v[ok])), 5)
        fila["l_swir1"], fila["l_swir2"] = round(float(np.mean(s1[ok])), 5), round(float(np.mean(s2[ok])), 5)
        if termica is not None:
            # Temperatura de superficie (Collection 2): K = DN · 0,00341802 + 149
            fila["l_temp_c"] = round(float(np.mean(termica[ok].astype("float32") * 0.00341802 + 149.0 - 273.15)), 2)
        filas.append(fila)
    return filas


def clima(potreros, desde, hasta):
    """Lluvia, temperatura y evapotranspiración diarias (Open-Meteo, reanálisis ERA5) en el centro de la finca."""
    c = unary_union([p["geom"] for p in potreros]).centroid
    r = requests.get("https://archive-api.open-meteo.com/v1/archive", timeout=60, params={
        "latitude": round(c.y, 3), "longitude": round(c.x, 3),
        "start_date": desde.isoformat(), "end_date": hasta.isoformat(),
        "daily": "precipitation_sum,temperature_2m_mean,temperature_2m_min,et0_fao_evapotranspiration,shortwave_radiation_sum",
        "timezone": "America/Bogota"})
    r.raise_for_status()
    d = r.json()["daily"]
    return {f: {k: d[k][i] for k in d if k != "time"} for i, f in enumerate(d["time"])}


def _cargar(ruta):
    if not ruta.exists():
        return []
    with ruta.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _guardar(ruta, filas):
    claves = []
    for f in filas:
        claves += [k for k in f if k not in claves]
    ruta.parent.mkdir(parents=True, exist_ok=True)
    with ruta.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=claves)
        w.writeheader()
        w.writerows(filas)


def main(slug, minutos):
    inicio = time.time()
    cb = json.loads((aforo.RAIZ / "config.json").read_text(encoding="utf-8"))
    f = aforo.fincas_de(cb)[slug]
    aforo.usar_finca(slug, f)
    config = aforo.config_de_finca(cb, f)
    potreros = aforo.cargar_potreros(config)
    exclusiones = aforo.cargar_exclusiones()
    carpeta = aforo.HISTORIAL.parent
    aforos = aforo.leer_aforos(potreros)
    fechas = sorted({dt.date.fromisoformat(a["fecha"]) for a in aforos})
    bbox = list(unary_union([p["geom"] for p in potreros]).bounds)
    cat = pystac_client.Client.open(aforo.STAC_URL, modifier=planetary_computer.sign_inplace)

    # --- Sentinel-2: solo escenas que ya se sabe que dejan algún potrero despejado (historial largo + reciente)
    vistos = {}
    for fuente in (aforo.HISTORIAL_LARGO, aforo.HISTORIAL):
        for r in _cargar(fuente):
            if r["sensor"] == "S2" and float(r["pct_valido"] or 0) >= 0.8:
                vistos.setdefault(r["escena"], r["fecha"])
    utiles = {e for e, fe in vistos.items()
              if any(0 <= (fa - dt.date.fromisoformat(fe)).days <= VENTANA["s2"] for fa in fechas)}

    def recorrer(nombre, ruta, ids_items, fn):
        hechas = _cargar(ruta)
        revisadas = {r["escena"] for r in hechas}
        marca = ruta.with_suffix(".revisadas")
        revisadas |= set(marca.read_text(encoding="utf-8").split()) if marca.exists() else set()
        pend = [it for it in ids_items if it.id not in revisadas]
        log(f"{nombre}: {len(ids_items)} escenas útiles, {len(pend)} por procesar")
        for it in pend:
            if time.time() - inicio > minutos * 60:
                log("  tiempo cumplido, sigue en la próxima corrida")
                break
            try:
                hechas += fn(it)
            except Exception as e:
                log(f"  ! {it.id}: {e}")
                continue
            revisadas.add(it.id)
        _guardar(ruta, hechas)
        marca.write_text("\n".join(sorted(revisadas)), encoding="utf-8")
        return hechas

    def items_por_rango(coleccion, ventana, filtro=None):
        out = {}
        for a, b in rangos_de(fechas, ventana):
            for it in aforo.buscar(cat, coleccion, bbox, a, b, filtro):
                out[it.id] = it
        return list(out.values())

    if utiles:
        s2_items = [it for it in items_por_rango("sentinel-2-l2a", VENTANA["s2"]) if it.id in utiles]
        recorrer("Sentinel-2 (10 bandas)", carpeta / "ml_s2.csv", s2_items,
                 lambda it: features_s2(it, potreros, config, exclusiones))
    s1_items = items_por_rango("sentinel-1-rtc", VENTANA["s1"])
    recorrer("Sentinel-1 radar", carpeta / "ml_s1.csv", s1_items,
             lambda it: [{"fecha": r["fecha"], "potrero": r["potrero"], "escena": r["escena"],
                          "s1_vv": r["vv_db"], "s1_vh": r["vh_db"]}
                         for r in aforo.procesar_s1(it, potreros, config, exclusiones, "")])
    l_items = items_por_rango("landsat-c2-l2", VENTANA["l"],
                              {"eo:cloud_cover": {"lt": 80}, "platform": {"in": ["landsat-8", "landsat-9"]}})
    recorrer("Landsat", carpeta / "ml_landsat.csv", l_items, lambda it: features_landsat(it, potreros, config, exclusiones))

    # --- Clima
    met = clima(potreros, fechas[0] - dt.timedelta(days=40), min(fechas[-1], dt.date.today() - dt.timedelta(days=6)))

    # --- Tabla final: por aforo, la observación más cercana ANTES (o el mismo día) de cada fuente
    def indexar(ruta):
        idx = {}
        for r in _cargar(ruta):
            idx.setdefault(r["potrero"], []).append(r)
        return idx
    fuentes = {"s2": indexar(carpeta / "ml_s2.csv"), "s1": indexar(carpeta / "ml_s1.csv"),
               "l": indexar(carpeta / "ml_landsat.csv")}
    tabla = []
    for a in aforos:
        fa = dt.date.fromisoformat(a["fecha"])
        fila = {"fecha": a["fecha"], "potrero": a["potrero"], "finca": a["potrero"].rsplit(" ", 1)[0],
                "aforo_kg_fv_m2": a["kg_fv_m2"], "mes": fa.month, "anio": fa.year}
        obs = {k: v for k, v in (x.split("=") for x in (a.get("observaciones") or "").split("; ") if "=" in x)}
        fila["dias_descanso"] = obs.get("dias_rotacion", "")
        for fuente, idx in fuentes.items():
            cands = [(fa - dt.date.fromisoformat(r["fecha"])).days for r in idx.get(a["potrero"], [])]
            mejor = None
            for r, d in zip(idx.get(a["potrero"], []), cands):
                if 0 <= d <= VENTANA[fuente] and (mejor is None or d < mejor[1]):
                    mejor = (r, d)
            if mejor:
                fila[f"{fuente}_dias_antes"] = mejor[1]
                fila.update({k: v for k, v in mejor[0].items() if k.startswith(fuente + "_")})
        if fila.get("s1_vv") not in (None, "") and fila.get("s1_vh") not in (None, ""):
            vv, vh = 10 ** (float(fila["s1_vv"]) / 10), 10 ** (float(fila["s1_vh"]) / 10)
            fila["s1_vh_menos_vv"] = round(float(fila["s1_vh"]) - float(fila["s1_vv"]), 3)
            fila["s1_rvi"] = round(4 * vh / (vv + vh), 4)
        for dias in (7, 15, 30):
            dd = [met.get((fa - dt.timedelta(days=i)).isoformat()) for i in range(1, dias + 1)]
            dd = [x for x in dd if x]
            if len(dd) == dias:
                fila[f"lluvia_{dias}d"] = round(sum(x["precipitation_sum"] or 0 for x in dd), 1)
                fila[f"temp_{dias}d"] = round(sum(x["temperature_2m_mean"] for x in dd) / dias, 2)
                fila[f"et0_{dias}d"] = round(sum(x["et0_fao_evapotranspiration"] or 0 for x in dd), 1)
                fila[f"radiacion_{dias}d"] = round(sum(x["shortwave_radiation_sum"] or 0 for x in dd), 1)
        fila["helada_30d"] = round(min((x["temperature_2m_min"] for x in [met.get((fa - dt.timedelta(days=i)).isoformat()) for i in range(1, 31)] if x), default=float("nan")), 1)
        tabla.append(fila)
    _guardar(carpeta / "dataset_ml.csv", tabla)
    con = {k: sum(1 for t in tabla if t.get(f"{k}_dias_antes") not in (None, "")) for k in fuentes}
    log(f"Tabla: {len(tabla)} aforos | con Sentinel-2 {con['s2']} | radar {con['s1']} | Landsat {con['l']} | "
        f"con clima {sum(1 for t in tabla if 'lluvia_30d' in t)}")


def rangos_de(fechas, antes):
    ventanas = sorted((f - dt.timedelta(days=antes), f) for f in fechas)
    unidos = []
    for a, b in ventanas:
        if unidos and a <= unidos[-1][1] + dt.timedelta(days=1):
            unidos[-1][1] = max(unidos[-1][1], b)
        else:
            unidos.append([a, b])
    return unidos


if __name__ == "__main__":
    with rasterio.Env(**aforo.OPCIONES_GDAL):
        main(sys.argv[1] if len(sys.argv) > 1 else "altos_de_medina", int(sys.argv[2]) if len(sys.argv) > 2 else 90)
