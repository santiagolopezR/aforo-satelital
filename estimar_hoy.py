"""
Estimación diaria del aforo de cada franja con el modelo v2 (para el mapa).

    python estimar_hoy.py altos_de_medina

  1. Entrena el modelo v2 con todos los aforos de campo (versión "con memoria" y "sin memoria").
  2. Para hoy toma lo último de cada satélite (Sentinel-2 hasta 20 días, radar 12, Landsat 16, con su edad),
     el clima de los últimos 3 meses hasta ayer y la memoria de aforos recientes si existe.
  3. Si la sub-finca tiene aforos en los últimos 30 días usa el modelo con memoria; si no, el sin memoria.
Salida: docs/<finca>/estimacion_hoy.json (el mapa la superpone a los datos de la finca).
"""
import csv
import datetime as dt
import json
import sys

import numpy as np
import planetary_computer
import pystac_client
import rasterio
import requests
from shapely.ops import unary_union
from sklearn.ensemble import HistGradientBoostingRegressor

import aforo
import dataset_ml
import variables_v2 as v2

VENTANA = dataset_ml.VENTANA   # s2: 20, s1: 12, l: 16 días


def clima_reciente(potreros, dias=92):   # 92 es el máximo que permite Open-Meteo
    """Clima diario de los últimos días hasta ayer (Open-Meteo, mismo formato que dataset_ml.clima)."""
    c = unary_union([p["geom"] for p in potreros]).centroid
    r = requests.get("https://api.open-meteo.com/v1/forecast", timeout=60, params={
        "latitude": round(c.y, 3), "longitude": round(c.x, 3), "past_days": dias, "forecast_days": 1,
        "daily": "precipitation_sum,temperature_2m_mean,temperature_2m_min,et0_fao_evapotranspiration,shortwave_radiation_sum",
        "timezone": "America/Bogota"})
    r.raise_for_status()
    d = r.json()["daily"]
    hoy = dt.date.today().isoformat()
    return {f: {k: d[k][i] for k in d if k != "time"} for i, f in enumerate(d["time"]) if f < hoy}


def ultimos_satelites(potreros, config, exclusiones, carpeta):
    """Último dato de cada sensor por potrero (con caché en data/<finca>/ml_hoy_*.csv)."""
    hoy = dt.date.today()
    bbox = list(unary_union([p["geom"] for p in potreros]).bounds)
    cat = pystac_client.Client.open(aforo.STAC_URL, modifier=planetary_computer.sign_inplace)
    fuentes = {
        "s2": ("sentinel-2-l2a", {"eo:cloud_cover": {"lt": config["max_nubes_escena"]}},
               lambda it: dataset_ml.features_s2(it, potreros, config, exclusiones)),
        "s1": ("sentinel-1-rtc", None,
               lambda it: [{"fecha": r["fecha"], "potrero": r["potrero"], "escena": r["escena"],
                            "s1_vv": r["vv_db"], "s1_vh": r["vh_db"]}
                           for r in aforo.procesar_s1(it, potreros, config, exclusiones, "")]),
        "l": ("landsat-c2-l2", {"eo:cloud_cover": {"lt": 80}, "platform": {"in": ["landsat-8", "landsat-9"]}},
              lambda it: dataset_ml.features_landsat(it, potreros, config, exclusiones)),
    }
    ultimo = {}
    for clave, (col, filtro, fn) in fuentes.items():
        ruta = carpeta / f"ml_hoy_{clave}.csv"
        cache = dataset_ml._cargar(ruta)
        hechas = {r["escena"] for r in cache}
        marca = ruta.with_suffix(".revisadas")
        hechas |= set(marca.read_text(encoding="utf-8").split()) if marca.exists() else set()
        items = aforo.buscar(cat, col, bbox, hoy - dt.timedelta(days=VENTANA[clave]), hoy, filtro)
        for it in sorted(items, key=lambda i: i.datetime):
            if it.id in hechas:
                continue
            try:
                cache += fn(it)
            except Exception as e:
                aforo.log(f"  ! {it.id}: {e}")
                continue
            hechas.add(it.id)
        dataset_ml._guardar(ruta, cache)
        marca.write_text("\n".join(sorted(hechas)), encoding="utf-8")
        for r in cache:
            edad = (hoy - dt.date.fromisoformat(r["fecha"])).days
            if 0 <= edad <= VENTANA[clave]:
                prev = ultimo.get((clave, r["potrero"]))
                if prev is None or edad < prev[0]:
                    ultimo[(clave, r["potrero"])] = (edad, r)
        aforo.log(f"{clave}: {len(items)} escenas recientes, {sum(1 for k in ultimo if k[0] == clave)} potreros con dato")
    return ultimo


def main(slug):
    cb = json.loads((aforo.RAIZ / "config.json").read_text(encoding="utf-8"))
    f = aforo.fincas_de(cb)[slug]
    aforo.usar_finca(slug, f)
    config = aforo.config_de_finca(cb, f)
    carpeta = aforo.HISTORIAL.parent
    potreros = aforo.cargar_potreros(config)
    exclusiones = aforo.cargar_exclusiones()

    # --- Entrenamiento con todos los aforos (mismas variables que la prueba del modelo v2)
    aforos = list(csv.DictReader((carpeta / "dataset_ml.csv").open(encoding="utf-8")))
    for r in aforos:
        r["_f"], r["_y"] = dt.date.fromisoformat(r["fecha"]), float(r["aforo_kg_fv_m2"])
    aforos.sort(key=lambda r: r["_f"])
    met_hist = dataset_ml.clima(potreros, aforos[0]["_f"] - dt.timedelta(days=90),
                                min(aforos[-1]["_f"], dt.date.today() - dt.timedelta(days=6)))
    v2.agregar_crecimiento(aforos, v2.balance_y_crecimiento(met_hist))
    v2.agregar_memoria(aforos, aforos)
    finca_de = {p["nombre"]: p.get("rotacion") or p["nombre"].rsplit(" ", 1)[0] for p in potreros}
    obs_s2 = dataset_ml._cargar(carpeta / "ml_s2.csv")
    v2.agregar_s2_finca(aforos, v2.escenas_por_finca(obs_s2, finca_de))

    # --- Filas de hoy, una por potrero
    hoy = dt.date.today()
    ultimo_aforo = {}
    for a in aforos:
        ultimo_aforo[a["potrero"]] = a["_f"]
    hoy_filas = []
    for p in potreros:
        r = {"fecha": hoy.isoformat(), "potrero": p["nombre"], "finca": finca_de[p["nombre"]], "_f": hoy, "mes": hoy.month}
        if p["nombre"] in ultimo_aforo:
            d = (hoy - ultimo_aforo[p["nombre"]]).days - 1   # días desde que salieron (ocupación de ~1 día)
            r["dias_descanso"] = d if 10 <= d <= 90 else ""
        hoy_filas.append(r)
    met = clima_reciente(potreros)
    v2.agregar_clima(hoy_filas, met)
    v2.agregar_crecimiento(hoy_filas, v2.balance_y_crecimiento(met))
    v2.agregar_memoria(hoy_filas, aforos)

    ultimo = ultimos_satelites(potreros, config, exclusiones, carpeta)
    for r in hoy_filas:
        for clave in ("s2", "s1", "l"):
            dato = ultimo.get((clave, r["potrero"]))
            if dato:
                r[f"{clave}_dias_antes"] = dato[0]
                r.update({k: v for k, v in dato[1].items() if k.startswith(clave + "_")})
        if r.get("s1_vv") not in (None, "") and r.get("s1_vh") not in (None, ""):
            vv, vh = 10 ** (float(r["s1_vv"]) / 10), 10 ** (float(r["s1_vh"]) / 10)
            r["s1_vh_menos_vv"] = round(float(r["s1_vh"]) - float(r["s1_vv"]), 3)
            r["s1_rvi"] = round(4 * vh / (vv + vh), 4)
    recientes_s2 = [r for (c, _), (e, r) in ultimo.items() if c == "s2"]
    v2.agregar_s2_finca(hoy_filas, v2.escenas_por_finca(recientes_s2, finca_de))

    # --- Dos modelos finales y predicción
    cols = sorted({k for r in aforos for k in r if not k.startswith("_")})
    g = v2.grupos(cols)
    y = np.array([r["_y"] for r in aforos])
    modelos = {}
    for clave in ("completo", "sin_memoria"):
        X = np.array([[v2.num(r.get(v)) for v in g[clave]] for r in aforos])
        modelos[clave] = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.04, max_leaf_nodes=15,
                                                       min_samples_leaf=30, l2_regularization=1.0, random_state=0).fit(X, y)
    precision = {}
    ruta_v2 = aforo.SALIDA_WEB.with_name("modelo_v2.json")
    if ruta_v2.exists():
        precision = {c["clave"]: {"franja": c["franja"], "finca_semana": c["finca_semana"]}
                     for c in json.loads(ruta_v2.read_text(encoding="utf-8"))["configuraciones"]}
    salida = {"finca": aforo.FINCA, "fecha": hoy.isoformat(),
              "generado": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"),
              "ultimo_aforo_campo": aforos[-1]["fecha"], "precision": precision, "potreros": {}, "subfincas": []}
    por_fin = {}
    for r in hoy_filas:
        clave = "completo" if np.isfinite(v2.num(r.get("nivel_finca_30d"))) else "sin_memoria"
        X = np.array([[v2.num(r.get(v)) for v in g[clave]]])
        kg = float(modelos[clave].predict(X)[0])
        edades = {c: r.get(f"{c}_dias_antes") for c in ("s2", "s1", "l") if r.get(f"{c}_dias_antes") not in (None, "")}
        salida["potreros"][r["potrero"]] = {"kg": round(kg, 2), "modelo": clave, "edades": edades,
                                            "dias_descanso": r.get("dias_descanso", "")}
        por_fin.setdefault(r["finca"], []).append((kg, clave))
    for fin, lst in sorted(por_fin.items()):
        salida["subfincas"].append({"finca": fin, "kg": round(float(np.mean([k for k, _ in lst])), 2), "n": len(lst),
                                    "modelo": max({c for _, c in lst}, key=lambda c: sum(1 for _, x in lst if x == c))})
    aforo.SALIDA_WEB.with_name("estimacion_hoy.json").write_text(json.dumps(salida, ensure_ascii=False), encoding="utf-8")
    aforo.log("Estimación de hoy: " + ", ".join(f"{s['finca']} {s['kg']} kg ({s['modelo']})" for s in salida["subfincas"]))


if __name__ == "__main__":
    with rasterio.Env(**aforo.OPCIONES_GDAL):
        main(sys.argv[1] if len(sys.argv) > 1 else "altos_de_medina")
