"""
Entrena y prueba modelos de machine learning para estimar el aforo (kg FV/m²) desde satélite.

    python entrenar_ml.py altos_de_medina

Prueba honesta: se entrena con los aforos hasta 2024 y se mide con los de 2025–2026 (el modelo no los vio).
También se valida dejando potreros completos por fuera (que no aprenda "este potrero siempre da X").

Compara grupos de variables para saber qué sirve de verdad:
  solo NDVI · todas las bandas de Sentinel-2 · + radar, Landsat y clima · + manejo · solo radar y clima ·
  y un control SIN satélite (solo clima y manejo): si ese control gana, el satélite no está aportando.

Salidas: docs/<finca>/modelo_ml.json (para la página docs/modelo/) y data/<finca>/modelo_*.joblib.
"""
import csv
import datetime as dt
import json
import math
import sys

import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
from sklearn.model_selection import GroupKFold

import aforo

CORTE = "2025-01-01"   # entrenar antes, probar desde esta fecha

NOMBRES = {
    "s2_ndvi": "NDVI", "s2_evi": "EVI", "s2_savi": "SAVI", "s2_msavi": "MSAVI", "s2_gndvi": "NDVI verde",
    "s2_ndre1": "Borde rojo 1 (NDRE)", "s2_ndre2": "Borde rojo 2", "s2_cire": "Clorofila borde rojo",
    "s2_ndmi": "Humedad (NDMI, onda corta)", "s2_nbr": "NBR (onda corta 2)", "s2_ndwi": "Agua (NDWI)",
    "s2_stvi": "Onda corta 1/2 (seco vs verde)", "s2_B02": "Azul", "s2_B03": "Verde", "s2_B04": "Rojo",
    "s2_B05": "Borde rojo B5", "s2_B06": "Borde rojo B6", "s2_B07": "Borde rojo B7", "s2_B08": "Infrarrojo B8",
    "s2_B8A": "Infrarrojo B8A", "s2_B11": "Onda corta B11", "s2_B12": "Onda corta B12",
    "s1_vv": "Radar VV", "s1_vh": "Radar VH", "s1_vh_menos_vv": "Radar VH−VV", "s1_rvi": "Radar RVI",
    "l_ndvi": "Landsat NDVI", "l_ndmi": "Landsat NDMI", "l_nbr": "Landsat NBR", "l_swir1": "Landsat onda corta 1",
    "l_swir2": "Landsat onda corta 2", "l_temp_c": "Landsat temperatura superficie",
    "dias_descanso": "Días de descanso", "mes": "Mes",
}


def nombre(v):
    if v in NOMBRES:
        return NOMBRES[v]
    for pre, txt in (("lluvia_", "Lluvia "), ("temp_", "Temperatura "), ("et0_", "Evapotranspiración "),
                     ("radiacion_", "Radiación ")):
        if v.startswith(pre):
            return txt + v[len(pre):].replace("d", " días")
    if v == "helada_30d":
        return "Mínima 30 días"
    for base, txt in NOMBRES.items():
        if v.startswith(base + "_"):
            return txt + " " + v[len(base) + 1:]
    return v


def metricas(y, p):
    y, p = np.asarray(y, float), np.asarray(p, float)
    pos = y > 0
    sse = float(((y - p) ** 2).sum())
    sst = float(((y - y.mean()) ** 2).sum()) or 1e-9
    return {"n": int(len(y)), "r2": round(1 - sse / sst, 3), "rmse": round(math.sqrt(sse / len(y)), 3),
            "mae": round(float(np.abs(y - p).mean()), 3), "mae_pct": round(float((np.abs(y[pos] - p[pos]) / y[pos]).mean() * 100), 1),
            "sesgo": round(float((p - y).mean()), 3)}


def modelo():
    return HistGradientBoostingRegressor(max_iter=400, learning_rate=0.04, max_leaf_nodes=15,
                                         min_samples_leaf=25, l2_regularization=1.0, random_state=0)


def main(slug):
    cb = json.loads((aforo.RAIZ / "config.json").read_text(encoding="utf-8"))
    f = aforo.fincas_de(cb)[slug]
    aforo.usar_finca(slug, f)
    ruta = aforo.HISTORIAL.parent / "dataset_ml.csv"
    filas = list(csv.DictReader(ruta.open(encoding="utf-8")))
    cols = sorted({k for r in filas for k in r})
    num = lambda v: float(v) if v not in (None, "") else np.nan

    s2 = [c for c in cols if c.startswith("s2_") and c not in ("s2_dias_antes", "s2_pixeles")]
    s1 = [c for c in cols if c.startswith("s1_") and c != "s1_dias_antes"]
    ls = [c for c in cols if c.startswith("l_") and c != "l_dias_antes"]
    met = [c for c in cols if c.split("_")[0] in ("lluvia", "temp", "et0", "radiacion") or c == "helada_30d"]
    man = ["dias_descanso", "mes"]
    tiene = lambda r, fuente: r.get(f"{fuente}_dias_antes") not in (None, "")
    grupos = [
        ("ndvi", "Solo NDVI", ["s2_ndvi"], "s2"),
        ("s2", "Sentinel-2: todas las bandas e índices", s2, "s2"),
        ("optico", "Sentinel-2 + radar + Landsat + clima", s2 + s1 + ls + met, "s2"),
        ("optico_manejo", "Todo lo anterior + manejo (descanso, mes)", s2 + s1 + ls + met + man, "s2"),
        ("radar", "Solo radar + clima (días nublados)", s1 + met, "s1"),
        ("radar_manejo", "Radar + clima + manejo", s1 + met + man, "s1"),
        ("control", "Control sin satélite: clima + manejo", met + man, "s1"),
    ]
    salida = {"finca": aforo.FINCA, "generado": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"),
              "corte": CORTE, "aforos": len(filas), "modelos": [], "importancia": {}, "dispersion": {}}
    for clave, titulo, vars_, req in grupos:
        datos = [r for r in filas if tiene(r, req)]
        X = np.array([[num(r.get(v)) for v in vars_] for r in datos])
        y = np.array([float(r["aforo_kg_fv_m2"]) for r in datos])
        entreno = np.array([r["fecha"] < CORTE for r in datos])
        info = {"clave": clave, "titulo": titulo, "variables": len(vars_), "filas": len(datos),
                "entreno": int(entreno.sum()), "prueba": int((~entreno).sum())}
        if entreno.sum() < 50 or (~entreno).sum() < 30:
            info["nota"] = "pocos datos"
            salida["modelos"].append(info)
            continue
        m = modelo().fit(X[entreno], y[entreno])
        pred = m.predict(X[~entreno])
        info["prueba_2025_2026"] = metricas(y[~entreno], pred)
        info["promedio_como_estimacion"] = metricas(y[~entreno], np.full((~entreno).sum(), y[entreno].mean()))
        # Validación por potreros: cada grupo de potreros se predice con un modelo que nunca los vio.
        grupos_pot = np.array([r["potrero"] for r in datos])
        pv = np.zeros(len(y))
        for tr, te in GroupKFold(n_splits=5).split(X, y, grupos_pot):
            pv[te] = modelo().fit(X[tr], y[tr]).predict(X[te])
        info["potreros_no_vistos"] = metricas(y, pv)
        salida["modelos"].append(info)
        aforo.log(f"{titulo}: prueba 2025-26 {info['prueba_2025_2026']} | potreros no vistos R² {info['potreros_no_vistos']['r2']}")
        if clave in ("optico", "radar", "optico_manejo"):
            imp = permutation_importance(m, X[~entreno], y[~entreno], n_repeats=5, random_state=0, scoring="r2")
            orden = np.argsort(-imp.importances_mean)[:15]
            salida["importancia"][clave] = [{"var": vars_[i], "nombre": nombre(vars_[i]),
                                              "valor": round(float(imp.importances_mean[i]), 4)} for i in orden]
            salida["dispersion"][clave] = [{"f": r["fecha"], "p": r["potrero"], "real": float(yy), "est": round(float(pp), 3)}
                                           for r, yy, pp in zip([d for d, e in zip(datos, entreno) if not e], y[~entreno], pred)]
        if clave in ("optico", "radar"):
            # Modelo final con todos los años, para usarlo en el mapa.
            final = modelo().fit(X, y)
            joblib.dump({"modelo": final, "variables": vars_, "titulo": titulo},
                        aforo.HISTORIAL.parent / f"modelo_{clave}.joblib")
    destino = aforo.SALIDA_WEB.with_name("modelo_ml.json")
    destino.write_text(json.dumps(salida, ensure_ascii=False), encoding="utf-8")
    aforo.log(f"Listo: {destino}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "altos_de_medina")
