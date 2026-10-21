"""
Modelo v2 del aforo (según la literatura: crecimiento + memoria + satélite + recalibración continua).

    python modelo_v2.py altos_de_medina

Ingredientes por aforo:
  1. Crecimiento del pasto durante el descanso: cada día crece según radiación × temperatura × humedad del suelo
     (balance lluvia − evapotranspiración). Kikuyo: base 8 °C, óptimo 22 °C.
  2. Memoria: aforo del ciclo anterior de esa franja y nivel reciente de la finca (otros potreros, 7 y 30 días).
  3. Satélite en dos escalas: la franja (último dato de cada sensor con su edad) y el promedio de su sub-finca
     en la última imagen despejada, más la diferencia franja − sub-finca.
  4. Recalibración mensual: cada mes se entrena con todo lo anterior y se predice ese mes (como en la finca real).
Se mide por franja y por sub-finca y semana. Resultado en docs/<finca>/modelo_v2.json (página docs/modelo/).
"""
import collections
import csv
import datetime as dt
import json
import math
import sys

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor

import aforo
import dataset_ml

DESDE_PRUEBA = "2024-01"
CLAVES_S2 = ["s2_ndvi", "s2_ndre1", "s2_ndre2", "s2_ndmi", "s2_evi", "s2_cire", "s2_B8A", "s2_B11", "s2_B05", "s2_B06", "s2_B07"]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def balance_y_crecimiento(met):
    """Índice diario de crecimiento del kikuyo a partir del clima (sin unidades; el modelo lo escala)."""
    dias = sorted(met)
    agua, cap, salida = 50.0, 100.0, {}
    for d in dias:
        m = met[d]
        t = m.get("temperature_2m_mean")
        lluvia, et0 = m.get("precipitation_sum") or 0.0, m.get("et0_fao_evapotranspiration") or 0.0
        rad = m.get("shortwave_radiation_sum") or 0.0
        agua = min(cap, max(0.0, agua + lluvia - et0))
        f_t = 0.0 if t is None else min(max((t - 8) / (22 - 8), 0.0), 1.0)
        f_w = min(agua / (0.5 * cap), 1.0)
        salida[d] = {"g": rad * f_t * f_w, "f_t": f_t, "f_w": f_w, "agua": agua}
    return salida


def main(slug):
    cb = json.loads((aforo.RAIZ / "config.json").read_text(encoding="utf-8"))
    f = aforo.fincas_de(cb)[slug]
    aforo.usar_finca(slug, f)
    config = aforo.config_de_finca(cb, f)
    carpeta = aforo.HISTORIAL.parent
    filas = list(csv.DictReader((carpeta / "dataset_ml.csv").open(encoding="utf-8")))
    for r in filas:
        r["_f"] = dt.date.fromisoformat(r["fecha"])
        r["_y"] = float(r["aforo_kg_fv_m2"])
    filas.sort(key=lambda r: r["_f"])

    # --- 1. Crecimiento con clima diario
    potreros = aforo.cargar_potreros(config)
    met = dataset_ml.clima(potreros, filas[0]["_f"] - dt.timedelta(days=90), min(filas[-1]["_f"], dt.date.today() - dt.timedelta(days=6)))
    crec = balance_y_crecimiento(met)

    def suma(fin, dias, k="g"):
        vals = [crec.get((fin - dt.timedelta(days=i)).isoformat(), {}).get(k) for i in range(1, dias + 1)]
        vals = [v for v in vals if v is not None]
        return sum(vals) if len(vals) >= 0.8 * dias else np.nan

    # --- 2. Memoria: ciclo anterior de la franja y nivel reciente de la sub-finca
    previo = {}
    por_finca = collections.defaultdict(list)
    for r in filas:
        ant = previo.get(r["potrero"])
        r["memoria_ciclo_anterior"] = ant[1] if ant else np.nan
        r["memoria_dias_desde_anterior"] = (r["_f"] - ant[0]).days if ant else np.nan
        previo[r["potrero"]] = (r["_f"], r["_y"])
        por_finca[r["finca"]].append(r)
    for fin, lst in por_finca.items():
        for r in lst:
            for dias in (7, 30):
                vals = [x["_y"] for x in lst if x["potrero"] != r["potrero"] and 1 <= (r["_f"] - x["_f"]).days <= dias]
                r[f"nivel_finca_{dias}d"] = float(np.mean(vals)) if len(vals) >= 3 else np.nan
    for r in filas:
        descanso = num(r.get("dias_descanso"))
        d = int(descanso) if np.isfinite(descanso) and 10 <= descanso <= 90 else 30
        r["crec_descanso"] = suma(r["_f"], d)
        r["crec_7d"] = suma(r["_f"], 7)
        r["crec_15d"] = suma(r["_f"], 15)
        r["f_temp_30d"] = suma(r["_f"], 30, "f_t") / 30
        r["f_agua_30d"] = suma(r["_f"], 30, "f_w") / 30
        hoy = crec.get((r["_f"] - dt.timedelta(days=1)).isoformat())
        r["agua_suelo"] = hoy["agua"] if hoy else np.nan

    # --- 3. Satélite a escala de sub-finca: promedio de sus franjas en cada imagen despejada
    s2 = list(csv.DictReader((carpeta / "ml_s2.csv").open(encoding="utf-8")))
    finca_de = {r["potrero"]: r["finca"] for r in filas}
    img = collections.defaultdict(list)
    for o in s2:
        fin = finca_de.get(o["potrero"]) or o["potrero"].rsplit(" ", 1)[0]
        img[(fin, o["fecha"])].append(o)
    escena_finca = {}
    for (fin, fe), lst in img.items():
        if len(lst) >= 5:
            escena_finca.setdefault(fin, []).append((dt.date.fromisoformat(fe), {k: np.nanmean([num(o.get(k)) for o in lst]) for k in CLAVES_S2}))
    for lst in escena_finca.values():
        lst.sort(key=lambda x: x[0])
    for r in filas:
        cands = [(r["_f"] - fe, v) for fe, v in escena_finca.get(r["finca"], []) if 0 <= (r["_f"] - fe).days <= 15]
        if cands:
            edad, v = min(cands, key=lambda c: c[0])
            r["s2finca_dias_antes"] = edad.days
            for k in CLAVES_S2:
                r[k.replace("s2_", "s2finca_")] = v[k]
                if r.get(k) not in (None, ""):
                    r[k.replace("s2_", "s2dif_")] = num(r[k]) - v[k]   # franja respecto a su sub-finca

    # --- Grupos de variables
    cols = sorted({k for r in filas for k in r if not k.startswith("_")})
    clima_v = [c for c in cols if c.split("_")[0] in ("lluvia", "temp", "et0", "radiacion") or c == "helada_30d"]
    crec_v = ["crec_descanso", "crec_7d", "crec_15d", "f_temp_30d", "f_agua_30d", "agua_suelo"]
    mem_v = ["memoria_ciclo_anterior", "memoria_dias_desde_anterior", "nivel_finca_7d", "nivel_finca_30d"]
    man_v = ["dias_descanso", "mes"]
    sat_v = [c for c in cols if c.startswith(("s2_", "s1_", "l_", "s2finca_", "s2dif_")) and c != "s2_pixeles"]
    configuraciones = [
        ("nivel", "Referencia: nivel reciente de la sub-finca (30 días)", None),
        ("sin_satelite", "Crecimiento + clima + memoria + manejo (sin satélite)", crec_v + clima_v + mem_v + man_v),
        ("completo", "Lo anterior + satélite (franja y sub-finca)", crec_v + clima_v + mem_v + man_v + sat_v),
        ("solo_satelite", "Solo satélite (franja y sub-finca)", sat_v),
    ]

    def modelo():
        return HistGradientBoostingRegressor(max_iter=300, learning_rate=0.04, max_leaf_nodes=15,
                                             min_samples_leaf=30, l2_regularization=1.0, random_state=0)

    def metr(o, p):
        o, p = np.asarray(o), np.asarray(p)
        return {"n": int(len(o)), "r2": round(float(1 - ((o - p) ** 2).sum() / ((o - o.mean()) ** 2).sum()), 3),
                "rmse": round(float(np.sqrt(((o - p) ** 2).mean())), 3),
                "mae_pct": round(float(np.mean(np.abs(o - p) / o) * 100), 1)}

    meses = sorted({r["fecha"][:7] for r in filas if r["fecha"][:7] >= DESDE_PRUEBA})
    salida = {"finca": aforo.FINCA, "generado": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"),
              "desde_prueba": DESDE_PRUEBA, "aforos": len(filas), "configuraciones": [], "serie_semanal": {}}
    for clave, titulo, vars_ in configuraciones:
        obs, pred, ref = [], [], []
        if vars_ is None:
            for r in filas:
                if r["fecha"][:7] >= DESDE_PRUEBA and np.isfinite(r["nivel_finca_30d"]):
                    obs.append(r["_y"]); pred.append(r["nivel_finca_30d"]); ref.append(r)
        else:
            X = np.array([[num(r.get(v)) for v in vars_] for r in filas])
            y = np.array([r["_y"] for r in filas])
            mes = np.array([r["fecha"][:7] for r in filas])
            for m in meses:
                tr, te = mes < m, mes == m
                if te.sum() == 0:
                    continue
                p = modelo().fit(X[tr], y[tr]).predict(X[te])
                obs += list(y[te]); pred += list(p); ref += [r for r, t in zip(filas, te) if t]
        franja = metr(obs, pred)
        # Por sub-finca y semana: promedio de aforos y de estimaciones de esa semana (3+ aforos)
        sem = collections.defaultdict(lambda: ([], []))
        for r, o, p in zip(ref, obs, pred):
            k = (r["finca"], r["_f"].isocalendar()[0], r["_f"].isocalendar()[1])
            sem[k][0].append(o); sem[k][1].append(p)
        so = [np.mean(a) for a, b in sem.values() if len(a) >= 3]
        sp = [np.mean(b) for a, b in sem.values() if len(a) >= 3]
        semana = metr(so, sp)
        salida["configuraciones"].append({"clave": clave, "titulo": titulo, "franja": franja, "finca_semana": semana})
        aforo.log(f"{titulo}: franja R² {franja['r2']} (error {franja['rmse']} kg, {franja['mae_pct']} %) | "
                  f"sub-finca/semana R² {semana['r2']} (error {semana['rmse']} kg, n={semana['n']})")
        if clave in ("completo", "sin_satelite"):
            puntos = sorted(((k[0], f"{k[1]}-{k[2]:02d}", float(np.mean(a)), float(np.mean(b)))
                             for k, (a, b) in sem.items() if len(a) >= 3))
            salida["serie_semanal"][clave] = [{"finca": fi, "semana": s, "real": round(o, 3), "est": round(p, 3)}
                                              for fi, s, o, p in puntos]
    aforo.SALIDA_WEB.with_name("modelo_v2.json").write_text(json.dumps(salida, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "altos_de_medina")
