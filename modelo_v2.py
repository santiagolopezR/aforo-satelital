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
import variables_v2 as v2

DESDE_PRUEBA = "2024-01"


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

    # --- Variables (ver variables_v2.py): crecimiento, memoria y satélite por sub-finca
    potreros = aforo.cargar_potreros(config)
    met = dataset_ml.clima(potreros, filas[0]["_f"] - dt.timedelta(days=90), min(filas[-1]["_f"], dt.date.today() - dt.timedelta(days=6)))
    v2.agregar_crecimiento(filas, v2.balance_y_crecimiento(met))
    v2.agregar_memoria(filas, filas)
    finca_de = {r["potrero"]: r["finca"] for r in filas}
    v2.agregar_s2_finca(filas, v2.escenas_por_finca(list(csv.DictReader((carpeta / "ml_s2.csv").open(encoding="utf-8"))), finca_de))
    cols = sorted({k for r in filas for k in r if not k.startswith("_")})
    g = v2.grupos(cols)
    configuraciones = [
        ("nivel", "Referencia: nivel reciente de la sub-finca (30 días)", None),
        ("sin_satelite", "Crecimiento + clima + memoria + manejo (sin satélite)", g["sin_satelite"]),
        ("completo", "Lo anterior + satélite (franja y sub-finca)", g["completo"]),
        ("sin_memoria", "Crecimiento + clima + manejo + satélite (sin aforos recientes)", g["sin_memoria"]),
        ("solo_satelite", "Solo satélite (franja y sub-finca)", g["solo_satelite"]),
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
            X = np.array([[v2.num(r.get(v)) for v in vars_] for r in filas])
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
