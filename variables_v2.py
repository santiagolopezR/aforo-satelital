"""
Variables del modelo v2, compartidas por la prueba (modelo_v2.py) y la estimación diaria (estimar_hoy.py).

Cada fila es un diccionario con al menos: fecha, potrero, finca (sub-finca), _f (date) y, si es un aforo, _y (kg FV/m²).
"""
import bisect
import collections
import datetime as dt

import numpy as np

CLAVES_S2 = ["s2_ndvi", "s2_ndre1", "s2_ndre2", "s2_ndmi", "s2_evi", "s2_cire", "s2_B8A", "s2_B11", "s2_B05", "s2_B06", "s2_B07"]
CREC_V = ["crec_descanso", "crec_7d", "crec_15d", "f_temp_30d", "f_agua_30d", "agua_suelo"]
MEM_V = ["memoria_ciclo_anterior", "memoria_dias_desde_anterior", "nivel_finca_7d", "nivel_finca_30d"]
MAN_V = ["dias_descanso", "mes"]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def balance_y_crecimiento(met):
    """Índice diario de crecimiento del kikuyo: radiación × f(temperatura) × f(agua del suelo)."""
    agua, cap, salida = 50.0, 100.0, {}
    for d in sorted(met):
        m = met[d]
        t = m.get("temperature_2m_mean")
        lluvia, et0 = m.get("precipitation_sum") or 0.0, m.get("et0_fao_evapotranspiration") or 0.0
        rad = m.get("shortwave_radiation_sum") or 0.0
        agua = min(cap, max(0.0, agua + lluvia - et0))
        f_t = 0.0 if t is None else min(max((t - 8) / (22 - 8), 0.0), 1.0)
        f_w = min(agua / (0.5 * cap), 1.0)
        salida[d] = {"g": rad * f_t * f_w, "f_t": f_t, "f_w": f_w, "agua": agua}
    return salida


def agregar_clima(filas, met):
    """Lluvia, temperatura, evapotranspiración y radiación de los 7/15/30 días previos (como dataset_ml.py)."""
    for r in filas:
        fa = r["_f"]
        for dias in (7, 15, 30):
            dd = [met.get((fa - dt.timedelta(days=i)).isoformat()) for i in range(1, dias + 1)]
            dd = [x for x in dd if x]
            if len(dd) == dias:
                r[f"lluvia_{dias}d"] = round(sum(x["precipitation_sum"] or 0 for x in dd), 1)
                r[f"temp_{dias}d"] = round(sum(x["temperature_2m_mean"] for x in dd) / dias, 2)
                r[f"et0_{dias}d"] = round(sum(x["et0_fao_evapotranspiration"] or 0 for x in dd), 1)
                r[f"radiacion_{dias}d"] = round(sum(x["shortwave_radiation_sum"] or 0 for x in dd), 1)
        mins = [met.get((fa - dt.timedelta(days=i)).isoformat()) for i in range(1, 31)]
        mins = [x["temperature_2m_min"] for x in mins if x and x.get("temperature_2m_min") is not None]
        r["helada_30d"] = round(min(mins), 1) if mins else np.nan


def agregar_crecimiento(filas, crec):
    def suma(fin, dias, k="g"):
        vals = [crec.get((fin - dt.timedelta(days=i)).isoformat(), {}).get(k) for i in range(1, dias + 1)]
        vals = [v for v in vals if v is not None]
        return sum(vals) if len(vals) >= 0.8 * dias else np.nan
    for r in filas:
        descanso = num(r.get("dias_descanso"))
        d = int(descanso) if np.isfinite(descanso) and 10 <= descanso <= 90 else 30
        r["crec_descanso"] = suma(r["_f"], d)
        r["crec_7d"] = suma(r["_f"], 7)
        r["crec_15d"] = suma(r["_f"], 15)
        r["f_temp_30d"] = suma(r["_f"], 30, "f_t") / 30
        r["f_agua_30d"] = suma(r["_f"], 30, "f_w") / 30
        ayer = crec.get((r["_f"] - dt.timedelta(days=1)).isoformat())
        r["agua_suelo"] = ayer["agua"] if ayer else np.nan


def agregar_memoria(filas, aforos):
    """
    Memoria de cada fila a partir de los aforos ANTERIORES a su fecha:
    aforo del ciclo anterior de la franja y nivel de los otros potreros de la sub-finca en 7 y 30 días.
    """
    por_pot, por_fin = collections.defaultdict(list), collections.defaultdict(list)
    for a in sorted(aforos, key=lambda a: a["_f"]):
        por_pot[a["potrero"]].append((a["_f"], a["_y"]))
        por_fin[a["finca"]].append((a["_f"], a["potrero"], a["_y"]))
    fechas_fin = {k: [x[0] for x in v] for k, v in por_fin.items()}
    for r in filas:
        previos = [x for x in por_pot.get(r["potrero"], []) if x[0] < r["_f"]]
        r["memoria_ciclo_anterior"] = previos[-1][1] if previos else np.nan
        r["memoria_dias_desde_anterior"] = (r["_f"] - previos[-1][0]).days if previos else np.nan
        lst, fechas = por_fin.get(r["finca"], []), fechas_fin.get(r["finca"], [])
        for dias in (7, 30):
            i0 = bisect.bisect_left(fechas, r["_f"] - dt.timedelta(days=dias))
            i1 = bisect.bisect_left(fechas, r["_f"])
            vals = [y for f, p, y in lst[i0:i1] if p != r["potrero"]]
            r[f"nivel_finca_{dias}d"] = float(np.mean(vals)) if len(vals) >= 3 else np.nan


def escenas_por_finca(obs_s2, finca_de):
    """Promedio de cada sub-finca en cada imagen Sentinel-2 (con 5+ franjas despejadas)."""
    img = collections.defaultdict(list)
    for o in obs_s2:
        fin = finca_de.get(o["potrero"]) or o["potrero"].rsplit(" ", 1)[0]
        img[(fin, o["fecha"])].append(o)
    escena = {}
    for (fin, fe), lst in img.items():
        if len(lst) >= 5:
            escena.setdefault(fin, []).append((dt.date.fromisoformat(fe), {k: np.nanmean([num(o.get(k)) for o in lst]) for k in CLAVES_S2}))
    for lst in escena.values():
        lst.sort(key=lambda x: x[0])
    return escena


def agregar_s2_finca(filas, escena):
    for r in filas:
        cands = [(r["_f"] - fe, v) for fe, v in escena.get(r["finca"], []) if 0 <= (r["_f"] - fe).days <= 15]
        if cands:
            edad, v = min(cands, key=lambda c: c[0])
            r["s2finca_dias_antes"] = edad.days
            for k in CLAVES_S2:
                r[k.replace("s2_", "s2finca_")] = v[k]
                if r.get(k) not in (None, ""):
                    r[k.replace("s2_", "s2dif_")] = num(r[k]) - v[k]   # franja respecto a su sub-finca


def grupos(cols):
    clima_v = [c for c in cols if c.split("_")[0] in ("lluvia", "temp", "et0", "radiacion") or c == "helada_30d"]
    sat_v = [c for c in cols if c.startswith(("s2_", "s1_", "l_", "s2finca_", "s2dif_")) and c != "s2_pixeles"]
    return {
        "sin_satelite": CREC_V + clima_v + MEM_V + MAN_V,
        "completo": CREC_V + clima_v + MEM_V + MAN_V + sat_v,
        "sin_memoria": CREC_V + clima_v + MAN_V + sat_v,
        "solo_satelite": sat_v,
    }
