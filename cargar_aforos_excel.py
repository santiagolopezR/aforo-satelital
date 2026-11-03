"""
Carga los aforos de platómetro desde el Excel de Altos de Medina (se corre en el PC, el Excel no se sube).

    python cargar_aforos_excel.py "D:/MAPA Altos de Medina Claude/SOLTAR EXCEL AFOROS"

Toma el .xlsx más reciente de la carpeta, lee la hoja ROTACION-AFORO-CONSUMO y saca SOLO los aforos
(fecha, finca, potrero, kg FV/m², área, vacas, días de ocupación y de rotación, colchón). Los suma a
aforos_campo/altos_de_medina_platometro.csv sin borrar nada (si una fecha y potrero ya estaban, manda el Excel)
y sube únicamente ese archivo a GitHub. El mapa se actualiza solo en unos minutos.
"""
import csv
import datetime as dt
import hashlib
import subprocess
import sys
from pathlib import Path

import openpyxl

RAIZ = Path(__file__).resolve().parent
DESTINO = RAIZ / "aforos_campo" / "altos_de_medina_platometro.csv"
HOJA = "ROTACION-AFORO-CONSUMO"
FINCAS = {"ARRIBA": "MANGA GRANDE", "LA POSADA": "LA POSADA", "PIONEROS": "PIONEROS", "ABAJO": "ABAJO"}
COLUMNAS = ["fecha", "potrero", "kg_fv_m2", "pct_ms", "metodo", "observaciones", "origen"]


def leer_excel(ruta):
    wb = openpyxl.load_workbook(ruta, read_only=True, data_only=True)
    if HOJA not in wb.sheetnames:
        sys.exit(f"El Excel no tiene la hoja '{HOJA}'. Hojas: {wb.sheetnames}")
    filas = wb[HOJA].iter_rows(values_only=True)
    cab = [str(c or "").strip().upper() for c in next(filas)]
    i = {n: cab.index(n) for n in ("FECHA", "FINCA", "POTRERO", "AFORO PLATOMETRO (KG/M2)", "NUMERO VACAS LOTE",
                                     "DIAS OCUPACION", "DIAS ROTACION", "COLCHON CM")}
    area = 4  # columna sin nombre claro ("3") con el área de la franja en m²
    salida, descartadas = [], 0
    for r in filas:
        fecha, finca, potrero, kg = r[i["FECHA"]], r[i["FINCA"]], r[i["POTRERO"]], r[i["AFORO PLATOMETRO (KG/M2)"]]
        if not isinstance(fecha, dt.datetime) or finca is None or potrero in (None, ""):
            continue
        if not isinstance(kg, (int, float)) or not 0.3 < kg <= 5:
            descartadas += 1
            continue
        finca = FINCAS.get(str(finca).strip().upper(), str(finca).strip().upper())
        pot = str(potrero).strip().upper()
        if pot.endswith(".0"):
            pot = pot[:-2]
        extra = {"area_m2": r[area], "vacas": r[i["NUMERO VACAS LOTE"]], "dias_ocupacion": r[i["DIAS OCUPACION"]],
                 "dias_rotacion": r[i["DIAS ROTACION"]], "colchon_cm": r[i["COLCHON CM"]]}
        obs = "; ".join(f"{k}={float(v):g}" for k, v in extra.items() if isinstance(v, (int, float)))
        salida.append({"fecha": fecha.date().isoformat(), "potrero": f"{finca} {pot}", "kg_fv_m2": round(float(kg), 3),
                       "pct_ms": "", "metodo": "Plato medidor", "observaciones": obs,
                       "origen": f"Excel Altos de Medina ({str(r[i['FINCA']]).strip().upper()})"})
    return salida, descartadas


def main(carpeta):
    carpeta = Path(carpeta)
    excels = sorted((p for p in carpeta.glob("*.xls*") if not p.name.startswith("~$")), key=lambda p: p.stat().st_mtime)
    if not excels:
        sys.exit(f"No hay ningún Excel en {carpeta}")
    ruta = excels[-1]
    print(f"Leyendo {ruta.name} ...")
    nuevas, descartadas = leer_excel(ruta)
    print(f"  {len(nuevas)} aforos válidos (fechas {min(n['fecha'] for n in nuevas)} a {max(n['fecha'] for n in nuevas)}); "
          f"{descartadas} lecturas descartadas por fuera de 0,3–5 kg/m².")

    actuales = list(csv.DictReader(DESTINO.open(encoding="utf-8"))) if DESTINO.exists() else []
    clave = lambda r: (r["fecha"], r["potrero"])
    union = {clave(r): r for r in actuales}
    agregadas = sum(1 for n in nuevas if clave(n) not in union)
    for n in nuevas:
        union[clave(n)] = n          # si ya estaba, manda el Excel; lo que no viene en el Excel se conserva
    filas = sorted(union.values(), key=lambda r: (r["fecha"], r["potrero"]))
    with DESTINO.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNAS)
        w.writeheader()
        w.writerows({c: r.get(c, "") for c in COLUMNAS} for r in filas)
    print(f"  Aforos nuevos: {agregadas}. Total en el sistema: {len(filas)} (último: {filas[-1]['fecha']}).")

    git = lambda *a: subprocess.run(["git", "-C", str(RAIZ), *a], capture_output=True, text=True)
    if not git("status", "--porcelain", str(DESTINO.relative_to(RAIZ))).stdout.strip():
        print("No hay cambios que subir: estos aforos ya estaban cargados.")
        return
    huella = hashlib.sha1(ruta.read_bytes()).hexdigest()[:8]
    git("pull", "--rebase", "--autostash")
    git("add", str(DESTINO.relative_to(RAIZ)))
    r = git("commit", "-m", f"Aforos de platómetro de Altos de Medina hasta {filas[-1]['fecha']} ({agregadas} nuevos, Excel {huella})")
    print(r.stdout.strip() or r.stderr.strip())
    r = git("push")
    if r.returncode:
        print("No se pudo subir a GitHub:\n" + r.stderr)
        sys.exit(1)
    print("Listo: aforos subidos. El mapa se actualiza solo en unos minutos.")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "D:/MAPA Altos de Medina Claude/SOLTAR EXCEL AFOROS")
