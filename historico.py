"""
Trae las imágenes Sentinel-2 viejas de las fechas en que hay aforos de campo, para calibrar.

    python historico.py altos_de_medina [minutos]

Para cada aforo busca imágenes desde 'dias_par_antes' días antes hasta 'dias_par_despues' días después,
las procesa igual que aforo.py (por píxel, sin árboles ni bordes) y las guarda en
data/<finca>/historial_largo.csv. Es incremental: las escenas ya revisadas no se vuelven a bajar, y se
detiene a los N minutos guardando lo hecho (la siguiente corrida sigue donde quedó).
"""
import csv
import datetime as dt
import json
import sys
import time

import planetary_computer
import pystac_client
import rasterio
from shapely.ops import unary_union

import aforo

aforo.guardar_pixeles = lambda *a, **k: None  # las imágenes viejas no van al mapa de detalle
log = aforo.log


def rangos(fechas, antes, despues):
    """Une las ventanas de todos los aforos en el menor número de rangos de fechas."""
    ventanas = sorted((f - dt.timedelta(days=antes), f + dt.timedelta(days=despues)) for f in fechas)
    unidos = []
    for a, b in ventanas:
        if unidos and a <= unidos[-1][1] + dt.timedelta(days=1):
            unidos[-1][1] = max(unidos[-1][1], b)
        else:
            unidos.append([a, b])
    return unidos


def main(slug, minutos):
    inicio = time.time()
    config_base = json.loads((aforo.RAIZ / "config.json").read_text(encoding="utf-8"))
    f = aforo.fincas_de(config_base)[slug]
    aforo.usar_finca(slug, f)
    config = aforo.config_de_finca(config_base, f)
    potreros = aforo.cargar_potreros(config)
    if not aforo.MASCARA_ARBOLES.exists():
        log("Falta la máscara de árboles: corre primero aforo.py")
        return
    exclusiones, huella = aforo.cargar_exclusiones(), aforo.huella_mascara(config)

    aforos = aforo.leer_aforos(potreros)
    fechas = sorted({dt.date.fromisoformat(a["fecha"]) for a in aforos})
    # Para el modelo de machine learning se traen también imágenes más viejas (el modelo usa su edad).
    antes = config.get("dias_historico_antes", config.get("dias_par_antes", config["dias_par_aforo_imagen"]))
    despues = config.get("dias_par_despues", config["dias_par_aforo_imagen"])
    lista_rangos = rangos(fechas, antes, despues)
    log(f"{f['nombre']}: {len(aforos)} aforos en {len(fechas)} fechas -> {len(lista_rangos)} rangos de búsqueda")

    # Lo ya hecho: filas guardadas y escenas revisadas (también las nubladas, que no dejan filas).
    todas = []
    if aforo.HISTORIAL_LARGO.exists():
        with aforo.HISTORIAL_LARGO.open(encoding="utf-8") as fh:
            todas = list(csv.DictReader(fh))
    previas = [r for r in todas if r.get("mascara") == huella]
    revisadas_ruta = aforo.HISTORIAL_LARGO.with_name("escenas_revisadas.txt")
    revisadas = set(revisadas_ruta.read_text(encoding="utf-8").split()) if revisadas_ruta.exists() else set()
    if len(previas) < len(todas):
        log("Cambió la máscara de árboles: se vuelven a revisar todas las escenas.")
        revisadas = set()

    bbox = list(unary_union([p["geom"] for p in potreros]).bounds)
    cat = pystac_client.Client.open(aforo.STAC_URL, modifier=planetary_computer.sign_inplace)
    items = {}
    for a, b in lista_rangos:
        for it in aforo.buscar(cat, "sentinel-2-l2a", bbox, a, b, {"eo:cloud_cover": {"lt": config["max_nubes_escena"]}}):
            items[it.id] = it
    pendientes = sorted((it for i, it in items.items() if i not in revisadas), key=lambda i: i.datetime)
    log(f"Escenas Sentinel-2: {len(items)} en total, {len(pendientes)} por revisar")

    nuevas, hechas = [], 0
    for it in pendientes:
        if time.time() - inicio > minutos * 60:
            log(f"Tiempo cumplido ({minutos} min): se guarda y la próxima corrida sigue.")
            break
        try:
            filas = aforo.procesar_s2(it, potreros, config, exclusiones, huella)
        except Exception as e:
            log(f"  ! {it.id}: {e}")
            continue
        nuevas += filas
        revisadas.add(it.id)
        hechas += 1
        if filas:
            log(f"  {it.datetime.date()} nubes {it.properties.get('eo:cloud_cover', 0):.0f}% -> {len(filas)} potreros")

    filas = previas + nuevas
    aforo.HISTORIAL_LARGO.parent.mkdir(parents=True, exist_ok=True)
    with aforo.HISTORIAL_LARGO.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=aforo.COLUMNAS)
        w.writeheader()
        for r in sorted(filas, key=lambda r: (r["fecha"], r["potrero"])):
            w.writerow({c: r.get(c, "") for c in aforo.COLUMNAS})
    revisadas_ruta.write_text("\n".join(sorted(revisadas)), encoding="utf-8")
    faltan = len(pendientes) - hechas
    log(f"Listo: {hechas} escenas revisadas ahora, {len(nuevas)} observaciones nuevas, {len(filas)} en total. Faltan {faltan}.")
    (aforo.HISTORIAL_LARGO.with_name("historico_estado.json")).write_text(json.dumps(
        {"escenas_total": len(items), "faltan": faltan, "observaciones": len(filas),
         "actualizado": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes")}, indent=1), encoding="utf-8")


if __name__ == "__main__":
    with rasterio.Env(**aforo.OPCIONES_GDAL):
        main(sys.argv[1] if len(sys.argv) > 1 else "altos_de_medina", int(sys.argv[2]) if len(sys.argv) > 2 else 200)
