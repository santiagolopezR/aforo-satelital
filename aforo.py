"""
Aforo satelital de forraje verde por potrero.

Lee los potreros de la carpeta potreros/ (GeoJSON, KML o KMZ), busca imágenes
Sentinel-2 (óptico) y Sentinel-1 (radar) gratuitas en Microsoft Planetary
Computer, calcula índices por potrero y estima kg de forraje verde por m2 con
un modelo provisional por especie (config.json).

Salidas:
  data/historial.csv  -> todas las observaciones por fecha y potrero
  docs/datos.json     -> lo que lee el tablero web (docs/index.html)
"""
import csv
import datetime as dt
import json
import math
import sys
import unicodedata
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

try:  # en Windows usa los certificados del sistema (antivirus/redes que inspeccionan HTTPS)
    import truststore
    truststore.inject_into_ssl()
except ImportError:
    pass

import numpy as np
import planetary_computer
import pystac_client
import rasterio
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.errors import WindowError
from rasterio.features import geometry_mask
from rasterio.windows import Window, from_bounds
from shapely.geometry import MultiPolygon, Polygon, box, mapping, shape
from shapely.ops import transform as shp_transform, unary_union

RAIZ = Path(__file__).resolve().parent
CARPETA_POTREROS = RAIZ / "potreros"
HISTORIAL = RAIZ / "data" / "historial.csv"
SALIDA_WEB = RAIZ / "docs" / "datos.json"
STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"

COLUMNAS = ["fecha", "potrero", "especie", "sensor", "escena", "pct_valido",
            "ndvi", "evi", "ndre", "msavi", "vv_db", "vh_db", "kg_ms_ha", "kg_fv_m2"]

ALIAS_ESPECIE = {
    "humidicola": "humidicola", "brachiaria humidicola": "humidicola", "humidicula": "humidicola",
    "decumbens": "decumbens", "cumbens": "decumbens", "brachiaria decumbens": "decumbens",
    "mombasa": "mombasa", "mombaza": "mombasa", "guinea": "mombasa", "panicum": "mombasa",
    "kikuyo": "kikuyo", "kikuyu": "kikuyo", "quicuyo": "kikuyo",
}

# Clases SCL de Sentinel-2 que se aceptan: 4 = vegetación, 5 = suelo desnudo.
SCL_VALIDAS = (4, 5)


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def log(*a):
    print(*a, flush=True)


def sin_tildes(t):
    return "".join(c for c in unicodedata.normalize("NFD", t) if unicodedata.category(c) != "Mn")


def normalizar_especie(valor, config):
    if not valor:
        return config["especie_por_defecto"]
    v = sin_tildes(str(valor)).lower().strip()
    if v in ALIAS_ESPECIE:
        return ALIAS_ESPECIE[v]
    for alias, esp in ALIAS_ESPECIE.items():
        if alias in v:
            return esp
    log(f"  ! Especie '{valor}' no reconocida, se usa {config['especie_por_defecto']}")
    return config["especie_por_defecto"]


# ---------------------------------------------------------------- potreros

def _prop(props, *claves):
    bajas = {sin_tildes(str(k)).lower(): v for k, v in props.items()}
    for c in claves:
        if bajas.get(c) not in (None, ""):
            return str(bajas[c]).strip()
    return None


def _leer_geojson(ruta):
    datos = json.loads(ruta.read_text(encoding="utf-8-sig"))
    feats = datos["features"] if datos.get("type") == "FeatureCollection" else [datos]
    for i, f in enumerate(feats):
        props = f.get("properties") or {}
        yield {
            "nombre": _prop(props, "nombre", "name", "potrero", "lote") or f"{ruta.stem}_{i + 1}",
            "especie": _prop(props, "especie", "pasto", "forraje"),
            "geom": shape(f["geometry"]),
        }


def _sin_ns(tag):
    return tag.split("}", 1)[-1]


def _coords_kml(texto):
    pts = []
    for par in texto.split():
        partes = par.split(",")
        if len(partes) >= 2:
            pts.append((float(partes[0]), float(partes[1])))
    return pts


def _leer_kml_texto(texto, nombre_archivo):
    raiz = ET.fromstring(texto)
    n = 0
    for pm in raiz.iter():
        if _sin_ns(pm.tag) != "Placemark":
            continue
        n += 1
        nombre, props, poligonos = None, {}, []
        for el in pm.iter():
            t = _sin_ns(el.tag)
            if t == "name" and nombre is None:
                nombre = (el.text or "").strip()
            elif t == "Data":
                for v in el:
                    if _sin_ns(v.tag) == "value":
                        props[el.get("name", "")] = (v.text or "").strip()
            elif t == "SimpleData":
                props[el.get("name", "")] = (el.text or "").strip()
            elif t == "Polygon":
                ext, huecos = None, []
                for b in el.iter():
                    bt = _sin_ns(b.tag)
                    if bt in ("outerBoundaryIs", "innerBoundaryIs"):
                        for c in b.iter():
                            if _sin_ns(c.tag) == "coordinates":
                                pts = _coords_kml(c.text or "")
                                if bt == "outerBoundaryIs":
                                    ext = pts
                                else:
                                    huecos.append(pts)
                if ext and len(ext) >= 3:
                    poligonos.append(Polygon(ext, huecos))
        if not poligonos:
            continue
        geom = poligonos[0] if len(poligonos) == 1 else MultiPolygon(poligonos)
        yield {
            "nombre": nombre or f"{nombre_archivo}_{n}",
            "especie": _prop(props, "especie", "pasto", "forraje"),
            "geom": geom,
        }


def _leer_kml(ruta):
    yield from _leer_kml_texto(ruta.read_bytes(), ruta.stem)


def _leer_kmz(ruta):
    with zipfile.ZipFile(ruta) as z:
        for n in z.namelist():
            if n.lower().endswith(".kml"):
                yield from _leer_kml_texto(z.read(n), ruta.stem)


def cargar_potreros(config):
    lectores = {".geojson": _leer_geojson, ".json": _leer_geojson, ".kml": _leer_kml, ".kmz": _leer_kmz}
    # especies.csv opcional: columnas nombre,especie (sirve cuando el KML no trae la especie)
    mapa_especie = {}
    ruta_csv = CARPETA_POTREROS / "especies.csv"
    if ruta_csv.exists():
        with ruta_csv.open(encoding="utf-8-sig") as f:
            for fila in csv.DictReader(f):
                fila = {sin_tildes(k).lower().strip(): (v or "").strip() for k, v in fila.items()}
                if fila.get("nombre"):
                    mapa_especie[fila["nombre"].lower()] = fila.get("especie")

    potreros, vistos = [], set()
    for ruta in sorted(CARPETA_POTREROS.iterdir()):
        lector = lectores.get(ruta.suffix.lower())
        if not lector:
            continue
        try:
            for p in lector(ruta):
                if p["geom"].is_empty:
                    continue
                if not p["geom"].is_valid:
                    p["geom"] = p["geom"].buffer(0)
                nombre = p["nombre"]
                base, k = nombre, 2
                while nombre.lower() in vistos:
                    nombre, k = f"{base} ({k})", k + 1
                vistos.add(nombre.lower())
                p["nombre"] = nombre
                p["especie"] = normalizar_especie(p["especie"] or mapa_especie.get(base.lower()), config)
                p["archivo"] = ruta.name
                potreros.append(p)
        except Exception as e:  # un archivo dañado no debe frenar a los demás
            log(f"  ! No se pudo leer {ruta.name}: {e}")
    return potreros


def zona_utm(lon, lat):
    zona = int((lon + 180) // 6) + 1
    return (32600 if lat >= 0 else 32700) + zona


def area_ha(geom):
    c = geom.centroid
    t = Transformer.from_crs(4326, zona_utm(c.x, c.y), always_xy=True).transform
    return shp_transform(t, geom).area / 10000


# ---------------------------------------------------------------- historial

def leer_historial():
    if not HISTORIAL.exists():
        return []
    with HISTORIAL.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def guardar_historial(filas):
    HISTORIAL.parent.mkdir(parents=True, exist_ok=True)
    filas.sort(key=lambda r: (r["fecha"], r["potrero"], r["sensor"]))
    with HISTORIAL.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNAS)
        w.writeheader()
        for r in filas:
            w.writerow({c: r.get(c, "") for c in COLUMNAS})


# ---------------------------------------------------------------- modelo

def estimar(ndvi, especie, config):
    """Modelo provisional sin calibrar. Devuelve (kg MS/ha, kg FV/m2)."""
    p = config["especies"][especie]
    if ndvi is None or ndvi < config["ndvi_minimo"]:
        return 0.0, 0.0
    ms = min(p["a"] * math.exp(p["b"] * ndvi), p["max_kg_ms_ha"])
    fv_m2 = ms / p["pct_ms"] / 10000
    return round(ms), round(fv_m2, 3)


# ---------------------------------------------------------------- satélite

def _num(x, dec=4):
    return "" if x is None or not np.isfinite(x) else round(float(x), dec)


def leer_banda(href, bounds, forma=None):
    """Lee solo el pedazo de la imagen que cubre 'bounds' (en el CRS de la imagen)."""
    with rasterio.open(href) as src:
        win = from_bounds(*bounds, transform=src.transform)
        if forma is None:
            win = win.round_offsets().round_lengths().intersection(Window(0, 0, src.width, src.height))
            datos = src.read(1, window=win)
            return datos, src.window_transform(win), src.window_bounds(win)
        # Bandas de otra resolución: se lee exactamente el mismo recuadro y se lleva a la malla de 10 m.
        datos = src.read(1, window=win, out_shape=forma, resampling=Resampling.nearest)
        return datos, None, None


def preparar_geoms(potreros, epsg, recorte_m):
    """Pasa cada potrero al CRS de la imagen y le quita el borde."""
    t = Transformer.from_crs(4326, epsg, always_xy=True).transform
    salida = {}
    for p in potreros:
        g = shp_transform(t, p["geom"]).buffer(-recorte_m)
        if g.is_empty:  # potrero muy pequeño: se usa sin recorte
            g = shp_transform(t, p["geom"])
        salida[p["nombre"]] = g
    return salida


def epsg_de(item, asset):
    """EPSG de la escena: proj:epsg (STAC viejo), proj:code (STAC nuevo) o leyendo la imagen."""
    for fuente in (item.properties, item.assets[asset].extra_fields):
        if fuente.get("proj:epsg"):
            return int(fuente["proj:epsg"])
        codigo = fuente.get("proj:code") or ""
        if codigo.upper().startswith("EPSG:"):
            return int(codigo.split(":")[1])
    with rasterio.open(item.assets[asset].href) as src:
        return src.crs.to_epsg()


def procesar_s2(item, potreros, config):
    epsg = epsg_de(item, "B04")
    geoms = preparar_geoms(potreros, epsg, config["recorte_borde_m"])
    bounds = unary_union(list(geoms.values())).bounds
    a = item.assets

    # Primero la máscara de nubes (es liviana): si ningún potrero está despejado no se baja nada más.
    try:
        scl20, transf20, caja20 = leer_banda(a["SCL"].href, bounds)
    except WindowError:
        return []  # la escena no cubre los potreros
    despejado20 = np.isin(scl20, SCL_VALIDAS)
    alguno = False
    for p in potreros:
        g = geoms[p["nombre"]]
        if not g.intersects(box(*caja20)):
            continue
        dentro = geometry_mask([mapping(g)], out_shape=scl20.shape, transform=transf20, invert=True, all_touched=True)
        if dentro.sum() and (dentro & despejado20).sum() / dentro.sum() >= config["min_pct_valido"]:
            alguno = True
            break
    if not alguno:
        return []

    try:
        rojo, transf, bounds_reales = leer_banda(a["B04"].href, bounds)
    except WindowError:
        return []
    forma = rojo.shape
    nir = leer_banda(a["B08"].href, bounds_reales)[0]
    azul = leer_banda(a["B02"].href, bounds_reales)[0]
    borde_rojo = leer_banda(a["B05"].href, bounds_reales, forma)[0]
    scl = leer_banda(a["SCL"].href, bounds_reales, forma)[0]
    if nir.shape != forma or azul.shape != forma:
        return []

    # Desde la línea base 04.00 (2022) los valores traen un desfase de +1000.
    linea_base = float(item.properties.get("s2:processing_baseline", "0") or 0)
    desfase = 1000 if linea_base >= 4.0 else 0

    def refl(x):
        return (x.astype("float32") - desfase) / 10000

    r, n, b, re = refl(rojo), refl(nir), refl(azul), refl(borde_rojo)
    validos_img = np.isin(scl, SCL_VALIDAS) & (rojo > 0) & (nir > 0)

    filas = []
    for p in potreros:
        g = geoms[p["nombre"]]
        if not g.intersects(box(*bounds_reales)):
            continue
        dentro = geometry_mask([mapping(g)], out_shape=forma, transform=transf, invert=True, all_touched=False)
        total = int(dentro.sum())
        if total < config["min_pixeles"]:
            continue
        ok = dentro & validos_img
        pct = ok.sum() / total
        if pct < config["min_pct_valido"]:
            continue
        rr, nn, bb, ee = r[ok], n[ok], b[ok], re[ok]
        with np.errstate(divide="ignore", invalid="ignore"):
            ndvi = np.nanmean((nn - rr) / (nn + rr))
            evi = np.nanmean(2.5 * (nn - rr) / (nn + 6 * rr - 7.5 * bb + 1))
            ndre = np.nanmean((nn - ee) / (nn + ee))
            msavi = np.nanmean((2 * nn + 1 - np.sqrt((2 * nn + 1) ** 2 - 8 * (nn - rr))) / 2)
        filas.append({
            "fecha": item.datetime.date().isoformat(), "potrero": p["nombre"], "especie": p["especie"],
            "sensor": "S2", "escena": item.id, "pct_valido": _num(pct, 3),
            "ndvi": _num(ndvi), "evi": _num(evi), "ndre": _num(ndre), "msavi": _num(msavi),
        })
    return filas


def procesar_s1(item, potreros, config):
    epsg = epsg_de(item, "vv")
    geoms = preparar_geoms(potreros, epsg, config["recorte_borde_m"])
    bounds = unary_union(list(geoms.values())).bounds
    try:
        vv, transf, bounds_reales = leer_banda(item.assets["vv"].href, bounds)
        vh = leer_banda(item.assets["vh"].href, bounds_reales, vv.shape)[0]
    except WindowError:
        return []
    forma = vv.shape
    validos_img = np.isfinite(vv) & np.isfinite(vh) & (vv > 0) & (vh > 0)

    filas = []
    for p in potreros:
        g = geoms[p["nombre"]]
        if not g.intersects(box(*bounds_reales)):
            continue
        dentro = geometry_mask([mapping(g)], out_shape=forma, transform=transf, invert=True)
        total = int(dentro.sum())
        if total < config["min_pixeles"]:
            continue
        ok = dentro & validos_img
        pct = ok.sum() / total
        if pct < config["min_pct_valido"]:
            continue
        # Se promedia en lineal y luego se pasa a decibeles.
        filas.append({
            "fecha": item.datetime.date().isoformat(), "potrero": p["nombre"], "especie": p["especie"],
            "sensor": "S1", "escena": item.id, "pct_valido": _num(pct, 3),
            "vv_db": _num(10 * np.log10(vv[ok].mean()), 2), "vh_db": _num(10 * np.log10(vh[ok].mean()), 2),
        })
    return filas


def buscar(catalogo, coleccion, bbox, desde, hasta, filtro=None):
    busqueda = catalogo.search(collections=[coleccion], bbox=bbox,
                               datetime=f"{desde.isoformat()}/{hasta.isoformat()}", query=filtro)
    return list(busqueda.items())


# ---------------------------------------------------------------- tablero

def exportar_web(potreros, filas, config):
    por_potrero = {}
    for r in filas:
        por_potrero.setdefault(r["potrero"], []).append(r)

    feats = []
    for p in potreros:
        serie = sorted(por_potrero.get(p["nombre"], []), key=lambda r: r["fecha"])
        s2 = [r for r in serie if r["sensor"] == "S2"]
        s1 = [r for r in serie if r["sensor"] == "S1"]
        ultimo = s2[-1] if s2 else None
        ha = area_ha(p["geom"])
        props = {
            "nombre": p["nombre"], "especie": p["especie"],
            "especie_nombre": config["especies"][p["especie"]]["nombre"],
            "area_ha": round(ha, 2),
            "ultima_fecha": ultimo["fecha"] if ultimo else None,
            "kg_fv_m2": float(ultimo["kg_fv_m2"]) if ultimo else None,
            "kg_ms_ha": float(ultimo["kg_ms_ha"]) if ultimo else None,
            "ndvi": float(ultimo["ndvi"]) if ultimo else None,
            "t_fv_potrero": round(float(ultimo["kg_fv_m2"]) * ha * 10, 1) if ultimo else None,
            "serie_s2": [{"f": r["fecha"], "fv": float(r["kg_fv_m2"]), "ms": float(r["kg_ms_ha"]),
                          "ndvi": float(r["ndvi"]), "v": float(r["pct_valido"])} for r in s2],
            "serie_s1": [{"f": r["fecha"], "vv": float(r["vv_db"]), "vh": float(r["vh_db"])}
                         for r in s1 if r["vv_db"] != "" and r["vh_db"] != ""],
        }
        feats.append({"type": "Feature", "properties": props, "geometry": mapping(p["geom"])})

    SALIDA_WEB.parent.mkdir(parents=True, exist_ok=True)
    SALIDA_WEB.write_text(json.dumps({
        "actualizado": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"),
        "modelo": "provisional sin calibrar",
        "type": "FeatureCollection", "features": feats,
    }, ensure_ascii=False), encoding="utf-8")


# ---------------------------------------------------------------- principal

def main():
    config = json.loads((RAIZ / "config.json").read_text(encoding="utf-8"))
    potreros = cargar_potreros(config)
    if not potreros:
        log("No hay potreros en la carpeta potreros/. Sube un .geojson, .kml o .kmz.")
        exportar_web([], [], config)
        return
    log(f"{len(potreros)} potreros: " + ", ".join(f"{p['nombre']} ({p['especie']})" for p in potreros))

    nombres = {p["nombre"] for p in potreros}
    # Se descartan filas de potreros que ya no están en la carpeta.
    historial = [r for r in leer_historial() if r["potrero"] in nombres]
    hechos = {(r["escena"], r["potrero"]) for r in historial}

    hoy = dt.date.today()
    inicio_nuevo = hoy - dt.timedelta(days=config["dias_historial_inicial"])
    ultimas = {}
    for r in historial:
        ultimas[r["potrero"]] = max(ultimas.get(r["potrero"], r["fecha"]), r["fecha"])
    # Desde la fecha más vieja que haga falta (potreros nuevos piden todo el historial inicial).
    desde = min(
        (dt.date.fromisoformat(ultimas[n]) - dt.timedelta(days=5)) if n in ultimas else inicio_nuevo
        for n in nombres
    )
    desde = max(desde, inicio_nuevo)
    bbox = list(unary_union([p["geom"] for p in potreros]).bounds)

    catalogo = pystac_client.Client.open(STAC_URL, modifier=planetary_computer.sign_inplace)
    nuevas = []

    s2 = buscar(catalogo, "sentinel-2-l2a", bbox, desde, hoy,
                {"eo:cloud_cover": {"lt": config["max_nubes_escena"]}})
    log(f"Sentinel-2: {len(s2)} escenas desde {desde}")
    for item in sorted(s2, key=lambda i: i.datetime):
        pendientes = [p for p in potreros if (item.id, p["nombre"]) not in hechos]
        if not pendientes:
            continue
        try:
            filas = procesar_s2(item, pendientes, config)
        except Exception as e:
            log(f"  ! {item.id}: {e}")
            continue
        log(f"  {item.datetime.date()} {item.id[:38]}  nubes {item.properties.get('eo:cloud_cover', 0):.0f}%  -> {len(filas)} potreros útiles")
        nuevas += filas
        # Escena revisada aunque estuviera nublada: se marca para no volver a bajarla.
        hechos.update((item.id, p["nombre"]) for p in pendientes)

    s1 = buscar(catalogo, "sentinel-1-rtc", bbox, desde, hoy)
    log(f"Sentinel-1: {len(s1)} escenas desde {desde}")
    for item in sorted(s1, key=lambda i: i.datetime):
        pendientes = [p for p in potreros if (item.id, p["nombre"]) not in hechos]
        if not pendientes:
            continue
        try:
            filas = procesar_s1(item, pendientes, config)
        except Exception as e:
            log(f"  ! {item.id}: {e}")
            continue
        log(f"  {item.datetime.date()} {item.id[:38]}  -> {len(filas)} potreros")
        nuevas += filas

    filas = historial + nuevas

    # Si un potrero cae en dos escenas el mismo día, se deja la de más píxeles válidos.
    mejor = {}
    for r in filas:
        k = (r["fecha"], r["potrero"], r["sensor"])
        if k not in mejor or float(r["pct_valido"] or 0) > float(mejor[k]["pct_valido"] or 0):
            mejor[k] = r
    filas = list(mejor.values())

    # El modelo se recalcula siempre, así un cambio en config.json se aplica a todo el historial.
    especie_de = {p["nombre"]: p["especie"] for p in potreros}
    for r in filas:
        r["especie"] = especie_de[r["potrero"]]
        if r["sensor"] == "S2" and r.get("ndvi") not in ("", None):
            r["kg_ms_ha"], r["kg_fv_m2"] = estimar(float(r["ndvi"]), r["especie"], config)

    guardar_historial(filas)
    exportar_web(potreros, filas, config)
    log(f"Listo: {len(nuevas)} observaciones nuevas, {len(filas)} en total.")


# Ajustes para leer imágenes por internet sin descargarlas completas.
OPCIONES_GDAL = dict(
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
    CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.TIF,.tiff",
    GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
    GDAL_HTTP_MULTIPLEX="YES",
    GDAL_HTTP_MAX_RETRY="4",
    GDAL_HTTP_RETRY_DELAY="3",
)

if __name__ == "__main__":
    with rasterio.Env(**OPCIONES_GDAL):
        sys.exit(main())
