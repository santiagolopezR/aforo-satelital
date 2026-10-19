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
from rasterio.errors import NotGeoreferencedWarning, WindowError
from rasterio.features import geometry_mask, rasterize
from rasterio.merge import merge as unir_rasters
from rasterio.warp import Resampling as ResamplingWarp, calculate_default_transform, reproject
import hashlib
from rasterio.windows import Window, from_bounds
from shapely.geometry import MultiPolygon, Polygon, box, mapping, shape
from shapely.ops import transform as shp_transform, unary_union

RAIZ = Path(__file__).resolve().parent
CARPETA_POTREROS = RAIZ / "potreros"
HISTORIAL = RAIZ / "data" / "historial.csv"
SALIDA_WEB = RAIZ / "docs" / "datos.json"
STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"

COLUMNAS = ["fecha", "potrero", "especie", "sensor", "escena", "pct_valido",
            "ndvi", "evi", "ndre", "msavi", "vv_db", "vh_db", "kg_ms_ha", "kg_fv_m2", "pct_pobre", "hist",
            "pct_arboles", "mascara"]
CARPETA_PIXELES = RAIZ / "docs" / "pixeles"
CARPETA_EXCLUSIONES = RAIZ / "exclusiones"
MASCARA_ARBOLES = CARPETA_EXCLUSIONES / "arboles.tif"         # la que se usa (Meta + ortomosaicos de dron)
MASCARA_META = CARPETA_EXCLUSIONES / "arboles_meta.tif"       # solo alturas de Meta, base y respaldo
KMZ_ARBOLES = RAIZ / "docs" / "arboles.kmz"
# Alturas de árboles de Meta / World Resources Institute (1 m, gratis), en cuadros "quadkey" de nivel 9.
URL_ALTURAS = "https://dataforgood-fb-data.s3.amazonaws.com/forests/v1/alsgedi_global_v6_float/chm/{}.tif"
COD_ARBOL = 255  # código del píxel excluido en la imagen del mapa de detalle

# Histograma del NDVI de los píxeles de cada potrero: clases de 0,01 entre -0,2 y 1,0.
HIST_MIN, HIST_PASO, HIST_N = -0.2, 0.01, 120

ALIAS_ESPECIE = {
    "humidicola": "humidicola", "brachiaria humidicola": "humidicola", "humidicula": "humidicola",
    "decumbens": "decumbens", "cumbens": "decumbens", "brachiaria decumbens": "decumbens",
    "mombasa": "mombasa", "mombaza": "mombasa", "guinea": "mombasa", "panicum": "mombasa",
    "kikuyo": "kikuyo", "kikuyu": "kikuyo", "quicuyo": "kikuyo",
}

# Clases SCL de Sentinel-2 que se aceptan: 4 = vegetación, 5 = suelo desnudo.
SCL_VALIDAS = (4, 5)


import warnings

# Los PNG del mapa de detalle no llevan georreferencia propia (va en pixeles/indice.json).
warnings.filterwarnings("ignore", category=NotGeoreferencedWarning)

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


def quadkey(lon, lat, nivel=9):
    n = 2 ** nivel
    x = int((lon + 180) / 360 * n)
    y = int((1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n)
    q = ""
    for i in range(nivel, 0, -1):
        m = 1 << (i - 1)
        q += str((1 if x & m else 0) + (2 if y & m else 0))
    return q


def generar_mascara_arboles(potreros, config):
    """
    Una sola vez: baja las alturas de árboles (1 m) sobre la finca y guarda 1 donde hay árbol
    (altura >= altura_arbol_m). Para regenerarla, se borra exclusiones/arboles_meta.tif.
    La máscara final (arboles.tif) arranca igual y arboles_dron.py la mejora con los ortomosaicos.
    """
    oeste, sur, este, norte = unary_union([p["geom"] for p in potreros]).buffer(0.001).bounds
    claves = sorted({quadkey(x, y) for x in (oeste, este) for y in (sur, norte)})
    log(f"Generando máscara de árboles (cuadros {', '.join(claves)})...")
    a3857 = Transformer.from_crs(4326, 3857, always_xy=True).transform
    x0, y0 = a3857(oeste, sur)
    x1, y1 = a3857(este, norte)
    fuentes = [rasterio.open("/vsicurl/" + URL_ALTURAS.format(c)) for c in claves]
    try:
        alturas, transf = unir_rasters(fuentes, bounds=(x0, y0, x1, y1), nodata=0)
    finally:
        for f in fuentes:
            f.close()
    arbol = (alturas[0] >= config["altura_arbol_m"]).astype("uint8")
    CARPETA_EXCLUSIONES.mkdir(parents=True, exist_ok=True)
    with rasterio.open(MASCARA_META, "w", driver="GTiff", width=arbol.shape[1], height=arbol.shape[0],
                       count=1, dtype="uint8", crs="EPSG:3857", transform=transf, compress="deflate",
                       nbits=1) as dst:
        dst.write(arbol, 1)
    log(f"  Máscara lista: {arbol.mean() * 100:.1f} % del área con árboles de {config['altura_arbol_m']} m o más")


def exportar_kmz_arboles():
    """Árboles detectados + bosques en un KMZ para revisarlos sobre Google Earth (docs/arboles.kmz)."""
    import zipfile
    from rasterio.features import shapes
    from xml.sax.saxutils import escape
    if not MASCARA_ARBOLES.exists():
        return
    with rasterio.open(MASCARA_ARBOLES) as src:
        m = src.read(1)
        a4326 = Transformer.from_crs(src.crs, 4326, always_xy=True).transform
        poligonos = [shape(g).simplify(0.6) for g, v in shapes(m, mask=m > 0, transform=src.transform) if v]

    def kml_poligono(g):
        partes = []
        for pol in (g.geoms if g.geom_type == "MultiPolygon" else [g]):
            def anillo(c):
                return " ".join(f"{x:.7f},{y:.7f},0" for x, y in c)
            huecos = "".join(f"<innerBoundaryIs><LinearRing><coordinates>{anillo(i.coords)}</coordinates></LinearRing></innerBoundaryIs>"
                             for i in pol.interiors)
            partes.append(f"<Polygon><outerBoundaryIs><LinearRing><coordinates>{anillo(pol.exterior.coords)}"
                          f"</coordinates></LinearRing></outerBoundaryIs>{huecos}</Polygon>")
        return "".join(partes)

    marcas = [f"<Placemark><styleUrl>#arbol</styleUrl>{kml_poligono(shp_transform(a4326, g))}</Placemark>"
              for g in poligonos if not g.is_empty]
    bosques = [f"<Placemark><name>{escape(str(z.get('nombre', '')))}</name><styleUrl>#bosque</styleUrl>{kml_poligono(z['geom'])}</Placemark>"
               for z in _leer_geojson(CARPETA_EXCLUSIONES / "bosques.geojson")] if (CARPETA_EXCLUSIONES / "bosques.geojson").exists() else []
    kml = ('<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document>'
           "<name>Árboles detectados (aforo satelital)</name>"
           '<Style id="arbol"><LineStyle><color>ff76e600</color><width>1.5</width></LineStyle><PolyStyle><color>5576e600</color></PolyStyle></Style>'
           '<Style id="bosque"><LineStyle><color>ff00ffff</color><width>2</width></LineStyle><PolyStyle><fill>0</fill></PolyStyle></Style>'
           f"<Folder><name>Árboles detectados ({len(marcas)})</name>{''.join(marcas)}</Folder>"
           f"<Folder><name>Bosques (capa QGIS)</name>{''.join(bosques)}</Folder>"
           "</Document></kml>")
    KMZ_ARBOLES.parent.mkdir(parents=True, exist_ok=True)
    nuevo = KMZ_ARBOLES.with_suffix(".tmp")
    with zipfile.ZipFile(nuevo, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(zipfile.ZipInfo("doc.kml", date_time=(2026, 1, 1, 0, 0, 0)), kml)  # fecha fija: sin cambios falsos en git
    if KMZ_ARBOLES.exists() and KMZ_ARBOLES.read_bytes() == nuevo.read_bytes():
        nuevo.unlink()
    else:
        nuevo.replace(KMZ_ARBOLES)


def cargar_exclusiones():
    """Zonas que no se cuentan como pasto: polígonos en exclusiones/ (bosques, árboles marcados a mano...)."""
    lectores = {".geojson": _leer_geojson, ".json": _leer_geojson, ".kml": _leer_kml, ".kmz": _leer_kmz}
    geoms = []
    if CARPETA_EXCLUSIONES.exists():
        for ruta in sorted(CARPETA_EXCLUSIONES.iterdir()):
            lector = lectores.get(ruta.suffix.lower())
            if lector:
                geoms += [z["geom"] if z["geom"].is_valid else z["geom"].buffer(0) for z in lector(ruta)]
    return geoms


def huella_mascara(config):
    """Identifica la máscara vigente; si cambia, las imágenes se vuelven a procesar."""
    h = hashlib.sha1(f"{config['altura_arbol_m']}|{config['fraccion_arbol_excluir']}|{config['recorte_borde_m']}".encode())
    if MASCARA_ARBOLES.exists():
        h.update(MASCARA_ARBOLES.read_bytes())
    if CARPETA_EXCLUSIONES.exists():
        for ruta in sorted(CARPETA_EXCLUSIONES.iterdir()):
            if ruta.suffix.lower() in (".geojson", ".json", ".kml", ".kmz"):
                h.update(ruta.name.encode() + ruta.read_bytes())
    # (arboles_meta.tif no entra: solo cuenta la máscara final arboles.tif)
    return h.hexdigest()[:10]


def mascara_excluida(forma, transf, epsg, exclusiones, config):
    """True en los píxeles de 10 m que son árbol o bosque (no se cuentan como pasto)."""
    fuera = np.zeros(forma, dtype=bool)
    if MASCARA_ARBOLES.exists():
        fraccion = np.zeros(forma, dtype="float32")
        with rasterio.open(MASCARA_ARBOLES) as src:
            reproject(src.read(1).astype("float32"), fraccion, src_transform=src.transform, src_crs=src.crs,
                      dst_transform=transf, dst_crs=epsg, resampling=ResamplingWarp.average)
        fuera |= fraccion >= config["fraccion_arbol_excluir"]
    if exclusiones:
        t = Transformer.from_crs(4326, epsg, always_xy=True).transform
        fuera |= geometry_mask([mapping(shp_transform(t, g)) for g in exclusiones],
                               out_shape=forma, transform=transf, invert=True)
    return fuera


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

def forma_curva(ndvi, p, config):
    """Parte de la curva que depende del NDVI; 0 en suelo desnudo."""
    return max(0.0, math.exp(p["b"] * (ndvi - config["ndvi_suelo"])) - 1)


def histograma(ndvi_pixeles):
    """Cuenta los píxeles por clase de NDVI; se guarda como texto 'clase:cantidad ...'."""
    idx = np.clip(((ndvi_pixeles - HIST_MIN) / HIST_PASO).astype(int), 0, HIST_N - 1)
    cuentas = np.bincount(idx, minlength=HIST_N)
    return " ".join(f"{i}:{c}" for i, c in enumerate(cuentas) if c)


def estimar_por_pixel(hist, especie, config):
    """
    Convierte cada píxel a kg con la curva de la especie y promedia (más exacto que promediar el NDVI).
    Devuelve (kg MS/ha, kg FV/m2, fracción del potrero con cobertura pobre).
    """
    total = ms_sum = fv_sum = pobre = 0
    umbral = config["umbral_cobertura_pobre_kg"]
    for par in hist.split():
        i, c = (int(v) for v in par.split(":"))
        ms, fv = estimar(HIST_MIN + (i + 0.5) * HIST_PASO, especie, config)
        total += c
        ms_sum += ms * c
        fv_sum += fv * c
        if fv < umbral:
            pobre += c
    if not total:
        return 0.0, 0.0, 0.0
    return round(ms_sum / total), round(fv_sum / total, 3), round(pobre / total, 3)


def estimar(ndvi, especie, config):
    """Devuelve (kg MS/ha, kg FV/m2) según la curva de la especie."""
    p = config["especies"][especie]
    if ndvi is None:
        return 0.0, 0.0
    ms = min(p["k"] * forma_curva(ndvi, p, config), p["max_kg_ms_ha"])
    fv_m2 = ms / p["pct_ms"] / 10000
    return round(ms), round(fv_m2, 3)


# ---------------------------------------------------------------- aforos de campo

COLUMNAS_AFOROS = ["fecha", "potrero", "kg_fv_m2", "pct_ms", "metodo", "observaciones", "origen"]
AFOROS = RAIZ / "aforos_campo" / "aforos.csv"
CALIBRACION = RAIZ / "data" / "calibracion.csv"
ISSUES_A_CERRAR = RAIZ / "data" / ".issues_a_cerrar.json"


def _clave(t):
    return sin_tildes(str(t)).upper().strip()


def leer_aforos(potreros):
    if not AFOROS.exists():
        return []
    por_clave = {_clave(p["nombre"]): p["nombre"] for p in potreros}
    salida = []
    with AFOROS.open(encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            try:
                r["kg_fv_m2"] = float(str(r["kg_fv_m2"]).replace(",", "."))
                dt.date.fromisoformat(r["fecha"])
            except (TypeError, ValueError):
                log(f"  ! Aforo con datos inválidos, se ignora: {r}")
                continue
            nombre = por_clave.get(_clave(r["potrero"]))
            if not nombre:
                continue  # potrero que ya no existe
            r["potrero"] = nombre
            pct = str(r.get("pct_ms") or "").replace(",", ".").strip()
            r["pct_ms"] = float(pct) / (100 if float(pct) > 1 else 1) if pct else None
            salida.append(r)
    return salida


def agregar_aforos(nuevos):
    existe = AFOROS.exists() and AFOROS.stat().st_size > 0
    AFOROS.parent.mkdir(parents=True, exist_ok=True)
    with AFOROS.open("a", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNAS_AFOROS)
        if not existe:
            w.writeheader()
        for r in nuevos:
            w.writerow({c: r.get(c, "") for c in COLUMNAS_AFOROS})


def _secciones_issue(cuerpo):
    """Convierte el cuerpo de un formulario de GitHub ('### Título\\n\\nvalor') en {título: valor}."""
    out = {}
    for bloque in ("\n" + (cuerpo or "")).split("\n### ")[1:]:
        titulo, _, valor = bloque.partition("\n")
        valor = valor.strip()
        out[_clave(titulo)] = "" if valor == "_No response_" else valor
    return out


def importar_issues(potreros):
    """Lee los formularios 'Registrar aforo' abiertos en GitHub y los pasa a aforos.csv."""
    import os
    import requests
    token, repo = os.environ.get("GITHUB_TOKEN"), os.environ.get("GITHUB_REPOSITORY")
    if not token or not repo:
        return
    cab = {"Authorization": f"Bearer {token}", "Accept": "application/vnd.github+json"}
    r = requests.get(f"https://api.github.com/repos/{repo}/issues", headers=cab,
                     params={"state": "open", "per_page": 100}, timeout=30)
    r.raise_for_status()
    por_clave = {_clave(p["nombre"]): p["nombre"] for p in potreros}
    nuevos, respuestas = [], []
    for issue in r.json():
        if "pull_request" in issue:
            continue
        s = _secciones_issue(issue.get("body"))
        if "FORRAJE VERDE (KG/M²)" not in s:
            continue  # no es un formulario de aforo
        errores = []
        # El repositorio es público: solo el dueño y los colaboradores pueden cargar aforos.
        if issue.get("author_association") not in ("OWNER", "MEMBER", "COLLABORATOR"):
            respuestas.append({"numero": issue["number"], "texto": "Solo el dueño del repositorio y sus "
                               "colaboradores pueden registrar aforos. Este no se tuvo en cuenta."})
            continue
        nombre = por_clave.get(_clave(s.get("POTRERO", "")))
        if not nombre:
            errores.append(f"no encontré el potrero «{s.get('POTRERO', '')}» (debe escribirse como sale en el mapa)")
        fecha = s.get("FECHA DEL AFORO") or issue["created_at"][:10]
        try:
            dt.date.fromisoformat(fecha)
        except ValueError:
            errores.append(f"la fecha «{fecha}» debe ser AAAA-MM-DD")
        try:
            kg = float(s.get("FORRAJE VERDE (KG/M²)", "").replace(",", "."))
            if not 0 <= kg <= 15:
                errores.append("los kg/m² deben estar entre 0 y 15")
        except ValueError:
            errores.append("los kg/m² deben ser un número, por ejemplo 1.8")
        pct = s.get("MATERIA SECA (%)", "").replace("%", "").replace(",", ".").strip()
        if pct:
            try:
                float(pct)
            except ValueError:
                errores.append("la materia seca debe ser un número, por ejemplo 25")
        if errores:
            respuestas.append({"numero": issue["number"], "texto": "No se pudo registrar: " + "; ".join(errores)
                               + ". Crea un aforo nuevo con los datos corregidos."})
            continue
        nuevos.append({"fecha": fecha, "potrero": nombre, "kg_fv_m2": kg, "pct_ms": pct,
                       "metodo": s.get("MÉTODO", s.get("METODO", "")),
                       "observaciones": " ".join(s.get("OBSERVACIONES", "").split()),
                       "origen": f"issue #{issue['number']}"})
        respuestas.append({"numero": issue["number"], "potrero": nombre, "fecha": fecha, "kg": kg})
    if nuevos:
        agregar_aforos(nuevos)
        log(f"Aforos nuevos desde GitHub: {len(nuevos)}")
    if respuestas:
        ISSUES_A_CERRAR.parent.mkdir(parents=True, exist_ok=True)
        ISSUES_A_CERRAR.write_text(json.dumps(respuestas, ensure_ascii=False), encoding="utf-8")


def emparejar(aforos, filas, config):
    """A cada aforo le busca la imagen Sentinel-2 más cercana del mismo potrero."""
    s2 = {}
    for r in filas:
        if r["sensor"] == "S2" and r.get("hist"):
            s2.setdefault(r["potrero"], []).append(r)
    pares = []
    for a in aforos:
        fa = dt.date.fromisoformat(a["fecha"])
        cands = [(abs((dt.date.fromisoformat(r["fecha"]) - fa).days), r) for r in s2.get(a["potrero"], [])]
        cands = [c for c in cands if c[0] <= config["dias_par_aforo_imagen"]]
        if cands:
            d, r = min(cands, key=lambda c: c[0])
            pares.append({"aforo": a, "obs": r, "dias": d})
    return pares


def calibrar(pares, especie_de, config):
    """Con suficientes aforos por especie recalcula k (y pct_ms si hay datos). Devuelve un resumen."""
    resumen = {}
    for esp, p in config["especies"].items():
        mios = [x for x in pares if especie_de[x["aforo"]["potrero"]] == esp]
        pcts = [x["aforo"]["pct_ms"] for x in mios if x["aforo"]["pct_ms"]]
        info = {"aforos": len(mios), "calibrada": False}
        if len(mios) >= config["min_aforos_calibrar"]:
            if len(pcts) >= 5:
                p["pct_ms"] = round(sum(pcts) / len(pcts), 3)
            ys = [x["aforo"]["kg_fv_m2"] * 10000 * (x["aforo"]["pct_ms"] or p["pct_ms"]) for x in mios]
            fs = [forma_curva(float(x["obs"]["ndvi"]), p, config) for x in mios]
            den = sum(f * f for f in fs)
            if den > 0:
                p["k"] = round(sum(y * f for y, f in zip(ys, fs)) / den, 1)
                info["calibrada"] = True
        info.update(k=p["k"], pct_ms=p["pct_ms"])
        resumen[esp] = info
    return resumen


def guardar_calibracion(pares, config, especie_de):
    CALIBRACION.parent.mkdir(parents=True, exist_ok=True)
    with CALIBRACION.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["fecha_aforo", "potrero", "especie", "kg_fv_m2_campo", "fecha_imagen", "ndvi",
                    "kg_fv_m2_satelite", "error_pct"])
        for x in sorted(pares, key=lambda x: (x["aforo"]["fecha"], x["aforo"]["potrero"])):
            a, o = x["aforo"], x["obs"]
            esp = especie_de[a["potrero"]]
            sat = estimar(float(o["ndvi"]), esp, config)[1]
            err = round((sat - a["kg_fv_m2"]) / a["kg_fv_m2"] * 100) if a["kg_fv_m2"] else ""
            w.writerow([a["fecha"], a["potrero"], esp, a["kg_fv_m2"], o["fecha"], o["ndvi"], sat, err])


def ajustar_con_aforos(serie, aforos_potrero, sat_en_aforo, config):
    """
    Corrige la serie satelital de un potrero con su último aforo de campo:
    el día del aforo vale lo medido y la diferencia se va desvaneciendo en 'dias_efecto_aforo'.
    """
    dias = config["dias_efecto_aforo"]
    for s in serie:
        fs = dt.date.fromisoformat(s["f"])
        previos = [a for a in aforos_potrero if dt.date.fromisoformat(a["fecha"]) <= fs]
        s["aj"] = None
        if not previos:
            continue
        a = max(previos, key=lambda a: a["fecha"])
        d = (fs - dt.date.fromisoformat(a["fecha"])).days
        base = sat_en_aforo.get(a["fecha"])
        if d >= dias or base is None:
            continue
        s["aj"] = round(max(0.0, s["fv"] + (a["kg_fv_m2"] - base) * (1 - d / dias)), 3)
    return serie


def cerrar_issues():
    """Responde y cierra en GitHub los formularios de aforo ya procesados (se llama tras guardar)."""
    import os
    import requests
    if not ISSUES_A_CERRAR.exists():
        return
    token, repo = os.environ.get("GITHUB_TOKEN"), os.environ.get("GITHUB_REPOSITORY")
    datos = json.loads(ISSUES_A_CERRAR.read_text(encoding="utf-8"))
    cab = {"Authorization": f"Bearer {token}", "Accept": "application/vnd.github+json"}
    web = json.loads(SALIDA_WEB.read_text(encoding="utf-8")) if SALIDA_WEB.exists() else {"features": []}
    sat = {f["properties"]["nombre"]: f["properties"] for f in web["features"]}
    for x in datos:
        texto = x.get("texto")
        if not texto:
            p = sat.get(x["potrero"], {})
            texto = (f"Registrado: **{x['potrero']}**, {x['fecha']}, **{x['kg']} kg FV/m²**.\n\n"
                     f"El satélite estimaba {p.get('kg_fv_m2_satelite', '—')} kg FV/m² "
                     f"(imagen del {p.get('ultima_fecha_satelite') or '—'}). El mapa ya usa el dato de campo.")
        u = f"https://api.github.com/repos/{repo}/issues/{x['numero']}"
        requests.post(u + "/comments", headers=cab, json={"body": texto}, timeout=30)
        requests.patch(u, headers=cab, json={"state": "closed"}, timeout=30)
    ISSUES_A_CERRAR.unlink()


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


def guardar_pixeles(fecha, ndvi, valido, excluido, potreros, geoms, transf, epsg, config, huella=""):
    """
    Guarda una imagen PNG de la finca en coordenadas GPS para el mapa de detalle:
    rojo = NDVI (0-250), verde = especie del potrero (1, 2, ...), alfa = píxel útil.
    El navegador convierte cada píxel a kg con la curva vigente de su especie.
    """
    especies = list(config["especies"])
    cod_especie = rasterize(
        [(mapping(geoms[p["nombre"]]), especies.index(p["especie"]) + 1) for p in potreros],
        out_shape=ndvi.shape, transform=transf, fill=0, dtype="uint8")
    rojo = (np.clip(np.nan_to_num(ndvi, nan=0), 0, 1) * 250).astype("uint8")
    alfa = np.where((valido | excluido) & (cod_especie > 0), 255, 0).astype("uint8")
    cod_especie = np.where(excluido & (cod_especie > 0), COD_ARBOL, cod_especie).astype("uint8")

    alto, ancho = ndvi.shape
    izq, abajo = transf * (0, alto)
    der, arriba = transf * (ancho, 0)
    dst_t, dst_w, dst_h = calculate_default_transform(epsg, 4326, ancho, alto, izq, abajo, der, arriba)
    bandas = []
    for banda in (rojo, cod_especie, np.zeros_like(rojo), alfa):
        destino = np.zeros((dst_h, dst_w), dtype="uint8")
        reproject(banda, destino, src_transform=transf, src_crs=epsg, dst_transform=dst_t, dst_crs=4326,
                  resampling=ResamplingWarp.nearest)
        bandas.append(destino)
    utiles = int((bandas[3] > 0).sum())
    if not utiles:
        return

    CARPETA_PIXELES.mkdir(parents=True, exist_ok=True)
    ruta_indice = CARPETA_PIXELES / "indice.json"
    indice = json.loads(ruta_indice.read_text(encoding="utf-8")) if ruta_indice.exists() else {}
    previa = indice.get(fecha)
    if previa and previa.get("mascara") == huella and previa["utiles"] >= utiles:
        return  # ya hay una imagen de ese día, con la misma máscara y más píxeles útiles
    with rasterio.open(CARPETA_PIXELES / f"{fecha}.png", "w", driver="PNG", width=dst_w, height=dst_h,
                       count=4, dtype="uint8") as dst:
        dst.write(np.stack(bandas))
    oeste, norte = dst_t * (0, 0)
    este, sur = dst_t * (dst_w, dst_h)
    indice[fecha] = {"limites": [[round(sur, 7), round(oeste, 7)], [round(norte, 7), round(este, 7)]],
                     "utiles": utiles, "mascara": huella}
    ruta_indice.write_text(json.dumps(dict(sorted(indice.items())), indent=1), encoding="utf-8")


def procesar_s2(item, potreros, config, exclusiones=(), huella=""):
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
    with np.errstate(divide="ignore", invalid="ignore"):
        ndvi_img = (n - r) / (n + r)
    validos_img &= np.isfinite(ndvi_img)
    excluido = mascara_excluida(forma, transf, epsg, exclusiones, config)

    filas = []
    for p in potreros:
        g = geoms[p["nombre"]]
        if not g.intersects(box(*bounds_reales)):
            continue
        dentro = geometry_mask([mapping(g)], out_shape=forma, transform=transf, invert=True, all_touched=False)
        n_dentro = int(dentro.sum())
        pasto = dentro & ~excluido  # sin árboles ni bosque
        total = int(pasto.sum())
        if total < config["min_pixeles"]:
            continue
        ok = pasto & validos_img
        pct = ok.sum() / total
        if pct < config["min_pct_valido"]:
            continue
        rr, nn, bb, ee = r[ok], n[ok], b[ok], re[ok]
        with np.errstate(divide="ignore", invalid="ignore"):
            ndvi = np.nanmean(ndvi_img[ok])
            evi = np.nanmean(2.5 * (nn - rr) / (nn + 6 * rr - 7.5 * bb + 1))
            ndre = np.nanmean((nn - ee) / (nn + ee))
            msavi = np.nanmean((2 * nn + 1 - np.sqrt((2 * nn + 1) ** 2 - 8 * (nn - rr))) / 2)
        filas.append({
            "fecha": item.datetime.date().isoformat(), "potrero": p["nombre"], "especie": p["especie"],
            "sensor": "S2", "escena": item.id, "pct_valido": _num(pct, 3),
            "ndvi": _num(ndvi), "evi": _num(evi), "ndre": _num(ndre), "msavi": _num(msavi),
            "hist": histograma(ndvi_img[ok]),
            "pct_arboles": _num(1 - total / n_dentro, 3), "mascara": huella,
        })
    if filas:
        guardar_pixeles(item.datetime.date().isoformat(), ndvi_img, validos_img, excluido, potreros, geoms,
                        transf, epsg, config, huella)
    return filas


def procesar_s1(item, potreros, config, exclusiones=(), huella=""):
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
    validos_img &= ~mascara_excluida(forma, transf, epsg, exclusiones, config)

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
            "sensor": "S1", "escena": item.id, "pct_valido": _num(pct, 3), "mascara": huella,
            "vv_db": _num(10 * np.log10(vv[ok].mean()), 2), "vh_db": _num(10 * np.log10(vh[ok].mean()), 2),
        })
    return filas


def buscar(catalogo, coleccion, bbox, desde, hasta, filtro=None):
    busqueda = catalogo.search(collections=[coleccion], bbox=bbox,
                               datetime=f"{desde.isoformat()}/{hasta.isoformat()}", query=filtro)
    return list(busqueda.items())


# ---------------------------------------------------------------- tablero

def exportar_web(potreros, filas, config, aforos=(), calibracion=None):
    por_potrero, aforos_de = {}, {}
    for r in filas:
        por_potrero.setdefault(r["potrero"], []).append(r)
    for a in aforos:
        aforos_de.setdefault(a["potrero"], []).append(a)
    ventana_base = 2 * config["dias_par_aforo_imagen"]

    feats = []
    for p in potreros:
        serie = sorted(por_potrero.get(p["nombre"], []), key=lambda r: r["fecha"])
        s2 = [{"f": r["fecha"], "fv": float(r["kg_fv_m2"]), "ms": float(r["kg_ms_ha"]),
               "ndvi": float(r["ndvi"]), "v": float(r["pct_valido"]), "pobre": float(r["pct_pobre"]),
               "arb": float(r.get("pct_arboles") or 0)}
              for r in serie if r["sensor"] == "S2" and r.get("hist")]
        s1 = [{"f": r["fecha"], "vv": float(r["vv_db"]), "vh": float(r["vh_db"])}
              for r in serie if r["sensor"] == "S1" and r["vv_db"] != "" and r["vh_db"] != ""]
        mis_aforos = sorted(aforos_de.get(p["nombre"], []), key=lambda a: a["fecha"])

        # Estimación satelital alrededor de cada aforo, para saber cuánto corregir.
        sat_en_aforo = {}
        for a in mis_aforos:
            fa = dt.date.fromisoformat(a["fecha"])
            cerca = [(abs((dt.date.fromisoformat(s["f"]) - fa).days), s["fv"]) for s in s2]
            cerca = [c for c in cerca if c[0] <= ventana_base]
            if cerca:
                sat_en_aforo[a["fecha"]] = min(cerca)[1]
        ajustar_con_aforos(s2, mis_aforos, sat_en_aforo, config)

        # Valor actual: lo más reciente entre la imagen (ajustada si aplica) y el aforo de campo.
        actual, fecha, fuente = None, None, None
        if s2:
            u = s2[-1]
            actual, fecha = (u["aj"], u["f"]) if u["aj"] is not None else (u["fv"], u["f"])
            fuente = "satélite ajustado con aforo" if u["aj"] is not None else "satélite"
        if mis_aforos and (fecha is None or mis_aforos[-1]["fecha"] >= fecha):
            actual, fecha, fuente = mis_aforos[-1]["kg_fv_m2"], mis_aforos[-1]["fecha"], "aforo de campo"

        ha = area_ha(p["geom"])
        frac_arboles = s2[-1]["arb"] if s2 else 0.0
        ha_pasto = ha * (1 - frac_arboles)
        props = {
            "nombre": p["nombre"], "especie": p["especie"],
            "especie_nombre": config["especies"][p["especie"]]["nombre"],
            "area_ha": round(ha, 2),
            "area_pasto_ha": round(ha_pasto, 2),
            "pct_arboles": round(frac_arboles, 3),
            "kg_fv_m2": actual, "ultima_fecha": fecha, "fuente": fuente,
            "kg_fv_m2_satelite": s2[-1]["fv"] if s2 else None,
            "ultima_fecha_satelite": s2[-1]["f"] if s2 else None,
            "kg_ms_ha": s2[-1]["ms"] if s2 else None,
            "ndvi": s2[-1]["ndvi"] if s2 else None,
            "pct_pobre": s2[-1]["pobre"] if s2 else None,
            "t_fv_potrero": round(actual * ha_pasto * 10, 1) if actual is not None else None,
            "serie_s2": s2, "serie_s1": s1,
            "aforos": [{"f": a["fecha"], "kg": a["kg_fv_m2"], "metodo": a.get("metodo", "")} for a in mis_aforos],
        }
        feats.append({"type": "Feature", "properties": props, "geometry": mapping(p["geom"])})

    calibradas = [config["especies"][e]["nombre"] for e, c in (calibracion or {}).items() if c["calibrada"]]
    SALIDA_WEB.parent.mkdir(parents=True, exist_ok=True)
    SALIDA_WEB.write_text(json.dumps({
        "actualizado": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"),
        "modelo": ("calibrado con aforos: " + ", ".join(calibradas)) if calibradas else "provisional sin calibrar",
        "calibracion": calibracion or {},
        "aforos_total": len(aforos),
        "curvas": [dict(id=e, **{c: p[c] for c in ("k", "b", "max_kg_ms_ha", "pct_ms")})
                   for e, p in config["especies"].items()],
        "ndvi_suelo": config["ndvi_suelo"],
        "umbral_pobre": config["umbral_cobertura_pobre_kg"],
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
    # Árboles y bosques: se generan una vez y cualquier cambio en la máscara obliga a reprocesar.
    if not MASCARA_ARBOLES.exists():
        try:
            if not MASCARA_META.exists():
                generar_mascara_arboles(potreros, config)
            MASCARA_ARBOLES.write_bytes(MASCARA_META.read_bytes())
        except Exception as e:
            log(f"  ! No se pudo generar la máscara de árboles (se usan solo los polígonos de exclusiones/): {e}")
    exportar_kmz_arboles()
    exclusiones = cargar_exclusiones()
    huella = huella_mascara(config)
    # Se reprocesan las filas hechas con otra máscara (o antes del cálculo por píxel).
    historial = [r for r in historial if r.get("mascara") == huella and (r["sensor"] != "S2" or r.get("hist"))]
    hechos = {(r["escena"], r["potrero"]) for r in historial}

    hoy = dt.date.today()
    inicio_nuevo = hoy - dt.timedelta(days=config["dias_historial_inicial"])

    def desde_para(sensor):
        """Fecha más vieja que haga falta buscar (potreros nuevos piden todo el historial inicial)."""
        ultimas = {}
        for r in historial:
            if r["sensor"] == sensor:
                ultimas[r["potrero"]] = max(ultimas.get(r["potrero"], r["fecha"]), r["fecha"])
        d = min((dt.date.fromisoformat(ultimas[n]) - dt.timedelta(days=5)) if n in ultimas else inicio_nuevo
                for n in nombres)
        return max(d, inicio_nuevo)

    bbox = list(unary_union([p["geom"] for p in potreros]).bounds)

    catalogo = pystac_client.Client.open(STAC_URL, modifier=planetary_computer.sign_inplace)
    nuevas = []

    desde = desde_para("S2")
    s2 = buscar(catalogo, "sentinel-2-l2a", bbox, desde, hoy,
                {"eo:cloud_cover": {"lt": config["max_nubes_escena"]}})
    log(f"Sentinel-2: {len(s2)} escenas desde {desde}")
    for item in sorted(s2, key=lambda i: i.datetime):
        pendientes = [p for p in potreros if (item.id, p["nombre"]) not in hechos]
        if not pendientes:
            continue
        try:
            filas = procesar_s2(item, pendientes, config, exclusiones, huella)
        except Exception as e:
            log(f"  ! {item.id}: {e}")
            continue
        log(f"  {item.datetime.date()} {item.id[:38]}  nubes {item.properties.get('eo:cloud_cover', 0):.0f}%  -> {len(filas)} potreros útiles")
        nuevas += filas
        # Escena revisada aunque estuviera nublada: se marca para no volver a bajarla.
        hechos.update((item.id, p["nombre"]) for p in pendientes)

    desde = desde_para("S1")
    s1 = buscar(catalogo, "sentinel-1-rtc", bbox, desde, hoy)
    log(f"Sentinel-1: {len(s1)} escenas desde {desde}")
    for item in sorted(s1, key=lambda i: i.datetime):
        pendientes = [p for p in potreros if (item.id, p["nombre"]) not in hechos]
        if not pendientes:
            continue
        try:
            filas = procesar_s1(item, pendientes, config, exclusiones, huella)
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

    # Aforos de campo: los formularios nuevos de GitHub pasan a aforos.csv y se recalibra por especie.
    especie_de = {p["nombre"]: p["especie"] for p in potreros}
    try:
        importar_issues(potreros)
    except Exception as e:
        log(f"  ! No se pudieron leer los aforos de GitHub: {e}")
    aforos = leer_aforos(potreros)
    pares = emparejar(aforos, filas, config)
    calibracion = calibrar(pares, especie_de, config)
    log(f"Aforos de campo: {len(aforos)} ({len(pares)} con imagen cercana). "
        + ", ".join(f"{e}: {c['aforos']}{' calibrada' if c['calibrada'] else ''}" for e, c in calibracion.items()))

    # El modelo se recalcula siempre, así un cambio de calibración se aplica a todo el historial.
    for r in filas:
        r["especie"] = especie_de[r["potrero"]]
        if r["sensor"] == "S2" and r.get("hist"):
            r["kg_ms_ha"], r["kg_fv_m2"], r["pct_pobre"] = estimar_por_pixel(r["hist"], r["especie"], config)

    guardar_historial(filas)
    guardar_calibracion(pares, config, especie_de)
    exportar_web(potreros, filas, config, aforos, calibracion)
    log(f"Listo: {len(nuevas)} observaciones nuevas, {len(filas)} en total.")


# Ajustes para leer imágenes por internet sin descargarlas completas.
OPCIONES_GDAL = dict(
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
    CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.TIF,.tiff",
    GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
    GDAL_HTTP_MULTIPLEX="YES",
    GDAL_HTTP_MAX_RETRY="4",
    GDAL_HTTP_RETRY_DELAY="3",
    GDAL_PAM_ENABLED="NO",  # no crear archivos .aux.xml junto a los PNG
)

if __name__ == "__main__":
    if "--cerrar-issues" in sys.argv:
        cerrar_issues()
    else:
        with rasterio.Env(**OPCIONES_GDAL):
            sys.exit(main())
