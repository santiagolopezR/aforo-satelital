"""
Detecta árboles con los ortomosaicos de dron de la finca y actualiza la máscara de árboles.

Se corre a mano en el PC donde están los ortomosaicos (las imágenes nunca se suben):

    python arboles_dron.py "D:/MAPA Puerto parra Claude/03_ortomosaicos"

Cómo funciona, por cada ortomosaico (el más reciente manda donde se traslapan):
  1. Lo lleva a una malla de ~0,6 m en la misma cuadrícula de la máscara de Meta.
  2. Toma ejemplos seguros de la máscara de Meta: centro de copas (árbol) y pasto lejos de árboles.
  3. Entrena un modelo (bosque aleatorio) con color, brillo y textura de SUS fotos y clasifica cada píxel.
  4. Limpia manchas sueltas y pasa a la malla de 1,2 m de la máscara.
Donde no hay ortomosaico (p. ej. Primavera) queda la máscara de Meta.

Salidas: exclusiones/arboles.tif (máscara final) y docs/arboles.kmz (para revisar en Google Earth).
"""
import re
import sys
import time
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.vrt import WarpedVRT
from scipy import ndimage
from sklearn.ensemble import RandomForestClassifier

import aforo

RES_LECTURA_MAX = 0.30  # m: se usa la vista reducida del ortomosaico más cercana a esto
MUESTRAS_POR_CLASE = 60000


def log(*a):
    print(*a, flush=True)


def fecha_de(ruta):
    m = re.search(r"(20\d\d)-(\d\d)", ruta.name)
    return m.group(0) if m else "0000-00"


def nivel_vista(src):
    """Elige la vista reducida (overview) con resolución <= RES_LECTURA_MAX."""
    res_m = src.res[0] * (111320 if src.crs.is_geographic else 1)
    mejor = None
    for i, factor in enumerate(src.overviews(1)):
        if res_m * factor <= RES_LECTURA_MAX:
            mejor = i
    return mejor


def leer_en_malla(ruta, transf, alto, ancho):
    """Lee el ortomosaico (RGBA) remuestreado a la malla fina; devuelve (rgb float32, cubierto bool)."""
    with rasterio.open(ruta) as base:
        nivel = nivel_vista(base)
    abrir = dict(overview_level=nivel) if nivel is not None else {}
    with rasterio.open(ruta, **abrir) as src:
        with WarpedVRT(src, crs="EPSG:3857", transform=transf, width=ancho, height=alto,
                       resampling=Resampling.average) as vrt:
            datos = vrt.read()
    rgb = datos[:3].astype("float32")
    cubierto = (datos[3] > 200) & (rgb.sum(axis=0) > 0)
    return rgb, cubierto


def caracteristicas(rgb):
    """Color normalizado, brillo, índice de verdor y textura en vecindarios de ~3 m y ~6,6 m."""
    R, G, B = rgb
    suma = R + G + B + 1
    r, g, b = R / suma, G / suma, B / suma
    brillo = (R + G + B) / 3
    exg = 2 * g - r - b
    capas = [r, g, b, brillo, exg]
    for tam in (5, 11):
        media = ndimage.uniform_filter(brillo, tam)
        cuad = ndimage.uniform_filter(brillo * brillo, tam)
        capas += [media, np.sqrt(np.maximum(cuad - media * media, 0)), ndimage.uniform_filter(exg, tam)]
    return np.stack(capas, axis=-1)


def limpiar(arbol, min_pixeles=6):
    """Quita ruido: apertura morfológica y manchas de menos de ~2 m²."""
    arbol = ndimage.binary_opening(arbol, iterations=1)
    etiquetas, n = ndimage.label(arbol)
    if n:
        tamanos = ndimage.sum(arbol, etiquetas, range(1, n + 1))
        arbol = np.isin(etiquetas, np.flatnonzero(tamanos >= min_pixeles) + 1)
    return arbol


def main(carpeta):
    t_ini = time.time()
    meta_ruta = aforo.CARPETA_EXCLUSIONES / "arboles_meta.tif"
    with rasterio.open(meta_ruta) as m:
        meta = m.read(1).astype(bool)
        t_meta, perfil = m.transform, m.profile

    # Malla fina: la mitad del píxel de la máscara, alineada con ella.
    t_fina = Affine(t_meta.a / 2, 0, t_meta.c, 0, t_meta.e / 2, t_meta.f)
    alto, ancho = meta.shape[0] * 2, meta.shape[1] * 2
    meta_fina = np.repeat(np.repeat(meta, 2, axis=0), 2, axis=1)
    # Ejemplos seguros: centro de copas de Meta y pasto a más de ~10 m de cualquier árbol de Meta.
    seguro_arbol = ndimage.binary_erosion(meta_fina, iterations=3)
    seguro_pasto = ~ndimage.binary_dilation(meta_fina, iterations=16)

    resultado = np.zeros((alto, ancho), dtype=bool)
    hecho = np.zeros((alto, ancho), dtype=bool)
    ortos = sorted(Path(carpeta).glob("*.tif"), key=fecha_de, reverse=True)
    rng = np.random.default_rng(0)
    for ruta in ortos:
        t0 = time.time()
        rgb, cubierto = leer_en_malla(ruta, t_fina, alto, ancho)
        usar = cubierto & ~hecho
        if usar.sum() < 1000:
            log(f"{ruta.name}: sin área nueva, se omite")
            continue
        X = caracteristicas(rgb)
        idx_a = np.flatnonzero((usar & seguro_arbol).ravel())
        idx_p = np.flatnonzero((usar & seguro_pasto).ravel())
        if len(idx_a) < 500 or len(idx_p) < 500:
            log(f"{ruta.name}: muy pocos ejemplos ({len(idx_a)} árbol / {len(idx_p)} pasto), se deja la máscara de Meta")
            hecho |= usar
            resultado[usar] = meta_fina[usar]
            continue
        idx_a = rng.choice(idx_a, min(MUESTRAS_POR_CLASE, len(idx_a)), replace=False)
        idx_p = rng.choice(idx_p, min(MUESTRAS_POR_CLASE, len(idx_p)), replace=False)
        Xf = X.reshape(-1, X.shape[-1])
        modelo = RandomForestClassifier(n_estimators=80, max_depth=14, min_samples_leaf=20, n_jobs=-1, random_state=0)
        modelo.fit(np.vstack([Xf[idx_a], Xf[idx_p]]), np.r_[np.ones(len(idx_a)), np.zeros(len(idx_p))])

        prob = np.zeros(alto * ancho, dtype="float32")
        idx_usar = np.flatnonzero(usar.ravel())
        for i in range(0, len(idx_usar), 2_000_000):
            bloque = idx_usar[i:i + 2_000_000]
            prob[bloque] = modelo.predict_proba(Xf[bloque])[:, 1]
        arbol = limpiar(prob.reshape(alto, ancho) >= 0.5) & usar
        resultado[usar] = arbol[usar]
        hecho |= usar
        log(f"{ruta.name}: {usar.sum() * (t_fina.a ** 2) / 10000:.0f} ha nuevas, "
            f"{arbol[usar].mean() * 100:.1f} % árbol (Meta: {meta_fina[usar].mean() * 100:.1f} %), {time.time() - t0:.0f} s")
        del X, Xf, prob, rgb

    # Donde no hubo ortomosaico queda Meta; se pasa de 0,6 m a la malla de 1,2 m (mayoría de 4 subpíxeles).
    resultado[~hecho] = meta_fina[~hecho]
    final = resultado.reshape(meta.shape[0], 2, meta.shape[1], 2).mean(axis=(1, 3)) >= 0.5
    perfil.update(nbits=1, compress="deflate")
    with rasterio.open(aforo.MASCARA_ARBOLES, "w", **perfil) as dst:
        dst.write(final.astype("uint8"), 1)
    cob = hecho.reshape(meta.shape[0], 2, meta.shape[1], 2).any(axis=(1, 3))
    log(f"Máscara final: {final.mean() * 100:.1f} % árbol (Meta sola: {meta.mean() * 100:.1f} %). "
        f"Ortomosaico en {cob.mean() * 100:.0f} % del área. {time.time() - t_ini:.0f} s")
    aforo.exportar_kmz_arboles()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
