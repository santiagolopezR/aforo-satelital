# Aforo satelital de potreros

Estima los **kg de forraje verde por m²** de cada potrero con imágenes satelitales gratuitas
(Sentinel-2 y Sentinel-1). Se actualiza solo todos los días y se ve en un mapa desde el celular.

> **Modelo provisional sin calibrar.** Los kilos son orientativos. Cuando haya aforos de campo,
> se ajustan los números de `config.json` y todo el historial se recalcula.

## Cómo se usa

### Subir o cambiar potreros
1. En GitHub, entrar a la carpeta **`potreros/`** → **Add file → Upload files**.
2. Subir el archivo `.kml`, `.kmz` (Google Earth) o `.geojson`. Puede traer muchos potreros en un solo archivo.
3. Clic en **Commit changes**. En unos minutos el proceso corre solo y el potrero aparece en el mapa
   con los últimos 60 días de historia.

**La especie de cada potrero** se toma de:
- una propiedad `especie` dentro del archivo (en Google Earth: Propiedades → datos extendidos), o
- el archivo **`potreros/especies.csv`**, con esta forma:

  ```
  nombre,especie
  Potrero 1,humidicola
  Potrero 2,mombasa
  La Loma,kikuyo
  El Bajo,decumbens
  ```

Si no se dice nada, se asume la especie de `especie_por_defecto` en `config.json`.

Para **quitar** un potrero, se borra su archivo o se quita del archivo, y su historial desaparece en la siguiente corrida.

### Ver los resultados
- **Mapa y tablero:** la página de GitHub Pages del repositorio (`https://<usuario>.github.io/<repositorio>/`).
- **Datos crudos:** `data/historial.csv` (se puede abrir en Excel).

### Correrlo a mano
Pestaña **Actions → Aforo satelital → Run workflow**.

## Qué hace por dentro
1. Lee los potreros y les quita 10 m de borde para no medir cercas ni caminos.
2. Busca imágenes nuevas en Microsoft Planetary Computer (gratis).
3. **Sentinel-2:** descarta nubes y sombras, y calcula NDVI, EVI, NDRE y MSAVI por potrero.
   Solo usa la imagen si al menos el 60 % del potrero está despejado.
4. **Sentinel-1 (radar):** guarda VV y VH. Por ahora no se usa en los kilos; queda listo para calibrar,
   sobre todo para Mombasa, donde el NDVI se satura.
5. Estima con la curva de abajo (con tope por especie) y `kg FV/m² = kg MS/ha ÷ %MS ÷ 10 000`.

## Aforos de campo (corrigen el satélite)
**Desde el celular:** en el mapa, toca el potrero → **Registrar aforo de este potrero**.
Se abre un formulario de GitHub con el nombre y la fecha ya puestos; solo se escriben los kg FV/m²
(0 si el potrero está sin pasto) y se envía. En 1–2 minutos el mapa se actualiza y el formulario
se cierra con un comentario que compara el dato de campo con lo que decía el satélite.

Solo cuentan los aforos del dueño del repositorio y de sus colaboradores
(Settings → Collaborators para agregar al mayordomo u otra persona con cuenta de GitHub).

También se pueden escribir directamente en `aforos_campo/aforos.csv`.

**Qué hace cada aforo:**
1. **Corrige ese potrero de inmediato.** El día del aforo el mapa muestra lo medido; en las imágenes
   siguientes la diferencia con el satélite se va desvaneciendo en 30 días (`dias_efecto_aforo`).
2. **Calibra la especie.** Cuando una especie junta 10 aforos con imagen despejada a ±5 días,
   el modelo de esa especie se recalcula solo (y el % de materia seca, si se anotó en 5 o más).
   La comparación aforo vs. satélite queda en `data/calibracion.csv`.

**Modelo:** `kg MS/ha = k · (e^(b · (NDVI − 0,25)) − 1)`; con NDVI de suelo desnudo (≤ 0,25) da 0.

## Pastoreo: con cuántos kilos va a entrar
Cada potrero se ubica en su ciclo y se proyecta su forraje a la fecha de entrada.

- **Salida de los animales:** la registrada; si no hay, la que ve el satélite (bajón fuerte entre dos
  imágenes despejadas). Si la salida conocida es de hace más de un ciclo, se marca como "salida vieja"
  y no se inventa la entrada.
- **Entrada:** salida + días de descanso de la rotación (`config.json` → `pastoreo.descanso_por_rotacion`, hoy 30 en todas).
- **Rebrote:** curva logística por especie aprendida de las imágenes de la finca (se recalcula en cada corrida).
  Especies con menos de 10 pares de imágenes usan la curva de las demás.
- **Proyección:** desde el último dato (satélite, aforo o residuo de 0,3 kg al salir) hasta la entrada,
  con un rango de ±40 % en la velocidad de rebrote.

**En el mapa:** selector "Colores: Hoy / En 7, 15, 30 días / A la entrada"; en cada potrero, la entrada
proyectada y una casilla "¿Cómo estará en N días?".

**Registrar entradas y salidas:** botón "Registrar entrada / salida de animales" en cada potrero
(mismo sistema que los aforos), o a mano en `pastoreos/pastoreos.csv`
(`fecha,potrero,movimiento,observaciones,origen`, movimiento = entrada o salida).

## Árboles y bosques (no cuentan como pasto)
Un cuadrito de 10 m con más de la mitad bajo árboles se saca del NDVI, de los kilos y del área en pasto.
En el mapa de detalle se pintan en verde intenso.

- **`exclusiones/arboles.tif`**: la máscara que se usa (1,2 m). Se arma con:
  - **Ortomosaicos de dron** (donde hay; el más reciente manda) con `arboles_dron.py`.
  - **Altura de árboles de Meta** (`exclusiones/arboles_meta.tif`, ≥ 3 m) donde no hay dron, p. ej. Primavera.
- **`exclusiones/bosques.geojson`**: los bosques de la capa QGIS.
- **Cualquier otro polígono** (KML/KMZ/GeoJSON) que se suba a `exclusiones/` también se excluye.

**Revisar en Google Earth:** en el mapa web, "Descargar árboles detectados para Google Earth (KMZ)".
Si sobra un árbol, avisar; si falta, dibujarlo en Google Earth y subir el KMZ a `exclusiones/`.

**Con un ortomosaico nuevo** (en el PC donde están las imágenes, que nunca se suben):

```
python arboles_dron.py "D:/MAPA Puerto parra Claude/03_ortomosaicos"
```

y subir `exclusiones/arboles.tif` y `docs/arboles.kmz`. El resto se recalcula solo.

## Varias fincas
Cada finca se define en `config.json` → `fincas` (carpeta de potreros, de datos, de la página, de exclusiones,
archivos de aforos y ajustes propios). En el mapa se cambia con el selector junto al título (`?finca=...`).

| Finca | Potreros | Pasto | Notas |
|---|---|---|---|
| La Holanda (Puerto Parra) | `potreros/` | humidícola, decumbens | rutas de siempre (`data/`, `docs/`, `exclusiones/`) |
| Altos de Medina | `potreros/altos_de_medina/` | kikuyo | franjas pequeñas: recorte de borde de 5 m |

## Calibración con aforos históricos (Altos de Medina)
`aforos_campo/altos_de_medina_platometro.csv`: 5.054 aforos de platómetro (nov 2022 – ago 2026), 4.792 cruzan con un potrero.
El platómetro se toma al entrar las vacas, así que cada aforo se empareja con la imagen despejada (≥ 80 %)
**de ese día o hasta 5 días antes**, nunca después.

`historico.py` trae las imágenes Sentinel-2 viejas de esas fechas (workflow "Historico para calibrar", corre cada 6 h
hasta terminar; es incremental). Con eso `aforo.py` ajusta la curva del kikuyo (escala `k` y forma `b`) y mide el error
con validación cruzada: `data/altos_de_medina/calibracion.csv` y el bloque `calibracion` de su `datos.json`.

## Modelo de machine learning
`dataset_ml.py` arma una fila por aforo con lo medido antes de él: Sentinel-2 (10 bandas y 12 índices, con percentiles),
Sentinel-1 radar (VV, VH, VH−VV, RVI), Landsat 8/9 (incluida la temperatura de superficie), clima diario
(Open-Meteo/ERA5: lluvia, temperatura, evapotranspiración, radiación) y manejo (días de descanso, mes).
`entrenar_ml.py` entrena modelos (gradient boosting) con los aforos hasta 2024 y los prueba con 2025–2026;
también deja potreros completos por fuera. Compara: solo NDVI, todas las bandas, + radar/Landsat/clima, + manejo,
solo radar (días nublados) y un control sin satélite. Resultados en la página `docs/modelo/` (R², error, real vs estimado
y qué variables pesan más). Todo corre en el workflow "Historico para calibrar".

## Limitaciones
- El píxel es de 10 m: potreros de menos de ~0,3 ha dan datos pobres.
- Con muchas semanas nubladas, el dato queda viejo (el tablero avisa si tiene más de 10 días).
- Maleza y pasto seco se confunden con forraje.
