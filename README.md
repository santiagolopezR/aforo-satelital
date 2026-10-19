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
5. Estima: `kg MS/ha = a · e^(b · NDVI)` (con tope por especie) y `kg FV/m² = kg MS/ha ÷ %MS ÷ 10 000`.

## Calibración (más adelante)
Anotar cada aforo de campo en `aforos_campo/aforos.csv`:

```
fecha,potrero,kg_fv_m2,pct_ms,observaciones
2026-11-02,Potrero 1,1.85,0.27,5 marcos de 0.25 m2
```

Con 30 o más aforos por especie se ajustan `a`, `b` y `pct_ms` en `config.json`.

## Limitaciones
- El píxel es de 10 m: potreros de menos de ~0,3 ha dan datos pobres.
- Con muchas semanas nubladas, el dato queda viejo (el tablero avisa si tiene más de 10 días).
- Maleza y pasto seco se confunden con forraje.
