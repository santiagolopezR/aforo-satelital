"""
Avisa por GitHub (correo y app del celular) cuando cambian los resultados del modelo de aforo.

    python notificar.py altos_de_medina

Reporta el modelo v2 (crecimiento + memoria + satélite, recalibración mensual). Solo avisa si las cifras
cambian de verdad (R² o error, a 2 decimales); la huella del último aviso queda en data/<finca>/notificado.txt.
"""
import json
import os
import sys

import requests

import aforo


def main(slug):
    cb = json.loads((aforo.RAIZ / "config.json").read_text(encoding="utf-8"))
    aforo.usar_finca(slug, aforo.fincas_de(cb)[slug])
    ruta = aforo.SALIDA_WEB.with_name("modelo_v2.json")
    marca = aforo.HISTORIAL.parent / "notificado.txt"
    if not ruta.exists():
        print("Todavía no hay resultados del modelo v2.")
        return
    d = json.loads(ruta.read_text(encoding="utf-8"))
    huella = ";".join(f"{c['clave']}:{c['franja']['r2']:.2f}/{c['finca_semana']['r2']:.2f}/{c['finca_semana']['rmse']:.2f}"
                      for c in d["configuraciones"])
    if marca.exists() and marca.read_text(encoding="utf-8").strip() == huella:
        print("Los resultados no cambiaron: no se avisa.")
        return
    token, repo = os.environ.get("GITHUB_TOKEN"), os.environ.get("GITHUB_REPOSITORY")
    if not token or not repo:
        print("Sin credenciales de GitHub: no se avisa.")
        return
    usuario, nombre_repo = repo.split("/")
    pagina = f"https://{usuario.lower()}.github.io/{nombre_repo}/modelo/?finca={slug}"
    filas = [f"| {c['titulo']} | {c['finca_semana']['r2']:.2f} | {c['finca_semana']['rmse']:.2f} | "
             f"{c['franja']['r2']:.2f} | {c['franja']['rmse']:.2f} |" for c in d["configuraciones"]]
    cuerpo = (f"@{usuario} cambiaron los resultados del modelo de aforo para **{d['finca']['nombre']}** "
              f"({d['aforos']} aforos). Prueba mes a mes desde {d['desde_prueba']}: cada mes se predice solo con lo anterior.\n\n"
              "| Modelo | R² sub-finca/semana | Error (kg) | R² por franja | Error (kg) |\n|---|---|---|---|---|\n"
              + "\n".join(filas) + f"\n\nVer: {pagina}")
    r = requests.post(f"https://api.github.com/repos/{repo}/issues", timeout=30,
                      headers={"Authorization": f"Bearer {token}", "Accept": "application/vnd.github+json"},
                      json={"title": f"Modelo de aforo actualizado: {d['finca']['nombre']}", "body": cuerpo})
    r.raise_for_status()
    marca.write_text(huella, encoding="utf-8")
    print("Aviso creado:", r.json()["html_url"])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "altos_de_medina")
