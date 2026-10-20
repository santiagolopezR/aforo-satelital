"""
Avisa por GitHub (correo y app del celular) cuando hay resultados nuevos del modelo de machine learning.

    python notificar.py altos_de_medina

Crea un "issue" mencionando al dueño del repositorio con el resumen de cada modelo. Solo avisa una vez por
cada resultado nuevo (recuerda el último en data/<finca>/notificado.txt).
"""
import json
import os
import sys

import requests

import aforo


def main(slug):
    cb = json.loads((aforo.RAIZ / "config.json").read_text(encoding="utf-8"))
    aforo.usar_finca(slug, aforo.fincas_de(cb)[slug])
    ruta = aforo.SALIDA_WEB.with_name("modelo_ml.json")
    marca = aforo.HISTORIAL.parent / "notificado.txt"
    if not ruta.exists():
        print("Todavía no hay resultados del modelo.")
        return
    d = json.loads(ruta.read_text(encoding="utf-8"))
    if marca.exists() and marca.read_text(encoding="utf-8").strip() == d["generado"]:
        print("Ya se avisó de estos resultados.")
        return
    token, repo = os.environ.get("GITHUB_TOKEN"), os.environ.get("GITHUB_REPOSITORY")
    if not token or not repo:
        print("Sin credenciales de GitHub: no se avisa.")
        return
    dueno = repo.split("/")[0]
    usuario, nombre_repo = repo.split("/")
    pagina = f"https://{usuario.lower()}.github.io/{nombre_repo}/modelo/?finca={slug}"
    filas = []
    for m in d["modelos"]:
        p = m.get("prueba_2025_2026")
        filas.append(f"| {m['titulo']} | {m['filas']} | "
                     + (f"{p['r2']:.2f} | {p['rmse']:.2f} | {p['mae_pct']:.0f} % |" if p else f"{m.get('nota', '')} | | |"))
    imp = d["importancia"].get("optico_manejo") or d["importancia"].get("optico") or []
    top = ", ".join(v["nombre"] for v in imp[:5]) or "—"
    cuerpo = (f"@{dueno} ya hay resultados del modelo de aforo para **{d['finca']['nombre']}** "
              f"({d['aforos']} aforos de campo).\n\n"
              "Prueba con aforos de 2025–2026 que el modelo no vio:\n\n"
              "| Datos usados | Aforos | R² | Error típico (kg) | Error medio |\n|---|---|---|---|---|\n"
              + "\n".join(filas)
              + f"\n\n**Variables que más pesan:** {top}.\n\nVer gráficos: {pagina}")
    r = requests.post(f"https://api.github.com/repos/{repo}/issues", timeout=30,
                      headers={"Authorization": f"Bearer {token}", "Accept": "application/vnd.github+json"},
                      json={"title": f"Resultados del modelo de aforo: {d['finca']['nombre']}", "body": cuerpo})
    r.raise_for_status()
    marca.write_text(d["generado"], encoding="utf-8")
    print("Aviso creado:", r.json()["html_url"])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "altos_de_medina")
