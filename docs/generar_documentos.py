"""Genera la documentación formal del programa sin dependencias externas.

Uso:
    python docs/generar_documentos.py

Salidas:
    docs/flujograma.pdf
    docs/flujograma.svg
    docs/flujograma.png
    docs/manual_usuario.pdf
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon


ROOT = Path(__file__).resolve().parent
BLUE = "#12355b"
CYAN = "#0e7490"
LIGHT = "#e8f1f5"
GREEN = "#dcfce7"
RED = "#fee2e2"
AMBER = "#fef3c7"
TEXT = "#17202a"


def page_base(title: str, page: int, total: int | None = None):
    fig = plt.figure(figsize=(11.69, 8.27), facecolor="white")
    fig.text(0.06, 0.95, "PROGRAMA DE DISEÑO DE COLUMNAS", color=BLUE, fontsize=9,
             fontweight="bold", va="top")
    fig.text(0.94, 0.95, "ACI 318-19", color=CYAN, fontsize=9, ha="right", va="top")
    fig.add_artist(plt.Line2D([0.06, 0.94], [0.935, 0.935], color=CYAN, lw=1.2))
    fig.text(0.06, 0.895, title, color=BLUE, fontsize=18, fontweight="bold", va="top")
    footer = f"Documento de usuario · Versión 1.0 · 23 de julio de 2026    |    {page}"
    fig.text(0.06, 0.035, footer, color="#64748b", fontsize=8, va="bottom")
    return fig


def rounded_box(ax, center, width, height, label, color=LIGHT, edge=BLUE, fontsize=9):
    x, y = center
    patch = FancyBboxPatch((x - width / 2, y - height / 2), width, height,
                           boxstyle="round,pad=0.012,rounding_size=0.025",
                           facecolor=color, edgecolor=edge, linewidth=1.4)
    ax.add_patch(patch)
    ax.text(x, y, label, ha="center", va="center", color=TEXT, fontsize=fontsize,
            wrap=True, linespacing=1.15)


def decision(ax, center, width, height, label, color=AMBER, fontsize=8.5):
    x, y = center
    vertices = [(x, y + height / 2), (x + width / 2, y),
                (x, y - height / 2), (x - width / 2, y)]
    ax.add_patch(Polygon(vertices, closed=True, facecolor=color, edgecolor=BLUE, linewidth=1.4))
    ax.text(x, y, label, ha="center", va="center", color=TEXT, fontsize=fontsize,
            wrap=True, linespacing=1.1)


def arrow(ax, start, end, label=None, color=BLUE, dy=0.0):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12,
                                 linewidth=1.3, color=color, connectionstyle="arc3,rad=0"))
    if label:
        mx, my = (start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + dy
        ax.text(mx, my, label, fontsize=8, color=color, ha="center", va="center",
                bbox=dict(facecolor="white", edgecolor="none", pad=1.5))


def flowchart_streamlit():
    fig = page_base("Flujograma principal · análisis Streamlit", 1)
    ax = fig.add_axes([0.045, 0.08, 0.91, 0.78])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    rounded_box(ax, (0.5, 0.95), 0.16, 0.055, "INICIO", color=GREEN)
    rounded_box(ax, (0.5, 0.85), 0.35, 0.075,
                "Ingresar B, H, recubrimiento, diámetro de estribo,\n"
                "barras longitudinales, f'c, fy y demandas Pu–Mux–Muy")
    rounded_box(ax, (0.5, 0.74), 0.20, 0.05, "Presionar CALCULAR", color="#dbeafe")
    decision(ax, (0.5, 0.62), 0.25, 0.10, "¿Geometría y\nmateriales válidos?")
    rounded_box(ax, (0.17, 0.62), 0.22, 0.065, "Mostrar error y\nsolicitar corrección", color=RED)
    rounded_box(ax, (0.5, 0.48), 0.31, 0.07,
                "Generar matriz perimetral\ny ubicar barras en la sección")
    decision(ax, (0.5, 0.36), 0.27, 0.105, "¿Detalle cumple?\nS libre ≥ 2.50 cm\nEstribo y ramas válidos")
    rounded_box(ax, (0.17, 0.36), 0.22, 0.065, "Solicitar B, recubrimiento\no diámetros válidos", color=RED)
    decision(ax, (0.5, 0.23), 0.23, 0.09, "¿1% ≤ ρ ≤ 8%?")
    rounded_box(ax, (0.17, 0.23), 0.22, 0.065, "Error si ρ > 8%\nAdvertencia si ρ < 1%", color=RED)
    rounded_box(ax, (0.5, 0.11), 0.34, 0.07,
                "Calcular P–Mx, P–My, contorno biaxial ACI\ny verificaciones de cortante", color="#dbeafe")
    rounded_box(ax, (0.82, 0.23), 0.23, 0.07,
                "Mostrar gráficos,\nD/C y memoria", color=GREEN)
    rounded_box(ax, (0.82, 0.11), 0.15, 0.05, "FIN", color=GREEN)

    arrow(ax, (0.5, 0.922), (0.5, 0.89))
    arrow(ax, (0.5, 0.812), (0.5, 0.765))
    arrow(ax, (0.5, 0.715), (0.5, 0.67))
    arrow(ax, (0.39, 0.62), (0.28, 0.62), "NO", dy=0.018)
    arrow(ax, (0.5, 0.57), (0.5, 0.515), "SÍ", dy=0.014)
    arrow(ax, (0.5, 0.445), (0.5, 0.415))
    arrow(ax, (0.39, 0.36), (0.28, 0.36), "NO", dy=0.018)
    arrow(ax, (0.5, 0.307), (0.5, 0.275), "SÍ", dy=0.014)
    arrow(ax, (0.39, 0.23), (0.28, 0.23), "NO", dy=0.018)
    arrow(ax, (0.5, 0.185), (0.5, 0.15), "SÍ", dy=0.014)
    arrow(ax, (0.67, 0.23), (0.75, 0.23))
    arrow(ax, (0.82, 0.195), (0.82, 0.14))
    arrow(ax, (0.28, 0.62), (0.28, 0.76), color="#b91c1c")
    arrow(ax, (0.28, 0.36), (0.28, 0.50), color="#b91c1c")
    arrow(ax, (0.28, 0.23), (0.28, 0.30), color="#b91c1c")
    return fig


def flowchart_wsgi():
    fig = page_base("Flujograma alternativo · servidor WSGI / Railway", 2)
    ax = fig.add_axes([0.06, 0.15, 0.88, 0.65])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    rounded_box(ax, (0.12, 0.5), 0.17, 0.13, "Solicitud HTTP\n/?b=40&h=40&…", color="#dbeafe")
    rounded_box(ax, (0.36, 0.5), 0.18, 0.13, "wsgi_app()\nLee QUERY_STRING")
    decision(ax, (0.57, 0.5), 0.16, 0.16, "¿/health?")
    rounded_box(ax, (0.78, 0.72), 0.18, 0.12, "200 OK\n{status: ok}", color=GREEN)
    rounded_box(ax, (0.78, 0.29), 0.24, 0.15,
                "Calcula el mismo núcleo\nque Streamlit y genera\nHTML + PNG base64")
    rounded_box(ax, (0.78, 0.08), 0.15, 0.08, "Respuesta HTML", color=GREEN)
    arrow(ax, (0.205, 0.5), (0.265, 0.5))
    arrow(ax, (0.45, 0.5), (0.49, 0.5))
    arrow(ax, (0.63, 0.56), (0.69, 0.68), "SÍ", dy=0.01)
    arrow(ax, (0.65, 0.45), (0.68, 0.34), "NO", dy=0.01)
    arrow(ax, (0.78, 0.215), (0.78, 0.13))
    fig.text(0.08, 0.10, "El modo WSGI comparte el motor numérico con Streamlit y permite una respuesta HTML autónoma.",
             fontsize=9, color="#475569")
    return fig


def text_block(fig, x, y, width, text, fontsize=10, color=TEXT, line_gap=0.024):
    lines = []
    for paragraph in text.split("\n"):
        if not paragraph.strip():
            lines.append("")
        else:
            lines.extend(textwrap.wrap(paragraph, width=width, break_long_words=False))
    for line in lines:
        fig.text(x, y, line, fontsize=fontsize, color=color, va="top")
        y -= line_gap
    return y


def manual_pdf():
    path = ROOT / "manual_usuario.pdf"
    with PdfPages(path) as pdf:
        fig = page_base("Manual de usuario", 1)
        fig.text(0.08, 0.69, "Programa de diseño de columnas\nde concreto armado", fontsize=28,
                 color=BLUE, fontweight="bold", va="top", linespacing=1.15)
        fig.text(0.08, 0.48, "Análisis de flexocompresión uniaxial, biaxial y cortante\nsegún ACI 318-19",
                 fontsize=15, color=CYAN, va="top", linespacing=1.4)
        fig.text(0.08, 0.28, "Manual formal de operación\nVersión 1.0 · 23 de julio de 2026",
                 fontsize=11, color="#475569", va="top", linespacing=1.5)
        fig.text(0.92, 0.12, "Documento de apoyo\npara revisión del ingeniero estructural",
                 fontsize=10, ha="right", color="#64748b", va="bottom")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        fig = page_base("1. Alcance y puesta en marcha", 2)
        y = 0.84
        y = text_block(fig, 0.08, y, 105,
                       "El programa analiza columnas rectangulares de concreto armado con estribos. "
                       "Calcula diagramas de interacción P–M, contorno biaxial ACI para una carga axial dada, "
                       "verificaciones de cuantía, disposición del acero y cortante. Los resultados son una ayuda "
                       "de diseño y deben ser revisados por un profesional responsable.", line_gap=0.026)
        fig.text(0.08, y - 0.02, "Ejecución local", fontsize=13, fontweight="bold", color=BLUE)
        text_block(fig, 0.08, y - 0.06, 105,
                   "1. Instale Python 3.11 o superior.\n"
                   "2. En la carpeta del proyecto ejecute: python -m pip install -r requirements.txt\n"
                   "3. Inicie la aplicación: streamlit run app.py\n"
                   "4. Abra http://localhost:8501 en el navegador.", line_gap=0.026)
        fig.text(0.55, 0.60, "Unidades de entrada", fontsize=13, fontweight="bold", color=BLUE)
        rows = [("B, H, rec, L, s", "cm"), ("Diámetros", "mm"), ("f'c, fy", "kgf/cm²"),
                ("Pu", "tonf"), ("Mux, Muy", "tonf·m")]
        table = fig.add_axes([0.55, 0.30, 0.34, 0.27])
        table.axis("off")
        table.table(cellText=rows, colLabels=["Entrada", "Unidad"], loc="center", cellLoc="left")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        fig = page_base("2. Datos de entrada", 3)
        y = 0.84
        sections = [
            ("Geometría", "B es la base, H el peralte, rec el recubrimiento libre y L la longitud libre de la columna. El recubrimiento debe permitir que barras y estribos quepan dentro de la sección."),
            ("Estribos", "Ingrese el diámetro real del estribo en milímetros. El programa calcula el área transversal y verifica las restricciones de confinamiento y separación."),
            ("Acero longitudinal", "Ingrese el número de barras en B y H, el diámetro de las barras de borde y el diámetro de las barras de esquina. La matriz se genera en el perímetro; el núcleo queda sin barras."),
            ("Materiales", "f'c es la resistencia especificada del concreto y fy la fluencia del acero. Se introducen en kgf/cm²."),
            ("Demandas", "Pu es la carga axial última; Mux y Muy son los momentos últimos alrededor de los ejes X y Y."),
        ]
        for title, body in sections:
            fig.text(0.08, y, title, fontsize=12, fontweight="bold", color=BLUE, va="top")
            y = text_block(fig, 0.08, y - 0.035, 110, body, fontsize=10, line_gap=0.024) - 0.02
        fig.text(0.08, 0.10, "Nota: el espaciamiento libre entre barras longitudinales se calcula automáticamente; no se solicita como dato independiente.",
                 fontsize=9, color="#92400e", bbox=dict(facecolor=AMBER, edgecolor="#f59e0b", pad=6))
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        fig = page_base("3. Lectura de resultados", 4)
        y = 0.84
        blocks = [
            ("Geometría y propiedades", "Ag, Ast, cuantía ρ, P0, Pn,max y φPn,max. La cuantía de referencia debe mantenerse entre 1% y 8%."),
            ("Diagramas P–M", "La curva negra representa resistencia nominal; la curva azul, resistencia de diseño φ; la línea naranja marca φPn,max y el punto rojo corresponde a la demanda."),
            ("Contorno biaxial ACI", "El contorno φMx–φMy se genera para el nivel Pu mediante compatibilidad de deformaciones y variación de la orientación del eje neutro. El punto de demanda debe quedar dentro de la región resistente; D/C ≤ 1.00 indica cumplimiento."),
            ("Verificación ACI", "Se revisan geometría, acero longitudinal, cuantía, acero transversal, longitud de confinamiento, ramales, separación de estribos y separación libre de barras."),
            ("Cortante", "Se muestran Vc, Vs, φVn y la separación máxima requerida. Verifique además los límites de detallado y las condiciones sísmicas que correspondan al proyecto."),
        ]
        for title, body in blocks:
            fig.text(0.08, y, title, fontsize=12, fontweight="bold", color=BLUE, va="top")
            y = text_block(fig, 0.08, y - 0.035, 110, body, fontsize=10, line_gap=0.024) - 0.025
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        fig = page_base("4. Errores y criterios de aceptación", 5)
        y = 0.84
        y = text_block(fig, 0.08, y, 110,
                       "El cálculo se detiene cuando los datos no son físicamente válidos o el armado no cumple los controles geométricos. Corrija los datos de entrada y vuelva a pulsar CALCULAR.", line_gap=0.026)
        fig.text(0.08, y - 0.02, "Situaciones habituales", fontsize=13, fontweight="bold", color=BLUE)
        items = [
            "B, H o recubrimiento incompatibles: aumente la sección o reduzca el recubrimiento dentro de los límites del proyecto.",
            "Espaciamiento libre menor a 2.50 cm: redistribuya las barras, aumente B/H o reduzca el diámetro.",
            "Cuantía menor a 1%: aumente el acero longitudinal; cuantía mayor a 8%: reduzca acero o aumente la sección.",
            "Pu fuera del contorno: aumente sección, resistencia del material o acero y revise esbeltez/segundo orden.",
            "D/C biaxial mayor que 1.00: la combinación Mux–Muy no es resistente para Pu; ajuste el diseño.",
            "Estribos insuficientes: aumente diámetro, número de ramales o reduzca la separación según el detalle requerido.",
        ]
        y -= 0.07
        for item in items:
            fig.text(0.10, y, "• " + item, fontsize=10, color=TEXT, va="top")
            y -= 0.055
        fig.text(0.08, 0.12, "Advertencia profesional: este software no sustituye la memoria de cálculo, la revisión de cargas, la esbeltez, el segundo orden, el detallado sísmico ni la firma del ingeniero estructural.",
                 fontsize=9, color="#991b1b", bbox=dict(facecolor=RED, edgecolor="#dc2626", pad=7), wrap=True)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def main():
    fig1 = flowchart_streamlit()
    fig1.savefig(ROOT / "flujograma.svg", format="svg", bbox_inches="tight")
    fig1.savefig(ROOT / "flujograma.png", format="png", dpi=220, bbox_inches="tight")
    with PdfPages(ROOT / "flujograma.pdf") as pdf:
        pdf.savefig(fig1, bbox_inches="tight")
        plt.close(fig1)
        fig2 = flowchart_wsgi()
        pdf.savefig(fig2, bbox_inches="tight")
        plt.close(fig2)
    manual_pdf()
    print("Documentos generados en", ROOT)


if __name__ == "__main__":
    main()
