"""Aplicación Streamlit para diseño preliminar de columnas de concreto armado.

Unidades internas: kgf y m (esfuerzos en kgf/m²). Las entradas y los resultados
mostrados usan cm, kgf/cm² y tonf para conservar la práctica habitual de diseño.
"""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from typing import Iterable
from urllib.parse import parse_qs

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

# Sistema base MKS solicitado (la unidad de fuerza usada es kgf).
m = 1.0
cm = 0.01 * m
mm = 0.001 * m
kgf = 1.0
ton = 1000.0 * kgf


class GeometriaError(ValueError):
    """Geometría físicamente imposible para la sección o el armado."""


class CuantiaError(ValueError):
    """Cuantía longitudinal fuera del intervalo de control de la aplicación."""


class BreslerError(ValueError):
    """No se puede resolver la comprobación biaxial para la demanda dada."""


class Rebar:
    """Catálogo de barras ASTM; áreas en m² y diámetros en m."""

    _diametros = {
        2: 6.4 * mm, 3: 9.5 * mm, 4: 12.7 * mm, 5: 15.9 * mm,
        6: 19.1 * mm, 7: 22.2 * mm, 8: 25.4 * mm, 10: 31.8 * mm,
    }
    _areas = {
        2: 0.32 * cm**2, 3: 0.71 * cm**2, 4: 1.27 * cm**2,
        5: 1.98 * cm**2, 6: 2.85 * cm**2, 7: 3.88 * cm**2,
        8: 5.10 * cm**2, 10: 7.92 * cm**2,
    }

    def __init__(self, matrix: Iterable[Iterable[int]]):
        self.matrix = np.asarray(matrix, dtype=int)
        if self.matrix.ndim != 2 or min(self.matrix.shape) < 1:
            raise GeometriaError("La matriz de barras debe tener por lo menos una fila y una columna.")
        invalidas = set(np.unique(self.matrix)) - {0, *self._areas}
        if invalidas:
            raise ValueError(f"Números de barra no admitidos: {sorted(invalidas)}.")

    def diametro(self, n: int) -> float:
        return self._diametros.get(int(n), 0.0)

    def area(self, n: int) -> float:
        return self._areas.get(int(n), 0.0)

    def define_acero(self) -> np.ndarray:
        return np.vectorize(self.area, otypes=[float])(self.matrix)


class VariablesGenerales:
    """Relaciones constitutivas de ACI 318-19 para concreto y acero."""

    ES = 2.0e6 * kgf / cm**2

    def esfuerzo_acero(self, fy: float, es: float) -> float:
        ey = fy / self.ES
        return float(np.clip(self.ES * es, -fy, fy))

    def beta1(self, fc: float) -> float:
        # fc está en kgf/m²: se convierte sólo para la ecuación tabulada ACI.
        fc_kgcm2 = fc * cm**2 / kgf
        if fc_kgcm2 <= 280:
            return 0.85
        if fc_kgcm2 >= 550:
            return 0.65
        return 0.85 - 0.05 * (fc_kgcm2 - 280) / 70.0


@dataclass(frozen=True)
class PuntoInteraccion:
    c: float
    pn: float
    mn: float
    phi: float
    et: float

    @property
    def pp(self) -> float:
        return self.phi * self.pn

    @property
    def mp(self) -> float:
        return self.phi * self.mn


@dataclass(frozen=True)
class ResultadoCortante:
    vc: float
    vs: float
    phi_vn: float
    s_requerido: float | None


class DiagramaInteraccion:
    """Compatibilidad de deformaciones de una sección rectangular armada.

    ``rebar_matrix`` representa la sección vista en planta: filas de arriba a
    abajo y columnas de izquierda a derecha. El cero indica una celda sin barra.
    """

    def __init__(self, b: float, h: float, rec: float, estribo: int,
                 fc: float, fy: float, rebar_matrix: Iterable[Iterable[int]]):
        self.b, self.h, self.rec = float(b), float(h), float(rec)
        self.fc, self.fy = float(fc), float(fy)
        self.rebar = Rebar(rebar_matrix)
        self.estribo = int(estribo)
        self.tie_diameter = self.rebar.diametro(self.estribo)
        self.vars = VariablesGenerales()
        self.centroide = self.h / 2.0
        self._validate()
        self._bars = self._bar_coordinates()

    def _validate(self) -> None:
        if self.b <= 0 or self.h <= 0 or self.rec <= 0:
            raise GeometriaError("B, H y el recubrimiento deben ser mayores que cero.")
        if 2 * self.rec >= self.b or 2 * self.rec >= self.h:
            raise GeometriaError("El recubrimiento deja sin núcleo útil a la sección.")
        if self.fc <= 0 or self.fy <= 0:
            raise ValueError("f'c y fy deben ser mayores que cero.")
        if self.tie_diameter <= 0:
            raise ValueError("El número de barra seleccionado para el estribo no es válido.")
        if not np.any(self.rebar.matrix > 0):
            raise GeometriaError("Debe existir al menos una barra longitudinal.")

    @property
    def ag(self) -> float:
        return self.b * self.h

    @property
    def ast(self) -> float:
        return float(self.rebar.define_acero().sum())

    @property
    def rho(self) -> float:
        return self.ast / self.ag

    @property
    def p0(self) -> float:
        return 0.85 * self.fc * (self.ag - self.ast) + self.fy * self.ast

    @property
    def pn_max(self) -> float:
        return 0.80 * self.p0

    @property
    def phi_p0(self) -> float:
        return 0.65 * self.p0

    def _bar_coordinates(self) -> list[tuple[float, float, float, int]]:
        """Devuelve (x, y, As, número), con x/y medidos desde caras izquierda/superior."""
        rows, cols = self.rebar.matrix.shape
        barras: list[tuple[float, float, float, int]] = []
        for i in range(rows):
            for j in range(cols):
                n = int(self.rebar.matrix[i, j])
                if n == 0:
                    continue
                db = self.rebar.diametro(n)
                # Se reparte la grilla entre centros extremos, respetando recubrimiento y estribo.
                x0 = self.rec + self.tie_diameter + db / 2
                y0 = self.rec + self.tie_diameter + db / 2
                x = self.b / 2 if cols == 1 else x0 + j * (self.b - 2 * x0) / (cols - 1)
                y = self.h / 2 if rows == 1 else y0 + i * (self.h - 2 * y0) / (rows - 1)
                if not (0 < x < self.b and 0 < y < self.h):
                    raise GeometriaError("El recubrimiento, estribo y barras no caben en la sección.")
                barras.append((x, y, self.rebar.area(n), n))
        return barras

    def capas_acero(self, eje: str, cara: str = "superior") -> list[tuple[float, float]]:
        """Agrupa As por profundidad desde la cara comprimida, para X o Y."""
        if eje not in {"x", "y"} or cara not in {"superior", "inferior"}:
            raise ValueError("Eje o cara no válido.")
        depth = self.h if eje == "x" else self.b
        layers: dict[float, float] = {}
        for x, y, area, _ in self._bars:
            coordinate = y if eje == "x" else x
            d = coordinate if cara == "superior" else depth - coordinate
            key = round(d, 9)
            layers[key] = layers.get(key, 0.0) + area
        return sorted(layers.items())

    def factor_phi(self, et: float) -> float:
        ey = self.fy / self.vars.ES
        if et <= ey:
            return 0.65
        if et >= 0.005:
            return 0.90
        return 0.65 + 0.25 * (et - ey) / (0.005 - ey)

    def calcular_punto(self, c: float, eje: str = "x", cara: str = "superior") -> PuntoInteraccion:
        """Calcula Pn y Mn por compatibilidad para una profundidad c dada."""
        if c <= 0:
            raise ValueError("La profundidad del eje neutro debe ser positiva.")
        depth = self.h if eje == "x" else self.b
        width = self.b if eje == "x" else self.h
        sign = 1.0 if cara == "superior" else -1.0
        beta = self.vars.beta1(self.fc)
        a = min(beta * c, depth)
        cc = 0.85 * self.fc * a * width
        pn = cc
        mn = sign * cc * (depth / 2 - a / 2)
        strains: list[float] = []
        for yi, asi in self.capas_acero(eje, cara):
            strain = 0.003 * (c - yi) / c
            fs = self.vars.esfuerzo_acero(self.fy, strain)
            force = asi * fs
            pn += force
            mn += sign * force * (depth / 2 - yi)
            if strain < 0:
                strains.append(abs(strain))
        et = max(strains, default=0.0)
        return PuntoInteraccion(c, pn, mn, self.factor_phi(et), et)

    def _branch(self, eje: str, cara: str, n: int = 160) -> list[PuntoInteraccion]:
        depth = self.h if eje == "x" else self.b
        # Logaritmo cerca de c=0 para representar correctamente tracción y transición.
        c_values = np.geomspace(max(depth * 1e-5, 1e-7), depth * 50.0, n)
        return [self.calcular_punto(c, eje, cara) for c in c_values]

    @staticmethod
    def _cut_curve(points: list[PuntoInteraccion], limit: float, reduced: bool) -> list[tuple[float, float]]:
        """Recorta P-M e incorpora por interpolación el punto exacto del límite."""
        result: list[tuple[float, float]] = []
        previous: PuntoInteraccion | None = None
        for point in points:
            p = point.pp if reduced else point.pn
            moment = point.mp if reduced else point.mn
            if p <= limit:
                result.append((moment, p))
            elif previous is not None:
                p0 = previous.pp if reduced else previous.pn
                m0 = previous.mp if reduced else previous.mn
                if p0 < limit < p:
                    t = (limit - p0) / (p - p0)
                    result.append((m0 + t * (moment - m0), limit))
                break
            previous = point
        return result

    def curva_interaccion(self, eje: str, reduced: bool = True) -> dict[str, list[tuple[float, float]]]:
        """Dos ramas del diagrama, cada una limitada a la compresión máxima ACI."""
        limit = (0.65 * self.pn_max) if reduced else self.pn_max
        return {
            "superior": self._cut_curve(self._branch(eje, "superior"), limit, reduced),
            "inferior": self._cut_curve(self._branch(eje, "inferior"), limit, reduced),
        }

    def capacidad_por_excentricidad(self, eje: str, eccentricidad: float) -> float:
        """Intersección de la curva φP-φM con la recta M=eP (compresión positiva)."""
        if abs(eccentricidad) < 1e-9:
            return 0.65 * self.pn_max
        curves = self.curva_interaccion(eje, reduced=True)
        candidates: list[float] = []
        for branch in curves.values():
            for (m1, p1), (m2, p2) in zip(branch, branch[1:]):
                g1, g2 = m1 - eccentricidad * p1, m2 - eccentricidad * p2
                if g1 == 0 and p1 > 0:
                    candidates.append(p1)
                elif g1 * g2 < 0:
                    t = -g1 / (g2 - g1)
                    p = p1 + t * (p2 - p1)
                    if p > 0:
                        candidates.append(p)
        if not candidates:
            raise BreslerError("La excentricidad no intersecta una capacidad de compresión válida.")
        return max(candidates)

    def verificar_biaxial(self, pu: float, mux: float, muy: float) -> dict[str, float | bool]:
        if pu <= 0:
            raise BreslerError("Bresler requiere una carga axial de compresión Pu mayor que cero.")
        ex, ey = muy / pu, mux / pu
        pnx = self.capacidad_por_excentricidad("x", ey)
        pny = self.capacidad_por_excentricidad("y", ex)
        pp0 = self.phi_p0
        denominator = 1 / pnx + 1 / pny - 1 / pp0
        if denominator <= 0 or not np.isfinite(denominator):
            raise BreslerError("La combinación de excentricidades no permite una solución de Bresler.")
        padm = 1 / denominator
        return {"pnx": pnx, "pny": pny, "p0": pp0, "padm": padm,
                "ratio": pu / padm, "cumple": pu <= padm, "ex": ex, "ey": ey}

    def verificar_cortante(self, pu: float, vu: float, s: float) -> ResultadoCortante:
        if s <= 0:
            raise ValueError("La separación de estribos debe ser mayor que cero.")
        # Se toma d en la dirección H; el programa muestra este supuesto explícitamente.
        avg_db = float(np.mean([self.rebar.diametro(n) for *_, n in self._bars]))
        d = self.h - self.rec - self.tie_diameter - avg_db / 2
        if d <= 0:
            raise GeometriaError("El peralte efectivo para cortante es negativo.")
        # La expresión ACI se evalúa en kgf, cm y kgf/cm²; Ag debe convertirse a cm².
        axial_factor = 1 + pu / (140 * (self.ag / cm**2) * kgf)
        vc = 0.53 * axial_factor * np.sqrt(self.fc * cm**2 / kgf) * (self.b / cm) * (d / cm)
        av = 2 * self.rebar.area(self.estribo) / cm**2
        vs = av * (self.fy * cm**2 / kgf) * (d / cm) / (s / cm)
        phi_vn = 0.75 * (vc + vs) * kgf
        needed_vs = vu / 0.75 - vc * kgf
        s_required = None if needed_vs <= 0 else av * (self.fy * cm**2 / kgf) * (d / cm) * cm / needed_vs
        return ResultadoCortante(vc * kgf, vs * kgf, phi_vn, s_required)


def fmt_force(value: float) -> str:
    return f"{value / ton:,.2f} tonf"


def plot_diagram(section: DiagramaInteraccion, axis: str, pu: float, mu: float):
    fig, ax = plt.subplots(figsize=(7, 5))
    for name, style, label in (("nominal", "--", "Capacidad nominal"), ("reducida", "-", "Capacidad φ")):
        for branch in section.curva_interaccion(axis, reduced=name == "reducida").values():
            values = np.asarray(branch)
            if len(values):
                ax.plot(values[:, 0] / ton, values[:, 1] / ton, style,
                        color="#172554" if name == "nominal" else "#0ea5e9",
                        alpha=0.8, label=label)
                label = None
    pmax = 0.65 * section.pn_max / ton
    ax.axhline(pmax, color="#f97316", lw=1.3, label="φPn,max")
    ax.scatter([mu / ton], [pu / ton], color="#dc2626", s=42, zorder=4, label="Demanda")
    ax.axhline(0, color="black", lw=0.7)
    ax.axvline(0, color="black", lw=0.7)
    ax.grid(color="#cbd5e1", linestyle="--", linewidth=0.6)
    title_axis = "Mₓ" if axis == "x" else "Mᵧ"
    ax.set(title=f"Diagrama P–{title_axis}", xlabel=f"{title_axis} (tonf·m)", ylabel="P (tonf)")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(dict(zip(labels, handles)).values(), dict(zip(labels, handles)).keys(), fontsize=8)
    fig.tight_layout()
    return fig


def plot_section(section: DiagramaInteraccion):
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.add_patch(plt.Rectangle((0, 0), section.b / cm, section.h / cm,
                               fill=False, lw=2, color="#1e293b"))
    for x, y, _, n in section._bars:
        ax.add_patch(plt.Circle((x / cm, y / cm), section.rebar.diametro(n) / (2 * cm),
                                color="#b91c1c", alpha=0.9))
        ax.text(x / cm, y / cm, f"#{n}", color="white", ha="center", va="center", fontsize=7)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set(xlabel="B (cm)", ylabel="H (cm)", title="Distribución de acero longitudinal")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    return fig


def default_matrix() -> pd.DataFrame:
    return pd.DataFrame([[8, 8, 8, 8], [8, 0, 0, 8], [8, 0, 0, 8], [8, 8, 8, 8]],
                        columns=["C1", "C2", "C3", "C4"])


def _vercel_number(values: dict[str, list[str]], name: str, default: float) -> float:
    """Obtiene un campo numérico de la calculadora WSGI de Vercel."""
    raw = values.get(name, [str(default)])[0].strip().replace(",", ".")
    return float(raw)


def _vercel_matrix(values: dict[str, list[str]]) -> np.ndarray:
    raw = values.get("armado", ["8,8,8,8;8,0,0,8;8,0,0,8;8,8,8,8"])[0]
    rows = [[int(item.strip()) for item in row.split(",")] for row in raw.split(";") if row.strip()]
    if not rows or any(len(row) != len(rows[0]) for row in rows):
        raise GeometriaError("La matriz debe usar filas separadas por ; y columnas separadas por ,.")
    return np.asarray(rows, dtype=int)


def _vercel_page(values: dict[str, list[str]]) -> str:
    """Vista HTML compacta para Vercel; reutiliza el mismo núcleo de cálculo."""
    defaults = {"b": 40.0, "h": 40.0, "rec": 4.0, "tie": 3.0, "fc": 280.0, "fy": 4200.0,
                "pu": 300.0, "mux": 25.0, "muy": 15.0, "vu": 35.0, "s": 15.0}
    current = {name: values.get(name, [str(default)])[0] for name, default in defaults.items()}
    armed = values.get("armado", ["8,8,8,8;8,0,0,8;8,0,0,8;8,8,8,8"])[0]
    results = ""
    try:
        data = {name: _vercel_number(values, name, default) for name, default in defaults.items()}
        section = DiagramaInteraccion(data["b"] * cm, data["h"] * cm, data["rec"] * cm, int(data["tie"]),
                                      data["fc"] * kgf / cm**2, data["fy"] * kgf / cm**2, _vercel_matrix(values))
        biaxial = section.verificar_biaxial(data["pu"] * ton, data["mux"] * ton * m, data["muy"] * ton * m)
        shear = section.verificar_cortante(data["pu"] * ton, data["vu"] * ton, data["s"] * cm)
        biaxial_class = "ok" if biaxial["cumple"] else "fail"
        shear_class = "ok" if data["vu"] * ton <= shear.phi_vn else "fail"
        spacing = "No requerido por resistencia" if shear.s_requerido is None else f"{shear.s_requerido / cm:.2f} cm"
        results = f"""
        <section class=\"results\"><h2>Resultados</h2>
        <div class=\"grid\"><article><h3>Sección</h3><p>Ag: <b>{section.ag / cm**2:.2f} cm²</b><br>Ast: <b>{section.ast / cm**2:.2f} cm²</b><br>ρ: <b>{section.rho:.3%}</b></p></article>
        <article class=\"{biaxial_class}\"><h3>Biaxial · Bresler</h3><p>φPnx: <b>{biaxial['pnx'] / ton:.2f} tonf</b><br>φPny: <b>{biaxial['pny'] / ton:.2f} tonf</b><br>φP admisible: <b>{biaxial['padm'] / ton:.2f} tonf</b><br>D/C: <b>{biaxial['ratio']:.3f}</b><br><strong>{'CUMPLE' if biaxial['cumple'] else 'FALLA'}</strong></p></article>
        <article class=\"{shear_class}\"><h3>Cortante</h3><p>Vc: <b>{shear.vc / ton:.2f} tonf</b><br>Vs: <b>{shear.vs / ton:.2f} tonf</b><br>φVn: <b>{shear.phi_vn / ton:.2f} tonf</b><br>Separación máxima: <b>{spacing}</b><br><strong>{'CUMPLE' if data['vu'] * ton <= shear.phi_vn else 'FALLA'}</strong></p></article></div></section>"""
    except (ValueError, TypeError, GeometriaError, CuantiaError, BreslerError) as error:
        results = f"<p class=\"error\"><strong>Datos no válidos:</strong> {escape(str(error))}</p>"

    fields = "".join(
        f'<label>{escape(name.upper())}<input name="{escape(name)}" type="number" step="any" value="{escape(str(current[name]))}"></label>'
        for name in defaults
    )
    return f"""<!doctype html><html lang=\"es\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>Columnas ACI 318-19</title>
    <style>body{{font-family:system-ui,sans-serif;background:#f1f5f9;color:#0f172a;margin:0}}main{{max-width:1000px;margin:2rem auto;padding:0 1rem}}header{{background:#172554;color:white;padding:1.5rem;border-radius:12px}}form,.results{{background:white;padding:1.25rem;border-radius:12px;margin-top:1rem;box-shadow:0 2px 8px #0f172a18}}.fields,.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:1rem}}label{{font-size:.82rem;font-weight:700}}input,textarea{{box-sizing:border-box;width:100%;margin-top:.3rem;padding:.55rem;border:1px solid #94a3b8;border-radius:6px}}textarea{{min-height:4.5rem}}button{{margin-top:1rem;background:#0284c7;color:white;border:0;border-radius:6px;padding:.65rem 1rem;font-weight:700;cursor:pointer}}article{{border-left:4px solid #64748b;padding:0 .85rem}}.ok{{border-color:#16a34a}}.fail,.error{{border-color:#dc2626;color:#991b1b}}.error{{background:#fee2e2;padding:1rem;border-radius:8px}}small{{color:#475569}}</style></head>
    <body><main><header><h1>Columnas de concreto armado · ACI 318-19</h1><p>Calculadora serverless: flexocompresión biaxial por Bresler y cortante.</p></header>
    <form method=\"get\"><h2>Datos de entrada</h2><div class=\"fields\">{fields}</div><label>Armado (0 = vacío; filas con ; y columnas con ,)<textarea name=\"armado\">{escape(armed)}</textarea></label><button type=\"submit\">Calcular</button><p><small>Unidades: B, H, r y s en cm; f'c/fy en kgf/cm²; Pu/Vu en tonf; Mux/Muy en tonf·m. Para diagramas P-M y editor gráfico completo, ejecute esta misma aplicación en Streamlit Cloud.</small></p></form>{results}</main></body></html>"""


def app(environ, start_response):
    """Entrada WSGI para Vercel; Vercel requiere esta variable en el módulo."""
    if environ.get("PATH_INFO", "/") == "/health":
        body = b'{"status":"ok"}'
        start_response("200 OK", [("Content-Type", "application/json; charset=utf-8"),
                                   ("Content-Length", str(len(body)))])
        return [body]
    body = _vercel_page(parse_qs(environ.get("QUERY_STRING", ""), keep_blank_values=True)).encode("utf-8")
    start_response("200 OK", [("Content-Type", "text/html; charset=utf-8"),
                               ("Content-Length", str(len(body)))])
    return [body]


def build_sidebar():
    with st.sidebar:
        st.header("Datos de entrada")
        st.caption("Unidades: cm, kgf/cm², tonf y tonf·m.")
        st.subheader("Geometría")
        b = st.number_input("B", min_value=1.0, value=40.0, step=1.0)
        h = st.number_input("H", min_value=1.0, value=40.0, step=1.0)
        rec = st.number_input("Recubrimiento r", min_value=0.1, value=4.0, step=0.5)
        tie = st.selectbox("Estribo", [2, 3, 4, 5], index=1, format_func=lambda n: f"Barra #{n}")
        st.subheader("Materiales")
        fc = st.number_input("f'c", min_value=100.0, value=280.0, step=10.0)
        fy = st.number_input("fy", min_value=2000.0, value=4200.0, step=100.0)
        st.subheader("Acero longitudinal")
        st.caption("Use 0 para una celda sin barra. Barras permitidas: #2–#8 y #10.")
        matrix = st.data_editor(default_matrix(), num_rows="dynamic", use_container_width=True,
                                key="armado", column_config={c: st.column_config.NumberColumn(c, min_value=0, step=1)
                                for c in default_matrix().columns})
        st.subheader("Demandas")
        pu = st.number_input("Pu", value=300.0, step=5.0)
        mux = st.number_input("Mux", value=25.0, step=1.0)
        muy = st.number_input("Muy", value=15.0, step=1.0)
        vu = st.number_input("Vu", min_value=0.0, value=35.0, step=1.0)
        s = st.number_input("Separación de estribos s", min_value=1.0, value=15.0, step=1.0)
        return b, h, rec, tie, fc, fy, matrix, pu, mux, muy, vu, s


def main() -> None:
    st.set_page_config(page_title="Columnas ACI 318-19", page_icon="🏗️", layout="wide")
    st.title("Columnas de concreto armado · ACI 318-19")
    st.caption("Análisis por compatibilidad de deformaciones, flexocompresión uniaxial/biaxial y cortante.")
    values = build_sidebar()
    b, h, rec, tie, fc, fy, matrix, pu, mux, muy, vu, s = values
    try:
        raw_matrix = matrix.fillna(0).astype(int).to_numpy()
        section = DiagramaInteraccion(b * cm, h * cm, rec * cm, tie,
                                      fc * kgf / cm**2, fy * kgf / cm**2, raw_matrix)
        if section.rho > 0.08:
            raise CuantiaError(f"Cuantía {section.rho:.2%}: excede el máximo de 8%.")
        if section.rho < 0.01:
            st.warning(f"Cuantía longitudinal {section.rho:.2%}: menor al mínimo de referencia de 1% ACI.")
    except (ValueError, TypeError, GeometriaError, CuantiaError) as error:
        st.error(f"Datos no válidos: {error}")
        st.stop()

    pu_i, mux_i, muy_i, vu_i, s_i = pu * ton, mux * ton * m, muy * ton * m, vu * ton, s * cm
    tabs = st.tabs(["Diagramas 2D", "Biaxial · Bresler", "Cortante", "Geometría y resumen"])
    with tabs[0]:
        col1, col2 = st.columns(2)
        with col1:
            st.pyplot(plot_diagram(section, "x", pu_i, mux_i), use_container_width=True)
        with col2:
            st.pyplot(plot_diagram(section, "y", pu_i, muy_i), use_container_width=True)
        st.caption("La línea naranja corresponde a φ·0.80·P₀. Las curvas se trazan para ambas caras comprimidas.")
    with tabs[1]:
        try:
            biaxial = section.verificar_biaxial(pu_i, mux_i, muy_i)
            cols = st.columns(5)
            for col, label, key in zip(cols, ["φPnx", "φPny", "φP0", "φP admisible", "D/C"],
                                       ["pnx", "pny", "p0", "padm", "ratio"]):
                col.metric(label, f"{biaxial[key] / ton:,.2f} tonf" if key != "ratio" else f"{biaxial[key]:.3f}")
            if biaxial["cumple"]:
                st.success("CUMPLE — Pu no supera la capacidad biaxial por carga recíproca de Bresler.")
            else:
                st.error("FALLA — Pu supera la capacidad biaxial por carga recíproca de Bresler.")
            st.caption(f"Excentricidades usadas: eₓ = {biaxial['ex'] / cm:.2f} cm; eᵧ = {biaxial['ey'] / cm:.2f} cm.")
        except BreslerError as error:
            st.error(f"No fue posible comprobar Bresler: {error}")
    with tabs[2]:
        try:
            shear = section.verificar_cortante(pu_i, vu_i, s_i)
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Vc", fmt_force(shear.vc))
            c2.metric("Vs", fmt_force(shear.vs))
            c3.metric("φVn", fmt_force(shear.phi_vn))
            c4.metric("Vu", fmt_force(vu_i))
            if vu_i <= shear.phi_vn:
                st.success("CUMPLE — Vu ≤ φVn.")
            else:
                st.error("FALLA — Vu > φVn.")
            if shear.s_requerido is None:
                st.info("El concreto aporta la resistencia requerida; el espaciamiento queda sujeto a los máximos normativos de detallado.")
            else:
                st.write(f"Separación máxima por resistencia: **{shear.s_requerido / cm:.2f} cm** (verificar además límites de espaciamiento ACI).")
            st.caption("Se emplea d en la dirección H y Av = dos ramas del estribo. Verifique requisitos de confinamiento y cortante sísmico aplicables.")
        except (ValueError, GeometriaError) as error:
            st.error(f"No fue posible comprobar el cortante: {error}")
    with tabs[3]:
        c1, c2 = st.columns([1, 1])
        with c1:
            st.pyplot(plot_section(section), use_container_width=True)
        with c2:
            st.subheader("Propiedades de la sección")
            summary = pd.DataFrame({
                "Propiedad": ["Ag", "Ast", "ρ", "P0 nominal", "Pn,max nominal", "φPn,max"],
                "Valor": [f"{section.ag / cm**2:.2f} cm²", f"{section.ast / cm**2:.2f} cm²",
                          f"{section.rho:.3%}", fmt_force(section.p0), fmt_force(section.pn_max),
                          fmt_force(0.65 * section.pn_max)],
            })
            st.table(summary)
            st.subheader("Matriz de barras")
            st.dataframe(pd.DataFrame(raw_matrix), use_container_width=True, hide_index=True)
            st.info("Supuestos: columnas con estribos (φ=0.65 en compresión), εcu=0.003, bloque rectangular de Whitney y acero elastoplástico perfecto.")


if __name__ == "__main__":
    main()
