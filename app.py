import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from io import BytesIO
from base64 import b64encode
from html import escape
from urllib.parse import parse_qs

m = 1.0
cm = 0.01 * m
mm = 0.001 * m
kgf = 1.0
ton = 1000.0 * kgf


class Rebar:
    def __init__(self, matrix):
        self.matrix = np.array(matrix)

    def diametro(self, n):
        tabla = {2: 6.4*mm, 3: 9.5*mm, 4: 12.7*mm, 5: 15.9*mm,
                 6: 19.1*mm, 7: 22.2*mm, 8: 25.4*mm, 10: 31.8*mm}
        return tabla.get(int(n), 0.0)

    def area(self, n):
        tabla = {2: 0.32*(cm**2), 3: 0.71*(cm**2), 4: 1.27*(cm**2), 5: 1.98*(cm**2),
                 6: 2.85*(cm**2), 7: 3.88*(cm**2), 8: 5.10*(cm**2), 10: 7.92*(cm**2)}
        return tabla.get(int(n), 0.0)


class VariablesGenerales:
    def esfuerzo_acero(self, fy, es):
        Es = 2.0e6 * (kgf / (cm**2))
        ey = fy / Es
        if es > ey:
            return fy
        elif es < -ey:
            return -fy
        else:
            return Es * es

    def beta1(self, fc):
        if 175 <= fc <= 280:
            return 0.85
        elif fc >= 550:
            return 0.65
        else:
            return 0.85 - 0.05 * (fc - 280) / 70.0


class DiagramaInteraccion:
    def __init__(self, b, h, rec, tie_diameter_m, fc, fy, rebar_matrix):
        self.b = b
        self.h = h
        self.rec = rec
        self.fc = fc
        self.fy = fy
        self.tie_diameter = tie_diameter_m
        self.rebar = Rebar(rebar_matrix)
        self.vg = VariablesGenerales()
        self.centroide = h / 2.0
        self._validate()
        self._bars = self._bar_coordinates()

    def _validate(self):
        if self.b <= 0 or self.h <= 0:
            raise ValueError("B y H deben ser mayores que cero.")
        if 2 * self.rec >= self.b or 2 * self.rec >= self.h:
            raise ValueError("El recubrimiento deja sin nucleo util a la seccion.")
        if self.fc <= 0 or self.fy <= 0:
            raise ValueError("f'c y fy deben ser mayores que cero.")

    def _bar_coordinates(self):
        rows, cols = self.rebar.matrix.shape
        barras = []
        x0 = self.rec + self.tie_diameter
        y0 = self.rec + self.tie_diameter
        denom_x = max(cols - 1, 1)
        denom_y = max(rows - 1, 1)
        dx = (self.b - 2 * x0) / denom_x
        dy = (self.h - 2 * y0) / denom_y
        if dx < 0 or dy < 0:
            raise ValueError("Recubrimiento y estribo no caben en la seccion.")
        for i in range(rows):
            for j in range(cols):
                n = int(self.rebar.matrix[i, j])
                if n == 0:
                    continue
                x = self.b / 2 if cols == 1 else x0 + j * dx
                y = self.h / 2 if rows == 1 else y0 + i * dy
                if not (0 < x < self.b and 0 < y < self.h):
                    raise ValueError("Barra fuera de la seccion. Revise recubrimiento y estribo.")
                barras.append((x, y, self.rebar.area(n), n, self.rebar.diametro(n)))
        if not barras:
            raise ValueError("Debe existir al menos una barra longitudinal.")
        return barras

    def ag(self):
        return self.b * self.h

    def ast(self):
        return sum(area for *_, area, _, _ in self._bars)

    def rho(self):
        return self.ast() / self.ag()

    def p0(self):
        return 0.85 * self.fc * (self.ag() - self.ast()) + self.fy * self.ast()

    def pn_max(self):
        return 0.80 * self.p0()

    def capas_acero(self, eje, cara="superior"):
        depth = self.h if eje == "x" else self.b
        layers = {}
        for x, y, area, n, db in self._bars:
            coord = y if eje == "x" else x
            d = coord if cara == "superior" else depth - coord
            key = round(d, 9)
            layers[key] = layers.get(key, 0.0) + area
        return sorted(layers.items())

    def factor_phi(self, et):
        Es = 2.0e6 * (kgf / (cm**2))
        ey = self.fy / Es
        if et <= ey:
            return 0.65
        elif et >= 0.005:
            return 0.90
        else:
            return 0.65 + 0.25 * (et - ey) / (0.005 - ey)

    def calcular_punto(self, c, eje="x", cara="superior"):
        if c <= 0:
            return 0, 0, 0, 0
        depth = self.h if eje == "x" else self.b
        width = self.b if eje == "x" else self.h
        beta = self.vg.beta1(self.fc)
        a = min(beta * c, depth)
        cc = 0.85 * self.fc * a * width
        pn = cc
        mn = cc * (depth / 2 - a / 2)
        if cara == "inferior":
            mn = -mn
        strains = []
        for yi, asi in self.capas_acero(eje, cara):
            strain = 0.003 * (c - yi) / c
            fs = self.vg.esfuerzo_acero(self.fy, strain)
            force = asi * fs
            pn += force
            mn += force * (depth / 2 - yi)
            if cara == "inferior":
                mn = -mn + 2 * (depth / 2 - yi) * force
            if strain < 0:
                strains.append(abs(strain))
        et = max(strains) if strains else 0.0
        phi = self.factor_phi(et)
        return pn, mn, phi, et

    def curva_interaccion(self, eje="x", n=60):
        depth = self.h if eje == "x" else self.b
        c_values = np.linspace(0.01 * cm, 3.0 * depth, n)
        p_list, m_list, pp_list, mp_list = [], [], [], []
        for c in c_values:
            pn, mn, phi, et = self.calcular_punto(c, eje, "superior")
            p_list.append(pn)
            m_list.append(mn)
            pp_list.append(phi * pn)
            mp_list.append(phi * mn)
        return p_list, m_list, pp_list, mp_list, self.pn_max()

    @staticmethod
    def _clip_polygon(polygon, nx, ny, limit):
        clipped = []
        for first, second in zip(polygon, polygon[1:] + polygon[:1]):
            q1 = nx * first[0] + ny * first[1]
            q2 = nx * second[0] + ny * second[1]
            inside1 = q1 >= limit - 1e-12
            inside2 = q2 >= limit - 1e-12
            if inside1:
                clipped.append(first)
            if inside1 != inside2:
                t = (limit - q1) / (q2 - q1) if q2 != q1 else 0
                clipped.append((first[0] + t * (second[0] - first[0]),
                                first[1] + t * (second[1] - first[1])))
        return clipped

    @staticmethod
    def _polygon_area_centroid(polygon):
        cross = [x1 * y2 - x2 * y1 for (x1, y1), (x2, y2) in zip(polygon, polygon[1:] + polygon[:1])]
        twice = sum(cross)
        if abs(twice) < 1e-14:
            return 0.0, 0.0, 0.0
        cx = sum((x1 + x2) * v for (x1, y1), (x2, y2), v in zip(polygon, polygon[1:] + polygon[:1], cross)) / (3 * twice)
        cy = sum((y1 + y2) * v for (x1, y1), (x2, y2), v in zip(polygon, polygon[1:] + polygon[:1], cross)) / (3 * twice)
        return abs(twice) / 2, cx, cy

    def calcular_punto_biaxial(self, c, theta):
        nx = float(np.cos(theta))
        ny = float(np.sin(theta))
        corners = [(0.0, 0.0), (self.b, 0.0), (self.b, self.h), (0.0, self.h)]
        q_max = max(nx * x + ny * y for x, y in corners)
        q_min = min(nx * x + ny * y for x, y in corners)
        projection = q_max - q_min
        beta = self.vg.beta1(self.fc)
        a = min(beta * c, projection)
        block = self._clip_polygon(corners, nx, ny, q_max - a)
        area, x_c, y_c = self._polygon_area_centroid(block)
        cc = 0.85 * self.fc * area
        pn = cc
        mx = cc * (self.h / 2 - y_c)
        my = cc * (self.b / 2 - x_c)
        strains = []
        for x, y, as_i, n, db in self._bars:
            depth = q_max - (nx * x + ny * y)
            strain = 0.003 * (c - depth) / c
            fs = self.vg.esfuerzo_acero(self.fy, strain)
            force = as_i * fs
            pn += force
            mx += force * (self.h / 2 - y)
            my += force * (self.b / 2 - x)
            if strain < 0:
                strains.append(abs(strain))
        et = max(strains) if strains else 0.0
        phi = self.factor_phi(et)
        return phi * pn, phi * mx, phi * my

    def contorno_aci(self, pu, n_theta=48):
        if pu <= 0:
            raise ValueError("Pu debe ser mayor que cero (compresion).")
        phi_pmax = 0.65 * self.pn_max()
        if pu > phi_pmax:
            raise ValueError(f"Pu excede \u03c6Pn,max ({phi_pmax / ton:.2f} tonf).")
        contour = []
        for theta in np.linspace(0, 2 * np.pi, n_theta, endpoint=False):
            projection = abs(np.cos(theta)) * self.b + abs(np.sin(theta)) * self.h
            c_values = np.geomspace(max(projection * 1e-5, 1e-7), projection * 50, 120)
            points = [self.calcular_punto_biaxial(c, theta) for c in c_values]
            candidates = []
            for (p1, mx1, my1), (p2, mx2, my2) in zip(points, points[1:]):
                if (p1 - pu) * (p2 - pu) <= 0 and p1 != p2:
                    t = (pu - p1) / (p2 - p1)
                    candidates.append((mx1 + t * (mx2 - mx1), my1 + t * (my2 - my1)))
            if candidates:
                contour.append(max(candidates, key=lambda pt: np.hypot(*pt)))
        if len(contour) < 8:
            raise ValueError("No se pudo formar el contorno; revise los datos.")
        return np.asarray(contour)


def plot_diagram(section, eje, pu, mu):
    pn, mn, pp, mp, pmax = section.curva_interaccion(eje)
    fig, ax = plt.subplots(figsize=(7, 5))
    mn_neg = [-x for x in mn]
    mp_neg = [-x for x in mp]
    ax.plot([x / ton for x in mn], [x / ton for x in pn], 'k--', lw=1.5, label='Nominal')
    ax.plot([x / ton for x in mn_neg], [x / ton for x in pn], 'k--', lw=1.5)
    ax.plot([x / ton for x in mp], [x / ton for x in pp], 'b-', lw=2, label='\u03c6 (dise\u00f1o)')
    ax.plot([x / ton for x in mp_neg], [x / ton for x in pp], 'b-', lw=2)
    ax.axhline(0.65 * pmax / ton, color='orange', lw=1.3, label='\u03c6Pn,max')
    ax.axhline(0, color='black', lw=0.7)
    ax.axvline(0, color='black', lw=0.7)
    ax.scatter([mu / ton], [pu / ton], color='red', s=50, zorder=5, label='Demanda')
    ax.grid(color='silver', linestyle='--', linewidth=0.5)
    title = 'Mx' if eje == 'x' else 'My'
    ax.set(xlabel=f'{title} (tonf\u00b7m)', ylabel='P (tonf)', title=f'Diagrama P\u2013{title}')
    ax.legend(fontsize=8, loc='lower right')
    fig.tight_layout()
    return fig


def plot_contorno_aci(contour, mux, muy):
    closed = np.vstack((contour, contour[:1]))
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.fill(closed[:, 0] / ton, closed[:, 1] / ton, color='#bae6fd', alpha=0.55, label='Region resistente')
    ax.plot(closed[:, 0] / ton, closed[:, 1] / ton, color='#0284c7', lw=2, label='Contorno ACI \u03c6P = Pu')
    ax.scatter([mux / ton], [muy / ton], color='red', s=48, zorder=5, label='Demanda')
    ax.axhline(0, color='black', lw=0.7)
    ax.axvline(0, color='black', lw=0.7)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(color='silver', linestyle='--', linewidth=0.5)
    ax.set(xlabel='Mx (tonf\u00b7m)', ylabel='My (tonf\u00b7m)', title='Contorno biaxial \u03c6Mx vs \u03c6My a Pu dado')
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def plot_section(section):
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.add_patch(plt.Rectangle((0, 0), section.b / cm, section.h / cm,
                                fill=False, lw=2.5, color='#1e293b'))
    for x, y, area, n, db in section._bars:
        r = db / (2 * cm)
        ax.add_patch(plt.Circle((x / cm, y / cm), r, color='#b91c1c', alpha=0.85))
        ax.text(x / cm, y / cm, f'#{int(n)}', color='white',
                ha='center', va='center', fontsize=9, fontweight='bold')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set(xlabel='B (cm)', ylabel='H (cm)', title='Distribucion de acero longitudinal')
    ax.grid(alpha=0.2)
    fig.tight_layout()
    return fig


def fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=130, bbox_inches='tight')
    plt.close(fig)
    return b64encode(buf.getvalue()).decode('ascii')


def _parse_matrix(raw):
    rows = [[int(x.strip()) for x in row.split(',')] for row in raw.split(';') if row.strip()]
    if not rows or any(len(r) != len(rows[0]) for r in rows):
        raise ValueError("Matriz invalida: use filas separadas por ; y columnas por ,")
    return np.array(rows, dtype=int)


def app(environ, start_response):
    qs = parse_qs(environ.get('QUERY_STRING', ''), keep_blank_values=True)
    if environ.get('PATH_INFO', '/') == '/health':
        body = b'{"status":"ok"}'
        start_response('200 OK', [('Content-Type', 'application/json; charset=utf-8')])
        return [body]
    html = ''
    try:
        b = float(qs.get('b', ['40'])[0])
        h = float(qs.get('h', ['40'])[0])
        rec = float(qs.get('rec', ['4'])[0])
        tie_mm = float(qs.get('tie', ['9.5'])[0])
        fc = float(qs.get('fc', ['280'])[0])
        fy = float(qs.get('fy', ['4200'])[0])
        pu = float(qs.get('pu', ['300'])[0])
        mux = float(qs.get('mux', ['25'])[0])
        muy = float(qs.get('muy', ['15'])[0])
        raw_mat = qs.get('armado', ['8,8,8,8;8,0,0,8;8,0,0,8;8,8,8,8'])[0]

        mat = _parse_matrix(raw_mat)
        sec = DiagramaInteraccion(b * cm, h * cm, rec * cm, tie_mm * mm,
                                  fc * kgf / (cm**2), fy * kgf / (cm**2), mat)

        fig1 = plot_diagram(sec, 'x', pu * ton, mux * ton * m)
        fig2 = plot_diagram(sec, 'y', pu * ton, muy * ton * m)
        fig3 = plot_section(sec)
        img1 = fig_to_b64(fig1)
        img2 = fig_to_b64(fig2)
        img3 = fig_to_b64(fig3)

        biaxial_html = ''
        if abs(muy) > 0.01:
            try:
                contour = sec.contorno_aci(pu * ton)
                fig4 = plot_contorno_aci(contour, mux * ton * m, muy * ton * m)
                img4 = fig_to_b64(fig4)
                mx_uni = float(np.max(np.abs(contour[:, 0]))) / ton
                my_uni = float(np.max(np.abs(contour[:, 1]))) / ton
                radial = np.hypot(mux, muy)
                cap_radial = np.hypot(mx_uni, my_uni)
                ratio = radial / cap_radial if cap_radial > 0 else 999
                status = 'CUMPLE' if ratio <= 1 else 'FALLA'
                biaxial_html = f'''
                <h2>Verificacion biaxial ACI</h2>
                <div class="grid">
                <div><img src="data:image/png;base64,{img4}" style="max-width:100%"></div>
                <div>
                <p><b>\u03c6Mnx (uniaxial):</b> {mx_uni:,.3f} tonf\u00b7m</p>
                <p><b>\u03c6Mny (uniaxial):</b> {my_uni:,.3f} tonf\u00b7m</p>
                <p><b>D/C radial:</b> {ratio:.3f}</p>
                <p class="{'ok' if ratio<=1 else 'fail'}"><b>{status}</b></p>
                <p><small>Contorno ACI 318-19 para la demanda Pu dada. La demanda debe caer dentro de la region azul.</small></p>
                </div></div>
                '''
            except Exception as e:
                biaxial_html = f'<p class="error">Biaxial: {escape(str(e))}</p>'

        ag = sec.ag() / cm**2
        ast = sec.ast() / cm**2
        rho = sec.rho()
        p0 = sec.p0() / ton
        pnmax = sec.pn_max() / ton
        phipmax = 0.65 * sec.pn_max() / ton

        html = f'''<!doctype html>
<html lang="es"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Columnas ACI 318-19</title>
<style>
body{{font-family:sans-serif;background:#f1f5f9;color:#0f172a;margin:0;padding:1rem}}
header{{background:#172554;color:#fff;padding:1rem;border-radius:10px;margin-bottom:1rem}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:1rem}}
.card{{background:#fff;padding:1rem;border-radius:10px;box-shadow:0 1px 4px #00000020}}
img{{max-width:100%;height:auto}}
table{{width:100%;border-collapse:collapse}}
td,th{{padding:.4rem;border:1px solid #cbd5e1;text-align:left}}
.ok{{color:#16a34a;font-weight:bold}}
.fail{{color:#dc2626;font-weight:bold}}
.error{{color:#dc2626;background:#fee2e2;padding:.5rem;border-radius:6px}}
form{{display:flex;flex-wrap:wrap;gap:.5rem;align-items:end;margin-bottom:1rem}}
form label{{font-size:.85rem}}
form input{{padding:.4rem;border:1px solid #94a3b8;border-radius:5px;width:90px}}
form button{{background:#0284c7;color:#fff;border:0;padding:.5rem 1rem;border-radius:5px;cursor:pointer}}
</style></head>
<body>
<header><h1>Columnas de concreto armado \u00b7 ACI 318-19</h1>
<p>Diagramas P\u2013M por compatibilidad de deformaciones y verificacion biaxial ACI.</p></header>
<form method="get">
<label>B (cm) <input name="b" value="{b}"></label>
<label>H (cm) <input name="h" value="{h}"></label>
<label>r (cm) <input name="rec" value="{rec}"></label>
<label>Estribo (mm) <input name="tie" value="{tie_mm}"></label>
<label>f\'c (kgf/cm\u00b2) <input name="fc" value="{fc}"></label>
<label>fy (kgf/cm\u00b2) <input name="fy" value="{fy}"></label>
<label>Pu (tonf) <input name="pu" value="{pu}"></label>
<label>Mux (tonf\u00b7m) <input name="mux" value="{mux}"></label>
<label>Muy (tonf\u00b7m) <input name="muy" value="{muy}"></label>
<label>Armado (; filas, , columnas)<br><textarea name="armado" rows="3" style="width:200px">{escape(raw_mat)}</textarea></label>
<button>Calcular</button>
</form>
<div class="grid">
<div class="card"><img src="data:image/png;base64,{img1}" style="width:100%"></div>
<div class="card"><img src="data:image/png;base64,{img2}" style="width:100%"></div>
</div>
{biaxial_html}
<div class="grid">
<div class="card"><img src="data:image/png;base64,{img3}" style="width:100%"></div>
<div class="card"><h2>Propiedades de la seccion</h2>
<table><tr><th>Propiedad</th><th>Valor</th><th>Descripcion</th></tr>
<tr><td>Ag</td><td>{ag:.2f} cm\u00b2</td><td>Area bruta de concreto</td></tr>
<tr><td>Ast</td><td>{ast:.2f} cm\u00b2</td><td>Area total de acero longitudinal</td></tr>
<tr><td>\u03c1</td><td>{rho:.3%}</td><td>Cuantia de acero (Ast/Ag)</td></tr>
<tr><td>P0 nominal</td><td>{p0:,.2f} tonf</td><td>Carga axial nominal en compresion pura</td></tr>
<tr><td>Pn,max nominal</td><td>{pnmax:,.2f} tonf</td><td>Max. nominal = 0.80\u00b7P0 (ACI)</td></tr>
<tr><td>\u03c6Pn,max</td><td>{phipmax:,.2f} tonf</td><td>Capacidad de dise\u00f1o = 0.65\u00b7Pn,max</td></tr>
</table>
<p><small>\u03b5cu=0.003, bloque de Whitney, acero elastoplastico. \u03c6=0.65 en compresion controlada.</small></p>
</div></div>
</body></html>'''
    except Exception as e:
        html = f'<p class="error">Error: {escape(str(e))}</p>'
    body = html.encode('utf-8')
    start_response('200 OK', [('Content-Type', 'text/html; charset=utf-8')])
    return [body]


if __name__ == "__main__":
    import streamlit as st

    st.set_page_config(page_title='Columnas ACI 318-19', page_icon='\U0001f3d7\ufe0f', layout='wide')
    st.title('Columnas de concreto armado \u00b7 ACI 318-19')
    st.caption(
        'Diagramas P\u2013M por compatibilidad de deformaciones (ACI 318-19). '
        'Si ingresa Mux y Muy ambos distintos de cero, se activa la verificaci\u00f3n biaxial.'
    )

    with st.sidebar:
        st.header('Datos de entrada')

        st.subheader('Geometr\u00eda')
        b = st.number_input('B (cm)', min_value=10.0, value=40.0, step=1.0,
                            help='Ancho de la seccion transversal')
        h = st.number_input('H (cm)', min_value=10.0, value=40.0, step=1.0,
                            help='Peralte de la seccion transversal')
        rec = st.number_input('Recubrimiento r (cm)', min_value=1.0, value=4.0, step=0.5,
                              help='Distancia desde la cara del concreto al borde exterior del estribo')
        tie_mm = st.number_input(
            'Estribo (mm)', min_value=4.0, max_value=25.0, value=9.5, step=0.5,
            help='Diametro del estribo en mm. Ej: #3 = 9.5 mm, #4 = 12.7 mm'
        )

        st.subheader('Materiales')
        fc = st.number_input("f'c (kgf/cm\u00b2)", min_value=100.0, value=280.0, step=10.0,
                             help='Resistencia especifica del concreto')
        fy = st.number_input('fy (kgf/cm\u00b2)', min_value=2000.0, value=4200.0, step=100.0,
                             help='Resistencia a la fluencia del acero')

        st.subheader('Acero longitudinal')
        st.caption('Use 0 para celdas sin barra. Barras disponibles: #2 (6.4 mm) a #10 (31.8 mm).')
        matrix = st.data_editor(
            pd.DataFrame([[8, 8, 8, 8], [8, 0, 0, 8], [8, 0, 0, 8], [8, 8, 8, 8]],
                         columns=['C1', 'C2', 'C3', 'C4']),
            num_rows='dynamic', use_container_width=True, key='armado',
            column_config={c: st.column_config.NumberColumn(c, min_value=0, step=1)
                           for c in ['C1', 'C2', 'C3', 'C4']})

        st.subheader('Demandas de dise\u00f1o')
        pu = st.number_input('Pu (tonf)', min_value=0.0, value=300.0, step=5.0,
                             help='Carga axial ultima de diseno (compresion positiva)')
        mux = st.number_input('Mux (tonf\u00b7m)', min_value=0.0, value=25.0, step=1.0,
                              help='Momento flector ultimo alrededor del eje X')
        muy = st.number_input('Muy (tonf\u00b7m)', min_value=0.0, value=15.0, step=1.0,
                              help='Momento flector ultimo alrededor del eje Y')

    try:
        raw_matrix = matrix.fillna(0).astype(int).to_numpy()
        section = DiagramaInteraccion(b * cm, h * cm, rec * cm, tie_mm * mm,
                                      fc * kgf / (cm**2), fy * kgf / (cm**2), raw_matrix)
        rho_val = section.rho()
        if rho_val > 0.08:
            st.error(f'Cuant\u00eda {rho_val:.2%}: excede el m\u00e1ximo de 8% (ACI 318-19 \u00a710.3).')
            st.stop()
        if rho_val < 0.01:
            st.warning(f'Cuant\u00eda longitudinal {rho_val:.2%}: menor al 1% m\u00ednimo (ACI 318-19 \u00a710.3).')
    except Exception as e:
        st.error(f'Error en datos: {e}')
        st.stop()

    pu_i = pu * ton
    mux_i = mux * ton * m
    muy_i = muy * ton * m

    st.subheader('Diagramas de interacci\u00f3n P\u2013M')
    st.caption(
        'Cada gr\u00e1fico muestra la capacidad de la columna para flexi\u00f3n en un solo eje. '
        'La curva negra discontinua es la capacidad nominal (P\u2099, M\u2099). '
        'La curva azul es la capacidad de dise\u00f1o reducida por \u03c6. '
        'La l\u00ednea naranja (\u03c6Pn,max) es el l\u00edmite m\u00e1ximo de compresi\u00f3n seg\u00fan ACI 318-19 '
        '(0.80\u00b7\u03c6\u00b7P\u2080 para columnas con estribos); toda demanda debe estar por debajo de esta l\u00ednea. '
        'El punto rojo es la demanda (Pu, Mu) ingresada.'
    )
    col1, col2 = st.columns(2)
    with col1:
        st.pyplot(plot_diagram(section, 'x', pu_i, mux_i), use_container_width=True)
    with col2:
        st.pyplot(plot_diagram(section, 'y', pu_i, muy_i), use_container_width=True)

    es_biaxial = abs(muy) > 0.01
    if es_biaxial:
        st.subheader('Verificaci\u00f3n biaxial ACI 318-19')
        st.caption(
            'Cuando hay momentos en ambos ejes (M\u2093 y M\u1d35), el contorno biaxial '
            'verifica que la combinaci\u00f3n (\u03c6M\u2093, \u03c6M\u1d35) a carga axial Pu est\u00e9 dentro de la '
            'capacidad de la secci\u00f3n. Los gr\u00e1ficos P\u2013M\u2093 y P\u2013M\u1d35 por s\u00ed solos no capturan '
            'la interacci\u00f3n biaxial. Este contorno es el m\u00e9todo exigido por ACI 318-19.'
        )
        try:
            contour = section.contorno_aci(pu_i)
            col1, col2 = st.columns([3, 2])
            with col1:
                st.pyplot(plot_contorno_aci(contour, mux_i, muy_i), use_container_width=True)
            with col2:
                mx_uni = float(np.max(np.abs(contour[:, 0]))) if len(contour) else 0
                my_uni = float(np.max(np.abs(contour[:, 1]))) if len(contour) else 0
                st.metric('\u03c6Mnx (uniaxial en Pu)', f'{mx_uni / ton:,.3f} tonf\u00b7m')
                st.metric('\u03c6Mny (uniaxial en Pu)', f'{my_uni / ton:,.3f} tonf\u00b7m')
                radial = np.hypot(mux_i, muy_i)
                cap_radial = np.hypot(mx_uni, my_uni)
                ratio = radial / cap_radial if cap_radial > 0 else 999
                st.metric('Relaci\u00f3n D/C radial', f'{ratio:.3f}')
                if ratio <= 1:
                    st.success('CUMPLE \u2014 La demanda cae dentro del contorno ACI.')
                else:
                    st.error('FALLA \u2014 La demanda excede el contorno ACI.')
                st.caption('D/C radial: distancia desde el origen hasta la demanda dividida por la distancia hasta el contorno en la misma direcci\u00f3n.')
        except Exception as e:
            st.error(f'No se pudo calcular el contorno biaxial: {e}')

    st.subheader('Geometr\u00eda de la secci\u00f3n')
    st.caption('Distribuci\u00f3n de las barras longitudinales en la secci\u00f3n transversal.')
    c1, c2 = st.columns([1, 1])
    with c1:
        st.pyplot(plot_section(section), use_container_width=True)
    with c2:
        desc = [
            '\u00c1rea bruta de la secci\u00f3n de concreto',
            '\u00c1rea total de acero longitudinal',
            'Cuant\u00eda de refuerzo = Ast / Ag',
            'Capacidad nominal en compresi\u00f3n pura sin excentricidad',
            'Capacidad nominal m\u00e1xima = 0.80\u00b7P\u2080 (ACI 318-19 \u00a722.4.2)',
            'Capacidad de dise\u00f1o = \u03c6\u00b7Pn,max, con \u03c6 = 0.65',
        ]
        data = {
            'Propiedad': ['Ag', 'Ast', '\u03c1', 'P\u2080 nominal', 'Pn,max nominal', '\u03c6Pn,max'],
            'Valor': [
                f'{section.ag() / cm**2:.2f} cm\u00b2',
                f'{section.ast() / cm**2:.2f} cm\u00b2',
                f'{section.rho():.3%}',
                f'{section.p0() / ton:,.2f} tonf',
                f'{section.pn_max() / ton:,.2f} tonf',
                f'{0.65 * section.pn_max() / ton:,.2f} tonf',
            ],
            'Descripci\u00f3n': desc,
        }
        st.table(pd.DataFrame(data))
        st.subheader('Matriz de barras')
        st.dataframe(pd.DataFrame(raw_matrix), use_container_width=True, hide_index=True)

    st.caption(
        'Supuestos: \u03b5cu = 0.003, bloque rectangular de Whitney (\u03b2\u2081 seg\u00fan ACI), '
        'acero elastopl\u00e1stico perfecto, columnas con estribos (\u03c6 = 0.65 en '
        'compresi\u00f3n controlada, \u03c6 = 0.90 en tracci\u00f3n controlada, transici\u00f3n lineal).'
    )
