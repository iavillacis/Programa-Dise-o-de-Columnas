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
    def __init__(self, b, h, rec, tie_diameter_m, fc, fy, diam_mm_matrix):
        self.b = b
        self.h = h
        self.rec = rec
        self.fc = fc
        self.fy = fy
        self.tie_diameter = tie_diameter_m
        self.vg = VariablesGenerales()
        self.centroide = h / 2.0
        self._validate()
        self._bars = self._bar_coordinates(diam_mm_matrix)

    def _validate(self):
        if self.b <= 0 or self.h <= 0:
            raise ValueError("B y H deben ser mayores que cero.")
        if 2 * self.rec >= self.b or 2 * self.rec >= self.h:
            raise ValueError("El recubrimiento deja sin nucleo util a la seccion.")
        if self.fc <= 0 or self.fy <= 0:
            raise ValueError("f'c y fy deben ser mayores que cero.")

    def _bar_coordinates(self, diam_mm_matrix):
        rows, cols = diam_mm_matrix.shape
        barras = []
        x0 = self.rec + self.tie_diameter
        y0 = self.rec + self.tie_diameter
        denom_x = max(cols - 1, 1)
        denom_y = max(rows - 1, 1)
        dx = (self.b - 2 * x0) / denom_x
        dy = (self.h - 2 * y0) / denom_y
        if dx < -1e-12 or dy < -1e-12:
            raise ValueError("Recubrimiento y estribo no caben en la seccion.")
        dx = max(dx, 0)
        dy = max(dy, 0)
        for i in range(rows):
            for j in range(cols):
                diam_mm_val = float(diam_mm_matrix[i, j])
                if diam_mm_val <= 0:
                    continue
                diam_m = diam_mm_val * mm
                area = np.pi * (diam_m) ** 2 / 4
                x = self.b / 2 if cols == 1 else x0 + j * dx
                y = self.h / 2 if rows == 1 else y0 + i * dy
                if not (0 < x < self.b and 0 < y < self.h):
                    raise ValueError("Barra fuera de la seccion. Revise recubrimiento y estribo.")
                barras.append((x, y, area, diam_mm_val, diam_m))
        if not barras:
            raise ValueError("Debe existir al menos una barra longitudinal.")
        return barras

    def verificar_espaciamiento(self):
        mensajes = []
        rows_dict = {}
        for x, y, area, dm, d in self._bars:
            yk = round(y, 9)
            rows_dict.setdefault(yk, []).append((x, d))
        for yk, bars in rows_dict.items():
            bars.sort(key=lambda b: b[0])
            for (x1, d1), (x2, d2) in zip(bars, bars[1:]):
                libre = (x2 - x1) - (d1 + d2) / 2
                if libre < 2.5 * cm:
                    mensajes.append(
                        f"Esp. libre horiz. = {libre / cm:.2f} cm en fila y={yk / cm:.1f} cm; minimo ACI = 2.50 cm"
                    )
                    break
        cols_dict = {}
        for x, y, area, dm, d in self._bars:
            xk = round(x, 9)
            cols_dict.setdefault(xk, []).append((y, d))
        for xk, bars in cols_dict.items():
            bars.sort(key=lambda b: b[0])
            for (y1, d1), (y2, d2) in zip(bars, bars[1:]):
                libre = (y2 - y1) - (d1 + d2) / 2
                if libre < 2.5 * cm:
                    mensajes.append(
                        f"Esp. libre vertical = {libre / cm:.2f} cm en columna x={xk / cm:.1f} cm; minimo ACI = 2.50 cm"
                    )
                    break
        return mensajes

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
        for x, y, area, dm, d in self._bars:
            coord = y if eje == "x" else x
            di = coord if cara == "superior" else depth - coord
            key = round(di, 9)
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
        for x, y, asi, dm, d in self._bars:
            depth = q_max - (nx * x + ny * y)
            strain = 0.003 * (c - depth) / c
            fs = self.vg.esfuerzo_acero(self.fy, strain)
            force = asi * fs
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
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.add_patch(plt.Rectangle((0, 0), section.b / cm, section.h / cm,
                                fill=False, lw=3, color='#1e293b'))
    for x, y, area, dm, d in section._bars:
        r = d / (2 * cm)
        ax.add_patch(plt.Circle((x / cm, y / cm), r, color='#b91c1c', alpha=0.85))
        ax.text(x / cm, y / cm, f'{dm:.1f}', color='white',
                ha='center', va='center', fontsize=8, fontweight='bold')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set(xlabel='B (cm)', ylabel='H (cm)', title='Distribucion de acero longitudinal (diametros en mm)')
    ax.grid(alpha=0.2)
    fig.tight_layout()
    return fig


def fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=130, bbox_inches='tight')
    plt.close(fig)
    return b64encode(buf.getvalue()).decode('ascii')


def _parse_matrix(raw):
    rows = [[float(x.strip()) for x in row.split(',')] for row in raw.split(';') if row.strip()]
    if not rows or any(len(r) != len(rows[0]) for r in rows):
        raise ValueError("Matriz invalida: use filas separadas por ; y columnas por ,")
    return np.array(rows, dtype=float)


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
        raw_mat = qs.get('armado', ['25.4,25.4,25.4,25.4;25.4,0,0,25.4;25.4,0,0,25.4;25.4,25.4,25.4,25.4'])[0]

        mat = _parse_matrix(raw_mat)
        sec = DiagramaInteraccion(b * cm, h * cm, rec * cm, tie_mm * mm,
                                  fc * kgf / (cm**2), fy * kgf / (cm**2), mat)

        esp_msgs = sec.verificar_espaciamiento()
        esp_html = ''
        if esp_msgs:
            esp_html = '<div class="error">' + '<br>'.join(escape(m) for m in esp_msgs) + '</div>'

        fig1 = plot_diagram(sec, 'x', pu * ton, mux * ton * m)
        fig2 = plot_diagram(sec, 'y', pu * ton, muy * ton * m)
        fig3 = plot_section(sec)
        img1 = fig_to_b64(fig1)
        img2 = fig_to_b64(fig2)
        img3 = fig_to_b64(fig3)

        biaxial_html = ''
        if abs(muy) > 0.01 and not esp_msgs:
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
                <p><b>\u03c6Mnx:</b> {mx_uni:,.3f} tonf\u00b7m</p>
                <p><b>\u03c6Mny:</b> {my_uni:,.3f} tonf\u00b7m</p>
                <p><b>D/C:</b> {ratio:.3f} <b>{status}</b></p>
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
form textarea{{padding:.4rem;border:1px solid #94a3b8;border-radius:5px}}
form button{{background:#0284c7;color:#fff;border:0;padding:.5rem 1rem;border-radius:5px;cursor:pointer}}
</style></head>
<body>
<header><h1>Columnas de concreto armado \u00b7 ACI 318-19</h1></header>
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
<label>Armado (mm; filas con ;, columnas con ,)<br><textarea name="armado" rows="4" style="width:260px">{escape(raw_mat)}</textarea></label>
<button style="padding:.6rem 2rem;font-size:1.1rem">CALCULAR</button>
</form>
{esp_html}
<div class="grid">
<div class="card"><img src="data:image/png;base64,{img1}" style="width:100%"></div>
<div class="card"><img src="data:image/png;base64,{img2}" style="width:100%"></div>
</div>
{biaxial_html}
<div class="grid">
<div class="card"><img src="data:image/png;base64,{img3}" style="width:100%"></div>
<div class="card"><h2>Propiedades de la seccion</h2>
<table><tr><th>Propiedad</th><th>Valor</th><th>Que es</th></tr>
<tr><td>Ag</td><td>{ag:.2f} cm\u00b2</td><td>Area bruta de concreto</td></tr>
<tr><td>Ast</td><td>{ast:.2f} cm\u00b2</td><td>Area total de acero longitudinal</td></tr>
<tr><td>\u03c1</td><td>{rho:.3%}</td><td>Cuantia = Ast/Ag</td></tr>
<tr><td>P0 nominal</td><td>{p0:,.2f} tonf</td><td>Capacidad en compresion pura</td></tr>
<tr><td>Pn,max nominal</td><td>{pnmax:,.2f} tonf</td><td>Max. nominal = 0.80\u00b7P0</td></tr>
<tr><td>\u03c6Pn,max</td><td>{phipmax:,.2f} tonf</td><td>Capacidad de dise\u00f1o = 0.65\u00b7Pn,max</td></tr>
</table>
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

    if "calcular" not in st.session_state:
        st.session_state.calcular = False

    with st.sidebar:
        st.header('Datos de entrada')
        st.subheader('Geometr\u00eda')
        b = st.number_input('B (cm)', min_value=10.0, value=40.0, step=1.0)
        h = st.number_input('H (cm)', min_value=10.0, value=40.0, step=1.0)
        rec = st.number_input('Rec. r (cm)', min_value=1.0, value=4.0, step=0.5)
        tie_mm = st.number_input('Estribo (mm)', min_value=4.0, max_value=25.0, value=9.5, step=0.5,
                                 help='Diametro del estribo. Referencia: #3=9.5mm, #4=12.7mm, #5=15.9mm')

        st.subheader('Materiales')
        fc = st.number_input("f'c (kgf/cm\u00b2)", min_value=100.0, value=280.0, step=10.0)
        fy = st.number_input('fy (kgf/cm\u00b2)', min_value=2000.0, value=4200.0, step=100.0)

        st.subheader('Armado longitudinal')
        st.caption('Ingrese el diametro en mm. Use 0 para celdas sin barra.')
        st.caption('Ej: 25.4 = barra #8, 12.7 = barra #4.')
        default_mat = pd.DataFrame(
            [[25.4, 25.4, 25.4, 25.4],
             [25.4, 0, 0, 25.4],
             [25.4, 0, 0, 25.4],
             [25.4, 25.4, 25.4, 25.4]],
            columns=['C1', 'C2', 'C3', 'C4']
        )
        matrix = st.data_editor(
            default_mat, num_rows='dynamic', height=320, use_container_width=True,
            key='armado',
            column_config={c: st.column_config.NumberColumn(c, min_value=0, step=1, format="%.1f")
                           for c in ['C1', 'C2', 'C3', 'C4']}
        )

        st.subheader('Demandas')
        pu = st.number_input('Pu (tonf)', min_value=0.0, value=300.0, step=5.0)
        mux = st.number_input('Mux (tonf\u00b7m)', min_value=0.0, value=25.0, step=1.0)
        muy = st.number_input('Muy (tonf\u00b7m)', min_value=0.0, value=15.0, step=1.0)

    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        if st.button('CALCULAR', type='primary', use_container_width=True):
            st.session_state.calcular = True
            st.rerun()

    if st.session_state.calcular:
        try:
            raw_matrix = matrix.fillna(0).astype(float).to_numpy()
            section = DiagramaInteraccion(
                b * cm, h * cm, rec * cm, tie_mm * mm,
                fc * kgf / (cm**2), fy * kgf / (cm**2), raw_matrix
            )

            esp_msgs = section.verificar_espaciamiento()
            for msg in esp_msgs:
                st.error(msg)

            rho_val = section.rho()
            if rho_val > 0.08:
                st.error(f'Cuantia {rho_val:.2%}: excede el maximo de 8% (ACI 318-19).')
                st.stop()
            if rho_val < 0.01:
                st.warning(f'Cuantia {rho_val:.2%}: menor al 1% minimo (ACI 318-19).')

            if esp_msgs:
                st.stop()
        except Exception as e:
            st.error(f'Error en datos: {e}')
            st.stop()

        pu_i = pu * ton
        mux_i = mux * ton * m
        muy_i = muy * ton * m

        st.subheader('Diagramas de interaccion P\u2013M')
        st.caption(
            'Curva nominal (negro discontinuo) y reducida por \u03c6 (azul). '
            'Linea naranja = \u03c6Pn,max (limite ACI 318-19). '
            'Punto rojo = demanda (Pu, Mu).'
        )
        col1, col2 = st.columns(2)
        with col1:
            st.pyplot(plot_diagram(section, 'x', pu_i, mux_i), use_container_width=True)
        with col2:
            st.pyplot(plot_diagram(section, 'y', pu_i, muy_i), use_container_width=True)

        if abs(muy) > 0.01:
            st.subheader('Verificacion biaxial ACI 318-19')
            st.caption(
                'Cuando hay momentos en ambos ejes, los diagramas P\u2013Mx y P\u2013My '
                'no son suficientes. El contorno biaxial verifica la interaccion '
                'combinada (Mx, My) a carga axial Pu.'
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
                    st.metric('D/C radial', f'{ratio:.3f}')
                    if ratio <= 1:
                        st.success('CUMPLE \u2014 La demanda esta dentro del contorno ACI.')
                    else:
                        st.error('FALLA \u2014 La demanda excede el contorno ACI.')
                    st.caption('D/C = distancia demanda / distancia contorno en misma direccion.')
            except Exception as e:
                st.error(f'Biaxial: {e}')

        st.subheader('Geometria de la seccion')
        col1, col2 = st.columns([3, 2])
        with col1:
            st.pyplot(plot_section(section), use_container_width=True)
        with col2:
            desc = [
                'Area bruta de concreto',
                'Area total de acero longitudinal',
                'Cuantia de refuerzo = Ast/Ag',
                'Capacidad en compresion pura (sin excentricidad)',
                'Maximo nominal admisible = 0.80 x P0',
                'Capacidad de diseno = 0.65 x Pn,max',
            ]
            data = {
                'Propiedad': ['Ag', 'Ast', '\u03c1', 'P0 nominal', 'Pn,max nominal', '\u03c6Pn,max'],
                'Valor': [
                    f'{section.ag() / cm**2:.2f} cm\u00b2',
                    f'{section.ast() / cm**2:.2f} cm\u00b2',
                    f'{section.rho():.3%}',
                    f'{section.p0() / ton:,.2f} tonf',
                    f'{section.pn_max() / ton:,.2f} tonf',
                    f'{0.65 * section.pn_max() / ton:,.2f} tonf',
                ],
                'Descripcion': desc,
            }
            st.table(pd.DataFrame(data))
            st.subheader('Matriz de barras (mm)')
            st.dataframe(pd.DataFrame(raw_matrix), use_container_width=True, hide_index=True)

        with st.expander('Supuestos tecnicos del calculo'):
            st.markdown(
                '- **\u03b5cu = 0.003**: Deformacion maxima del concreto en compresion segun ACI 318-19.\n'
                '- **Bloque rectangular de Whitney**: Distribucion simplificada del esfuerzo de compresion '
                'en el concreto (0.85\u00b7f\'c en un area de profundidad a = \u03b2\u2081\u00b7c).\n'
                '- **Acero elastoplastico perfecto**: El acero tiene comportamiento elastico lineal hasta '
                'la fluencia (fy), luego plastico sin endurecimiento.\n'
                '- **\u03c6 variable**: 0.65 en compresion controlada, 0.90 en traccion controlada, '
                'transicion lineal segun ACI 318-19 \u00a721.2.\n'
                '- **Columnas con estribos**: Factor de confinamiento 0.80 para Pn,max (ACI \u00a722.4.2).'
            )

        if esp_msgs:
            with st.expander('Errores de espaciamiento', expanded=True):
                for msg in esp_msgs:
                    st.error(msg)
