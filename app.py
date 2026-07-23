import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from io import BytesIO

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

    def define_acero(self):
        filas, columnas = self.matrix.shape
        matriz_areas = np.zeros((filas, columnas))
        for i in range(filas):
            for j in range(columnas):
                matriz_areas[i, j] = self.area(self.matrix[i, j])
        return matriz_areas


class VariablesGenerales:
    def esfuerzo_acero(self, fy, es):
        Es = 2.0 * 10**6 * (kgf / (cm**2))
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
            return 0.85 - (0.05 * (fc - 280) / 70.0)


class DiagramaInteraccion:
    def __init__(self, b, h, rec, estribo, fc, fy, rebar_matrix):
        self.b = b
        self.h = h
        self.rec = rec
        self.fc = fc
        self.fy = fy
        self.rebar = Rebar(rebar_matrix)
        self.vg = VariablesGenerales()
        self.estribo_num = int(estribo)
        self.tie_diameter = self.rebar.diametro(self.estribo_num)
        self.centroide = h / 2.0
        self._bars = self._bar_coordinates()

    def _bar_coordinates(self):
        rows, cols = self.rebar.matrix.shape
        barras = []
        x0 = self.rec + self.tie_diameter
        y0 = self.rec + self.tie_diameter
        dx = (self.b - 2 * x0) / max(cols - 1, 1)
        dy = (self.h - 2 * y0) / max(rows - 1, 1)
        for i in range(rows):
            for j in range(cols):
                n = int(self.rebar.matrix[i, j])
                if n == 0:
                    continue
                x = self.b / 2 if cols == 1 else x0 + j * dx
                y = self.h / 2 if rows == 1 else y0 + i * dy
                barras.append((x, y, self.rebar.area(n), n, self.rebar.diametro(n)))
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
        Es = 2.0 * 10**6 * (kgf / (cm**2))
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
        pmax = self.pn_max()
        return p_list, m_list, pp_list, mp_list, pmax

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

    ax.plot([x / ton for x in mn], [x / ton for x in pn], 'k--', lw=1.5, label='Capacidad nominal')
    ax.plot([x / ton for x in mn_neg], [x / ton for x in pn], 'k--', lw=1.5)
    ax.plot([x / ton for x in mp], [x / ton for x in pp], 'b-', lw=2, label='Capacidad \u03c6')
    ax.plot([x / ton for x in mp_neg], [x / ton for x in pp], 'b-', lw=2)

    phi_pmax = 0.65 * pmax
    ax.axhline(phi_pmax / ton, color='orange', lw=1.3, label='\u03c6Pn,max')
    ax.axhline(0, color='black', lw=0.7)
    ax.axvline(0, color='black', lw=0.7)

    ax.scatter([mu / ton], [pu / ton], color='red', s=50, zorder=5, label='Demanda')

    ax.grid(color='silver', linestyle='--', linewidth=0.5)
    title = 'Mx' if eje == 'x' else 'My'
    ax.set(xlabel=f'{title} (tonf\u00b7m)', ylabel='P (tonf)',
           title=f'Diagrama P\u2013{title}')
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def plot_contorno_aci(contour, mux, muy):
    closed = np.vstack((contour, contour[:1]))
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    ax.fill(closed[:, 0] / ton, closed[:, 1] / ton, color='#bae6fd', alpha=0.55, label='Regi\u00f3n resistente')
    ax.plot(closed[:, 0] / ton, closed[:, 1] / ton, color='#0284c7', lw=2, label='Contorno ACI \u03c6P = Pu')
    ax.scatter([mux / ton], [muy / ton], color='red', s=48, zorder=5, label='Demanda')
    ax.axhline(0, color='black', lw=0.7)
    ax.axvline(0, color='black', lw=0.7)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(color='silver', linestyle='--', linewidth=0.5)
    ax.set(xlabel='Mx (tonf\u00b7m)', ylabel='My (tonf\u00b7m)',
           title='Contorno biaxial ACI 318-19')
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def plot_section(section):
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.add_patch(plt.Rectangle((0, 0), section.b / cm, section.h / cm,
                                fill=False, lw=2, color='#1e293b'))
    for x, y, area, n, db in section._bars:
        ax.add_patch(plt.Circle((x / cm, y / cm), db / (2 * cm),
                                color='#b91c1c', alpha=0.9))
        ax.text(x / cm, y / cm, f'#{int(n)}', color='white',
                ha='center', va='center', fontsize=7)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set(xlabel='B (cm)', ylabel='H (cm)', title='Distribuci\u00f3n de acero')
    ax.grid(alpha=0.25)
    fig.tight_layout()
    return fig


def figure_to_base64(fig):
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=130, bbox_inches='tight')
    plt.close(fig)
    from base64 import b64encode
    return b64encode(buf.getvalue()).decode('ascii')


def default_matrix():
    return pd.DataFrame([[8, 8, 8, 8], [8, 0, 0, 8], [8, 0, 0, 8], [8, 8, 8, 8]],
                        columns=['C1', 'C2', 'C3', 'C4'])


st.set_page_config(page_title='Columnas ACI 318-19', page_icon='\U0001f3d7\ufe0f', layout='wide')
st.title('Columnas de concreto armado \u00b7 ACI 318-19')
st.caption('Diagramas P\u2013M por compatibilidad de deformaciones y verificaci\u00f3n biaxial ACI.')

with st.sidebar:
    st.header('Datos de entrada')
    st.caption('Unidades: cm, kgf/cm\u00b2, tonf, tonf\u00b7m.')
    b = st.number_input('B', min_value=10.0, value=40.0, step=1.0)
    h = st.number_input('H', min_value=10.0, value=40.0, step=1.0)
    rec = st.number_input('Recubrimiento r', min_value=1.0, value=4.0, step=0.5)
    tie = st.selectbox('Estribo', [2, 3, 4, 5], index=1,
                       format_func=lambda n: f'Barra #{n}')
    fc = st.number_input("f\'c", min_value=100.0, value=280.0, step=10.0)
    fy = st.number_input('fy', min_value=2000.0, value=4200.0, step=100.0)
    st.subheader('Acero longitudinal')
    st.caption('0 = vac\u00edo. Barras: #2\u2013#8, #10.')
    matrix = st.data_editor(default_matrix(), num_rows='dynamic',
                            use_container_width=True, key='armado',
                            column_config={c: st.column_config.NumberColumn(c, min_value=0, step=1)
                                           for c in default_matrix().columns})
    st.subheader('Demandas')
    pu = st.number_input('Pu', value=300.0, step=5.0)
    mux = st.number_input('Mux', value=25.0, step=1.0)
    muy = st.number_input('Muy', value=15.0, step=1.0)

try:
    raw_matrix = matrix.fillna(0).astype(int).to_numpy()
    section = DiagramaInteraccion(b * cm, h * cm, rec * cm, tie,
                                  fc * kgf / (cm**2), fy * kgf / (cm**2), raw_matrix)
    rho_val = section.rho()
    if rho_val > 0.08:
        st.error(f'Cuant\u00eda {rho_val:.2%}: excede el m\u00e1ximo de 8%.')
        st.stop()
    if rho_val < 0.01:
        st.warning(f'Cuant\u00eda longitudinal {rho_val:.2%}: menor al 1% m\u00ednimo ACI.')
except Exception as e:
    st.error(f'Error en datos: {e}')
    st.stop()

pu_i = pu * ton
mux_i = mux * ton * m
muy_i = muy * ton * m

st.subheader('Diagramas de interacci\u00f3n P\u2013M')
col1, col2 = st.columns(2)
with col1:
    st.pyplot(plot_diagram(section, 'x', pu_i, mux_i), use_container_width=True)
with col2:
    st.pyplot(plot_diagram(section, 'y', pu_i, muy_i), use_container_width=True)
st.caption('Curva nominal (negro discontinuo), curva \u03c6 (azul), l\u00edmite \u03c6Pn,max (naranja) y punto de demanda (rojo).')

es_biaxial = abs(muy) > 0.01
if es_biaxial:
    st.subheader('Verificaci\u00f3n biaxial ACI 318-19')
    try:
        contour = section.contorno_aci(pu_i)
        col1, col2 = st.columns([3, 2])
        with col1:
            fig = plot_contorno_aci(contour, mux_i, muy_i)
            st.pyplot(fig, use_container_width=True)
        with col2:
            mx_uni = float(np.max(np.abs(contour[:, 0]))) if len(contour) else 0
            my_uni = float(np.max(np.abs(contour[:, 1]))) if len(contour) else 0
            st.metric('\u03c6Mnx (uniaxial en Pu)', f'{mx_uni / ton:,.3f} tonf\u00b7m')
            st.metric('\u03c6Mny (uniaxial en Pu)', f'{my_uni / ton:,.3f} tonf\u00b7m')
            dx, dy = mux_i, muy_i
            radial = np.hypot(dx, dy)
            cap_radial = np.hypot(mx_uni, my_uni)
            ratio = radial / cap_radial if cap_radial > 0 else 999
            st.metric('Relaci\u00f3n D/C (radial)', f'{ratio:.3f}')
            if ratio <= 1:
                st.success('CUMPLE \u2014 La demanda est\u00e1 dentro del contorno ACI.')
            else:
                st.error('FALLA \u2014 La demanda excede el contorno ACI.')
            st.caption('Contorno de carga ACI 318-19 por compatibilidad de deformaciones.')
    except Exception as e:
        st.error(f'No se pudo calcular el contorno biaxial: {e}')

st.subheader('Geometr\u00eda de la secci\u00f3n')
c1, c2 = st.columns([1, 1])
with c1:
    st.pyplot(plot_section(section), use_container_width=True)
with c2:
    data = {
        'Propiedad': ['Ag', 'Ast', '\u03c1', 'P0 nominal', 'Pn,max nominal', '\u03c6Pn,max'],
        'Valor': [f'{section.ag() / cm**2:.2f} cm\u00b2',
                  f'{section.ast() / cm**2:.2f} cm\u00b2',
                  f'{section.rho():.3%}',
                  f'{section.p0() / ton:,.2f} tonf',
                  f'{section.pn_max() / ton:,.2f} tonf',
                  f'{0.65 * section.pn_max() / ton:,.2f} tonf'],
    }
    st.table(pd.DataFrame(data))
    st.subheader('Matriz de barras')
    st.dataframe(pd.DataFrame(raw_matrix), use_container_width=True, hide_index=True)

st.caption('\u03b5cu = 0.003, bloque rectangular de Whitney, acero elastopl\u00e1stico perfecto. Columnas con estribos (\u03c6 = 0.65 en compresi\u00f3n controlada).')
