import numpy as np
from io import BytesIO
from base64 import b64encode
from html import escape
from urllib.parse import parse_qs
from string import Template
import warnings
import os
import gc
import sys

import verification

_plt_mod = None

def _plt():
    global _plt_mod
    if _plt_mod is None:
        if 'matplotlib.pyplot' in sys.modules:
            _plt_mod = sys.modules['matplotlib.pyplot']
        else:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            _plt_mod = plt
    return _plt_mod

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

    def curva_interaccion(self, eje="x", n=40):
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

    def contorno_aci(self, pu, n_theta=24):
        if pu <= 0:
            raise ValueError("Pu debe ser mayor que cero (compresion).")
        phi_pmax = 0.65 * self.pn_max()
        if pu > phi_pmax:
            raise ValueError(f"Pu excede \u03c6Pn,max ({phi_pmax / ton:.2f} tonf).")
        contour = []
        for theta in np.linspace(0, 2 * np.pi, n_theta, endpoint=False):
            projection = abs(np.cos(theta)) * self.b + abs(np.sin(theta)) * self.h
            c_values = np.geomspace(max(projection * 1e-5, 1e-7), projection * 50, 60)
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
    fig, ax = _plt().subplots(figsize=(7, 5.2))
    mn_neg = [-x for x in mn]
    mp_neg = [-x for x in mp]
    ax.plot([x / ton for x in mn], [x / ton for x in pn], 'k--', lw=1.5, label='Nominal')
    ax.plot([x / ton for x in mn_neg], [x / ton for x in pn], 'k--', lw=1.5)
    ax.plot([x / ton for x in mp], [x / ton for x in pp], 'b-', lw=2, label='\u03c6 (dise\u00f1o)')
    ax.plot([x / ton for x in mp_neg], [x / ton for x in pp], 'b-', lw=2)
    ax.axhline(0.65 * pmax / ton, color='orange', lw=1.3, label='\u03c6Pn,max')
    ax.axhline(0, color='black', lw=0.7)
    ax.axvline(0, color='black', lw=0.7)
    ax.scatter([mu / ton], [pu / ton], color='red', s=56, zorder=5, label='Demanda')
    ax.grid(color='silver', linestyle='--', linewidth=0.5)
    title = 'Mx' if eje == 'x' else 'My'
    ax.set(xlabel=f'{title} (tonf\u00b7m)', ylabel='P (tonf)', title=f'Diagrama P\u2013{title}')
    ax.legend(fontsize=8, loc='lower right')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fig.tight_layout()
    return fig


def plot_contorno_aci(contour, mux, muy):
    closed = np.vstack((contour, contour[:1]))
    fig, ax = _plt().subplots(figsize=(6.5, 5.5))
    ax.fill(closed[:, 0] / ton, closed[:, 1] / ton, color='#bae6fd', alpha=0.55, label='Region resistente')
    ax.plot(closed[:, 0] / ton, closed[:, 1] / ton, color='#0284c7', lw=2, label='Contorno ACI \u03c6P = Pu')
    ax.scatter([mux / ton], [muy / ton], color='red', s=52, zorder=5, label='Demanda')
    ax.axhline(0, color='black', lw=0.7)
    ax.axvline(0, color='black', lw=0.7)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(color='silver', linestyle='--', linewidth=0.5)
    ax.set(xlabel='Mx (tonf\u00b7m)', ylabel='My (tonf\u00b7m)', title='Contorno biaxial \u03c6Mx vs \u03c6My a Pu dado')
    ax.legend(fontsize=8)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fig.tight_layout()
    return fig


def plot_section(section):
    _p = _plt()
    fig, ax = _p.subplots(figsize=(10, 8))
    ax.add_patch(_p.Rectangle((0, 0), section.b / cm, section.h / cm,
                                fill=False, lw=3, color='#1e293b'))
    for x, y, area, dm, d in section._bars:
        r = d / (2 * cm)
        ax.add_patch(_p.Circle((x / cm, y / cm), r, color='#b91c1c', alpha=0.85))
        ax.text(x / cm, y / cm, f'{dm:.1f}', color='white',
                ha='center', va='center', fontsize=9, fontweight='bold')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set(xlabel='B (cm)', ylabel='H (cm)', title='Distribucion de acero longitudinal (diametros en mm)')
    ax.grid(alpha=0.2)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fig.tight_layout()
    return fig


def fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=85, bbox_inches='tight')
    _plt().close(fig)
    return b64encode(buf.getvalue()).decode('ascii')


def generar_matriz_barras(n_B, n_H, diam_long, diam_corner):
    mat = np.full((n_H, n_B), diam_long)
    if n_B >= 1 and n_H >= 1:
        mat[0, 0] = diam_corner
    if n_B >= 2:
        mat[0, -1] = diam_corner
    if n_H >= 2:
        mat[-1, 0] = diam_corner
    if n_B >= 2 and n_H >= 2:
        mat[-1, -1] = diam_corner
    return mat


def matriz_to_html(mat):
    rows, cols = mat.shape
    html = '<table class="matrix">'
    for i in range(rows):
        html += '<tr>'
        for j in range(cols):
            val = float(mat[i, j])
            cls = 'corner' if val > 0 and (
                (i == 0 and j == 0) or
                (i == 0 and j == cols - 1) or
                (i == rows - 1 and j == 0) or
                (i == rows - 1 and j == cols - 1)
            ) else ''
            html += f'<td class="{cls}">{val:.1f}</td>'
        html += '</tr>'
    html += '</table>'
    return html


def _parse_matrix(raw):
    rows = [[float(x.strip()) for x in row.split(',')] for row in raw.split(';') if row.strip()]
    if not rows or any(len(r) != len(rows[0]) for r in rows):
        raise ValueError("Matriz invalida: use filas separadas por ; y columnas por ,")
    return np.array(rows, dtype=float)


def _verif_to_html(vr):
    geo = vr["geometria"]
    la = vr["acero_longitudinal"]
    cm_check = vr["cuantia_minima"]
    at = vr["acero_transversal"]
    ram = vr["ramales"]
    se = vr["separacion_estribos"]

    html = '''
<div class="card">
<h2>Verificación ACI 318-19</h2>
<div class="summary-grid">
<div class="summary-item"><div class="label">Ag</div><div class="value">{ag:.2f} cm²</div></div>
<div class="summary-item"><div class="label">bc</div><div class="value">{bc:.2f} cm</div></div>
<div class="summary-item"><div class="label">Ac</div><div class="value">{ac:.2f} cm²</div></div>
<div class="summary-item"><div class="label">N° varillas</div><div class="value">{n_total}</div></div>
<div class="summary-item"><div class="label">As total</div><div class="value">{as_col:.2f} cm²</div></div>
<div class="summary-item"><div class="label">ρ</div><div class="value">{pcol:.2%}</div><div class="status {cm_cls}">{cm_txt}</div></div>
<div class="summary-item"><div class="label">Ash req.</div><div class="value">{ash:.2f} cm²</div></div>
<div class="summary-item"><div class="label">Ramales</div><div class="value">{n_ram}</div></div>
<div class="summary-item"><div class="label">Lo</div><div class="value">{lo:.1f} cm</div></div>
</div>
<h3>Separación de estribos</h3>
<p>Dentro de Lo: s &lt; {sdentro:.2f} cm &nbsp;|&nbsp; Fuera de Lo: s &lt; {sfuera:.2f} cm</p>
<h3>Separación acero longitudinal</h3>
<p>Dirección B: {sepB:.2f} cm &nbsp;|&nbsp; Dirección H: {sepH:.2f} cm</p>
</div>'''.format(
        ag=geo["Ag_real"], bc=geo["bc"], ac=geo["Ac"],
        n_total=la["n_total"], as_col=la["As_col"],
        pcol=cm_check["p_col"],
        cm_cls="ok" if cm_check["cumple"] else "fail",
        cm_txt="CUMPLE ρ≥1%" if cm_check["cumple"] else "NO CUMPLE ρ<1%",
        ash=at["Ash"], n_ram=ram["n_usar"],
        lo=vr["Lo"],
        sdentro=se["dentro_Lo"], sfuera=se["fuera_Lo"],
        sepB=vr["separacion_long_B"], sepH=vr["separacion_long_H"]
    )
    return html


TEMPLATE_PATH = os.path.join(os.path.dirname(__file__), 'templates', 'report.html')


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
        n_B = int(qs.get('n_B', ['4'])[0])
        n_H = int(qs.get('n_H', ['4'])[0])
        diam_long = float(qs.get('diam_long', ['22'])[0])
        diam_corner = float(qs.get('diam_corner', ['25'])[0])
        L = float(qs.get('L', ['300'])[0])
        s_estribo = float(qs.get('s_estribo', ['10'])[0])

        raw_mat = qs.get('armado', [None])[0]
        if raw_mat:
            mat = _parse_matrix(raw_mat)
        else:
            mat = generar_matriz_barras(n_B, n_H, diam_long, diam_corner)

        sec = DiagramaInteraccion(b * cm, h * cm, rec * cm, tie_mm * mm,
                                  fc * kgf / (cm**2), fy * kgf / (cm**2), mat)

        esp_msgs = sec.verificar_espaciamiento()
        esp_html = ''
        if esp_msgs:
            esp_html = '<div class="error">' + '<br>'.join(escape(m) for m in esp_msgs) + '</div>'

        fig1 = plot_diagram(sec, 'x', pu * ton, mux * ton * m)
        fig2 = plot_diagram(sec, 'y', pu * ton, muy * ton * m)
        fig3 = plot_section(sec)
        img_x = fig_to_b64(fig1)
        img_y = fig_to_b64(fig2)
        img_sec = fig_to_b64(fig3)

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
                status_cls = 'ok' if ratio <= 1 else 'fail'
                biaxial_html = f'''
                <div class="card">
                <h2>Verificación biaxial ACI</h2>
                <div class="grid-2">
                <div><img src="data:image/png;base64,{img4}" style="width:100%"></div>
                <div>
                <div class="summary-grid">
                <div class="summary-item"><div class="label">ϕMnx</div><div class="value">{mx_uni:,.3f} tonf·m</div></div>
                <div class="summary-item"><div class="label">ϕMny</div><div class="value">{my_uni:,.3f} tonf·m</div></div>
                <div class="summary-item"><div class="label">D/C radial</div><div class="value">{ratio:.3f}</div><div class="status {status_cls}">{status}</div></div>
                </div></div></div></div>
                '''
            except Exception as e:
                biaxial_html = f'<div class="error">Biaxial: {escape(str(e))}</div>'

        ag_val = sec.ag() / cm**2
        ast_val = sec.ast() / cm**2
        rho_val = sec.rho()
        p0_val = sec.p0() / ton
        pnmax_val = sec.pn_max() / ton
        phipmax_val = 0.65 * sec.pn_max() / ton

        matrix_html = matriz_to_html(mat)

        vr = verification.verificar_columna(
            B=b, H=h, rec=rec,
            n_var_B=n_B, n_var_H=n_H,
            phi_long_mm=diam_long, phi_esq_mm=diam_corner,
            phi_estribo_mm=tie_mm,
            fc=fc, fy=fy,
            L=L, s_estribo=s_estribo
        )
        verificacion_html = _verif_to_html(vr)

        with open(TEMPLATE_PATH, 'r', encoding='utf-8') as f:
            tmpl = Template(f.read())

        html = tmpl.substitute(
            b_cm=b, h_cm=h, rec_cm=rec, tie_mm=tie_mm,
            fc_val=fc, fy_val=fy,
            n_B=n_B, n_H=n_H, diam_long=diam_long, diam_corner=diam_corner,
            L_val=L, s_estribo=s_estribo,
            pu_val=pu, mux_val=mux, muy_val=muy,
            espaciado_html=esp_html,
            img_x=img_x, img_y=img_y, img_sec=img_sec,
            ag_val=f'{ag_val:.2f}', ast_val=f'{ast_val:.2f}',
            rho_val=f'{rho_val:.3%}',
            p0_val=f'{p0_val:,.2f}', pnmax_val=f'{pnmax_val:,.2f}',
            phipmax_val=f'{phipmax_val:,.2f}',
            biaxial_html=biaxial_html,
            verificacion_html=verificacion_html,
            matrix_html=matrix_html,
        )
    except Exception as e:
        html = f'<p class="error">Error: {escape(str(e))}</p>'
    _plt().close('all')
    gc.collect()
    body = html.encode('utf-8')
    start_response('200 OK', [('Content-Type', 'text/html; charset=utf-8')])
    return [body]


if __name__ == "__main__":
    import streamlit as st

    _plt()

    st.set_page_config(page_title='Columnas ACI 318-19', page_icon='\U0001f3d7\ufe0f', layout='wide')

    st.markdown("""
    <style>
    div[data-testid="column"]:nth-of-type(2) .stButton button {
        font-size: 1.4rem !important;
        padding: 0.75rem 2rem !important;
        font-weight: 700 !important;
        border-radius: 10px !important;
        box-shadow: 0 4px 14px rgba(2,132,199,0.35) !important;
        border: 2px solid #0284c7 !important;
        transition: all 0.2s ease !important;
    }
    div[data-testid="column"]:nth-of-type(2) .stButton button:hover {
        box-shadow: 0 6px 20px rgba(2,132,199,0.5) !important;
        transform: scale(1.02) !important;
    }
    .verif-card {
        background: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 8px;
        padding: 0.8rem 1rem;
        text-align: center;
    }
    .verif-card .label {
        font-size: 0.75rem;
        text-transform: uppercase;
        color: #64748b;
        letter-spacing: 0.03em;
    }
    .verif-card .value {
        font-size: 1.1rem;
        font-weight: 700;
        color: #0f172a;
        margin-top: 0.1rem;
    }
    .verif-card .status-pass {
        color: #16a34a;
        font-weight: 700;
        font-size: 0.85rem;
    }
    .verif-card .status-fail {
        color: #dc2626;
        font-weight: 700;
        font-size: 0.85rem;
    }
    .matrix-table {
        border-collapse: collapse;
        margin: 0 auto;
        font-size: 0.9rem;
    }
    .matrix-table td {
        width: 2.8rem;
        height: 2.8rem;
        text-align: center;
        border: 1px solid #94a3b8;
        font-weight: 600;
    }
    .matrix-table td.corner {
        background: #dbeafe;
        color: #1e40af;
    }
    .matrix-table td.interior {
        background: #fefce8;
    }
    .st-emotion-cache-1r4qj8v {
        font-size: 1rem;
    }
    </style>
    """, unsafe_allow_html=True)

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
                                 help='Diametro del estribo. #3=9.5mm, #4=12.7mm, #5=15.9mm')
        L = st.number_input('L (cm) - Longitud libre columna', min_value=100.0, value=300.0, step=10.0)
        s_estribo = st.number_input('s (cm) - Separacion estribos propuesta', min_value=1.0, value=10.0, step=0.5)

        st.subheader('Materiales')
        fc = st.number_input("f'c (kgf/cm\u00b2)", min_value=100.0, value=280.0, step=10.0)
        fy = st.number_input('fy (kgf/cm\u00b2)', min_value=2000.0, value=4200.0, step=100.0)

        st.subheader('Armado longitudinal')
        st.caption('Las esquinas usan el diametro de esquina; el resto usa el diametro longitudinal.')
        col_a, col_b = st.columns(2)
        with col_a:
            n_B = st.number_input('N\u00b0 varillas en B', min_value=2, max_value=12, value=4, step=1)
        with col_b:
            n_H = st.number_input('N\u00b0 varillas en H', min_value=2, max_value=12, value=4, step=1)
        diam_long = st.number_input('\u00d8 acero longitudinal (mm)', min_value=6.0, max_value=50.0, value=22.0, step=1.0)
        diam_corner = st.number_input('\u00d8 acero esquinas (mm)', min_value=6.0, max_value=50.0, value=25.0, step=1.0)

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
            raw_matrix = generar_matriz_barras(n_B, n_H, diam_long, diam_corner)

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

        # --- Expander 1: Geometria y Propiedades ---
        with st.expander('Geometr\u00eda y Propiedades de la secci\u00f3n', expanded=True):
            col1, col2 = st.columns([7, 3])
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
                st.table(data)

                st.markdown('**Matriz de barras (mm)**')
                st.dataframe(raw_matrix, use_container_width=True, hide_index=True)

        # --- Expander 2: Diagramas P-M ---
        with st.expander('Diagramas de interacci\u00f3n P\u2013M', expanded=True):
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

        # --- Expander 3: Verificacion ACI 318-19 ---
        with st.expander('Verificaci\u00f3n ACI 318-19', expanded=True):
            vr = verification.verificar_columna(
                B=b, H=h, rec=rec,
                n_var_B=n_B, n_var_H=n_H,
                phi_long_mm=diam_long, phi_esq_mm=diam_corner,
                phi_estribo_mm=tie_mm,
                fc=fc, fy=fy,
                L=L, s_estribo=s_estribo
            )
            geo = vr["geometria"]
            la = vr["acero_longitudinal"]
            cm_check = vr["cuantia_minima"]
            at = vr["acero_transversal"]
            ram = vr["ramales"]
            se = vr["separacion_estribos"]

            st.markdown('**1. Geometr\u00eda de la columna**')
            c1, c2, c3 = st.columns(3)
            c1.metric('Ag (area bruta)', f'{geo["Ag_real"]:.2f} cm\u00b2')
            c2.metric('bc (ancho confinado)', f'{geo["bc"]:.2f} cm')
            c3.metric('Ac (area confinada)', f'{geo["Ac"]:.2f} cm\u00b2')

            st.markdown('**2. Acero longitudinal**')
            c1, c2 = st.columns(2)
            c1.metric('N\u00b0 total de varillas', str(la["n_total"]))
            c2.metric('As total', f'{la["As_col"]:.2f} cm\u00b2')

            st.markdown('**3. Cuant\u00eda m\u00ednima**')
            c1, c2, c3 = st.columns(3)
            c1.metric('\u03c1 (Ast/Ag)', f'{cm_check["p_col"]:.3%}')
            if cm_check["cumple"]:
                c2.success('CUMPLE \u03c1 \u2265 1%')
            else:
                c2.error(f'NO CUMPLE (min {cm_check["As_min"]:.2f} cm\u00b2)')
            c3.metric('Estado', 'OK' if cm_check["cumple"] else 'Aumentar As')

            st.markdown('**4. Acero transversal (estribos) \u2014 Ash**')
            c1, c2, c3 = st.columns(3)
            c1.metric('Ash1', f'{at["Ash1"]:.3f} cm\u00b2')
            c2.metric('Ash2', f'{at["Ash2"]:.3f} cm\u00b2')
            c3.metric('Ash (mayor)', f'{at["Ash"]:.3f} cm\u00b2')

            st.markdown('**5. Zona protegida Lo**')
            st.metric('Lo', f'{vr["Lo"]:.1f} cm')

            st.markdown('**6. N\u00famero de ramales del estribo**')
            c1, c2 = st.columns(2)
            c1.metric('N\u00b0 exacto', f'{ram["n_exacto"]:.2f}')
            c2.metric('N\u00b0 a usar', str(ram["n_usar"]))

            st.markdown('**7. Separaci\u00f3n m\u00e1xima de estribos**')
            c1, c2 = st.columns(2)
            c1.metric('Dentro de Lo', f's < {se["dentro_Lo"]:.2f} cm')
            c2.metric('Fuera de Lo', f's < {se["fuera_Lo"]:.2f} cm')

            st.markdown('**8. Separaci\u00f3n del acero longitudinal**')
            c1, c2 = st.columns(2)
            c1.metric('Direcci\u00f3n B', f'{vr["separacion_long_B"]:.2f} cm')
            c2.metric('Direcci\u00f3n H', f'{vr["separacion_long_H"]:.2f} cm')

        # --- Expander 4: Biaxial ---
        if abs(muy) > 0.01:
            with st.expander('Verificaci\u00f3n biaxial ACI 318-19', expanded=True):
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

        # --- Expander 5: Supuestos ---
        if esp_msgs:
            with st.expander('Errores de espaciamiento', expanded=True):
                for msg in esp_msgs:
                    st.error(msg)

        with st.expander('Supuestos t\u00e9cnicos del c\u00e1lculo', expanded=False):
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
