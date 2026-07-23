import math

FA = math.pi / 4 / 100


def area_barra(diametro_mm):
    return FA * diametro_mm ** 2


def geometria_columna(B, H, rec):
    Ag_real = B * H
    bc = B - 2 * rec
    Ac = bc * (H - 2 * rec)
    return {"Ag_real": Ag_real, "bc": bc, "Ac": Ac}


def acero_longitudinal(n_var_B, n_var_H, phi_long_mm, phi_esq_mm):
    n_total = n_var_B * 2 + (n_var_H - 2) * 2
    As_col = 4 * area_barra(phi_esq_mm) + (n_total - 4) * area_barra(phi_long_mm)
    return {"n_total": n_total, "As_col": As_col}


def verificar_cuantia_minima(As_col, B, H, pmin=0.01):
    Ag_real = B * H
    p_col = As_col / Ag_real
    cumple = p_col >= pmin
    As_min = pmin * B * H
    return {"p_col": p_col, "cumple": cumple, "As_min": As_min}


def acero_transversal(bc, s, fc, fy, Ag, Ac):
    Ash1 = 0.3 * bc * s * (fc / fy) * (Ag / Ac - 1)
    Ash2 = 0.09 * bc * s * (fc / fy)
    return {"Ash1": Ash1, "Ash2": Ash2, "Ash": max(Ash1, Ash2)}


def longitud_zona_protegida(B, H, L):
    b = max(B, H)
    return max(45, b, L / 6)


def numero_ramales(Ash, phi_estribo_mm):
    area_est = area_barra(phi_estribo_mm)
    n_exacto = Ash / area_est
    return {"n_exacto": n_exacto, "n_usar": math.ceil(n_exacto)}


def separacion_estribos(diametro_menor_varilla_mm):
    dv_cm = diametro_menor_varilla_mm / 10
    dentro_Lo = min(6 * dv_cm, 10)
    fuera_Lo = min(6 * dv_cm, 15)
    return {"dentro_Lo": dentro_Lo, "fuera_Lo": fuera_Lo}


def separacion_acero_longitudinal(bc, phi_estribo_mm, phi_esq_mm, phi_long_mm, n_var_lado):
    num = (10 * bc - 2 * phi_estribo_mm - 2 * phi_esq_mm
           - (n_var_lado - 2) * phi_long_mm)
    den = 10 * (n_var_lado - 1)
    return num / den


def verificar_columna(B, H, rec, n_var_B, n_var_H, phi_long_mm, phi_esq_mm,
                       phi_estribo_mm, fc, fy, L, s_estribo, pmin=0.01):
    geo = geometria_columna(B, H, rec)
    long_ = acero_longitudinal(n_var_B, n_var_H, phi_long_mm, phi_esq_mm)
    cuantia = verificar_cuantia_minima(long_["As_col"], B, H, pmin)
    transv = acero_transversal(geo["bc"], s_estribo, fc, fy,
                                geo["Ag_real"], geo["Ac"])
    Lo = longitud_zona_protegida(B, H, L)
    ramales = numero_ramales(transv["Ash"], phi_estribo_mm)
    sep_estr = separacion_estribos(min(phi_long_mm, phi_esq_mm))
    sep_long_B = separacion_acero_longitudinal(geo["bc"], phi_estribo_mm,
                                                phi_esq_mm, phi_long_mm, n_var_B)
    sep_long_H = separacion_acero_longitudinal(geo["bc"], phi_estribo_mm,
                                                phi_esq_mm, phi_long_mm, n_var_H)

    return {
        "geometria": geo,
        "acero_longitudinal": long_,
        "cuantia_minima": cuantia,
        "acero_transversal": transv,
        "Lo": Lo,
        "ramales": ramales,
        "separacion_estribos": sep_estr,
        "separacion_long_B": sep_long_B,
        "separacion_long_H": sep_long_H,
    }
