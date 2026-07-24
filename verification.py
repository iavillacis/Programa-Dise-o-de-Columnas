# =============================================================================
# MODULO DE VERIFICACION ACI 318-19
# =============================================================================
# Este modulo contiene las funciones que verifican si una columna
# rectangular cumple con los requisitos del codigo ACI 318-19.
# Cada funcion revisa un aspecto especifico y devuelve un diccionario
# con los resultados.
# =============================================================================

import math

# FA = factor de area para barras: pi/4 / 100
# (pi/4 * diametro_mm^2) produce area en mm2; dividir /100 da cm2
FA = math.pi / 4 / 100


# =============================================================================
# FUNCION AUXILIAR: area_barra
# =============================================================================
# Calcula el area de una barra de acero en cm2 dado su diametro en mm.
# Area = (pi/4) * diametro^2, convertido de mm2 a cm2.
# =============================================================================
def area_barra(diametro_mm):
    return FA * diametro_mm ** 2


# =============================================================================
# 1. GEOMETRIA DE LA COLUMNA
# =============================================================================
# Calcula las propiedades geometricas basicas:
#   - Ag_real = Area bruta de concreto (B x H)
#   - bc = Ancho confinado (B - 2*recubrimiento)
#   - Ac = Area confinada (bc x (H - 2*rec))
# =============================================================================
def geometria_columna(B, H, rec):
    Ag_real = B * H
    bc = B - 2 * rec
    Ac = bc * (H - 2 * rec)
    return {"Ag_real": Ag_real, "bc": bc, "Ac": Ac}


# =============================================================================
# 2. ACERO LONGITUDINAL
# =============================================================================
# Cuenta el numero total de barras y el area total de acero longitudinal.
# n_var_B = barras en direccion B, n_var_H = barras en direccion H.
# Las 4 esquinas usan el diametro de esquina; el resto usa el longitudinal.
# =============================================================================
def acero_longitudinal(n_var_B, n_var_H, phi_long_mm, phi_esq_mm):
    n_total = n_var_B * 2 + (n_var_H - 2) * 2
    As_col = 4 * area_barra(phi_esq_mm) + (n_total - 4) * area_barra(phi_long_mm)
    return {"n_total": n_total, "As_col": As_col}


# =============================================================================
# 3. CUANTIA MINIMA DE ACERO
# =============================================================================
# Verifica que la cuantia de refuerzo longitudinal (rho = As/Ag) sea
# mayor o igual al minimo especificado por ACI 318-19 (1% por defecto).
# Si cumple -> OK, si no -> se muestra el area de acero minima requerida.
# =============================================================================
def verificar_cuantia_minima(As_col, B, H, pmin=0.01):
    Ag_real = B * H
    p_col = As_col / Ag_real
    cumple = p_col >= pmin
    As_min = pmin * B * H
    return {"p_col": p_col, "cumple": cumple, "As_min": As_min}


# =============================================================================
# 4. ACERO TRANSVERSAL (ESTRIBOS) - Ash REQUERIDO
# =============================================================================
# Calcula el area de acero transversal (Ash) necesaria segun ACI 318-19
# Ecuaciones:
#   Ash1 = 0.30 * bc * s * (fc/fy) * (Ag/Ac - 1)
#   Ash2 = 0.09 * bc * s * (fc/fy)
# El Ash requerido es el mayor entre Ash1 y Ash2.
# bc = ancho confinado, s = separacion propuesta, fc = resistencia
# del concreto, fy = fluencia del acero, Ag = area bruta,
# Ac = area confinada.
# =============================================================================
def acero_transversal(bc, s, fc, fy, Ag, Ac):
    Ash1 = 0.3 * bc * s * (fc / fy) * (Ag / Ac - 1)
    Ash2 = 0.09 * bc * s * (fc / fy)
    return {"Ash1": Ash1, "Ash2": Ash2, "Ash": max(Ash1, Ash2)}


# =============================================================================
# 5. ZONA PROTEGIDA Lo
# =============================================================================
# Calcula la longitud de la zona protegida (Lo) segun ACI 318-19.
# Es la zona cerca de la union columna-viga donde los estribos deben
# estar mas cerrados. Lo es el mayor entre:
#   - 45 cm (minimo ACI)
#   - b (la mayor dimension de la columna)
#   - L/6 (un sexto de la longitud libre)
# =============================================================================
def longitud_zona_protegida(B, H, L):
    b = max(B, H)
    return max(45, b, L / 6)


# =============================================================================
# 6. NUMERO DE RAMALES DEL ESTRIBO
# =============================================================================
# Calcula cuantos ramales (patas) de estribo se necesitan para
# proporcionar el Ash requerido. n_exacto = Ash / area de 1 estribo.
# n_usar es el entero superior (se redondea hacia arriba).
# =============================================================================
def numero_ramales(Ash, phi_estribo_mm):
    area_est = area_barra(phi_estribo_mm)
    n_exacto = Ash / area_est
    return {"n_exacto": n_exacto, "n_usar": math.ceil(n_exacto)}


# =============================================================================
# 7. SEPARACION MAXIMA DE ESTRIBOS
# =============================================================================
# Calcula la separacion maxima permitida entre estribos segun ACI:
#   - Dentro de Lo: s <= 6*dv (donde dv = diametro menor de la
#     varilla longitudinal) y no mas de 10 cm
#   - Fuera de Lo: s <= 6*dv y no mas de 15 cm
# =============================================================================
def separacion_estribos(diametro_menor_varilla_mm):
    dv_cm = diametro_menor_varilla_mm / 10
    dentro_Lo = min(6 * dv_cm, 10)
    fuera_Lo = min(6 * dv_cm, 15)
    return {"dentro_Lo": dentro_Lo, "fuera_Lo": fuera_Lo}


# =============================================================================
# 8. SEPARACION DEL ACERO LONGITUDINAL
# =============================================================================
# Calcula la separacion libre entre barras longitudinales en una
# direccion dada. Se usa para verificar que la separacion entre
# barras cumpla con los requisitos ACI de espaciamiento minimo.
# =============================================================================
def separacion_acero_longitudinal(bc, phi_estribo_mm, phi_esq_mm, phi_long_mm, n_var_lado):
    num = (10 * bc - 2 * phi_estribo_mm - 2 * phi_esq_mm
           - (n_var_lado - 2) * phi_long_mm)
    den = 10 * (n_var_lado - 1)
    return num / den


# =============================================================================
# FUNCION PRINCIPAL: verificar_columna
# =============================================================================
# Ejecuta TODAS las verificaciones ACI 318-19 de una sola vez:
# geometria, acero longitudinal, cuantia minima, acero transversal,
# zona protegida Lo, numero de ramales, separacion de estribos, y
# separacion del acero longitudinal.
#
# Parametros de entrada:
#   B, H = dimensiones de la columna (cm)
#   rec = recubrimiento (cm)
#   n_var_B, n_var_H = numero de barras en cada direccion
#   phi_long_mm, phi_esq_mm = diametros de acero longitudinal y esquinas (mm)
#   phi_estribo_mm = diametro del estribo (mm)
#   fc = resistencia del concreto (kgf/cm2)
#   fy = fluencia del acero (kgf/cm2)
#   L = longitud libre de la columna (cm)
#   s_estribo = separacion propuesta entre estribos (cm)
#   pmin = cuantia minima (0.01 = 1% por defecto)
#
# Devuelve: diccionario con todos los resultados de verificacion
# =============================================================================
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
