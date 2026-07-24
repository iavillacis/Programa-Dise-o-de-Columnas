# Manual de Usuario — Columnas de Concreto Armado ACI 318-19

## 1. ¿Qué hace este programa?

Calcula la resistencia de **columnas rectangulares de concreto armado** con estribos, siguiendo el código **ACI 318-19**. Genera:

- Diagramas de interacción **P–M** (carga axial vs momento flector)
- Verificación **biaxial** (Mx + My combinados)
- Revisión de todos los requisitos del **ACI 318-19** (cuantías, estribos, separaciones, zona protegida Lo)

---

## 2. Cómo ejecutar el programa

### Opción 1: Local (en tu computadora)

```bash
python -m streamlit run app.py
```

Se abrirá una ventana en tu navegador en `http://localhost:8501`.

### Opción 2: Railway (en la nube)

Accede a la URL pública que te proporcionó Railway (ej: `https://tu-app.railway.app`). No necesitas instalar nada.

---

## 3. Pantalla principal

La interfaz tiene dos zonas:

### Barra lateral izquierda (ingreso de datos)

Todos los parámetros se ingresan aquí, organizados en secciones:

#### Geometría
| Parámetro | Descripción | Ejemplo |
|-----------|-------------|---------|
| **B** | Base de la columna (cm) | 40 |
| **H** | Peralte de la columna (cm) | 40 |
| **Rec. r** | Recubrimiento libre (cm) | 4.0 |
| **Estribo** | Diámetro del estribo (mm) | 9.5 (#3) |
| **L** | Longitud libre de la columna (cm) | 300 |
| **s** | Separación propuesta entre estribos (cm) | 10 |

#### Materiales
| Parámetro | Descripción | Ejemplo |
|-----------|-------------|---------|
| **f'c** | Resistencia del concreto (kgf/cm²) | 280 |
| **fy** | Fluencia del acero (kgf/cm²) | 4200 |

#### Armado longitudinal
| Parámetro | Descripción | Ejemplo |
|-----------|-------------|---------|
| **N° varillas en B** | Barras en la dirección horizontal | 4 |
| **N° varillas en H** | Barras en la dirección vertical | 4 |
| **Ø acero longitudinal** | Diámetro de las barras de borde (mm) | 22 |
| **Ø acero esquinas** | Diámetro de las barras de esquina (mm) | 25 |

> **Importante:** Solo el perímetro tiene barras. El centro queda sin acero para confinar el núcleo de concreto.

#### Demandas
| Parámetro | Descripción | Ejemplo |
|-----------|-------------|---------|
| **Pu** | Carga axial última (tonf) | 300 |
| **Mux** | Momento flector último en eje X (tonf·m) | 25 |
| **Muy** | Momento flector último en eje Y (tonf·m) | 15 |

### Botón CALCULAR

Una vez ingresados todos los datos, presiona **CALCULAR** en el centro de la pantalla.

---

## 4. Resultados

Los resultados se muestran en 5 secciones expandibles:

### Expander 1: Geometría y Propiedades de la sección

Muestra una tabla con:

| Propiedad | Significado |
|-----------|-------------|
| **Ag** | Área bruta de concreto (cm²) |
| **Ast** | Área total de acero longitudinal (cm²) |
| **ρ** | Cuantía de refuerzo = Ast/Ag |
| **P0 nominal** | Capacidad en compresión pura (tonf) |
| **Pn,max nominal** | Máximo nominal admisible = 0.80 × P0 |
| **φPn,max** | Capacidad de diseño = 0.65 × Pn,max |

Debajo se muestra la **matriz de barras** con los diámetros en mm.

### Expander 2: Diagramas P–M

Dos gráficos (uno para cada eje):
- **Curva negra punteada**: resistencia nominal (Pn, Mn)
- **Curva azul**: resistencia de diseño (φPn, φMn)
- **Línea naranja**: φPn,max (límite ACI)
- **Punto rojo**: demanda (Pu, Mu)

La columna **cumple** si el punto rojo está dentro de la curva azul.

### Expander 3: Verificación ACI 318-19

8 revisiones numeradas:

| # | Revisión | Qué significa |
|---|----------|---------------|
| 1 | Geometría | Ag, bc, Ac |
| 2 | Acero longitudinal | N° varillas y área total |
| 3 | Cuantía mínima | ρ ≥ 1% |
| 4 | Ash | Área de estribos requerida |
| 5 | Lo | Longitud de zona protegida |
| 6 | Ramales | Número de patas del estribo |
| 7 | Separación de estribos | Máximo permitido dentro/fuera de Lo |
| 8 | Separación acero longitudinal | Distancia libre entre barras |

### Expander 4: Verificación biaxial

Aparece solo cuando `Muy > 0.01`. Muestra:
- **Gráfico de contorno** (Mx vs My a Pu constante)
- **φMnx** y **φMny** uniaxiales
- **Relación D/C radial**: si ≤ 1, la columna **CUMPLE**; si > 1, **FALLA**

### Expander 5: Supuestos técnicos

Explica las hipótesis de cálculo del ACI 318-19:
- εcu = 0.003 (deformación máxima del concreto)
- Bloque rectangular de Whitney
- Acero elastoplástico perfecto
- φ variable según deformación
- Factor de confinamiento 0.80

---

## 5. Interpretación de resultados

### La columna es adecuada si:
- El punto rojo está **dentro** de las curvas azules en los diagramas P–M
- La cuantía está entre **1% y 8%**
- La verificación ACI muestra **"CUMPLE"** en cuantía mínima
- **D/C radial ≤ 1** en la verificación biaxial (si aplica)

### La columna NO es adecuada si:
- El punto rojo está **fuera** de alguna curva azul
- La cuantía es **menor a 1%** o **mayor a 8%**
- **D/C radial > 1** en biaxial
- Hay errores de espaciamiento entre barras

---

## 6. Consejos prácticos

- Empieza con **B = H = 40 cm**, **4 barras por lado**, **φ22 mm** y **φ25 mm** en esquinas
- Para columnas más cargadas, aumenta el número de barras o el diámetro
- La separación de estribos **s = 10 cm** suele funcionar bien
- Si el D/C biaxial es > 1, prueba aumentar la sección o el acero
- La **zona protegida Lo** se calcula automáticamente; los estribos deben estar más cerrados dentro de ella

---

## 7. Solución de problemas

| Problema | Posible causa | Solución |
|----------|---------------|----------|
| "Debe existir al menos una barra" | n_B o n_H muy pequeños | Usa al menos 2×2 |
| "Recubrimiento deja sin núcleo" | rec demasiado grande | Reduce rec o aumenta B/H |
| "Pu excede φPn,max" | Columna muy pequeña para la carga | Aumenta B, H o f'c |
| "No se pudo formar el contorno" | Muy inconsistente con Pu | Revisa datos de demanda |
| Gráficos no se ven | Error de Pillow o memoria | Haz "Refresh" en el navegador |
