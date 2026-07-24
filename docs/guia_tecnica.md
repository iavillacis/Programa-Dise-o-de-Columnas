# Guía Técnica — Programa de Diseño de Columnas ACI 318-19

## 1. Arquitectura del programa

```
ProgramaColumnas/
├── app.py                    # Código principal (Streamlit + WSGI)
├── verification.py           # Módulo de verificación ACI 318-19
├── requirements.txt          # Dependencias Python
├── Procfile                  # Comando de inicio para Railway
├── railway.json              # Configuración de Railway
├── .gitignore                # Archivos ignorados por git
├── .streamlit/
│   └── config.toml           # Tema visual de Streamlit
├── templates/
│   └── report.html           # Plantilla HTML para modo WSGI
└── docs/
    ├── flujograma.md          # Diagrama de flujo
    ├── manual_usuario.md      # Manual de usuario
    └── guia_tecnica.md        # Esta guía
```

---

## 2. Modos de ejecución

El programa tiene **dos modos** que comparten el mismo núcleo de cálculo:

### Modo Streamlit (interfaz gráfica)

Se activa al ejecutar: `streamlit run app.py`

En `app.py`, la sección `if __name__ == "__main__":` (línea ~760) contiene toda la lógica de la interfaz. Usa la librería **Streamlit** para crear formularios, botones y gráficos interactivos.

### Modo WSGI (servidor HTTP)

Se activa cuando un servidor WSGI (como gunicorn) llama a la función `wsgi_app()` (línea ~550). Esta función:

1. Lee parámetros desde la URL (`QUERY_STRING`)
2. Realiza los mismos cálculos que el modo Streamlit
3. Genera una página HTML completa usando la plantilla `templates/report.html`
4. Devuelve el HTML como respuesta HTTP

---

## 3. Flujo de cálculo

### 3.1 Ingreso de datos

Los parámetros de entrada se convierten a unidades internas:

| Entrada (usuario) | Unidad | Conversión | Variable interna |
|-------------------|--------|------------|------------------|
| B, H, rec | cm | × 0.01 | m |
| Ø estribo | mm | × 0.001 | m |
| f'c, fy | kgf/cm² | × 1 | kgf/m² |
| Pu | tonf | × 1000 | kgf |
| Mux, Muy | tonf·m | × 1000 | kgf·m |

### 3.2 Generación de matriz de barras

`generar_matriz_barras(n_B, n_H, diam_long, diam_corner)` crea una matriz `n_H × n_B` donde:

- **Esquinas** (4 posiciones) → `diam_corner`
- **Bordes** (perímetro no esquina) → `diam_long`
- **Centro** → 0 (sin barra)

Esto asegura que el acero solo confine el perímetro, dejando el núcleo de concreto sin barras interiores.

### 3.3 Clase DiagramaInteraccion

Esta clase (aproximadamente líneas 78–392) representa una columna rectangular con armado. Métodos principales:

| Método | Descripción |
|--------|-------------|
| `__init__` | Recibe geometría, materiales y matriz de barras. Calcula coordenadas de cada barra. |
| `ag()` | Área bruta = b × h |
| `ast()` | Área total de acero = suma de todas las barras |
| `rho()` | Cuantía = Ast / Ag |
| `p0()` | Compresión pura = 0.85·f'c·(Ag − Ast) + fy·Ast |
| `pn_max()` | Máximo nominal = 0.80 × P0 |
| `calcular_punto(c)` | Para una profundidad de eje neutro `c`, calcula (Pn, Mn, φ, εt) |
| `curva_interaccion()` | Barre 40 valores de `c` y genera toda la curva P–M |
| `contorno_aci(pu)` | Genera el contorno biaxial para la carga axial Pu |
| `verificar_espaciamiento()` | Revisa que la separación libre entre barras ≥ 2.5 cm |

### 3.4 Cálculo de la curva P–M (uniaxial)

Para cada profundidad de eje neutro `c`:

1. **Deformación en el concreto**: εcu = 0.003 (máxima)
2. **Bloque de Whitney**: a = β₁ · c, esfuerzo uniforme 0.85·f'c
3. **Deformación en el acero**: εs = 0.003 · (c − d) / c (por secciones planas)
4. **Esfuerzo en el acero**: elastoplástico perfecto (lineal hasta fy, luego constante)
5. **Equilibrio**: Σ fuerzas = Pn, Σ momentos = Mn
6. **Factor φ**: según la deformación neta de tracción εt:
   - εt ≤ εy → φ = 0.65 (compresión controlada)
   - εt ≥ 0.005 → φ = 0.90 (tracción controlada)
   - Intermedio → interpolación lineal

### 3.5 Contorno biaxial

Cuando hay momentos en ambos ejes (Mx y My), se usa el método del **contorno de interacción ACI**:

1. Para cada ángulo θ (24 direcciones), se calcula el par (φ·Mnx, φ·Mny) donde φ·Pn = Pu
2. Se busca en 60 valores de profundidad de eje neutro (escala logarítmica)
3. Se interpola linealmente entre los dos puntos que cruzan Pu
4. El contorno resultante tiene 24 puntos que definen la región resistente
5. Si la demanda (Mux, Muy) cae dentro del contorno → cumple

### 3.6 Verificación ACI (`verification.py`)

Módulo independiente con 8 funciones de verificación:

| Función | Fórmula ACI | Verifica |
|---------|-------------|----------|
| `geometria_columna` | Ag, bc, Ac | Dimensiones básicas |
| `acero_longitudinal` | n_total = n_B·2 + (n_H−2)·2 | N° barras y área total |
| `verificar_cuantia_minima` | ρ ≥ 1% | Cuantía mínima |
| `acero_transversal` | Ash = max(Ash1, Ash2) | Estribos requeridos |
| `longitud_zona_protegida` | Lo = max(45 cm, b, L/6) | Zona de confinamiento |
| `numero_ramales` | n = Ash / A_estribo | Ramales necesarios |
| `separacion_estribos` | Dentro Lo: s ≤ 6d·min(10); Fuera: s ≤ 6d·min(15) | Separación máxima |
| `separacion_acero_longitudinal` | (10·bc − 2φ_est − 2φ_esq − (n−2)·φ_long) / (10·(n−1)) | Separación entre barras |

---

## 4. Conexión con la plantilla HTML (`templates/report.html`)

En el modo WSGI, los resultados se renderizan en HTML usando la plantilla `templates/report.html`.

### Flujo:

```
wsgi_app() 
  → calcula todo (igual que Streamlit)
  → fig_to_b64(fig) convierte cada gráfico a una cadena base64
  → _verif_to_html(vr) formatea la verificación ACI como HTML
  → matriz_to_html(mat) convierte la matriz de barras a tabla HTML
  → Template('report.html').substitute({...}) reemplaza variables
  → devuelve página HTML completa
```

### Variables reemplazadas en `report.html`:

| Variable | Contenido |
|----------|-----------|
| `{{b_cm}}`, `{{h_cm}}` | Dimensiones de la columna |
| `{{fc_val}}`, `{{fy_val}}` | Materiales |
| `{{img_x}}`, `{{img_y}}` | Diagramas P–Mx y P–My (base64) |
| `{{img_sec}}` | Sección transversal (base64) |
| `{{ag_val}}`, `{{ast_val}}`, `{{rho_val}}` | Propiedades |
| `{{p0_val}}`, `{{pnmax_val}}`, `{{phipmax_val}}` | Capacidades |
| `{{biaxial_html}}` | Verificación biaxial (HTML completo) |
| `{{verificacion_html}}` | Verificación ACI (HTML completo) |
| `{{matrix_html}}` | Matriz de barras (tabla HTML) |
| `{{espaciado_html}}` | Errores de espaciamiento |

Las imágenes se incrustan directamente en el HTML como `data:image/png;base64,...`, por lo que el reporte es un **archivo HTML autónomo** sin dependencias externas.

---

## 5. Railway — Plataforma de despliegue

### ¿Qué es Railway?

Railway es una plataforma **Platform as a Service (PaaS)** que permite desplegar aplicaciones web en la nube sin necesidad de configurar servidores. Similar a Heroku, pero más moderno y con facturación por uso.

### Características principales:

- **Despliegue desde GitHub**: cada `git push` dispara un build y deploy automático
- **Escalado automático**: Railway maneja la memoria, CPU y tráfico
- **Domino público**: cada proyecto recibe una URL `*.railway.app`
- **Buildpacks**: detecta automáticamente el lenguaje (Python, Node.js, etc.)
- **Variables de entorno**: Railway asigna `$PORT` dinámicamente

### Cómo funciona el despliegue

1. El repositorio se conecta a Railway (via `railway.json` y `Procfile`)
2. Railway lee `requirements.txt` y ejecuta `pip install`
3. Railway ejecuta: `streamlit run app.py --server.port $PORT --server.address 0.0.0.0`
4. La aplicación queda disponible en la URL pública

### Archivos de configuración

**`railway.json`**:
```json
{
  "build": { "builder": "NIXPACKS" },
  "deploy": {
    "startCommand": "streamlit run app.py --server.port $PORT --server.address 0.0.0.0",
    "healthcheckPath": "/health",
    "restartPolicyType": "ON_FAILURE"
  }
}
```

- `healthcheckPath`: Railway revisa `GET /health` para saber si la app está viva
- `startCommand`: comando exacto para iniciar la app

**`Procfile`** (respaldo):
```
web: streamlit run app.py --server.port $PORT --server.address 0.0.0.0
```

**`.streamlit/config.toml`**:
```toml
[theme]
primaryColor = "#0ea5e9"
backgroundColor = "#ffffff"
secondaryBackgroundColor = "#f1f5f9"
textColor = "#0f172a"
font = "sans serif"

[server]
headless = true
enableCORS = false
enableXsrfProtection = true
```

- `headless = true`: necesario para Railway (no hay navegador)
- `enableCORS = false`: evita problemas de seguridad en producción

---

## 6. Dependencias (`requirements.txt`)

| Librería | Propósito |
|----------|-----------|
| `streamlit` | Interfaz web interactiva |
| `numpy` | Cálculos numéricos y matrices |
| `matplotlib` | Gráficos (diagramas P–M, contorno, sección) |
| `pandas` | Usada internamente por Streamlit (no se importa explícitamente) |
| `Pillow` | Procesamiento de imágenes PNG |

---

## 7. Gestión de memoria

Para evitar consumir memoria innecesaria:

- **matplotlib se carga bajo demanda** con la función `_plt()` — solo se importa cuando se necesita generar un gráfico
- **matplotlib.use('Agg')** solo se activa dentro de `wsgi_app()`, no en el módulo global, para no interferir con Streamlit
- **Las figuras se cierran** después de renderizarse (`plt.close(fig)`)
- **gc.collect()** se llama al final de cada petición WSGI para liberar memoria
- **PIL.MAX_IMAGE_PIXELS = None** evita el error "decompression bomb" en imágenes grandes

---

## 8. Unidades

| Sistema | Longitud | Fuerza | Esfuerzo |
|---------|----------|--------|----------|
| Interno | metro (m) | kgf | kgf/m² |
| Entrada | cm | tonf | kgf/cm² |
| Gráficos | cm | tonf | — |

Constantes en `app.py`:

```python
m = 1.0        # 1 metro
cm = 0.01 * m  # 1 centímetro
mm = 0.001 * m # 1 milímetro
kgf = 1.0      # 1 kilogramo-fuerza
ton = 1000.0 * kgf  # 1 tonelada
```
