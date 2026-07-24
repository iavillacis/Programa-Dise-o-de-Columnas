# Diagrama de Flujo — Programa de Diseño de Columnas ACI 318-19

```mermaid
flowchart TD
    A["INICIO"] --> B["Usuario ingresa datos: 
    B, H, rec, estribo, L, s_estribo
    f'c, fy
    n_B, n_H, diam_long, diam_corner
    Pu, Mux, Muy"]
    B --> C["Presiona CALCULAR"]
    C --> D{"¿Datos válidos?
    (B > 0, H > 0, 
     rec < B/2, rec < H/2,
     fc > 0, fy > 0)"}
    D -->|No| E["Mostrar error"]
    E --> C
    D -->|Sí| F["Generar matriz de barras
    (solo perímetro: 
     esquinas = diam_corner,
     bordes = diam_long,
     centro = 0)"]
    F --> G["Crear objeto DiagramaInteraccion
    (columna con su armado)"]
    G --> H["Verificar espaciamiento
    entre barras ≥ 2.5 cm"]
    H --> I{"¿Cuantía ρ = Ast/Ag?"}
    I -->|"ρ > 8%"| J["ERROR: excede máximo ACI"]
    J --> K["DETENER"]
    I -->|"ρ < 1%"| L["ADVERTENCIA: mínimo 1%"]
    L --> M
    I -->|"1% ≤ ρ ≤ 8%"| M["Continuar"]
    K --> M
    M --> N["Convertir demandas:
    Pu, Mux, Muy → unidades internas"]
    
    N --> O["EXPANDER 1:
    Tabla de propiedades
    (Ag, Ast, ρ, P0, Pn,max, φPn,max)
    Matriz de barras"]
    
    O --> P["EXPANDER 2:
    Diagramas P–Mx y P–My
    (curva nominal + φ + φPn,max + demanda)"]
    
    P --> Q["EXPANDER 3:
    Verificación ACI 318-19
    Geometría → Acero long. → Cuantía mínima
    → Ash → Lo → Ramales → Separaciones"]
    
    Q --> R{"¿|Muy| > 0.01?"}
    R -->|Sí| S["EXPANDER 4:
    Contorno biaxial ACI
    φMnx vs φMny a Pu dado
    D/C radial"]
    R -->|No| T["(Saltar biaxial)"]
    S --> U
    T --> U
    
    U --> V["EXPANDER 5:
    Supuestos técnicos
    (εcu, Whitney, φ variable, etc.)"]
    
    V --> W["FIN"]
```

## Flujo alternativo — Modo WSGI (Railway / servidor)

```mermaid
flowchart LR
    A["Petición HTTP
    GET /?b=40&h=40&...
    (parámetros en URL)"] 
    A --> B["wsgi_app()
    Lee QUERY_STRING"]
    B --> C{"¿/health?"}
    C -->|Sí| D["200 OK
    {'status':'ok'}"]
    C -->|No| E["Calcula todo
    (mismo núcleo que Streamlit)"]
    E --> F["Genera HTML con
    imágenes PNG en base64
    + tabla matriz + verificaciones"]
    F --> G["200 OK
    text/html"]
```

## Leyenda de formas

| Forma | Significado |
|-------|-------------|
| Rectángulo | Proceso / acción |
| Rombo | Decisión / condición |
| Rectángulo redondeado | Inicio / Fin |
| Rectángulo con doble línea | Subproceso / expansor |
