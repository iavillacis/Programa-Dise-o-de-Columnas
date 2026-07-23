# Programa de columnas ACI 318-19

Aplicación Streamlit para revisar columnas rectangulares de concreto armado a flexocompresión uniaxial, biaxial por carga recíproca de Bresler y cortante.

## Ejecución local

```powershell
python -m pip install -r requirements.txt
streamlit run app.py
```

## Despliegue en Streamlit Community Cloud

1. Cree un repositorio GitHub y suba `app.py`, `requirements.txt` y `.streamlit/config.toml`.
2. En [Streamlit Community Cloud](https://share.streamlit.io/), seleccione el repositorio, rama y el archivo principal `app.py`.
3. El servicio instalará automáticamente las dependencias de `requirements.txt`.

## Despliegue en Vercel

El archivo `app.py` también expone una variable WSGI superior llamada `app`, que Vercel detecta como una Python Function. El despliegue en Vercel ofrece una calculadora serverless de Bresler y cortante con el mismo núcleo de cálculo. La versión Streamlit conserva el editor visual de barras y los diagramas P-M, por lo que sigue siendo la opción recomendada para el análisis gráfico completo.

1. Importe el repositorio en Vercel sin definir un *Build Command* ni un *Output Directory*.
2. Vercel instalará `requirements.txt` y detectará automáticamente `app` como entrada WSGI.
3. Compruebe la función con la ruta `/health` después del despliegue.

> Nota de ingeniería: el resultado es una ayuda de diseño. Antes de emitir planos, valide cargas, esbeltez, segundo orden, detallado, confinamiento y disposiciones sísmicas aplicables del ACI 318-19.
