# Programa de columnas ACI 318-19

Aplicación Streamlit para revisar columnas rectangulares de concreto armado a flexocompresión uniaxial, biaxial por carga recíproca de Bresler y cortante.

## Ejecución local

```powershell
python -m pip install -r requirements.txt
streamlit run app.py
```

## Despliegue

1. Cree un repositorio GitHub y suba `app.py`, `requirements.txt` y `.streamlit/config.toml`.
2. En [Streamlit Community Cloud](https://share.streamlit.io/), seleccione el repositorio, rama y el archivo principal `app.py`.
3. El servicio instalará automáticamente las dependencias de `requirements.txt`.

Streamlit mantiene una sesión y un servidor Python persistente; por ello Streamlit Community Cloud es el destino compatible recomendado. Vercel Serverless no admite de forma nativa una aplicación Streamlit/WSGI persistente. Para usar un dominio administrado en Vercel, configure un *rewrite* o proxy hacia el despliegue de Streamlit, en vez de intentar ejecutar `app.py` como una función serverless.

> Nota de ingeniería: el resultado es una ayuda de diseño. Antes de emitir planos, valide cargas, esbeltez, segundo orden, detallado, confinamiento y disposiciones sísmicas aplicables del ACI 318-19.
