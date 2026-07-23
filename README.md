# Programa de columnas ACI 318-19

Aplicacion Streamlit para revision de columnas rectangulares de concreto armado a flexocompresion uniaxial y biaxial por contorno de carga ACI 318-19.

## Ejecucion local

```powershell
python -m pip install -r requirements.txt
streamlit run app.py
```

## Despliegue en Streamlit Community Cloud

1. Cree un repositorio GitHub y suba `app.py`, `requirements.txt` y `.streamlit/config.toml`.
2. En [Streamlit Community Cloud](https://share.streamlit.io/), seleccione el repositorio, rama y el archivo principal `app.py`.
3. El servicio instalara automaticamente las dependencias de `requirements.txt`.

> Nota de ingenieria: el resultado es una ayuda de diseno. Antes de emitir planos, valide cargas, esbeltez, segundo orden, detallado, confinamiento y disposiciones sismicas aplicables del ACI 318-19.
