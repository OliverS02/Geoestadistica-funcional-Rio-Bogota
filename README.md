# Análisis Funcional Espacial del Índice de Calidad del Agua en el Río Bogotá

Este repositorio contiene el código y los datos asociados al proyecto de investigación sobre el **Índice de Calidad del Agua (ICA)** en el Río Bogotá, utilizando herramientas de análisis funcional y geoestadística en redes hidrográficas.

## Metodología

El proyecto combina técnicas de **Functional Data Analysis (FDA)** y **geoestadística de redes**:

- **Descomposición de Karhunen–Loève (KL)** para reconstruir curvas funcionales a partir de observaciones dispersas.
- **Distancia de Mahalanobis funcional** para medir desviaciones multivariadas en los scores de ICA.
- **Cartas de control T² y Q** para monitoreo estadístico funcional con correlación espacial.
- **Kriging funcional sobre la red hidrográfica** para predicciones espaciotemporales, considerando la topología del río y la conectividad entre afluentes.

## Estructura del repositorio

- `R/` : Funciones principales para análisis funcional y espacial.
- `scripts/` : Ejemplos de análisis, ajuste de curvas, visualización y predicciones.
- `data/` : Datos de ejemplo del ICA y de la red hidrográfica del Río Bogotá.
- `figures/` : Gráficos generados durante el análisis (matriz de covarianza, curvas medias, mapas de predicción).

## Requisitos

- R (>= 4.2)
- Paquetes: `funData`, `fdapace`, `SSN2`, `SSNbler`, `ggplot2`, `sf`, `viridis`, `dplyr`, `lubridate`, `zoo`

## Cómo usar

1. Clonar el repositorio:
    ```bash
    git clone https://github.com/tu_usuario/nombre_repositorio.git
    ```
2. Abrir R y ejecutar los scripts de `scripts/` siguiendo el orden de análisis.
3. Ajustar rutas de archivos si es necesario.

## Resultados

Los scripts permiten:

- Ajustar y visualizar curvas funcionales del ICA para múltiples sitios.
- Estimar la función media ponderada espacialmente y la matriz de covarianza entre ubicaciones.
- Realizar predicciones espaciotemporales del ICA en fechas representativas a lo largo de la red hidrográfica.
- Evaluar la incertidumbre de las predicciones mediante la varianza del error de predicción.

## Licencia

Este proyecto está bajo licencia MIT. Consulta el archivo `LICENSE` para más detalles.

## Referencia

Si utilizas este código, cita el repositorio:

Oliver Danilo Sierra Suarez, 2025. *Análisis Funcional Espacial del Índice de Calidad del Agua en el Río Bogotá*. Repositorio de GitHub: [https://github.com/tu_usuario/nombre_repositorio](https://github.com/tu_usuario/nombre_repositorio)
