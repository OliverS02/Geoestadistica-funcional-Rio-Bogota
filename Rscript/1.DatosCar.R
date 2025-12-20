###############################################################################
# Script: Lectura y preprocesamiento de datos de calidad del agua
#
# Descripción:
# Este script realiza la lectura de datos consolidados de calidad del agua,
# ajusta los tipos de datos, calcula el Índice de Calidad del Agua (ICA)
# a partir de cinco variables fisicoquímicas y agrupa espacialmente las
# observaciones mediante k-means con el fin de incrementar el número de
# observaciones temporales por ubicación para análisis funcional posterior.
#
# El script está pensado para ejecutarse en un entorno local. Las rutas
# deben ajustarse según la ubicación de los archivos en cada sistema.
#
# Entradas:
# - DatosCarConsolidados.csv : archivo CSV con variables fisicoquímicas y
#   coordenadas espaciales.
#
# Salidas:
# - Objeto `datosFinales`: data.frame con ICA calculado y coordenadas
#   agrupadas mediante k-means.
###############################################################################

## ---------------------------------------------------------------------------
## Configuración del entorno (ruta local)
## ---------------------------------------------------------------------------

# Ruta local de trabajo (ajustar según el sistema del usuario)
setwd("C:/Users/sierr/OneDrive/Documentos/Calidad del agua")

## ---------------------------------------------------------------------------
## Librerías
## ---------------------------------------------------------------------------

library(sf)
library(ggplot2)
library(dplyr)
library(readr)

## ---------------------------------------------------------------------------
## Lectura de datos
## ---------------------------------------------------------------------------

# Datos consolidados de calidad del agua
datosCAR <- read.csv("DatosCarConsolidados.csv")

## ---------------------------------------------------------------------------
## Análisis descriptivo inicial
## ---------------------------------------------------------------------------

# Visualización básica de la distribución espacial de los puntos de muestreo
plot(datosCAR$Este.X., datosCAR$Norte..Y.)

# Inspección visual del conjunto de datos
View(datosCAR)

## ---------------------------------------------------------------------------
## Ajuste de tipos de datos
## ---------------------------------------------------------------------------

# Conversión de variables a formato numérico:
# - Reemplazo de separadores decimales
# - Conversión de valores faltantes codificados como texto
for (j in 2:dim(datosCAR)[2]) {
  datosCAR[[j]] <- as.numeric(
    gsub(",", ".", gsub("#N/A", NA, as.character(datosCAR[[j]])))
  )
}

## ---------------------------------------------------------------------------
## Cálculo del Índice de Calidad del Agua (ICA)
## ---------------------------------------------------------------------------

# El ICA se calcula como una combinación ponderada de cinco subíndices:
# ICA = 0.2 * (I_OD + I_SST + I_DQO + I_CE + I_PH)

## ---- I_OD: Oxígeno disuelto -----------------------------------------------

# Función para la concentración de saturación de oxígeno disuelto
# Dependiente de temperatura y altitud (Benson & Krause)
Cp <- function(Temp, A) {
  (14.652 - 0.41022 * Temp + 0.007991 * Temp^2 - 0.000077774 * Temp^3) /
    (1 + 0.000022 * A)
}

# Se asume una altitud promedio de Bogotá de 2640 m
PS_OD <- datosCAR$Oxígeno.disuelto / Cp(datosCAR$Temperatura.agua, 2640)

# Subíndice de oxígeno disuelto
I_OD <- ifelse(PS_OD > 1, 2 - PS_OD, PS_OD)

## ---- I_SST: Sólidos suspendidos totales ------------------------------------

I_SST <- ifelse(
  datosCAR$Sólidos.suspendidos.totales >= 6.67,
  ifelse(
    datosCAR$Sólidos.suspendidos.totales >= 340,
    0,
    1.02 - 0.003 * datosCAR$Sólidos.suspendidos.totales
  ),
  1
)

## ---- I_DQO: Demanda química de oxígeno -------------------------------------

I_DQO <- case_when(
  datosCAR$DQO <= 20                    ~ 0.91,
  datosCAR$DQO > 20 & datosCAR$DQO <= 25 ~ 0.71,
  datosCAR$DQO > 25 & datosCAR$DQO <= 40 ~ 0.51,
  datosCAR$DQO > 40 & datosCAR$DQO <= 80 ~ 0.26,
  datosCAR$DQO > 80                     ~ 0.71
)

## ---- I_CE: Conductividad eléctrica -----------------------------------------

I_CE <- (1 - 10^(-3.26 + 1.34 * log10(datosCAR$Conductividad)))
I_CE <- ifelse(I_CE < 0, 0, I_CE)

## ---- I_PH: Potencial de hidrógeno ------------------------------------------

I_PH <- case_when(
  datosCAR$pH < 4                      ~ 0.1,
  datosCAR$pH >= 4 & datosCAR$pH <= 7  ~ 0.02628419 * exp(datosCAR$pH * 0.520025),
  datosCAR$pH > 7 & datosCAR$pH <= 8   ~ 1,
  datosCAR$pH > 8 & datosCAR$pH <= 11  ~ exp((8 - datosCAR$pH) * 0.5187742),
  datosCAR$DQO > 80                    ~ 0.71
)

## ---- Construcción del ICA --------------------------------------------------

# Matriz de subíndices
indices <- cbind(I_OD, I_SST, I_DQO, I_CE, I_PH)

# Cálculo del índice de calidad del agua
ICA <- 0.2 * indices %*% cbind(rep(1, dim(indices)[2]))

# Adición del ICA al conjunto de datos

datosCAR$ICA<-ICA

datosCAR <- datosCAR %>%
  filter(!is.na(ICA))


## ---------------------------------------------------------------------------
## Exploración del número de observaciones por sitio
## ---------------------------------------------------------------------------

obs_por_sitio <- datosCAR %>%
  group_by(Este.X., Norte..Y.) %>%
  count()

obs_por_sitio

## ---------------------------------------------------------------------------
## Agrupación espacial mediante k-means
## ---------------------------------------------------------------------------

# Matriz de coordenadas espaciales
coords <- as.matrix(datosCAR %>% select(Este.X., Norte..Y.))

# Evaluación del número de clusters mediante WSS
size <- seq(10, 200, 10)
wss <- numeric(20)

for (k in 1:20) {
  km <- kmeans(coords, centers = size[k])
  wss[k] <- km$tot.withinss
}

# Visualización del criterio WSS

plot(size, wss, type = "b",
     xlab = "Número de clusters k",
     ylab = "WSS")

# Selección de k = 50 clusters (zona de estabilización del WSS)
km <- kmeans(coords, centers = 40)

# Coordenadas representativas de cada cluster
new_coords <- km$centers

# Asignación de nuevas coordenadas a cada observación
datosCAR$X <- new_coords[, 1][km$cluster]
datosCAR$Y <- new_coords[, 2][km$cluster]




