############################################################
# Kriging funcional del ICA sobre la red del Río Bogotá
############################################################

## ---------------------------------------------------------------------------
## 1) Configuración del entorno
## ---------------------------------------------------------------------------
setwd("C:/Users/sierr/OneDrive/Documentos/Calidad del agua")

library(SSN2)
library(SSNbler)
library(funData)
library(ggplot2)
library(sf)
library(viridis)

# Cargar modelos ajustados
modelos_ajustados <- readRDS(file = 'modelos_ssn_rio_bogota.rds')

############################################################
# 2) Predicciones de scores para cada componente
############################################################
e1_hat <- predict(modelos_ajustados$e1$modelo, newdata = 'preds', interval = 'prediction')
e2_hat <- predict(modelos_ajustados$e2$modelo, newdata = 'preds', interval = 'prediction')
e3_hat <- predict(modelos_ajustados$e3$modelo, newdata = 'preds', interval = 'prediction')
e4_hat <- predict(modelos_ajustados$e4$modelo, newdata = 'preds', interval = 'prediction')

# Varianza total del error de predicción (suma de los 4 scores)
var_ep <- (e1_hat[,3]-e1_hat[,2])^2/(4*qnorm(0.975)^2) +
  (e2_hat[,3]-e2_hat[,2])^2/(4*qnorm(0.975)^2) +
  (e3_hat[,3]-e3_hat[,2])^2/(4*qnorm(0.975)^2) +
  (e4_hat[,3]-e4_hat[,2])^2/(4*qnorm(0.975)^2)

# Matriz de predicciones de scores
E_hat <- cbind(e1_hat[,1], e2_hat[,1], e3_hat[,1], e4_hat[,1])

# Reconstruir curvas funcionales con KL
pred_ICA <- funData(
  argvals = curvas@argvals,
  X       = KL_ICA$mu + E_hat %*% t(KL_ICA$phi)
)

############################################################
# 3) Selección de fechas de interés
############################################################
fechas_objetivo <- as.Date(c("2017-12-15", "2018-03-15", "2019-06-15"))
fechas_num <- as.numeric(fechas_objetivo - fecha0)

# Encontrar índices más cercanos en KL_ICA$workGrid
indices_fechas <- sapply(fechas_num, function(fn){
  which.min(abs(KL_ICA$workGrid - fn))
})
names(indices_fechas) <- c("diciembre de 2017", "marzo de 2018", "junio de 2019")

############################################################
# 4) Función para generar mapas de predicción
############################################################
mapa_prediccion <- function(idx, score_label = "ICA") {
  ggplot(data = NULL) +
    # Departamentos
    geom_sf(data = cundi_utm, fill = "forestgreen", alpha = 0.4, color = NA) +
    # Red de ríos
    geom_sf(data = red_hidri_bogota_snn$edges, color = "lightblue", size = 0.6) +
    # Puntos de predicción coloreados por score
    geom_sf(
      data = red_hidri_bogota_snn$preds$preds,
      aes(color = pred_ICA@X[, idx]),
      size = 0.3,
      alpha = 0.8
    ) +
    # Puntos observados en blanco
    geom_sf(
      data = red_hidri_bogota_snn$obs,
      color = "white",
      size = 0.3,
      alpha = 0.9
    ) +
    # Escala de color continua
    scale_color_viridis_c(option = "viridis") +
    coord_sf(datum = st_crs(red_hidri_bogota_snn$edges)) +
    labs(
      color    = score_label,
      title    = paste0("Predicciones de ", score_label, " sobre la red del Río Bogotá"),
      subtitle = paste("Fecha:", names(idx))
    ) +
    theme_minimal() +
    theme(
      legend.text  = element_text(size = 8),
      legend.title = element_text(size = 10)
    )
}

# =========================
# 5) Generar mapas para las fechas
# =========================
map_dic2017 <- mapa_prediccion(indices_fechas["diciembre de 2017"])
map_mar2018 <- mapa_prediccion(indices_fechas["marzo de 2018"])
map_jun2019 <- mapa_prediccion(indices_fechas["junio de 2019"])

# Mostrar mapas
map_dic2017
map_mar2018
map_jun2019

############################################################
# 6) Mapa de varianzas del error de predicción
############################################################
ggplot(data = NULL) +
  geom_sf(data = cundi_utm, fill = "forestgreen", alpha = 0.4, color = NA) +
  geom_sf(data = red_hidri_bogota_snn$edges, color = "lightblue", size = 0.6) +
  geom_sf(
    data = red_hidri_bogota_snn$preds$preds,
    aes(color = var_ep),
    size = 0.3,
    alpha = 0.8
  ) +
  geom_sf(
    data = red_hidri_bogota_snn$obs,
    color = "white",
    size = 0.3,
    alpha = 0.9
  ) +
  scale_color_viridis_c(option = "plasma") +
  coord_sf(datum = st_crs(red_hidri_bogota_snn$edges)) +
  labs(
    color = "Varianza del error",
    title = "Varianza de los errores de predicción en la red del Río Bogotá"
  ) +
  theme_minimal() +
  theme(
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 10)
  )
