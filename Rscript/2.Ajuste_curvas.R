############################################################
# Ajuste de curvas funcionales con Sparse FDA
#
# Descripción:
#   Este script ajusta curvas funcionales del Índice de Calidad
#   del Agua (ICA) para múltiples ubicaciones espaciales bajo
#   un esquema de muestreo irregular en el tiempo.
#
#   Se utiliza Sparse Functional Data Analysis (FPCA) para
#   reconstruir trayectorias temporales suavizadas a partir de
#   observaciones escasas y no uniformes en cada sitio.
#
#   Las curvas se visualizan únicamente en intervalos temporales
#   donde existe información observada en al menos una ubicación
#   espacial. Los gaps globales (más de 4 meses sin datos) no
#   se muestran.
#
# Entradas:
#   - datosCAR: data.frame con columnas
#       * X, Y                  : coordenadas espaciales
#       * Fecha..mes.y.año.      : fecha de observación (Date)
#       * ICA                   : índice de calidad del agua
#
# Salidas:
#   - Objeto FPCA (KL_ICA)
#   - Curvas funcionales reconstruidas (curvas)
############################################################

# Librerías
library(funData)
library(fdapace)
library(dplyr)
library(lubridate)
library(zoo)

############################################################
# Preparación y ordenamiento de los datos
############################################################


datosFinales <- datosCAR %>%
  select(X, Y, Fecha..mes.y.año., ICA) %>%
  # Convertir la fecha a Date (día/mes/año)
  mutate(
    Fecha..mes.y.año. = as.Date(Fecha..mes.y.año., format = "%d/%m/%Y"),
    ICA = as.numeric(ICA)
  ) %>%
  arrange(X, Y, Fecha..mes.y.año.)

############################################################
# Agregación temporal por sitio
############################################################

datosAgregados <- datosFinales %>%
  group_by(X, Y, Fecha..mes.y.año.) %>%
  summarise(median_ICA = median(ICA), .groups = "drop")

############################################################
# Escala temporal numérica
############################################################

fecha0 <- min(datosAgregados$Fecha..mes.y.año.)
datosAgregados <- datosAgregados %>%
  mutate(time_num = as.numeric(Fecha..mes.y.año. - fecha0))

############################################################
# Identificación de réplicas espaciales
############################################################

n <- datosAgregados %>%
  group_by(X, Y) %>%
  count() %>%
  pull(n)


datosAgregados$site_id <- rep(seq_along(n), n)

# Función para construir índices de inicio y fin por réplica
index_matrix <- function(n) {
  cbind(
    start = c(1, cumsum(n)[-length(n)] + 1),
    end   = cumsum(n)
  )
}

index <- index_matrix(n)

############################################################
# VISUALIZACIÓN DE LA IRREGULARIDAD DEL MUESTREO
############################################################

# Crear un identificador de sitio si no existe
if (!"site_id" %in% names(datosFinales)) {
  n_sites <- datosFinales %>%
    group_by(X, Y) %>%
    count() %>%
    pull(n)
  datosFinales$site_id <- rep(seq_along(n_sites), n_sites)
}

par(mar = c(6, 4, 2, 1))

fechas <- datosFinales$Fecha..mes.y.año.

# Colores por sitio
sites <- unique(datosFinales$site_id)
col_sites <- setNames(
  viridisLite::viridis(length(sites), option = "turbo"),
  sites
)

cols <- col_sites[as.character(datosFinales$site_id)]

plot(
  fechas,
  datosFinales$site_id,
  xaxt = "n",
  xlab = "Fecha",
  ylab = "Punto de muestreo",
  pch  = 19,
  cex  = 0.7,
  col  = cols
)

axis.Date(
  side = 1,
  at = seq(min(fechas), max(fechas), by = "6 months"),
  labels = FALSE
)

text(
  x = seq(min(fechas), max(fechas), by = "6 months"),
  y = par("usr")[3] - 0.05 * diff(par("usr")[3:4]),
  labels = format(seq(min(fechas), max(fechas), by = "6 months"), "%Y-%m"),
  srt = 45,
  adj = 1,
  xpd = TRUE,
  cex = 0.8
)

box()




############################################################
# Construcción de objetos para Sparse FPCA
############################################################

ICA_list   <- vector("list", length(n))
times_list <- vector("list", length(n))

for (j in seq_along(n)) {
  ICA_list[[j]]   <- datosAgregados$median_ICA[index[j, 1]:index[j, 2]]
  times_list[[j]] <- datosAgregados$time_num[index[j, 1]:index[j, 2]]
}

############################################################
# Sparse FPCA
############################################################

KL_ICA <- FPCA(
  Ly = ICA_list,
  Lt = times_list
)

############################################################
# Reconstrucción de curvas funcionales
############################################################

curvas <- funData(
  argvals = list(KL_ICA$workGrid),
  X = KL_ICA$mu + KL_ICA$xiEst %*% t(KL_ICA$phi)
)


############################################################
# Identificación de gaps globales (>6 meses sin datos)
############################################################

# Ordenar todas las fechas con observaciones
fechas_globales <- sort(unique(datosAgregados$Fecha..mes.y.año.))
# Diferencias entre fechas consecutivas
diffs <- diff(fechas_globales)

# Umbral de gap: 120 días
umbral_gap <- 180


# Detectar índices donde la diferencia excede el umbral
gap_idx <- which(diffs > umbral_gap)

# Construir tabla de intervalos
intervalos_gaps <- tibble(
  inicio = fechas_globales[gap_idx] + 1,
  fin    = fechas_globales[gap_idx + 1] - 1
)

# Convertir a escala numérica de días desde fecha0
fecha0 <- min(datosAgregados$Fecha..mes.y.año.)
intervalos_gaps <- intervalos_gaps %>%
  mutate(
    inicio_num = as.numeric(inicio - fecha0),
    fin_num    = as.numeric(fin - fecha0)
  )

############################################################
# GRÁFICO DE CURVAS FUNCIONALES CON FRANJAS ICA Y GAPS
############################################################

# ==============================================
# CONFIGURACIÓN DEL PLOT
# ==============================================
par(mar = c(7, 4, 2, 1))

# Extraer datos del funData
X  <- curvas@X
tt <- curvas@argvals[[1]]

ylim <- c(0,1)

# ==============================================
# PLOT BASE (vacío)
# ==============================================
plot(
  NA,
  xlim = range(tt),
  ylim = ylim,
  xlab = "Fecha",
  ylab = "ICA",
  xaxt = "n"
)

# ==============================================
# FRANJAS ICA (IDEAM)
# ==============================================
rect(min(tt), 0.00, max(tt), 0.25,
     col = rgb(0.6, 0.0, 0.0, 0.25), border = NA) # Muy Mala
rect(min(tt), 0.25, max(tt), 0.50,
     col = rgb(1.0, 0.4, 0.0, 0.25), border = NA) # Mala
rect(min(tt), 0.50, max(tt), 0.70,
     col = rgb(1.0, 0.8, 0.0, 0.25), border = NA) # Media
rect(min(tt), 0.70, max(tt), 0.90,
     col = rgb(0.4, 0.8, 0.0, 0.25), border = NA) # Buena
rect(min(tt), 0.90, max(tt), 1.00,
     col = rgb(0.0, 0.6, 0.4, 0.25), border = NA) # Excelente

# ==============================================
# CURVAS AJUSTADAS SEGÚN GAP
# ==============================================

# Crear vector lógico de validez para cada punto de tt
is_valid <- rep(TRUE, length(tt))
if (!is.null(intervalos_gaps) && nrow(intervalos_gaps) > 0) {
  for (i in 1:nrow(intervalos_gaps)) {
    # Puntos dentro del gap se marcan como no válidos
    is_valid[tt >= intervalos_gaps$inicio_num[i] & tt < intervalos_gaps$fin_num[i]] <- FALSE
  }
}

# Dibujar curvas por tramos consecutivos
for (i in 1:nrow(X)) {
  r <- rle(is_valid)
  idx_end <- cumsum(r$lengths)
  idx_start <- c(1, head(idx_end + 1, -1))
  
  for (j in seq_along(idx_start)) {
    # Tramos válidos (con datos)
    if (r$values[j]) {
      lines(tt[idx_start[j]:idx_end[j]], X[i, idx_start[j]:idx_end[j]],
            col = "gray40", lty = 1, lwd = 1.2)
    } else {
      # Tramos dentro de gaps
      lines(tt[idx_start[j]:idx_end[j]], X[i, idx_start[j]:idx_end[j]],
            col = rgb(0.7, 0.7, 0.7, 0.2), lty = 1, lwd = 1)
    }
  }
}

# ==============================================
# BARRA INFERIOR: DATOS DISPONIBLES Y GAPS
# ==============================================
y_bottom <- par("usr")[3]        # borde inferior del plot
height   <- 0.02 * diff(par("usr")[3:4])  # altura de la barra

# Dibujar toda la barra de verde claro (intervalos con datos)
rect(min(tt), y_bottom, max(tt), y_bottom + height,
     col = rgb(0.0, 0.6, 0.0, 0.5), border = NA)

# Dibujar gaps en rojo claro (intervalos sin datos)
if (!is.null(intervalos_gaps) && nrow(intervalos_gaps) > 0) {
  for (i in 1:nrow(intervalos_gaps)) {
    rect(
      xleft   = intervalos_gaps$inicio_num[i],
      xright  = intervalos_gaps$fin_num[i],
      ybottom = y_bottom,
      ytop    = y_bottom + height,
      col     = rgb(1, 0, 0, 0.5),  # rojo claro semitransparente
      border  = NA
    )
  }
}

# ==============================================
# LÍNEAS VERTICALES VERDES EN EXTREMOS DE INTERVALOS CON DATOS
# ==============================================
# Calcular intervalos con datos: espacios entre los gaps
info_starts <- c(min(tt), intervalos_gaps$fin_num + 1)
info_ends   <- c(intervalos_gaps$inicio_num - 1, max(tt))

for (i in seq_along(info_starts)) {
  abline(v = info_starts[i], col = "darkgreen", lty = 2, lwd = 1.5) # punteada
  abline(v = info_ends[i],   col = "darkgreen", lty = 2, lwd = 1.5)
}

# ==============================================
# EJE X EN FECHAS
# ==============================================
fecha0 <- min(datosAgregados$Fecha..mes.y.año.)

fechas_ticks <- seq(
  from = min(datosAgregados$Fecha..mes.y.año.),
  to   = max(datosAgregados$Fecha..mes.y.año.),
  by   = "3 months"
)

ticks <- as.numeric(fechas_ticks - fecha0)

axis(1, at = ticks, labels = FALSE)

text(
  x = ticks,
  y = par("usr")[3] - 0.08 * diff(par("usr")[3:4]),
  labels = format(fechas_ticks, "%Y-%m"),
  srt = 45,
  adj = 1,
  xpd = TRUE,
  cex = 0.75
)


