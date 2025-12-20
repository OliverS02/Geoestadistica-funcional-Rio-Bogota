#############################################################
# Script: Creación de red hidrológica y SSN para Río Bogotá
# Objetivo: Generar LSN y SSN a partir de HydroRIVERS y shapefiles locales
#            para los municipios de Cundinamarca intersectados por el Río Bogotá
############################################################

# ==========================
# Librerías
# ==========================
library(sf)
library(dplyr)
library(SSNbler)   # CRAN
library(ggplot2)
library(tidyr)

# ==========================
# Directorio de trabajo
# ==========================
setwd("C:/Users/sierr/OneDrive/Documentos/Calidad del agua")

# ==========================
# Shapefile de departamentos y municipios
# ==========================
Dept <- st_read("Departamentos_Abril_2025_shp/Departamento.shp") %>%
  st_transform(3116)   # CRS proyectado (metros)

cundi_utm <- Dept %>% filter(DeNombre == "Cundinamarca")

Muni <- st_read("Servicios_P%C3%BAblicos_-_Municipios/Servicios_Públicos_-_Municipios.shp") %>%
  st_transform(3116) %>%
  filter(DEPTO == "CUNDINAMARCA")

# ==========================
# Río Bogotá
# ==========================
rio <- st_read('Bogota/CuerpoAgua.shp') %>% st_transform(3116)
rio_bogota <- rio %>% filter(NOMBRE == "Río Bogotá")

# ==========================
# Municipios que intersectan el río
# ==========================
municipios_con_rio <- Muni %>%
  filter(st_intersects(geometry, rio_bogota, sparse = FALSE) %>% rowSums() > 0)

# ==========================
# Visualización inicial
# ==========================
ggplot(data = NULL) +
  geom_sf(data = cundi_utm, fill = "grey90", color = "white") +
  geom_sf(data = municipios_con_rio, fill = "lightgreen", alpha = 0.4, color = NA) +
  geom_sf(data = rio_bogota, color = "blue", size = 1) +
  theme_minimal() +
  labs(title = "Río Bogotá")

# ==========================
# Cargar HydroRIVERS
# ==========================
ruta_gdb <- "HydroRIVERS_v10_sa.gdb/HydroRIVERS_v10_sa.gdb"
rivers <- st_read(ruta_gdb, layer = "HydroRIVERS_v10_sa") %>%
  st_transform(3116)

# Seleccionar ríos dentro de municipios con río
Cuerpo_Agua <- rivers[st_intersects(rivers, municipios_con_rio, sparse = FALSE) %>% rowSums() > 0, ]

# ==========================
# Preparar edges para SSN
# ==========================
edges <- Cuerpo_Agua %>%
  st_transform(3116) %>%
  mutate(
    seg_id = as.character(HYRIV_ID),
    next_down_seg_id = ifelse(NEXT_DOWN == 0, NA, as.character(NEXT_DOWN))
  ) %>%
  select(seg_id, next_down_seg_id, LENGTH_KM, DIST_DN_KM, everything())

# Asegurar geometrías LINESTRING
if(!all(sf::st_geometry_type(edges) == "LINESTRING")){
  edges <- sf::st_cast(edges, "LINESTRING")
}

# ==========================
# Crear carpeta LSN
# ==========================
lsn_path <- "lsn_cuerpo_agua"
dir.create(lsn_path, showWarnings = FALSE)

# ==========================
# Crear LSN
# ==========================
lines_to_lsn(
  streams = edges,
  lsn_path = lsn_path,
  check_topology = TRUE,
  snap_tolerance = 0,
  topo_tolerance = 0,
  remove_ZM = TRUE,
  overwrite = TRUE
)

# Leer edges generados por lines_to_lsn
edges_lsn <- st_read(file.path(lsn_path, "edges.gpkg"), quiet = TRUE)

# ==========================
# Calcular upstream distance
# ==========================
edges_lsn <- updist_edges(
  edges = edges_lsn,
  lsn_path = lsn_path,
  calc_length = TRUE,
  save_local = TRUE,
  overwrite = TRUE
)

# ==========================
# Resolver confluencias problemáticas (>2 entradas)
# ==========================
rels <- read.csv(file.path(lsn_path, "relationships.csv"))
confluencias <- table(rels$toedge)
problem_nodes <- as.numeric(names(confluencias[confluencias > 2]))

rels_fixed <- rels
for(n in problem_nodes){
  from_edges <- rels$fromedge[rels$toedge == n]
  if(length(from_edges) > 2){
    dis <- edges_lsn$DIS_AV_CMS[match(from_edges, edges_lsn$rid)]
    keep <- from_edges[order(dis, decreasing = TRUE)[1:2]]
    remove <- setdiff(from_edges, keep)
    rels_fixed <- rels_fixed[!(rels_fixed$fromedge %in% remove & rels_fixed$toedge == n), ]
  }
}

write.csv(rels_fixed, file.path(lsn_path, "relationships_fixed.csv"), row.names = FALSE)
file.copy(file.path(lsn_path, "relationships_fixed.csv"), file.path(lsn_path, "relationships.csv"), overwrite = TRUE)

# ==========================
# Ensamblar SSN
# ==========================
ssn_outdir <- "SSN_Cuerpo_Agua.ssn"
ssn <- SSNbler::ssn_assemble(
  edges = edges_lsn,
  lsn_path = lsn_path,
  ssn_path = ssn_outdir,
  overwrite = TRUE
)

# ==========================
# Calcular AFV usando CATCH_SKM
# ==========================
SSNbler::afv_edges(
  edges = edges_lsn,
  lsn_path = lsn_path,
  infl_col = "UPLAND_SKM",
  segpi_col = "segPI",
  afv_col = "afv"
)


# ------------------------------------------------------------
# Añadir puntos observados, variables y ubicaciones a predecir
# ------------------------------------------------------------

# --------------------------------------------------
# 1) Snap de puntos observados a la red (LSN)
# --------------------------------------------------
# puntos_sf: puntos de muestreo en sf
# edges_lsn: edges ya procesados con lines_to_lsn()
# lsn_path: carpeta donde se guardó la LSN

puntos_sf<- datosAgregados %>% 
            select(X,Y) %>% 
            distinct() %>% 
            st_as_sf(coords = c("X", "Y"), crs = 3116)

obs_snapped <- sites_to_lsn(
  sites = puntos_sf,
  edges = edges_lsn,
  lsn_path = lsn_path,
  file_name = "obs",
  snap_tolerance = 50000,  # tolerancia en metros para snap
  save_local = TRUE,
  overwrite = TRUE
)

# --------------------------------------------------
# 2) Generar puntos de predicción a lo largo de los edges
# --------------------------------------------------

# --------------------------------------------------
#  Generar puntos de predicción sobre la red
# --------------------------------------------------
pred_geom <- st_line_sample(
  edges_lsn,
  density = 1 / 1000
)

# Eliminar geometrías vacías si existen
pred_geom <- pred_geom[!st_is_empty(pred_geom)]

pred_sf <- st_cast(pred_geom, "POINT") |>
  st_sf(crs = st_crs(edges_lsn))

# --------------------------------------------------
#  Eliminar duplicados geométricos (FORMA CORRECTA)
# --------------------------------------------------
coords <- st_coordinates(pred_sf)

pred_sf <- pred_sf |>
  mutate(
    x = coords[, 1],
    y = coords[, 2]
  ) |>
  distinct(x, y, .keep_all = TRUE) |>
  select(-x, -y)

# --------------------------------------------------
# 3) Snapear puntos de predicción a la LSN
# --------------------------------------------------
preds <- sites_to_lsn(
  sites = pred_sf,
  edges = edges_lsn,
  lsn_path = lsn_path,
  file_name = "preds",
  snap_tolerance = 200,
  overwrite = TRUE
)

# --------------------------------------------------
# 4) Calcular upstream distance para todos los sitios
# --------------------------------------------------
site.list <- updist_sites(
  sites = list(obs = obs_snapped, preds = preds),
  edges = edges_lsn,
  length_col = "LENGTH_KM",
  save_local = TRUE,
  lsn_path = lsn_path
)

# --------------------------------------------------
# 5) Calcular AFV en los edges
# --------------------------------------------------
# Asegurarse de que la columna "pi" exista

edges_lsn <- afv_edges(
  edges      = edges_lsn,
  lsn_path   = lsn_path,
  infl_col   = "UPLAND_SKM",  # variable de influencia
  segpi_col  = "pi",
  afv_col    = "afv",
  save_local = TRUE
)

# Calcular AFV para los puntos de muestreo
site.list <- afv_sites(
  sites      = site.list,
  edges      = edges_lsn,
  afv_col    = "afv",
  save_local = TRUE,
  lsn_path   = lsn_path
)


# --------------------------------------------------
# 6) Añadir scores de KL (ejemplo)
# --------------------------------------------------
# Asegurarse que KL_ICA$xiEst tenga la misma longitud que los puntos observados
site.list$obs <- site.list$obs %>%
  mutate(
    e_1 = KL_ICA$xiEst[,1],
    e_2 = KL_ICA$xiEst[,2],
    e_3 = KL_ICA$xiEst[,3],
    e_4 = KL_ICA$xiEst[,4]
  )

# --------------------------------------------------
# 7) Ensamblar el SSN
# --------------------------------------------------
ssn_outdir <- file.path(lsn_path, "red_hidri_bogota.ssn")

red_hidri_bogota_snn <- ssn_assemble(
  edges       = edges_lsn,
  lsn_path    = lsn_path,
  obs_sites   = site.list$obs,
  preds_list  = site.list["preds"],
  ssn_path    = ssn_outdir,
  import      = TRUE,
  check       = TRUE,
  afv_col     = "afv",
  overwrite   = TRUE
)

# --------------------------------------------------
# 8) Visualización básica
# --------------------------------------------------
ggplot(data = NULL) +
  geom_sf(data = cundi_utm, fill = "forestgreen", alpha = 0.4) +
  geom_sf(
    data = red_hidri_bogota_snn$edges,
    aes(linewidth = afv),
    color = "lightblue"
  ) +
  scale_linewidth(range = c(0.1, 2.5)) +
  geom_sf(
    data = red_hidri_bogota_snn$preds$preds,
    size = 0.5,
    shape = 21,
    fill = "blue",
    color = "blue"
  ) +
  geom_sf(
    data = red_hidri_bogota_snn$obs,
    size = 1.7,
    aes(color = e_1)  # variable de interés
  ) +
  coord_sf(datum = st_crs(red_hidri_bogota_snn$edges)) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    color = "Score 1",
    linewidth = "Área AFV",
    title = "Upstream Distance y AFV en la red del Río Bogotá"
  ) +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

# Guardar objeto SSN
saveRDS(red_hidri_bogota_snn, file = file.path(lsn_path, "red_hidri_bogota_ssn.rds"))


