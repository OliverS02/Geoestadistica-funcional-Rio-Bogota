############################################################
## MODELADO SSN – RÍO BOGOTÁ
## Selección de estructuras de covarianza y ajuste final
############################################################

library(SSN2)

############################################################
# 1. Crear matrices de distancia hidrológica
#    (requisito previo para ssn_lm)
############################################################
ssn_create_distmat(red_hidri_bogota_snn)

############################################################
# 2. Espacio de búsqueda de estructuras de covarianza
############################################################

# Tail-up
TU <- c("linear", "spherical", "exponential", "mariah", "epa", "none")

# Tail-down
TD <- c("linear", "spherical", "exponential", "mariah", "epa", "none")

# Euclidiana
EU <- c(
  "spherical", "exponential", "gaussian", "cosine", "cubic",
  "pentaspherical", "wave", "jbessel", "gravity", "rquad",
  "magnetic", "none"
)

############################################################
# 3. Inicialización de almacenamiento de resultados
############################################################

best_model1 <- best_model2 <- best_model3 <- best_model4 <- c(NA, NA, NA)

best_ll1 <- best_ll2 <- best_ll3 <- best_ll4 <- -Inf

n_fail <- 0   # número de fallos de ajuste
n_fit  <- 0   # número total de intentos

############################################################
# 4. Función segura para ajuste SSN
#    (evita que errores detengan el grid search)
############################################################

safe_ssn_lm <- function(formula, TU, TD, EU) {
  tryCatch(
    {
      mod <- ssn_lm(
        formula = formula,
        ssn.object = red_hidri_bogota_snn,
        tailup_type = TU,
        taildown_type = TD,
        euclid_type = EU,
        additive = "afv"
      )
      list(ok = TRUE, model = mod)
    },
    error = function(e) {
      list(ok = FALSE, msg = e$message)
    }
  )
}

############################################################
# 5. Grid search sobre estructuras (criterio: logLik)
############################################################

for (model_TU in TU) {
  for (model_TD in TD) {
    for (model_EU in EU) {
      
      n_fit <- n_fit + 1
      
      # Ajustes para cada score funcional
      m1 <- safe_ssn_lm(e_1 ~ 1, model_TU, model_TD, model_EU)
      m2 <- safe_ssn_lm(e_2 ~ 1, model_TU, model_TD, model_EU)
      m3 <- safe_ssn_lm(e_3 ~ 1, model_TU, model_TD, model_EU)
      m4 <- safe_ssn_lm(e_4 ~ 1, model_TU, model_TD, model_EU)
      
      # Score 1
      if (m1$ok) {
        ll1 <- glances(m1$model)$logLik
        if (ll1 > best_ll1) {
          best_ll1 <- ll1
          best_model1 <- c(model_TU, model_TD, model_EU)
        }
      } else n_fail <- n_fail + 1
      
      # Score 2
      if (m2$ok) {
        ll2 <- glances(m2$model)$logLik
        if (ll2 > best_ll2) {
          best_ll2 <- ll2
          best_model2 <- c(model_TU, model_TD, model_EU)
        }
      } else n_fail <- n_fail + 1
      
      # Score 3
      if (m3$ok) {
        ll3 <- glances(m3$model)$logLik
        if (ll3 > best_ll3) {
          best_ll3 <- ll3
          best_model3 <- c(model_TU, model_TD, model_EU)
        }
      } else n_fail <- n_fail + 1
      
      # Score 4
      if (m4$ok) {
        ll4 <- glances(m4$model)$logLik
        if (ll4 > best_ll4) {
          best_ll4 <- ll4
          best_model4 <- c(model_TU, model_TD, model_EU)
        }
      } else n_fail <- n_fail + 1
      
      # Progreso
      if (n_fit %% 2 == 0) {
        cat("Iteraciones:", n_fit,
            "| Fallos:", n_fail, "\n")
      }
    }
  }
}

############################################################
# 6. Resumen de mejores estructuras encontradas
############################################################

list(
  best_e1 = list(model = best_model1, logLik = best_ll1),
  best_e2 = list(model = best_model2, logLik = best_ll2),
  best_e3 = list(model = best_model3, logLik = best_ll3),
  best_e4 = list(model = best_model4, logLik = best_ll4),
  total_fits = n_fit,
  total_failures = n_fail
)

############################################################
# 7. Ajuste final de modelos con estructuras óptimas
############################################################

ssn_mod1 <- ssn_lm(
  e_1 ~ 1,
  ssn.object = red_hidri_bogota_snn,
  tailup_type = best_model1[1],
  taildown_type = best_model1[2],
  euclid_type = best_model1[3],
  additive = "afv"
)

ssn_mod2 <- ssn_lm(
  e_2 ~ 1,
  ssn.object = red_hidri_bogota_snn,
  tailup_type = best_model2[1],
  taildown_type = best_model2[2],
  euclid_type = best_model2[3],
  additive = "afv"
)

ssn_mod3 <- ssn_lm(
  e_3 ~ 1,
  ssn.object = red_hidri_bogota_snn,
  tailup_type = best_model3[1],
  taildown_type = best_model3[2],
  euclid_type = best_model3[3],
  additive = "afv"
)

ssn_mod4 <- ssn_lm(
  e_4 ~ 1,
  ssn.object = red_hidri_bogota_snn,
  tailup_type = best_model4[1],
  taildown_type = best_model4[2],
  euclid_type = best_model4[3],
  additive = "afv"
)

############################################################
# 8. Almacenar modelos finales y metadatos
############################################################

modelos_ssn <- list(
  e1 = list(modelo = ssn_mod1,
            logLik = as.numeric(logLik(ssn_mod1)),
            estructura = best_model1),
  e2 = list(modelo = ssn_mod2,
            logLik = as.numeric(logLik(ssn_mod2)),
            estructura = best_model2),
  e3 = list(modelo = ssn_mod3,
            logLik = as.numeric(logLik(ssn_mod3)),
            estructura = best_model3),
  e4 = list(modelo = ssn_mod4,
            logLik = as.numeric(logLik(ssn_mod4)),
            estructura = best_model4)
)

############################################################
# 9. Guardar resultados
############################################################

saveRDS(modelos_ssn, "modelos_ssn_rio_bogota.rds")


modelos_ajustados
