#Brincando
library(nortest)
library(ggplot2)
library(dplyr)
library(car)
library(moments)
library(lmtest)


ns <- c(10, 100, 500, 1000)
beta0 <- 3
beta1 <- 5
beta2 <- 4
beta3 <- 6
B <- 1000

# --------------------------
# Modelo 1: Y ~ X
simular_modelo1 <- function(n) {
  resultados <- replicate(B, {
    X <- runif(n)
    e <- rnorm(n)
    Y <- beta0 + beta1 * X + e
    fit <- lm(Y ~ X)
    c(AIC = AIC(fit), BIC = BIC(fit))
  })
  list(AIC = resultados["AIC", ], BIC = resultados["BIC", ])
}

# --------------------------
# Modelo 2: Y ~ X + Z
simular_modelo2 <- function(n) {
  resultados <- replicate(B, {
    X <- runif(n)
    Z <- runif(n)
    e <- rnorm(n)
    Y <- beta0 + beta1 * X + beta2 * Z + e
    fit <- lm(Y ~ X + Z)
    c(AIC = AIC(fit), BIC = BIC(fit))
  })
  list(AIC = resultados["AIC", ], BIC = resultados["BIC", ])
}

# --------------------------
# Modelo 3: Y ~ X + Z + W
simular_modelo3 <- function(n) {
  resultados <- replicate(B, {
    X <- runif(n)
    Z <- runif(n)
    W <- runif(n)
    e <- rnorm(n)
    Y <- beta0 + beta1 * X + beta2 * Z + beta3 * W + e
    fit <- lm(Y ~ X + Z + W)
    c(AIC = AIC(fit), BIC = BIC(fit))
  })
  list(AIC = resultados["AIC", ], BIC = resultados["BIC", ])
}

# --------------------------
# Rodar simulações para cada tamanho de amostra
resultados_modelo1 <- lapply(ns, simular_modelo1)
resultados_modelo2 <- lapply(ns, simular_modelo2)
resultados_modelo3 <- lapply(ns, simular_modelo3)

names(resultados_modelo1) <- as.character(ns)
names(resultados_modelo2) <- as.character(ns)
names(resultados_modelo3) <- as.character(ns)

# Função auxiliar para extrair as métricas
processar_resultados <- function(resultados) {
  df_final <- data.frame()
  
  for (n_str in names(resultados)) {
    n_val <- as.integer(n_str)
    res <- resultados[[n_str]]
    
    # Criar data.frame para AIC
    df_aic <- data.frame(
      n = n_val,
      criterio = "AIC",
      media = mean(res$AIC),
      variancia = var(res$AIC),
      curtose = kurtosis(res$AIC) - 3,
      assimetria = skewness(res$AIC),
      stringsAsFactors = FALSE
    )
    
    # Criar data.frame para BIC
    df_bic <- data.frame(
      n = n_val,
      criterio = "BIC",
      media = mean(res$BIC),
      variancia = var(res$BIC),
      curtose = kurtosis(res$BIC) - 3,
      assimetria = skewness(res$BIC),
      stringsAsFactors = FALSE
    )
    
    # Combinar os resultados
    df_final <- rbind(df_final, df_aic, df_bic)
  }
  
  return(df_final)
}

# Gerar o dataframe de resultados
resultados_finais <- processar_resultados(resultados_modelo1)


#--------------------------------------------------------------
#Histogramas

# Combinar os resultados em um único data frame
df_plot <- do.call(rbind, lapply(seq_along(ns), function(i) {
  data.frame(
    n = ns[i],
    AIC = resultados_modelo1[[i]]$AIC,
    BIC = resultados_modelo1[[i]]$BIC
  )
}))

# 1. Histograma para o AIC
ggplot(df_plot, aes(x = AIC)) +
  geom_histogram(fill = "#81A1C1", color = "white", bins = 30, alpha = 0.8) +
  facet_wrap(~n, scales = "free", labeller = label_both) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribuição do Critério AIC por tamanho da amostra",
    x = "Valor do AIC",
    y = "Frequência"
  ) +
  theme(
    strip.text = element_text(face = "bold")
  )

# 2. Histograma para o BIC
ggplot(df_plot, aes(x = BIC)) +
  geom_histogram(fill = "#BF616A", color = "white", bins = 30, alpha = 0.8) +
  facet_wrap(~n, scales = "free", labeller = label_both) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribuição do Critério BIC por tamanho da amostra",
    x = "Valor do BIC",
    y = "Frequência"
  ) +
  theme(
    strip.text = element_text(face = "bold")
  )

#-----------------------------------------------------------------------------
# Testando independência só para AIC

durbin_resultado <- function(valores, criterio, n) {
  lag_1 <- c(NA, head(valores, -1))
  df <- na.omit(data.frame(v = valores, v_lag = lag_1))
  modelo <- lm(v ~ v_lag, data = df)
  resultado_dw <- dwtest(modelo)
  
  data.frame(
    Tamanho_amostral = as.integer(n),
    DW = resultado_dw$statistic,
    p_value = resultado_dw$p.value
  )
}

resultados_dw <- do.call(rbind, lapply(names(resultados_modelo1), function(n) {
  aic_vals <- resultados_modelo1[[n]]$AIC
  durbin_resultado(aic_vals, "AIC", n)
}))

print(resultados_dw)

# ----------------------------------------------------------------------------
#Testes de normalidade

# Função atualizada com Anderson-Darling
aplicar_testes_normalidade <- function(vetor) {
  ad <- ad.test(vetor)
  lf <- lillie.test(vetor)
  return(data.frame(
    Anderson_Darling_p = ad$p.value,
    Lilliefors_p       = lf$p.value
  ))
}

# Testes para AIC
cat("---- Testes de normalidade para AIC ----\n")
for (i in seq_along(ns)) {
  n <- ns[i]
  aic_vals <- resultados_modelo1[[i]]$AIC
  cat(paste0("n = ", n, "\n"))
  print(aplicar_testes_normalidade(aic_vals))
  cat("\n")
}

# Testes para BIC
cat("---- Testes de normalidade para BIC ----\n")
for (i in seq_along(ns)) {
  n <- ns[i]
  bic_vals <- resultados_modelo1[[i]]$BIC
  cat(paste0("n = ", n, "\n"))
  print(aplicar_testes_normalidade(bic_vals))
  cat("\n")
}


#----------------------------------------------------------------
#QQ-Plot
  
# QQ-Plot para AIC
cat("==== QQ-Plots para AIC ====\n")
par(mfrow = c(2, 2))  # Layout 2x2
for (i in seq_along(ns)) {
  n <- ns[i]
  aic_vals <- resultados_modelo1[[i]]$AIC
  
  qqPlot(aic_vals,
         main = paste("QQ-Plot do AIC - n =", n),
         ylab = "Quantis do AIC",
         xlab = "Quantis teóricos normais",
         col = "black",      # pontos pretos
         pch = 20,
         envelope = 0.95)
}

# QQ-Plot para BIC
cat("==== QQ-Plots para BIC ====\n")
par(mfrow = c(2, 2))  # Reset layout 2x2
for (i in seq_along(ns)) {
  n <- ns[i]
  bic_vals <- resultados_modelo1[[i]]$BIC
  
  qqPlot(bic_vals,
         main = paste("QQ-Plot do BIC - n =", n),
         ylab = "Quantis do BIC",
         xlab = "Quantis teóricos normais",
         col = "black",     # pontos pretos
         pch = 20,
         envelope = 0.95)
}

# Reset layout para 1 gráfico por vez (opcional)
par(mfrow = c(1, 1))


# ------------------------------------------------------------

#-----------------------
