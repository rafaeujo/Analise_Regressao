RNGkind("Marsaglia")

library(stats)
library(strucchange)
library(car)
library(tictoc)
library(nortest)
library(e1071)

tic("Questão 1")

beta0 <- 3
beta1 <- 5
ns <- c(10, 50, 100, 1000)
B <- 10000

simulacao_normal <- function(n) {
  set.seed(1305)
  X <- rnorm(n)
  
  replicate(B, {
    Y <- beta0 + beta1 * X + rnorm(n)
    fit <- lm(Y ~ X)
    list(
      fit = fit,
      ordinary = residuals(fit),
      standard = rstandard(fit),
      student = rstudent(fit),
      recursive = recresid(fit),
      fitted = fitted(fit),
      coefficients = coef(fit)
    )
  }, simplify = FALSE)
}

simulacoes_normal <- lapply(ns, simulacao_normal)

resumo_residuos <- function(simulacoes) {
  tipos <- c("ordinary", "standard", "student", "recursive")
  
  sapply(tipos, function(tipo) {
    medias <- sapply(simulacoes, function(sim) mean(sim[[tipo]]))
    variancias <- sapply(simulacoes, function(sim) var(sim[[tipo]]))
    
    c(
      media_das_medias = mean(medias),
      media_das_variancias = mean(variancias),
    )
  })
}

# Extrair e organizar os dados
resultados <- lapply(simulacoes_normal, function(sim) {
  sapply(sim, function(run) {
    tipos <- c("ordinary", "standard", "student", "recursive")
    sapply(tipos, function(tipo) {
      c(mean = mean(run[[tipo]]), var = var(run[[tipo]]))
    })
  }, simplify = "array")
})

# Preparar dados para plotagem
library(ggplot2)
library(dplyr)

dados_plot <- data.frame()

for (i in seq_along(ns)) {
  n <- ns[i]
  sim <- simulacoes_normal[[i]]
  
  for (tipo in c("ordinary", "standard", "student", "recursive")) {
    # Calcular médias e variâncias para todas as replicações
    medias <- sapply(sim, function(run) mean(run[[tipo]]))
    variancias <- sapply(sim, function(run) var(run[[tipo]]))
    
    dados_plot <- rbind(dados_plot,
                        data.frame(
                          n = n,
                          tipo = tipo,
                          estatistica = "Média",
                          valor = mean(medias)
                        ),
                        data.frame(
                          n = n,
                          tipo = tipo,
                          estatistica = "Variância",
                          valor = mean(variancias)
                        ))
  }
}


ggplot(dados_plot, aes(x = n, y = valor, color = tipo, shape = tipo)) +
  geom_line(aes(linetype = tipo), size = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = ns) +
  facet_wrap(~estatistica, scales = "free_y") +
  labs(title = NULL,
       x = "Tamanho da amostra (n)",
       y = NULL,
       color = "Tipo de resíduo",
       shape = "Tipo de resíduo",
       linetype = "Tipo de resíduo") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

###############################################################################
#Normalidade
###############################################################################
calcular_proporcoes_normalidade_multiteste <- function(simulacoes, ns, alpha = 0.05) {
  # Carregar pacotes necessários para os testes
  if (!require(nortest)) install.packages("nortest"); library(nortest) # Para Anderson-Darling
  if (!require(tseries)) install.packages("tseries"); library(tseries) # Para Jarque-Bera
  
  resultados <- list()
  
  for (i in seq_along(ns)) {
    n <- ns[i]
    sims <- simulacoes[[i]]
    
    # Inicializar contadores para cada teste e tipo de resíduo
    rej_shapiro <- rej_ad <- rej_jb <- 
      matrix(0, nrow = 4, ncol = 1, 
             dimnames = list(c("ordinary", "standard", "student", "recursive"), NULL))
    
    for (j in 1:B) {
      # Extrair resíduos
      res <- sims[[j]]
      
      # Teste de Shapiro-Wilk
      if (n <= 5000) {  # Shapiro-Wilk tem limite de 5000 observações
        rej_shapiro[1] <- rej_shapiro[1] + (shapiro.test(res$ordinary)$p.value < alpha)
        rej_shapiro[2] <- rej_shapiro[2] + (shapiro.test(res$standard)$p.value < alpha)
        rej_shapiro[3] <- rej_shapiro[3] + (shapiro.test(res$student)$p.value < alpha)
        rej_shapiro[4] <- rej_shapiro[4] + (shapiro.test(res$recursive)$p.value < alpha)
      }
      
      # Teste de Lilliefors (Kolmogorov-Smirnov com parâmetros estimados)
      rej_lillie[1] <- rej_lillie[1] + (lillie.test(res$ordinary)$p.value < alpha)
      rej_lillie[2] <- rej_lillie[2] + (lillie.test(res$standard)$p.value < alpha)
      rej_lillie[3] <- rej_lillie[3] + (lillie.test(res$student)$p.value < alpha)
      rej_lillie[4] <- rej_lillie[4] + (lillie.test(res$recursive)$p.value < alpha)
      
      
      # Teste Jarque-Bera
      rej_jb[1] <- rej_jb[1] + (jarque.bera.test(res$ordinary)$p.value < alpha)
      rej_jb[2] <- rej_jb[2] + (jarque.bera.test(res$standard)$p.value < alpha)
      rej_jb[3] <- rej_jb[3] + (jarque.bera.test(res$student)$p.value < alpha)
      rej_jb[4] <- rej_jb[4] + (jarque.bera.test(res$recursive)$p.value < alpha)
    }
    
    # Calcular proporções
    proporcoes <- data.frame(
      Residuo = c("ordinary", "standard", "student", "recursive"),
      Shapiro_Wilk = rej_shapiro/B,
      Anderson_Darling = rej_ad/B,
      Jarque_Bera = rej_jb/B,
      n = n
    )
    
    resultados[[i]] <- proporcoes
  }
  
  # Combinar todos os resultados em um único data.frame
  do.call(rbind, resultados)
}

# Primeiro, vamos obter os resultados de normalidade para a simulação normal
proporcoes_normalidade <- calcular_proporcoes_normalidade_multiteste(simulacoes_normal, ns)

# Carregar pacote ggplot2 para gráficos
library(ggplot2)
library(tidyr)  # Para pivot_longer

# Transformar os dados para formato longo (melhor para ggplot)
dados_grafico <- proporcoes_normalidade %>%
  pivot_longer(
    cols = c(Shapiro_Wilk, Anderson_Darling, Jarque_Bera),
    names_to = "Teste",
    values_to = "Proporcao_Rejeicao"
  )

ggplot(dados_grafico, aes(x = n, y = Proporcao_Rejeicao, 
                          color = Residuo, linetype = Residuo)) +
  geom_line(size = 1) +
  geom_point(size = 2, shape = 21, fill = "white") +
  facet_wrap(~ Teste, ncol = 1, scales = "fixed") +
  scale_x_continuous(breaks = unique(dados_grafico$n)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Proporção de Rejeição por Teste de Normalidade",
    subtitle = "Simulação com Erros Normais",
    x = "Tamanho da Amostra (n)",
    y = "Proporção de Rejeição",
    color = "Tipo de Resíduo",
    linetype = "Tipo de Resíduo"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )



###############################################################################
# Autocorrelacao
###############################################################################

calcular_autocorrelacao_media_multilag <- function(simulacoes, ns, max_lag = 5) {
  resultados <- data.frame()
  tipos_residuos <- c("ordinary", "standard", "student", "recursive")
  
  for (i in seq_along(ns)) {
    n <- ns[i]
    sims <- simulacoes[[i]]
    
    for (tipo in tipos_residuos) {
      # Calcular autocorrelações para cada simulação
      autocorrs_lags <- sapply(sims, function(sim) {
        acf(sim[[tipo]], plot = FALSE, lag.max = max_lag)$acf[2:(max_lag + 1)]  # Lags 1 a 5
      })
      
      # Se autocorrs_lags for vetor (B = 1), transformar em matriz
      if (is.vector(autocorrs_lags)) {
        autocorrs_lags <- matrix(autocorrs_lags, nrow = max_lag, ncol = 1)
      }
      
      # Calcular média por lag
      medias_por_lag <- rowMeans(autocorrs_lags)
      
      # Montar data frame
      df <- data.frame(
        n = n,
        tipo = tipo,
        lag = 1:max_lag,
        autocorrelacao_media = medias_por_lag
      )
      
      resultados <- rbind(resultados, df)
    }
  }
  
  return(resultados)
}

# Calcular autocorrelações médias para lags de 1 a 5
autocorrelacoes_medias_lags <- calcular_autocorrelacao_media_multilag(simulacoes_normal, ns)

# Plotar
ggplot(autocorrelacoes_medias_lags, aes(x = lag, y = autocorrelacao_media, 
                                        color = tipo, group = tipo)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  facet_wrap(~ n, ncol = 2, labeller = label_both) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Autocorrelação média dos resíduos (lags 1 a 5)",
       x = "Lag",
       y = "Autocorrelação média",
       color = "Tipo de resíduo") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

################################################################################
gerar_acfplots <- function(simulacao, n, max_lag = 20) {
  par(mfrow = c(2, 2))  # Layout 2x2
  tipos <- c("ordinary", "standard", "student", "recursive")
  
  for (tipo in tipos) {
    acf(simulacao[[tipo]],
        main = paste("ACF -", tipo, "\nn =", n),
        lag.max = max_lag)
  }
}

# Gerar ACF plots para as 4 simulações aleatórias de cada tamanho de amostra
for (i in seq_along(ns)) {
  n <- ns[i]
  sims <- simulacoes_aleatorias[[i]]
  
  for (j in 1:1) {
    cat("Simulação", j, "para n =", n, "\n")
    gerar_acfplots(sims[[j]], n, max_lag = 20)
    readline(prompt = "Pressione [Enter] para continuar...")  # Pausa entre os plots
  }
}


################################################################################
if (!require(e1071)) install.packages("e1071")


calcular_assimetria_curtose <- function(simulacoes, ns) {
  tipos_residuos <- c("ordinary", "standard", "student", "recursive")
  resultados <- data.frame()
  
  for (i in seq_along(ns)) {
    n <- ns[i]
    sims <- simulacoes[[i]]
    
    for (tipo in tipos_residuos) {
      skewness_vals <- sapply(sims, function(sim) skewness(sim[[tipo]]))
      kurtosis_vals <- sapply(sims, function(sim) kurtosis(sim[[tipo]]))
      
      resultados <- rbind(resultados,
                          data.frame(
                            n = n,
                            tipo = tipo,
                            media_assimetria = mean(skewness_vals),
                            dp_assimetria = sd(skewness_vals),
                            media_curtose = mean(kurtosis_vals),
                            dp_curtose = sd(kurtosis_vals)
                          ))
    }
  }
  
  return(resultados)
}

resultados_momentos <- calcular_assimetria_curtose(simulacoes_normal, ns)

library(tidyr)
library(ggplot2)

# Transformar para formato longo
dados_longos <- resultados_momentos %>%
  pivot_longer(cols = c(media_assimetria, media_curtose),
               names_to = "momento", values_to = "valor")

ggplot(dados_longos, aes(x = factor(n), y = valor, color = tipo, group = tipo)) +
  geom_line(linewidth = 1, linetype = "solid") +
  geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1.2) +
  facet_wrap(~ momento, scales = "free_y", labeller = label_both) +
  labs(
    title = "Média da Assimetria e Curtose dos Resíduos",
    x = "Tamanho da Amostra (n)",
    y = "Valor Médio",
    color = "Tipo de Resíduo"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "#e0e0e0"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed")
  )



ljung_box_residuos <- function(simulacoes, lags = 5) {
  tipos <- c("ordinary", "standard", "student", "recursive")
  resultados <- data.frame()
  
  for (tipo in tipos) {
    pvals <- sapply(simulacoes, function(sim) {
      bt <- Box.test(sim[[tipo]], lag = lags, type = "Ljung-Box")
      return(bt$p.value)
    })
    
    resultados <- rbind(resultados,
                        data.frame(
                          tipo = tipo,
                          prop_rejeicao = mean(pvals < 0.05)
                        ))
  }
  return(resultados)
}

ljung_todos <- lapply(simulacoes_normal, ljung_box_residuos, lags = 5)

# Adiciona informação de n
for (i in seq_along(ns)) {
  ljung_todos[[i]]$n <- ns[i]
}

ljung_df <- do.call(rbind, ljung_todos)

library(ggplot2)

ggplot(ljung_df, aes(x = factor(n), y = prop_rejeicao, color = tipo, group = tipo)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(title = "Proporção de rejeição do Teste de Ljung-Box (lag = 5)",
       x = "Tamanho da amostra (n)",
       y = "Proporção de rejeição",
       color = "Tipo de Resíduo") +
  theme_minimal()
toc()