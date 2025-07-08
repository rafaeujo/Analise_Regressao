RNGkind("Marsaglia")

library(stats)
library(strucchange)
library(car)
library(tictoc)
library(nortest)
library(e1071)
library(ggplot2)
library(dplyr)
library(tidyr)

tic("Questão 1 T")

beta0 <- 3
beta1 <- 5
ns <- c(10, 50, 100, 1000)
B <- 10000
gl <-  5

simulacao_normal <- function(n) {
  set.seed(1305)
  X <- rnorm(n)
  
  replicate(B, {
    Y <- beta0 + beta1 * X + rt(n, gl)
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
# Configuração de cores para distribuição T-Student (escala de vermelhos)
cores_tstudent <- c("#990000", "#CC0000", "#FF3333", "#FF6666")

# Gráfico de Médias e Variâncias dos Resíduos (T-Student)
ggplot(dados_plot, aes(x = factor(n), y = valor, color = tipo, group = tipo)) +
  geom_line(size = 1) +
  geom_point(aes(shape = tipo), size = 3) +
  scale_color_manual(values = cores_tstudent) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(
    title = "Médias e Variâncias(Erros ~ T-Student)",
    x = "Tamanho da amostra (n)",
    y = "Valor",
    color = "Resíduo",
    shape = "Resíduo"
  ) +
  facet_wrap(~estatistica, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
    
  )

###############################################################################
#Normalidade
###############################################################################
calcular_proporcoes_normalidade_multiteste <- function(simulacoes, ns, alpha = 0.05, B = length(simulacoes[[1]])) {
  # Carregar pacotes necessários
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
      res <- sims[[j]]
      
      # Teste de Shapiro-Wilk (apenas para n <= 5000)
      if (n <= 5000) {
        rej_shapiro[1] <- rej_shapiro[1] + (shapiro.test(res$ordinary)$p.value < alpha)
        rej_shapiro[2] <- rej_shapiro[2] + (shapiro.test(res$standard)$p.value < alpha)
        rej_shapiro[3] <- rej_shapiro[3] + (shapiro.test(res$student)$p.value < alpha)
        rej_shapiro[4] <- rej_shapiro[4] + (shapiro.test(res$recursive)$p.value < alpha)
      }
      
      # Teste de Anderson-Darling
      rej_ad[1] <- rej_ad[1] + (ad.test(res$ordinary)$p.value < alpha)
      rej_ad[2] <- rej_ad[2] + (ad.test(res$standard)$p.value < alpha)
      rej_ad[3] <- rej_ad[3] + (ad.test(res$student)$p.value < alpha)
      rej_ad[4] <- rej_ad[4] + (ad.test(res$recursive)$p.value < alpha)
      
      # Teste Jarque-Bera
      rej_jb[1] <- rej_jb[1] + (jarque.bera.test(res$ordinary)$p.value < alpha)
      rej_jb[2] <- rej_jb[2] + (jarque.bera.test(res$standard)$p.value < alpha)
      rej_jb[3] <- rej_jb[3] + (jarque.bera.test(res$student)$p.value < alpha)
      rej_jb[4] <- rej_jb[4] + (jarque.bera.test(res$recursive)$p.value < alpha)
    }
    
    # Calcular proporções de rejeição
    proporcoes <- data.frame(
      Residuo = c("ordinary", "standard", "student", "recursive"),
      Shapiro_Wilk = rej_shapiro/B,
      Anderson_Darling = rej_ad/B,
      Jarque_Bera = rej_jb/B,
      n = n
    )
    
    resultados[[i]] <- proporcoes
  }
  
  do.call(rbind, resultados)
}
# Primeiro, vamos obter os resultados de normalidade para a simulação normal
proporcoes_normalidade <- calcular_proporcoes_normalidade_multiteste(simulacoes_normal, ns)

# Transformar os dados para formato longo (melhor para ggplot)
library(tidyr)
library(ggplot2)

# Primeiro, certifique-se que os dados estão no formato longo
dados_grafico <- proporcoes_normalidade %>%
  pivot_longer(
    cols = c("Shapiro_Wilk", "Anderson_Darling", "Jarque_Bera"),  # Nomes entre aspas
    names_to = "Teste",
    values_to = "Proporcao_Rejeicao"
  )

# Criar o gráfico com ajustes
# Gráfico de Testes de Normalidade (T-Student)
ggplot(dados_grafico, aes(x = factor(n), y = Proporcao_Rejeicao, color = Residuo, group = Residuo)) +
  geom_line(size = 1) +
  geom_point(aes(shape = Residuo), size = 3) +
  scale_color_manual(values = cores_tstudent) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(
    title = "Testes de Normalidade (Erros ~ T-Student)",
    x = "Tamanho da Amostra (n)",
    y = "Proporção de Rejeição",
    color = "Resíduo",
    shape = "Resíduo"
  ) +
  facet_wrap(~ Teste, ncol = 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
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
ggplot(autocorrelacoes_medias_lags, aes(x = lag, y = autocorrelacao_media, color = tipo)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = cores_tstudent) +
  labs(
    title = "Autocorrelação(Erros ~ T-Student)",
    subtitle = "Lags 1 a 5",
    x = "Lag",
    y = "Autocorrelação média",
    color = "Tipo de resíduo"
  ) +
  facet_wrap(~ n, ncol = 2) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

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

# Transformar para formato longo
dados_longos <- resultados_momentos %>%
  pivot_longer(cols = c(media_assimetria, media_curtose),
               names_to = "momento", values_to = "valor")

# Gráfico de Assimetria e Curtose (T-Student)
ggplot(dados_longos, aes(x = factor(n), y = valor, color = tipo, group = tipo)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = cores_tstudent) +
  labs(
    title = "Assimetria e Curtose (Erros ~ T-Student)",
    x = "Tamanho da Amostra (n)",
    y = "Valor",
    color = "Tipo de Resíduo"
  ) +
  facet_wrap(~ momento, scales = "free_y") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

toc()
