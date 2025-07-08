##Estimando R-quadrado ajustado no modelo de regressão linear e calculando a moda

#Bibliotecas
library(car)
library(MASS)
library(tictoc)

tic("Questão 2")
#Estabelecimento dos números aleatórios
RNGkind('Marsaglia')

#Definindo parâmetros iniciais e variáveis
REP <- 10000
beta0 <- 3
beta1 <- 5
beta2 <- 7
n <- c(10,50,100,1000)

#Definindo a função
Rquadrados <- function(n){
  erro <- rnorm(n)
  x <- runif(n)
  y <- beta0 + beta1*x + erro
  fit <- lm(y~x)
  return(summary(fit)$adj.r.squared)}

#Simulando para diferentes tamanhos de amostras
set.seed(1305)

rquadreps <- data.frame(
  "n10" = replicate(REP, Rquadrados(n[1])),
  "n50" = replicate(REP, Rquadrados(n[2])),
  "n100" = replicate(REP, Rquadrados(n[3])),
  "n1000" = replicate(REP, Rquadrados(n[4]))
)

#Testando a distribuição de R-quadrado
# Distribuição Beta ajustada aos R² simulados

# Função de verossimilhança negativa para a distribuição Beta
Log_veross_neg <- function(par, x) {
  a <- par[1]
  b <- par[2]
  if (a <= 0 || b <= 0) return(Inf)
  -sum(dbeta(x, shape1 = a, shape2 = b, log = TRUE))
}

# Ajustando a distribuição Beta aos R² simulados
ajustar_beta_para_r2 <- function(r2) {
  r2 <- r2[r2 > 0 & r2 < 1]  
  
  mu <- mean(r2)
  sigma2 <- var(r2)
  a0 <- mu * ((mu * (1 - mu)) / sigma2 - 1)
  b0 <- (1 - mu) * ((mu * (1 - mu)) / sigma2 - 1)
  
  ajuste <- optim(
    par = c(a0, b0),
    fn = Log_veross_neg,
    x = r2,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    hessian = TRUE
  )
  
  a_hat <- ajuste$par[1]
  b_hat <- ajuste$par[2]
  
  hess_inv <- solve(ajuste$hessian)
  se <- sqrt(diag(hess_inv))
  IC_a <- a_hat + c(-1.96, 1.96) * se[1]
  IC_b <- b_hat + c(-1.96, 1.96) * se[2]
  
  return(list(a = a_hat, b = b_hat, IC_a = IC_a, IC_b = IC_b))
}

# Guardando os ajustes em uma lista nomeada
resultados_ajuste <- list(
  "10" = ajustar_beta_para_r2(rquadreps$n10),
  "50" = ajustar_beta_para_r2(rquadreps$n50),
  "100" = ajustar_beta_para_r2(rquadreps$n100),
  "1000" = ajustar_beta_para_r2(rquadreps$n1000)
)

# Montando tabela com os resultados
tabela_resultados <- do.call(rbind, lapply(names(resultados_ajuste), function(n) {
  ajuste <- resultados_ajuste[[n]]
  data.frame(
    n = as.integer(n),
    a = ajuste$a,
    a_lower = ajuste$IC_a[1],
    a_upper = ajuste$IC_a[2],
    b = ajuste$b,
    b_lower = ajuste$IC_b[1],
    b_upper = ajuste$IC_b[2]
  )
}))

# Exibindo tabela
print(tabela_resultados, row.names = FALSE)

# # Função para gerar QQ-Plot para um dado tamanho amostral
# plot_qq_r2 <- function(nome_coluna, resultados_ajuste) {
#   # Extraindo os dados de R² ajustado da coluna correspondente
#   r2_simulado <- rquadreps[[nome_coluna]]
#   r2_simulado <- r2_simulado[r2_simulado > 0 & r2_simulado < 1]  # beta só definida entre 0 e 1
#   
#   # Pegando os parâmetros ajustados
#   parametros_beta <- resultados_ajuste[[sub("n", "", nome_coluna)]]  # nome sem 'n'
#   a_hat <- parametros_beta$a
#   b_hat <- parametros_beta$b
#   
#   # Gerando o QQ-plot
#   qqPlot(
#     r2_simulado,
#     distribution = "beta",
#     shape1 = a_hat,
#     shape2 = b_hat,
#     main = paste("QQ Plot para R² simulado (", nome_coluna, ")"),
#     xlab = "Quantis Teóricos da Beta",
#     ylab = "Quantis Amostrais do R²",
#     col = "black",
#     pch = 19
#   )
# }
# 
# # Nomes das colunas em rquadreps
# colunas_r2 <- c("n10", "n50", "n100", "n1000")
# 
# # Gráfico 2x2 com os QQ-Plots
# par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
# 
# for (col in colunas_r2) {
#   plot_qq_r2(col, resultados_ajuste)
# }

# Função para gerar QQ-Plot para um dado tamanho amostral
plot_qq_r2 <- function(nome_coluna, resultados_ajuste) {
  # Extraindo os dados de R² ajustado da coluna correspondente
  r2_simulado <- rquadreps[[nome_coluna]]
  r2_simulado <- r2_simulado[r2_simulado > 0 & r2_simulado < 1]  # beta só definida entre 0 e 1
  
  # Pegando os parâmetros ajustados
  parametros_beta <- resultados_ajuste[[sub("n", "", nome_coluna)]]  # nome sem 'n'
  a_hat <- parametros_beta$a
  b_hat <- parametros_beta$b
  
  # Criando título com tamanho amostral
  n <- sub("n", "", nome_coluna)
  plot_title <- paste("n =", n, "| Erros normais")
  
  # Gerando o QQ-plot
  qqPlot(
    r2_simulado,
    distribution = "beta",
    shape1 = a_hat,
    shape2 = b_hat,
    main = plot_title,  # Usando main para o título
    xlab = "Quantis Teóricos da Beta",
    ylab = "Quantis Amostrais do R²",
    col = "black",
    pch = 19
  )
}

# Configuração da área de plotagem
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0), mar = c(4, 4, 2, 1))

# Gerar apenas n10 e n1000
plot_qq_r2("n10", resultados_ajuste)
plot_qq_r2("n1000", resultados_ajuste)

# Título principal reduzido e ajustado
mtext("QQ-Plots do R² Ajustado (n = 10 vs n = 1000)", 
      outer = TRUE, cex = 0.9, font = 2, line = 0.5)
mtext("Distribuição Beta Teórica com Erros Normais", 
      outer = TRUE, cex = 0.7, font = 1, line = -0.8)

# Função da moda teórica da Beta
moda_teorica <- function(a, b) {
  if (a > 1 && b > 1) {
    return((a - 1) / (a + b - 2))
  } else {
    return(NA)
  }
}

# Função para calcular a moda teórica e intervalo simétrico baseado nos quantis
calcular_moda_e_intervalo_quantil <- function(a, b, conf = 0.95) {
  # Calcula os quantis simétricos da Beta
  alfa <- (1 - conf) / 2
  IC <- qbeta(c(alfa, 1 - alfa), shape1 = a, shape2 = b)
  
  # Moda teórica
  moda <- moda_teorica(a, b)
  
  return(list(
    moda = moda,
    IC = IC
  ))
}

# Aplica a função para cada linha da tabela_resultados
modas_com_intervalo <- mapply(
  FUN = calcular_moda_e_intervalo_quantil,
  a = tabela_resultados$a,
  b = tabela_resultados$b,
  SIMPLIFY = FALSE
)

# Monta o data.frame final com os resultados
modas_beta_quantil <- data.frame(
  n = tabela_resultados$n,
  moda_teorica = sapply(modas_com_intervalo, function(x) x$moda),
  IC_lower = sapply(modas_com_intervalo, function(x) x$IC[1]),
  IC_upper = sapply(modas_com_intervalo, function(x) x$IC[2])
)

# Visualiza a tabela
print(modas_beta_quantil, row.names = FALSE)

#####################################################

#Para erros com distribuição T-studente g.l = 5

rquadrepsT <- read.csv('rquadrepsT')

resultados_ajusteT <- list(
  "10" = ajustar_beta_para_r2(rquadrepsT$n10),
  "50" = ajustar_beta_para_r2(rquadrepsT$n50),
  "100" = ajustar_beta_para_r2(rquadrepsT$n100),
  "1000" = ajustar_beta_para_r2(rquadrepsT$n1000)
)

# Montando tabela com os resultados
tabela_resultadosT <- do.call(rbind, lapply(names(resultados_ajusteT), function(n) {
  ajuste <- resultados_ajusteT[[n]]
  data.frame(
    n = as.integer(n),
    a = ajuste$a,
    a_lower = ajuste$IC_a[1],
    a_upper = ajuste$IC_a[2],
    b = ajuste$b,
    b_lower = ajuste$IC_b[1],
    b_upper = ajuste$IC_b[2]
  )
}))

# Exibindo tabela
print(tabela_resultadosT, row.names = FALSE)

#Modas T
# Aplica a função para cada linha da tabela_resultados
modas_com_intervaloT <- mapply(
  FUN = calcular_moda_e_intervalo_quantil,
  a = tabela_resultadosT$a,
  b = tabela_resultadosT$b,
  SIMPLIFY = FALSE
)

# Monta o data.frame final com os resultados
modas_beta_quantilT <- data.frame(
  n = tabela_resultados$n,
  moda_teorica = sapply(modas_com_intervaloT, function(x) x$moda),
  IC_lower = sapply(modas_com_intervaloT, function(x) x$IC[1]),
  IC_upper = sapply(modas_com_intervaloT, function(x) x$IC[2])
)

# Visualiza a tabela
print(modas_beta_quantil, row.names = FALSE)

#####################################

#Fazendo um modelo com  3 betas

#Definindo a função
Rquadrados2 <- function(n){
  erro <- rnorm(n)
  x1 <- runif(n)
  x2 <- runif(n)
  y <- beta0 + beta1*x1 + beta2*x2 + erro
  fit <- lm(y~ x1 + x2)
  return(summary(fit)$adj.r.squared)}

set.seed(1305)

rquadreps2 <- data.frame(
  "n10" = replicate(REP, Rquadrados2(n[1])),
  "n50" = replicate(REP, Rquadrados2(n[2])),
  "n100" = replicate(REP, Rquadrados2(n[3])),
  "n1000" = replicate(REP, Rquadrados2(n[4]))
)

# Guardando os ajustes em uma lista nomeada
resultados_ajuste2 <- list(
  "10" = ajustar_beta_para_r2(rquadreps2$n10),
  "50" = ajustar_beta_para_r2(rquadreps2$n50),
  "100" = ajustar_beta_para_r2(rquadreps2$n100),
  "1000" = ajustar_beta_para_r2(rquadreps2$n1000)
)

# Montando tabela com os resultados
tabela_resultados2 <- do.call(rbind, lapply(names(resultados_ajuste2), function(n) {
  ajuste <- resultados_ajuste2[[n]]
  data.frame(
    n = as.integer(n),
    a = ajuste$a,
    a_lower = ajuste$IC_a[1],
    a_upper = ajuste$IC_a[2],
    b = ajuste$b,
    b_lower = ajuste$IC_b[1],
    b_upper = ajuste$IC_b[2]
  )
}))

# Nomes das colunas em rquadreps
colunas_r2 <- c("n10", "n50", "n100", "n1000")

# Gráfico 2x2 com os QQ-Plots
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))

for (col in colunas_r2) {
  plot_qq_r2(col, resultados_ajuste2)
}

#Modas com 3 valores de beta
# Aplica a função para cada linha da tabela_resultados
modas_com_intervalo2 <- mapply(
  FUN = calcular_moda_e_intervalo_quantil,
  a = tabela_resultados2$a,
  b = tabela_resultados2$b,
  SIMPLIFY = FALSE
)

# Monta o data.frame final com os resultados
modas_beta_quantil2 <- data.frame(
  n = tabela_resultados$n,
  moda_teorica = sapply(modas_com_intervalo2, function(x) x$moda),
  IC_lower = sapply(modas_com_intervalo2, function(x) x$IC[1]),
  IC_upper = sapply(modas_com_intervalo2, function(x) x$IC[2])
)
toc()

write.csv(modas_beta_quantil, "Moda_Base.csv")
write.csv(modas_beta_quantil2, "Moda_Base2.csv")
