##Estimando R-quadrado ajustado no modelo de regressão linear

#Bibliotecas
library(car)
library(MASS)

#Estabelecimento dos números aleatórios
RNGkind('Marsaglia')

#Definindo parâmetros iniciais e variáveis
REP <- 10000
beta0 <- 3
beta1 <- 5
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

system.time(rquadreps <- data.frame(
  "n10" = replicate(REP, Rquadrados(n[1])),
  "n50" = replicate(REP, Rquadrados(n[2])),
  "n100" = replicate(REP, Rquadrados(n[3])),
  "n1000" = replicate(REP, Rquadrados(n[4]))
)
)
write.csv(rquadreps, file = "resultados_r2_ajustado.csv", row.names = TRUE)
rquadreps <- read.csv("resultados_r2_ajustado.csv")


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
system.time(resultados_ajuste <- list(
  "10" = ajustar_beta_para_r2(rquadreps$n10),
  "50" = ajustar_beta_para_r2(rquadreps$n50),
  "100" = ajustar_beta_para_r2(rquadreps$n100),
  "1000" = ajustar_beta_para_r2(rquadreps$n1000)
))

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

# Função para gerar QQ-Plot para um dado tamanho amostral
plot_qq_r2 <- function(nome_coluna, resultados_ajuste) {
  # Extraindo os dados de R² ajustado da coluna correspondente
  r2_simulado <- rquadreps[[nome_coluna]]
  r2_simulado <- r2_simulado[r2_simulado > 0 & r2_simulado < 1]  # beta só definida entre 0 e 1
  
  # Pegando os parâmetros ajustados
  parametros_beta <- resultados_ajuste[[sub("n", "", nome_coluna)]]  # nome sem 'n'
  a_hat <- parametros_beta$a
  b_hat <- parametros_beta$b
  
  # Gerando o QQ-plot
  qqPlot(
    r2_simulado,
    distribution = "beta",
    shape1 = a_hat,
    shape2 = b_hat,
    main = paste("QQ Plot para R² simulado (", nome_coluna, ")"),
    xlab = "Quantis Teóricos da Beta",
    ylab = "Quantis Amostrais do R²",
    col = "black",
    pch = 19
  )
}

# Nomes das colunas em rquadreps
colunas_r2 <- c("n10", "n50", "n100", "n1000")

# Gráfico 2x2 com os QQ-Plots
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))

for (col in colunas_r2) {
  plot_qq_r2(col, resultados_ajuste)
}
