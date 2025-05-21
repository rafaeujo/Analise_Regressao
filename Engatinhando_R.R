install.packages("car")
library(car)
RNGkind('Marsaglia')
set.seed(1234)

#Simulando Regressão Linear

#Definindo parâmetros e variáveis

REP <- 10000
n <- 50

beta0 <- 3
beta1 <- 5

#Definindo a função

Simulando <- function(){
  erro <- rnorm(n)
  x <- runif(n)
  y <- beta0 + beta1*x + erro
  fit <- lm(y~x)
  return(fit$coefficients)
}

repli <- replicate(REP, Simulando())

#Análises

mean(repli[1,])
mean(repli[2,])
hist(repli[1,], main = "Histograma de beta0")
hist(repli[2,], main = "Histograma de beta1")

#Testando a normalidade dos erros

par(mfrow = c(1,2))
qqPlot(repli[1,], envelope = 0.95)
qqPlot(repli[2,], envelope = 0.95)

#Mudando a distribuição dos erros

beta0_T <- 3
beta1_T <- 5

Simulando_T <- function(){
  erro2 <- rt(n, df = 3)
  x <- runif(n)
  y <- beta0_T + beta1_T*x + erro2
  fit <- lm(y~x)
  return(fit$coefficients)
}

repli_T <- replicate(REP, Simulando_T())

par(mfrow = c(2,2))
qqPlot(repli_T[1,], envelope = 0.95)
qqPlot(repli_T[2,], envelope = 0.95)
qqPlot(repli[1,], envelope = 0.95)
qqPlot(repli[2,], envelope = 0.95)

var(repli_T[1,])
var(repli_T[2,])

#Teorema Central do Limite

Cauchy <- function(){
  n2 <-  1000
  x_cauchy <-  rcauchy(n2)
  return(mean(x_cauchy))
}
par(mfrow  = c(1,1))
repli2 <- replicate(REP, Cauchy())
qqPlot(repli2, distribution = 't', df = 1)
hist(repli2)

# Função modificada para testar consistência
Simulando_consistencia <- function(n) {
  erro <- rnorm(n)
  x <- runif(n)
  y <- beta0 + beta1*x + erro
  fit <- lm(y~x)
  return(fit$coefficients)
}

# Testando para diferentes tamanhos de amostra
tamanhos <- c(50, 100, 500, 1000, 5000)
resultados <- list()

for (i in seq_along(tamanhos)) {
  repli <- replicate(REP, Simulando_consistencia(tamanhos[i]))
  resultados[[i]] <- list(
    n = tamanhos[i],
    media_beta0 = mean(repli[1,]),
    media_beta1 = mean(repli[2,]),
    var_beta0 = var(repli[1,]),
    var_beta1 = var(repli[2,])
  )
}

# Exibindo os resultados
for (res in resultados) {
  cat("\nTamanho da amostra:", res$n,
      "\nMédia de beta0:", res$media_beta0,
      "\nMédia de beta1:", res$media_beta1,
      "\nVariância de beta0:", res$var_beta0,
      "\nVariância de beta1:", res$var_beta1, "\n")
}

# Função de simulação com erros t-Student
Simulando_t <- function(n) {
  erro <- rt(n, df = 3)  # Erros com distribuição t(3)
  x <- runif(n)          # Variável explicativa uniforme
  y <- beta0 + beta1*x + erro
  fit <- lm(y ~ x)
  return(fit$coefficients)
}

# Tamanhos de amostra para testar consistência
tamanhos <- c(50, 100, 500, 1000, 5000)
resultados <- list()

# Executando a simulação para cada tamanho de amostra
for (i in seq_along(tamanhos)) {
  repli <- replicate(REP, Simulando_t(tamanhos[i]))
  resultados[[i]] <- list(
    n = tamanhos[i],
    media_beta0 = mean(repli[1,]),
    media_beta1 = mean(repli[2,]),
    var_beta0 = var(repli[1,]),
    var_beta1 = var(repli[2,])
  )
}

# Exibindo os resultados
for (res in resultados) {
  cat("\nTamanho da amostra:", res$n,
      "\nMédia de beta0:", round(res$media_beta0, 4),
      "\nMédia de beta1:", round(res$media_beta1, 4),
      "\nVariância de beta0:", round(res$var_beta0, 6),
      "\nVariância de beta1:", round(res$var_beta1, 6), "\n")
}

correlacao_empirica <- cor(t(repli)[,1], t(repli)[,2])
