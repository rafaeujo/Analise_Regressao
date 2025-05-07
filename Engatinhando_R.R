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

# Tarefa: Criar uma função que receba um vetor de letras e retorne cada letra e o número que estará associado à letra

~/Rafael Jordane/Analise_Regressao
