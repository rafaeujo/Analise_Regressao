install.packages("car")
library(car)
RNGkind('Marsaglia')
set.seed(1234)

#Simulando Regressão Linear

#Definindo parâmetros e variáveis

REP <- 10000
n <- 10

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

qqPlot(erro, envelope = 0.99)
