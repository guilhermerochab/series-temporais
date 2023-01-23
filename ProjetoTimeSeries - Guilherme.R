#Carregamento das bibliotecas
library("urca")
library("tseries")
library("forecast")
library(astsa)
library(MLmetrics)

#Carregamento dos dados
setwd('/Users/guilhermerocha/Mestrado/Séries Temporais e Previsão/projeto/dados')
dados = read.csv('aod_-30_40.txt', sep = '\ ')

#Verificação dos dados nulos
sapply(dados, function(dados) sum(is.na(dados)))

#Correção do dado da coluna 220
dados[220,]$X.2 = dados[220,]$X.3

#Filtragem das colunas necessárias (Data e AOC)
dados = dados[,c(2,4)]

#Renomear as colunas
colnames(dados) = c('Mês', 'AOD')

#Criação da série temporal
serie = ts(dados$AOD, start=2003, frequency = 12)

#Gráfico da Série Temporal
plot(serie,xlab="Tempo",ylab="AOD",main="Profundidade Óptica do Aerossol (AOD)", lw=2)

#Histograma
hist(serie, main="Distribuição dos Dados", xlab = "AOD", ylab='Frequência')

#Teste Q-Q
qqnorm(serie, main="Gráfico Q-Q Normal", xlab="Quantis Teóricos", ylab="Quantis da Amostra")
qqline(serie)

#Teste de normalidade de Shapiro-Wilk
shapiro.test(serie)

#Teste de estacionaridade KPSS (Kwiatkowski-Phillips-Schmidt-Shin)
kpss.test(serie)

#Teste de estacionaridade Philips-Perron
pp.test(serie)

#Teste de estacionaridade Dickey Fuller
adf.test(serie)

#Gráfico da Função Autocorrelação (ACF)
acf(serie, main="Função de Autocorrelação (ACF)")

#Gráfico da Função Autocorrelação Parcial (PACF)
pacf(serie, main="Função de Autocorrelação Parcial (PACF)", ylab="PACF")

#Teste de Autocorrelação (Ljung-Box)
Box.test(serie, type = "Ljung-Box")

#Gráfico de decompisção classica
plot(decompose(serie))

#Gráfico de decomposição multiplicativa
plot(decompose(serie, type="multiplicative"))

#Decomposição STL
stl.11=stl(serie,s.window=11)
autoplot(stl.11, xlab="Tempo")

#Componente sazonal da decomposição STL
monthplot(stl.11,col=4,lwd=2,main="Componente Sazonal da Decomposição STL", ylab="Sazonalidade", xlab="Meses")

#Teste de normalidade dos resíduos da decomposição STL
residuos = stl.11$time.series[,3]
shapiro.test(residuos)

#Gráficos dos resíduso da decomposição STL
hist(residuos, main = "Distribuição dos Resíduos da Decomposição STL", xlab = "Resíduos", ylab="Frequência")
qqnorm(residuos, main="Q-Q Normal dos Resíduos", xlab="Quantis Teóricos", ylab="Quantis da Amostra")
qqline(residuos)

#Separação da base de treinamento e teste
base_treinamento = window(serie, start=c(2003,1), end=c(2019,12))
base_teste = window(serie, start=c(2020,1), end=c(2021,12))
plot(serie, lw=2, col=4, main="Divisão série de treinamento (em azul) e teste (em vermelho)", xlab="Tempo", ylab="AOD")
lines(base_teste,lw=2, col="red")

# MODELO AR
modelo_ar = arima(base_treinamento, order = c(2,0,0))
summary(modelo_ar)
checkresiduals(modelo_ar)
shapiro.test(resid(modelo_ar))
qqnorm(resid(modelo_ar))
qqline(resid(modelo_ar))
acf(resid(modelo_ar))
pacf(resid(modelo_ar))
previsao_modelo_ar = forecast(modelo_ar, h = 24)
plot(previsao_modelo_ar)

# MODELO MA
modelo_ma = arima(base_treinamento, order = c(0,0,1))
summary(modelo_ma)
checkresiduals(modelo_ma)
shapiro.test(resid(modelo_ma))
qqnorm(resid(modelo_ma))
qqline(resid(modelo_ma))
acf(resid(modelo_ma))
pacf(resid(modelo_ma))
previsao_modelo_ma = forecast(modelo_ma, h = 36)
plot(previsao_modelo_ma)

# MODELO ARMA
modelo_arma = arima(base_treinamento, order = c(1, 0, 2))
summary(modelo_arma)
checkresiduals(modelo_arma)
shapiro.test(resid(modelo_arma))
qqnorm(resid(modelo_arma))
qqline(resid(modelo_arma))
acf(resid(modelo_arma))
pacf(resid(modelo_arma))
previsao_modelo_arma = forecast(modelo_arma, h = 36)
plot(previsao_modelo_arma)

# MODELO ARIMA
modelo_arima = arima(base_treinamento, order = c(3, 1, 3))
summary(modelo_arima)
checkresiduals(modelo_arima)
shapiro.test(resid(modelo_arima))
qqnorm(resid(modelo_arima))
qqline(resid(modelo_arima))
acf(resid(modelo_arima))
pacf(resid(modelo_arima))
previsao_modelo_arima = forecast(modelo_arima, h = 36)
plot(previsao_modelo_arima)

# MODELOS SARIMA
sarima1 = arima(base_treinamento, order = c(0,0,0), seasonal=c(2,1,0))
summary(sarima1)
checkresiduals(sarima1)
Box.test(sarima1$residuals, type = "Ljung-Box")

sarima2 = arima(base_treinamento, order = c(0,0,0), seasonal=c(1,1,0))
summary(sarima2)
checkresiduals(sarima2)
Box.test(sarima2$residuals, type = "Ljung-Box")

sarima3 = arima(base_treinamento, order = c(0,0,0), seasonal=c(0,1,0))
summary(sarima3)
checkresiduals(sarima3)
Box.test(sarima3$residuals, type = "Ljung-Box")

sarima4 = arima(base_treinamento, order = c(0,0,1), seasonal=c(2,1,0))
summary(sarima4)
checkresiduals(sarima4)
Box.test(sarima4$residuals, type = "Ljung-Box")

sarima5 = arima(base_treinamento, order = c(0,0,1), seasonal=c(1,1,0))
summary(sarima5)
checkresiduals(sarima5)
Box.test(sarima5$residuals, type = "Ljung-Box")

modelo_sarima = arima(base_treinamento, order = c(0,0,1), seasonal=c(0,1,1))
summary(modelo_sarima)
checkresiduals(modelo_sarima)
shapiro.test(resid(modelo_sarima))
qqnorm(resid(modelo_sarima))
qqline(resid(modelo_sarima))
acf(resid(modelo_sarima))
pacf(resid(modelo_sarima))
previsao_modelo_sarima = forecast(modelo_sarima, h = 24)
plot(base_treinamento)
lines(base_treinamento-modelo_sarima$residuals, col="red")
plot(previsao_modelo_sarima)
lines(base_teste, col="red")

# FUNÇÃO AUTO-ARIMA
modelo_auto = auto.arima(base_treinamento, trace = T, stepwise = F, approximation = F)
summary(modelo_auto)
checkresiduals(modelo_auto)
Box.test(modelo_auto$residuals, type = "Ljung-Box")
plot(residuals(modelo_auto), main="Resíduos do Modelo SARIMA (0,0,0)(2,1,0)[12]", xlab = "Tempo", ylab="AOD", lw=2)
hist(residuals(modelo_auto), main="Distribuição dos Resíduos do Modelo SARIMA (0,0,0)(2,1,0)[12]", xlab = "Resíduos", ylab="Frequência")
qqnorm(residuals(modelo_auto), main="Gráfico Q-Q dos Resíduos do Modelo SARIMA (0,0,0)(2,1,0)[12]", xlab="Quantis Teóricos", ylab="Quantis da Amostra")
qqline(residuals(modelo_auto))
acf(residuals(modelo_auto), main="Função Auto Correlação dos Resíduos")
pacf(residuals(modelo_auto), main="Função Auto Correlação Parcial dos Resíduos", ylab="PACF")
shapiro.test(residuals(modelo_auto))
plot(base_treinamento, main="Modelagem do Modelo SARIMA (0,0,0)(2,1,0)[12]", xlab="Tempo", ylab="AOD")
lines(base_treinamento-modelo_auto$residuals, col="blue")

# Previsões
previsoes1 = forecast(modelo_auto, h = 24)
plot(previsoes1, main="Previsões do Modelo SARIMA(0,0,0)(2,1,0)[12]", xlab="Tempo", ylab="AOD")
previsoes1 = as.data.frame(previsoes1)
previsoes1 = previsoes1[ ,1]

previsoes2 = forecast(sarima2, h = 24)
plot(previsoes2, main="Previsões do Modelo SARIMA(0,0,0)(0,1,1)[12]", xlab="Tempo", ylab="AOD")
previsoes2 = as.data.frame(previsoes2)
previsoes2 = previsoes2[ ,1]
previsoes2

previsoes3 = forecast(sarima3, h = 24)
plot(previsoes3, main="Previsões do Modelo SARIMA(0,0,0)(0,1,0)[12]", xlab="Tempo", ylab="AOD")
previsoes3 = as.data.frame(previsoes3)
previsoes3 = previsoes3[ ,1]
previsoes3

previsoes4 = forecast(sarima4, h = 24)
plot(previsoes4, main="Previsões do Modelo SARIMA(0,0,1)(2,1,0)[12]", xlab="Tempo", ylab="AOD")
previsoes4 = as.data.frame(previsoes4)
previsoes4 = previsoes4[ ,1]
previsoes4

previsoes5 = forecast(sarima5, h = 24)
plot(previsoes5, main="Previsões do Modelo SARIMA(0,0,1)(1,1,0)[12]", xlab="Tempo", ylab="AOD")
previsoes5 = as.data.frame(previsoes5)
previsoes5 = previsoes5[ ,1]
previsoes5

#TESTES
teste = as.vector(base_teste)
teste

#Erro médio
me1 = mean(previsoes1 - teste)
me2 = mean(previsoes2 - teste)
me3 = mean(previsoes3 - teste)
me4 = mean(previsoes4 - teste)
me5 = mean(previsoes5 - teste)
me1
me2
me3
me4
me5

#Erro absoluto médio
mae1 = mean(abs(previsoes1 - teste))
mae2 = mean(abs(previsoes2 - teste))
mae3 = mean(abs(previsoes3 - teste))
mae4 = mean(abs(previsoes4 - teste))
mae5 = mean(abs(previsoes5 - teste))
mae1
mae2
mae3
mae4
mae5

#Raiz do erro quadrático médio
rmse1 = sqrt(mean((previsoes1 - teste) ** 2))
rmse2 = sqrt(mean((previsoes2 - teste) ** 2))
rmse3 = sqrt(mean((previsoes3 - teste) ** 2))
rmse4 = sqrt(mean((previsoes4 - teste) ** 2))
rmse5 = sqrt(mean((previsoes5 - teste) ** 2))
rmse1
rmse2
rmse3
rmse4
rmse5

#MAPE
mape1 = MAPE(previsoes1,teste)
mape2 = MAPE(previsoes2,teste)
mape3 = MAPE(previsoes3,teste)
mape4 = MAPE(previsoes4,teste)
mape5 = MAPE(previsoes5,teste)
mape1
mape2
mape3
mape4
mape5

#Previsão do melhor modelo
previsao_vencedor = arima(serie, order = c(0,0,0), seasonal=c(2,1,0))
previsao_vencedor = forecast(previsao_vencedor, h=36)
plot(previsao_vencedor, main="Previsão do Modelo SARIMA (0,0,0)(2,1,0)[12] para 36 meses", xlab="Tempo", ylab="AOD")
