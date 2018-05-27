rm(list = ls())
library(TSA)
library(forecast)
library(lmtest)
library(fGarch)
library(readr)

# This  function sort the AIC and BIC accoring to their score
sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
} 

# This function produce the out put for residual analysis
residual.analysis <- function(model, std = TRUE){
  library(TSA)
  library(FitAR)
  if (std == TRUE){
    res.model = rstandard(model)
  }else{
    res.model = residuals(model)
  }
  par(mfrow=c(3,2))
  plot(res.model,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
  abline(h=0)
  hist(res.model,main="Histogram of standardised residuals")
  qqnorm(res.model,main="QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  acf(res.model,main="ACF of standardised residuals")
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = length(model$residuals)-1 , StartLag = k + 1, k = 0, SquaredQ = FALSE)
}

# This code read the data set
data <- read_csv("C:/Users/ndnad/Desktop/Time Series Analysis/Project/sunspotnumbers.csv", 
                               col_types = cols(Apr = col_number(), 
                                                           Aug = col_number(), Dec = col_number(), 
                                                           Feb = col_number(), Jan = col_number(), 
                                                           Jul = col_number(), Jun = col_number(), 
                                                           Mar = col_number(), May = col_number(), 
                                                           Nov = col_number(), Oct = col_number(), 
                                                          Sep = col_number()))

# This code change the data frame to time seris 
data.ts = as.vector(t(data[,-1]))
data.ts = ts(data.ts,start=c(1998,1), end=c(2016,12), frequency=12)
class(data.ts)

# plot  of time series plot
plot(data.ts,ylab='sun spots no',xlab='Year',type='o', 
  
                                            main = "Time series plot of monthy sun spot number ")


# seasonility check (can you help me here)
plot(data.ts, type = "l", ylab='sun spot',main = "Time series plot.")
points(y=data.ts,x=time(data.ts), pch=as.vector(season(data.ts)))  # this code is not working 

# scatter plot
plot(y=data.ts,x=zlag(data.ts),ylab='sun spot', xlab='Previous Year sun spot' , main = "Scatter plot of neighboring sun spots")
# there is ovious highly positive autocorrelations in the series 

# this is not working 
y=data.ts
x = zlag(data.ts)        # Generate first lag of the Spawners series
index = 2:length(x)          
cor(y[index],x[index]) 
# thus the series  is higly auto correlated with the previous year 

par(mfrow=c(1,2))
acf(data.ts, main="The sample ACF of sun spot number series")

pacf(data.ts, main="The sample PACF of sun spot number series")

# The Dickey-Fuller Unit-Root test (ADF test)
ar(diff(data.ts)) # To find the value of lag in the following adfTest() function
#adfTest(data.ts, lags = 11, type = "ct", title = NULL,description = NULL)
adf.test(data.ts, k = 11)
# The p value of 0.7063 > 0.5, tells us we can not reject the null hypothesis of data is non stationary. 
# with this unit roots test we are conforming that the series is *seasonal nonstationary*. We need to look for 
 # Seasonality and existence of trend are apparent from the ACF and PACF plots

#--------------- seasonal arima model----------------------------

# First fit a plain model with only the first seasonal difference with order D = 1 
# and see if we can get rid of the seasonal trend effect
# by inspecting the autocorrelation structure of the residuals.
m1.sunspot = arima(data.ts,order=c(0,0,0),seasonal=list(order=c(0,1,0), period=4))
res.m1 = residuals(m1.sunspot);  
plot(res.m1,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m1, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m1, lag.max = 36, main = "The sample PACF of the residuals")
# from the above time series plot we have not seen any clear trends, however seasonal lags
# are still significant in ACF and lags are in decreasing in PACF.so, we add SARMA(0,1) 
#  see if we get rid of seasonal component.

# adding MA(1) in sasoanl part
m2.sunspot = arima(data.ts,order=c(0,0,0),seasonal=list(order=c(0,1,1), period=4))
res.m2 = residuals(m2.sunspot);  
plot(res.m2,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m2, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m2, lag.max = 36, main = "The sample PACF of the residuals")
# first seasonal lags in ACF is resolved while others higher order are still significant 
# and decresing  patterns of lags in PACF. We think this is due to ordinary series 

# So, we will apply differentiating on the ordinary seris and see if we can see the trend more clearly.
m4.sunspot = arima(data.ts,order=c(0,1,0),seasonal=list(order=c(0,1,1), period=4))
res.m4 = residuals(m4.sunspot);  
plot(res.m4,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res.m4, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res.m4, lag.max = 36, main = "The sample PACF of the residuals")
# we get rid of trends and seaonal effects. Not a single seasonal lags in ACF and PACF are significant 
# for lower lags 

# We are going to see the posible sets of models using eacf and BIC of residulas of the 
# above models

eacf(res.m4)
res = armasubsets(y=data.ts, nar=14,nma=14,y.name='test',ar.method='ols') 
plot(res)


# From the EACF, we will include AR(1) order as well.
# SARIMA(0,1,2)x(0,1,1)_4
# SARIMA(0,1,3)x(0,1,1)_4 
# SARIMA(1,1,2)x(0,1,1)_4 and
# SARIMA(2,1,1)x(0,1,1)_4 will be fitted

# from BIC table we wil include AR(3)
# SARIMA(1,1,0)x(0,1,1)_4
# SARIMA(3,1,0)x(0,1,1)_4

m4_012.sunspot = arima(data.ts ,order=c(0,1,2),seasonal=list(order=c(0,1,1), period=4),method = "ML")
res_012 = residuals(m4_012.sunspot);  
plot(res_012,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res_012, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res_012, lag.max = 36, main = "The sample PACF of the residuals")
# still have some significant lags in ACF and PACF

m4_013.sunspot = arima(data.ts ,order=c(0,1,3),seasonal=list(order=c(0,1,1), period=4),method = "ML")
res_013 = residuals(m4_013.sunspot);  
plot(res_013,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res_013, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res_013, lag.max = 36, main = "The sample PACF of the residuals")
# still have some significant lags in ACF and PACF

m4_112.sunspot = arima(data.ts ,order=c(1,1,2),seasonal=list(order=c(0,1,1), period=4),method = "ML")
res_112 = residuals(m4_112.sunspot);  
plot(res_112,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res_112, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res_112, lag.max = 36, main = "The sample PACF of the residuals")

m4_211.sunspot = arima(data.ts ,order=c(2,1,1),seasonal=list(order=c(0,1,1), period=4),method = "ML")
res_211 = residuals(m4_211.sunspot);  
plot(res_211,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res_211, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res_211, lag.max = 36, main = "The sample PACF of the residuals")

# from BIC table ma(3)
m4_110.sunspot = arima(data.ts ,order=c(1,1,0),seasonal=list(order=c(0,1,1), period=4),method = "ML")
res_110 = residuals(m4_110.sunspot);  
plot(res_110,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res_110, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res_110, lag.max = 36, main = "The sample PACF of the residuals")

m4_310.sunspot = arima(data.ts ,order=c(3,1,0),seasonal=list(order=c(0,1,1), period=4),method = "ML")
res_310 = residuals(m4_310.sunspot);  
plot(res_310,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
par(mfrow=c(1,2))
acf(res_310, lag.max = 36, main = "The sample ACF of the residuals")
pacf(res_310, lag.max = 36, main = "The sample PACF of the residuals")


coeftest(m4_012.sunspot) # all significant 
coeftest(m4_013.sunspot) # MA(3) is not significant 
coeftest(m4_112.sunspot) # MA(1) and AR(1) is not significant
coeftest(m4_110.sunspot) # all significant
coeftest(m4_211.sunspot) # AR(1) is not significant
coeftest(m4_310.sunspot) # all are signifcant 


sc.AIC=AIC(m4_012.sunspot, m4_013.sunspot, m4_112.sunspot,m4_110.sunspot, m4_211.sunspot, m4_310.sunspot)
sc.BIC=BIC(m4_012.sunspot, m4_013.sunspot, m4_112.sunspot,m4_110.sunspot, m4_211.sunspot, m4_310.sunspot)

sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "bic")
# using coef test, AIC and BIC m4_012.sunpost or [sarima(0,1,2)X(0,,1)_4] is best model 

residual.analysis(model = m4_012.sunspot) # good but still have some lags signicant
residual.analysis(model = m4_013.sunspot) # still good but proble with lag 12
residual.analysis(model = m4_112.sunspot) # still good but problem with lag 12
residual.analysis(model = m4_110.sunspot) # good
residual.analysis(model = m4_211.sunspot) # still good but slight significan line in ACF
residual.analysis(model = m4_310.sunspot) # still good but slight significan line in ACF

m4_012.sunspot = Arima(data.ts, order=c(0,1,2), seasonal=list(order=c(0,1,1), period=4), method = "ML")
preds1 = forecast(m4_012.sunspot, h = 24)
plot(preds1)

m4_310.sunspot = Arima(data.ts, order=c(3,1,0), seasonal=list(order=c(0,1,1), period=4), method = "ML")
preds2 = forecast(m4_310.sunspot, h = 24)
plot(preds2)

# choose either of above, I think (m4_012.sunspot) is better model because this modle is quite small than 
# the (m4_310.sunspot)


