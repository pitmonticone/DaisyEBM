
# Aumentiamo la risoluzione spaziale latitudinale (x5)
library(tidyverse)
library(plotly)

Func <-  function(x){
  0.7768699*cos(0.0164348*x)^2+0.4617747
}
gauss <-  function(x,m,sd,b){
  ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
}
gauss2 <- function(x){
  1-((85)/(sqrt(2*pi*35^2)))*exp(-(x)^2/(35^2))
}

#EBM <-function(S,A,B,K,a_ice,a,T_c
#Input Data
S <- 1370
K <- 3.86
A <- 204
B <- 2.17
sigma <- 5.67*10^(-8)
ai <- 0.62
ab <- 0.25
Temp <- rep(0,360)

Zones <- seq(-89.75, 89.75, by = 0.5)
cosZones <- abs(cospi(Zones/180))
SunWt <- Func(Zones)
Ti <- gauss(Zones,0,50,31.6)
T <- gauss(Zones,0,50,31.6)
Tl <- -10
Tu <- 20
m <- gauss2(Zones)
Rin <- sigma*T^4*(1-m*tanh(19*T^6*10^(-16)))
alb <- function(x,a,b){
  ifelse(x>=Tl & x<Tu,a-((x-Tl)/(Tu-Tl))*(a-b),ifelse(x<Tl,a,b))}
a <- alb(T,ai,ab)

#Test Run

for(i in c(1:360)) {Tcos <- cosZones*T
Tm <- sum(Tcos)/sum(cosZones)
T <- (Rin*(1-a)+K*Tm-A) / (B+K)
a <- alb(T,ai,ab)
Temp[i] <-  Tm
} 

