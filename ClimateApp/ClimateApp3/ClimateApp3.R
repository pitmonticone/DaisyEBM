library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(broom)
library(scales)
library(shinythemes)

# INITIALIZATION ----

### Random 
gauss <-  function(x,m,sd,b){
  ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
}

### Solar Luminosity (Latitude)
Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
### Solar Luminosity  (Time)
Sun1 <- function(x,a){a*x/100} 
Sun2 <- function(x){1370*(sinpi((x+90)/180))^2}
Sun3 <- function(x){1370-(((1370)/(sqrt(2*pi*1^2)))*exp(-(x-50)^2/(2*1^2)))}
Sun4 <- function(x){ifelse(x==50,1370/3,1370)}
Sun5 <- function(x){(1/100) * (abs(x-150)+25)}
Sun6 <- function(x){1370*(1+0.1*cospi(x/180))}
Incident <- function(x,y){x*y/4}

# Albedo
Step <- function(x,c){ifelse(x<c, 0.6, 0.3)}
alb <- function(x,a,b,c,d){
  (exp(c*(x+d)) / (exp(c*(x+d))+1)) * (b-a) + a 
}

# EBM
ebm01 <- function(cycles,S,A,B,K,ai,ab,gamma,delta) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  t <- c(1:cycles)
  Temperature <- rep(0,cycles)
  Ti <- gauss(Zones,0,50,31.6)
  SunWt <- Func(Zones)
  Rin <- Incident(S,SunWt)
  T <- Ti
  a <- alb(T,ai,ab,gamma,delta)
  
  for(i in t)  { Tcos <- cosZones*T
  Tm <- sum(Tcos)/sum(cosZones)
  T <- (Rin*(1-a)+K*Tm-A) / (B+K)
  a <- alb(T,ai,ab,gamma,delta)
  Temperature[i] <-  Tm  }
  
  return( data.frame(Zones,T,a,Ti)) } 
ebm02 <- function(cycles,S,A,B,K,ai,ab,gamma,delta) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  t <- c(1:cycles)
  Temperature <- rep(0,cycles)
  Ti <- gauss(Zones,0,50,31.6)
  SunWt <- Func(Zones)
  Rin <- Incident(S,SunWt)
  T <- Ti
  a <- alb(T,ai,ab,gamma,delta)
  for(i in t) 
  {Tcos <- cosZones*T
  Tm <- sum(Tcos)/sum(cosZones)
  T <- (Rin*(1-a)+K*Tm-A) / (B+K)
  a <- alb(T,ai,ab,gamma,delta)
  Temperature[i] <-  Tm  }
  return( data.frame(t,Temperature) ) } 

ebm11 <- function(cycles1,cycles2,S,A,B,K,ai,ab,gamma,delta) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun1 <- function(x,a){a*x/100} 
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  J <- rep(0,cycles1)
  TEMP1 <- matrix(NA, nrow=90, ncol=cycles1)
  Temp <- rep(0,cycles1)
  SunWt <- Func(Zones)
  
  for(j in c(1:cycles1)){
    T <- gauss(Zones,0,50,31.6)
    a <- alb(T,ai,ab,gamma,delta)
    S1 <- Sun1(S,j)
    Rin <- Incident(S1,SunWt) 
    
    for(i in c(1:cycles2))
    {Tcos <- cosZones*T
    Tm <- sum(Tcos)/sum(cosZones)
    T <- (Rin*(1-a)+K*Tm-A) / (B+K)
    a <- alb(T,ai,ab,gamma,delta)
    TEMP1[,j] <- T
    J[j] <- j
    Temp[j] <- Tm  } } 
  
  return( data.frame(J,Temp) )   } 
ebm12 <- function(cycles1,cycles2,S,A,B,K,ai,ab,gamma,delta) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun1 <- function(x,a){a*x/100}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  J <- rep(0,cycles1)
  TEMP1 <- matrix(NA, nrow=90, ncol=cycles1)
  Temp <- rep(0,cycles1)
  
  
  for(j in c(1:cycles1)){
    T <- gauss(Zones,0,50,31.6)
    a <- alb(T,ai,ab,gamma,delta)
    S1 <- Sun1(S,j)
    Rin <- Incident(S1,SunWt)
    
    for(i in c(1:cycles2))
    {Tcos <- cosZones*T
    Tm <- sum(Tcos)/sum(cosZones)
    T <- (Rin*(1-a)+K*Tm-A) / (B+K)
    a <- alb(T,ai,ab,gamma,delta)
    TEMP1[,j] <- T
    J[j] <- j
    Temp[j] <- Tm  } } 
  
  return( TEMP1 ) }

ebm21 <- function(cycles1,cycles2,S,A,B,K,ai,ab,gamma,delta ) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun2 <- function(x){S*(sinpi((x+90)/180))^2}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  J <- rep(0,cycles1)
  Temp2 <- rep(0,cycles1)
  T2 <- gauss(Zones,0,50,31.6)
  TEMP2 <- matrix(NA, nrow=90, ncol=cycles1)
  a <- alb(T2,ai,ab,gamma,delta)
  Sarr <- rep(0,cycles1) 
  
  for(j in c(1:cycles1)){
    S2 <- Sun2(j)
    Rin <- Incident(S2,SunWt)
    for(i in c(1:cycles2))
    {Tcos <- cosZones*T2
    Tm <- sum(Tcos)/sum(cosZones)
    T2 <- (Rin*(1-a)+K*Tm-A) / (B+K)
    a <- alb(T2,ai,ab,gamma,delta)
    }
    Sarr[j] <- Sun2(j)
    TEMP2[,j] <- T2
    J[j] <- j
    Temp2[j] <- Tm  }  
  
  return( data.frame(J,Temp2) ) }
ebm22 <- function(cycles1,cycles2,S,A,B,K,ai,ab,gamma,delta ) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun2 <- function(x){S*(sinpi((x+90)/180))^2}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  J <- rep(0,cycles1)
  Temp2 <- rep(0,cycles1)
  T2 <- gauss(Zones,0,50,31.6)
  TEMP2 <- matrix(NA, nrow=90, ncol=cycles1)
  a <- alb(T2,ai,ab,gamma,delta)
  Sarr <- rep(0,cycles1) 
  SunWt <- Func(Zones)
  
  for(j in c(1:cycles1)){
    S2 <- Sun2(j)
    Rin <- Incident(S2,SunWt)
    for(i in c(1:cycles2))
    {Tcos <- cosZones*T2
    Tm <- sum(Tcos)/sum(cosZones)
    T2 <- (Rin*(1-a)+K*Tm-A) / (B+K)
    a <- alb(T2,ai,ab,gamma,delta)
    }
    Sarr[j] <- Sun2(j)
    TEMP2[,j] <- T2
    J[j] <- j
    Temp2[j] <- Tm  }  
  
  return( TEMP2 )  }
ebm23 <- function(cycles1,cycles2,S,A,B,K,ai,ab,gamma,delta ) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun2 <- function(x){S*(sinpi((x+90)/180))^2}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  J <- rep(0,cycles1)
  Temp2 <- rep(0,cycles1)
  T2 <- gauss(Zones,0,50,31.6)
  TEMP2 <- matrix(NA, nrow=90, ncol=cycles1)
  a <- alb(T2,ai,ab,gamma,delta)
  Sarr <- rep(0,cycles1) 
  SunWt <- Func(Zones)
  
  for(j in c(1:cycles1)){
    S2 <- Sun2(j)
    Rin <- Incident(S2,SunWt)
    for(i in c(1:cycles2))
    {Tcos <- cosZones*T2
    Tm <- sum(Tcos)/sum(cosZones)
    T2 <- (Rin*(1-a)+K*Tm-A) / (B+K)
    a <- alb(T2,ai,ab,gamma,delta)
    }
    Sarr[j] <- Sun2(j)
    TEMP2[,j] <- T2
    J[j] <- j
    Temp2[j] <- Tm  }  
  
  return( data.frame(Sarr,Temp2) ) }

ebm31 <- function(cycles1,cycles2,S,A,B,K,ai,ab,gamma,delta ) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun3 <- function(x){S-(((S)/(sqrt(2*pi*1^2)))*exp(-(x-50)^2/(2*1^2)))}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  J <- rep(0,cycles1)
  TEMP3 <- matrix(NA, nrow=90, ncol=cycles1)
  Temp <- rep(0,cycles1)
  T3 <- gauss(Zones,0,50,31.6)
  a <- alb(T3,ai,ab,gamma,delta)
  
  for(j in c(1:cycles1)){
    S3 <- Sun3(j)
    Rin <- Incident(S3,SunWt)
    for(i in c(1:cycles2))
    {Tcos <- cosZones*T3
    Tm <- sum(Tcos)/sum(cosZones)
    T3 <- (Rin*(1-a)+K*Tm-A) / (B+K)
    a <- alb(T3,ai,ab,gamma,delta)
    }
    TEMP3[,j] <- T3
    J[j] <- j
    Temp[j] <- Tm
  }
  
  return( data.frame(J,Temp) ) }
ebm32 <- function(cycles1,cycles2,S,A,B,K,ai,ab,gamma,delta ) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun3 <- function(x){S-(((S)/(sqrt(2*pi*1^2)))*exp(-(x-50)^2/(2*1^2)))}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  J <- rep(0,cycles1)
  TEMP3 <- matrix(NA, nrow=90, ncol=cycles1)
  Temp <- rep(0,cycles1)
  T3 <- gauss(Zones,0,50,31.6)
  a <- alb(T3,ai,ab,gamma,delta)
  
  for(j in c(1:cycles1)){
    S3 <- Sun3(j)
    Rin <- Incident(S3,SunWt)
    for(i in c(1:cycles2))
    {Tcos <- cosZones*T3
    Tm <- sum(Tcos)/sum(cosZones)
    T3 <- (Rin*(1-a)+K*Tm-A) / (B+K)
    a <- alb(T3,ai,ab,gamma,delta)
    }
    TEMP3[,j] <- T3
    J[j] <- j
    Temp[j] <- Tm
  }
  
  return( TEMP3 ) }

ebm41 <- function(cycles1,cycles2,S,A,B,K,ai,ab,gamma,delta ) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun4 <- function(x){ifelse(x==50,S/3,S)}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  J <- rep(0,cycles1)
  Temp <- rep(0,cycles1)
  T4 <- gauss(Zones,0,50,31.6)
  TEMP4 <- matrix(NA, nrow=90, ncol=cycles1)
  a <- alb(T4,ai,ab,gamma,delta)
  
  for(j in c(1:cycles1)){
    S4 <- Sun4(j)
    Rin <- Incident(S4,SunWt)
    for(i in c(1:cycles2))
    {Tcos <- cosZones*T4
    Tm <- sum(Tcos)/sum(cosZones)
    T4 <- (Rin*(1-a)+K*Tm-A) / (B+K)
    a <- alb(T4,ai,ab,gamma,delta)
    }
    TEMP4[,j] <- T4
    J[j] <- j
    Temp[j] <- Tm
  }
  
  return( data.frame(J,Temp) ) }
ebm42 <- function(cycles1,cycles2,S,A,B,K,ai,ab,gamma,delta ) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun4 <- function(x){ifelse(x==50,S/3,S)}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  J <- rep(0,cycles1)
  Temp <- rep(0,cycles1)
  T4 <- gauss(Zones,0,50,31.6)
  TEMP4 <- matrix(NA, nrow=90, ncol=cycles1)
  a <- alb(T4,ai,ab,gamma,delta)
  
  for(j in c(1:cycles1)){
    S4 <- Sun4(j)
    Rin <- Incident(S4,SunWt)
    for(i in c(1:cycles2))
    {Tcos <- cosZones*T4
    Tm <- sum(Tcos)/sum(cosZones)
    T4 <- (Rin*(1-a)+K*Tm-A) / (B+K)
    a <- alb(T4,ai,ab,gamma,delta)
    }
    TEMP4[,j] <- T4
    J[j] <- j
    Temp[j] <- Tm
  }
  return( TEMP4 ) }

ebm51 <- function(cycles1,cycles2,A,B,K,ai,ab,gamma,delta ) {
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun5 <- function(x){(1/100) * (abs(x-150)+25)}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  J <- rep(0,cycles1)
  Temp5 <- rep(0,cycles1)
  T5 <- gauss(Zones,0,50,31.6)
  TEMP5 <- matrix(NA, nrow=90, ncol=cycles1)
  a <- alb(T5,ai,ab,gamma,delta)
  Sarr <- rep(0,cycles1) 
  
  for(j in c(0:cycles1)){
    S <- Sun5(j)
    Rin <- 1370*Incident(S,SunWt)
    for(i in c(1:cycles2))
    {Tcos <- cosZones*T5
    Tm <- sum(Tcos)/sum(cosZones)
    T5 <- (Rin*(1-a)+K*Tm-A) / (B+K)
    a <- alb(T5,ai,ab,gamma,delta)
    }
    Sarr[j] <- Sun5(j)
    J[j] <- j
    Temp5[j] <- Tm
  }
  return( data.frame(Sarr,Temp5) ) }

ebm_ND1 <- function(cycles,S,w0,b0,A,B,K,ai,ab,aW,aB,gamma,delta){
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun6 <- function(x){S*(1+0.1*cospi(x/180))}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  T <- gauss(Zones,0,50,31.6)-6
  
  w <- rep(w0,length(Zones)) #0.5
  b <- rep(b0,length(Zones)) #0.2
  u <- rep(1-w0-b0,length(Zones))
  
  a <- w*aW+b*aB+u*alb(T,ai,ab,gamma,delta)
  
  Barr <- rep(0,cycles)
  Warr <- rep(0,cycles)
  Uarr <- rep(0,cycles)
  Tarr <- rep(0,cycles)
  I <- rep(0,cycles)
  
  TEMP <- matrix(NA, nrow=90, ncol=cycles)
  
  for(i in c(1:cycles)) {
    S6 <- Sun6(i)  # oppure costante S <- 1370
    Rin <- Incident(S6,SunWt)
    Tcos <- cosZones*T
    Tm <- sum(Tcos)/sum(cosZones)
    T <- (Rin*(1-a)+K*Tm-A) / (B+K)
    TEMP[,i] <- T
    a <- alb(T,ai,ab,gamma,delta)
    I[i] <- i
    Tarr[i] <- T[45]
  } 
  
  return(data.frame(Zones,w,b,u,T) )
}
ebm_ND2 <- function(cycles,S,w0,b0,A,B,K,ai,ab,aW,aB,gamma,delta){
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun6 <- function(x){S*(1+0.1*cospi(x/180))}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  T <- gauss(Zones,0,50,31.6)-6
  
  w <- rep(w0,length(Zones)) #0.5
  b <- rep(b0,length(Zones)) #0.2
  u <- rep(1-w0-b0,length(Zones))
  
  a <- w*aW+b*aB+u*alb(T,ai,ab,gamma,delta)
  
  Barr <- rep(0,cycles)
  Warr <- rep(0,cycles)
  Uarr <- rep(0,cycles)
  Tarr <- rep(0,cycles)
  I <- rep(0,cycles)
  
  TEMP <- matrix(NA, nrow=90, ncol=cycles)
  
  for(i in c(1:cycles)) {
    S6 <- Sun6(i)  # oppure costante S <- 1370
    Rin <- Incident(S6,SunWt)
    Tcos <- cosZones*T
    Tm <- sum(Tcos)/sum(cosZones)
    T <- (Rin*(1-a)+K*Tm-A) / (B+K)
    TEMP[,i] <- T
    a <- alb(T,ai,ab,gamma,delta)
    I[i] <- i
    Tarr[i] <- T[45]
  } 
  return( data.frame(I,Barr,Warr,Uarr,Tarr) )}
ebm_ND3 <- function(cycles,S,w0,b0,A,B,K,ai,ab,aW,aB,gamma,delta){
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun6 <- function(x){S*(1+0.1*cospi(x/180))}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  T <- gauss(Zones,0,50,31.6)-6
  
  w <- rep(w0,length(Zones)) #0.5
  b <- rep(b0,length(Zones)) #0.2
  u <- rep(1-w0-b0,length(Zones))
  
  a <- w*aW+b*aB+u*alb(T,ai,ab,gamma,delta)
  
  Barr <- rep(0,cycles)
  Warr <- rep(0,cycles)
  Uarr <- rep(0,cycles)
  Tarr <- rep(0,cycles)
  I <- rep(0,cycles)
  
  TEMP <- matrix(NA, nrow=90, ncol=cycles)
  
  for(i in c(1:cycles)) {
    S6 <- Sun6(i)  # oppure costante S <- 1370
    Rin <- Incident(S6,SunWt)
    Tcos <- cosZones*T
    Tm <- sum(Tcos)/sum(cosZones)
    T <- (Rin*(1-a)+K*Tm-A) / (B+K)
    TEMP[,i] <- T
    a <- alb(T,ai,ab,gamma,delta)
    I[i] <- i
    Tarr[i] <- T[45]
  } 
  return( TEMP )}

ebm_D1 <- function(cycles,S,w0,b0,T0,c,k,A,B,K,ai,ab,aW,aB,gamma,delta){
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun6 <- function(x){S*(1+0.1*cospi(x/180))}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  T <- gauss(Zones,0,50,31.6)-6
  
  w <- rep(w0,length(Zones)) #0.5
  b <- rep(b0,length(Zones)) #0.2
  u <- rep(1-w0-b0,length(Zones))
  
  a <- w*aW+b*aB+u*alb(T,ai,ab,gamma,delta)
  
  Barr <- rep(0,cycles)
  Warr <- rep(0,cycles)
  Uarr <- rep(0,cycles)
  Tarr <- rep(0,cycles)
  I <- rep(0,cycles)
  
  TEMP <- matrix(NA, nrow=90, ncol=cycles)
  
  for(i in c(1:cycles)) {
    S6 <- Sun6(i)  # oppure costante S <- 1370
    Rin <- Incident(S6,SunWt)
    Tcos <- cosZones*T
    Tm <- sum(Tcos)/sum(cosZones)
    T <- (Rin*(1-a)+K*Tm-A) / (B+K)
    TEMP[,i] <- T
    Tw <- T+c*(a-aW)
    Tb <- T+c*(a-aB)
    Fw <- 1-k*(T0-Tw)^2
    Fb <- 1-k*(T0-Tb)^2
    for(j in c(1:length(Zones))){
      if(Fw[j]<0){Fw[j]=0}
      if(Fb[j]<0){Fb[j]=0}  }
    w <- w+w*(u*Fw-D)
    b <- b+b*(u*Fb-D)
    for(j in c(1:length(Zones))){
      if(w[j]<0.001){w[j]=0.001}
      if(b[j]<0.001){b[j]=0.001}  }
    u <- 1-w-b
    a <- w*aW+b*aB+u*alb(T,ai,ab,gamma,delta)
    Barr[i] <- b[45]
    Warr[i] <- w[45]
    Uarr[i] <- u[45]
    I[i] <- i
    Tarr[i] <- T[45]
  } 
  return ( data.frame(Zones,w,b,u,T) )
}
ebm_D2 <- function(cycles,S,w0,b0,T0,c,k,A,B,K,ai,ab,aW,aB,gamma,delta){
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }
  Sun6 <- function(x){S*(1+0.1*cospi(x/180))}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  T <- gauss(Zones,0,50,31.6)-6
  
  w <- rep(w0,length(Zones)) #0.5
  b <- rep(b0,length(Zones)) #0.2
  u <- rep(1-w0-b0,length(Zones))
  
  a <- w*aW+b*aB+u*alb(T,ai,ab,gamma,delta)
  
  Barr <- rep(0,cycles)
  Warr <- rep(0,cycles)
  Uarr <- rep(0,cycles)
  Tarr <- rep(0,cycles)
  I <- rep(0,cycles)
  
  TEMP <- matrix(NA, nrow=90, ncol=cycles)
  
  for(i in c(1:cycles)) {
    S6 <- Sun6(i)  # oppure costante S <- 1370
    Rin <- Incident(S6,SunWt)
    Tcos <- cosZones*T
    Tm <- sum(Tcos)/sum(cosZones)
    T <- (Rin*(1-a)+K*Tm-A) / (B+K)
    TEMP[,i] <- T
    Tw <- T+c*(a-aW)
    Tb <- T+c*(a-aB)
    Fw <- 1-k*(T0-Tw)^2
    Fb <- 1-k*(T0-Tb)^2
    for(j in c(1:length(Zones))){
      if(Fw[j]<0){Fw[j]=0}
      if(Fb[j]<0){Fb[j]=0}  }
    w <- w+w*(u*Fw-D)
    b <- b+b*(u*Fb-D)
    for(j in c(1:length(Zones))){
      if(w[j]<0.001){w[j]=0.001}
      if(b[j]<0.001){b[j]=0.001}  }
    u <- 1-w-b
    a <- w*aW+b*aB+u*alb(T,ai,ab,gamma,delta)
    Barr[i] <- b[45]
    Warr[i] <- w[45]
    Uarr[i] <- u[45]
    I[i] <- i
    Tarr[i] <- T[45]
  } 
  return ( data.frame(I,Barr,Warr,Uarr,Tarr) )
}
ebm_D3 <- function(cycles,S,w0,b0,T0,c,k,A,B,K,ai,ab,aW,aB,gamma,delta){
  Incident <- function(x,y){ x*y/4 }
  Func <-  function(x){ 0.7768699*cos(0.0164348*x)^2+0.4617747 }  
  Sun6 <- function(x){S*(1+0.1*cospi(x/180))}
  gauss <-  function(x,m,sd,b){
    ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
  }
  
  Zones <- seq(-89, 89, by = 2)
  cosZones <- abs(cospi(Zones/180))
  SunWt <- Func(Zones)
  T <- gauss(Zones,0,50,31.6)-6
  
  w <- rep(w0,length(Zones)) #0.5
  b <- rep(b0,length(Zones)) #0.2
  u <- rep(1-w0-b0,length(Zones))
  
  a <- w*aW+b*aB+u*alb(T,ai,ab,gamma,delta)
  
  Barr <- rep(0,cycles)
  Warr <- rep(0,cycles)
  Uarr <- rep(0,cycles)
  Tarr <- rep(0,cycles)
  I <- rep(0,cycles)
  
  TEMP <- matrix(NA, nrow=90, ncol=cycles)
  
  for(i in c(1:cycles)) {
    S6 <- Sun6(i) 
    Rin <- Incident(S6,SunWt)
    Tcos <- cosZones*T
    Tm <- sum(Tcos)/sum(cosZones)
    T <- (Rin*(1-a)+K*Tm-A) / (B+K)
    TEMP[,i] <- T
    Tw <- T+c*(a-aW)
    Tb <- T+c*(a-aB)
    Fw <- 1-k*(T0-Tw)^2
    Fb <- 1-k*(T0-Tb)^2
    for(j in c(1:length(Zones))){
      if(Fw[j]<0){Fw[j]=0}
      if(Fb[j]<0){Fb[j]=0}  }
    w <- w+w*(u*Fw-D)
    b <- b+b*(u*Fb-D)
    for(j in c(1:length(Zones))){
      if(w[j]<0.001){w[j]=0.001}
      if(b[j]<0.001){b[j]=0.001}  }
    u <- 1-w-b
    a <- w*aW+b*aB+u*alb(T,ai,ab,gamma,delta)
    Barr[i] <- b[45]
    Warr[i] <- w[45]
    Uarr[i] <- u[45]
    I[i] <- i
    Tarr[i] <- T[45]
  } 
  return ( TEMP )
}


# UI ----
ui <- dashboardPage(  
  dashboardHeader(title = "Climate Modelling"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Widgets", tabName = "widgets", icon = icon("th"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
              fluidRow(
                box(
                  title = "Input Parameters",
                  numericInput("S",
                               label = "Solar Luminosity (\\(\\ S \\))", 
                               value = 1370,
                               min = NA, max = NA, step = NA),
                  numericInput("A",
                               label = "\\(\\ A \\)", 
                               value = 204,
                               min = NA, max = NA, step = NA),
                  numericInput("B",
                               label = "\\(\\ B \\)", 
                               value = 2.17,
                               min = NA, max = NA, step = NA),
                  numericInput("K",
                               label = "\\(\\ K \\)", 
                               value = 3.86,
                               min = NA, max = NA, step = NA),
                  numericInput("ai",
                               label = "Ice Albedo (\\(\\alpha\\))", 
                               value = 0.62,
                               min = NA, max = NA, step = NA),
                  numericInput("ab",
                               label = "Bare Ground Albedo (\\(\\beta\\))", 
                               value = 0.25,
                               min = NA, max = NA, step = NA),
                  numericInput("gamma",
                               label = "Steepness Albedo (\\(\\gamma\\))", 
                               value = 2.2,
                               min = NA, max = NA, step = NA),
                  numericInput("delta",
                               label = "Offset Albedo (\\(\\delta\\))", 
                               value = 10/2.2,
                               min = NA, max = NA, step = NA),
                  numericInput("aW",
                               label = "White Albedo (\\(\\ a_{w} \\))", 
                               value = 0.75,
                               min = NA, max = NA, step = NA),
                  numericInput("aB",
                               label = "Black Albedo (\\(\\ a_{b} \\))", 
                               value = 0.25,
                               min = NA, max = NA, step = NA),
                  numericInput("w0",
                               label = "% White (\\(\\ w_{0} \\))", 
                               value = 0.5,
                               min = NA, max = NA, step = NA),
                  numericInput("b0",
                               label = "% Black (\\(\\ b_{0} \\))", 
                               value = 0.25,
                               min = NA, max = NA, step = NA),
                  numericInput("u0",
                               label = "% Bare Ground (\\(\\ 1-w_0-b_0 \\))", 
                               value = 1-0.5-0.25,
                               min = NA, max = NA, step = NA),
                  numericInput("c",
                               label = "\\(\\ c \\)", 
                               value = 7,
                               min = NA, max = NA, step = NA),
                  numericInput("k",
                               label = "\\(\\ k \\)", 
                               value = 0.003265*0.75,
                               min = NA, max = NA, step = NA),
                  numericInput("T0",
                               label = "\\(\\ T_0 \\)", 
                               value = 20,
                               min = NA, max = NA, step = NA),
                  numericInput("D",
                               label = "Death Rate (\\(\\ D \\))", 
                               value = 0.3,
                               min = NA, max = NA, step = NA),
                  br(),
                  br(),
                  h3("Source :"),
                  a("Climate Modelling Project", href="https://pitmonticone.github.io/Climate-Physics/"),
                  h3("Authors :"),
                  tags$ul( tags$li( a("Pietro Monticone", href="https://github.com/pitmonticone") ), 
                           tags$li( a("Davide Orsenigo", href="https://github.com/dadorse") )
                  )
                )
              )
      ),
      
      # Second tab content
      tabItem(tabName = "widgets",
              h2("Widgets tab content")
      )
    )
  )
)

# SERVER ----
server <- function(input, output) { }

shinyApp(ui, server)