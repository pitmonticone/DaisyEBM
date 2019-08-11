library(shiny)
library(tidyverse)
library(plotly)
library(broom)
library(scales)
library(shinythemes) #---> themeSelector()
# risonanza stocastica
# a(input$T,input$a,input$b,input$c) 
# LATEX FORMULAS
# COSIN^2 --> COSIN 
# EBM GENERALIZZATO CON FUNZIONI (COMPATTIFICAZIONE)
# GRID PLOTLY WITH DAISIES


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
ui <- fluidPage(
  theme = shinytheme("spacelab"),
  titlePanel("Climate Modelling"),
  withMathJax(),
  
  sidebarLayout(
    sidebarPanel(
      h3("Input Parameters"),
      
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
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Albedo",
                 uiOutput("Albedo"),
                 plotOutput("plot_albedo")
        ),
        tabPanel("Run 0", 
                 uiOutput("Run0"),
                 plotOutput("plot1_Run0"),plotOutput("plot2_Run0")),
        tabPanel("Run 1", 
                 uiOutput("Run1"),
                 plotOutput("plot1_Run1"),plotlyOutput("plot2_Run1")),
        tabPanel("Run 2", 
                 uiOutput("Run2"),
                 plotOutput("plot1_Run2"),plotlyOutput("plot2_Run2")),
        tabPanel("Run 3 & 4",
                 uiOutput("Run34"),
                 plotOutput("plot1_Run3"),plotOutput("plot1_Run4"),plotlyOutput("plot2_Run3"),plotlyOutput("plot2_Run4")),
        tabPanel("Hysteresis Cycles",
                 uiOutput("Hysteresis"),
                 plotOutput("plot1_Hysteresis"), plotOutput("plot2_Hysteresis")),
        tabPanel("Without Daisies",
                 uiOutput("NoDaisy"),
                 plotOutput("plot1_noDaisy"), plotOutput("plot2_noDaisy"), plotlyOutput("plot3_noDaisy")),
        tabPanel("With Daisies",
                 uiOutput("Daisy"),
                 plotOutput("plot1_Daisy"), plotOutput("plot2_Daisy"), plotlyOutput("plot3_Daisy") )
      )
    )
  )
)



# SERVER ----
server <- function(input, output,session){
  
  output$Albedo <- renderUI({
    withMathJax(
      p('$$a(T)=\\frac{e^{\\gamma(T+\\delta)}}{e^{\\gamma(T+\\delta)}+1} (\\beta-\\alpha) + \\alpha$$')
    )
  })
  
  output$Run0 <- renderUI({
    withMathJax(
      p('$$a(T)=\\frac{e^{\\gamma(T+\\delta)}}{e^{\\gamma(T+\\delta)}+1} (\\beta-\\alpha) + \\alpha$$'),
      p('$$ T =\\frac{ R_{in}(1-a(T))+K \\bar{T}-A }{ B+K }$$')
      
    )
  })
  
  output$Run1 <- renderUI({
    withMathJax(
      p('$$ S_1(t)= \\frac{S}{100}t $$'),
      p('$$ a(T)=\\frac{e^{\\gamma(T+\\delta)}}{e^{\\gamma(T+\\delta)}+1} (\\beta-\\alpha) + \\alpha$$'),
      p('$$ T =\\frac{ R_{in}(1-a(T))+K \\bar{T}-A }{ B+K }$$')
    )
  })
  
  output$Run2 <- renderUI({
    withMathJax(
      p('$$ S_2(t)= S(\\sin^2(t+90)) $$'),
      p('$$ a(T)=\\frac{e^{\\gamma(T+\\delta)}}{e^{\\gamma(T+\\delta)}+1} (\\beta-\\alpha) + \\alpha$$'),
      p('$$ T =\\frac{ R_{in}(1-a(T))+K \\bar{T}-A }{ B+K }$$')
    )
  })
  
  output$Run34 <- renderUI({
    withMathJax(
      p('$$ S_3(t)= S\\Big(1-\\frac{1}{\\sqrt{2\\pi}}e^{-\\frac{(t-50)^2}{2}}\\Big) $$'),
      p('$$ S_4(t)= S\\Big(1-\\frac{1}{3}\\delta(t-50)\\Big) $$'),
      p('$$ a(T)=\\frac{e^{\\gamma(T+\\delta)}}{e^{\\gamma(T+\\delta)}+1} (\\beta-\\alpha) + \\alpha$$'),
      p('$$ T =\\frac{ R_{in}(1-a(T))+K \\bar{T}-A }{ B+K }$$')
    )
  })
  
  output$Hysteresis <- renderUI({
    withMathJax(
      p('$$ S_2(t)= S(\\sin^2(t+90)) $$'),
      p('$$ S_5(t)=\\frac{1}{100}(|t-150|+25) $$'),
      p('$$ a(T)=\\frac{e^{\\gamma(T+\\delta)}}{e^{\\gamma(T+\\delta)}+1} (\\beta-\\alpha) + \\alpha$$'),
      p('$$ T =\\frac{ R_{in}(1-a(T))+K \\bar{T}-A }{ B+K }$$')
    )
  })
  
  output$NoDaisy <- renderUI({
    withMathJax(
      p('$$ S_6(t)=S\\Big(1+\\frac{1}{10} cos(t)\\Big)  $$'),
      p('$$ a(T)=\\frac{e^{\\gamma(T+\\delta)}}{e^{\\gamma(T+\\delta)}+1} (\\beta-\\alpha) + \\alpha$$'),
      p('$$ T =\\frac{ R_{in}(1-a(T))+K \\bar{T}-A }{ B+K }$$')
    )
  })
  
  output$Daisy <- renderUI({
    withMathJax(
      p('$$ S_6(t)=S\\Big(1+\\frac{1}{10} cos(t)\\Big)  $$'),
      p('$$ a(T)=wa_w+ba_b+u\\Big(\\frac{e^{\\gamma(T+\\delta)}}{e^{\\gamma(T+\\delta)}+1} (\\beta-\\alpha) + \\alpha\\Big)$$'),
      p("$$ T_w=T+c(a-a_w) $$"),
      p("$$ T_b=T+c(a-a_b) $$"),
      p("$$ F_w=1-k(T_0-T_w)^2 $$"),
      p("$$ F_b=1-k(T_0-T_b)^2 $$"),
      p("$$ w'=w+w(uF_w-D) $$"),
      p("$$ b'=b+b(uF_b-D) $$"),
      p('$$ T =\\frac{ R_{in}(1-a(T))+K \\bar{T}-A }{ B+K }$$')
    )
  })
  
  output$plot_albedo <- renderPlot({ 
    plot_albedo <- ggplot(data.frame(x= seq(-15,15, by=0.1),y=alb(seq(-15,15, by=0.1),input$ai,input$ab,input$gamma,input$delta)))+geom_line(aes(x,y),colour="blue")+xlab("Temperature")+ylab("Albedo")
    plot_albedo  })
  
  output$plot1_Run0 <- renderPlot({ 
    
    data01 <- ebm01(300,input$S,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta)
    
    plot01 <- ggplot(data01,aes(Zones)) +
      geom_line(aes(y=T, colour = "Temperature",linetype="Temperature")) +
      geom_line(aes(y=a*25, colour = "Albedo",linetype="Albedo"))+
      geom_line(aes(y=Ti, colour = "Initial Temperature",linetype="Initial Temperature"))+xlab("Latitude")+ylab("Temperature")+scale_colour_manual(name="Legend",values=c("blue","red","red"))+scale_linetype_manual(name="Legend",values=c("Temperature"=1, "Albedo"=1, "Initial Temperature"=2))
    
    plot01
   })
  output$plot2_Run0 <- renderPlot({ 
    data02 <-ebm02(300,input$S,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta)
    
    plot02 <- ggplot(data02) +
      geom_line(aes(t, Temperature),colour = 'green')+xlab("Time")+ylab("Mean Temperature")
    plot02
  })
  
  output$plot1_Run1 <- renderPlot({
    data11 <- ebm11(100,300,input$S,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta)
   
    plot11 <- ggplot(data11,aes(J,Temp))+geom_line(aes(J, Temp),colour = 'red')+xlab("Time")+ylab("Mean Temperature")+ggtitle("Non Dynamic Mean Temperature (T(t) independent from T(t-1))")
    
    plot11  })
  output$plot2_Run1 <- renderPlotly({
    TEMP1 <- ebm12(100,300,input$S,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta)
    
    plot_ly(z=~TEMP1)%>% add_surface() %>%   layout(
      title = "With Initialization", scene = list(
        xaxis = list(title = "S/100"),
        yaxis = list(title = "Latitude"),
        zaxis = list(title = "T")  )) 
  })
  
  output$plot1_Run2 <- renderPlot({
    data21 <- ebm21(180,300,input$S,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta )
    
    plot21 <- ggplot(data21,aes(J,Temp2))+geom_line(aes(J, Temp2),colour = 'red')+xlab("Time")+ylab("Mean Temperature")+ggtitle("Sinusoidal Perturbation")
    
    plot21 })
  output$plot2_Run2 <- renderPlotly({
    TEMP2 <- ebm22(180,300,input$S,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta )
    
    plot_ly(z=~TEMP2) %>% add_surface() %>% layout(
      title = "Without Initialization",scene = list(
        xaxis = list(title = "S/100"),
        yaxis = list(title = "Latitude"),
        zaxis = list(title = "T") ))
    })
  
  output$plot1_Run3 <- renderPlot({
    
    data31 <- ebm31(200,300,input$S,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta )

    plot31 <- ggplot(data31,aes(J,Temp))+geom_line(aes(J, Temp),colour = 'red')+xlab("Time")+ylab("Mean Temperature")+ggtitle("Gaussian Perturbation")
    
    plot31 })
  output$plot1_Run4 <- renderPlot({
    data41 <- ebm41(200,300,input$S,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta )
    plot41 <- ggplot(data41,aes(J,Temp))+geom_line(aes(J, Temp),colour = 'red')+xlab("Time")+ylab("Mean Temperature")+ggtitle("Delta Perturbation")
    
    plot41  })
  output$plot2_Run3 <- renderPlotly({
    TEMP3 <- ebm32(200,300,input$S,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta )

    plot_ly(z=~TEMP3) %>% add_surface()%>% layout(
      title = "Without Initialization", scene = list(
        xaxis = list(title = "S/100"),
        yaxis = list(title = "Latitude"),
        zaxis = list(title = "T")   ))
    })
  output$plot2_Run4 <- renderPlotly({
    TEMP4 <- ebm42(200,300,input$S,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta )

    plot_ly(z=~TEMP4) %>% add_surface()%>% layout(
      title = "Without Initialization",scene = list(
        xaxis = list(title = "S/100"),
        yaxis = list(title = "Latitude"),
        zaxis = list(title = "T")    ))
  })
  
  output$plot1_Hysteresis <- renderPlot({
    data23 <- ebm23(180,300,input$S,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta )
    
    plot23 <- ggplot(data23,aes(Sarr,Temp2))+geom_point(aes(Sarr, Temp2),colour = 'yellow')+ylab("Mean Temperature")+xlab("S_2(t)")+ggtitle("Mean temperature vs. Sinusoidally Fluctuating Incoming Radiation")
    plot23  })
  output$plot2_Hysteresis <- renderPlot({
    data51 <- ebm51(300,60,input$A,input$B,input$K,input$ai,input$ab,input$gamma,input$delta )
    
    plot51 <- ggplot(data51,aes(Sarr,Temp5,group=1))+geom_point(aes(Sarr, Temp5),colour = 'yellow')+geom_segment(aes(x = Sarr[89], y = Temp5[89], xend = Sarr[90], yend = Temp5[90]),colour = 'yellow')+geom_segment(aes(x = Sarr[258], y = Temp5[258], xend = Sarr[259], yend = Temp5[259]),colour = 'yellow')+ylab("Mean Temperature")+xlab("S_5(t)")+ggtitle("Mean temperature vs. Linearly Fluctuating Incoming Radiation")
    plot51  })
  
  output$plot1_noDaisy <- renderPlot({
    
    data_ND1 <- ebm_ND1(500,input$S,input$w0,input$b0,input$A,input$B,input$K,input$ai,input$ab,input$aW,input$aB,input$gamma,input$delta)
    
    plot_ND1 <- ggplot(data_ND1,aes(Zones))+geom_line(aes(y=w, colour = "% White"))+geom_line(aes(y=b, colour = "% Black"))+geom_line(aes(y=u, colour="% Bare Ground"))+geom_line(aes(y=T/5, colour="Temperature"))+ylab("Temperature & Albedo")+xlab("Latitude")+ggtitle("Without Daisies")+scale_colour_manual(name="Legend",values=c("brown","black","white", "red"))
    plot_ND1  
    })
  output$plot2_noDaisy <- renderPlot({
    data_ND2 <- ebm_ND2(500,input$S,input$w0,input$b0,input$A,input$B,input$K,input$ai,input$ab,input$aW,input$aB,input$gamma,input$delta)

    plot_ND2 <- ggplot(data_ND2,aes(I))+geom_line(aes(y=Barr,colour="% Black"))+ geom_line(aes(y=Warr,colour="% White"))+ geom_line(aes(y=Uarr, colour="% Bare Ground"))+geom_line(aes(y=Tarr/25,colour="Temperature"))+ylab("Temperature & Albedo")+xlab("Time")+ggtitle("Without Daisies")+scale_colour_manual(name="Legend",values=c("brown","black","white", "red"))
    
    plot_ND2 })
  output$plot3_noDaisy <- renderPlotly({

    TEMP <- ebm_ND3(500,input$S,input$w0,input$b0,input$A,input$B,input$K,input$ai,input$ab,input$aW,input$aB,input$gamma,input$delta)
    
    plot_ly(z=~TEMP)%>% add_surface() %>% layout( title="Without Daisies",
                                                  scene = list(
                                                    xaxis = list(title = "Time"),
                                                    yaxis = list(title = "Latitude"),
                                                    zaxis = list(title = "T")  )) 
    })
  
  output$plot1_Daisy <- renderPlot({
    data_D1 <- ebm_D1(500,input$S,input$w0,input$b0,input$T0,input$c,input$k,input$A,input$B,input$K,input$ai,input$ab,input$aW,input$aB,input$gamma,input$delta)

    plot_D1 <- ggplot(data_D1 ,aes(Zones))+geom_line(aes(y=w, colour = "% White"))+geom_line(aes(y=b, colour = "% Black"))+geom_line(aes(y=u, colour="% Bare Ground"))+geom_line(aes(y=T/5, colour="Temperature"))+ylab("Temperature & Albedo")+xlab("Latitude")+ggtitle("With Daisies")+scale_colour_manual(name="Legend",values=c("brown","black","white", "red"))

    plot_D1  }) 
  output$plot2_Daisy <- renderPlot({
    data_D2 <- ebm_D2(500,input$S,input$w0,input$b0,input$T0,input$c,input$k,input$A,input$B,input$K,input$ai,input$ab,input$aW,input$aB,input$gamma,input$delta)

    plot_D2 <- ggplot(data_D2,aes(I))+geom_line(aes(y=Barr,colour="% Black"))+ geom_line(aes(y=Warr,colour="% White"))+ geom_line(aes(y=Uarr, colour="% Bare Ground"))+geom_line(aes(y=Tarr/25,colour="Temperature"))+ylab("Temperature & Albedo")+xlab("Time")+ggtitle("With Daisies")+scale_colour_manual(name="Legend",values=c("brown","black","white", "red"))
    
    plot_D2
    })
  output$plot3_Daisy <- renderPlotly({
    
    TEMP <- ebm_D3(500,input$S,input$w0,input$b0,input$T0,input$c,input$k,input$A,input$B,input$K,input$ai,input$ab,input$aW,input$aB,input$gamma,input$delta)

    plot_ly(z=~TEMP)%>% add_surface() %>% layout(title="With Daisies",
                                                 scene = list(
                                                   xaxis = list(title = "Time"),
                                                   yaxis = list(title = "Latitude"),
                                                   zaxis = list(title = "T")  )) 
    })
}

# APP ----
shinyApp(ui = ui, server = server)







