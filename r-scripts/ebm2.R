
# Aumentiamo la risoluzione spaziale latitudinale (x5)
library(tidyverse)
library(plotly)

Func <-  function(x){
  0.7768699*cos(0.0164348*x)^2+0.4617747
}
gauss <-  function(x,m,sd,b){
  ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
}
Sun1 <- function(x,a){a*x/100}
Sun2 <- function(x){1370*(sinpi((x+90)/180))^2}
Sun3 <- function(x){1370-(((1370)/(sqrt(2*pi*1^2)))*exp(-(x-50)^2/(2*1^2)))} #gauss
Sun4 <- function(x){ifelse(x==50,1370/3,1370)}

Incident <- function(x,y){x*y/4}
#EBM <-function(S,A,B,K,a_ice,a,T_c
#Input Data
alpha <- 1
S <- 1370
A <- 204
B <- 2.17
K <- 3.86
a_ice <- 0.6
a <- 0.3
T_c <- -10
Temperature <- rep(0,300)
Temp <- rep(0,100)
Temp2 <- rep(0,180)
Tempw <- rep(0,200)
J <- rep(0,100)
Albedo <- rep(0,300)
t <- c(1:300)

# 
Zones <- seq(-89, 89, by = 2)
cosZones <- abs(cospi(Zones/180))
SunWt <- Func(Zones)
Rin <- Incident(S,SunWt)
Ti <- gauss(Zones,0,50,31.6)
T <- gauss(Zones,0,50,31.6)
a_i <- ifelse(T<T_c, a_ice, a)

#Test Run

for(i in c(1:300)) {Tcos <- cosZones*T
Tm <- sum(Tcos)/sum(cosZones)
print(Tm)
T <- (Rin*(1-a_i)+K*Tm-A) / (B+K)
a_i <- ifelse(T<T_c, a_ice, a)
Temperature[i] <-  Tm
Albedo[i] <- a_i 
} 
M <- matrix(NA, nrow=90, ncol=100)
data1 <- data.frame(Zones,T,a_i,Ti)
plot1 <- ggplot(data1) +
  geom_line(aes(Zones, T),colour = 'red') +
  geom_line(aes(Zones, a_i*25), colour = 'blue')+
  geom_line(aes(Zones, Ti), colour = 'red', linetype = "dashed")+
  scale_y_continuous(sec.axis = sec_axis(~./25,name = "a"))
plot1

data2 <- data.frame(t,Temperature)
plot2 <- ggplot(data2) +
  geom_line(aes(t, Temperature),colour = 'green')+ylab("<T>")
plot2

#Run1
for(j in c(1:100)){
  T <- gauss(Zones,0,50,31.6)
  a_i <- ifelse(T<T_c, a_ice, a)
  S <- Sun1(1370,j)
  Rin <- Incident(S,SunWt)
  
  for(i in c(1:300))
  {Tcos <- cosZones*T
  Tm <- sum(Tcos)/sum(cosZones)
  T <- (Rin*(1-a_i)+K*Tm-A) / (B+K)
  a_i <- ifelse(T<T_c, a_ice, a)
  #Temperature[i] <-  Tm
  #Albedo[i] <- a_i 
  }
  M[,j] <- T
  J[j] <- j
  Temp[j] <- Tm
}

data3 <- data.frame(J,Temp)
plot3 <- ggplot(data3,aes(J,Temp))+geom_line(aes(J, Temp),colour = 'red')+ylab("<T>")
plot3

#Run2: senza inizializzazione di (T,a)

N <- matrix(NA, nrow=90, ncol=180)
T2 <- gauss(Zones,0,50,31.6)
a_i <- ifelse(T2<T_c, a_ice, a)
Sarr <- rep(0,180) 
for(j in c(1:180)){
  S <- Sun2(j)
  Rin <- Incident(S,SunWt)
  for(i in c(1:300))
  {Tcos <- cosZones*T2
  Tm <- sum(Tcos)/sum(cosZones)
  #print(Tm)
  T2 <- (Rin*(1-a_i)+K*Tm-A) / (B+K)
  a_i <- ifelse(T2<T_c, a_ice, a)
  #Temperature[i] <-  Tm
  #Albedo[i] <- a_i 
  }
  Sarr[j] <- Sun2(j)
  N[,j] <- T2
  J[j] <- j
  Temp2[j] <- Tm
  #print(Tm)
}

data3 <- data.frame(J,Temp2)
plot3 <- ggplot(data3,aes(J,Temp2))+geom_line(aes(J, Temp2),colour = 'red')+ylab("<T>")
plot3

data <- data.frame(Sarr,Temp2)
plot <- ggplot(data,aes(Sarr,Temp2))+geom_point(aes(Sarr, Temp2),colour = 'yellow')+ylab("<T>")
plot

#Run3: senza inizializzazione di (T,a)

H <- matrix(NA, nrow=90, ncol=200)
T3 <- gauss(Zones,0,50,31.6)
a_i <- ifelse(T2<T_c, a_ice, a)

for(j in c(1:200)){
  S <- Sun3(j)
  Rin <- Incident(S,SunWt)
  for(i in c(1:300))
  {Tcos <- cosZones*T3
  Tm <- sum(Tcos)/sum(cosZones)
  #print(Tm)
  T3 <- (Rin*(1-a_i)+K*Tm-A) / (B+K)
  a_i <- ifelse(T3<T_c, a_ice, a)
  #Temperature[i] <-  Tm
  #Albedo[i] <- a_i 
  }
  H[,j] <- T3
  J[j] <- j
  Temp[j] <- Tm
  #print(Tm)
}

data3 <- data.frame(J,Temp)
plot3 <- ggplot(data3,aes(J,Temp))+geom_line(aes(J, Temp),colour = 'red')+ylab("<T>")
plot3

#Run4: senza inizializzazione di (T,a)

U <- matrix(NA, nrow=90, ncol=200)
T4 <- gauss(Zones,0,50,31.6)
a_i <- ifelse(T2<T_c, a_ice, a)

for(j in c(1:200)){
  S <- Sun4(j)
  Rin <- Incident(S,SunWt)
  for(i in c(1:3))
  {Tcos <- cosZones*T4
  Tm <- sum(Tcos)/sum(cosZones)
  #print(Tm)
  T4 <- (Rin*(1-a_i)+K*Tm-A) / (B+K)
  a_i <- ifelse(T4<T_c, a_ice, a)
  #Temperature[i] <-  Tm
  #Albedo[i] <- a_i 
  }
  U[,j] <- T4
  J[j] <- j
  Temp[j] <- Tm
  #print(Tm)
}

data3 <- data.frame(J,Temp)
plot3 <- ggplot(data3,aes(J,Temp))+geom_line(aes(J, Temp),colour = 'red')+ylab("<T>")
plot3



plot_ly(z=~M) %>% add_surface()
plot_ly(z=~N) %>% add_surface()
plot_ly(z=~H) %>% add_surface()
plot_ly(z=~U) %>% add_surface()


# CICLO ISTERESI
library(tidyverse)
library(plotly)

#functions
Func <-  function(x){
  0.7768699*cos(0.0164348*x)^2+0.4617747
}
gauss <-  function(x,m,sd,b){
  ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
}
Sun1 <- function(x,a){a*x/100}
Sun2 <- function(x){1370*(sinpi((x+90)/180))^2}
Sun3 <- function(x){1370-(((1370)/(sqrt(2*pi*1^2)))*exp(-(x-50)^2/(2*1^2)))} #gauss
Sun4 <- function(x){ifelse(x==50,1370/3,1370)}
Sunw <- function(x){(1/100) * abs(x-101)}
Incident <- function(x,y){x*y/4}

W <- matrix(NA, nrow=90, ncol=200)
Tw <- gauss(Zones,0,50,31.6)
a_i <- ifelse(Tw<T_c, a_ice, a)
Sarr <- rep(0,200) 
for(j in c(1:200)){
  S <- Sunw(j)
  Rin <- 1370*Incident(S,SunWt)
  for(i in c(1:300))
  {Tcos <- cosZones*Tw
  Tm <- sum(Tcos)/sum(cosZones)
  #print(Tm)
  Tw <- (Rin*(1-a_i)+K*Tm-A) / (B+K)
  a_i <- ifelse(Tw<T_c, a_ice, a)
  #Temperature[i] <-  Tm
  #Albedo[i] <- a_i 
  }
  Sarr[j] <- Sunw(j)
  W[,j] <- Tw
  J[j] <- j
  Tempw[j] <- Tm
  #print(Tm)
}


dataw <- data.frame(Sarr,Tempw)
plotw <- ggplot(dataw,aes(Sarr,Tempw))+geom_point(aes(Sarr, Tempw),colour = 'yellow')+ylab("<T>")
plotw






