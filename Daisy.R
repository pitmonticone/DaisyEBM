#Packages
library(tidyverse)
library(plotly)



#Sun <- function(x){1370*(1+0.1*cospi(x/180))}




#Inpit Parameters
S <- 1370
A <- 204
B <- 2.17
K <- 3.86
ai <- 0.62
ab <- 0.25
aW <- 0.75
aB <- 0.25
#aU <- 0.4
c <- 7
k <- 0.003265*0.75
T0 <- 20
D <- 0.3 # Death Rate

Zones <- seq(-89, 89, by = 2)

#Functions
Func <-  function(x){
  0.7768699*cos(0.0164348*x)^2+0.4617747
}
gauss <-  function(x,m,sd,b){
  ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
}
Incident <- function(x,y){x*y/4}
alb <- function(x,a,b){ (-exp(1.6*x+10)/(exp(1.6*x+10)+1))*(a-b)+a} #naked

#Initialization
cosZones <- abs(cospi(Zones/180))
SunWt <- Func(Zones)
Rin <- Incident(S,SunWt)
T <- gauss(Zones,0,50,31.6)-6

w <- rep(0.5,length(Zones))
b <- rep(0.2,length(Zones))
u <- rep(0.3,length(Zones))
a <- w*aW+b*aB+u*alb(T,0.62,0.25)
Barr <- rep(0,500)
Warr <- rep(0,500)
Uarr <- rep(0,500)
Tarr <- rep(0,500)
I <- rep(0,500)

TEMP <- matrix(NA, nrow=90, ncol=500)
for(i in c(1:500)) {
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
a <- w*aW+b*aB+u*alb(T,0.62,0.25)
Barr[i] <- b[45]
Warr[i] <- w[45]
Uarr[i] <- u[45]
I[i] <- i
Tarr[i] <- T[45]
} 

plot <- ggplot(data.frame(Zones,w,b,u,T))+geom_line(aes(Zones,w,color="white"))+geom_line(aes(Zones,b,color="Black"),color="black")+geom_line(aes(Zones,u),color="brown")+geom_line(aes(Zones,T/5),color="green")
plot
plot2 <- ggplot(data.frame(I,Barr,Warr,Uarr,Tarr))+geom_line(aes(I,Barr),color="black")+ geom_line(aes(I,Warr),color="white")+ geom_line(aes(I,Uarr),color="brown")+geom_line(aes(I,Tarr/25),color="pink")
plot2

plot_ly(z=~TEMP)%>% add_surface() %>% layout( 
    scene = list(
    xaxis = list(title = "Time"),
    yaxis = list(title = "Latitude"),
    zaxis = list(title = "T")  )) 

