# Simple Energy Balance Model [Climate Modelling Primer]
library(tidyverse)

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
Albedo <- rep(0,300)

# 
Zones <- c(85,75,65,55,45,35,25,15,6)
cosZones <- cospi(Zones/180)
SunWt <- c(0.5,0.531,0.624,0.77,0.892,1.021,1.12,1.189,1.219)
data3 <- data.frame(Zones,SunWt)
Func <-  function(Zones,d,n,k){
         d*cos(n*Zones)+k
}
plot(Zones,SunWt)
curve(Func(x,1,0.015,0.3),add=TRUE)
Fit <-  nls(SunWt ~ Func(Zones,d,n,k),data= data3, start=list(d=1,n=0.015,k=0.3))

plot(Zones,SunWt,main="ciao")
curve(Func(x,0.7768699,0.0164348,0.4617747))
ggplot(data3, aes(Zones,SunWt))+geom_point()+geom_smooth(aes(y=Func(Zones,0.7768699,0.0164348,0.4617747)))

Rin <- S*SunWt/4
T <- c(-15,-15,-5,5,10,15,18,22,24)
a_i <- ifelse(T<T_c, a_ice, a)

for(i in c(1:300))
  {Tcos <- cosZones*T
   Tm <- sum(Tcos)/sum(cosZones)
   #print(Tm)
   T <- (Rin*(1-a_i)+K*Tm-A) / (B+K)
   a_i <- ifelse(T<T_c, a_ice, a)
   Temperature[i] <-  Tm
   Albedo[i] <- a_i 
   }


data1 <- data.frame(Zones,T,a_i)
plot1 <- ggplot(data1) +
             geom_line(aes(Zones, T),colour = 'red') +
             geom_line(aes(Zones, a_i*25), colour = 'blue')+
             scale_y_continuous(sec.axis = sec_axis(~./25,name = "a"))
plot1
t <- c(1:300)
data2 <- data.frame(t,Temperature)
plot2 <- ggplot(data2) +
          geom_line(aes(t, Temperature),colour = 'green')
plot2

# Aumentare length(Zones) ?
# Complessificare? 
# Plottare Tm=Tm(S) variando S in (0,1)
# Glaciated Mode
# Daisy : margerite distribuite per Zone + Growth Rate 







