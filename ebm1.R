# Simple Energy Balance Model [Climate Modelling Primer]
library(tidyverse)
#get.params <- function(plot, layer = 1){
  
  # return empty string if the specified geom layer doesn't use stat = "smooth"
  if(!"StatSmooth" %in% class(plot$layers[[layer]]$stat)){
    message("No smoothing function was used in this geom layer.")
    return("")
  }
  
  # recreate data used by this layer, in the format expected by StatSmooth
  # (this code chunk takes heavy reference from ggplot2:::ggplot_build.ggplot)
  layer.data <- plot$layers[[layer]]$layer_data(plot$data)
  layout <- ggplot2:::create_layout(plot$facet, plot$coordinates)
  data <- layout$setup(list(layer.data), plot$data, plot$plot_env)
  data[[1]] <- plot$layers[[layer]]$compute_aesthetics(data[[1]], plot)
  scales <- plot$scales
  data[[1]] <- ggplot2:::scales_transform_df(scales = scales, df = data[[1]])
  layout$train_position(data, scales$get_scales("x"), scales$get_scales("y"))
  data <- layout$map_position(data)[[1]]
  
  # set up stat params (e.g. replace "auto" with actual method / formula)
  stat.params <- suppressMessages(
    plot$layers[[layer]]$stat$setup_params(data = data, 
                                           params = plot$layers[[layer]]$stat_params)
  )
  
  # reverse the last step in setup_params; we don't need the actual function
  # for mgcv::gam, just the name
  if(identical(stat.params$method, mgcv::gam)) stat.params$method <- "gam"
  
  print(stat.params)
  
  return(paste0("Method: ", stat.params$method, ", Formula: ", deparse(stat.params$formula)))
}


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
         d*cos(n*Zones)^2+k
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
   print(Tm)
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







