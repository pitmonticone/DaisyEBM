library(shiny)
library(tidyverse)
library(plotly)
library(broom)
library(scales)
library(shinythemes) #---> themeSelector()


# INITIALIZATION ----

Func <-  function(x){
  0.7768699*cos(0.0164348*x)^2+0.4617747
}
gauss <-  function(x,m,sd,b){
  ((24+b)/(0.00798*sqrt(2*pi*sd^2)))*exp(-(x-m)^2/(2*sd^2))-b
}
Sun1 <- function(x,a){a*x/100}
Sun2 <- function(x){1370*(sinpi((x+90)/180))^2}
Sun3 <- function(x){1370-(((1370)/(sqrt(2*pi*1^2)))*exp(-(x-50)^2/(2*1^2)))}
Sun4 <- function(x){ifelse(x==50,1370/3,1370)}
Sun5 <- function(x){(1/100) * (abs(x-150)+25)}
Sun6 <- function(x){1370*(1+0.1*cospi(x/180))}
Incident <- function(x,y){x*y/4}
Step <- function(x,c){ifelse(x<c, 0.6, 0.3)}
alb <- function(x,a,b){
  (-exp(2.2*x+10)/(exp(2.2*x+10)+1))*(a-b)+a}
Temperature <- rep(0,300)
Temp <- rep(0,100)
Temp2 <- rep(0,180)
Temp5 <- rep(0,300)
J <- rep(0,100)
Albedo <- rep(0,300)
t <- c(1:300)
Zones <- seq(-89, 89, by = 2)
Ti <- gauss(Zones,0,50,31.6)
cosZones <- abs(cospi(Zones/180))
SunWt <- Func(Zones)
Barr <- rep(0,500)
Warr <- rep(0,500)
Uarr <- rep(0,500)
Tarr <- rep(0,500)
w <- rep(0.5,length(Zones))
b <- rep(0.2,length(Zones))
u <- rep(0.3,length(Zones))
I <- rep(0,500)
TEMP <- matrix(NA, nrow=90, ncol=500)
TEMP1 <- matrix(NA, nrow=90, ncol=100)
T2 <- gauss(Zones,0,50,31.6)
TEMP2 <- matrix(NA, nrow=90, ncol=180)
Sarr <- rep(0,180) 
TEMP4 <- matrix(NA, nrow=90, ncol=200)
T4 <- gauss(Zones,0,50,31.6)
TEMP3 <- matrix(NA, nrow=90, ncol=200)
T3 <- gauss(Zones,0,50,31.6)



# UI ----
ui <- fluidPage(
  theme = shinytheme("spacelab"),
  titlePanel("Climate Modelling"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Input Parameters"),
      numericInput("S",
                   label = "Solar Constant", 
                   value = 1370,
                   min = NA, max = NA, step = NA),
      numericInput("A",
                   label = "A", 
                   value = 204,
                   min = NA, max = NA, step = NA),
      numericInput("B",
                   label = "B", 
                   value = 2.17,
                   min = NA, max = NA, step = NA),
      numericInput("K",
                   label = "K", 
                   value = 3.86,
                   min = NA, max = NA, step = NA),
      numericInput("ai",
                   label = "Ice Albedo", 
                   value = 0.62,
                   min = NA, max = NA, step = NA),
      numericInput("ab",
                   label = "Bare Ground Albedo", 
                   value = 0.25,
                   min = NA, max = NA, step = NA),
      numericInput("aW",
                   label = "White Albedo", 
                   value = 0.75,
                   min = NA, max = NA, step = NA),
      numericInput("aB",
                   label = "Black Albedo", 
                   value = 0.25,
                   min = NA, max = NA, step = NA),
      numericInput("c",
                   label = "c", 
                   value = 7,
                   min = NA, max = NA, step = NA),
      numericInput("k",
                   label = "k", 
                   value = 0.003265*0.75,
                   min = NA, max = NA, step = NA),
      numericInput("T0",
                   label = "T_0", 
                   value = 20,
                   min = NA, max = NA, step = NA),
      numericInput("D",
                   label = "Death Rate", 
                   value = 0.3,
                   min = NA, max = NA, step = NA)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Albedo", plotOutput("plot_albedo")),
        tabPanel("Run 0", plotOutput("plot1_Run0"),plotOutput("plot2_Run0")),
        tabPanel("Run 1", plotOutput("plot1_Run1"),plotlyOutput("plot2_Run1")),
        tabPanel("Run 2", plotOutput("plot1_Run2"),plotlyOutput("plot2_Run2")),
        tabPanel("Run 3 & 4", plotOutput("plot1_Run3"),plotOutput("plot1_Run4"),plotlyOutput("plot2_Run3"),plotlyOutput("plot2_Run4")),
        tabPanel("Hysteresis Cycles", plotOutput("plot1_Hysteresis"), plotOutput("plot2_Hysteresis")),
        tabPanel("Without Daisies", plotOutput("plot1_noDaisy"), plotOutput("plot2_noDaisy"), plotlyOutput("plot3_noDaisy")),
        tabPanel("With Daisies", plotOutput("plot1_Daisy"), plotOutput("plot2_Daisy"), plotlyOutput("plot3_Daisy") )
      )
    )
  )
)

# SERVER ----
server <- function(input, output,session) {
  
  output$plot_albedo <- renderPlot({ 
    ggplot(data.frame(x= seq(-15,15, by=0.1),y=alb(seq(-15,15, by=0.1),input$ai,input$ab)))+geom_line(aes(x,y),colour="blue")+xlab("Temperature")+ylab("Albedo")
  })
  
  output$plot1_Run0 <- renderPlot({ 
    
    Rin <- Incident(input$S,SunWt)
    T <- gauss(Zones,0,50,31.6)
    a <- alb(T,input$ai,input$ab)
    
    for(i in c(1:300)) {Tcos <- cosZones*T
    Tm <- sum(Tcos)/sum(cosZones)
    T <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
    a <- alb(T,input$ai,input$ab)
    Temperature[i] <-  Tm
    } 
    
    data1 <- data.frame(Zones,T,a,Ti)
    ggplot(data1,aes(Zones)) +
      geom_line(aes(y=T, colour = "Temperature",linetype="Temperature")) +
      geom_line(aes(y=a*25, colour = "Albedo",linetype="Albedo"))+
      geom_line(aes(y=Ti, colour = "Initial Temperature",linetype="Initial Temperature"))+xlab("Latitude")+ylab("Temperature")+scale_colour_manual(name="Legend",values=c("blue","red","red"))+scale_linetype_manual(name="Legend",values=c("Temperature"=1, "Albedo"=1, "Initial Temperature"=2))
  })
  
  output$plot2_Run0 <- renderPlot({ 
    Rin <- Incident(input$S,SunWt)
    T <- gauss(Zones,0,50,31.6)
    a <- alb(T,input$ai,input$ab)
    
    for(i in c(1:300)) {Tcos <- cosZones*T
    Tm <- sum(Tcos)/sum(cosZones)
    T <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
    a <- alb(T,input$ai,input$ab)
    Temperature[i] <-  Tm
    } 
    data2 <- data.frame(t,Temperature)
    ggplot(data2) +
      geom_line(aes(t, Temperature),colour = 'green')+xlab("Time")+ylab("Mean Temperature")
  })
  
  output$plot1_Run1 <- renderPlot({
    for(j in c(1:100)){
      T <- gauss(Zones,0,50,31.6)
      a <- alb(T,input$ai,input$ab)
      S <- Sun1(1370,j)
      Rin <- Incident(S,SunWt)
      
      for(i in c(1:300))
      {Tcos <- cosZones*T
      Tm <- sum(Tcos)/sum(cosZones)
      T <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      a <- alb(T,input$ai,input$ab)}
      TEMP1[,j] <- T
      J[j] <- j
      Temp[j] <- Tm
    }
    
    data3 <- data.frame(J,Temp)
    ggplot(data3,aes(J,Temp))+geom_line(aes(J, Temp),colour = 'red')+xlab("Time")+ylab("Mean Temperature")+ggtitle("Non Dynamic Mean Temperature (T(t) independent from T(t-1))")
    
  })
  
  output$plot2_Run1 <- renderPlotly({
    for(j in c(1:100)){
      T <- gauss(Zones,0,50,31.6)
      a <- alb(T,input$ai,input$ab)
      S <- Sun1(1370,j)
      Rin <- Incident(S,SunWt)
      
      for(i in c(1:300))
      {Tcos <- cosZones*T
      Tm <- sum(Tcos)/sum(cosZones)
      T <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      a <- alb(T,input$ai,input$ab)}
      TEMP1[,j] <- T
      J[j] <- j
      Temp[j] <- Tm
    }
    
    plot_ly(z=~TEMP1)%>% add_surface() %>%   layout(
      title = "With Initialization", scene = list(
        xaxis = list(title = "S/100"),
        yaxis = list(title = "Latitude"),
        zaxis = list(title = "T")  )) 
  })
  
  output$plot1_Run2 <- renderPlot({
    a <- alb(T2,input$ai,input$ab)

    for(j in c(1:180)){
      S <- Sun2(j)
      Rin <- Incident(S,SunWt)
      for(i in c(1:300))
      {Tcos <- cosZones*T2
      Tm <- sum(Tcos)/sum(cosZones)
      T2 <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      a <- alb(T2,input$ai,input$ab)
      }
      Sarr[j] <- Sun2(j)
      TEMP2[,j] <- T2
      J[j] <- j
      Temp2[j] <- Tm
    }
    
    data3 <- data.frame(J,Temp2)
    ggplot(data3,aes(J,Temp2))+geom_line(aes(J, Temp2),colour = 'red')+xlab("Time")+ylab("Mean Temperature")+ggtitle("Sinusoidal Perturbation")
  })
  
  output$plot2_Run2 <- renderPlotly({
    a <- alb(T2,input$ai,input$ab)
    
    for(j in c(1:180)){
      S <- Sun2(j)
      Rin <- Incident(S,SunWt)
      for(i in c(1:300))
      {Tcos <- cosZones*T2
      Tm <- sum(Tcos)/sum(cosZones)
      T2 <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      a <- alb(T2,input$ai,input$ab)
      }
      Sarr[j] <- Sun2(j)
      TEMP2[,j] <- T2
      J[j] <- j
      Temp2[j] <- Tm
    }
    
    plot_ly(z=~TEMP2) %>% add_surface() %>% layout(
      title = "Without Initialization",scene = list(
        xaxis = list(title = "S/100"),
        yaxis = list(title = "Latitude"),
        zaxis = list(title = "T") ))
  })
  
  output$plot1_Run3 <- renderPlot({
    
    TEMP3 <- matrix(NA, nrow=90, ncol=200)
    T3 <- gauss(Zones,0,50,31.6)
    a <- alb(T3,input$ai,input$ab)
    
    for(j in c(1:200)){
      S <- Sun3(j)
      Rin <- Incident(S,SunWt)
      for(i in c(1:300))
      {Tcos <- cosZones*T3
      Tm <- sum(Tcos)/sum(cosZones)
      T3 <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      a <- alb(T3,input$ai,input$ab)
      }
      TEMP3[,j] <- T3
      J[j] <- j
      Temp[j] <- Tm
    }
    
    data3 <- data.frame(J,Temp)
    ggplot(data3,aes(J,Temp))+geom_line(aes(J, Temp),colour = 'red')+xlab("Time")+ylab("Mean Temperature")+ggtitle("Gaussian Perturbation")
  })
  
  output$plot1_Run4 <- renderPlot({
    a <- alb(T4,input$ai,input$ab)
    
    for(j in c(1:200)){
      S <- Sun4(j)
      Rin <- Incident(S,SunWt)
      for(i in c(1:300))
      {Tcos <- cosZones*T4
      Tm <- sum(Tcos)/sum(cosZones)
      T4 <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      a <- alb(T4,input$ai,input$ab)
      }
      TEMP4[,j] <- T4
      J[j] <- j
      Temp[j] <- Tm
    }
    
    data3 <- data.frame(J,Temp)
    ggplot(data3,aes(J,Temp))+geom_line(aes(J, Temp),colour = 'red')+xlab("Time")+ylab("Mean Temperature")+ggtitle("Delta Perturbation")
  })
  
  output$plot2_Run3 <- renderPlotly({
    a <- alb(T3,input$ai,input$ab)
    
    for(j in c(1:200)){
      S <- Sun3(j)
      Rin <- Incident(S,SunWt)
      for(i in c(1:300))
      {Tcos <- cosZones*T3
      Tm <- sum(Tcos)/sum(cosZones)
      T3 <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      a <- alb(T3,input$ai,input$ab)
      }
      TEMP3[,j] <- T3
      J[j] <- j
      Temp[j] <- Tm
    }
    
    plot_ly(z=~TEMP3) %>% add_surface()%>% layout(
      title = "Without Initialization", scene = list(
        xaxis = list(title = "S/100"),
        yaxis = list(title = "Latitude"),
        zaxis = list(title = "T")   ))
  })
  
  output$plot2_Run4 <- renderPlotly({
    a <- alb(T4,input$ai,input$ab)
    
    for(j in c(1:200)){
      S <- Sun4(j)
      Rin <- Incident(S,SunWt)
      for(i in c(1:300))
      {Tcos <- cosZones*T4
      Tm <- sum(Tcos)/sum(cosZones)
      T4 <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      a <- alb(T4,input$ai,input$ab)
      }
      TEMP4[,j] <- T4
      J[j] <- j
      Temp[j] <- Tm
    }
    plot_ly(z=~TEMP4) %>% add_surface()%>% layout(
      title = "Without Initialization",scene = list(
        xaxis = list(title = "S/100"),
        yaxis = list(title = "Latitude"),
        zaxis = list(title = "T")    ))
  })
  
  output$plot1_Hysteresis <- renderPlot({
    T2 <- gauss(Zones,0,50,31.6)
    TEMP2 <- matrix(NA, nrow=90, ncol=180)
    a <- alb(T2,input$ai,input$ab)
    Sarr <- rep(0,180) 
    for(j in c(1:180)){
      S <- Sun2(j)
      Rin <- Incident(S,SunWt)
      for(i in c(1:300))
      {Tcos <- cosZones*T2
      Tm <- sum(Tcos)/sum(cosZones)
      T2 <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      a <- alb(T2,input$ai,input$ab)
      }
      Sarr[j] <- Sun2(j)
      TEMP2[,j] <- T2
      J[j] <- j
      Temp2[j] <- Tm
    }
    
    data <- data.frame(Sarr,Temp2)
    ggplot(data,aes(Sarr,Temp2))+geom_point(aes(Sarr, Temp2),colour = 'yellow')+ylab("Mean Temperature")+xlab("S_2(t)")+ggtitle("Mean temperature vs. Sinusoidally Fluctuating Incoming Radiation")
  })
  
  output$plot2_Hysteresis <- renderPlot({
    
    Temp5 <- rep(0,300)
    T5 <- gauss(Zones,0,50,31.6)
    a <- alb(T5,input$ai,input$ab)
    Sarr <- rep(0,300) 
    for(j in c(0:300)){
      S <- Sun5(j)
      Rin <- 1370*Incident(S,SunWt)
      for(i in c(1:60))
      {Tcos <- cosZones*T5
      Tm <- sum(Tcos)/sum(cosZones)
      T5 <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      a <- alb(T5,input$ai,input$ab)
      }
      Sarr[j] <- Sun5(j)
      J[j] <- j
      Temp5[j] <- Tm
    }
    
    dataw <- data.frame(Sarr,Temp5)
    ggplot(dataw,aes(Sarr,Temp5,group=1))+geom_point(aes(Sarr, Temp5),colour = 'yellow')+geom_segment(aes(x = Sarr[89], y = Temp5[89], xend = Sarr[90], yend = Temp5[90]),colour = 'yellow')+geom_segment(aes(x = Sarr[258], y = Temp5[258], xend = Sarr[259], yend = Temp5[259]),colour = 'yellow')+ylab("Mean Temperature")+xlab("S_5(t)")+ggtitle("Mean temperature vs. Linearly Fluctuating Incoming Radiation")
  })
  
  output$plot1_noDaisy <- renderPlot({
    
    SunWt <- Func(Zones)
    Rin <- Incident(input$S,SunWt)
    T <- gauss(Zones,0,50,31.6)-6
    
    SunWt <- Func(Zones)
    Rin <- Incident(input$S,SunWt)
    T <- gauss(Zones,0,50,31.6)-6
    a <- w*input$aW+b*input$aB+u*alb(T,input$ai,input$ab)
    
    for(i in c(1:500)) {
      S <- Sun6(i)  # oppure costante S <- 1370
      Rin <- Incident(S,SunWt)
      Tcos <- cosZones*T
      Tm <- sum(Tcos)/sum(cosZones)
      T <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      TEMP[,i] <- T
      a <- alb(T,input$ai,input$ab)
      I[i] <- i
      Tarr[i] <- T[45]
    } 
    
    ggplot(data.frame(Zones,w,b,u,T),aes(Zones))+geom_line(aes(y=w, colour = "% White"))+geom_line(aes(y=b, colour = "% Black"))+geom_line(aes(y=u, colour="% Bare Ground"))+geom_line(aes(y=T/5, colour="Temperature"))+ylab("Temperature & Albedo")+xlab("Latitude")+ggtitle("Without Daisies")+scale_colour_manual(name="Legend",values=c("brown","black","white", "red"))
  })
  
  output$plot2_noDaisy <- renderPlot({
    SunWt <- Func(Zones)
    Rin <- Incident(input$S,SunWt)
    T <- gauss(Zones,0,50,31.6)-6
    
    SunWt <- Func(Zones)
    Rin <- Incident(input$S,SunWt)
    T <- gauss(Zones,0,50,31.6)-6
    a <- w*input$aW+b*input$aB+u*alb(T,input$ai,input$ab)
    
    for(i in c(1:500)) {
      S <- Sun6(i)  # oppure costante S <- 1370
      Rin <- Incident(S,SunWt)
      Tcos <- cosZones*T
      Tm <- sum(Tcos)/sum(cosZones)
      T <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      TEMP[,i] <- T
      a <- alb(T,input$ai,input$ab)
      I[i] <- i
      Tarr[i] <- T[45]
    } 
    
    ggplot(data.frame(I,Barr,Warr,Uarr,Tarr),aes(I))+geom_line(aes(y=Barr,colour="% Black"))+ geom_line(aes(y=Warr,colour="% White"))+ geom_line(aes(y=Uarr, colour="% Bare Ground"))+geom_line(aes(y=Tarr/25,colour="Temperature"))+ylab("Temperature & Albedo")+xlab("Time")+ggtitle("Without Daisies")+scale_colour_manual(name="Legend",values=c("brown","black","white", "red"))
  })
  
  output$plot3_noDaisy <- renderPlotly({
    SunWt <- Func(Zones)
    Rin <- Incident(input$S,SunWt)
    T <- gauss(Zones,0,50,31.6)-6
    
    a <- w*input$aW+b*input$aB+u*alb(T,input$ai,input$ab)
    
    for(i in c(1:500)) {
      S <- Sun6(i)  # oppure costante S <- 1370
      Rin <- Incident(S,SunWt)
      Tcos <- cosZones*T
      Tm <- sum(Tcos)/sum(cosZones)
      T <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      TEMP[,i] <- T
      a <- alb(T,input$ai,input$ab)
      I[i] <- i
      Tarr[i] <- T[45]
    } 
    
    plot_ly(z=~TEMP)%>% add_surface() %>% layout( title="Without Daisies",
                                                  scene = list(
                                                    xaxis = list(title = "Time"),
                                                    yaxis = list(title = "Latitude"),
                                                    zaxis = list(title = "T")  )) 
  })
  
  output$plot1_Daisy <- renderPlot({
    
    SunWt <- Func(Zones)
    Rin <- Incident(input$S,SunWt)
    T <- gauss(Zones,0,50,31.6)-6
    a <- w*input$aW+b*input$aB+u*alb(T,input$ai,input$ab)
    
    for(i in c(1:500)) {
      S <- Sun6(i)  # oppure costante S <- 1370
      Rin <- Incident(S,SunWt)
      Tcos <- cosZones*T
      Tm <- sum(Tcos)/sum(cosZones)
      T <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      TEMP[,i] <- T
      Tw <- T+input$c*(a-input$aW)
      Tb <- T+input$c*(a-input$aB)
      Fw <- 1-input$k*(input$T0-Tw)^2
      Fb <- 1-input$k*(input$T0-Tb)^2
      for(j in c(1:length(Zones))){
        if(Fw[j]<0){Fw[j]=0}
        if(Fb[j]<0){Fb[j]=0}  }
      w <- w+w*(u*Fw-input$D)
      b <- b+b*(u*Fb-input$D)
      for(j in c(1:length(Zones))){
        if(w[j]<0.001){w[j]=0.001}
        if(b[j]<0.001){b[j]=0.001}  }
      u <- 1-w-b
      a <- w*input$aW+b*input$aB+u*alb(T,input$ai,input$ab)
      Barr[i] <- b[45]
      Warr[i] <- w[45]
      Uarr[i] <- u[45]
      I[i] <- i
      Tarr[i] <- T[45]
    } 
    
    ggplot(data.frame(Zones,w,b,u,T),aes(Zones))+geom_line(aes(y=w, colour = "% White"))+geom_line(aes(y=b, colour = "% Black"))+geom_line(aes(y=u, colour="% Bare Ground"))+geom_line(aes(y=T/5, colour="Temperature"))+ylab("Temperature & Albedo")+xlab("Latitude")+ggtitle("With Daisies")+scale_colour_manual(name="Legend",values=c("brown","black","white", "red"))
  })
  
  output$plot2_Daisy <- renderPlot({
    
    a <- w*input$aW+b*input$aB+u*alb(T,0.62,0.25)
    
    for(i in c(1:500)) {
      S <- Sun6(i)  # oppure costante S <- 1370
      Rin <- Incident(S,SunWt)
      Tcos <- cosZones*T
      Tm <- sum(Tcos)/sum(cosZones)
      T <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      TEMP[,i] <- T
      Tw <- T+input$c*(a-input$aW)
      Tb <- T+input$c*(a-input$aB)
      Fw <- 1-input$k*(input$T0-Tw)^2
      Fb <- 1-input$k*(input$T0-Tb)^2
      for(j in c(1:length(Zones))){
        if(Fw[j]<0){Fw[j]=0}
        if(Fb[j]<0){Fb[j]=0}  }
      w <- w+w*(u*Fw-input$D)
      b <- b+b*(u*Fb-input$D)
      for(j in c(1:length(Zones))){
        if(w[j]<0.001){w[j]=0.001}
        if(b[j]<0.001){b[j]=0.001}  }
      u <- 1-w-b
      a <- w*input$aW+b*input$aB+u*alb(T,input$ai,input$ab)
      Barr[i] <- b[45]
      Warr[i] <- w[45]
      Uarr[i] <- u[45]
      I[i] <- i
      Tarr[i] <- T[45]
    } 
    
    ggplot(data.frame(I,Barr,Warr,Uarr,Tarr),aes(I))+geom_line(aes(y=Barr,colour="% Black"))+ geom_line(aes(y=Warr,colour="% White"))+ geom_line(aes(y=Uarr, colour="% Bare Ground"))+geom_line(aes(y=Tarr/25,colour="Temperature"))+ylab("Temperature & Albedo")+xlab("Time")+ggtitle("With Daisies")+scale_colour_manual(name="Legend",values=c("brown","black","white", "red"))
  })
  
  output$plot3_Daisy <- renderPlotly({
    
    a <- w*input$aW+b*input$aB+u*alb(T,input$ai,input$ab)

    for(i in c(1:500)) {
      S <- Sun6(i)  # oppure costante S <- 1370
      Rin <- Incident(S,SunWt)
      Tcos <- cosZones*T
      Tm <- sum(Tcos)/sum(cosZones)
      T <- (Rin*(1-a)+input$K*Tm-input$A) / (input$B+input$K)
      TEMP[,i] <- T
      Tw <- T+input$c*(a-input$aW)
      Tb <- T+input$c*(a-input$aB)
      Fw <- 1-input$k*(input$T0-Tw)^2
      Fb <- 1-input$k*(input$T0-Tb)^2
      for(j in c(1:length(Zones))){
        if(Fw[j]<0){Fw[j]=0}
        if(Fb[j]<0){Fb[j]=0}  }
      w <- w+w*(u*Fw-input$D)
      b <- b+b*(u*Fb-input$D)
      for(j in c(1:length(Zones))){
        if(w[j]<0.001){w[j]=0.001}
        if(b[j]<0.001){b[j]=0.001}  }
      u <- 1-w-b
      a <- w*input$aW+b*input$aB+u*alb(T,0.62,0.25)
      Barr[i] <- b[45]
      Warr[i] <- w[45]
      Uarr[i] <- u[45]
      I[i] <- i
      Tarr[i] <- T[45]
    } 
    
    plot_ly(z=~TEMP)%>% add_surface() %>% layout(title="With Daisies",
                                                 scene = list(
                                                   xaxis = list(title = "Time"),
                                                   yaxis = list(title = "Latitude"),
                                                   zaxis = list(title = "T")  )) 
  })
  
}


# APP ----
shinyApp(ui = ui, server = server)