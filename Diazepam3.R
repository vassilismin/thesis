library(deSolve)  

##############
# User input #
##############

weight <- 70 #in kg
gender <- 0 #0 for male | 1 for female
dose <- 10  # in mg
infusion_time <- 5/60 #infusion time in hours
#C0_MU <- 0
##C0_AD <- 0
#C0_GO <- 0
#C0_SK <- 0
#C0_HT <- 0
#C0_BR <- 0
#C0_KI <- 0
#C0_RE <- 0
#C0_ST <- 0
#C0_IN <- 0
#C0_LU <- 0
#C0_LI <- 0
#C0_ART <- 0
#C0_VEN <- 0
user_input <-data.frame(weight ,gender, dose, infusion_time)


##########################################
# Function for creating parameter vector #
##########################################

create.params <- function(input){
    w <- input$weight
    gender <- input$gender
    dose <- input$dose
    Tinf <- input$infusion_time

    k1<-c(9.61e-02,-4.88e-06,3.05e-10,-3.62e-15,1.22e-20,0,0.17)
    k2<-c(3.95e-02,1.59e-05,-6.99e-10,1.09e-14,-5.26e-20,0,0.05)
    k3<-c(1.67e-04,6.2e-10,-6.54e-13,2.48e-17,-2.85e-22,1.03e-27,0.001)
    k4<-c(1.07e-01,-3.26e-06,6.11e-11,-5.43e-16,1.83e-21,0,0.05)
    k5<-c(8.53e-03,-4.07e-07,1.4e-11,-1.9e-16,1.05e-21,-1.94e-27,0.04)
    k6<-c(1.19e-01,-3.51e-06,4.28e-11,-1.82e-16,0,0,0.12)
    k7<-c(7.31e-03,-8.29e-08,5.71e-13,0,0,0,0.19)
    k8<-rep(0,6)
    k9<-c(1.88e-03,8.76e-08,-2.52e-12,1.86e-17,0,0,0.01)
    k10<-c(1.74e-02,-5.3e-07,1.18e-11,-6.74e-17,0,0,0.14)
    k11<-c(1.67e-02,-9.96e-08,-1.09e-13,1.13e-17,0,0,1)
    k12<-c(3.49e-02,-3.23e-07,2.13e-12,0,0,0,0.065)
    k13<-c(3.66e-02,-3.44e-07,5.00e-12, -2.59e-17,0.0,0.0)
    k14<-c(5.49e-02,-5.15e-07, 7.50e-12,-3.87e-17,0.0,0.0)
  
    k1_W<-c(1.17e-01,-3.59e-06,3.19e-10,-3.55e-15,-7.58e-22,0.0)
    k2_W<-c(5.91e-02,1.20e-05,-5.80e-10,1.12e-14,-6.36e-20,0.0)
    k3_W<-c(1.94e-04,-8.32e-09,3.15e-13,0.0,0.0,0.0)
    k4_W<-c(9.54e-02,-1.7e-06,-1.64e-13,2.64e-16,-1.49e-21,0.0)
    k5_W<-c(5.72e-03,-1.02e-07,2.53e-12,-2.71e-17,9.29e-23,0.0)
    k6_W<-c(1.12e-01,-3.33e-06,4.04e-11,-1.70e-16,0.0,0.0)
    k7_W<-c(8.04e-03,-1.38e-07,2.19e-12,-1.34e-17,0.0,0.0)
    k8_W<-rep(0,6)
    k9_W<-c(1.88e-03,8.76e-08,-2.52e-12,1.86e-17,0.0,0.0)
    k10_W<-c(1.89e-02,-6.62e-07,1.56e-11,-9.87e-17,0.0,0.0)
    k11_W<-c(1.74e-02,-7.14e-08,-6.78e-14,0.0,0.0,0.0)
    k12_W<-c(3.59e-02,-4.76e-07,8.50e-12,-5.45e-17,0.0,0.0)
    k13_W<-c(3.66e-02,-3.44e-07,5.00e-12, -2.59e-17,0.0,0.0)
    k14_W<-c(5.49e-02,-5.15e-07, 7.50e-12,-3.87e-17,0.0,0.0)
  
    #density<-c(1.041,0.916,1,1,1.03,1.035,1.05,1,1.05,1.042,1.05,1,1,1)#in kg/L
    density <- rep(1,14)
    flow_frac<-c(0.17,0.05,0.001,0.05,0.04,0.12,0.19,0,0.01,0.14,1,0.065)
  
    const<-list(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14)
    const_W<-list(k1_W,k2_W,k3_W,k4_W,k5_W,k6_W,k7_W,k8_W,k9_W,k10_W,k11_W,k12_W,k13_W,k14_W)
    vol<-rep(0,14)
    flow<-rep(0,12)
    TF<-(187*(w)^0.81)*60/1000 #Total flow
    weight <- w
  
    if (gender==0){
      if (w>75){
        w<-75*1000
      }else{
        w<-w*1000
      }
      co<-const
      for (i in 1:14){
        vol[i]<-(co[[i]][1]+co[[i]][2]*w+co[[i]][3]*w^2+co[[i]][4]*w^3+co[[i]][5]*w^4+co[[i]][6]*w^5)
        vol[i]<-vol[i]*weight/density[i]
      }
    }else if (gender==1){
      if (w>65){
        w<-65*1000
      }else{
        w<-w*1000
      }
      co<-const_W
      for (i in 1:14){
        vol[i]<-(co[[i]][1]+co[[i]][2]*w+co[[i]][3]*w^2+co[[i]][4]*w^3+co[[i]][5]*w^4+co[[i]][6]*w^5)
        vol[i]<-vol[i]*weight/density[i]
      }
    }
    for (i in 1:12){
      flow[i]<-TF*flow_frac[i]
    }
    flow[8]<-TF-sum(flow)+TF
    vol[8]<-weight-sum(vol)
    combine<-c(flow,vol)
  
    return(  c("Q_MU" = combine[1], "Q_AD" = combine[2], "Q_TE" = combine[3], "Q_SK" = combine[4],
               "Q_HT" = combine[5],"Q_BR" = combine[6], "Q_KI" = combine[7],  "Q_RE" = combine[8], 
               "Q_ST" = combine[9], "Q_SPL" = combine[10], "Q_LU" = combine[11], "Q_LI"=combine[12], 
               "Q_ART" = combine[11], "Q_VEN" = combine[11], "V_MU" = combine[13], "V_AD" = combine[14], 
               "V_TE" = combine[15], "V_SK" = combine[16], "V_HT" = combine[17],"V_BR" = combine[18], 
               "V_KI" = combine[19], "V_RE"=combine[20], "V_ST"=combine[21],"V_SPL"=combine[22],
               "V_LU"=combine[23], "V_LI"=combine[24], "V_ART"=combine[25], "V_VEN"=combine[26], 
             "dose" = dose*1e06, "Tinf" = Tinf, "weight" = w, "gender" = gender))
  
}

### store them once
params <- create.params(user_input)

#################################################
# Function for creating initial values for ODEs #
#################################################

create.inits <- function(parameters){
  with( as.list(parameters),{
    Mu <- 0; Ad <- 0; Go <- 0; Sk <- 0; Ht <- 0; Br <- 0; Ki <- 0; Re <- 0; St <- 0; Spl <- 0; Lu <- 0; 
    Li <- 0; Art <- 0; Ven <- 0; 
    
    return(c("Mu" = Mu, "Ad" =  Ad,  "Go" = Go, "Sk" = Sk, "Ht" = Ht, "Br" = Br, "Ki" = Ki, "Re" = Re,
             "St" = St, "Spl" =Spl, "Lu" = Lu, "Li" = Li, "Art" = Art, "Ven" = Ven))
  }) 
}
##store the values
inits <- create.inits(params)

#################################################
# Function for creating events #
#################################################

create.events<- function(parameters){
    with( as.list(parameters),{
        parts <- 5
        addition <- dose/(V_VEN*1000*parts)
        addition2 <- dose/(V_VEN*1000)
        interval <- (100-70)/10
        start <- 70
        end <- 100
        events <- list(data = rbind(data.frame(var = c("Ven"),  time = seq(0, Tinf , by=Tinf/parts), 
                                         value = addition, method = c("add")), 
                                    
                                    data.frame(var = c("Ven"),  time = seq(start, end , by=interval), 
                                                            value = addition2, method = c("add")) ) )
        return(events)
  }) 
}

events <- create.events(params)

###################
# Custom function #
###################

custom.func <- function(x){
  if (x == 1){
    y = 64.5
  } else{
    y=100
  }
  return(y)
}
#################
# ODEs system   #
#################

ode.func <- function(time, Initial.values, Parameters, custom.func){
  with( as.list(c(Initial.values, Parameters)),{
    
  Kp_init <- c(0.56, 3.34, 0.76, 0.46, 0.86, 0.32, 0.73, 0.5, 0.75, 0.53, 0.71, 1.35) 
  a <- 0.12 # scaling factor
  Kp <- Kp_init * a # updated Kp vector
  Kp[2] <- 3.32 # Kp of adipose compartment
  
  KP_MU <- Kp[1]
  KP_AD <- Kp[2]
  KP_TE <- Kp[3]
  KP_SK <- Kp[4]
  KP_HT <- Kp[5]
  KP_BR <- Kp[6]
  KP_KI <- Kp[7]
  KP_RE <- Kp[8]
  KP_ST <- Kp[9]
  KP_SPL <- Kp[10]
  KP_LU <- Kp[11]
  KP_LI <- Kp[12]
  CL <- custom.func(1)
  Q_HA <-Q_LI
  Q_HP <-Q_ST + Q_SPL
  Q_H <- Q_HA+Q_HP
  f_UB <- 0.015/0.65 
  
  dMu <- (-Q_MU*Mu/(KP_MU*V_MU))+(Q_MU*Art/V_MU)
  dAd <- (-Q_AD*Ad/(KP_AD*V_AD))+(Q_AD*Art/V_AD)
  dGo <- (-Q_TE*Go/(KP_TE*V_TE))+(Q_TE*Art/V_TE)
  dSk <- (-Q_SK*Sk/(KP_SK*V_SK))+(Q_SK*Art/V_SK)
  dHt <- (-Q_HT*Ht/(KP_HT*V_HT))+(Q_HT*Art/V_HT)
  dBr <- (-Q_BR*Br/(KP_BR*V_BR))+(Q_BR*Art/V_BR)
  dKi <- (-Q_KI*Ki/(KP_KI*V_KI))+(Q_KI*Art/V_KI)
  dRe <- (-Q_RE*Re/(KP_RE*V_RE))+(Q_RE*Art/V_RE)
  dSt <- (-Q_ST*St/(KP_ST*V_ST))+(Q_ST*Art/V_ST)
  dSpl <- (-Q_SPL*Spl/(KP_SPL*V_SPL))+(Q_SPL*Art/V_SPL)
  dLu <- (-Q_LU*Lu/(KP_LU*V_LU))+(Q_VEN*Ven/V_LU)
  dLi <- (Q_HA*Art/V_LI)+(Q_ST*St/(KP_ST*V_LI))+(Q_SPL*Spl/(KP_SPL*V_LI))-(((Q_H+CL*f_UB)/(KP_LI*V_LI))*Li)
  dArt <- (-Q_ART*Art/V_ART)+(Q_LU*Lu/(V_ART*KP_LU))
  dVen <- (-Q_LU*Ven/V_VEN)+(Q_H*Li/(KP_LI*V_VEN))+(Q_KI*Ki/(KP_KI*V_VEN))+
      (Q_MU*Mu/(KP_MU*V_VEN))+(Q_AD*Ad/(KP_AD*V_VEN))+(Q_SK*Sk/(KP_SK*V_VEN))+
      (Q_TE*Go/(KP_TE*V_VEN))+(Q_HT*Ht/(KP_HT*V_VEN))+(Q_BR*Br/(KP_BR*V_VEN))+
      (Q_RE*Re/(KP_RE*V_VEN))
  
  
    
  list(c(dMu= dMu, dAd = dAd, dGo = dGo, dSk = dSk, dHt = dHt, dBr = dBr, dKi = dKi, dRe = dRe, dSt = dSt,
         dSpl = dSpl, dLu = dLu, dLi = dLi, dArt = dArt, dVen = dVen))       
  }) }

##############################################

sample_time <- c(0, 5/60, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 24, 36, 48, 72, 100) # in hours
sample_time <- seq(0,100, by=1)
solution <- ode(times = sample_time,  func = ode.func, y = inits, parms = params, 
                custom.func = custom.func, method="lsodes",  events = events)
solution

plot(solution[,1],solution[,15], type = "l")
#deploy.pbpk(user_input, predicted.feats, create.params, create.inits, create.events, 
 #            custom.func, ode.func)