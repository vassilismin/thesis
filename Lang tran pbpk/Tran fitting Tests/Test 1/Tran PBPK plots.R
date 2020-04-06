library(ggplot2)
library(rstan)
library(reshape2)
library(dplyr)
library(deSolve)
library(MASS)
library(openxlsx)
options(max.print=1000000)

setwd("C:\\Users\\vassi\\Documents\\Diploma Thesis\\Lang Tran Files\\R  Stan codes\\Tran fitting Tests\\Test 1")

load("Results.RData")

exp_fit <- extract(fit)

dose <- 1.7e03 #micro grams
times <- c(1e-20) #infusion time in hours



##########################################
# Function for creating parameter vector #
##########################################

    ##############
    # Parameters #
    ##############
    # In our code
    #1: Liver, 2:gut, 3: kidney, 4:heart, 5:spleen, 6:brain, 7:other tissues
    # In Tran's cod7
    # 3: liver, 4: gut, 5: kidney, 6:heart, 7:spleen, 8:brain, 9:other tissues
    
    # Rat weight
    Wrat <- 250 #mL
    # Tissue volumes in mL
    Vtis <- c(10.3, 6.0, 1.2, 1.2, 0.6, 1.2, NA)
    Vven <- 5.6 
    Vart <- 11.3
    Vtis[length(Vtis)] <- Wrat - sum(Vtis, na.rm = TRUE)-Vven-Vart
    #Fraction of capillary blood
    fcap <- c(0.06, 0.0265, 0.13, 0.1, 0.1, 0.033, 0.1)
    # Volume of capillary blood per tissue
    Vcap <- fcap*Vtis
    
    
    # Regional blood flows in mL/day
    QC <- 120120 # Cardiac output
    # The first flow is from hepatic artery to liver (liver input)
    #Q   <- c(2523, 16704, 13248, 5616, 864, 1872, NA) 
    Q   <- c(2523, 16704, 13248, 5616, 864, 1872, NA)  
    Q[length(Q)]<- QC - sum(Q, na.rm = TRUE)
    Q_liver <- Q[1]+Q[2]+Q[5] # Liver output
    QTotal <- QC
    Q_bile <- 20 # Daily bile production by the liver 
    
    Depo <- 0 # Olfactory deposition fraction
    Depu <- 0.44 # Upper airways deposition fraction
    Depa <- 0.3 # Alveolar deposition fraction
    conc <- 0 # External concentration
    VI <- 0.18 # Volume inhaled (L/min)
    CONV <- 1.44 # Conversion factor
    
    Da <- Depo * conc * VI * CONV # Olfactory deposition rate [micro grams/day]
    Du <- Depu * conc * VI * CONV # Upper airways deposition rate [micro grams/day]
    Do <- Depa * conc * VI * CONV # Alveolar deposition rate [micro grams/day]
    
    #####################
    # Priors #
    #####################
    
    # lamdai_1: plasma to tissue PC [unitless]
    # lamdai_2: fractional translocation from capillaries to venous blood [unitless]
    # lamdai_3: sequestration rate [d^-1]  
    # lamda3_5: bile:liver PC [unitless]
    # lamda4_5: GI to liver PC [unitless]
    # lamda4_6 GI absorption rate [d^-1]  
    # lamda4_4: fecal elimination rate [d^-1]  
    
    lambda3 <- c(4.39, 0.0001, 0.24, 1e-5, 1.73)  #lambda3[5] taken from code for semmler
    lambda4 <- c(2, 0.29185, 0.07, 4, 1, 1.44) #lambda4[6] taken from code for semmler
    lambda5 <- c(0.59, 0.05, 0.73, 1e-5)
    lambda6 <- c(0.79, 0.9, 0.69)
    lambda7 <- c(0.37, 0.003, 0.12)
    lambda8 <- c(1.33, 0.9, 0.58)
    lambda9 <- c(6.62, 0.9, 0.18)
    
    
    k_B  <- 1e-5 # Translocation from Olfactory to brain [d^-1]
    kb  <-  0.29 # Translocation from interstitium to blood [d^-1] 
    ko  <- 1e-5 # Clearance from Olfactory to GI [d^-1]
    ku  <- 109 # Clearance from upper airwais to GI [d^-1]
    kt  <- 0.015 # Clearance from upper airwais to GI [d^-1]
    kr  <- 4 # Macrophage phagocytosis rate [d^-1]
    kd  <- 0.033 # Macrophage death rate [d^-1]
    ki  <- 3.5 # interstitialization rate [d^-1]
    kl  <- 0.07 # Translocation from interstitium to Lymph nodes [d^-1]
    lambdaI  <- 0.15 # fraction of translocation from blood to interstitium [unitless]
    
    ####Parameters after fitting
    for (i in 1:37) {
      if (i <=5){
        lambda3[i] = mean(exp_fit$theta[,i])
      } else if ((i>=6) &  (i<=11)){
        lambda4[i-5] = mean(exp_fit$theta[,i])
      } else if ((i>=12) & (i<=15)){
        lambda5[i-11] = mean(exp_fit$theta[,i])
      } else if ((i>=16) & (i<=18)){
        lambda6[i-15] = mean(exp_fit$theta[,i])
      } else if ((i>=19) & (i<=21)){
        lambda7[i-18] = mean(exp_fit$theta[,i])
      } else if ((i>=22) & (i<=24)){
        lambda8[i-21] = mean(exp_fit$theta[,i])
      } else if ((i>=25) & (i<= 27)){
        lambda9[i-24] = mean(exp_fit$theta[,i])
      }
      
      k_B <- mean(exp_fit$theta[,28])
      kb <- mean(exp_fit$theta[,29])
      ko <- mean(exp_fit$theta[,30])
      ku <- mean(exp_fit$theta[,31])
      kt <- mean(exp_fit$theta[,32])
      kr <- mean(exp_fit$theta[,33])
      kd <- mean(exp_fit$theta[,34])
      ki <- mean(exp_fit$theta[,35])
      kl <- mean(exp_fit$theta[,36])
      lambdaI <- mean(exp_fit$theta[,37])
    }
    
    params <- c(lambda3[], lambda4[], lambda5[], lambda6[], lambda7[], lambda8[], lambda9[], k_B, kb, ko, ku, kt, kr, kd, ki, kl, lambdaI, 
                Vtis[],Vcap[], Q[], QTotal, Da, Do, Du, Vven, Vart)
    

#########################
### Initial conditions###
#########################
    
M <- c(rep(0,29), dose, 0, 0)
inits <- as.vector(M)

###############
# ODEs system #
###############

ode.func <- function(time, M, params){
  with( as.list(params),{
    
    dMdt<-rep(0,32)
    
    dMdt[1]  = Do - (ko*M[1]) - (k_B*M[1]) ;   # Olfactory
    dMdt[2]  = Du - (ku*M[2])  ;             # Upper airways
    dMdt[3]  = Da - (kr*M[3])  + (kd*M[4]) - (ki*M[3]) ;  # Alveolar free
    dMdt[4]  = (kr*M[3]) - (kd*M[4]) - (kt*M[4]) ;          # Alveolar Mac
    dMdt[5]  = (ki*M[3]) - (kl*M[5]) - (kb*M[5]) + ((M[30]/Vven)*lambdaI*QTotal);    # I
    dMdt[6]  = (kl*M[5]) ;  # Lymph       to arterial  from ven
    
    #Liver
    dMdt[7]  = ((M[29]/Vart)*Q[1]) - ((M[7]/Vcap[1])*Q[1]) + ((M[8]/Vtis[1])*lambda3[1]*Q[1]) - ((M[7]/Vcap[1])*Q_liver*lambda3[2]) - 
      ((M[7]/Vcap[1])*Q_bile*lambda3[5]) + ((M[10]/Vcap[2])*Q[2]*lambda4[5]) + ((M[20]/Vcap[5])*Q[5]) ;  #capillary
    dMdt[8] =  ((M[7]/Vcap[1])*Q[1]) - ((M[8]/Vtis[1])*lambda3[1]*Q[1]) - (lambda3[3]*M[8])  ; #tissue
    dMdt[9]  = (lambda3[3]*M[8]) ; #seq
    
    #GI
    dMdt[10]  = ((M[29]/Vart)*Q[2]) - ((M[10]/Vcap[2])*Q[2]) + ((M[11]/Vtis[2])*lambda4[1]*Q[2])  - ((M[10]/Vcap[2])*Q[2]*lambda4[5])  ;  #capillary
    dMdt[11]  = ((M[10]/Vcap[2])*Q[2]) - ((M[11]/Vtis[2])*lambda4[1]*Q[2]) - (lambda4[3]*M[11]) + (lambda4[6]*M[13]) ; #tissue
    dMdt[12]  = (lambda4[3]*M[11]); #seq
    dMdt[13]  = (ko*M[1]) + (ku*M[2]) + (kt*M[4]) - (lambda4[6]*M[13]) - (lambda4[4]*M[13])+ ((M[7]/Vcap[1])*Q_bile*lambda3[5]) ; # A4
    
    #Kidney
    dMdt[14]  = ((M[29]/Vart)*Q[3]) - ((M[14]/Vcap[3])*Q[3]) + ((M[15]/Vtis[3])*lambda5[1]*Q[3]) - ((M[14]/Vcap[3])*Q[3]*lambda5[2]) - 
      (M[14]*lambda5[4]);  #capillary
    dMdt[15] =  ((M[14]/Vcap[3])*Q[3]) - ((M[15]/Vtis[3])*lambda5[1]*Q[3]) - (lambda5[3]*M[15])  ; #tissue
    dMdt[16]  = (lambda5[3]*M[15]) ; #seq
    
    #Heart
    dMdt[17]  = ((M[29]/Vart)*Q[4]) - ((M[17]/Vcap[4])*Q[4]) + ((M[18]/Vtis[4])*lambda6[1]*Q[4]) - ((M[17]/Vcap[4])*Q[4]*lambda6[2]) ;  #capillary
    dMdt[18] =  ((M[17]/Vcap[4])*Q[4]) - ((M[18]/Vtis[4])*lambda6[1]*Q[4]) - (lambda6[3]*M[18]) ; #tissue
    dMdt[19]  = (lambda6[3]*M[18]) ; #seq
    
    #Spleen
    dMdt[20]  = ((M[29]/Vart)*Q[5]) - ((M[20]/Vcap[5])*Q[5]) + ((M[21]/Vtis[5])*lambda7[1]*Q[5]) - ((M[20]/Vcap[5])*Q[5]) ;  #capillary
    dMdt[21] =  ((M[20]/Vcap[5])*Q[5]) - ((M[21]/Vtis[5])*lambda7[1]*Q[5]) - (lambda7[3]*M[21]) ; #tissue
    dMdt[22]  = (lambda7[3]*M[21]) ; #seq
    
    # Brain
    dMdt[23]  =  ((M[29]/Vart)*Q[6]) - ((M[23]/Vcap[6])*Q[6]) + ((M[24]/Vtis[6])*lambda8[1]*Q[6]) - ((M[23]/Vcap[6])*Q[6]*lambda8[2]) ;  #capillary
    dMdt[24]  =  (k_B*M[1]) + ((M[23]/Vcap[6])*Q[6]) - ((M[24]/Vtis[6])*lambda8[1]*Q[6]) - (lambda8[3]*M[24]) ; #tissue
    dMdt[25]  =  (lambda8[3]*M[24]) ; #seq
    
    #Others
    dMdt[26]  = ((M[29]/Vart)*Q[7]) - ((M[26]/Vcap[7])*Q[7]) + ((M[27]/Vtis[7])*lambda9[1]*Q[7]) - ((M[26]/Vcap[7])*Q[7]*lambda9[2]) ;  #capillary
    dMdt[27] =  ((M[26]/Vcap[7])*Q[7]) - ((M[27]/Vtis[7])*lambda9[1]*Q[7]) - (lambda9[3]*M[27]) ; #tissue
    dMdt[28]  = (lambda9[3]*M[27]) ; #seq
    
    #Blood
    dMdt[29] =  (kb*M[5]) - ((M[29]/Vart)*Q[7])- ((M[29]/Vart)*Q[6])- ((M[29]/Vart)*Q[5])- ((M[29]/Vart)*Q[4])- ((M[29]/Vart)*Q[3])-
      ((M[29]/Vart)*Q[2])- ((M[29]/Vart)*Q[1]);   #art
    dMdt[30] =  ((M[23]/Vcap[6])*Q[6]*lambda8[2])  + ((M[7]/Vcap[1])*Q_liver*lambda3[2]) + ((M[17]/Vcap[4])*Q[4]*lambda6[2]) + 
      ((M[26]/Vcap[7])*Q[7]*lambda9[2]) + ((M[14]/Vcap[3])*Q[3]*lambda5[2]) - ((M[30]/Vven)*lambdaI*QTotal) ; #ven
    
    #Feces
    dMdt[31] = (lambda4[4]*M[13]) ;
    
    #Urine
    dMdt[32] =  (M[14]*lambda5[4]) ;
    
    return(list(dMdt))
    
  })}

sample_time <- c(10/(60*24), 1/24, 1, 7, 28, 56) #in days
#sample_time <- seq(0, 56, 0.001)
ltime <- length(sample_time)

solution <- ode(times = sample_time, func = ode.func, y = inits, parms = params, method = "bdf", maxsteps = 1000000)
solution <- solution[,2:33]

Total_amount <- matrix(0, nrow = ltime, ncol = 6)

for (t in 1:ltime) { # t->time indice 
  ###Total amount of NPs in each organ
  
  # Amount in Lungs
  Total_amount[t,1] <- solution[t,1] + solution[t,2] + solution[t,3] + solution[t,4] + solution[t,5] + solution[t,6]
  
  # Amount in Brain
  Total_amount[t,2] <- solution[t,23] + solution[t,24] + solution[t,25]
  
  # Amount in Kidneys
  Total_amount[t,3] <- solution[t,14] + solution[t,15] + solution[t,16]
  
  # Amount in Liver
  Total_amount[t,4] <- solution[t,7] + solution[t,8] + solution[t,9]
  
  # Amount in Spleen
  Total_amount[t,5] <- solution[t,20] + solution[t,21] + solution[t,22]
  
  # Amount in Blood
  Total_amount[t,6] <- solution[t,29] +solution[t,30] 
}

M_Total <- as.data.frame(cbind(sample_time,Total_amount[,]))
colnames(M_Total) <- c("Time", "Lungs", "Brain", "Kidneys", "Liver", "Spleen", "Blood")

M_Capillaries <- as.data.frame(cbind(sample_time, solution[,23], solution[,14], solution[,7], solution[,20]))
colnames(M_Capillaries) <- c("Time", "Brain", "Kidneys", "Liver", "Spleen")

M_Tissue <- as.data.frame(cbind(sample_time, Total_amount[,1], solution[,24], solution[,15], solution[,8], solution[,21], solution[,29], solution[,30]))
colnames(M_Tissue) <- c("Time", "Lungs", "Brain", "Kidneys", "Liver", "Spleen", "Art", "Ven")

M_Sequestered <- as.data.frame(cbind(sample_time, solution[,25], solution[,16], solution[,9], solution[t,22]))
colnames(M_Sequestered) <- c("Time", "Brain", "Kidneys", "Liver", "Spleen")

ggplot(M_Total, aes(x=Time, y=Liver, colour = "Total"))+
  geom_line(size=1.2) +
  geom_line(data=M_Capillaries, aes(x=Time, y=Liver, colour = "Capillaries"),size=1.2)+
  geom_line(data = M_Tissue, aes(x = Time, y = Liver, colour = "Tissue"),size=1.2)+
  geom_line(data = M_Sequestered, aes(x = Time, y= Liver, colour = "Sequestered"),size=1.2)+
  
  labs(x = expression("Time (days)"),y = expression("Amount(micrograms)"))+
  scale_colour_manual(name="Subcompartment", values=c("Sequestered"=4 ,"Tissue"=3, "Capillaries"=2, "Total"=1))+
  ggtitle("Liver")+
  theme(plot.title = element_text(hjust = 0.5,size=26), axis.title.y =element_text(hjust = 0.5,size=20,face="bold"),
        axis.text.y=element_text(size=18),
        axis.title.x =element_text(hjust = 0.5,size=20,face="bold"),
        axis.text.x=element_text(size=18),
        legend.title=element_text(hjust = 0.5,size=20), 
        legend.text=element_text(size=18))





