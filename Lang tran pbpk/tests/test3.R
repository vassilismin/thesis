###### 1)Sto script xrisimopoiountai autousies oi exisoseis tou Tran kathos kai oi times tvn parametron. Ginetai elegxos ton exiseon me deSolve
###### 2)Allagh stis diaforikes gia ta tissues: prosthiki tvn lambda...2 ston oro eisodou


library(deSolve)
dose <- 1250
dE <- matrix(0, ncol=30)
E <- c(rep(0,29),dose)
inits <- as.vector(E)


##############
# Parameters #
##############

Vtis <- c(10.3, 6.0, 3.7, 1.2, 0.6, 1.2, (250- 10.3 - 6.0 - 3.7 - 1.2 - 0.6 - 1.2))   # From published studies -- Volume of tissue
Vcap <- c((0.06*Vtis[1]), (0.0265*Vtis[2]), (0.13*Vtis[3]), (0.1*Vtis[4]), (0.1*Vtis[5]), (0.033*Vtis[6]), (0.1*Vtis[7])) # From published studies -- Volume of capillary
Q   <- c(16992, 16704, 13248, 5616, 864, 1872, (120120 - 16992 - 16704 - 13248 - 5616 - 864 - 1872)) # From published studies -- Blood flow to tissue
QTotal <- sum(Q[])

lambdav <- 1
Da <- 0
Du <- 0
Do <- 0
Vven <- 5.6 
Vart <- 11.3
lambda3 <- c(0.09,0.0001,0.01,0)
lambda4 <- c(1,0,0.29185,0.5,0.9)
lambda5 <- c(0.9,0.05,0.09,0.0)
lambda6 <- c(0.07,0.9,0.003)
lambda7 <- c(0.08,0.003,0.08)
lambda8 <- c(0.0006,0.9,0.0011)
lambda9 <- c(0.9,0.9,0.003)
lambdaI  <- 1.0

fo <- 0
fu <- 0
ft <- 0

kB  <- 0
kb  <- 500
ko  <- 0
ku  <- 0
kt  <- 0
kr  <- 0
kd  <- 0
ki  <- 0
kl  <- 0.3

params <- c(Vtis[],Vcap[], Q[], QTotal, lambdav, Da, Do, Du, Vven, Vart, lambda3[], lambda4[], lambda5[], lambda6[], lambda7[], lambda8[], lambda9[], lambdaI, fo, fu, ft, kB, kb, ko, ku, kt, kr, kd, ki, kl)

###############
# ODEs system #
###############
ode.func <- function(time, inits, params){
  with( as.list(c(inits, params)),{
    
    dE[1]  = Do - (ko*E[1]) - (kB*E[1]);   # Olfactory
    dE[2]  = Du - (ku*E[2])  ;             # Upper airways
    dE[3]  = Da - (kr*E[3])  + (kd*E[4]) - (ki*E[3]);  # Alveolar free
    dE[4]  = (kr*E[3]) - (kd*E[4]) - (kt*E[4]) ;          # Alveolar Mac
    dE[5]  = (ki*E[3]) - (kl*E[5]) - (kb*E[5]) + ((E[30]/Vven)*lambdaI*QTotal*lambdav);    # I
    dE[6]  = kl*E[5];  # Lymph       to arterial  from ven
    
    
    #Liver
    dE[7]  = ((E[29]/Vart)*Q[3]) - ((E[7]/Vcap[3])*Q[3]) + ((E[8]/Vtis[3])*lambda3[1]*Q[3]) - ((E[7]/Vcap[3])*Q[3]*lambda3[2]) - 
      ((E[7]/Vcap[3])*Q[3]*lambda3[5]) + ((E[10]/Vcap[4])*Q[4]*lambda4[5]) ;  #capillary
    dE[8] =  ((E[7]/Vcap[3])*Q[3])*lambda3[2] - ((E[8]/Vtis[3])*lambda3[1])*Q[3] - (lambda3[3]*E[8])  ; #tissue
    dE[9]  = lambda3[3]*E[8]; #seq
    
    #GI
    dE[10]  = ((E[29]/Vart)*Q[4]) - ((E[10]/Vcap[4])*Q[4]) + ((E[11]/Vtis[4])*lambda4[1]*Q[4]) + ((E[7]/Vcap[3])*Q[3]*lambda3[5]) - ((E[10]/Vcap[4])*Q[4]*lambda4[5])  ;  #capillary
    dE[11]  = ((E[10]/Vcap[4])*Q[4]*lambda4[5]) - ((E[11]/Vtis[4])*lambda4[1]*Q[4]) - (lambda4[3]*E[11]) + (lambda4[6]*E[13]) ; #tissue
    dE[12]  = lambda4[3]*E[11]; #seq
    dE[13]  = (fo*ko*E[1]) + (fu*ku*E[2]) + (ft*kt*E[4]) - (lambda4[6]*E[13]) - (lambda4[4]*E[13]) ; # A4
    
    #Kidney
    dE[14]  = ((E[29]/Vart)*Q[5]) - ((E[14]/Vcap[5])*Q[5]) + ((E[15]/Vtis[5])*lambda5[1]*Q[5]) - ((E[14]/Vcap[5])*Q[5]*lambda5[2]) - (E[14]*lambda5[4]);  #capillary
    dE[15] =  ((E[14]/Vcap[5])*Q[5]*lambda5[2]) - ((E[15]/Vtis[5])*lambda5[1]*Q[5]) - (lambda5[3]*E[15])  ; #tissue
    dE[16]  = lambda5[3]*E[15]; #seq
    
    #Heart
    dE[17]  = ((E[29]/Vart)*Q[6]) - ((E[17]/Vcap[6])*Q[6]) + ((E[18]/Vtis[6])*lambda6[1]*Q[6]) - ((E[17]/Vcap[6])*Q[6]*lambda6[2]);  #capillary
    dE[18] =  ((E[17]/Vcap[6])*Q[6]*lambda6[2]) - ((E[18]/Vtis[6])*lambda6[1]*Q[6]) - (lambda6[3]*E[18]) ; #tissue
    dE[19]  = lambda6[3]*E[18]; #seq
    
    #Spleen
    dE[20]  = ((E[29]/Vart)*Q[7]) - ((E[20]/Vcap[7])*Q[7]) + ((E[21]/Vtis[7])*lambda7[1]*Q[7]) - ((E[20]/Vcap[7])*Q[7]*lambda7[2]);  #capillary
    dE[21] =  ((E[20]/Vcap[7])*Q[7]*lambda7[2]) - ((E[21]/Vtis[7])*lambda7[1]*Q[7]) - (lambda7[3]*E[21]) ; #tissue
    dE[22]  = lambda7[3]*E[21]; #seq
    
    # Brain
    dE[23]  =  ((E[29]/Vart)*Q[8]) - ((E[23]/Vcap[8])*Q[8]) + ((E[24]/Vtis[8])*lambda8[1]*Q[8]) - ((E[23]/Vcap[8])*Q[8]*lambda8[2]);  #capillary
    dE[24]  =  kB*E[1] + ((E[23]/Vcap[8])*Q[8]*lambda8[2]) - ((E[24]/Vtis[8])*lambda8[1]*Q[8]) - (lambda8[3]*E[24]) ; #tissue
    dE[25]  =  lambda8[3]*E[24]; #seq
    
    #Others
    dE[26]  = ((E[29]/Vart)*Q[9]) - ((E[26]/Vcap[9])*Q[9]) + ((E[27]/Vtis[9])*lambda9[1]*Q[9]) - ((E[26]/Vcap[9])*Q[9]*lambda9[2]);  #capillary
    dE[27] =  ((E[26]/Vcap[9])*Q[9]*lambda9[2]) - ((E[27]/Vtis[9])*lambda9[1]*Q[9]) - (lambda9[3]*E[27]); #tissue
    dE[28]  = lambda9[3]*E[27]; #seq
    
    #Blood
    dE[29] =  kb*E[5] - ((E[29]/Vart)*QTotal);   #art
    dE[30] =  ((E[23]/Vcap[8])*Q[8])  + ((E[7]/Vcap[3])*Q[3]) + ((E[17]/Vcap[6])*Q[6]) + ((E[20]/Vcap[7])*Q[7]) + 
      ((E[26]/Vcap[9])*Q[9]) + ((E[14]/Vcap[5])*Q[5]) - ((E[30]/Vven)*lambdaI*QTotal) ; #ven
    
    return(list(c(dE[1] , dE[2], dE[3], dE[4], dE[5], dE[6], dE[7], dE[8], dE[9], dE[10], dE[11], dE[12], dE[13], 
                  dE[14], dE[15], dE[16], dE[17], dE[18] ,dE[19], dE[20], dE[21], dE[22], dE[23], dE[24], dE[25],
                  dE[26], dE[27], dE[28], dE[29], dE[30])))
    
  })}

sample_time <- c(0, 1, 14, 28) #in days
#sample_time <- seq(0,1,0.01)
solution <- ode(times = sample_time, func = ode.func, y = inits, parms = params, method = "daspk")

solution