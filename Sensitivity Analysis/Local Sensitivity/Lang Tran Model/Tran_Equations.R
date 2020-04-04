#Sourced by "Tran Local Sensitivity.R"


dose <- 100 #micro grams

# Initial conditions. In this scenario we study the impact of an intravenous injection
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
solution <- ode(times = sample_time, func = ode.func, y = inits, parms = params, method = "bdf")

Total_amount <- matrix(0, nrow = 6, ncol = 6)

for (t in 1:6) { # t->time indice 
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

#solution
# Check that mass balance is satisfied
#rowSums(solution[,2:dim(solution)[2]])