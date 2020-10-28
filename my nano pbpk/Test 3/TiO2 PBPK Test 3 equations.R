library(openxlsx)
library(deSolve)
library(ggplot2)

setwd("C:\\Users\\vassi\\Documents\\Diploma Thesis\\My nano PBPK\\Final Model\\Fitting\\Test 3")

####################################################################
###################             DATA             ###################
####################################################################

doses <- c(18.15, 11.29, 16.53, 108.5, 46.92)

### Important!!! each compartment has a specific index vectors Tissue_fractions, Regional_flow_fractions, Capillary_fractions and cannot be changed
# The index of each compartment:
#Rest of Body (rob) --> 1
#Heart (ht) --> 2
#Kidneys (ki) --> 3
#Brain (br) --> 4
#Spleen (spl) --> 5
#Lungs (lu) --> 6
#Liver (li) --> 7
#Uterus (ut) --> 8
#Skeleton (skel) --> 9
#Adipose (ad) --> 10
#Skin (skin) --> 11
#Muscles (mu) --> 12


####################
### User's INPUT ###
####################
#### If any of these compartments don not exist in pbpk, just give it the value NA in compartments vector, example: "Heart" = NA and it will remove it 
#### from the equilibriums and the corresponding V_tis, V_cap, Q will be equal to NA.


compartments <- list( "RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"="Brain", "Spleen"="Spleen",
                      "Lungs"="Lungs", "Liver"="Liver", "Uterus"="Uterus", "Skeleton"="Skeleton", "Adipose"="Adipose", "Skin"="Skin", "Muscles"="Muscles") #used as input in function, compartments that are used in pbpk

BW <- 263 # Total Body weight of rat in g

#####################################
### Function to create Parameters ###
#####################################

create.params <- function(comp_names, w){
  
  # List with names of all possible compartments
  all_comps <- list("RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"="Brain", "Spleen"="Spleen",
                    "Lungs"="Lungs", "Liver"="Liver", "Uterus"="Uterus", "Skeleton"="Skeleton", "Adipose"="Adipose", "Skin"="Skin", "Muscles"="Muscles") # List with names of all possible compartments
  
  ### Density of tissues/organs
  d_tissue <- 1 #g/ml
  d_skeleton <- 1.92 #g/ml
  d_adipose <- 0.940 #g/ml
  
  Q_total <- (1.54*w^0.75)*60 # Total Cardiac Output (ml/h)
  
  Total_Blood <- 0.06*w+0.77 # Total blood volume (ml)
  
  #Arterial blood volume
  Vart <- 1.2905*w/100 #0.15*Total_Blood #(ml)
  
  #Veins blood volume
  Vven <- 2.968*w/100 #0.64*Total_Blood #(ml)
  
  fr_ad <- 0.0199*w + 1.644 # w in g,  Brown et al.1997 p.420. This equation gives the  adipose % of body weight 
  
  #read data from excel
  fractions <- openxlsx::read.xlsx("Rat physiological parameters.xlsx", sheet = 1, colNames = T, rowNames = T)
  fractions <- as.matrix(sapply(fractions, as.numeric))
  
  #Tissue weight fraction 
  Tissue_fractions <- fractions[,1]/100 # % of BW. Na values refers to the volume of the rest organs(RoB)
  Tissue_fractions[10] <- fr_ad/100
  #Regional blood flow fraction
  Regional_flow_fractions <- fractions[,2]/100 # % of total cardiac output
  #Capillary volume fractions (fractions of tissue volume)
  Capillary_fractions <- fractions[,3] # of tissue volume
  #Macrophage content as fraction tissue volume for each tissue/organ
  Macrophage_fractions <- fractions[,4] 
  
  W_tis <- rep(0,length(comp_names))
  V_tis <- rep(0,length(comp_names))
  V_cap <- rep(0,length(comp_names))
  V_macro <- rep(0,length(comp_names))  #one more for blood compartment
  Q <- rep(0,length(comp_names))
  
  
  for (i in 1:length(comp_names)) {
    control <- comp_names[i]
    
    Tissue_fractions[i] <- ifelse(is.na(control), NA, Tissue_fractions[i])
    Regional_flow_fractions[i] <- ifelse(is.na(control), NA, Regional_flow_fractions[i])
    Capillary_fractions[i] <- ifelse(is.na(control), NA, Capillary_fractions[i])
    Macrophage_fractions[i] <- ifelse(is.na(control), NA, Macrophage_fractions[i])
    
    ### Calculation of tissue weights  
    W_tis[i] <- w*Tissue_fractions[i]
    
    
    ###Calculation of tissue volumes
    
    if (i==9){
      V_tis[i] <- W_tis[i]/d_skeleton
    } else if(i==10){
      V_tis[i] <- W_tis[i]/d_adipose
    } else{
      V_tis[i] <- W_tis[i]/d_tissue 
    }
    
    ###Calculation of capillary volumes
    V_cap[i] <- V_tis[i]*Capillary_fractions[i]
    
    ###Volume of macrophage contents
    V_macro[i] <- V_tis[i]*Macrophage_fractions[i]
    
    ###Calculation of regional blood flows
    Q[i] <- Q_total*Regional_flow_fractions[i]
  }
  
  Vm_ven <- 0.02*Vven #macrophage content in veins
  Vm_art <- 0.01*Vart #0.02*Vart #macrophage content in arteries
  
  ### Calculations for "Rest of Body" compartment
  W_tis[1] <- w - sum(W_tis[2:length(W_tis)], na.rm = TRUE)
  V_tis[1] <- W_tis[1]     #(considering that the density of the rest tissues is 1 g/ml)
  Q[1] <- Q_total - sum(Q[2:length(Q)],na.rm = TRUE) + Q[6]
  V_cap[1] <- Total_Blood - Vven - Vart - sum(V_cap[2:length(V_cap)], na.rm = TRUE)
  V_macro[1] <- V_tis[1]*Macrophage_fractions[1]
  #Capillary_fractions[1] <- V_cap[1]/V_tis[1]
  
  
  parameters <- matrix(c(W_tis[],V_tis[],V_cap[],Q[],V_macro[]), ncol = 5)
  colnames(parameters) <- c("W_tis", "V_tis", "V_cap", "Q", "V_macro")
  rownames(parameters) <- all_comps
  
  
  return(list(parameters,
              "Q_total"=Q_total, "V_blood"=Total_Blood, "Vven"=Vven, "Vart"=Vart, "Vm_ven"=Vm_ven, "Vm_art"=Vm_art,
              
              "W_rob"=parameters[1,1], "W_ht"=parameters[2,1], "W_ki"=parameters[3,1], "W_br"=parameters[4,1], "W_spl"=parameters[5,1], "W_lu"=parameters[6,1], "W_li"=parameters[7,1], "W_ut"=parameters[8,1], "W_skel"=parameters[9,1],
              "W_ad"=parameters[10,1], "W_skin"=parameters[11,1], "W_mu"=parameters[12,1], 
              
              "Vtis_rob"=parameters[1,2], "Vtis_ht"=parameters[2,2], "Vtis_ki"=parameters[3,2], "Vtis_br"=parameters[4,2], "Vtis_spl"=parameters[5,2], "Vtis_lu"=parameters[6,2], "Vtis_li"=parameters[7,2], "Vtis_ut"=parameters[8,2], "Vtis_skel"=parameters[9,2],
              "Vtis_ad"=parameters[10,2], "Vtis_skin"=parameters[11,2], "Vtis_mu"=parameters[12,2],
              
              "Vcap_rob"=parameters[1,3], "Vcap_ht"=parameters[2,3], "Vcap_ki"=parameters[3,3], "Vcap_br"=parameters[4,3], "Vcap_spl"=parameters[5,3], "Vcap_lu"=parameters[6,3], "Vcap_li"=parameters[7,3], "Vcap_ut"=parameters[8,3], "Vcap_skel"=parameters[9,3],
              "Vcap_ad"=parameters[10,3], "Vcap_skin"=parameters[11,3], "Vcap_mu"=parameters[12,3],
              
              "Vm_rob"=parameters[1,5], "Vm_ht"=parameters[2,5], "Vm_ki"=parameters[3,5], "Vm_br"=parameters[4,5], "Vm_spl"=parameters[5,5], "Vm_lu"=parameters[6,5], "Vm_li"=parameters[7,5], "Vm_ut"=parameters[8,5], "Vm_skel"=parameters[9,5],
              "Vm_ad"=parameters[10,5], "Vm_skin"=parameters[11,5], "Vm_mu"=parameters[12,5],
              "Q_rob"=parameters[1,4]+parameters[6,4], "Q_ht"=parameters[2,4], "Q_ki"=parameters[3,4], "Q_br"=parameters[4,4], "Q_spl"=parameters[5,4], "Q_lu"=parameters[6,4], "Q_li"=parameters[7,4], "Q_ut"=parameters[8,4], "Q_skel"=parameters[9,4],
              "Q_ad"=parameters[10,4], "Q_skin"=parameters[11,4], "Q_mu"=parameters[12,4],
              
              
              P_1=6, #liver
              P_2=0.55, #spleen
              P_3=0.01, #kidneys, lungs, blood
              P_4=0.1, #heart, uterus
              P_5=0.005, #Skeleton, Soft tissue 
              P_6=1e-10, #brain
              
              x_1=1000, #liver
              x_2=1000, #spleen
              x_3=0.0001, #kidneys, lungs, blood
              x_4=0.001, #heart, uterus
              x_5=10, #Skeleton, Soft tissue
              x_6=1e-10, #brain,
              k_de=0.0001,
              Km=5, #micro-g/ml 
              Pup=0.995, #ml/h/(ml pcs)
              
             
              CLE_f=0.0008, #1/h
              CLE_u=1.7 #1/h
              
              
  ))
}

params<-create.params(compartments,BW)
#print(params)
dose <- doses[1]
inits <- c(Mcap_lu=0, Mcap_spl=0, Mcap_li=0, Mcap_ki=0, Mcap_ht=0, Mcap_br=0, Mcap_ut=0, Mcap_skel=0, Mcap_st=0, Mcap_rob=0,
           Mtis_lu=0, Mtis_spl=0, Mtis_li=0, Mtis_ki=0, Mtis_ht=0, Mtis_br=0, Mtis_ut=0, Mtis_skel=0, Mtis_st=0, Mtis_rob=0,
           Mm_lu=0, Mm_spl=0, Mm_li=0, Mm_ki=0, Mm_ht=0, Mm_br=0, Mm_ut=0, Mm_skel=0, Mm_st=0, Mm_rob=0, 
           M_ven=dose, M_art=0, Mm_ven=0, Mm_art=0, M_feces=0, M_urine=0)


###############
# ODEs system #
###############
ode.func <- function(time, inits, params){
  with( as.list(c(inits,params)),{
    
    #In Kreyling data, soft tissue compartment consists of muscles, adipose and skin. As a result these 3 compartments will be represented as 1 (Soft tissue) here.
    W_st <- W_ad + W_skin + W_mu
    Vtis_st <- Vtis_ad + Vtis_skin + Vtis_mu
    Vcap_st <- Vcap_ad + Vcap_skin + Vcap_mu
    Q_st <- Q_ad + Q_skin + Q_mu
    Vm_st <- Vm_ad + Vm_skin + Vm_mu
    
    #Capillary concentrations
    Ccap_lu <- Mcap_lu/Vcap_lu
    Ccap_spl <- Mcap_spl/Vcap_spl
    Ccap_li <- Mcap_li/Vcap_li
    Ccap_ki <- Mcap_ki/Vcap_ki
    Ccap_ht <- Mcap_ht/Vcap_ht
    Ccap_br <- Mcap_br/Vcap_br
    Ccap_ut <- Mcap_ut/Vcap_ut
    Ccap_skel <- Mcap_skel/Vcap_skel
    Ccap_st <- Mcap_st/Vcap_st
    Ccap_rob <- Mcap_rob/Vcap_rob
    
    #Tissue concentration
    Ctis_lu <- Mtis_lu/Vtis_lu
    Ctis_spl <- Mtis_spl/Vtis_spl
    Ctis_li <- Mtis_li/Vtis_li
    Ctis_ki <- Mtis_ki/Vtis_ki
    Ctis_ht <- Mtis_ht/Vtis_ht
    Ctis_br <- Mtis_br/Vtis_br
    Ctis_ut <- Mtis_ut/Vtis_ut
    Ctis_skel <- Mtis_skel/Vtis_skel
    Ctis_st <- Mtis_st/Vtis_st
    Ctis_rob <- Mtis_rob/Vtis_rob
    
    #Concentration in Macrophages cells sub-compartment
    Cm_lu <- Mm_lu/Vm_lu
    Cm_spl <- Mm_spl/Vm_spl
    Cm_li <- Mm_li/Vm_li
    Cm_ki <- Mm_ki/Vm_ki
    Cm_ht <- Mm_ht/Vm_ht
    Cm_br <- Mm_br/Vm_br
    Cm_ut <- Mm_ut/Vm_ut
    Cm_skel <- Mm_skel/Vm_skel
    Cm_st <- Mm_st/Vm_st
    Cm_rob <- Mm_rob/Vm_rob
    
    Cm_art <- Mm_art/Vm_art
    Cm_ven <- Mm_ven/Vm_ven
    
    C_ven <- M_ven/Vven
    C_art <- M_art/Vart
    
    
    PA_lu<-x_3*Q_total
    PA_spl<-x_2*Q_spl
    PA_li<-x_1*Q_li
    PA_ki<-x_3*Q_ki
    PA_ht<-x_4*Q_ht
    PA_br<-x_6*Q_br
    PA_ut<-x_4*Q_ut
    PA_skel<-x_5*Q_skel
    PA_st<-x_5*Q_st
    PA_rob<-x_4*Q_rob
    
    P_lu <- P_3
    P_spl <- P_2
    P_li <- P_1
    P_ki <- P_3
    P_ht <- P_4
    P_br <- P_6
    P_ut <- P_4
    P_skel <- P_5
    P_st <- P_5
    P_rob <- P_4
    
    Pup_lu<-Pup*Vm_lu*(1-(Cm_lu/(Km+Cm_lu))) #ml/h
    Pup_spl<-Pup*Vm_spl*(1-(Cm_spl/(Km+Cm_spl))) #ml/h
    Pup_li<-Pup*Vm_li*(1-(Cm_li/(Km+Cm_li))) #ml/h
    Pup_ki<-Pup*Vm_ki*(1-(Cm_ki/(Km+Cm_ki))) #ml/h
    Pup_ht<-Pup*Vm_ht*(1-(Cm_ht/(Km+Cm_ht))) #ml/h
    Pup_br<-Pup*Vm_br*(1-(Cm_br/(Km+Cm_br))) #ml/h
    Pup_ut<-Pup*Vm_ut*(1-(Cm_ut/(Km+Cm_ut))) #ml/h
    Pup_skel<-Pup*Vm_skel*(1-(Cm_skel/(Km+Cm_skel))) #ml/h
    Pup_st<-Pup*Vm_st*(1-(Cm_st/(Km+Cm_st))) #ml/h
    Pup_rob<-Pup*Vm_rob*(1-(Cm_rob/(Km+Cm_rob))) #ml/h
    
    Pup_ven<-Pup*Vm_ven*(1-(Cm_ven/(Km+Cm_ven))) #ml/h
    Pup_art<-Pup*Vm_art*(1-(Cm_art/(Km+Cm_art))) #ml/h
    
    #Lungs
    dMcap_lu <- Q_total*C_ven - Q_total*Ccap_lu - PA_lu*Ccap_lu + PA_lu*Ctis_lu/P_lu
    dMtis_lu <- PA_lu*Ccap_lu - PA_lu*Ctis_lu/P_lu - Pup_lu*Ctis_lu + k_de*Mm_lu
    dMm_lu   <- Pup_lu*Ctis_lu - k_de*Mm_lu
    
    #Spleen
    dMcap_spl <- Q_spl*C_art - Q_spl*Ccap_spl - PA_spl*Ccap_spl + PA_spl*Ctis_spl/P_spl
    dMtis_spl <- PA_spl*Ccap_spl - PA_spl*Ctis_spl/P_spl - Pup_spl*Ctis_spl + k_de*Mm_spl
    dMm_spl   <- Pup_spl*Ctis_spl - k_de*Mm_spl
    
    #Liver
    dMcap_li <- Q_li*C_art + Q_spl*Ccap_spl - (Q_li+Q_spl)*Ccap_li - PA_li*Ccap_li + PA_li*Ctis_li/P_li
    dMtis_li <- PA_li*Ccap_li - PA_li*Ctis_li/P_li - Pup_li*Ctis_li + k_de*Mm_li - CLE_f*Mtis_li
    dMm_li   <- Pup_li*Ctis_li - k_de*Mm_li
    dM_feces <- CLE_f*Mtis_li
    
    #Kidneys
    dMcap_ki <- Q_ki*C_art - Q_ki*Ccap_ki - PA_ki*Ccap_ki + PA_ki*Ctis_ki/P_ki
    dMtis_ki <- PA_ki*Ccap_ki - PA_ki*Ctis_ki/P_ki - Pup_ki*Ctis_ki + k_de*Mm_ki - CLE_u*Mtis_ki
    dMm_ki   <- Pup_ki*Ctis_ki - k_de*Mm_ki
    dM_urine <- CLE_u*Mtis_ki
    
    #Heart
    dMcap_ht <- Q_ht*C_art - Q_ht*Ccap_ht - PA_ht*Ccap_ht + PA_ht*Ctis_ht/P_ht
    dMtis_ht <- PA_ht*Ccap_ht - PA_ht*Ctis_ht/P_ht - Pup_ht*Ctis_ht + k_de*Mm_ht
    dMm_ht   <- Pup_ht*Ctis_ht - k_de*Mm_ht
    
    #Brain
    dMcap_br <- Q_br*C_art - Q_br*Ccap_br - PA_br*Ccap_br + PA_br*Ctis_br/P_br
    dMtis_br <- PA_br*Ccap_br - PA_br*Ctis_br/P_br - Pup_br*Ctis_br + k_de*Mm_br
    dMm_br   <- Pup_br*Ctis_br - k_de*Mm_br
    
    #Uterus
    dMcap_ut <- Q_ut*C_art - Q_ut*Ccap_ut - PA_ut*Ccap_ut + PA_ut*Ctis_ut/P_ut
    dMtis_ut <- PA_ut*Ccap_ut - PA_ut*Ctis_ut/P_ut - Pup_ut*Ctis_ut + k_de*Mm_ut
    dMm_ut   <- Pup_ut*Ctis_ut - k_de*Mm_ut
    
    #Skeleton
    dMcap_skel <- Q_skel*C_art - Q_skel*Ccap_skel - PA_skel*Ccap_skel + PA_skel*Ctis_skel/P_skel
    dMtis_skel <- PA_skel*Ccap_skel - PA_skel*Ctis_skel/P_skel - Pup_skel*Ctis_skel + k_de*Mm_skel
    dMm_skel   <- Pup_skel*Ctis_skel - k_de*Mm_skel
    
    #Soft tissue
    dMcap_st <- Q_st*C_art - Q_st*Ccap_st - PA_st*Ccap_st + PA_st*Ctis_st/P_st
    dMtis_st <- PA_st*Ccap_st - PA_st*Ctis_st/P_st - Pup_st*Ctis_st + k_de*Mm_st
    dMm_st   <- Pup_st*Ctis_st - k_de*Mm_st
    
    #Rest of Body
    dMcap_rob <- Q_rob*C_art - Q_rob*Ccap_rob - PA_rob*Ccap_rob + PA_rob*Ctis_rob/P_rob
    dMtis_rob <- PA_rob*Ccap_rob - PA_rob*Ctis_rob/P_rob - Pup_rob*Ctis_rob + k_de*Mm_rob
    dMm_rob   <- Pup_rob*Ctis_rob - k_de*Mm_rob
    
    #Veins
    dM_ven <- - Q_total*C_ven + (Q_li+Q_spl)*Ccap_li + Q_ki*Ccap_ki + Q_ht*Ccap_ht + Q_br*Ccap_br + Q_ut*Ccap_ut + Q_skel*Ccap_skel + 
      Q_st*Ccap_st + Q_rob*Ccap_rob - Pup_ven*C_ven + k_de*Mm_ven
    dMm_ven <- Pup_ven*C_ven - k_de*Mm_ven
    
    #Arteries
    dM_art <- Q_total*Ccap_lu - Q_spl*C_art - Q_li*C_art - Q_ki*C_art - Q_ht*C_art - Q_br*C_art - Q_ut*C_art - Q_skel*C_art - Q_st*C_art - Q_rob*C_art - 
      Pup_art*C_art + k_de*Mm_art
    dMm_art <-Pup_art*C_art - k_de*Mm_art
    
    
    #Total amounts in each compartment
    #Lungs
    Lu_total <- Mcap_lu + Mtis_lu + Mm_lu
    
    #Spleen
    Spl_total <- Mcap_spl + Mtis_spl + Mm_spl
    
    #Liver
    Li_total <- Mcap_li + Mtis_li + Mm_li
    
    #Kidneys
    Ki_total <- Mcap_ki + Mtis_ki + Mm_ki
    
    #Heart
    Ht_total <- Mcap_ht + Mtis_ht + Mm_ht
    
    #Brain
    Br_total <- Mcap_br + Mtis_br + Mm_br
    
    #Uterus
    Ut_total <- Mcap_ut + Mtis_ut + Mm_ut
    
    #Skeleton
    Skel_total <- Mcap_skel + Mtis_skel + Mm_skel 
    
    #Soft tissue
    St_total <- Mcap_st + Mtis_st + Mm_st
    
    #Rest of Body
    Rob_total <- Mcap_rob + Mtis_rob + Mm_rob
    
    #Blood
    Blood_total <- M_ven + Mm_ven + M_art + Mm_art
    
    Feces_total = M_feces
    Urine_total = M_urine
    
    list(c(dMcap_lu=dMcap_lu, dMcap_spl=dMcap_spl, dMcap_li=dMcap_li, dMcap_ki=dMcap_ki, dMcap_ht=dMcap_ht, dMcap_br=dMcap_br, dMcap_ut=dMcap_ut, dMcap_skel=dMcap_skel, dMcap_st=dMcap_st, dMcap_rob=dMcap_rob,
           dMtis_lu=dMtis_lu, dMtis_spl=dMtis_spl, dMtis_li=dMtis_li, dMtis_ki=dMtis_ki, dMtis_ht=dMtis_ht, dMtis_br=dMtis_br, dMtis_ut=dMtis_ut, dMtis_skel=dMtis_skel, dMtis_st=dMtis_st, dMtis_rob=dMtis_rob,
           dMm_lu=dMm_lu, dMm_spl=dMm_spl, dMm_li=dMm_li, dMm_ki=dMm_ki, dMm_ht=dMm_ht, dMm_br=dMm_br, dMm_ut=dMm_ut, dMm_skel=dMm_skel, dMm_st=dMm_st, dMm_rob=dMm_rob,
           dM_ven=dM_ven, dM_art=dM_art, dMm_ven=dMm_ven, dMm_art=dMm_art, dM_feces=dM_feces, dM_urine=dM_urine), 
         Lu_total=Lu_total, Spl_total=Spl_total, Li_total=Li_total, Ki_total=Ki_total, Ht_total=Ht_total, Br_total=Br_total, Ut_total=Ut_total, 
         Skel_total=Skel_total, St_total=St_total, Rob_total=Rob_total, Blood_total=Blood_total, Feces_total=Feces_total, Urine_total=Urine_total)
    
    
  })}

#sample_time <- c(0, 10/60, 1, 1*24, 7*24, 28*24) #in hours
sample_time <- seq(0, 28*24, 0.1) #in hours

solution <- ode(times = sample_time, func = ode.func, y = inits, parms = params, method = "bdf")
#print(solution)
#rowSums(solution[,2:(length(inits)+1)])

Total_amounts <- solution[,38:50]


solution <- solution[,2:(length(inits)+1)]
ltime <- length(sample_time)
Capillaries_amount <- matrix(0, nrow = ltime, ncol = 11)

for (t in 1:ltime) { # t->time indice 
  ### Amount of NPs in capillaries of each organ
  
  # Amount in Lungs
  Capillaries_amount[t,1] <- solution[t,1]
  
  # Amount in Spleen
  Capillaries_amount[t,2] <- solution[t,2] 
  
  # Amount in Liver
  Capillaries_amount[t,3] <- solution[t,3] 
  
  # Amount in Kidneys
  Capillaries_amount[t,4] <- solution[t,4] 
  
  # Amount in Heart
  Capillaries_amount[t,5] <- solution[t,5] 
  
  # Amount in Brain
  Capillaries_amount[t,6] <- solution[t,6] 
  
  # Amount in Uterus
  Capillaries_amount[t,7] <- solution[t,7]  
  
  # Amount in Skeleton
  Capillaries_amount[t,8] <- solution[t,8] 
  
  # Amount in Soft tissue
  Capillaries_amount[t,9] <- solution[t,9] 
  
  # Amount in Rest of Body
  Capillaries_amount[t,10] <- solution[t,10] 
  
  # Amount in Blood
  Capillaries_amount[t,11] <- 0
}


Tissue_amount <- matrix(0, nrow = ltime, ncol = 11)
for (t in 1:ltime) { # t->time indice 
  ### Amount of NPs in tissue of each organ
  
  # Amount in Lungs
  Tissue_amount[t,1] <- solution[t,11]
  
  # Amount in Spleen
  Tissue_amount[t,2] <- solution[t,12] 
  
  # Amount in Liver
  Tissue_amount[t,3] <- solution[t,13] 
  
  # Amount in Kidneys
  Tissue_amount[t,4] <- solution[t,14] 
  
  # Amount in Heart
  Tissue_amount[t,5] <- solution[t,15] 
  
  # Amount in Brain
  Tissue_amount[t,6] <- solution[t,16] 
  
  # Amount in Uterus
  Tissue_amount[t,7] <- solution[t,17]  
  
  # Amount in Skeleton
  Tissue_amount[t,8] <- solution[t,18] 
  
  # Amount in Soft tissue
  Tissue_amount[t,9] <- solution[t,19] 
  
  # Amount in Rest of Body
  Tissue_amount[t,10] <- solution[t,20] 
  
  # Amount in Blood
  Tissue_amount[t,11] <- solution[t,31] + solution[t,32]
}

PCs_amount <- matrix(0, nrow = ltime, ncol = 11)
for (t in 1:ltime) { # t->time indice 
  ### Amount of NPs in tissue of each organ
  
  # Amount in Lungs
  PCs_amount[t,1] <- solution[t,21]
  
  # Amount in Spleen
  PCs_amount[t,2] <- solution[t,22] 
  
  # Amount in Liver
  PCs_amount[t,3] <- solution[t,23] 
  
  # Amount in Kidneys
  PCs_amount[t,4] <- solution[t,24] 
  
  # Amount in Heart
  PCs_amount[t,5] <- solution[t,25] 
  
  # Amount in Brain
  PCs_amount[t,6] <- solution[t,26] 
  
  # Amount in Uterus
  PCs_amount[t,7] <- solution[t,27]  
  
  # Amount in Skeleton
  PCs_amount[t,8] <- solution[t,28] 
  
  # Amount in Soft tissue
  PCs_amount[t,9] <- solution[t,29] 
  
  # Amount in Rest of Body
  PCs_amount[t,10] <- solution[t,30] 
  
  # Amount in Blood
  PCs_amount[t,11] <- solution[t,33] + solution[t,34]
}

###Prepare experimental data
#read data
data <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 1, colNames = T, rowNames = T)
sd <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 2, colNames = T, rowNames = T)
feces_exp <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 3, colNames = T, rowNames = F)
urine_exp <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 4, colNames = T, rowNames = F)
colnames(feces_exp) <- c("Time", "Feces") # Time in days, Feces in micro_g
colnames(urine_exp) <- c("Time", "Urine_excretion_rate", "Cumulative_ID") # Time in days
#Transform data
Transformed_data <- data
for (d in 1:length(doses)) {
  Transformed_data[d,] <- ((data[d,]/100)*doses[1])  #results give TiO2 in micro grams
  sd[d,] <- sd[d,]*doses[1]/100
}
#Drop Carcass column
Transformed_data <- subset(Transformed_data, select = -(Carcass))
sd <- subset(sd, select = -(Carcass))
feces_exp[,1] <- feces_exp[,1]*24 #transform time to hours
urine_exp[,1] <- urine_exp[,1]*24 #transform time to hours
feces_exp[,2] <- feces_exp[,2]*doses[1]/100
urine_exp <- subset(urine_exp, select = -(Urine_excretion_rate))
urine_exp[,2] <- urine_exp[,2]*doses[1]/100


# Time of sampling
time <-  c(1, 4, 24, 7*24, 28*24) #in hours
urine_time <- urine_exp[,1]*24 #in hours
feces_time <- feces_exp[,1]*24 #in hours


Liver_data <- as.data.frame(cbind(time,Transformed_data[,1], sd[,1]))
Spleen_data <- as.data.frame(cbind(time,Transformed_data[,2], sd[,2]))
Kidneys_data <- as.data.frame(cbind(time,Transformed_data[,3], sd[,3]))
Lungs_data <- as.data.frame(cbind(time,Transformed_data[,4], sd[,4]))
Heart_data <- as.data.frame(cbind(time,Transformed_data[,5], sd[,5]))
Brain_data <- as.data.frame(cbind(time,Transformed_data[,6], sd[,6]))
Uterus_data <- as.data.frame(cbind(time,Transformed_data[,7], sd[,7]))
Blood_data <- as.data.frame(cbind(time,Transformed_data[,8], sd[,8]))
Skeleton_data <- as.data.frame(cbind(time,Transformed_data[,9], sd[,9]))
Soft_tissue_data <- as.data.frame(cbind(time,Transformed_data[,10], sd[,10]))
Feces_data <- feces_exp
Urine_data <- urine_exp


colnames(Liver_data) <- c("Time", "Mass_data", "SD")
colnames(Spleen_data) <- colnames(Liver_data)
colnames(Kidneys_data) <- colnames(Liver_data)
colnames(Lungs_data) <- colnames(Liver_data)
colnames(Heart_data) <- colnames(Liver_data)
colnames(Brain_data) <- colnames(Liver_data)
colnames(Uterus_data) <- colnames(Liver_data)
colnames(Blood_data) <- colnames(Liver_data)
colnames(Skeleton_data) <- colnames(Liver_data)
colnames(Soft_tissue_data) <- colnames(Liver_data)
colnames(Feces_data) <- c("time", "Feces")
colnames(Urine_data) <- c("time", "Urine")



observed_data <- list(Lungs_data, Spleen_data, Liver_data, Kidneys_data, Heart_data, Brain_data, Uterus_data,
                      Skeleton_data, Soft_tissue_data, Blood_data, Feces_data, Urine_data)

#Prepare results for plots
Lungs <- as.data.frame(cbind(sample_time, Total_amounts[,1], Capillaries_amount[,1], Tissue_amount[,1], PCs_amount[,1]))
Spleen <- as.data.frame(cbind(sample_time, Total_amounts[,2], Capillaries_amount[,2], Tissue_amount[,2], PCs_amount[,2]))
Liver <- as.data.frame(cbind(sample_time, Total_amounts[,3], Capillaries_amount[,3], Tissue_amount[,3], PCs_amount[,3]))
Kidneys <- as.data.frame(cbind(sample_time, Total_amounts[,4], Capillaries_amount[,4], Tissue_amount[,4], PCs_amount[,4]))
Heart <- as.data.frame(cbind(sample_time, Total_amounts[,5], Capillaries_amount[,5], Tissue_amount[,5], PCs_amount[,5]))
Brain <- as.data.frame(cbind(sample_time, Total_amounts[,6], Capillaries_amount[,6], Tissue_amount[,6], PCs_amount[,6]))
Uterus <- as.data.frame(cbind(sample_time, Total_amounts[,7], Capillaries_amount[,7], Tissue_amount[,7], PCs_amount[,7]))
Skeleton <- as.data.frame(cbind(sample_time, Total_amounts[,8], Capillaries_amount[,8], Tissue_amount[,8], PCs_amount[,8]))
Soft_tissue <- as.data.frame(cbind(sample_time, Total_amounts[,9], Capillaries_amount[,9], Tissue_amount[,9], PCs_amount[,9]))
#Rest_of_body <- as.data.frame(cbind(sample_time, Total_amounts[,10], Capillaries_amount[,10], Tissue_amount[,10], PCs_amount[,10]))
Blood <- as.data.frame(cbind(sample_time, Total_amounts[,11], Capillaries_amount[,11], Tissue_amount[,11], PCs_amount[,11]))
Feces <- as.data.frame(cbind(sample_time, Total_amounts[,12]))
Urine <- as.data.frame(cbind(sample_time, Total_amounts[,13]))

colnames(Lungs) <- c("Time", "Total", "Capillaries", "Tissue", "Macrophage")
colnames(Spleen) <- colnames(Lungs)
colnames(Liver) <- colnames(Lungs)
colnames(Kidneys) <- colnames(Lungs)
colnames(Heart) <- colnames(Lungs)
colnames(Brain) <- colnames(Lungs)
colnames(Uterus) <- colnames(Lungs)
colnames(Skeleton) <- colnames(Lungs)
colnames(Soft_tissue) <- colnames(Lungs)
#colnames(Rest_of_body) <- colnames(Lungs)
colnames(Blood) <- colnames(Lungs)
colnames(Feces) <- c("Time", "Total")
colnames(Urine) <- colnames(Feces)

bag_of_data <- list(Lungs, Spleen, Liver, Kidneys, Heart, Brain, Uterus, Skeleton, Soft_tissue, Blood)
names(bag_of_data) <- c( "Lungs", "Spleen", "Liver", "Kidneys", "Heart", "Brain", "Uterus", "Skeleton", "Soft_tissue",  "Blood")
comp_names <- c( "Lungs", "Spleen", "Liver", "Kidneys", "Heart", "Brain", "Uterus", "Skeleton", "Soft_tissue",  "Blood")
counter <-1

setwd("C:\\Users\\vassi\\Documents\\Diploma Thesis\\My nano PBPK\\Final Model\\Fitting\\Test 3\\Test Plots")

for (dat in bag_of_data) {
  comp_name<-comp_names[counter]
  save_name<-paste0(comp_name, ".png", sep="")
  data_to_plot <- dat
  experimental <- observed_data[[counter]]
  
  my_plot <- ggplot(data_to_plot, aes(x=Time, y=Total, colour = "Total"))+
    geom_line(size=1.2) +
    geom_line(data = data_to_plot, aes(x=Time, y=Capillaries, colour = "Capillaries"),size=1.2)+
    geom_line(data = data_to_plot, aes(x = Time, y = Tissue, colour = "Tissue"),size=1.2)+
    geom_line(data = data_to_plot, aes(x = Time, y= Macrophage, colour = "Macrophage"),size=1.2)+
    geom_point(data = experimental, aes(x = Time, y=Mass_data,  colour = "Mass_data"),size=3)+
    #geom_errorbar(data = experimental, aes(ymin = ifelse((Mass_data - SD)>0, Mass_data - SD, 0), ymax = Mass_data + SD ), width=1)+
    
    labs(title = rlang::expr(!!comp_name),  y = "TiO2 (ug)", x = "Time (in hours)") +
    
    #labs(x = expression("Time (days)"),y = expression("Amount(micrograms)"))+
    #scale_colour_manual(name="Subcompartment", values=c("Total"=1, "Capillaries"=2, "Tissue"=3, "Macrophage"=4, "Mass_data"=5))+
    
    theme(plot.title = element_text(hjust = 0.5,size=26), axis.title.y =element_text(hjust = 0.5,size=20,face="bold"),
          axis.text.y=element_text(size=18),
          axis.title.x =element_text(hjust = 0.5,size=20,face="bold"),
          axis.text.x=element_text(size=18),
          legend.title=element_text(hjust = 0.5,size=20), 
          legend.text=element_text(size=18))
  
  png(rlang::expr(!!save_name), width = 15, height = 10, units = 'in', res = 500)
  print(my_plot)
  
  dev.off()
  counter <- counter +1
}

feces_plot <- ggplot(feces_exp, aes(Time, Feces, colour = "Mass_data"))+
  geom_point(size = 3) +
  geom_line(data = Feces, aes(x=Time, y=Total, colour="Total"))+
  
  labs(title = "Feces", y = "TiO2 (micro-grams)", x = "Time (hours)" )+
  theme(plot.title = element_text(hjust = 0.5,size=26),
        axis.title.y =element_text(hjust = 0.5,size=20,face="bold"),
        axis.text.y=element_text(size=18),
        axis.title.x =element_text(hjust = 0.5,size=20,face="bold"),
        axis.text.x=element_text(size=18),
        legend.title=element_text(hjust = 0.5,size=20), 
        legend.text=element_text(size=18))

ggsave("Feces.png", width = 15, height = 10, units = 'in')
dev.off()


urine_plot <- ggplot(urine_exp, aes(Time, Cumulative_ID, colour = "Mass_data"))+
  geom_point(size = 3) +
  geom_line(data = Urine, aes(x=Time, y=Total, colour="Total"))+
  
  
  labs(title = "Urine", y = "TiO2 (micro-grams)", x = "Time (days)" )+
  theme(plot.title = element_text(hjust = 0.5,size=26),
        axis.title.y =element_text(hjust = 0.5,size=20,face="bold"),
        axis.text.y=element_text(size=18),
        axis.title.x =element_text(hjust = 0.5,size=20,face="bold"),
        axis.text.x=element_text(size=18),
        legend.title=element_text(hjust = 0.5,size=20), 
        legend.text=element_text(size=18))

ggsave("Urine.png", width = 15, height = 10, units = 'in')
dev.off()

#print(urine_plot)



