library(openxlsx)
library(deSolve)
setwd("C:\\Users\\vassi\\Documents\\Diploma Thesis\\Parameters Function") #excel with fractions
dose=20 #micro g

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
  Q[1] <- Q_total - sum(Q[2:length(Q)],na.rm = TRUE)
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
              
              #Extra parameters (to be fitted)
              #Carlander
              P=3.8,
              x_fast=111,
              x_rest=0.2,
              x_brain=21.1,
              k_de=0.001,
              #Aborig
              Km=50, #micro-g/ml (not fitted)
              Pup=0.05, #ml/h/(ml pcs)
              
              #Li 2014
              CLE_f=0.000118, #1/h
              CLE_u=0.000001 #1/h
              
              
  ))
}

params<-create.params(compartments,BW)
print(params)

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
    
    
    PA_lu<-x_rest*Q_total
    PA_spl<-x_fast*Q_spl
    PA_li<-x_fast*Q_li
    PA_ki<-x_rest*Q_ki
    PA_ht<-x_rest*Q_ht
    PA_br<-x_brain*Q_br
    PA_ut<-x_rest*Q_ut
    PA_skel<-x_rest*Q_skel
    PA_st<-x_rest*Q_st
    PA_rob<-x_rest*Q_rob
    
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
    dMcap_lu <- Q_total*C_ven - Q_total*Ccap_lu - PA_lu*Ccap_lu + PA_lu*Ctis_lu/P
    dMtis_lu <- PA_lu*Ccap_lu - PA_lu*Ctis_lu/P - Pup_lu*Ctis_lu + k_de*Mm_lu
    dMm_lu   <- Pup_lu*Ctis_lu - k_de*Mm_lu
    
    #Spleen
    dMcap_spl <- Q_spl*C_art - Q_spl*Ccap_spl - PA_spl*Ccap_spl + PA_spl*Ctis_spl/P
    dMtis_spl <- PA_spl*Ccap_spl - PA_spl*Ctis_spl/P - Pup_spl*Ctis_spl + k_de*Mm_spl
    dMm_spl   <- Pup_spl*Ctis_spl - k_de*Mm_spl
    
    #Liver
    dMcap_li <- Q_li*C_art + Q_spl*Ccap_spl - (Q_li+Q_spl)*Ccap_li - PA_li*Ccap_li + PA_li*Ctis_li/P
    dMtis_li <- PA_li*Ccap_li - PA_li*Ctis_li/P - Pup_li*Ctis_li + k_de*Mm_li - CLE_f*Mtis_li
    dMm_li   <- Pup_li*Ctis_li - k_de*Mm_li
    dM_feces <- CLE_f*Mtis_li
    
    #Kidneys
    dMcap_ki <- Q_ki*C_art - Q_ki*Ccap_ki - PA_ki*Ccap_ki + PA_ki*Ctis_ki/P
    dMtis_ki <- PA_ki*Ccap_ki - PA_ki*Ctis_ki/P - Pup_ki*Ctis_ki + k_de*Mm_ki - CLE_u*Mtis_ki
    dMm_ki   <- Pup_ki*Ctis_ki - k_de*Mm_ki
    dM_urine <- CLE_u*Mtis_ki
    
    #Heart
    dMcap_ht <- Q_ht*C_art - Q_ht*Ccap_ht - PA_ht*Ccap_ht + PA_ht*Ctis_ht/P
    dMtis_ht <- PA_ht*Ccap_ht - PA_ht*Ctis_ht/P - Pup_ht*Ctis_ht + k_de*Mm_ht
    dMm_ht   <- Pup_ht*Ctis_ht - k_de*Mm_ht
    
    #Brain
    dMcap_br <- Q_br*C_art - Q_br*Ccap_br - PA_br*Ccap_br + PA_br*Ctis_br/P
    dMtis_br <- PA_br*Ccap_br - PA_br*Ctis_br/P - Pup_br*Ctis_br + k_de*Mm_br
    dMm_br   <- Pup_br*Ctis_br - k_de*Mm_br
    
    #Uterus
    dMcap_ut <- Q_ut*C_art - Q_ut*Ccap_ut - PA_ut*Ccap_ut + PA_ut*Ctis_ut/P
    dMtis_ut <- PA_ut*Ccap_ut - PA_ut*Ctis_ut/P - Pup_ut*Ctis_ut + k_de*Mm_ut
    dMm_ut   <- Pup_ut*Ctis_ut - k_de*Mm_ut
    
    #Skeleton
    dMcap_skel <- Q_skel*C_art - Q_skel*Ccap_skel - PA_skel*Ccap_skel + PA_skel*Ctis_skel/P
    dMtis_skel <- PA_skel*Ccap_skel - PA_skel*Ctis_skel/P - Pup_skel*Ctis_skel + k_de*Mm_skel
    dMm_skel   <- Pup_skel*Ctis_skel - k_de*Mm_skel
    
    #Soft tissue
    dMcap_st <- Q_st*C_art - Q_st*Ccap_st - PA_st*Ccap_st + PA_st*Ctis_st/P
    dMtis_st <- PA_st*Ccap_st - PA_st*Ctis_st/P - Pup_st*Ctis_st + k_de*Mm_st
    dMm_st   <- Pup_st*Ctis_st - k_de*Mm_st
    
    #Rest of Body
    dMcap_rob <- Q_rob*C_art - Q_rob*Ccap_rob - PA_rob*Ccap_rob + PA_rob*Ctis_rob/P
    dMtis_rob <- PA_rob*Ccap_rob - PA_rob*Ctis_rob/P - Pup_rob*Ctis_rob + k_de*Mm_rob
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
    
    list(c(dMcap_lu=dMcap_lu, dMcap_spl=dMcap_spl, dMcap_li=dMcap_li, dMcap_ki=dMcap_ki, dMcap_ht=dMcap_ht, dMcap_br=dMcap_br, dMcap_ut=dMcap_ut, dMcap_skel=dMcap_skel, dMcap_st=dMcap_st, dMcap_rob=dMcap_rob,
           dMtis_lu=dMtis_lu, dMtis_spl=dMtis_spl, dMtis_li=dMtis_li, dMtis_ki=dMtis_ki, dMtis_ht=dMtis_ht, dMtis_br=dMtis_br, dMtis_ut=dMtis_ut, dMtis_skel=dMtis_skel, dMtis_st=dMtis_st, dMtis_rob=dMtis_rob,
           dMm_lu=dMm_lu, dMm_spl=dMm_spl, dMm_li=dMm_li, dMm_ki=dMm_ki, dMm_ht=dMm_ht, dMm_br=dMm_br, dMm_ut=dMm_ut, dMm_skel=dMm_skel, dMm_st=dMm_st, dMm_rob=dMm_rob,
           dM_ven=dM_ven, dM_art=dM_art, dMm_ven=dMm_ven, dMm_art=dMm_art, dM_feces=dM_feces, dM_urine=dM_urine), 
         Lu_total=Lu_total, Spl_total=Spl_total, Li_total=Li_total, Ki_total=Ki_total, Ht_total=Ht_total, Br_total=Br_total, Ut_total=Ut_total, 
         Skel_total=Skel_total, St_total=St_total, Rob_total=Rob_total, Blood_total=Blood_total)
    
    
  })}

sample_time <- c(0, 10/60, 1, 1*24, 7*24, 28*24, 56*24) #in hours
solution <- ode(times = sample_time, func = ode.func, y = inits, parms = params, method = "bdf")
print(solution)
rowSums(solution[,2:(length(inits)+1)])


