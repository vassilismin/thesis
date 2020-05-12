library(openxlsx)

### Iportant!!! each compartment has a specific index vectors Tissue_fractions, Regional_flow_fractions, Capillary_fractions and cannot be changed
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

####################
### User's INPUT ###
####################
#### If any of these compartments don not exist in pbpk, just give it the value NA in compartments vector, example: "Heart" = NA and it will remove it 
#### from the equilibriums and the corresponding V_tis, V_cap, Q will be equal to NA.


compartments <- list( "RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"="Brain", "Spleen"="Spleen",
                      "Lungs"="Lungs", "Liver"="Liver", "Uterus"="Uterus", "Skeleton"="Skeleton", "Adipose"="Adipose", "Skin"="Skin") #used as input in function, compartments that are used in pbpk
BW <- 263 # Total Body weight of rat in g

#####################################
### Function to create Parameters ###
#####################################

create.params <- function(comp_names, w){
  
  # List with names of all possible compartments
  all_comps <- list("RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"="Brain", "Spleen"="Spleen",
                    "Lungs"="Lungs", "Liver"="Liver", "Uterus"="Uterus", "Skeleton"="Skeleton", "Adipose"="Adipose", "Skin"="Skin") # List with names of all possible compartments
  
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
  
  setwd("C:\\Users\\vassi\\Documents\\Diploma Thesis\\Parameters Function")
  #read data
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
  
  
  parms <- matrix(c(W_tis[],V_tis[],V_cap[],Q[],V_macro[]), ncol = 5)
  colnames(parms) <- c("W_tis", "V_tis", "V_cap", "Q", "V_macro")
  rownames(parms) <- all_comps
  
  
  return(list(parms,
              "Q_total"=Q_total, "V_blood"=Total_Blood, "Vven"=Vven, "Vart"=Vart, "Vm_ven"=Vm_ven, "Vm_art"=Vm_art,
              
              "W_rob"=parms[1,1], "W_ht"=parms[2,1], "W_ki"=parms[3,1], "W_br"=parms[4,1], "W_spl"=parms[5,1], "W_lu"=parms[6,1], "W_li"=parms[7,1], "W_ut"=parms[8,1], "W_skel"=parms[9,1],
              "W_ad"=parms[10,1], "W_skin"=parms[11,1],  
              
              "Vtis_rob"=parms[1,2], "Vtis_ht"=parms[2,2], "Vtis_ki"=parms[3,2], "Vtis_br"=parms[4,2], "Vtis_spl"=parms[5,2], "Vtis_lu"=parms[6,2], "Vtis_li"=parms[7,2], "Vtis_ut"=parms[8,2], "Vtis_skel"=parms[9,2],
              "Vtis_ad"=parms[10,2], "Vtis_skin"=parms[11,2],
              
              "Vcap_rob"=parms[1,3], "Vcap_ht"=parms[2,3], "Vcap_ki"=parms[3,3], "Vcap_br"=parms[4,3], "Vcap_spl"=parms[5,3], "Vcap_lu"=parms[6,3], "Vcap_li"=parms[7,3], "Vcap_ut"=parms[8,3], "Vcap_skel"=parms[9,3],
              "Vcap_ad"=parms[10,3], "Vcap_skin"=parms[11,3], 
              
              "Vm_rob"=parms[1,5], "Vm_ht"=parms[2,5], "Vm_ki"=parms[3,5], "Vm_br"=parms[4,5], "Vm_spl"=parms[5,5], "Vm_lu"=parms[6,5], "Vm_li"=parms[7,5], "Vm_ut"=parms[8,5], "Vm_skel"=parms[9,5],
              "Vm_ad"=parms[10,5], "Vm_skin"=parms[11,5],
              "Q_rob"=parms[1,4], "Q_ht"=parms[2,4], "Q_ki"=parms[3,4], "Q_br"=parms[4,4], "Q_spl"=parms[5,4], "Q_lu"=parms[6,4], "Q_li"=parms[7,4], "Q_ut"=parms[8,4], "Q_skel"=parms[9,4],
              "Q_ad"=parms[10,4], "Q_skin"=parms[11,4]
              
  ))
}

params<-create.params(compartments,BW)
print(params)