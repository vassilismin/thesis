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
  
  Q_total <- 1.54*w^0.75 # Total Cardiac Output (ml/min)
  
  Total_Blood <- 0.06*w+0.77 # Total blood volume (ml)
  
  #Arterial blood volume
  Vart <- 1.2905*w/100 #0.15*Total_Blood #(ml)
  
  #Veins blood volume
  Vven <- 2.968*w/100 #0.64*Total_Blood #(ml)
  
  
  fr_ad <- 0.0199*w + 1.644 # w in g,  Brown et al.1997 p.420. This equation gives the  adipose % of body weight 
  
  #Tissue weight fraction 
  Tissue_fractions <- c(NA, 0.33 , 0.73, 0.57, 0.20, 0.5, 3.66, 0.011, 7, fr_ad, 19.03)/100 # % of BW. Na values refers to the volume of the rest organs(RoB)
  
  #Regional blood flow fraction
  Regional_flow_fractions <- c(NA, 4.9, 14.1, 2.0, 1.22, 2.1, 17.4, 1.11, 12.2, 7.0, 5.8)/100 # % of total cardiac output
  
  #Capillary volume fractions (fractions of tissue volume)
  Capillary_fractions <- c(NA, 0.26, 0.16, 0.03, 0.22, 0.36, 0.21, 0.077, 0.04, 0.0055, 0.02) # of tissue volume
                                                                                       # Where NA, it is the same value as Rest of Body due to luck od data for uterus and adipose
  
  W_tis <- rep(0,length(comp_names))
  V_tis <- rep(0,length(comp_names))
  V_cap <- rep(0,length(comp_names))
  Q <- rep(0,length(comp_names))
  
  
  for (i in 1:length(comp_names)) {
    control <- comp_names[i]
    
    Tissue_fractions[i] <- ifelse(is.na(control), NA, Tissue_fractions[i])
    Regional_flow_fractions[i] <- ifelse(is.na(control), NA, Regional_flow_fractions[i])
    Capillary_fractions[i] <- ifelse(is.na(control), NA, Capillary_fractions[i])
      
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
    
    ###Calculation of regional blood flows
    Q[i] <- Q_total*Regional_flow_fractions[i]
  }

  ### Calculations for "Rest of Body" compartment
  W_tis[1] <- w - sum(W_tis[2:length(W_tis)], na.rm = TRUE)
  V_tis[1] <- W_tis[1]     #(considering that the density of the rest tissues is 1 g/ml)
  Q[1] <- Q_total - sum(Q[2:length(Q)],na.rm = TRUE)
  
  a<- -(- Vven - Vart - sum(V_cap[2:length(V_cap)], na.rm = TRUE))
  V_cap[1] <- Total_Blood - Vven - Vart - sum(V_cap[2:length(V_cap)], na.rm = TRUE)
  #Capillary_fractions[1] <- V_cap[1]/V_tis[1]
  

  parms <- matrix(c(W_tis[],V_tis[],V_cap[],Q[]), ncol = 4)
  colnames(parms) <- c("W_tis", "V_tis", "V_cap", "Q")
  rownames(parms) <- all_comps
  
  
  return(list(parms, "a"=a,
              "Q_total"=Q_total, "V_blood"=Total_Blood, "Vven"=Vven, "Vart"=Vart,
              
              "W_rob"=parms[1,1], "W_ht"=parms[2,1], "W_ki"=parms[3,1], "W_br"=parms[4,1], "W_spl"=parms[5,1], "W_lu"=parms[6,1], "W_li"=parms[7,1], "W_ut"=parms[8,1], "W_skel"=parms[9,1],
              "W_ad"=parms[10,1], "W_skin"=parms[11,1],  
              
              "V_rob"=parms[1,2], "V_ht"=parms[2,2], "V_ki"=parms[3,2], "V_br"=parms[4,2], "V_spl"=parms[5,2], "V_lu"=parms[6,2], "V_li"=parms[7,2], "V_ut"=parms[8,2], "V_skel"=parms[9,2],
              "V_ad"=parms[10,2], "W_skin"=parms[11,2],
              
              "Vcap_rob"=parms[1,3], "Vcap_ht"=parms[2,3], "Vcap_ki"=parms[3,3], "Vcap_br"=parms[4,3], "Vcap_spl"=parms[5,3], "Vcap_lu"=parms[6,3], "Vcap_li"=parms[7,3], "Vcap_ut"=parms[8,3], "Vcap_skel"=parms[9,3],
              "Vcap_ad"=parms[10,3], "Vcap_skin"=parms[11,3], 
              
              "Q_rob"=parms[1,4], "Q_ht"=parms[2,4], "Q_ki"=parms[3,4], "Q_br"=parms[4,4], "Q_spl"=parms[5,4], "Q_lu"=parms[6,4], "Q_li"=parms[7,4], "Q_ut"=parms[8,4], "Q_skel"=parms[9,4],
              "Q_ad"=parms[10,4], "Q_skin"=parms[11,4]))
}

params <- create.params(compartments,BW)
print(params)