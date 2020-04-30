### Iportant!!! each compartment has a specific index vectors Tissue_rates, Regional_flow_rates, Capillary_rates and cannot be changed
# The index of each compartment:
#Rest of Body (rob) --> 1
#Heart (ht) --> 2
#Kidneys (ki) --> 3
#Brain (br) --> 4
#Spleen (spl) --> 5
#Lungs (lu) --> 6
#Liver (li) --> 7

####################
### User's INPUT ###
####################
#### If any of these compartments don not exist in pbpk, just give it the value NA in compartments vector, example: "Heart" = NA and it will remove it 
#### from the equilibriums and the corresponding V_tis, V_cap, Q will be equal to 0.
  

compartments <- list( "RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"="Brain", "Spleen"="Spleen",
                     "Lungs"="Lungs", "Liver"="Liver") #used as input in function, compartments that are used in pbpk
BW <- 250 # Total Body weight of rat in g

#####################################
### Function to create Parameters ###
#####################################

create.params <- function(comp_names, w){
  
  all_comps <- list("RoB"="RoB", "Heart"="Heart", "Kidneys"="Kidneys", "Brain"="Brain", "Spleen"="Spleen",
                    "Lungs"="Lungs", "Liver"="Liver") # List with names of all possible compartments
  
  Q_total <- 1.54*BW^0.75 # Total Cardiac Output (ml/min)
  
  Total_Blood <- 0.06*BW+0.77 # Total blood volume (ml)
  
  #Arterial blood volume
  Vart <- 0.15*Total_Blood #(ml)
  
  #Veins blood volume
  Vven <- 0.64*Total_Blood #(ml)
  
  
  #Tissue weight or volumes rates 
  Tissue_rates <- c(NA, 0.33 , 0.73, 0.60, 0.20, 0.5, 3.66)/100 # % of BW. Na values refers to the volume of the rest organs(RoB)
  
  #Regional blood flow rates
  Regional_flow_rates <- c(NA, 4.9, 14.1, 2.0, 1.22, 2.1, 17.4)/100 # % of total cardiac output
  
  #Capillary volume rates (rates of tissue volume)
  Capillary_rates <- c(NA, 0.26, 0.16, 0.03, 0.22, 0.36, 0.21)/100 # % of tissue volume
  
  #Tissue_rates[1] <- 100 - sum(Tissue_rates[2:length(Tissue_rates)])
  #Regional_flow_rates[1] <- 100 - sum(Regional_flow_rates[2:length(Regional_flow_rates)])
  #Capillary_rates[1] <- 100 - sum(Capillary_rates[2:length(Capillary_rates)]) - ((Vven + Vart)/Total_Blood)*100
  
  #Arterial blood volume
  Vart <- 0.15*Total_Blood
  
  #Veins blood volume
  Vven <- 0.64*Total_Blood
  
  V_tis <- rep(0,length(comp_names))
  V_cap <- rep(0,length(comp_names))
  Q <- rep(0,length(comp_names))
  
  
  for (i in 1:length(comp_names)) {
    control <- comp_names[i]
    
    Tissue_rates[i] <- ifelse(is.na(control), 0, Tissue_rates[i])
    Regional_flow_rates[i] <- ifelse(is.na(control), 0, Regional_flow_rates[i])
    Capillary_rates[i] <- ifelse(is.na(control), 0, Capillary_rates[i])
      
      
    
    V_tis[i] <- w*Tissue_rates[i]
    V_cap[i] <- V_tis[i]*Capillary_rates[i]
    Q[i] <- Q_total*Regional_flow_rates[i]
    }

  ### Calculations for "Rest of Body" compartment
  V_tis[1] <- w - sum(V_tis[2:length(V_tis)])
  Q[1] <- Q_total - sum(Q[2:length(Q)])
  V_cap[1] <- Total_Blood - Vven - Vart - sum(V_cap[2:length(V_cap)])
  
  parms <- matrix(c(V_tis[],V_cap[],Q[]), ncol = 3)
  colnames(parms) <- c("V_tis", "V_cap", "Q")
  rownames(parms) <- all_comps
  
  
  return(list(parms,
              "Q_total"=Q_total, "V_blood"=Total_Blood,
              "V_rob"=parms[1,1], "V_ht"=parms[2,1], "V_ki"=parms[3,1], "V_br"=parms[4,1], "V_spl"=parms[5,1], "V_lu"=parms[6,1], "V_li"=parms[7,1], "Vven"=Vven, "Vart"=Vart,
              "Vcap_rob"=parms[1,2],  "Vcap_ht"=parms[2,2], "Vcap_ki"=parms[3,2], "Vcap_br"=parms[4,2], "Vcap_spl"=parms[5,2], "Vcap_lu"=parms[6,2], "Vcap_li"=parms[7,2],
              "Q_rob"=parms[1,3], "Q_ht"=parms[2,3], "Q_ki"=parms[3,3], "Q_br"=parms[4,3], "Q_spl"=parms[5,3], "Q_lu"=parms[6,3], "Q_li"=parms[7,3]))
}

params <- create.params(compartments,BW)
print(params)