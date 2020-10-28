library(openxlsx)
library(rstan)
library(parallel)
options(max.print=1000000)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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
  
  
  return(c("Q_total"=Q_total, "V_blood"=Total_Blood, "Vven"=Vven, "Vart"=Vart, "Vm_ven"=Vm_ven, "Vm_art"=Vm_art,
              
              "W_rob"=parameters[1,1], "W_ht"=parameters[2,1], "W_ki"=parameters[3,1], "W_br"=parameters[4,1], "W_spl"=parameters[5,1], "W_lu"=parameters[6,1], "W_li"=parameters[7,1], "W_ut"=parameters[8,1], "W_skel"=parameters[9,1],
              "W_st"=parameters[10,1] + parameters[11,1] + parameters[12,1], 
              
              "Vtis_rob"=parameters[1,2], "Vtis_ht"=parameters[2,2], "Vtis_ki"=parameters[3,2], "Vtis_br"=parameters[4,2], "Vtis_spl"=parameters[5,2], "Vtis_lu"=parameters[6,2], "Vtis_li"=parameters[7,2], "Vtis_ut"=parameters[8,2], "Vtis_skel"=parameters[9,2],
              "Vtis_st"=parameters[10,2] + parameters[11,2] + parameters[12,2],
              
              "Vcap_rob"=parameters[1,3], "Vcap_ht"=parameters[2,3], "Vcap_ki"=parameters[3,3], "Vcap_br"=parameters[4,3], "Vcap_spl"=parameters[5,3], "Vcap_lu"=parameters[6,3], "Vcap_li"=parameters[7,3], "Vcap_ut"=parameters[8,3], "Vcap_skel"=parameters[9,3],
              "Vcap_st"=parameters[10,3] + parameters[11,3] + parameters[12,3],
              
              "Vm_rob"=parameters[1,5], "Vm_ht"=parameters[2,5], "Vm_ki"=parameters[3,5], "Vm_br"=parameters[4,5], "Vm_spl"=parameters[5,5], "Vm_lu"=parameters[6,5], "Vm_li"=parameters[7,5], "Vm_ut"=parameters[8,5], "Vm_skel"=parameters[9,5],
              "Vm_st"=parameters[10,5] + parameters[11,5] + parameters[12,5],
              
              "Q_rob"=parameters[1,4]+parameters[6,4], "Q_ht"=parameters[2,4], "Q_ki"=parameters[3,4], "Q_br"=parameters[4,4], "Q_spl"=parameters[5,4], "Q_lu"=parameters[6,4], "Q_li"=parameters[7,4], "Q_ut"=parameters[8,4], "Q_skel"=parameters[9,4],
              "Q_st"=parameters[10,4] + parameters[11,4] + parameters[12,4]
              
  ))
}

params<-create.params(compartments,BW)
#print(params)

dose <- doses[1]

##############
### Priors ###
##############
P_1=6 #liver
P_2=0.55 #spleen
P_3=0.01 #kidneys lungs blood
P_4=0.1 #heart uterus
P_5=0.005 #Skeleton Soft tissue 
P_6=1e-10 #brain
x_1=1000 #liver
x_2=1000 #spleen
x_3=0.0001 #kidneys lungs blood
x_4=0.001 #heart uterus
x_5=10 #Skeleton Soft tissue
x_6=1e-10 #brain
k_de=0.0001
Km=5 #micro-g/ml 
Pup=0.995 #ml/h/(ml pcs)

#Li 2014
CLE_f=0.0008 #1/h
CLE_u=1.7 #1/h

eta_tr <- c(P_1, P_2, P_3, P_4, P_5, P_6, x_1, x_2, x_3, x_4, x_5, x_6,  k_de, Km, Pup, CLE_f, CLE_u)

#sampling time
time <- c(10/60, 1, 1*24, 7*24, 28*24) #in hours



#######################
###Experimental data###
#######################
data <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 1, colNames = T, rowNames = T)
sd <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 2, colNames = T, rowNames = T)
feces <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 3, colNames = T, rowNames = F)
urine <- openxlsx::read.xlsx("Kreyling IV data.xlsx", sheet = 4, colNames = T, rowNames = F)
colnames(feces) <- c("Time", "Feces") # Time in days, Feces in micro_g
colnames(urine) <- c("Time", "Urine_excretion_rate", "Cumulative_ID") # Time in days
Transformed_data <- data
for (d in 1:length(doses)) {
  Transformed_data[d,] <- ((data[d,]/100)*doses[1])  #results give TiO2 in micro grams
  sd[d,] <- sd[d,]*doses[1]/100
}
#Drop Carcass column
Transformed_data <- subset(Transformed_data, select = -(Carcass))
sd <- subset(sd, select = -(Carcass))
feces[,2] <- feces[,2]*doses[1]/100
feces_time <- feces[,1]*24 #transform time to hours
feces <- feces[,2]
urine <- subset(urine, select = -(Urine_excretion_rate))
urine_time <- urine[,1]*24 #transform time to hours
urine <- urine[,2]*doses[1]/100


#############################################################################################

DataList <- list(
  eta_tr = eta_tr, 
  N_diff = 36,
  N_compart= 10,
  N_param = length(eta_tr),
  params = params,
  time = time,
  feces_time = feces_time,
  urine_time = urine_time,
  mass = Transformed_data, #MassData,
  feces = feces,
  urine = urine,
  t_init = 0,
  m0 = c(rep(0,30),doses[1], rep(0,5)), #micrograms
  rel_tol = 1e-06, 
  abs_tol = 1e-06,
  max_num_steps = 1e05)

tic = proc.time()

fit <- stan(file = 'TiO2 PBPK Test 3 Fitting.stan', data = DataList, iter = 1000, warmup=400, chains=4)
options(max.print=5.5E5) 

#check_hmc_diagnostics(fit)
#check_divergences(fit)
#check_treedepth(fit)
#check_energy(fit)
#if I had set  control=list(max_treedepth=15) the the correct command would be check_treedepth(fit, 15) 

#library(shinystan)
#launch_shinystan(fit)
#print(fit)


#clock<-proc.time() - tic
#print(clock)
#pairs(fit, pars=c("par"))
#traceplot(fit,c("epsilon"));
#stan_dens(fit,("par[1]"),separate_chains = TRUE)
#exp_fit <- extract(fit)
#mean(exp_fit$theta[,1,2])
