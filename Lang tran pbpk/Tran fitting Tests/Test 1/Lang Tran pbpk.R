library(rstan)
library(ggplot2)
library(parallel)

setwd("C:\\Users\\vassi\\Documents\\Diploma Thesis\\Lang Tran Files\\R  Stan codes\\Tran fitting Tests\\Test 1")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####################################################################
###################             DATA             ###################
####################################################################

dose <- 1.7e03 #micro grams

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

params <- c(Vtis[],Vcap[], Q[], QTotal, Da, Do, Du, Vven, Vart, Q_liver, Q_bile)


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

eta_tr <- c(lambda3[], lambda4[], lambda5[], lambda6[], lambda7[], lambda8[], lambda9[], lambdaI,  k_B, kb, ko, ku, kt, kr, kd, ki, kl)

#Time of sampling
time <- c(10/(60*24), 1/24, 1, 7,28, 56) #in days

#read data
data <- openxlsx::read.xlsx("elgrabli_data.xlsx", sheet = 1, colNames = T, rowNames = T)
# Discard control, urine (for the time) and the first row which is the compartment names 
# and also transpose the data so that each column represents a compartment
TransformedData <- t(as.matrix(sapply(data, as.numeric))[1:8,1:7])
# Retrieve the compartment names
colnames(TransformedData) <- rownames(data[1:8,])
# Subtract the control values
for (i in 2:7){
  for (j in 1:8){
    TransformedData[i,j] <- TransformedData[i,j] -TransformedData[1,j]
    if( TransformedData[i,j]<=0){
      TransformedData[i,j] <- TransformedData[i,j]+1e-25
    }
  }
}
#Discard the control measurments
TransformedData <- TransformedData[2:7,]
# Discard testis and lymph nodes for the first test
TransformedData <- TransformedData[,1:6]
#create an array to store the number of samples of each patient
Nsamp <- apply(TransformedData, 2, function(x) length(na.omit(x)))

# Create a vectorized form of the data and a corresponding time vector
time_vec <- list()
MassData <- list()
for (i in 1:dim(TransformedData)[2]){
  if (Nsamp[i] != length(time)){
    time_vec[[i]] <- time[-which(is.na(TransformedData[,i]))]
    MassData[[i]] <- TransformedData[-which(is.na(TransformedData[,i])),i]
  } else {
    time_vec[[i]] <- time
    MassData[[i]] <- TransformedData[,i]
  }
}
# Unlist the lists and convert them to vectors
time_vec <- as.vector(unlist(time_vec))
MassData <- as.vector(unlist(MassData))

#############################################################################################

DataList <- list( eta_tr = eta_tr,
                  N_diff = 32,
                  N_compart = 6,
                  N_param = length(eta_tr),
                  params = params,
                  time = time,
                  mass = TransformedData,
                  samp = Nsamp,
                  N_obs = length(MassData),
                  t_init = 0,
                  m0 = c(rep(0,29), dose, 0, 0),
                  rel_tol = 1e-06, 
                  abs_tol = 1e-06,
                  max_num_steps = 1e05)

tic = proc.time()

fit <- stan(file = 'Lang Tran pbpk.stan', data = DataList, iter = 1000, warmup=400, chains=4)
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