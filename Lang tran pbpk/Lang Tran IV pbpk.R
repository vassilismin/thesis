library(rstan)
library(ggplot2)
library(parallel)

setwd("C:\\Users\\vassi\\Documents\\Diploma Thesis\\Lang Tran Files\\R  Stan codes")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####################################################################
###################             DATA             ###################
####################################################################
dose <- 1.7e03 #micro grams


# tissue volumes (in ml, assuming density = 1 g/ml according to Brown et al., 1997 page 26) 
V_tot <-250 # ;body weight, experimental data - ml
V_li <- 0.037*V_tot # volume of liver, literature (Brown et al., 1997) - ml 
V_gi <- 0.0248*V_tot # volume of gastrointestinal track, literature (Brown et al., 1997) - ml 
V_ki <- 0.0073*V_tot # volume of kidneys, literature (Brown et al., 1997) - ml
V_ht <- 0.0033*V_tot # volume of heart, literature (Brown et al., 1997) - ml
V_spl <- 0.0020*V_tot # volume of spleen, literature (Brown et al., 1997) - ml
V_br <- 0.0060*V_tot # volume of brain, literature (Brown et al., 1997) - ml
V_re <- V_tot - V_br - V_ht - V_ki - V_li - V_spl - V_gi # volume of rest of the body, experimental data - ml

#blood volumes
Vven <- 5.6  #psaximo
Vart <- 11.3 #psaximo

# volume of capillary blood assuming density = 1
# (capillary blood is given as a percentage of tissue volume)
Vb_li <-  0.21 * V_li # volume of blood in liver, literature (Brown et al., 1997) - ml
Vb_gi <-  0.0265 * V_br  # volume of blood in brain, literature (Sweeney et al ( mice)) - ml
Vb_ki <- 0.16 * V_ki # volume of blood in kidneys, literature (Brown et al., 1997) - ml
Vb_ht <- 0.26 * V_ht # volume of blood in heart, literature (Brown et al., 1997) - ml
Vb_spl <- 0.22 * V_spl # volume of blood in spleen, literature (Brown et al., 1997) - ml
Vb_br <- 0.03 * V_br # volume of blood in brain, literature (Brown et al., 1997) - ml
Vb_re <- 0.04 * V_re # volume of blood in rest of the body, literature (Brown et al., 1997) - ml 
# volume of arterial and venous blood, literature (Law et al., 2017) - ml
V_blood <- 0.0816 * V_tot - (Vb_spl + Vb_li + Vb_gi + Vb_ht + Vb_ki + Vb_br + Vb_re) 

#Regional blood flows (in mL per day)
fQs = 0.0086 # fraction of cardiac output to spleen, literature (Sweeney et al., 2014) - unitless
fQli = 0.174 # fraction of cardiac output to liver, literature (Brown et al., 1997) - unitless
fQgi = 0.14	# fraction of cardiac output to gastrointestinal track, literature (Brown et al., 1997) - unitless
fQh = 0.0448 # fraction of cardiac output to heart, literature (Brown et al., 1997) - unitless
fQk = 0.0948 # fraction of cardiac output to kidneys, literature (Brown et al., 1997) - unitless
fQbr = 0.02	# fraction of cardiac output to brain, literature (Brown et al., 1997) - unitless
fQrest = 1-fQs-fQli-fQgi-fQh-fQk-fQbr # fraction of cardiac output to rest of the body, fitted value - unitless

Q_tot <- 58705*(V_tot/1000)^0.75 # cardiac output, literature (Brown et al., 1997) - mL/day, 
Q_li <- fQli*Q_tot	# blood flow to liver - mL/day
Q_gi <- fQgi*Q_tot # blood flow to brain - mL/day		
Q_ki <- fQk*Q_tot # blood flow to kidneys - mL/day
Q_ht <- fQh*Q_tot	# blood flow to heart - mL/day
Q_spl <- fQs*Q_tot # blood flow to spleen - mL/day		
Q_br <- fQbr*Q_tot # blood flow to brain - mL/day		
Q_re <- fQrest*Q_tot # blood flow to rest of the body - mL/day

Da <- 0
Du <- 0
Do <- 0
kB <- 0

ko <- 0
ku <- 0
kt <- 0
kr <- 0
kd <- 0
ki <- 0


fo <- 0
fu <- 0
ft <- 0

lamdav <- 1  #psaximo

params <- c(V_gi, V_br, V_ht, V_ki, V_li, V_spl, V_re, Vven, Vart, Vb_gi, Vb_br,Vb_ht,
            Vb_ki, Vb_li, Vb_spl, Vb_re, V_blood, Q_tot, Q_gi, Q_br, Q_ht, Q_ki, Q_spl, Q_li, Q_re, Da, Du, Do, kB, ko, ku, kt,
            kr, kd, ki, fo, fu, ft, lamdav)

#####################
# Priors #
#####################
lamda31 <- 3.8 # partition coefficient liver tissue:capillary blood, fitted value - unitless
lamda41 <- 3.8 # partition coefficient GI tissue:capillary blood, fitted value - unitless
lamda51 <- 3.8 # partition coefficient Kidneys tissue:capillary blood, fitted value - unitless
lamda61 <- 3.8 # partition coefficient Heart tissue:capillaryblood, fitted value - unitless
lamda71 <- 3.8 # partition coefficient Spleen tissue:capillary blood, fitted value - unitless
lamda81 <- 3.8 # partition coefficient Brain tissue:capillary blood, fitted value - unitless
lamda91 <- 3.8 # partition coefficient Other tissues:capillary blood, fitted value - unitless

lamda32 <- 111 # permeability coefficient from blood to liver, fitted value - unitless
lamda42 <- 111 # permeability coefficient from blood to GI, fitted value - unitless
lamda52 <- 0.2 # permeability coefficient from blood to kidneys, fitted value - unitless
lamda62 <- 0.2 # permeability coefficient from blood to heart, fitted value - unitless
lamda72 <- 111 # permeability coefficient from blood to spleen, fitted value - unitless
#lamda82 <- 0   # permeability coefficient from blood to brain, fitted value - unitless
lamda82 <- 1e-20  # permeability coefficient from blood to brain, fitted value - unitless
lamda92 <- 0.2 # permeability coefficient from blood to other organs, fitted value - unitless

lamda33 <- 1968 #Liver sequestration rate, fitted value - 1/day
lamda43 <- 1968 #GI sequestration rate, fitted value - 1/day
lamda53 <- 1968 #Kidney sequestration rate, fitted value - 1/day
lamda63 <- 1968 #Heart sequestration rate, fitted value - 1/day
lamda73 <- 1368 #Spleen sequestration rate, fitted value - 1/day
lamda83 <- 1968 #Brain sequestration rate, fitted value - 1/day
lamda93 <- 1968 #Other organs sequestration rate, fitted value - 1/day

lamda35 <- 50    #bile:liver partition coefficient
lamda44 <- 1e-03 #fecal elimination from GI - 1/day
#lamda54 <- 0     #urinary elimination from kidneys - 1/day
lamda54 <- 1e-20     #urinary elimination from kidneys - 1/day

lamdaI <- 1 #psaximo
kb <- 500   #psaximo
kl <- 0.3   #psaximo

eta_tr <- c(lamda31, lamda32, lamda33, lamda41, lamda42, lamda43, lamda51, lamda52, lamda53, lamda61, lamda62, lamda63, lamda71,
            lamda72, lamda73, lamda81,lamda82, lamda83, lamda91, lamda92, lamda93, lamda35, lamda44, lamda54, lamdaI, kb, kl)


# Time of sampling
time <-  c(10/60, 1, 24, 24*7, 24*28, 24*56)
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
data <- openxlsx::read.xlsx("elgrabli_data.xlsx", sheet = 1, colNames = T, rowNames = T)
DataList <- list(
  eta_tr = eta_tr, #cl, a (kp coeff), kp_fat,kp_rest
  N_diff = 30,
  N_compart= 6,
  N_param = length(eta_tr),
  params = params,
  time = time,
  mass = TransformedData, #MassData,
  samp = Nsamp, 
  N_obs = length(MassData),
  t_init = 0,
  m0 = c(rep(0,29),dose),
  rel_tol = 1e-06, 
  abs_tol = 1e-06,
  max_num_steps = 1e05
)

tic = proc.time()

fit <- stan(file = 'lang tran stan model.stan', data = DataList, iter = 1000, warmup=400, chains=4)
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