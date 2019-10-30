library(rstan)
library(ggplot2)
library(parallel)

#setwd("C:/Users/user/Desktop/paper_revision/numerical/multivariate")
#setwd("/home/periklis/Desktop/test")


####################################################################
###################             DATA             ###################
####################################################################

dose <- 1.7e03 #micro grams

# tissue weights (in g)
W_tot <-253.4 # ;body weight, experimental data - g
W_lu <-0.005*W_tot # weight of lungs, experimental data - g
W_bm <- 0.0559*W_tot  # weight of bone marrow, literature (Travlos 2006) - g
#   W_bm_exp <- 0.02 # weight of extracted bone marrow 
W_br <- 0.0060*W_tot # weight of brain, experimental data - g
W_ht <- 0.0033*W_tot # weight of heart, experimental data - g
W_ki <- 0.0073*W_tot # weight of kidneys, experimental data - g
W_li <- 0.037*W_tot # weight of liver, experimental data - g 
W_spl <- 0.0020*W_tot # weight of spleen, experimental data - g
W_re <- W_tot - W_lu - W_bm - W_br - W_ht - W_ki - W_li - W_spl # weight of rest of the body, experimental data - g

# Weight of capillary blood assuming density = 1
# (capillary blood is given as a percentage of tissue volume)
Wb_lu <- 0.36 * W_lu # weight of blood in lungs, literature (Brown et al., 1997) - g
Wb_bm <- 0.1* W_bm # weight of blood in bone marrow, estimated - g
Wb_br <-  0.03 * W_br  # weight of blood in brain, literature (Brown et al., 1997) - g
Wb_ht <- 0.26 * W_ht # weight of blood in heart, literature (Brown et al., 1997) - g   
Wb_ki <- 0.16 * W_ki # weight of blood in kidneys, literature (Brown et al., 1997) - g 
Wb_li <-  0.21 * W_li # weight of blood in liver, literature (Brown et al., 1997) - g
Wb_spl <- 0.22 * W_spl # weight of blood in spleen, literature (Brown et al., 1997) - g
Wb_re <- 0.04 * W_re # weight of blood in rest of the body, literature (Brown et al., 1997) - g 
# weight of arterial and venous blood, experimental data and literature - g
W_blood <- 0.0816 * (W_tot - (Wb_spl + Wb_li + Wb_lu + Wb_br + Wb_ht + Wb_ki + Wb_bm + Wb_re))

#Regional blood flows (in mL per hour)
fQs = 0.0086 # fraction of cardiac output to spleen, literature (Bernareggi and Rowland, 1991) - unitless
fQl = 0.174 # fraction of cardiac output to liver, literature (chapter 5) - unitless
fQbr = 0.02	# fraction of cardiac output to brain, literature (chapter 5) - unitless
fQh = 0.0448 # fraction of cardiac output to heart, literature (chapter 5) - unitless
fQk = 0.0948 # fraction of cardiac output to kidneys, literature (chapter 5) - unitless
fQbm = 0.0267 # fraction of cardiac output to bone marrow, literature (Brookes, 1967) - unitless
fQrest = 1-fQs-fQl-fQbr-fQh-fQk-fQbm # fraction of cardiac output to rest of the body, fitted value - unitless
Q_tot <- 4980 # cardiac output, literature (chapter 5) - mL/h
Q_bm <- fQbm*Q_tot # blood flow to bone marrow - mL/h
Q_br <- fQbr*Q_tot # blood flow to brain - mL/h		
Q_ht <- fQh*Q_tot	# blood flow to heart - mL/h
Q_ki <- fQk*Q_tot # blood flow to kidneys - mL/h
Q_spl <- fQs*Q_tot # blood flow to spleen - mL/h		
Q_li <- fQl*Q_tot	# blood flow to liver - mL/h
Q_re <- fQrest*Q_tot # blood flow to rest of the body - mL/h

params <- c(W_tot,W_lu,  W_bm, W_br, W_ht, W_ki, W_li, W_spl, W_re, Wb_lu, Wb_bm, Wb_br,Wb_ht,
            Wb_ki, Wb_li, Wb_spl, Wb_re, W_blood, Q_tot, Q_bm, Q_br, Q_ht, Q_ki, Q_spl, Q_li, Q_re,
            dose)

#####################
# Priors #
#####################
#PCs uptake capacity per organ weight (in micrograms per g of tissue)
M_lu_cap <- 25.5 # maximum phagocytizing cells uptake per lung weight, fitted value - ug/g
M_bm_cap <- 41.2 # maximum phagocytizing cells uptake per bone marrow weight, fitted value - ug/g
M_br_cap <- 0.0827 # maximum phagocytizing cells uptake per brain weight, fitted value - ug/g
M_ht_cap <- 5.03 # maximum phagocytizing cells uptake per heart weight, fitted value - ug/g
M_ki_cap <- 1.08 # maximum phagocytizing cells uptake per kidney weight, fitted value - ug/g
M_li_cap <- 74.8 # maximum phagocytizing cells uptake per liver weight, fitted value - ug/g  
M_spl_cap <- 631 # maximum phagocytizing cells uptake per spleen weight, fitted value - ug/g
M_blood_cap <- 0.0396 # maximum phagocytizing cells uptake per blood weight, fitted value - ug/g
M_re_cap <- 17.6 # maximum phagocytizing cells uptake per slowly perfused tissue weight, fitted value - ug/g

x_fast <-1.06e-03 # permeability coefficient from blood to fast perfused tissue, fitted value - unitless
x_br <-0 # permeability coefficient from blood to brain tissue, fitted value - unitless
x_re <-8.25e-05 # permeability coefficient from blood to rest of the body, fitted value - unitless
P <-0.147 # partition coefficient tissue:blood, fitted value - unitless
k_ab0 <- 16.1 # maximum uptake rate by phagocytizing cells, fitted value - 1/h
k_ab0_spl  <- 0.112 # maximum uptake rate by phagocytizing cells in spleen, fitted value - 1/h
k_de <- 4.9e-19 # desorption rate by phagocytizing cells, fitted value - 1/h   
CLE_f <- 1.18e-02 # clearance rate to feces from liver, fitted value - mL/h
CLE_u <- 6.56e-03 # clearance rate to urine from blood in kidneys, fitted value - 1/h
frbr <- 0.346 # fraction of capillary blood remained in brain when measured - unitless
fro <- 0.177 # fraction of capillary blood remained in other organs when measured - unitless

# Time of sampling
time <-  c(10/60, 1, 24, 24*7, 24*28, 24*56)
#read data
data <- openxlsx::read.xlsx("elgrabli_data.xlsx", sheet = 1, colNames = T, rowNames = T)
# Discard control, urine (for the time) and the first row which is the compartment names 
# and also transpose the data so that each column represents a compartment
TransformedData <- t(as.matrix(sapply(data, as.numeric))[1:8,2:7])
# Retrieve the compartment names
colnames(TransformedData) <- rownames(data[1:8,])
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

DataList <- list(
                   eta_tr = c(70.8,1,3.34), #cl, a (kp coeff), kp_fat,kp_rest
                   eta_tr_std = c(62.5,1,0.44), 
                   R = 0.65,
                   N_compart = 14.0,
                   N_subj = 23,
                   N_param = 3.0,
                   params = params,
                   time = time_vec,
                   con = concen,
                   samp = samp, 
                   N_obs = length(concen),
                   t_init = 0,
                   y0 = rep(0,14),
                   rel_tol = 1e-04, 
                   abs_tol = 1e-02,
                   max_num_steps = 1e05,
                   a=15
                 
)

tic = proc.time()

fit <- stan(file = 'test1.stan', data = DataList, 
            iter = 5400, warmup=400, chains=4)
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


