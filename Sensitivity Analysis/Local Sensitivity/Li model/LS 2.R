setwd("C:\\Users\\vassi\\Documents\\Diploma Thesis\\Li_model_tests\\Sensitivity analysis")
results <- matrix(0, nrow = 21, ncol = 16)
Dx <- 1.2   # to pososto metavolis ton parametron


###Epanaliptiki diadikasia ypologismou apotelesmaton gia kathe nea parametro
for (i in 1:16) {

####Parameters


#PCs uptake capacity per organ weight (in micrograms per g of tissue)
M_lu_cap <- 25.5/2 # maximum phagocytizing cells uptake per lung weight, fitted value - ug/g
M_bm_cap <- 41.2 # maximum phagocytizing cells uptake per bone marrow weight, fitted value - ug/g
M_br_cap <- 0.0827 # maximum phagocytizing cells uptake per brain weight, fitted value - ug/g
M_ht_cap <- 5.03/2 # maximum phagocytizing cells uptake per heart weight, fitted value - ug/g
M_ki_cap <- 1.08/2 # maximum phagocytizing cells uptake per kidney weight, fitted value - ug/g
M_li_cap <- 74.8 # maximum phagocytizing cells uptake per liver weight, fitted value - ug/g  
M_spl_cap <- 631 # maximum phagocytizing cells uptake per spleen weight, fitted value - ug/g
M_blood_cap <- 0.0396 # maximum phagocytizing cells uptake per blood weight, fitted value - ug/g
M_re_cap <- 17.6/2 # maximum phagocytizing cells uptake per slowly perfused tissue weight, fitted value - ug/g

x_fast <- 111 # permeability coefficient from blood to fast perfused tissue, fitted value - unitless
x_re <- 0.2 # permeability coefficient from blood to rest of the body, fitted value - unitless
P <-3.8 # partition coefficient tissue:blood, fitted value - unitless
k_ab0 <- 82 # maximum uptake rate by phagocytizing cells, fitted value - 1/h
k_ab0_spl  <- 57 # maximum uptake rate by phagocytizing cells in spleen, fitted value - 1/h
k_de <- 4.9e-19 # desorption rate by phagocytizing cells, fitted value - 1/h   
CLE_f <-  1e-04# clearance rate to feces from liver, fitted value - mL/h
#frbr <- 0.346 # fraction of capillary blood remained in brain when measured - unitless

##Tissue mass (g)
W_tot <-250 # ;body weight, experimental data - g
W_lu <-0.005*W_tot # weight of lungs, literature (Brown et al., 1997) - g
W_bm <- 0.0559*W_tot  # weight of bone marrow, literature (Brown et al., 1997) - g
W_br <- 0.0060*W_tot # weight of brain, literature (Brown et al., 1997) - g
W_ht <- 0.0033*W_tot # weight of heart, literature (Brown et al., 1997) - g
W_ki <- 0.0073*W_tot # weight of kidneys, literature (Brown et al., 1997) - g
W_li <- 0.037*W_tot # weight of liver, literature (Brown et al., 1997) - g 
W_spl <- 0.0020*W_tot # weight of spleen, literature (Brown et al., 1997) - g
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
# weight of arterial and venous blood, literature (Law et al., 2017) - g
W_blood <- 0.0816 * W_tot - (Wb_spl + Wb_li + Wb_lu + Wb_br + Wb_ht + Wb_ki + Wb_bm + Wb_re) 

#Regional blood flows (in mL per hour)
fQs = 0.0086 # fraction of cardiac output to spleen, literature (Sweeney et al., 2014) - unitless
fQl = 0.174 # fraction of cardiac output to liver, literature (Brown et al., 1997) - unitless
fQbr = 0.02	# fraction of cardiac output to brain, literature (Brown et al., 1997) - unitless
fQh = 0.0448 # fraction of cardiac output to heart, literature (Brown et al., 1997) - unitless
fQk = 0.0948 # fraction of cardiac output to kidneys, literature (Brown et al., 1997) - unitless
fQbm = 0.0267 # fraction of cardiac output to bone marrow, literature (Brookes, 1967) - unitless
fQrest = 1-fQs-fQl-fQbr-fQh-fQk-fQbm # fraction of cardiac output to rest of the body, fitted value - unitless
Q_tot <- 14100*(W_tot/1000)^0.75 # cardiac output, literature (Brown et al., 1997) - mL/h, 
Q_bm <- fQbm*Q_tot # blood flow to bone marrow - mL/h
Q_br <- fQbr*Q_tot # blood flow to brain - mL/h		
Q_ht <- fQh*Q_tot	# blood flow to heart - mL/h
Q_ki <- fQk*Q_tot # blood flow to kidneys - mL/h
Q_spl <- fQs*Q_tot # blood flow to spleen - mL/h		
Q_li <- fQl*Q_tot	# blood flow to liver - mL/h
Q_re <- fQrest*Q_tot # blood flow to rest of the body - mL/h

CLE_u <- 0 # clearance rate to urine from blood in kidneys, fitted value - 1/h
x_br <- 0 # permeability coefficient from blood to brain tissue, fitted value - unitless


params <- c(M_lu_cap, M_bm_cap, M_br_cap, M_ht_cap, M_ki_cap, M_li_cap, M_spl_cap,
            M_blood_cap, M_re_cap, x_fast, x_re, P, k_ab0, k_ab0_spl, k_de, CLE_f, #frbr,  #mexri edo einai oi metavlites stis opoies ginetai sensitivity analysis
            W_lu,  W_bm, W_br, W_ht, W_ki, W_li, W_spl, W_re, Wb_lu, Wb_bm, Wb_br,Wb_ht,
            Wb_ki, Wb_li, Wb_spl, Wb_re, W_blood, Q_tot, Q_bm, Q_br, Q_ht, Q_ki, Q_spl, Q_li, Q_re,
            CLE_u, x_br)

  init_params <- params        # xrisimopoieitai pio kato gia ton ypologismo arxikon apotelesmaton xoris kapoia metavoli stis parametrous
  params[i] <- Dx*params[i]    # epivoli metavolis tis parametrou i kata Dx (oristike stin arxi) 
  
  
  
  ###Sto parakato kommati ginetai antistoixisi twn newn timwn twn parametrwn me tis antistoixes onomasies tous 
  M_lu_cap <- params[1] # maximum phagocytizing cells uptake per lung weight, fitted value - ug/g
  M_bm_cap <- params[2] # maximum phagocytizing cells uptake per bone marrow weight, fitted value - ug/g
  M_br_cap <- params[3] # maximum phagocytizing cells uptake per brain weight, fitted value - ug/g
  M_ht_cap <- params[4] # maximum phagocytizing cells uptake per heart weight, fitted value - ug/g
  M_ki_cap <- params[5] # maximum phagocytizing cells uptake per kidney weight, fitted value - ug/g
  M_li_cap <- params[6] # maximum phagocytizing cells uptake per liver weight, fitted value - ug/g  
  M_spl_cap <- params[7] # maximum phagocytizing cells uptake per spleen weight, fitted value - ug/g
  M_blood_cap <- params[8] # maximum phagocytizing cells uptake per blood weight, fitted value - ug/g
  M_re_cap <- params[9] # maximum phagocytizing cells uptake per slowly perfused tissue weight, fitted value - ug/g
  
  x_fast <- params[10] # permeability coefficient from blood to fast perfused tissue, fitted value - unitless
  x_re <- params[11] # permeability coefficient from blood to rest of the body, fitted value - unitless
  P <- params[12] # partition coefficient tissue:blood, fitted value - unitless
  k_ab0 <- params[13] # maximum uptake rate by phagocytizing cells, fitted value - 1/h
  k_ab0_spl  <- params[14] # maximum uptake rate by phagocytizing cells in spleen, fitted value - 1/h
  k_de <- params[15] # desorption rate by phagocytizing cells, fitted value - 1/h   
  CLE_f <-  params[16]
  
  
  
  source("Li-equations.R")     # epilisi toy montelou gia ti nea timi tis parametrou i

  solution <- solution[2,2:22] # kratao ta apotelesmata mono gia 1 xroniki stigmh epilisis kai petao ta ipoloipa apotelesmata kathos kai tin timi tou xronou oloklirosis
  
  
  results[,i] <- solution
  
}

#results

###Ypologismos apotelesmaton xoris kapoia metavoli stis parametrous
params <- init_params
source("Li-equations.R")          # epilisi toy montelou gia ti arxikes times ton parametron
init_solution <- solution[2,2:22] # kratao ta apotelesmata mono gia 1 xroniki stigmh epilisis kai petao ta ipoloipa apotelesmata kathos kai tin timi tou xronou oloklirosis


###Ypologismos metavolis twn diaforikwn dydt ws pros tin metavoli tis parametrou poy elegxetai
for (y in 1:16) {
  for (x in 1:21) {
  results[x,y] <- abs((results[x,y] - init_solution[x]))/(Dx*params[y])     
  }
  
}


results


###Sensitivity Index (SI) calculation
SI <- matrix(0, nrow = 21, ncol = 16)
for (k in 1:21) {
  SI[k,] <- results[k,]/init_solution[k]
}
SI