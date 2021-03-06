setwd("C:\\Users\\vassi\\Documents\\Diploma Thesis\\Li_model_tests\\Sensitivity analysis")
options(max.print=1000000)
library(ggplot2)

Dx <- 1.2   # to pososto metavolis ton parametron
ar <- array(0, c(6, 21, 16))
ar2 <- array(0, c(6, 16, 21))


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

  ar[,,i] <- solution[,2:22] #dimiourgia 3d array (6x21x16), diladi ta apotelesmata twn 21 diaforikwn se 6 xronikes stigmes epilisis gia metavoli se kathe mia apo tis 16 parametrous
  
}

###Metatropi tou pinaka ar se diastaseis 6x16x21
for (z in 1:16) {
  for (w in 1:21) {
      ar2[,z,w] <- ar[,w,z]  
      }
} 

#ar2

###Ypologismos apotelesmaton xoris kapoia metavoli stis parametrous
params <- init_params

for (i in 1:16) {
  #Sto parakato kommati ginetai antistoixisi twn newn timwn twn parametrwn me tis antistoixes onomasies tous 
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
  
}

source("Li-equations.R")          # epilisi toy montelou gia ti arxikes times ton parametron
init_solution <- solution[,2:22] # kratao ta apotelesmata mono gia oles tis xronikes stigmes epilisis



###Ypologismos metavolis twn diaforikwn dydt ws pros tin metavoli tis parametrou poy elegxetai
results <- array(0, c(6, 16, 21)) 
for (w in 1:21) {      ### deiktis diamerismatos
    for (t in 1:6) {   ### deiktis xronikis stigmis
      results[t,,w] <- abs(ar2[t,,w]-init_solution[t,w])
      }
}

for (p in 1:16) { # opou p --> deiktis gia kathe allagmeni parametro
  results[,p,] <- results[,p,]/(Dx*params[p])
}
#results


###Sensitivity Index (SI) calculation
SI <- array(0, c(6, 16, 21))
eps <- 1e-10
for (k in 1:21) {
  for (t in 1:6) {
    SI[t,,k] <- results[t,,k]/(init_solution[t,k]+eps)
  }
}
#SI


#Dimiourgia data frame gia kathe diamerisma gia ta SI olon ton parametrwn      ## Nanoparticles in tissue
data_comp1 <- as.data.frame(rbind(cbind(sample_time, SI[,,1])))                # lungs
data_comp2 <- as.data.frame(rbind(cbind(sample_time, SI[,,2])))                # rest of body
data_comp3 <- as.data.frame(rbind(cbind(sample_time, SI[,,3])))                # bone marrow
data_comp4 <- as.data.frame(rbind(cbind(sample_time, SI[,,4])))                # brain
data_comp5 <- as.data.frame(rbind(cbind(sample_time, SI[,,5])))                # heart
data_comp6 <- as.data.frame(rbind(cbind(sample_time, SI[,,6])))                # kidneys
data_comp7 <- as.data.frame(rbind(cbind(sample_time, SI[,,7])))                # liver
data_comp8 <- as.data.frame(rbind(cbind(sample_time, SI[,,8])))                # spleen
data_comp9 <- as.data.frame(rbind(cbind(sample_time, SI[,,9])))                # arterial blood
data_comp10 <- as.data.frame(rbind(cbind(sample_time, SI[,,10])))              # venous blood

                                                                               ## Nanoparticles uptaken in PCs
data_comp11 <- as.data.frame(rbind(cbind(sample_time, SI[,,11])))              # lungs
data_comp12 <- as.data.frame(rbind(cbind(sample_time, SI[,,12])))              # rest of the body
data_comp13 <- as.data.frame(rbind(cbind(sample_time, SI[,,13])))              # bone marrow
data_comp14 <- as.data.frame(rbind(cbind(sample_time, SI[,,14])))              # brain
data_comp15 <- as.data.frame(rbind(cbind(sample_time, SI[,,15])))              # heart
data_comp16 <- as.data.frame(rbind(cbind(sample_time, SI[,,16])))              # kidneys
data_comp17 <- as.data.frame(rbind(cbind(sample_time, SI[,,17])))              # liver
data_comp18 <- as.data.frame(rbind(cbind(sample_time, SI[,,18])))              # spleen
data_comp19 <- as.data.frame(rbind(cbind(sample_time, SI[,,19])))              # blood
data_comp20 <- as.data.frame(rbind(cbind(sample_time, SI[,,20])))              #feces
data_comp21 <- as.data.frame(rbind(cbind(sample_time, SI[,,21])))              #urine

colnames(data_comp1) <- c("Time", "M_lu_cap", "M_bm_cap", "M_br_cap", "M_ht_cap", "M_ki_cap", "M_li_cap", "M_spl_cap", "M_blood_cap", "M_re_cap", "x_fast", "x_re", "P", "k_ab0", "k_ab0_spl", "k_de", "CLE_f")
colnames(data_comp2) <- colnames(data_comp1)
colnames(data_comp3) <- colnames(data_comp1)
colnames(data_comp4) <- colnames(data_comp1)
colnames(data_comp5) <- colnames(data_comp1)
colnames(data_comp6) <- colnames(data_comp1)
colnames(data_comp7) <- colnames(data_comp1)
colnames(data_comp8) <- colnames(data_comp1)
colnames(data_comp9) <- colnames(data_comp1)
colnames(data_comp10) <- colnames(data_comp1)
colnames(data_comp11) <- colnames(data_comp1)
colnames(data_comp12) <- colnames(data_comp1)
colnames(data_comp13) <- colnames(data_comp1)
colnames(data_comp14) <- colnames(data_comp1)
colnames(data_comp15) <- colnames(data_comp1)
colnames(data_comp16) <- colnames(data_comp1)
colnames(data_comp17) <- colnames(data_comp1)
colnames(data_comp18) <- colnames(data_comp1)
colnames(data_comp19) <- colnames(data_comp1)
colnames(data_comp20) <- colnames(data_comp1)
colnames(data_comp21) <- colnames(data_comp1)


ggplot(data_comp1, aes(x=Time)) +
  geom_point(data = data_comp1, aes(x=Time, y=M_lu_cap, color = "M_lu_cap"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=M_bm_cap, color = "M_bm_cap"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=M_br_cap,  color = "M_br_cap"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=M_ht_cap,  color = "M_ht_cap"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=M_ki_cap,  color = "M_ki_cap"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=M_li_cap,  color = "M_li_cap"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=M_spl_cap,  color = "M_spl_cap"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=M_blood_cap,  color = "M_blood_cap"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=M_re_cap,  color = "M_re_cap"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=x_fast,  color = "x_fast"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=x_re,  color = "x_re"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=P,  color = "P"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=k_ab0,  color = "k_ab0"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=k_ab0_spl,  color = "k_ab0_spl"),  size = 3) +
  #geom_point(data = data_comp1, aes(x=Time, y=k_de, color = "k_de"),  size = 3) +
  geom_point(data = data_comp1, aes(x=Time, y=CLE_f, color = "CLE_f"),  size = 3) +
  
  labs(title = "SI vs Time", subtitle = "Lungs compartment", y = "SI", x = "Time (in minutes)") +
  scale_colour_manual(name = "Parameters",
                     breaks = c("M_lu_cap", "M_bm_cap", "M_br_cap", "M_ht_cap", "M_ki_cap", "M_li_cap", "M_spl_cap", "M_blood_cap", "M_re_cap", "x_fast", "x_re", "P", "k_ab0", "k_ab0_spl",  "CLE_f"),
                     values = c("grey", "red", "royalblue", "pink", "navy", "maroon", "orange", "yellow", "violetred", "rosybrown", "khaki", "hotpink", "cyan", "salmon", "black")) +
  theme(legend.title=element_text(hjust = 0.5,size=17), 
        legend.text=element_text(size=14))
