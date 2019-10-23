library(deSolve)        
# 1:lungs (lu), 2:rest of the body (re), 3:bone marrow (bm), 4:brain (br), 5:heart (ht)
# 6:kidneys (ki), 7:liver (li), 8:spleen (spl), 9:arterial blood (art), 10:venous blood (ven)
options(max.print=999999)

##############
# User input #
##############

dose <- c(7000)  # in mg
times <- c(0.000001) #infusion time in hours
compartment <- "Ven_blood"

user_input <-list("dose" = dose, "times" = times, "compartment" = compartment)
predicted.feats <- c("Lu_tissue", "Rob_tissue", "Bm_tissue", "Br_tissue", "Ht_tissue", "Ki_tissue", 
                     "Li_tissue", "Spl_tissue", "Art_blood", "Ven_blood", "Lu_pc", "Rob_pc", 
                     "Bm_pc", "Br_pc", "Ht_pc", "Ki_pc","Li_pc", "Spl_pc", "Blood_pc","Cle_f", 
                     "Cle_u","Lu_total", "Rob_total", "Bm_total", "Br_total", "Ht_total", "Ki_total",
                     "Li_total", "Spl_total", "Blood_total")

##########################################
# Function for creating parameter vector #
##########################################

create.params <- function(input){
  with( as.list(input),{
    
  ############################
  # Physiological parameters #
  ############################
  # tissue weights (in g)
  W_tot <-253.4 # ;body weight, experimental data - g
  W_lu <-1.2 # weight of lungs, experimental data - g
  W_bm <- W_tot*0.03  # weight of bone marrow, literature (Travlos 2006) - g
  W_bm_exp <- 0.02 # weight of extracted bone marrow 
  W_br <- 1.35 # weight of brain, experimental data - g
  W_ht <- 0.87 # weight of heart, experimental data - g
  W_ki <- 2.37 # weight of kidneys, experimental data - g
  W_li <- 10.03 # weight of liver, experimental data - g 
  W_spl <- 0.78 # weight of spleen, experimental data - g
  W_re <- 212.8 # weight of rest of the body, experimental data - g
  
  # Weight of capillary blood assuming density = 1
  # (capillary blood is given as a percentage of tissue volume)
  Wb_lu <- 0.36 * W_lu # weight of blood in lungs, literature (Brown et al., 1997) - g
  Wb_bm <- 0.1* W_bm # weight of blood in bone marrow, estimated - g
  Wb_br <-  0.07 * W_br  # weight of blood in brain, literature (Brown et al., 1997) - g
  Wb_ht <- 0.26 * W_ht # weight of blood in heart, literature (Brown et al., 1997) - g   
  Wb_ki <- 0.16 * W_ki # weight of blood in kidneys, literature (Brown et al., 1997) - g 
  Wb_li <-  0.21 * W_li # weight of blood in liver, literature (Brown et al., 1997) - g
  Wb_spl <- 0.22 * W_spl # weight of blood in spleen, literature (Brown et al., 1997) - g
  Wb_re <- 0.04 * W_re # weight of blood in rest of the body, literature (Brown et al., 1997) - g 
  # weight of arterial and venous blood, experimental data and literature - g
  W_blood <- 0.065 * W_tot - (Wb_spl + Wb_li + Wb_lu + Wb_br + Wb_ht + Wb_ki + Wb_bm + Wb_re)
  
  #Regional blood flows (in mL per hour)
  fQs = 0.0146 # fraction of cardiac output to spleen, literature (Bernareggi and Rowland, 1991) - unitless
  fQl = 0.183 # fraction of cardiac output to liver, literature (chapter 5) - unitless
  fQbr = 0.02	# fraction of cardiac output to brain, literature (chapter 5) - unitless
  fQh = 0.051 # fraction of cardiac output to heart, literature (chapter 5) - unitless
  fQk = 0.141 # fraction of cardiac output to kidneys, literature (chapter 5) - unitless
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
  
  #####################
  # Fitted parameters #
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
  
  return(list(  "W_tot" = W_tot, "W_lu" = W_lu, "W_bm" = W_bm, "W_bm_exp" = W_bm_exp, "W_br" = W_br, 
                "W_ht" = W_ht, "W_ki" = W_ki, "W_li" = W_li, "W_spl" = W_spl, "W_re" = W_re, "Wb_lu" = Wb_lu,
                "Wb_bm" = Wb_bm, "Wb_br" = Wb_br, "Wb_ht" = Wb_ht, "Wb_ki" = Wb_ki, "Wb_li" = Wb_li, 
                "Wb_spl" = Wb_spl, "Wb_re" = Wb_re, "W_blood" = W_blood, "fQs" = fQs, "fQl" = fQl, 
                "fQbr" = fQbr, "fQh" = fQh, "fQk" = fQk, "fQbm" = fQbm, "fQrest" = fQrest,"Q_tot" =  Q_tot,
                "Q_bm" = Q_bm, "Q_br" = Q_br, "Q_ht" = Q_ht, "Q_ki" = Q_ki, "Q_spl" = Q_spl, "Q_li" = Q_li,
                "Q_re" = Q_re, "M_lu_cap" = M_lu_cap, "M_bm_cap" = M_bm_cap, "M_br_cap" = M_br_cap, 
                "M_ht_cap" = M_ht_cap, "M_ki_cap" = M_ki_cap, "M_li_cap" = M_li_cap, "M_spl_cap" = M_spl_cap,
                "M_blood_cap" = M_blood_cap, "M_re_cap" = M_re_cap, "x_fast" = x_fast, "x_br" = x_br,
                "x_re" = x_re, "P" = P, "k_ab0" = k_ab0, "k_ab0_spl" = k_ab0_spl, "k_de" = k_de, 
                "CLE_f" = CLE_f, "CLE_u" = CLE_u, "frbr" = frbr, "fro" = fro, "dose" = dose,
                "times" = times, compartment))
  }) 
}

### store them once
params <- create.params(user_input)

#################################################
# Function for creating initial values for ODEs #
#################################################

create.inits <- function(parameters){
  with( as.list(parameters),{
    Lu_tissue <- 0; Rob_tissue <- 0; Bm_tissue <- 0; Br_tissue <- 0; Ht_tissue <- 0; Ki_tissue <- 0; 
    Li_tissue <- 0; Spl_tissue <- 0; Art_blood <- 0; Ven_blood <- 0; Lu_pc<- 0; Rob_pc <- 0; 
    Bm_pc <- 0; Br_pc <- 0; Ht_pc <- 0; Ki_pc <- 0;Li_pc <- 0; Spl_pc <- 0; Blood_pc <- 0; feces<- 0; 
    urine<- 0; 
    
    return(c("Lu_tissue" = Lu_tissue, "Rob_tissue" = Rob_tissue, "Bm_tissue" = Bm_tissue, 
             "Br_tissue" = Br_tissue, "Ht_tissue" = Ht_tissue, "Ki_tissue" = Ki_tissue, 
              "Li_tissue" = Li_tissue, "Spl_tissue" = Spl_tissue, "Art_blood" = Art_blood, 
             "Ven_blood" = Ven_blood, "Lu_pc" = Lu_pc, "Rob_pc" = Rob_pc, 
              "Bm_pc" = Bm_pc, "Br_pc" = Br_pc, "Ht_pc" = Ht_pc, "Ki_pc" = Ki_pc, "Li_pc" = Li_pc,
             "Spl_pc" = Spl_pc, "Blood_pc" = Blood_pc, "feces" = feces, 
              "urine" = urine))
  }) 
}
##store the values
inits <- create.inits(params)

#################################################
# Function for creating events #
#################################################
create.events<- function(parameters){
  with( as.list(parameters),{
    
    ldose <- length(dose)
    ltimes <- length(times)
    
    addition <- dose
    if (ltimes == ldose){
      events <- list(data = rbind(data.frame(var = c(compartment),  time = times, 
                                             value = addition, method = c("add")) ))
    }else{
      stop("The times when the drug is injected should be equal in number to the doses")
    }
    
    
    return(events)
  }) 
}

events <- create.events(params)

###################
# Custom function #
###################

custom.func <- function(){
  return()
}

#################
# ODEs system #
#################

ode.func <- function(time, Initial.values, Parameters, custom.func){
  with( as.list(c(Initial.values, Parameters)),{
 
  # concentrations in tissues
  C_lu <- Lu_tissue/W_lu
  C_re  <-  Rob_tissue/W_re
  C_bm  <- Bm_tissue/W_bm
  C_br  <-  Br_tissue/W_br
  C_ht  <-  Ht_tissue/W_ht
  C_ki  <-  Ki_tissue/W_ki
  C_li  <-  Li_tissue/W_li
  C_spl <- Spl_tissue/ W_spl
  C_art <- Art_blood/(0.2*(W_blood + Wb_spl + Wb_li + Wb_lu + Wb_br + Wb_ht + Wb_ki + Wb_bm + Wb_re))
  C_ven <- Ven_blood/(0.8*(W_blood + Wb_spl + Wb_li + Wb_lu + Wb_br + Wb_ht + Wb_ki + Wb_bm + Wb_re))
  
  # Uptake rates by phagocytizing cells
  kluab <- k_ab0*(1-(Lu_pc/(M_lu_cap*W_lu)))
  kreab <- k_ab0*(1-(Rob_pc/(M_re_cap*W_re)))
  kbmab <- k_ab0*(1-(Bm_pc/(M_bm_cap*W_bm)))
  kbrab <- k_ab0*(1-(Br_pc/(M_br_cap*W_br)))
  khtab <- k_ab0*(1-(Ht_pc/(M_ht_cap*W_ht)))
  kkiab <- k_ab0*(1-(Ki_pc/(M_ki_cap*W_ki)))
  kliab <- k_ab0*(1-(Li_pc/(M_li_cap*W_li)))
  ksplab <- k_ab0_spl*(1-(Spl_pc/(M_spl_cap*W_spl)))
  kbloodab <- k_ab0*(1-(Blood_pc/(M_blood_cap*W_blood)))
  
  # Nanoparticles in tissue
  # lungs
  dLu_tissue <- ((x_fast*Q_tot)/(1+x_fast)) * (C_ven - C_lu/P) - (W_lu*C_lu*kluab - Lu_pc*k_de)
  # rest of the body
  dRob_tissue <- ((x_re*Q_re)/(1+x_re)) * (C_art - C_re/P) - (W_re*C_re*kreab - Rob_pc*k_de)
  # bone marrow
  dBm_tissue <- ((x_fast*Q_bm)/(1+x_fast)) * (C_art - C_bm/P) - (W_bm*C_bm*kbmab - Bm_pc*k_de)
  # brain
  dBr_tissue <- ((x_br*Q_br)/(1+x_br)) * (C_art - C_br/P) - (W_br*C_br*kbrab - Br_pc*k_de )
  # heart
  dHt_tissue <- ((x_fast*Q_ht)/(1+x_fast)) * (C_art - C_ht/P) - (W_ht*C_ht*khtab - Ht_pc*k_de )
  # Kidneys
  dKi_tissue <- ((x_fast*Q_ki)/(1+x_fast)) * ( C_art - C_ki/P ) -(W_ki*C_ki*kkiab - Ki_pc*k_de )-
             C_art*CLE_u*x_fast/(1+x_fast) 
  # Liver
  dLi_tissue <- C_art*x_fast*Q_li/(x_fast+1) + Q_spl*x_fast*(C_art+x_fast*C_spl/P)/((1+x_fast)*(1+x_fast))-
            (((C_li/P)*x_fast*(Q_li+Q_spl))/(x_fast+1))-(W_li*C_li*kliab- Li_pc*k_de) - Li_tissue*CLE_f
  # spleen
  dSpl_tissue <- ((x_fast*Q_spl)/(1+x_fast)) * (C_art - C_spl/P) - (W_spl*C_spl*ksplab - Spl_pc*k_de)
  # arterial blood
  dArt_blood <- Q_tot *((C_ven + x_fast*C_lu/P)/(1+x_fast)-C_art)-((0.2*W_blood*C_art)*kbloodab-0.2*Blood_pc*k_de)
  # venous blood
  dVen_blood <- C_li*x_fast*(Q_li+Q_spl)/(P*(1+x_fast)) + Q_spl*x_fast*C_spl/(P*(1+x_fast)*
                (1+x_fast)) + (C_br*Q_br*x_br)/(P*(1+x_br)) + (C_ht*x_fast*Q_ht)/(P*(1+x_fast))+
               (C_ki*Q_ki*x_fast)/(P*(1+x_fast)) + (C_bm*Q_bm*x_fast)/(P*(1+x_fast))+
                (Q_re*C_re*x_re)/(P*(1+x_re)) + (Q_li/(1+x_fast) + Q_spl/((1+x_fast)*(1+x_fast))+
                Q_br/(1+x_br) +  Q_ht/(1+x_fast) +  Q_ki*(1-CLE_u/Q_ki)/(1+x_fast) + Q_bm/(1+x_fast)+
                Q_re/(1+x_re))*C_art - Q_tot*C_ven - ((0.8*W_blood*C_ven)* kbloodab -0.8* Blood_pc * k_de) 
  
  
  # Nanoparticles uptaken in PCs
  # lungs
  dLu_pc <- W_lu*C_lu*kluab - Lu_pc*k_de                                   
  # rest of the body                                    
  dRob_pc <- W_re*C_re*kreab - Rob_pc*k_de    
  # bone marrow                                    
  dBm_pc <- W_bm*C_bm*kbmab - Bm_pc*k_de  
  # brain                                  
  dBr_pc <- W_br*C_br*kbrab - Br_pc*k_de
  # heart                                  
  dHt_pc <- W_ht*C_ht*khtab - Ht_pc*k_de 
  # kidneys                                
  dKi_pc <- W_ki*C_ki*kkiab - Ki_pc*k_de
  # liver                               
  dLi_pc <- W_li*C_li*kliab - Li_pc*k_de
  # spleen                               
  dSpl_pc <- W_spl*C_spl*ksplab - Spl_pc*k_de 
  # blood
  dBlood_pc <- (0.2*W_blood*C_art + 0.8*W_blood*C_ven)*kbloodab - Blood_pc*k_de 
  
  #Nanoparticles in excreta
  dfeces <- Li_tissue*CLE_f # Rate of amount in feces from liver - ug/h
  durine<- C_art*CLE_u  # Rate of amount in urine from kidney -ug/h
  
  # Total amount of NPs in each organ
  # Amount in lungs
  Lu_total <- Lu_tissue + Lu_pc + fro*(C_ven+x_fast*C_lu/P)/(1+x_fast)*Wb_lu
 # Amount in bone marrow
  Bm_total <- Bm_tissue + Bm_pc + fro*(C_art+x_fast*C_bm/P)/(1+x_fast)*Wb_bm
  # Amount in rest of the body
  Rob_total <- Rob_tissue + Rob_pc + fro*(C_art+x_re*C_re/P)/(1+x_re)*Wb_re + Bm_total*(1-(W_bm_exp/W_bm))
  # Amount in brain
  Br_total <- Br_tissue + Br_pc + frbr*(C_art+x_br*C_br/P)/(1+x_br)*Wb_br
  # Amount in heart
  Ht_total <- Ht_tissue + Ht_pc + fro*(C_art+x_fast*C_ht/P)/(1+x_fast)*Wb_ht
  # Amount in kidneys
  Ki_total <- Ki_tissue + Ki_pc + fro*(C_art+x_fast*C_ki/P)/(1+x_fast)*Wb_ki
  # Amount in liver
  Li_total <- Li_tissue + Li_pc + fro*(C_art+x_fast*C_li/P)/(1+x_fast)*Wb_li
  # Amount in spleen
  Spl_total <- Spl_tissue + Spl_pc + fro*(C_art+x_fast*C_spl/P)/(1+x_fast)*Wb_spl
  # Amount in blood
  Blood_total <- Art_blood+Ven_blood+Blood_pc
  
  list(c(dLu_tissue = dLu_tissue, dRob_tissue = dRob_tissue, dBm_tissue = dBm_tissue, dBr_tissue = dBr_tissue, 
         dHt_tissue = dHt_tissue, dKi_tissue = dKi_tissue, dLi_tissue = dLi_tissue, dSpl_tissue = dSpl_tissue,
         dArt_blood = dArt_blood, dVen_blood = dVen_blood, dLu_pc = dLu_pc, dRob_pc = dRob_pc, 
         dBm_pc = dBm_pc, dBr_pc = dBr_pc, dHt_pc = dHt_pc, dKi_pc = dKi_pc, dLi_pc = dLi_pc, dSpl_pc = dSpl_pc,
         dBlood_pc = dBlood_pc, dfeces = dfeces, durine= durine), Lu_total = Lu_total, Rob_total = Rob_total, 
         Bm_total = Bm_total, Br_total = Br_total, Ht_total = Ht_total, Ki_total = Ki_total,
       Li_total = Li_total, Spl_total = Spl_total, Blood_total = Blood_total)
  })
}

##############################################

sample_time <- seq(0,120,0.1) # in hours
solution <- ode(times = sample_time,  func = ode.func, y = inits, parms = params, 
                custom.func = custom.func, method="lsodes",  events = events)
#solution

#library(xlsx)
write.xlsx(solution, file="Li_data2.xlsx", sheetName="sheet1", row.names=FALSE)

rowSums(solution[,23:31])
