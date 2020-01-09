library(deSolve)
dE <- matrix(0, ncol=30)

##############
# User Input #
##############

weight <- 250 #in g
dose <- 1.7e-03 #in grams
infusion_time <- 0 


##############
# Parameters #
##############

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


params <- c(V_gi, V_br, V_ht, V_ki, V_li, V_spl, V_re, Vven, Vart, Vb_gi, Vb_br,Vb_ht,
            Vb_ki, Vb_li, Vb_spl, Vb_re, V_blood, Q_tot, Q_gi, Q_br, Q_ht, Q_ki, Q_spl, Q_li, Q_re, Da, Du, Do, kB, ko, ku, kt,
            kr, kd, ki, fo, fu, ft, lamdav, lamda31, lamda32, lamda33, lamda41, lamda42, lamda43, lamda51, lamda52, lamda53, lamda61, lamda62, lamda63, lamda71,
            lamda72, lamda73, lamda81,lamda82, lamda83, lamda91, lamda92, lamda93, lamda35, lamda44, lamda54, lamdaI, kb, kl)

#################################################
# Function for creating initial values for ODEs #
#################################################

#E <- matrix(0, ncol=30)
E <- c(rep(0,29),dose)
inits <- as.vector(E)


###############
# ODEs system #
###############

func <- function(time=0, inits, params){
  with( as.list(c(inits, params)),{
    
    dE[1]  = Do - (ko*E[1]) - (kB*E[1])   # Olfactory
    dE[2]  = Du - (ku*E[2]) # Upper airways
    dE[3]  = Da - (kr*E[3])  + (kd*E[4]) - (ki*E[3]) #Alveolar free
    dE[4]  = (kr*E[3]) - (kd*E[4]) - (kt*E[4]) #Alveolar Mac
    dE[5]  = (ki*E[3]) - (kl*E[5]) - (kb*E[5]) + ((E[30]/Vven)*lamdaI*Q_tot*lamdav) #I
    dE[6]  = kl*E[5];  #Lymph to arterial  from venous Blood
    
    #Liver
    dE[7]  = ((E[29]/Vart)*Q_li) - ((E[7]/Vb_li)*Q_li) + ((E[8]/V_li)*lamda31*Q_li) - ((E[7]/Vb_li)*Q_li*lamda32) - 
              ((E[7]/Vb_li)*Q_li*lamda35) + ((E[10]/Vb_gi)*Q_gi) #capillary
    dE[8] =  ((E[7]/Vb_li)*Q_li) - ((E[8]/V_li)*lamda31*Q_li) - (lamda33*E[8]) #tissue
    dE[9]  = lamda33*E[8] #seq
    
    #GI
    dE[10]  = (fo*ko*E[1]) + (fu*ku*E[2]) + (ft*kt*E[4]) + (E[7]/Vb_li)*Q_li*lamda35 + ((E[29]/Vart)*Q_gi) - 
              ((E[10]/Vb_gi)*Q_gi) + ((E[11]/V_gi)*lamda41*Q_gi)  - (lamda44*E[10])  #capillary
    dE[11]  = ((E[10]/Vb_gi)*Q_gi) - ((E[11]/V_gi)*lamda41*Q_gi) - (lamda43*E[11]) #tissue
    dE[12]  = lamda43*E[11] #seq
    dE[13]  = 0;
    
    #Kidneys
    dE[14]  = ((E[29]/Vart)*Q_ki) - ((E[14]/Vb_ki)*Q_ki) + ((E[15]/V_ki)*lamda51*Q_ki) - 
               ((E[14]/Vb_ki)*Q_ki*lamda52) - (E[14]*lamda54) #capillary
    dE[15] =  ((E[14]/Vb_ki)*Q_ki) - ((E[15]/V_ki)*lamda51*Q_ki) - (lamda53*E[15])  #tissue
    dE[16]  = lamda53*E[15] #seq
    
    #Heart
    dE[17]  = ((E[29]/Vart)*Q_ht) - ((E[17]/Vb_ht)*Q_ht) + ((E[18]/V_ht)*lamda61*Q_ht) -
              ((E[17]/Vb_ht)*Q_ht*lamda62) #capillary
    dE[18] =  ((E[17]/Vb_ht)*Q_ht) - ((E[18]/V_ht)*lamda61*Q_ht) - (lamda63*E[18]) #tissue
    dE[19]  = lamda63*E[18] #seq
    
    #Spleen
    dE[20]  = ((E[29]/Vart)*Q_spl) - ((E[20]/Vb_spl)*Q_spl) + ((E[21]/V_spl)*lamda71*Q_spl) - 
              ((E[20]/Vb_spl)*Q_spl*lamda72) #capillary
    dE[21] =  ((E[20]/Vb_spl)*Q_spl) - ((E[21]/V_spl)*lamda71*Q_spl) - (lamda73*E[21]) #tissue
    dE[22]  = lamda73*E[21] #seq
    
    #Brain
    dE[23]  =  ((E[29]/Vart)*Q_br) - ((E[23]/Vb_br)*Q_br) + ((E[24]/V_br)*lamda81*Q_br) - 
               ((E[23]/Vb_br)*Q_br*lamda82) #capillary
    dE[24]  =  kB*E[1] + ((E[23]/Vb_br)*Q_br) - ((E[24]/V_br)*lamda81*Q_br) - (lamda83*E[24]) #tissue
    dE[25]  =  lamda83*E[24] #seq
    
    #Others
    dE[26]  = ((E[29]/Vart)*Q_re) - ((E[26]/Vb_re)*Q_re) + ((E[27]/V_re)*lamda91*Q_re) -
              ((E[26]/Vb_re)*Q_re*lamda92) #capillary
    dE[27] =  ((E[26]/Vb_re)*Q_re) - ((E[27]/V_re)*lamda91*Q_re) - (lamda93*E[27]) #tissue
    dE[28]  = lamda93*E[27] #seq
    
    #Blood
    dE[29] =  kb*E[5] - ((E[29]/Vart)*Q_tot) #art
    dE[30] =  ((E[23]/Vb_br)*Q_br*lamda82)  + ((E[7]/Vb_li)*Q_li*lamda32) + ((E[17]/Vb_ht)*Q_ht*lamda62) + 
              ((E[20]/Vb_spl)*Q_spl*lamda72) +
              ((E[26]/Vb_re)*Q_re*lamda92) + ((E[14]/Vb_ki)*Q_ki*lamda52) - ((E[30]/Vven)*lamdaI*Q_tot) #ven
    
    return(list(c(dE[1] , dE[2], dE[3], dE[4], dE[5], dE[6], dE[7], dE[8], dE[9], dE[10], dE[11], dE[12], dE[13], 
                  dE[14], dE[15], dE[16], dE[17], dE[18] ,dE[19], dE[20], dE[21], dE[22], dE[23], dE[24], dE[25],
                  dE[26], dE[27], dE[28], dE[29], dE[30]))) 
    
  })
}
  


###########
sample_time <- c(0, 10/(60*24), 1/24, 1, 7, 28, 56) #in days
solution <- ode(times = sample_time, func = func, y = inits, parms = params, method = "lsodes" )

solution
