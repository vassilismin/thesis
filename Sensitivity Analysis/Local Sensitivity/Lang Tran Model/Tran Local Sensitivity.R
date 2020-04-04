setwd("C:\\Users\\vassi\\Documents\\Diploma Thesis\\Lang Tran Files\\Sensitivity Analysis\\Local Sensitivity")
options(max.print=1000000)
library(deSolve)
library(ggplot2)

init <- 1
Dx <- 0.2   # to pososto metavolis ton parametron
mult <- init + Dx
ar <- array(0, c(6, 6, 37))  # 6 xronikes stigmes epilisis, 6 diamerismata, 37 parametroi pros analisi 
ar2 <- array(0, c(6, 37, 6))

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
  Depa <- 03 # Alveolar deposition fraction
  conc <- 0 # External concentration
  VI <- 0.18 # Volume inhaled (L/min)
  CONV <- 1.44 # Conversion factor
  
  Da <- Depo * conc * VI * CONV # Olfactory deposition rate [micro grams/day]
  Du <- Depu * conc * VI * CONV # Upper airways deposition rate [micro grams/day]
  Do <- Depa * conc * VI * CONV # Alveolar deposition rate [micro grams/day]
  
  # lamdai_1: plasma to tissue PC [unitless]
  # lamdai_2: fractional translocation from capillaries to venous blood [unitless]
  # lamdai_3: sequestration rate [d^-1]  
  # lamda3_5: bile:liver PC [unitless]
  # lamda4_5: GI to liver PC [unitless]
  # lamda4_6 GI absorption rate [d^-1]  
  # lamda4_4: fecal elimination rate [d^-1]  
  
  lambda3 <- c(4.39, 0.0001, 0.24, 0, 1.73)  #lambda3[5] taken from code for semmler
  lambda4 <- c(2, 0.29185, 0.07, 4, 1, 1.44) #lambda4[6] taken from code for semmler
  lambda5 <- c(0.59, 0.05, 0.73, 0.0)
  lambda6 <- c(0.79, 0.9, 0.69)
  lambda7 <- c(0.37, 0.003, 0.12)
  lambda8 <- c(1.33, 0.9, 0.58)
  lambda9 <- c(6.62, 0.9, 0.18)
  
  
  k_B  <- 0 # Translocation from Olfactory to brain [d^-1]
  kb  <-  0.29 # Translocation from interstitium to blood [d^-1] 
  ko  <- 0 # Clearance from Olfactory to GI [d^-1]
  ku  <- 109 # Clearance from upper airwais to GI [d^-1]
  kt  <- 0.015 # Clearance from upper airwais to GI [d^-1]
  kr  <- 4 # Macrophage phagocytosis rate [d^-1]
  kd  <- 0.033 # Macrophage death rate [d^-1]
  ki  <- 3.5 # interstitialization rate [d^-1]
  kl  <- 0.07 # Translocation from interstitium to Lymph nodes [d^-1]
  lambdaI  <- 0.15 # fraction of translocation from blood to interstitium [unitless]
  
  params <- c(lambda3[], lambda4[], lambda5[], lambda6[], lambda7[], lambda8[], lambda9[], k_B, kb, ko, ku, kt, kr, kd, ki, kl, lambdaI, # se autes ginetai sensitivity analysis
              Vtis[],Vcap[], Q[], QTotal, Da, Do, Du, Vven, Vart)
  
  
  source("Tran_Equations.R")          # epilisi toy montelou gia ti arxikes times ton parametron
  init_solution <- Total_amount
  
  init_params <- params        # arxikes times parametrwn
  
  
  
  
  
  
  ###Epanaliptiki diadikasia ypologismou apotelesmaton gia kathe nea parametro
  for (i in 1:37) {
  params[i] <- mult*params[i]    # epivoli metavolis tis parametrou i kata mult 
  
  ###Sto parakato kommati ginetai antistoixisi twn newn timwn twn parametrwn me tis antistoixes onomasies tous
  if (i <=5){
    lambda3[i] = params[i]
  } else if ((i>=6) &  (i<=11)){
    lambda4[i-5] = params[i]
  } else if ((i>=12) & (i<=15)){
    lambda5[i-11] = params[i]
  } else if ((i>=16) & (i<=18)){
    lambda6[i-15] = params[i]
  } else if ((i>=19) & (i<=21)){
    lambda7[i-18] = params[i]
  } else if ((i>=22) & (i<=24)){
    lambda8[i-21] = params[i]
  } else if ((i>=25) & (i<= 27)){
    lambda9[i-24] = params[i]
  }
  
  k_B <- params[28]
  kb <- params[29]
  ko <- params[30]
  ku <- params[31]
  kt <- params[32]
  kr <- params[33]
  kd <- params[34]
  ki <- params[35]
  kl <- params[36]
  lambdaI <- params[37]
  
  source("Tran_Equations.R")   # epilisi toy montelou gia ti nea timi tis parametrou i
  
  ar[,,i] <- Total_amount   #dimiourgia 3d array (6x6x37), diladi ta apotelesmata twn 6 diamerismatwn se 6 xronikes stigmes epilisis 
                               #gia metavoli se kathe mia apo tis 37 parametrous
  params <- init_params 
}

###Metatropi tou pinaka ar se diastaseis 6x37x6
for (z in 1:37) {
  for (w in 1:6) {
    ar2[,z,w] <- ar[,w,z]
  }
  
}
  ###Sensitivity Index (SI) calculation
  eps <- 1e-15
  SI <- array(0, c(6, 37, 6)) 
  for (w in 1:6) {      ### deiktis diamerismatos
    for (t in 1:6) {   ### deiktis xronikis stigmis
      for (p in 1:37){
        SI[t,p,w] <- (abs(ar2[t,p,w]-init_solution[t,w]) / (Dx*init_params[p] + eps)) /(init_solution[t,w] + eps)
      }
    }	
  }

  #Dimiourgia data frame gia kathe diamerisma gia ta SI olon ton parametrwn      
  data_comp1 <- as.data.frame( cbind(sample_time, SI[,,1]))              # Total amount in lungs
  data_comp2 <- as.data.frame( cbind(sample_time, SI[,,2]))              # Total amount in brain
  data_comp3 <- as.data.frame( cbind(sample_time, SI[,,3]))              # Total amount in kidneys
  data_comp4 <- as.data.frame( cbind(sample_time, SI[,,4]))              # Total amount in liver
  data_comp5 <- as.data.frame( cbind(sample_time, SI[,,5]))              # Total amount in spleen
  data_comp6 <- as.data.frame( cbind(sample_time, SI[,,6]))              # Total amount in blood

colnames(data_comp1) <- c("Time","lambda31","lambda32","lambda33","lambda34","lambda35","lambda41","lambda42","lambda43","lambda44","lambda45","lambda46",
                          "lambda51","lambda52","lambda53","lambda54","lambda61","lambda62","lambda63","lambda71","lambda72","lambda73","lambda81","lambda82","lambda83",
                          "lambda91","lambda92","lambda93","k_B","kb","ko","ku","kt","kr","kd","ki","kl","lambdaI")
colnames(data_comp2) <- colnames(data_comp1)
colnames(data_comp3) <- colnames(data_comp1)
colnames(data_comp4) <- colnames(data_comp1)
colnames(data_comp5) <- colnames(data_comp1)
colnames(data_comp6) <- colnames(data_comp1)


ggplot(data_comp1, aes(x=Time, y=lambda31,  color = "lambda31")) +
  geom_line(size=1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda32,  color = "lambda32"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda33,  color = "lambda33"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda35,  color = "lambda35"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda41,  color = "lambda41"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda42,  color = "lambda42"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda43,  color = "lambda43"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda44,  color = "lambda44"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda45,  color = "lambda45"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda46,  color = "lambda46"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda51,  color = "lambda51"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda52,  color = "lambda52"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda53,  color = "lambda53"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda61,  color = "lambda61"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda62,  color = "lambda62"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda63,  color = "lambda63"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda71,  color = "lambda71"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda73,  color = "lambda73"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda81,  color = "lambda81"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda82,  color = "lambda82"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda83,  color = "lambda83"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda91,  color = "lambda91"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda92,  color = "lambda92"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambda93,  color = "lambda93"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=k_B,  color = "k_B"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=kb,  color = "kb"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=kl,  color = "kl"),  size = 1.2) +
  geom_line(data = data_comp1, aes(x=Time, y=lambdaI,  color = "lambdaI"),  size = 1.2) +
  
  labs(title = "SI vs Time", subtitle = "Compartment 5", y = "SI", x = "Time (in days)") +
  theme(legend.title=element_text(hjust = 0.5,size=17), 
        legend.text=element_text(size=14))
  