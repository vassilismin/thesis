library(deSolve)
dose <- 1.7e03 #micro grams
y <- c(rep(0,9),dose, rep(0,11))
inits <- as.vector(y)



  ode.func <- function(time, y, params){
    with( as.list( params),{
      
      dydt<-rep(0,21)
      
      # concentrations in tissues
      C_lu = y[1]/W_lu
      C_re = y[2]/W_re
      C_bm = y[3]/W_bm
      C_br = y[4]/W_br
      C_ht = y[5]/W_ht
      C_ki = y[6]/W_ki
      C_li = y[7]/W_li
      C_spl = y[8]/W_spl
      C_art = y[9]/(0.2*(W_blood + Wb_spl + Wb_li + Wb_lu + Wb_br + Wb_ht + Wb_ki + 
                           Wb_bm +   Wb_re))
      C_ven = y[10]/(0.8*(W_blood + Wb_spl + Wb_li + Wb_lu + Wb_br + Wb_ht + Wb_ki + 
                            Wb_bm + Wb_re))
      
      # Uptake rates by phagocytizing cells
      kluab = k_ab0*(1-(y[11]/(M_lu_cap*W_lu)))
      kreab = k_ab0*(1-(y[12]/(M_re_cap*W_re)))
      kbmab = k_ab0*(1-(y[13]/(M_bm_cap*W_bm)))
      kbrab = k_ab0*(1-(y[14]/(M_br_cap*W_br))) 
      khtab = k_ab0*(1-(y[15]/(M_ht_cap*W_ht)))
      kkiab = k_ab0*(1-(y[16]/(M_ki_cap*W_ki)))
      kliab = k_ab0*(1-(y[17]/(M_li_cap*W_li)))
      ksplab = k_ab0_spl*(1-(y[18]/(M_spl_cap*W_spl)))
      kbloodab = k_ab0*(1-(y[19]/(M_blood_cap*W_blood)))
      
      ## Nanoparticles in tissue
      # lungs
      dydt[1] = ((x_re*Q_tot)/(1+x_re))*(C_ven-C_lu/P)-(W_lu*C_lu*kluab-y[11]*k_de)
      # rest of body
      dydt[2] = ((x_re*Q_re)/(1+x_re))*(C_art-C_re/P)-(W_re*C_re*kreab-y[12]*k_de)
      # bone marrow
      dydt[3] = ((x_fast*Q_bm)/(1+x_fast))*(C_art-C_bm/P)-(W_bm*C_bm*kbmab-y[13]*k_de)
      # brain
      dydt[4] = ((x_br*Q_br)/(1+x_br))*(C_art-C_br/P)-(W_br*C_br*kbrab-y[14]*k_de)
      # heart
      dydt[5] = ((x_re*Q_ht)/(1+x_re))*(C_art-C_ht/P)-(W_ht*C_ht*khtab-y[15]*k_de)
      # kidneys
      dydt[6] = ((x_re*Q_ki)/(1+x_re))*(C_art-C_ki/P)-(W_ki*C_ki*kkiab-y[16]*k_de)-C_art*y[21]*x_fast/(1+x_re)
      # liver
      dydt[7] = C_art*x_fast*Q_li/(x_fast+1) + Q_spl*x_fast*(C_art+x_fast*C_spl/P)/((1+x_fast) * (1 +x_fast))-
        (((C_li/P)*x_fast*(Q_li+Q_spl))/(x_fast+1))-(W_li*C_li*kliab- y[17]*k_de)-y[7]*y[20]
      # spleen
      dydt[8] = ((x_fast*Q_spl)/(1+x_fast))*(C_art - C_spl/P) - (W_spl*C_spl*ksplab - y[18]*k_de)
      # arterial blood
      dydt[9] = Q_tot *((C_ven + x_re*C_lu/P)/(1+x_re)-C_art)- ((0.2*W_blood*C_art)*kbloodab-0.2*y[19]*k_de)
      # venous blood
      dydt[10] = C_li*x_fast*(Q_li+Q_spl)/(P*(1+x_fast)) +  Q_spl*x_fast*C_spl/(P*(1+x_fast)*(1+x_fast))+
        (C_br*Q_br*x_br)/(P*(1+x_br)) + (C_ht*x_re*Q_ht)/(P*(1+x_re))+ (C_ki*Q_ki*x_re)/(P*(1+x_re))+ 
        (C_bm*Q_bm*x_fast)/(P*(1+x_fast))+ (Q_re*C_re*x_re)/(P*(1+x_re))+(Q_li/(1+x_fast)+ 
                                                                            Q_spl/((1+x_fast)*(1+x_fast)) +Q_br/(1+x_br)+Q_ht/(1+x_re)+Q_ki*(1-y[21]/Q_ki)/(1+x_re)+
                                                                            Q_bm/(1+x_fast)+Q_re/(1+x_re))*C_art - Q_tot*C_ven - ((0.8*W_blood*C_ven)* kbloodab -0.8* y[19] * k_de)
      
      ## Nanoparticles uptaken in PCs
      # lungs
      dydt[11] = W_lu*C_lu*kluab - y[11]*k_de
      # rest of the body
      dydt[12] = W_re*C_re*kreab - y[12]*k_de
      # bone marrow
      dydt[13] = W_bm*C_bm*kbmab - y[13]*k_de
      # brain
      dydt[14] = W_br*C_br*kbrab - y[14]*k_de
      # heart
      dydt[15] = W_ht*C_ht*khtab - y[15]*k_de
      # kidneys
      dydt[16] = W_ki*C_ki*kkiab - y[16]*k_de 
      # liver
      dydt[17] = W_li*C_li*kliab - y[17]*k_de
      # spleen
      dydt[18] = W_spl*C_spl*ksplab - y[18]*k_de
      # blood
      dydt[19] = (0.2*W_blood*C_art + 0.8*W_blood*C_ven)*kbloodab - y[19]*k_de
      
      ##Nanoparticles in excreta (CLE_f --> y[20], CLE_u --> y[21])
      # Rate of amount in feces from liver - ug/h
      dydt[20] = y[7]*CLE_f
      # Rate of amount in urine from kidney -ug/h
      dydt[21] = C_art*CLE_u
      
  return(list(dydt))
      
    })}
  
  sample_time <- c(10/60, 1, 24, 24*7, 24*28, 24*56)
  
  solution <- ode(times = sample_time, func = ode.func, y = inits, parms = params, method = "lsodes")

  #as.vector(solution)
  