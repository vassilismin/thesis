functions{
    real [] pbpk(real t,
                real[] M,
                real[] theta,
                real[] rdata,
                int[] idata) {

    real dMdt[36] ;

    real Q_total; real V_blood; real Vven; real Vart; real Vm_ven; real Vm_art;

    real Vtis_rob; real Vcap_rob; real Vm_rob; real Q_rob;
    real Vtis_ht; real Vcap_ht; real Vm_ht; real Q_ht;
    real Vtis_ki; real Vcap_ki; real Vm_ki; real Q_ki;
    real Vtis_br; real Vcap_br; real Vm_br; real Q_br;
    real Vtis_spl; real Vcap_spl; real Vm_spl; real Q_spl;
    real Vtis_lu; real Vcap_lu; real Vm_lu; real Q_lu;
    real Vtis_li; real Vcap_li; real Vm_li; real Q_li;
    real Vtis_ut; real Vcap_ut; real Vm_ut; real Q_ut;
    real Vtis_skel; real Vcap_skel; real Vm_skel; real Q_skel;
    real Vtis_st; real Vcap_st; real Vm_st; real Q_st;

    real P_1; real P_2; real P_3; real P_4; real P_5; real P_6;
    real x_1; real x_2; real x_3; real x_4; real x_5; real x_6;
    real k_de; real Km; real Pup; real CLE_f; real CLE_u;

    real P_lu; real P_spl; real P_li; real P_ki; real P_ht;
    real P_br; real P_ut; real P_skel; real P_st; real P_rob;

    real Ccap_lu; real Ctis_lu; real Cm_lu;
    real Ccap_spl; real Ctis_spl; real Cm_spl;
    real Ccap_li; real Ctis_li; real Cm_li;
    real Ccap_ki; real Ctis_ki; real Cm_ki;
    real Ccap_ht; real Ctis_ht; real Cm_ht;
    real Ccap_br; real Ctis_br; real Cm_br;
    real Ccap_ut; real Ctis_ut; real Cm_ut;
    real Ccap_skel; real Ctis_skel; real Cm_skel;
    real Ccap_st; real Ctis_st; real Cm_st;
    real Ccap_rob; real Ctis_rob; real Cm_rob;

    real C_art; real Cm_art;
    real C_ven; real Cm_ven;

    real PA_lu; real PA_spl; real PA_li; real PA_ki; real PA_ht;
    real PA_br; real PA_ut; real PA_skel; real PA_st; real PA_rob;

    real Pup_lu; real Pup_spl; real Pup_li; real Pup_ki; real Pup_ht;
    real Pup_br; real Pup_ut; real Pup_skel; real Pup_st; real Pup_rob;
    real Pup_ven; real Pup_art;

    Vcap_lu = rdata[1];
    Vcap_spl = rdata[2];
    Vcap_li = rdata[3];
    Vcap_ki = rdata[4];
    Vcap_ht = rdata[5];
    Vcap_br = rdata[6];
    Vcap_ut = rdata[7];
    Vcap_skel = rdata[8];
    Vcap_st = rdata[9];
    Vcap_rob = rdata[10];

    Vtis_lu = rdata[11];
    Vtis_spl = rdata[12];
    Vtis_li = rdata[13];
    Vtis_ki = rdata[14];
    Vtis_ht = rdata[15];
    Vtis_br = rdata[16];
    Vtis_ut = rdata[17];
    Vtis_skel = rdata[18];
    Vtis_st = rdata[19];
    Vtis_rob = rdata[20];

    Vm_lu = rdata[21];
    Vm_spl = rdata[22];
    Vm_li = rdata[23];
    Vm_ki = rdata[24];
    Vm_ht = rdata[25];
    Vm_br = rdata[26];
    Vm_ut = rdata[27];
    Vm_skel = rdata[28];
    Vm_st = rdata[29];
    Vm_rob = rdata[30];

    Q_lu = rdata[31];
    Q_spl = rdata[32];
    Q_li = rdata[33];
    Q_ki = rdata[34];
    Q_ht = rdata[35];
    Q_br = rdata[36];
    Q_ut = rdata[37];
    Q_skel = rdata[38];
    Q_st = rdata[39];
    Q_rob = rdata[40];

    Vven = rdata[41];
    Vm_ven = rdata[42];
    Vart = rdata[43];
    Vm_art = rdata[44];
    Q_total = rdata[45];

    ///////////params///////////
    P_1 = exp(theta[1]);
    P_2 = exp(theta[2]);
    P_3 = exp(theta[3]);
    P_4 = exp(theta[4]);
    P_5 = exp(theta[5]);
    P_6 = exp(theta[6]);
    x_1 = exp(theta[7]);
    x_2 = exp(theta[8]);
    x_3 = exp(theta[9]);
    x_4 = exp(theta[10]);
    x_5 = exp(theta[11]);
    x_6 = exp(theta[12]);
    k_de = exp(theta[13]);
    Km = exp(theta[14]);
    Pup = exp(theta[15]);
    CLE_f = exp(theta[16]);
    CLE_u = exp(theta[17]);

    /////////////////
    // ODEs system //
    /////////////////

   //Capillary concentrations
   Ccap_lu = M[1]/Vcap_lu;
   Ccap_spl = M[4]/Vcap_spl;
   Ccap_li = M[7]/Vcap_li;
   Ccap_ki = M[10]/Vcap_ki;
   Ccap_ht = M[13]/Vcap_ht;
   Ccap_br = M[16]/Vcap_br;
   Ccap_ut = M[19]/Vcap_ut;
   Ccap_skel = M[22]/Vcap_skel;
   Ccap_st = M[25]/Vcap_st;
   Ccap_rob = M[28]/Vcap_rob;

   //Tissue concentration
   Ctis_lu = M[2]/Vtis_lu;
   Ctis_spl = M[5]/Vtis_spl;
   Ctis_li = M[8]/Vtis_li;
   Ctis_ki = M[11]/Vtis_ki;
   Ctis_ht = M[14]/Vtis_ht;
   Ctis_br = M[17]/Vtis_br;
   Ctis_ut = M[20]/Vtis_ut;
   Ctis_skel = M[23]/Vtis_skel;
   Ctis_st = M[26]/Vtis_st;
   Ctis_rob = M[29]/Vtis_rob;

   //Concentration in Macrophages cells sub-compartment
   Cm_lu = M[3]/Vm_lu;
   Cm_spl = M[6]/Vm_spl;
   Cm_li = M[9]/Vm_li;
   Cm_ki = M[12]/Vm_ki;
   Cm_ht = M[15]/Vm_ht;
   Cm_br = M[18]/Vm_br;
   Cm_ut = M[21]/Vm_ut;
   Cm_skel = M[24]/Vm_skel;
   Cm_st = M[27]/Vm_st;
   Cm_rob = M[30]/Vm_rob;

   Cm_art = M[34]/Vm_art;
   Cm_ven = M[32]/Vm_ven;

   C_ven = M[31]/Vven;
   C_art = M[33]/Vart;


    PA_lu=x_3*Q_total;
    PA_spl=x_2*Q_spl;
    PA_li=x_1*Q_li;
    PA_ki=x_3*Q_ki;
    PA_ht=x_4*Q_ht;
    PA_br=x_6*Q_br;
    PA_ut=x_4*Q_ut;
    PA_skel=x_5*Q_skel;
    PA_st=x_5*Q_st;
    PA_rob=x_4*Q_rob;

    P_lu = P_3;
    P_spl = P_2;
    P_li = P_1;
    P_ki = P_3;
    P_ht = P_4;
    P_br = P_6;
    P_ut = P_4;
    P_skel = P_5;
    P_st = P_5;
    P_rob = P_6;

   Pup_lu=Pup*Vm_lu*(1-(Cm_lu/(Km+Cm_lu))); //ml/h
   Pup_spl=Pup*Vm_spl*(1-(Cm_spl/(Km+Cm_spl))); //ml/h
   Pup_li=Pup*Vm_li*(1-(Cm_li/(Km+Cm_li))); //ml/h
   Pup_ki=Pup*Vm_ki*(1-(Cm_ki/(Km+Cm_ki))); //ml/h
   Pup_ht=Pup*Vm_ht*(1-(Cm_ht/(Km+Cm_ht))); //ml/h
   Pup_br=Pup*Vm_br*(1-(Cm_br/(Km+Cm_br))); //ml/h
   Pup_ut=Pup*Vm_ut*(1-(Cm_ut/(Km+Cm_ut))); //ml/h
   Pup_skel=Pup*Vm_skel*(1-(Cm_skel/(Km+Cm_skel))); //ml/h
   Pup_st=Pup*Vm_st*(1-(Cm_st/(Km+Cm_st))); //ml/h
   Pup_rob=Pup*Vm_rob*(1-(Cm_rob/(Km+Cm_rob))); //ml/h

   Pup_ven=Pup*Vm_ven*(1-(Cm_ven/(Km+Cm_ven))); //ml/h
   Pup_art=Pup*Vm_art*(1-(Cm_art/(Km+Cm_art))); //ml/h

   //Lungs
   dMdt[1] = Q_total*C_ven - Q_total*Ccap_lu - PA_lu*Ccap_lu + PA_lu*Ctis_lu/P_lu;
   dMdt[2] = PA_lu*Ccap_lu - PA_lu*Ctis_lu/P_lu - Pup_lu*Ctis_lu + k_de*M[3];
   dMdt[3]   = Pup_lu*Ctis_lu - k_de*M[3];

   //Spleen
   dMdt[4] = Q_spl*C_art - Q_spl*Ccap_spl - PA_spl*Ccap_spl + PA_spl*Ctis_spl/P_spl;
   dMdt[5] = PA_spl*Ccap_spl - PA_spl*Ctis_spl/P_spl - Pup_spl*Ctis_spl + k_de*M[6];
   dMdt[6]   = Pup_spl*Ctis_spl - k_de*M[6];

   //Liver
   dMdt[7] = Q_li*C_art + Q_spl*Ccap_spl - (Q_li+Q_spl)*Ccap_li - PA_li*Ccap_li + PA_li*Ctis_li/P_li;
   dMdt[8] = PA_li*Ccap_li - PA_li*Ctis_li/P_li - Pup_li*Ctis_li + k_de*M[9] - CLE_f*M[8];
   dMdt[9]   = Pup_li*Ctis_li - k_de*M[9];
   dMdt[35] = CLE_f*M[8];

   //Kidneys
   dMdt[10] = Q_ki*C_art - Q_ki*Ccap_ki - PA_ki*Ccap_ki + PA_ki*Ctis_ki/P_ki;
   dMdt[11] = PA_ki*Ccap_ki - PA_ki*Ctis_ki/P_ki - Pup_ki*Ctis_ki + k_de*M[12] - CLE_u*M[11];
   dMdt[12]   = Pup_ki*Ctis_ki - k_de*M[12];
   dMdt[36] = CLE_u*M[11];

   //Heart
   dMdt[13] = Q_ht*C_art - Q_ht*Ccap_ht - PA_ht*Ccap_ht + PA_ht*Ctis_ht/P_ht;
   dMdt[14] = PA_ht*Ccap_ht - PA_ht*Ctis_ht/P_ht - Pup_ht*Ctis_ht + k_de*M[15];
   dMdt[15] = Pup_ht*Ctis_ht - k_de*M[15];

   //Brain
   dMdt[16] = Q_br*C_art - Q_br*Ccap_br - PA_br*Ccap_br + PA_br*Ctis_br/P_br;
   dMdt[17] = PA_br*Ccap_br - PA_br*Ctis_br/P_br - Pup_br*Ctis_br + k_de*M[18];
   dMdt[18] = Pup_br*Ctis_br - k_de*M[18];

   //Uterus
   dMdt[19] = Q_ut*C_art - Q_ut*Ccap_ut - PA_ut*Ccap_ut + PA_ut*Ctis_ut/P_ut;
   dMdt[20] = PA_ut*Ccap_ut - PA_ut*Ctis_ut/P_ut - Pup_ut*Ctis_ut + k_de*M[21];
   dMdt[21] = Pup_ut*Ctis_ut - k_de*M[21];

   //Skeleton
   dMdt[22] = Q_skel*C_art - Q_skel*Ccap_skel - PA_skel*Ccap_skel + PA_skel*Ctis_skel/P_skel;
   dMdt[23] = PA_skel*Ccap_skel - PA_skel*Ctis_skel/P_skel - Pup_skel*Ctis_skel + k_de*M[24];
   dMdt[24] = Pup_skel*Ctis_skel - k_de*M[24];

   //Soft tissue
   dMdt[25] = Q_st*C_art - Q_st*Ccap_st - PA_st*Ccap_st + PA_st*Ctis_st/P_st;
   dMdt[26] = PA_st*Ccap_st - PA_st*Ctis_st/P_st - Pup_st*Ctis_st + k_de*M[27];
   dMdt[27] = Pup_st*Ctis_st - k_de*M[27];

   //Rest of Body
   dMdt[28] = Q_rob*C_art - Q_rob*Ccap_rob - PA_rob*Ccap_rob + PA_rob*Ctis_rob/P_rob;
   dMdt[29] = PA_rob*Ccap_rob - PA_rob*Ctis_rob/P_rob - Pup_rob*Ctis_rob + k_de*M[30];
   dMdt[30] = Pup_rob*Ctis_rob - k_de*M[30];

   //Veins
   dMdt[31] = - Q_total*C_ven + (Q_li+Q_spl)*Ccap_li + Q_ki*Ccap_ki + Q_ht*Ccap_ht + Q_br*Ccap_br + Q_ut*Ccap_ut + Q_skel*Ccap_skel +
             Q_st*Ccap_st + Q_rob*Ccap_rob - Pup_ven*C_ven + k_de*M[32];
   dMdt[32] = Pup_ven*C_ven - k_de*M[32];

   //Arteries
   dMdt[33] = Q_total*Ccap_lu - Q_spl*C_art - Q_li*C_art - Q_ki*C_art - Q_ht*C_art - Q_br*C_art - Q_ut*C_art - Q_skel*C_art - Q_st*C_art - Q_rob*C_art -
             Pup_art*C_art + k_de*M[34];
   dMdt[34] =Pup_art*C_art - k_de*M[34];

   return dMdt;
   }
}

//////////////////////////////////////////////////////////////////////////

data{

        int<lower=0> N_param;                // Number of parameters to be estimated
        int<lower=0> N_compart;              //number of observed compartments
        int<lower=0> N_diff;                 // number of differential equations
        real  time[5];
        real  urine_time[11];
        real  feces_time[4];
        real  mass[5,N_compart];
        real  feces[4];
        real  urine[11];
        real  m0[N_diff];           // Initial concentration in compartments
        real  t_init;                  // Initial time
        real  eta_tr[N_param];
        real  params[56];      // Matrix containing the individual parameters
        real  rel_tol;
        real  abs_tol;
        real  max_num_steps;

}
//////////////////////////////////////////////////////////////////
transformed data {
      real rdata[0];
      int idata[0];
      vector[N_param]  eta_tr_std ;
      vector[N_param]  eta_std ;
      vector[N_param]  eta ;
      vector [N_param] H;                //covariance matrix



      for (i in 1:N_param){
              eta_tr_std[i] = eta_tr[i];
              eta_std[i]= sqrt(log(((eta_tr_std[i]^2)/(eta_tr[i])^2)+1));
              eta[i]=log(((eta_tr[i])^2)/sqrt((eta_tr_std[i]^2)+(eta_tr[i])^2));
              H[i] = eta_std[i];
      }
}
//////////////////////////////////////////////////////////////////



parameters{

        real<lower=0>  sigma1;
        real<lower=0>  sigma2;
        real<lower=0>  sigma3;

            vector [N_param] theta_tr;
            }

////////////////////////////////////////////////////////////////////
model{
real m_hat[5,N_diff];
real total_m_hat[5,N_compart];
real feces_excr[4,N_diff];
real urine_excr[11,N_diff];
real feces_hat[4];
real urine_hat[11];



//priors
sigma1 ~ normal(0,1);
sigma2 ~ normal(0,1);
sigma3 ~ normal(0,1);


theta_tr[:] ~normal(eta[:],H[:]);

//likelihood

m_hat[:,:] = integrate_ode_bdf(pbpk, m0, t_init, time,
            to_array_1d(theta_tr[:]), params[:], idata,
            rel_tol, abs_tol, max_num_steps);

feces_excr[:,:] = integrate_ode_bdf(pbpk, m0, t_init, feces_time,
            to_array_1d(theta_tr[:]), params[:], idata,
            rel_tol, abs_tol, max_num_steps);

urine_excr[:,:] = integrate_ode_bdf(pbpk, m0, t_init, urine_time,
            to_array_1d(theta_tr[:]), params[:], idata,
            rel_tol, abs_tol, max_num_steps);

for (i in 1:5){
          //Total amount of NPs in each organ

          //Amount in Liver
          total_m_hat[i,1] = m_hat[i,7] + m_hat[i,8] + m_hat[i,9];

          //Amount in Spleen
          total_m_hat[i,2] = m_hat[i,4] + m_hat[i,5] + m_hat[i,6];

          //Amount in Kidneys
          total_m_hat[i,3] = m_hat[i,10] + m_hat[i,11] + m_hat[i,12];

          //Ammount in Lungs
          total_m_hat[i,4] = m_hat[i,1] + m_hat[i,2] + m_hat[i,3];

          //Amount in Heart
          total_m_hat[i,5] = m_hat[i,13] + m_hat[i,14] + m_hat[i,15];

          //Amount in Brain
          total_m_hat[i,6] = m_hat[i,16] + m_hat[i,17] + m_hat[i,18];

          //Amount in Uterus
          total_m_hat[i,7] = m_hat[i,19] + m_hat[i,20] + m_hat[i,21];

          //Amount in Blood
          total_m_hat[i,8] = m_hat[i,31] + m_hat[i,32] + m_hat[i,33] + m_hat[i,34];

          //Amount in Skeleton
          total_m_hat[i,9] = m_hat[i,22] + m_hat[i,23] + m_hat[i,24];

          //Amount in Soft Tissue
          total_m_hat[i,10] = m_hat[i,25] + m_hat[i,26] + m_hat[i,27];
          }

for (j in 1:N_compart){
  if(j==1){
    to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma1);
  }
  else if(j==2){
    to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma1);
  }
  else if(j==3){
    to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma2);
  }
  else if(j==4){
    to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma2);
  }
  else if(j==5){
    to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma3);
  }
  else if(j==6){
    to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma3);
  }
  else if(j==7){
    to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma3);
  }
  else if(j==8){
    to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma2);
  }
  else if(j==9){
    to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma1);
  }
  else if(j==10){
    to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma2);
  }
}

feces_hat[:] = feces_excr[:,35];
urine_hat[:] = urine_excr[:,36];
to_vector(feces[:]) ~ normal(to_vector(feces_hat[:]),sigma2);
to_vector(urine[:]) ~ normal(to_vector(urine_hat[:]),sigma2);
}

generated quantities{
  vector [N_param] theta;

  theta = exp(theta_tr);
}
