functions{
 
        real [] pbpk(real t,
                    real[] y,
                    real[] theta,
                    real[] rdata,
                    int[] idata) {
                
              real dydt[21] ;

              real C_lu; real C_re; real C_bm; real C_br; real C_ht; real C_ki; real C_li; real C_spl;            
              real C_art; real C_ven;
              
              real W_lu; real W_re; real W_bm; real W_br; real W_ht; real W_ki; real W_li; real W_spl;real W_blood;            
              
              real Wb_lu; real Wb_re; real Wb_bm; real Wb_br; real Wb_ht; real Wb_ki; real Wb_li; real Wb_spl;   
              
              real Q_tot; real Q_bm; real Q_br; real Q_ht; real Q_ki; real Q_spl; real Q_li; real Q_re; 
 
              real CLE_u; real x_br; 
              
              real M_lu_cap; real M_bm_cap; real M_br_cap;real M_ht_cap;
              real M_ki_cap; real M_li_cap; real M_spl_cap; real M_blood_cap;
              real M_re_cap; real x_fast; real x_re; real P;
              real k_ab0; real k_ab0_spl; real k_de; real CLE_f;

              real kluab; real kreab; real kbmab; real kbrab; real khtab; real kkiab; real kliab; 
              real ksplab; real kbloodab; 

              
              // tissue weights (in g)
              W_lu = rdata[1];
              W_bm = rdata[2];
              W_br = rdata[3];
              W_ht = rdata[4];
              W_ki = rdata[5];
              W_li = rdata[6];
              W_spl =  rdata[7];
              W_re =  rdata[8] ;

              // Weight of capillary blood assuming density = 1
              Wb_lu = rdata[9]; 
              Wb_bm = rdata[10]; 
              Wb_br = rdata[11];
              Wb_ht = rdata[12]; 
              Wb_ki = rdata[13]; 
              Wb_li = rdata[14]; 
              Wb_spl = rdata[15]; 
              Wb_re = rdata[16]; 
              W_blood = rdata[17]; 

              //Regional blood flows (in mL per hour)

             Q_tot = rdata[18]; 
             Q_bm = rdata[19]; 
             Q_br = rdata[20]; 
             Q_ht = rdata[21]; 
             Q_ki = rdata[22]; 
             Q_spl = rdata[23]; 	
             Q_li = rdata[24]; 
             Q_re = rdata[25]; 

             CLE_u = rdata[26]; 
             x_br = rdata[27]; 
             
            ///////////params///////////
            M_lu_cap = exp(theta[1]);
            M_bm_cap = exp(theta[2]);
            M_br_cap = exp(theta[3]);
            M_ht_cap = exp(theta[4]);
            M_ki_cap = exp(theta[5]);
            M_li_cap = exp(theta[6]);
            M_spl_cap = exp(theta[7]);
            M_blood_cap = exp(theta[8]);
            M_re_cap = exp(theta[9]);

            x_fast = exp(theta[10]);
            x_re = exp(theta[11]);
            P = exp(theta[12]);
            k_ab0 = exp(theta[13]);
            k_ab0_spl = exp(theta[14]);
            k_de = exp(theta[15]);
            CLE_f = exp(theta[16]);


              // concentrations in tissues
              C_lu = y[1]/W_lu;
              C_re = y[2]/W_re;
              C_bm = y[3]/W_bm;
              C_br = y[4]/W_br;
              C_ht = y[5]/W_ht;
              C_ki = y[6]/W_ki;
              C_li = y[7]/W_li;
              C_spl = y[8]/W_spl;
              C_art = y[9]/(0.2*(W_blood + Wb_spl + Wb_li + Wb_lu + Wb_br + Wb_ht + Wb_ki + 
              Wb_bm +   Wb_re));
              C_ven = y[10]/(0.8*(W_blood + Wb_spl + Wb_li + Wb_lu + Wb_br + Wb_ht + Wb_ki + 
              Wb_bm + Wb_re));

              // Uptake rates by phagocytizing cells
              kluab = k_ab0*(1-(y[11]/(M_lu_cap*W_lu)));
              kreab = k_ab0*(1-(y[12]/(M_re_cap*W_re)));
              kbmab = k_ab0*(1-(y[13]/(M_bm_cap*W_bm)));
              kbrab = k_ab0*(1-(y[14]/(M_br_cap*W_br))); 
              khtab = k_ab0*(1-(y[15]/(M_ht_cap*W_ht)));
              kkiab = k_ab0*(1-(y[16]/(M_ki_cap*W_ki)));
              kliab = k_ab0*(1-(y[17]/(M_li_cap*W_li)));
              ksplab = k_ab0_spl*(1-(y[18]/(M_spl_cap*W_spl)));
              kbloodab = k_ab0*(1-(y[19]/(M_blood_cap*W_blood)));

             // Nanoparticles in tissue
             // lungs
             dydt[1] = ((x_re*Q_tot)/(1+x_re))*(C_ven-C_lu/P)-(W_lu*C_lu*kluab-y[11]*k_de);
             // rest of body
             dydt[2] = ((x_re*Q_re)/(1+x_re))*(C_art-C_re/P)-(W_re*C_re*kreab-y[12]*k_de);
             // bone marrow
             dydt[3] = ((x_fast*Q_bm)/(1+x_fast))*(C_art-C_bm/P)-(W_bm*C_bm*kbmab-y[13]*k_de);
             // brain
             dydt[4] = ((x_br*Q_br)/(1+x_br))*(C_art-C_br/P)-(W_br*C_br*kbrab-y[14]*k_de);
             // heart
             dydt[5] = ((x_re*Q_ht)/(1+x_re))*(C_art-C_ht/P)-(W_ht*C_ht*khtab-y[15]*k_de);
             // kidneys
             dydt[6] = ((x_re*Q_ki)/(1+x_re))*(C_art-C_ki/P)-(W_ki*C_ki*kkiab-y[16]*k_de)-C_art*CLE_u*x_fast/(1+x_re) ;
             // liver
             dydt[7] = C_art*x_fast*Q_li/(x_fast+1) + Q_spl*x_fast*(C_art+x_fast*C_spl/P)/((1+x_fast) * (1 +x_fast))-(((C_li/P)*x_fast*(Q_li+Q_spl))/(x_fast+1))-(W_li*C_li*kliab- y[17]*k_de)-y[7]*CLE_f;
             // spleen
             dydt[8] = ((x_fast*Q_spl)/(1+x_fast))*(C_art - C_spl/P) - (W_spl*C_spl*ksplab - y[18]*k_de);
             // arterial blood
             dydt[9] = Q_tot *((C_ven + x_re*C_lu/P)/(1+x_re)-C_art)- ((0.2*W_blood*C_art)*kbloodab-0.2*y[19]*k_de);
             // venous blood
             dydt[10] = C_li*x_fast*(Q_li+Q_spl)/(P*(1+x_fast)) +  Q_spl*x_fast*C_spl/(P*(1+x_fast)*(1+x_fast))+ (C_br*Q_br*x_br)/(P*(1+x_br)) + (C_ht*x_re*Q_ht)/(P*(1+x_re))+ (C_ki*Q_ki*x_re)/(P*(1+x_re))+ (C_bm*Q_bm*x_fast)/(P*(1+x_fast))+ (Q_re*C_re*x_re)/(P*(1+x_re))+(Q_li/(1+x_fast) +Q_spl/((1+x_fast)*(1+x_fast)) +Q_br/(1+x_br)+Q_ht/(1+x_re)+Q_ki*(1-CLE_u/Q_ki)/(1+x_re)+ Q_bm/(1+x_fast)+Q_re/(1+x_re))*C_art - Q_tot*C_ven - ((0.8*W_blood*C_ven)* kbloodab -0.8* y[19] * k_de);

             // Nanoparticles uptaken in PCs
             // lungs
             dydt[11] = W_lu*C_lu*kluab - y[11]*k_de;
             // rest of the body
             dydt[12] = W_re*C_re*kreab - y[12]*k_de;
             // bone marrow
             dydt[13] = W_bm*C_bm*kbmab - y[13]*k_de;
             // brain
             dydt[14] = W_br*C_br*kbrab - y[14]*k_de;
             // heart
             dydt[15] = W_ht*C_ht*khtab - y[15]*k_de;
             // kidneys
             dydt[16] = W_ki*C_ki*kkiab - y[16]*k_de; 
             // liver
             dydt[17] = W_li*C_li*kliab - y[17]*k_de;
             // spleen
             dydt[18] = W_spl*C_spl*ksplab - y[18]*k_de;
             // blood
             dydt[19] = (0.2*W_blood*C_art + 0.8*W_blood*C_ven)*kbloodab - y[19]*k_de;

             //Nanoparticles in excreta (CLE_f --> y[20], CLE_u --> y[21])
             // Rate of amount in feces from liver - ug/h
             dydt[20] = y[7]*CLE_f;
             // Rate of amount in urine from kidney -ug/h
             dydt[21] = C_art*CLE_u;
            
             return dydt;       

        } 
}
//////////////////////////////////////////////////////////////////////////
        
data{
                
        int<lower=0> N_param;                 // Number of parameters to be estimated
        int<lower=0> N_compart;               //number of observed compartments
        int<lower=0> N_diff;              // number of differential equations
        int<lower=0>  N_obs;                 // Total number of observations
        real  time[6];
        real  mass[6,6];
        int   samp[N_compart];                   // Number of samples of each compartment
        real m0[N_diff];           // Initial concentration in compartments
        real t_init;                  // Initial time
        real eta_tr[N_param];
        real  params[27];      // Matrix containing the individual parameters
        real rel_tol;
        real abs_tol;
        real max_num_steps;

}

////////////////////////////////////////////////////////////////////////////
        
transformed data {
                
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
        
 
          
   //print(solution);
   //print(mass);
   //print(time)

}
//////////////////////////////////////////////////////////////////
        
parameters{
                
    real<lower=0>  sigma1;
    real<lower=0>  sigma2;
    vector [N_param] theta_tr;

                
}
////////////////////////////////////////////////////
transformed parameters{
        
     
}
////////////////////////////////////////////////////////////////////
        
model{

real m_hat[6,21];
real total_m_hat[6,N_compart];
//svise apo edw
real Wb_ki; real Wb_br; real Wb_lu; real Wb_li; real W_blood;real Wb_spl ;real Wb_ht; real Wb_bm ; real Wb_re;
real W_lu; real W_br;real W_ki; real W_li; real W_spl;
real C_art[6]; real C_ven[6]; real C_ki[6]; real C_br[6];real C_spl[6];real C_lu[6];real C_li[6];
real x_re ; real x_fast;  real P; real fro;real frbr;  
//mexri edw
int pos; // uxiliary variable for segmentation
pos = 1;


//priors
 sigma1 ~ normal(0,10);
 sigma2 ~ normal(0,50);

 theta_tr[:] ~normal(eta[:],H[:]);
 
 //svise apo edw
 x_re = exp(theta_tr[11]);
 x_fast = exp(theta_tr[10]);
 P = exp(theta_tr[12]);
 fro =  exp(theta_tr[18]);
 frbr =  exp(theta_tr[17]);
 W_lu = params[1];
 W_br = params[3];
 W_ki = params[5];
 W_li = params[6];
 W_spl =  params[7];
 Wb_lu = params[9]; 
 Wb_bm = params[10]; 
 Wb_br = params[11];
 Wb_ht = params[12]; 
 Wb_ki = params[13]; 
 Wb_li = params[14]; 
 Wb_spl = params[15]; 
 Wb_re = params[16]; 
 W_blood = params[17];
 //mexri edw
 
 //likelihood~  

       m_hat[:,:] = integrate_ode_bdf(pbpk,m0,t_init, time,
                                                 to_array_1d(theta_tr[:]),params[:],idata,
                                                 rel_tol, abs_tol,max_num_steps);
       for (i in 1:6){
           // svise edw
          C_lu[i] = m_hat[i,1]/W_lu;
          C_br[i] = m_hat[i,4]/W_br;
          C_ki[i] = m_hat[i,6]/W_ki;
          C_li[i] = m_hat[i,7]/W_li;
          C_spl[i] = m_hat[i,8]/W_spl;
          C_art[i] = m_hat[i,9]/(0.2*(W_blood + Wb_spl + Wb_li + Wb_lu + Wb_br + Wb_ht + Wb_ki + 
              Wb_bm +   Wb_re));
          C_ven[i] =m_hat[i,10]/(0.8*(W_blood + Wb_spl + Wb_li + Wb_lu + Wb_br + Wb_ht + Wb_ki + 
              Wb_bm + Wb_re));
        //mexri edw
                                      
       // Total amount of NPs in each organ
       // Amount in kidneys
       
       //svise apo to fro kai meta 
       total_m_hat[i,1] = m_hat[i,6] + m_hat[i,16] +  fro*(C_art[i]+x_re*C_ki[i]/P)/(1+x_re)*Wb_ki;
       // Amount in brain
       total_m_hat[i,2] = m_hat[i,4] + m_hat[i,14] + frbr*(C_art[i]+0*C_br[i]/P)/(1+0)*Wb_br;
       // Amount in spleen
       total_m_hat[i,3] = m_hat[i,8] + m_hat[i,18]+  fro*(C_art[i]+x_fast*C_spl[i]/P)/(1+x_fast)*Wb_spl;
       // Amount in lungs
       total_m_hat[i,4] = m_hat[i,1] + m_hat[i,11] +  fro*(C_ven[i]+x_re*C_lu[i]/P)/(1+x_re)*Wb_lu ;
       // Amount in liver
       total_m_hat[i,5] = m_hat[i,7] + m_hat[i,17]+  fro*(C_art[i]+x_fast*C_li[i]/P)/(1+x_fast)*Wb_li;
       // Amount in blood
       total_m_hat[i,6] = m_hat[i,9]+m_hat[i,10]+m_hat[i,19];
       }
    

       for (j in 1:N_compart){
         if (j == 1 || j==2 || j==4 || j==6){
           //log(mass[i,j]) ~ normal(log(total_m_hat[i,j]+1e-15),sigma);
             to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma1);
         } else{
             to_vector(mass[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma2);

         }
         
       }

}


generated quantities{
      vector [N_param] theta;
      theta = exp(theta_tr);


}
