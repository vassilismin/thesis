functions{
        real []  pbpk(real t,
                      real[] M,
                      real[] theta,
                      real[] rdata,
                      int[] idata) {
        real dMdt[32] ;

        real lambda3[5]; real lambda4[6]; real lambda5[4]; real lambda6[3];
        real lambda7[3]; real lambda8[3]; real lambda9[3];

        real Vtis[7]; real Vven; real Vart; real Vcap[7];

        real Q[7]; real Q_liver; real QTotal; real Q_bile;
        real Da; real Do; real Du;

        real k_B; real kb; real ko; real ku; real kt; real kr;
        real kd; real ki; real kl;

        real lambdaI;

        Vtis[1] = rdata[1];
        Vtis[2] = rdata[2];
        Vtis[3] = rdata[3];
        Vtis[4] = rdata[4];
        Vtis[5] = rdata[5];
        Vtis[6] = rdata[6];
        Vtis[7] = rdata[7];
        Vven = rdata[8];
        Vart = rdata[9];
        Vcap[1] = rdata[10];
        Vcap[2] = rdata[11];
        Vcap[3] = rdata[12];
        Vcap[4] = rdata[13];
        Vcap[5] = rdata[14];
        Vcap[6] = rdata[15];
        Vcap[7] = rdata[16];
        Q[1] = rdata[17];
        Q[2] = rdata[18];
        Q[3] = rdata[19];
        Q[4] = rdata[20];
        Q[5] = rdata[21];
        Q[6] = rdata[22];
        Q[7] = rdata[23];
        Q_liver = rdata[24];
        QTotal = rdata[25];
        Q_bile = rdata[26];
        Da = rdata[27];
        Do = rdata[28];
        Du = rdata[29];


        ///////////params///////////
        lambda3[1] = exp(theta[1]);
        lambda3[2] = exp(theta[2]);
        lambda3[3] = exp(theta[3]);
        lambda3[4] = exp(theta[4]);
        lambda3[5] = exp(theta[5]);

        lambda4[1] = exp(theta[6]);
        lambda4[2] = exp(theta[7]);
        lambda4[3] = exp(theta[8]);
        lambda4[4] = exp(theta[9]);
        lambda4[5] = exp(theta[10]);
        lambda4[6] = exp(theta[11]);

        lambda5[1] = exp(theta[12]);
        lambda5[2] = exp(theta[13]);
        lambda5[3] = exp(theta[14]);
        lambda5[4] = exp(theta[15]);

        lambda6[1] = exp(theta[16]);
        lambda6[2] = exp(theta[17]);
        lambda6[3] = exp(theta[18]);

        lambda7[1] = exp(theta[19]);
        lambda7[2] = exp(theta[20]);
        lambda7[3] = exp(theta[21]);

        lambda8[1] = exp(theta[22]);
        lambda8[2] = exp(theta[23]);
        lambda8[3] = exp(theta[24]);

        lambda9[1] = exp(theta[25]);
        lambda9[2] = exp(theta[26]);
        lambda9[3] = exp(theta[27]);

        k_B = exp(theta[28]);
        kb = exp(theta[29]);
        ko = exp(theta[30]);
        ku = exp(theta[31]);
        kt = exp(theta[32]);
        kr = exp(theta[33]);
        kd = exp(theta[34]);
        ki = exp(theta[35]);
        kl = exp(theta[36]);
        lambdaI = exp(theta[37]);

        ///////////////
        // ODEs system //
        ///////////////

        dMdt[1]  = Do - (ko*M[1]) - (k_B*M[1]) ;   // Olfactory
        dMdt[2]  = Du - (ku*M[2])  ;             // Upper airways
        dMdt[3]  = Da - (kr*M[3])  + (kd*M[4]) - (ki*M[3]) ;  // Alveolar free
        dMdt[4]  = (kr*M[3]) - (kd*M[4]) - (kt*M[4]) ;          // Alveolar Mac
        dMdt[5]  = (ki*M[3]) - (kl*M[5]) - (kb*M[5]) + ((M[30]/Vven)*lambdaI*QTotal);    // I
        dMdt[6]  = (kl*M[5]) ;  // Lymph       to arterial  from ven

        //Liver
        dMdt[7]  = ((M[29]/Vart)*Q[1]) - ((M[7]/Vcap[1])*Q[1]) + ((M[8]/Vtis[1])*lambda3[1]*Q[1]) - ((M[7]/Vcap[1])*Q_liver*lambda3[2]) -
                 ((M[7]/Vcap[1])*Q_bile*lambda3[5]) + ((M[10]/Vcap[2])*Q[2]*lambda4[5]) + ((M[20]/Vcap[5])*Q[5]) ;  //capillary
        dMdt[8] =  ((M[7]/Vcap[1])*Q[1]) - ((M[8]/Vtis[1])*lambda3[1]*Q[1]) - (lambda3[3]*M[8])  ; //tissue
        dMdt[9]  = (lambda3[3]*M[8]) ; //seq

        //GI
        dMdt[10]  = ((M[29]/Vart)*Q[2]) - ((M[10]/Vcap[2])*Q[2]) + ((M[11]/Vtis[2])*lambda4[1]*Q[2])  - ((M[10]/Vcap[2])*Q[2]*lambda4[5])  ;  //capillary
        dMdt[11]  = ((M[10]/Vcap[2])*Q[2]) - ((M[11]/Vtis[2])*lambda4[1]*Q[2]) - (lambda4[3]*M[11]) + (lambda4[6]*M[13]) ; //tissue
        dMdt[12]  = (lambda4[3]*M[11]); //seq
        dMdt[13]  = (ko*M[1]) + (ku*M[2]) + (kt*M[4]) - (lambda4[6]*M[13]) - (lambda4[4]*M[13])+ ((M[7]/Vcap[1])*Q_bile*lambda3[5]) ; // A4

        //Kidney
        dMdt[14]  = ((M[29]/Vart)*Q[3]) - ((M[14]/Vcap[3])*Q[3]) + ((M[15]/Vtis[3])*lambda5[1]*Q[3]) - ((M[14]/Vcap[3])*Q[3]*lambda5[2]) -
                  (M[14]*lambda5[4]);  //capillary
        dMdt[15] =  ((M[14]/Vcap[3])*Q[3]) - ((M[15]/Vtis[3])*lambda5[1]*Q[3]) - (lambda5[3]*M[15])  ; //tissue
        dMdt[16]  = (lambda5[3]*M[15]) ; //seq

        //Heart
        dMdt[17]  = ((M[29]/Vart)*Q[4]) - ((M[17]/Vcap[4])*Q[4]) + ((M[18]/Vtis[4])*lambda6[1]*Q[4]) - ((M[17]/Vcap[4])*Q[4]*lambda6[2]) ;  //capillary
        dMdt[18] =  ((M[17]/Vcap[4])*Q[4]) - ((M[18]/Vtis[4])*lambda6[1]*Q[4]) - (lambda6[3]*M[18]) ; //tissue
        dMdt[19]  = (lambda6[3]*M[18]) ; //seq

        //Spleen
        dMdt[20]  = ((M[29]/Vart)*Q[5]) - ((M[20]/Vcap[5])*Q[5]) + ((M[21]/Vtis[5])*lambda7[1]*Q[5]) - ((M[20]/Vcap[5])*Q[5]) ;  //capillary
        dMdt[21] =  ((M[20]/Vcap[5])*Q[5]) - ((M[21]/Vtis[5])*lambda7[1]*Q[5]) - (lambda7[3]*M[21]) ; //tissue
        dMdt[22]  = (lambda7[3]*M[21]) ; //seq

        // Brain
        dMdt[23]  =  ((M[29]/Vart)*Q[6]) - ((M[23]/Vcap[6])*Q[6]) + ((M[24]/Vtis[6])*lambda8[1]*Q[6]) - ((M[23]/Vcap[6])*Q[6]*lambda8[2]) ;  //capillary
        dMdt[24]  =  (k_B*M[1]) + ((M[23]/Vcap[6])*Q[6]) - ((M[24]/Vtis[6])*lambda8[1]*Q[6]) - (lambda8[3]*M[24]) ; //tissue
        dMdt[25]  =  (lambda8[3]*M[24]) ; //seq

        //Others
        dMdt[26]  = ((M[29]/Vart)*Q[7]) - ((M[26]/Vcap[7])*Q[7]) + ((M[27]/Vtis[7])*lambda9[1]*Q[7]) - ((M[26]/Vcap[7])*Q[7]*lambda9[2]) ;  //capillary
        dMdt[27] =  ((M[26]/Vcap[7])*Q[7]) - ((M[27]/Vtis[7])*lambda9[1]*Q[7]) - (lambda9[3]*M[27]) ; //tissue
        dMdt[28]  = (lambda9[3]*M[27]) ; //seq

        //Blood
        dMdt[29] =  (kb*M[5]) - ((M[29]/Vart)*Q[7])- ((M[29]/Vart)*Q[6])- ((M[29]/Vart)*Q[5])- ((M[29]/Vart)*Q[4])- ((M[29]/Vart)*Q[3])-
                    ((M[29]/Vart)*Q[2])- ((M[29]/Vart)*Q[1]);   //art
        dMdt[30] =  ((M[23]/Vcap[6])*Q[6]*lambda8[2])  + ((M[7]/Vcap[1])*Q_liver*lambda3[2]) + ((M[17]/Vcap[4])*Q[4]*lambda6[2]) +
                  ((M[26]/Vcap[7])*Q[7]*lambda9[2]) + ((M[14]/Vcap[3])*Q[3]*lambda5[2]) - ((M[30]/Vven)*lambdaI*QTotal) ; //ven

        //Feces
        dMdt[31] = (lambda4[4]*M[13]) ;

        //Urine
        dMdt[32] =  (M[14]*lambda5[4]) ;

        return dMdt;
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
      real  params[29];      // Matrix containing the individual parameters
      real rel_tol;
      real abs_tol;
      real max_num_steps;
}

////////////////////////////////////////////////////////////////////////////
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
    vector [N_param] theta_tr;
    }

////////////////////////////////////////////////////////////////////

model{
real m_hat[6,32];
real total_m_hat[6,N_compart];


//priors
sigma1 ~ normal(0,10);
sigma2 ~ normal(0,50);

theta_tr[:] ~normal(eta[:],H[:]);


//likelihood

m_hat[:,:] = integrate_ode_bdf(pbpk, m0, t_init, time,
            to_array_1d(theta_tr[:]), params[:], idata,
            rel_tol, abs_tol, max_num_steps);


for (i in 1:6){
//Total amount of NPs in each organ

//Amount in kidneys
total_m_hat[i,1] = m_hat[i,14] + m_hat[i,15] + m_hat[i,16];

//Amount in brain
total_m_hat[i,2] = m_hat[i,23] + m_hat[i,24] + m_hat[i,25];

//Amount in spleen
total_m_hat[i,3] = m_hat[i,20] + m_hat[i,21] + m_hat[i,22];

//Ammount in lungs
total_m_hat[i,4] = m_hat[i,1] + m_hat[i,2] + m_hat[i,3] + m_hat[i,4] + m_hat[i,5] + m_hat[i,6];

//Amount in liver
total_m_hat[i,5] = m_hat[i,7] + m_hat[i,8] + m_hat[i,9];

//Amount in blood
total_m_hat[i,6] = m_hat[i,29] + m_hat[i,30];
}

for (j in 1:N_compart){
      to_vector(m_hat[:,j]) ~ normal(to_vector(total_m_hat[:,j]),sigma1);
  }
}

generated quantities{
      vector [N_param] theta;
      theta = exp(theta_tr);


}
