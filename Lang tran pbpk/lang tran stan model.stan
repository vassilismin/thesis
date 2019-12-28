functions{
        real[] pbpk(real t,
                    real[] E,
                    real[] theta,
                    real[] rdata,
                    int[] idata) {

          real dEdt[30] ;

          real V_gi; real V_br; real V_ht; real V_ki; real V_li; real V_spl; real V_re; real V_blood; real Vven; real Vart;    //volume of tissues

          real Vb_gi; real Vb_br; real Vb_ht; real Vb_ki; real Vb_li; real Vb_spl; real Vb_re;          //volume of capillary in each organ

          real Q_tot; real Q_gi; real Q_br; real Q_ht; real Q_ki; real Q_spl; real Q_li; real Q_re;     //blood flows

          real Da; real Do; real Du; real kB; real kb; real ko; real ku; real kt; real kr; real kd; real ki; real kl; real fo; real fu; real ft;

          real lamda31; real lamda32; real lamda33; real lamda41; real lamda42; real lamda43; real lamda51; real lamda52; real lamda53; real lamda61; real lamda62; real lamda63;
          real lamda71; real lamda72; real lamda73; real lamda81; real lamda82; real lamda83; real lamda91; real lamda92; real lamda93; real lamda35; real lamda44; real lamda54;
          real lamdaI; real lamdav;


          // tissue volumes (in ml)
          V_li = rdata[1];
          V_gi = rdata[2];
          V_ki = rdata[3];
          V_ht = rdata[4];
          V_spl = rdata[5];
          V_br = rdata[6];
          V_re = rdata[7];

          // volume of capillary blood (in ml)
          Vb_li = rdata[8];
          Vb_gi = rdata[9];
          Vb_ki = rdata[10];
          Vb_ht = rdata[10];
          Vb_spl = rdata[11];
          Vb_br = rdata[11];
          Vb_re = rdata[12];
          V_blood = rdata[13];

          //Regional blood flows (in mL per day)
          Q_tot = rdata[14];
          Q_gi = rdata[15];
          Q_ki = rdata[16];
          Q_ht = rdata[17];
          Q_spl = rdata[18];
          Q_br = rdata[19];
          Q_re = rdata[20];

          Da = rdata[21];
          Du = rdata[21];
          Do = rdata[22];

          kB = rdata[23];
          kb = rdata[24];
          ko = rdata[25];
          ku = rdata[26];
          kt = rdata[27];
          kr = rdata[28];
          kd = rdata[29];
          ki = rdata[30];
          kl = rdata[31];

          fo = rdata[32];
          fu = rdata[33];
          ft = rdata[34];

          Vven = rdata[35];
          Vart = rdata[36];


        /////params/////
        lamda31 = exp(theta[1]);
        lamda41 = exp(theta[2]);
        lamda51 = exp(theta[3]);
        lamda61 = exp(theta[4]);
        lamda71 = exp(theta[5]);
        lamda81 = exp(theta[6]);
        lamda91 = exp(theta[7]);

        lamda32 = exp(theta[8]);
        lamda42 = exp(theta[9]);
        lamda52 = exp(theta[10]);
        lamda62 = exp(theta[11]);
        lamda72 = exp(theta[12]);
        lamda82 = exp(theta[13]);
        lamda92 = exp(theta[14]);

        lamda33 = exp(theta[15]);
        lamda43 = exp(theta[16]);
        lamda53 = exp(theta[17]);
        lamda63 = exp(theta[18]);
        lamda73 = exp(theta[19]);
        lamda83 = exp(theta[20]);
        lamda93 = exp(theta[21]);

        lamda35 = exp(theta[22]);
        lamda44 = exp(theta[23]);
        lamda54 = exp(theta[24]);

        lamdaI = exp(theta[25]);
        kb = exp(theta[26]);
        kl = exp(theta[27]);

//Differential equations

          dEdt[1]  = Do - (ko*E[1]) - (kB*E[1]);   // Olfactory
          dEdt[2]  = Du - (ku*E[2])  ;             // Upper airways
          dEdt[3]  = Da - (kr*E[3])  + (kd*E[4]) - (ki*E[3]);  // Alveolar free
          dEdt[4]  = (kr*E[3]) - (kd*E[4]) - (kt*E[4]) ;          // Alveolar Mac
          dEdt[5]  = (ki*E[3]) - (kl*E[5]) - (kb*E[5]) + ((E[30]/Vven)*lamdaI*Q_tot*lamdav);    // I
          dEdt[6]  = kl*E[5];  // Lymph to arterial  from venous Blood

          //Liver
          dEdt[7]  = ((E[29]/Vart)*Q_li) - ((E[7]/Vb_li)*Q_li) + ((E[8]/V_li)*lamda31*Q_li) - ((E[7]/Vb_li)*Q_li*lamda32) - ((E[7]/Vb_li)*Q_li*lamda35) + ((E[10]/Vb_gi)*Q_gi) ;  //capillary
          dEdt[8] =  ((E[7]/Vb_li)*Q_li) - ((E[8]/V_li)*lamda31*Q_li) - (lamda33*E[8])  ; //tissue
          dEdt[9]  = lamda33*E[8]; //seq

          //GI
          dEdt[10]  = (fo*ko*E[1]) + (fu*ku*E[2]) + (ft*kt*E[4]) + (E[7]/Vb_li)*Q_li*lamda35 + ((E[29]/Vart)*Q_gi) - ((E[10]/Vb_gi)*Q_gi) + ((E[11]/V_gi)*lamda41*Q_gi)  - (lamda44*E[10])  ;  //capillary
          dEdt[11]  = ((E[10]/Vb_gi)*Q_gi) - ((E[11]/V_gi)*lamda41*Q_gi) - (lamda43*E[11]) ; //tissue
          dEdt[12]  = lamda43*E[11]; //seq
          dEdt[13]  = 0;

          //Kidneys
          dEdt[14]  = ((E[29]/Vart)*Q_ki) - ((E[14]/Vb_ki)*Q_ki) + ((E[15]/V_ki)*lamda51*Q_ki) - ((E[14]/Vb_ki)*Q_ki*lamda52) - (E[14]*lamda54);  //capillary
          dEdt[15] =  ((E[14]/Vb_ki)*Q_ki) - ((E[15]/V_ki)*lamda51*Q_ki) - (lamda53*E[15])  ; //tissue
          dEdt[16]  = lamda53*E[15]; //seq

          //Heart
          dEdt[17]  = ((E[29]/Vart)*Q_ht) - ((E[17]/Vb_ht)*Q_ht) + ((E[18]/V_ht)*lamda61*Q_ht) - ((E[17]/Vb_ht)*Q_ht*lamda62);  //capillary
          dEdt[18] =  ((E[17]/Vb_ht)*Q_ht) - ((E[18]/V_ht)*lamda61*Q_ht) - (lamda63*E[18]) ; //tissue
          dEdt[19]  = lamda63*E[18]; //seq

          //Spleen
          dEdt[20]  = ((E[29]/Vart)*Q_spl) - ((E[20]/Vb_spl)*Q_spl) + ((E[21]/V_spl)*lamda71*Q_spl) - ((E[20]/Vb_spl)*Q_spl*lamda72);  //capillary
          dEdt[21] =  ((E[20]/Vb_spl)*Q_spl) - ((E[21]/V_spl)*lamda71*Q_spl) - (lamda73*E[21]) ; //tissue
          dEdt[22]  = lamda73*E[21]; //seq

          //Brain
          dEdt[23]  =  ((E[29]/Vart)*Q_br) - ((E[23]/Vb_br)*Q_br) + ((E[24]/V_br)*lamda81*Q_br) - ((E[23]/Vb_br)*Q_br*lamda82);  //capillary
          dEdt[24]  =  kB*E[1] + ((E[23]/Vb_br)*Q_br) - ((E[24]/V_br)*lamda81*Q_br) - (lamda83*E[24]) ; //tissue
          dEdt[25]  =  lamda83*E[24]; //seq

          //Others
          dEdt[26]  = ((E[29]/Vart)*Q_re) - ((E[26]/Vb_re)*Q_re) + ((E[27]/V_re)*lamda91*Q_re) - ((E[26]/Vb_re)*Q_re*lamda92);  //capillary
          dEdt[27] =  ((E[26]/Vb_re)*Q_re) - ((E[27]/V_re)*lamda91*Q_re) - (lamda93*E[27]); //tissue
          dEdt[28]  = lamda93*E[27]; //seq

          //Blood
          dEdt[29] =  kb*E[5] - ((E[29]/Vart)*Q_tot);   //art
          dEdt[30] =  ((E[23]/Vb_br)*Q_br*lamda82)  + ((E[7]/Vb_li)*Q_li*lamda32) + ((E[17]/Vb_ht)*Q_ht*lamda62) + ((E[20]/Vb_spl)*Q_spl*lamda72) + ((E[26]/Vb_re)*Q_re*lamda92) + ((E[14]/Vb_ki)*Q_ki*lamda52) - ((E[30]/Vven)*lamdaI*Q_tot) ;  //ven

          return dEdt;
          }
}

//////////////////////////////////////////////////////////////////////////

data{

        int<lower=0> N_param;                 // Number of parameters to be estimated
        int<lower=0> N_compart;               //number of observed compartments
        int<lower=0> N_diff;              // number of differential equations
        int<lower=0>  N_obs;                 // Total number of observations
        real time[6];
        real mass[6,6];
        int  samp[N_compart];                   // Number of samples of each compartment
        real m0[N_diff];           // Initial concentration in compartments
        real t_init;                  // Initial time
        real eta_tr[N_param];
        real  params[39];      // Matrix containing the individual parameters
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

real m_hat[6,30];
real total_m_hat[6,N_compart];

int pos; // uxiliary variable for segmentation
pos = 1;


//priors
 sigma1 ~ normal(0,10);
 sigma2 ~ normal(0,50);

 theta_tr[:] ~normal(eta[:],H[:]);

 //likelihood~

        m_hat[:,:] = integrate_ode_bdf(pbpk,m0,t_init, time,
                                                  to_array_1d(theta_tr[:]),params[:],idata,
                                                  rel_tol, abs_tol,max_num_steps);
        for (i in 1:6){

        // Total amount of NPs in each organ
        // Amount in kidneys
        total_m_hat[i,1] = m_hat[i,14] + m_hat[i,15] + m_hat[i,16];
        // Amount in brain
        total_m_hat[i,2] = m_hat[i,23] + m_hat[i,24] + m_hat[i,25];
        // Amount in spleen
        total_m_hat[i,3] = m_hat[i,20] + m_hat[i,21] + m_hat[i,22];
        // Amount in Alveolar (lungs)
        total_m_hat[i,4] = m_hat[i,1] + m_hat[i,2] + m_hat[i,3] + m_hat[i,4] + m_hat[i,5] + m_hat[i,6];
        // Amount in liver
        total_m_hat[i,5] = m_hat[i,7] + m_hat[i,8] + m_hat[i,9];
        // Amount in blood
        total_m_hat[i,6] = m_hat[i,29]+m_hat[i,30];
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
