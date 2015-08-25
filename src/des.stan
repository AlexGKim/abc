functions{

  real[] spline(real[] x, real[] y, int n, real yp1, real ypn){
    real y2[n];
    real u[n-1];
    real sig;
    real p;
    real qn;
    real un;
    int k;

    if (yp1 > .99e30){
      y2[1]<-0;
      u[1]<-0;
    } else {
      y2[1] <- -0.5;
      u[1] <- (3./(x[2]-x[1])) * ((y[2]-y[1])/(x[2]-x[1])-yp1);
    }

    for (i in 2:n-1){
      sig<-(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p <-sig*y2[i-1]+2;
      y2[i]<-(sig-1)/p;
      u[i] <-(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]<-(6*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }

    if (ypn > .99e30){
      qn <- 0;
      un <-0;
    } else {
      qn <- 0.5;
      un <- (3./(x[n]-x[n-1])) * (ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }

    y2[n] <- (un-qn*u[n-1])/(qn*y2[n-1]+1);
    //for (k in n-1:1:-1){
    for (k_ in 1:n-1){
      k<- n-k_;
      y2[k] <- y2[k]*y2[k+1]*u[k];
    }
    return y2;
  }

  real splint(real[] xa, real[] ya, real[] y2a, int n, real x){
    real y;
    real h;
    real a;
    real b;
    int klo;
    int khi;
    int k;
    klo<-1;
    khi <-n;

    while (khi-klo > 1){
      k<-(khi+klo)/2;
      if (xa[k] > x){
        khi <-k;
      } else{
        klo <- k;
      }
    }
    h <- xa[khi]-xa[klo];
    a<-(xa[khi]-x)/h;
    b<-(x-xa[klo])/h;
    y<-a*ya[klo]+b*ya[khi]+((a^3-a)*y2a[klo]+(b^3-b)*y2a[khi])*(h^2)/6;
    return y; 
  }

  // differential equation we want to solve is d(z+1)/H = dr/sqrt(1-kr2)
  // H^2/H_0^2 = Omega_R a-4 + Omega_M a-3 + Omega_k a-2 + Omega_L
  real[] friedmann(real ainv,
          real[] r,
          real[] theta,
          real[] x_r,
          int[] i_r){
    real drdainv[1];
 //   real omega_k;
 //   omega_k <- (1-theta[1]-theta[2]);
 //   drdainv[1] <- theta[1]*ainv^3 + theta[2]+ omega_k*ainv^2;
 //   drdainv[1] <- sqrt((1-fabs(omega_k)/omega_k*r[1]*r[1])/drdainv[1]);

    //do flat for simplicity
    drdainv[1] <- theta[1]*ainv^3 + (1- theta[1])*ainv^(3*(1+theta[3]));
    drdainv[1] <- sqrt(1/drdainv[1]);
    return drdainv;
  }
}

data{
  int<lower=1> N_sn;
  int<lower=0> N_obs;

  int<lower=1> N_adu_max;

  int<lower=1> N_SNIa;

  real ADU0;

  real<lower=0> zmin;
  real<lower=zmin> zmax;

  vector<lower = ADU0>[N_adu_max] adu_obs[N_obs];
  vector<lower = ADU0>[N_adu_max] adu_mis[N_sn-N_obs];

  vector[N_obs] trans_ainv_obs;
  int<lower=0, upper=1> snIa_obs[N_obs];

  vector[N_obs] host_zs_obs;  # not used for now
  vector[N_sn-N_obs] host_zs_mis; 
  vector[N_sn-N_obs] host2_zs_mis; 

//  int n_int;
}

transformed data {

  int N_mis;

//  real ainv_int[n_int];

  //oft used numbers calculated once
  real ln10d25;
  real loggalaxyProb;
  real lognotgalaxyProb;

  //conditions needed for the integration
  real x_r[0];
  real ainv0;

  # real zmin3;
  # real zmax3;

  // containers for subsets of typed Ia and non-Ia
  vector[N_SNIa] adu_SNIa;
  vector[N_obs-N_SNIa] adu_nonIa;
  int index_SNIa[N_SNIa];
  int index_nonIa[N_obs-N_SNIa];

  //containers for all redshifts used for calculation of distances.  More efficient if sorted.
  real ainv_all[N_obs+ 2*(N_sn-N_obs)];  
  int ainv_all_ind_obs[N_obs] ;  
  int ainv_all_ind_mis[N_sn-N_obs];  
  int ainv2_all_ind_mis[N_sn-N_obs]; 

  real galaxyProb;

//  vector[N_sn-N_obs] ainv_zs_mis;  // need to do this way in case N_mis =0, transformed data is used
//  vector[N_sn-N_obs] ainv2_zs_mis;  // need to do this way in case N_mis =0, transformed data is used

//  ainv_zs_mis <- 1+host_zs_mis;
//  ainv2_zs_mis <- 1+host2_zs_mis;

//  real logvolumedensity;

//  vector[N_sn-N_obs] lp_gal_mis_2;

  N_mis <- N_sn-N_obs;

  galaxyProb <- 0.98; 
  ln10d25 <- log(10.)/2.5;
  loggalaxyProb <- log(galaxyProb);
  lognotgalaxyProb <- log(1-galaxyProb);
  

  # zmin3 <-(zmin*0.5)^3;
  # zmax3 <-(zmax*1.5)^3;

  ainv0 <- 1;

  // redshift nodes for interpolation
  /*
  for (i in 1:n_int){
    ainv_int[i] <- (1+zmin*.1) + (zmax*1.5-zmin*.1)/(n_int-1)*(i-1);
  }
*/

  /*  Was for uniform volume probability for second redshift
  logvolumedensity <- 4./3*pi()*((zmax*1.5)^3 - (zmin*.5)^3);
  logvolumedensity <- -log(logvolumedensity) + log(4*pi()); //extra piece from r jacobian
  for (i in 1:N_mis){
    lp_gal_mis_2[i] <-  logvolumedensity + 2*log(host_zs_mis[i]) + lognotgalaxyProb;  // logvolumedensity has 4pi piece of jacaboian
  }
  */

  // vectors for observed SNIa and nonIa used fo efficient calculation
  {
    int iaindex;
    int noniaindex;
    iaindex <- 1;
    noniaindex <- 1;
    for (i in 1:N_obs){
      if (snIa_obs[i] == 1){
        adu_SNIa[iaindex] <- adu_obs[i][1];
        index_SNIa[iaindex] <- i;
        iaindex <- iaindex+1;
      } else{
        adu_nonIa[noniaindex] <- adu_obs[i][1];
        index_nonIa[noniaindex] <- i;
        noniaindex <- noniaindex+1;
      }
    }
  }

  //sorting of all redshifts and identification with original arrays for efficient integration
  {
    int so[N_obs+2*(N_mis)];
    for (i in 1:N_obs){
      ainv_all[i] <- trans_ainv_obs[i];
    }
    for (i in 1:N_mis){
      ainv_all[N_obs+i] <- 1+host_zs_mis[i];
      ainv_all[N_obs+N_mis+i] <- 1+host2_zs_mis[i];
    } 
    so <- sort_indices_asc(ainv_all);
    ainv_all <- sort_asc(ainv_all);

    for (i in 1:N_obs+2*(N_mis)){
      if (so[i] <= N_obs){
        ainv_all_ind_obs[so[i]]<-i;
      } else if (so[i] <= (N_obs+N_mis)){
        ainv_all_ind_mis[so[i]-N_obs]<-i;
      } else {
        ainv2_all_ind_mis[so[i]-(N_obs+N_mis)]<-i;
      }
    }
  }
}

parameters{

  // transient parameters
  real <lower=0> alpha_Ia;
  real <lower=0> alpha_nonIa;
  real <lower=0> sigma_Ia;  //these sigmas are in log space
  real <lower=0> sigma_nonIa;


  // cosmology parameters
  real <lower=0.0, upper=1> Omega_M;
//  real <lower=0.0, upper=1> Omega_L;
  real <lower=-2, upper=0> w;

  // relative rate parameter
  simplex[2] snIa_rate_0;
  simplex[2] snIa_rate_1;

  // true redshifts
  //in principle observed guys have redshift uncertainty but for efficiency ignore for the moment
//  vector<lower=1+zmin*0.1, upper=1+zmax*1.5>[N_obs] ainv_true_obs;
//  vector<lower=1+zmin*0.5, upper=1+zmax*1.5>[N_mis] ainv_true_mis;
 }

transformed parameters{

  // log probability
  // P(ADU|theta_t,g)
  vector[4] lp_mis[N_mis];

  // log probability
  // P(\theta_g|g)
//  vector[2] lp_gal_mis[N_mis];

//  vector[N_obs] logainv_true_obs;

  vector[N_SNIa] adu_true_SNIa;
  vector[N_obs-N_SNIa] adu_true_nonIa;


//  logainv_true_obs <- log(ainv_true_obs);


  {
//    real luminosity_distance_int_s[n_int];
    real luminosity_distance[N_obs+2*N_mis,1];
//    real y2[n_int];

    real adu_true_mis;
    real adu_true_mis2;
    real rate0;
    real rate1;
    real rate2;
    real rate3;
    # vector[2] logsnIa_rate_0;
    # vector[2] logsnIa_rate_1;

    real renorm;

#    vector[N_mis] ainv_true_mis;
#    vector[N_mis] logainv_true_mis;

    // get luminosity distance at nodes
    {
      real theta[3];

      int x_i[0];
      real r0[1];
      //real luminosity_distance_int[n_int,1];
       // the initial condition for solving the cosmology ODE
      r0[1] <- 0;

      theta[1] <- Omega_M;
      theta[2] <- (1-Omega_M);
      theta[3] <- w;

      /*
      old integration
      luminosity_distance_int <- integrate_ode(friedmann, r0, ainv0, ainv_int, theta, x_r, x_i);
      for (m in 1:n_int){
        luminosity_distance_int_s[m] <- ainv_int[m]*luminosity_distance_int[m,1];
      }
      */

      //new direct integration
      luminosity_distance <- integrate_ode(friedmann, r0, ainv0, ainv_all, theta, x_r, x_i);
      for (m in 1:N_obs+2*N_mis){
        luminosity_distance[m,1] <- ainv_all[m]*luminosity_distance[m,1];
      }
    }

    for (s in 1:N_SNIa){
//      adu_true_SNIa[s] <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, trans_ainv_obs[index_SNIa[s]])^(-2);
      adu_true_SNIa[s] <- luminosity_distance[ainv_all_ind_obs[index_SNIa[s]],1]^(-2);  //powers do not seem to be vectorized in Stan
    }

    for (s in 1:(N_obs-N_SNIa)){
//      adu_true_nonIa[s] <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, trans_ainv_obs[index_nonIa[s]])^(-2);
      adu_true_nonIa[s] <- luminosity_distance[ainv_all_ind_obs[index_nonIa[s]],1]^(-2);
    }
 
/*
  old integration
    // get spline at the desired coordinates
    // should work well in 1+z vs d_L since the relationship is linear
    y2 <- spline(ainv_int, luminosity_distance_int_s,n_int, 0., 0.);

    for (s in 1:N_SNIa){
//      adu_true_SNIa[s] <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, ainv_true_obs[index_SNIa[s]])^(-2);
      adu_true_SNIa[s] <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, trans_ainv_obs[index_SNIa[s]])^(-2);
    }

    for (s in 1:(N_obs-N_SNIa)){
      adu_true_nonIa[s] <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, trans_ainv_obs[index_nonIa[s]])^(-2);
    }
 */
    /*
    logsnIa_rate_0 <- log(snIa_rate_0);
    logsnIa_rate_1 <- log(snIa_rate_1);
    */

    for (s in 1:N_mis){

      //the case for correct galaxy attribution
#      adu_true_mis  <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, ainv_zs_mis[s])^(-2); 
#     adu_true_mis  <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, ainv2_zs_mis[s])^(-2); 
 
      adu_true_mis  <- luminosity_distance[ainv_all_ind_mis[s],1]^(-2);
      adu_true_mis2  <- luminosity_distance[ainv2_all_ind_mis[s],1]^(-2); 

#   snIa_rate_0 rate at z=0, snIa_rate_1 rate at zmax
#   [1] for sn Ia, [2] for nonIa
      rate0 <- snIa_rate_0[1]+ (ainv_all_ind_mis[s]-1)/zmax*(snIa_rate_1[1]-snIa_rate_0[1]); //SN Ia good z
      rate1 <- snIa_rate_0[2]+ (ainv_all_ind_mis[s]-1)/zmax*(snIa_rate_1[2]-snIa_rate_0[2]);  //non-Ia good z
      rate2 <- snIa_rate_0[1]+ (ainv2_all_ind_mis[s]-1)/zmax*(snIa_rate_1[1]-snIa_rate_0[1]);
      rate3 <- snIa_rate_0[2]+ (ainv2_all_ind_mis[s]-1)/zmax*(snIa_rate_1[2]-snIa_rate_0[2]);
      // explicit handling of normalization of truncated distribution
      renorm <-  galaxyProb*
        (rate0*erfc((ADU0-adu_true_mis*alpha_Ia)/sqrt(2)/(adu_true_mis*alpha_Ia*sigma_Ia*ln10d25))
          + rate1*erfc((ADU0-adu_true_mis*alpha_nonIa)/sqrt(2)/(adu_true_mis*alpha_nonIa*sigma_nonIa*ln10d25)));
      renorm <-  renorm + (1-galaxyProb)*
        (rate2*erfc((ADU0-adu_true_mis2*alpha_Ia)/sqrt(2)/(adu_true_mis2*alpha_Ia*sigma_Ia*ln10d25))
          + rate3*erfc((ADU0-adu_true_mis2*alpha_nonIa)/sqrt(2)/(adu_true_mis2*alpha_nonIa*sigma_nonIa*ln10d25)));
      renorm <- renorm/2;
      renorm <- -log(renorm);

     // marginalize over type
      lp_mis[s][1] <- normal_log(adu_mis[s][1], adu_true_mis*alpha_Ia,  adu_true_mis*alpha_Ia*sigma_Ia*ln10d25) + log(rate0) + loggalaxyProb + renorm;
      lp_mis[s][2] <- normal_log(adu_mis[s][1], adu_true_mis*alpha_nonIa, adu_true_mis*alpha_nonIa*sigma_nonIa*ln10d25)+ log(rate1) + loggalaxyProb + renorm;

      lp_mis[s][3] <- normal_log(adu_mis[s][1], adu_true_mis2*alpha_Ia,  adu_true_mis2*alpha_Ia*sigma_Ia*ln10d25) + log(rate2) + lognotgalaxyProb + renorm;
      lp_mis[s][4] <- normal_log(adu_mis[s][1], adu_true_mis2*alpha_nonIa, adu_true_mis2*alpha_nonIa*sigma_nonIa*ln10d25)+ log(rate3) + lognotgalaxyProb + renorm;

      // marginalize over possible hosts
      # lp_gal_mis[s][1] <- lognormal_log(1+host_zs_mis[s], logainv_true_mis[s], 0.001) + loggalaxyProb;
      # lp_gal_mis[s][2] <- lognormal_log(1+host2_zs_mis[s], logainv_true_mis[s], 0.001) + lognotgalaxyProb;

//      lp_gal_mis[s][2] <- lp_gal_mis_2[s] ;
    }
  }
}

model{
  // magnitude zeropoint constrained by a prior of nearby SNe
  alpha_Ia ~ lognormal(log(2.),0.02*ln10d25);

  // P(z_s|g)
  //trans_ainv_obs ~ lognormal(logainv_true_obs,0.001);

  //collect classified SNe and vectorize accordingly
  // P(ADU, Ts | rate) = P(ADU| SNIa) P(SNIa|rate) for guys with Ts=1
  for (s in 1:N_SNIa){
    adu_SNIa[s] ~ normal(adu_true_SNIa[s]*alpha_Ia, adu_true_SNIa[s]*alpha_Ia*sigma_Ia*ln10d25) T[ADU0,]; //trancated only works on univariate
  }
  increment_log_prob(N_SNIa*log(snIa_rate_0[1]));

  if (N_obs-N_SNIa > 0){
    for (s in 1:N_obs-N_SNIa){ 
      adu_nonIa[s] ~ normal(adu_true_nonIa[s]*alpha_nonIa, adu_true_nonIa[s]*alpha_nonIa*sigma_nonIa*ln10d25) T[ADU0,];
    }
    increment_log_prob((N_obs-N_SNIa)*log(snIa_rate_0[2]));
  }
  
  # from Stan manual 30.1 no speed benefit from vectorization of increment_log_prob 

  for (s in 1:N_mis){
    increment_log_prob(log_sum_exp(lp_mis[s]));
#    increment_log_prob(log_sum_exp(lp_gal_mis[s]));
  }

}