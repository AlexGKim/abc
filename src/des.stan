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

  real<lower=0> zmin;
  real<lower=zmin> zmax;

  vector[N_adu_max] adu_obs[N_obs];
  vector[N_adu_max] adu_mis[N_sn-N_obs];

  vector[N_obs] trans_ainv_obs;
  int<lower=0, upper=1> snIa_obs[N_obs];

  vector[N_obs] host_zs_obs;

  vector[N_sn-N_obs] host_zs_mis;  // need to do this way in case N_mis =0, transformed data is used
  
  int n_int;

  int N_SNIa;
}

transformed data {

  int N_mis;

  real ainv_int[n_int];

  # real snIa_logit;
  # real nonIa_logit;

  int index;

  real ln10d25;

  real x_r[0];
  real ainv0;

  real zmin3;
  real zmax3;

  vector[N_obs] host_zs_obs3;

  vector[N_sn-N_obs] host_zs_mis3;  // need to do this way in case N_mis =0, transformed data is used

  vector[N_SNIa] adu_SNIa;
  vector[N_obs-N_SNIa] adu_nonIa;
  int index_SNIa[N_SNIa];
  int index_nonIa[N_obs-N_SNIa];


  zmin3 <-(zmin*0.5)^3;
  zmax3 <-(zmax*1.5)^3;

  ainv0 <- 1;

  ln10d25 <- log(10.)/2.5;

  # snIa_logit <- logit(1-1e-6);
  # nonIa_logit <- logit(1e-6);

  N_mis <- N_sn-N_obs;

  for (i in 1:N_obs){
    host_zs_obs3[i]<-pow(host_zs_obs[i],3.);
  }
  for (i in 1:N_mis){
    host_zs_mis3[i]<-pow(host_zs_mis[i],3.);
  }

  for (i in 1:n_int){
    ainv_int[i] <- (1+zmin*.1) + (zmax*1.5-zmin*.1)/(n_int-1)*(i-1);
  }

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
}

parameters{

  // transient parameters
  real <lower=0> alpha_Ia;
  real <lower=0> alpha_nonIa;
  real <lower=0> sigma_Ia;  //these sigmas are in log space
  real <lower=0> sigma_nonIa;


  // cosmology parameters
  real <lower=0.0, upper=1> Omega_M;
  real <lower=0.0, upper=1> Omega_L;
  real <lower=-2, upper=0> w;

  // relative rate parameter
#  real<lower=0, upper=1> snIa_rate;
  simplex[2] snIa_rate;

  // true redshifts spectroscopy
  vector<lower=1+zmin*0.1, upper=1+zmax*1.5>[N_obs] ainv_true_obs;
  vector<lower=1+zmin*0.1, upper=1+zmax*1.5>[N_mis] ainv_true_mis;
 }

transformed parameters{


  // log probability
  // P(ADU|theta_t,g)
//  vector[2] lp_obs[N_obs];
  vector[2] lp_mis[N_mis];

  // log probability
  // P(\theta_g|g)
  vector[2] lp_gal_obs[N_obs];
  vector[2] lp_gal_mis[N_mis];


  vector[N_obs] logainv_true_obs;
  vector[N_mis] logainv_true_mis;

  vector[N_SNIa] adu_true_SNIa;
  vector[N_obs-N_SNIa] adu_true_nonIa;

  vector[2] snIa_rate_logit;

  snIa_rate_logit[1] <- logit(snIa_rate[1]);
  snIa_rate_logit[2] <- logit(snIa_rate[2]);


  logainv_true_obs <- log(ainv_true_obs);
  logainv_true_mis <- log(ainv_true_mis);



  {
    real luminosity_distance_int_s[n_int];
    real y2[n_int];

    real adu_true_mis;

    // get luminosity distance at nodes
    {
      real theta[3];

      int x_i[0];
      real r0[1];
      real luminosity_distance_int[n_int,1];
       // the initial condition for solving the cosmology ODE
      r0[1] <- 0;

      theta[1] <- Omega_M;
      theta[2] <- Omega_L;
      theta[3] <- w;
      luminosity_distance_int <- integrate_ode(friedmann, r0, ainv0, ainv_int, theta, x_r, x_i);
      for (m in 1:n_int){
        luminosity_distance_int_s[m] <- ainv_int[m]*luminosity_distance_int[m,1];
      }
    }

    // get spline at the desired coordinates
    // should work well in 1+z vs d_L since the relationship is linear
    y2 <- spline(ainv_int, luminosity_distance_int_s,n_int, 0., 0.);

    
    for (s in 1:N_SNIa){
      adu_true_SNIa[s] <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, ainv_true_obs[index_SNIa[s]])^(-2);
    }

    for (s in 1:(N_obs-N_SNIa)){
      adu_true_nonIa[s] <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, ainv_true_obs[index_nonIa[s]])^(-2);
    }
 
    for (s in 1:N_obs) {
        //predicted adu

      //adu_true_obs[s]  <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, ainv_true_obs[s])^(-2);
      // marginalize over type
      //lp_obs[s][1] <- normal_log(adu_obs[s][1], adu_true_obs*alpha_Ia, adu_true_obs*alpha_Ia*sigma_Ia*ln10d25) + log(snIa_rate[1]) +  bernoulli_logit_log(snIa_obs[s],snIa_logit);
      //lp_obs[s][2] <- normal_log(adu_obs[s][1], adu_true_obs*alpha_nonIa, adu_true_obs*alpha_nonIa* sigma_nonIa*ln10d25)+ log(snIa_rate[2])  + bernoulli_logit_log(snIa_obs[s], nonIa_logit) ;

      // marginalize over possible hosts
      // *********  P(gals|z) \propto P(z|gals)/ P(z) assume flat P(z) ********
      /*
      for (t in 1:2){
        lp_gal_obs[s][t] <- lognormal_log(1+host_zs_obs[s][t][1], logainv_true_obs[s], 0.001) + bernoulli_logit_log(1, host_zs_obs[s][t][2]);
      }
      */
      lp_gal_obs[s][1] <- lognormal_log(1+host_zs_obs[s], logainv_true_obs[s], 0.001) + log(0.98);
      lp_gal_obs[s][2] <- uniform_log(host_zs_obs3[s], zmin3, zmax3) + log(0.02);
    }
    
    for (s in 1:N_mis){
      adu_true_mis  <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, ainv_true_mis[s])^(-2);
     // marginalize over type
      lp_mis[s][1] <- normal_log(adu_mis[s][1], adu_true_mis*alpha_Ia,  adu_true_mis*alpha_Ia*sigma_Ia*ln10d25) + log(snIa_rate[1]) ;
      lp_mis[s][2] <- normal_log(adu_mis[s][1], adu_true_mis*alpha_nonIa, adu_true_mis*alpha_nonIa*sigma_nonIa*ln10d25)+ log(snIa_rate[2]);

      // marginalize over possible hosts
      // *********  P(gals|z) \propto P(z|gals)/ P(z) assume flat P(z) ********
      /*
      for (t in 1:2){
        lp_gal_mis[s][t] <- lognormal_log(1+host_zs_mis[s][t][1], logainv_true_mis[s], 0.001) + bernoulli_logit_log(1, host_zs_mis[s][t][2]);
      }
      */
      lp_gal_mis[s][1] <- lognormal_log(1+host_zs_mis[s], logainv_true_mis[s], 0.001) + log(0.98);
      lp_gal_mis[s][2] <- uniform_log(host_zs_mis3[s], zmin3,zmax3) + log(0.02);
    }
  
  }
}

model{

  // magnitude zeropoint constrained by a prior of nearby SNe
  alpha_Ia ~ lognormal(log(2.),0.02*ln10d25);

  // P(z_s|g)
  trans_ainv_obs ~ lognormal(logainv_true_obs,0.001);

  //collect classified SNe and vectorize accordingly
  adu_SNIa ~ normal(adu_true_SNIa*alpha_Ia, adu_true_SNIa*alpha_Ia*sigma_Ia*ln10d25);
  increment_log_prob(N_SNIa*log(snIa_rate[1]));

  if (N_obs-N_SNIa > 0){
    adu_nonIa ~ normal(adu_true_nonIa*alpha_nonIa, adu_true_nonIa*alpha_nonIa*sigma_nonIa*ln10d25);
    increment_log_prob((N_obs-N_SNIa)*log(snIa_rate[2]));
  }
  
  # from Stan manual 30.1 no speed benefit from vectorization of increment_log_prob 
  for (s in 1:N_obs){
//    increment_log_prob(log_sum_exp(lp_obs[s]));
    increment_log_prob(log_sum_exp(lp_gal_obs[s]));
  }

  for (s in 1:N_mis){
    increment_log_prob(log_sum_exp(lp_mis[s]));
    increment_log_prob(log_sum_exp(lp_gal_mis[s]));
  }

}