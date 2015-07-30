functions{

  real[] mu (real[] redshifts, real Om, int n_sne, int nzadd){
    real zmin;
    real zmax;
    real zstep;

    real redshifts_sort_fill [2*(n_sne + nzadd) - 1];
    int unsort_inds[n_sne + nzadd];

    real Hinv_sort_fill [2*(n_sne + nzadd) - 1];
    real r_com_sort[n_sne + nzadd];


    real model_mu[n_sne];

    zmin <- 0.;
    zmax <- 2.51;
    zstep <- 0.1;

    for (i in 1: 2*(n_sne + nzadd) - 1) {
        Hinv_sort_fill[i] <- 1./sqrt( Om*pow(1. + redshifts_sort_fill[i], 3) + (1. - Om) );
    }

    r_com_sort[1] <- 0.; // Redshift = 0 should be first element!
    for (i in 2:(n_sne + nzadd)) {
        r_com_sort[i] <- r_com_sort[i - 1] + (Hinv_sort_fill[2*i - 3] + 4.*Hinv_sort_fill[2*i - 2] + Hinv_sort_fill[2*i - 1])*(redshifts_sort_fill[2*i - 1] - redshifts_sort_fill[2*i - 3])/6.;
    }

    for (i in 1:n_sne) {
        model_mu[i] <- 5.*log10((1. + redshifts[i])*r_com_sort[unsort_inds[i] + 1]) + 43.1586133146;
    }
    return model_mu;
  }

  real[] spline(real[] x, real[] y, int n, real yp1, real ypn, int nmax){
    real y2[n];
    real u[nmax];
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

  // differential equation we want to solve is dz/H = dr/sqrt(1-kr2)
  // H^2/H_0^2 = Omega_R a-4 + Omega_M a-3 + Omega_k a-2 + Omega_L
  real[] friedmann(real z,
          real[] r,
          real[] theta,
          real[] x_r,
          int[] i_r){
    real H2;
    real drdz[1];
    real omega_k;
    omega_k <- (1-theta[1]-theta[2]);
    H2 <- theta[1]*(1+z)^3 + theta[2]+ omega_k*(1+z)^2;
    drdz[1] <- sqrt((1-fabs(omega_k)/omega_k*r[1])/H2);
    return drdz;
  }

/*
  real[,] luminosity_distance (real[] redshifts, real[] theta, int n_sne, int bin){

    real omega_k;

    real redshifts_sort_fill [bin*n_sne+1];

    int unsort_inds[n_sne];
    real sorted_redshifts[n_sne];

    real Hinv_sort_fill [bin*n_sne+1];

    real r_com_sort[bin*n_sne+1];

    real ans[n_sne, 1];

    unsort_inds <- sort_indices_asc(redshifts);
    sorted_redshifts <- sort_asc(redshifts);

    redshifts_sort_fill[1] <- 0;
    for (j in 2:bin){
      redshifts_sort_fill[j]<-sorted_redshifts[1]/bin*(j-1);
    }
    for (i in 2:n_sne)
      for (j in 1: bin)
        redshifts_sort_fill[(i-1)*bin+j]<-sorted_redshifts[i-1]+(sorted_redshifts[i]-sorted_redshifts[i-1])/bin*(j-1);
    redshifts_sort_fill[bin*n_sne+1]<-sorted_redshifts[n_sne];

    print(redshifts);
    print(redshifts_sort_fill);
    for (i in 1:n_sne)
      print(redshifts[(unsort_inds[i]-1)*bin+1]);


    omega_k <- (1-theta[1]-theta[2]);
    for (i in 1: bin*n_sne+1) {
        Hinv_sort_fill[i] <- 1./sqrt(  theta[1]*(1+redshifts_sort_fill[i])^3 + theta[2]+ omega_k*(1+redshifts_sort_fill[i])^2 );
    }

    r_com_sort[1] <- 0.; // Redshift = 0 should be first element!
    for (i in 2:bin*n_sne+1) {
        r_com_sort[i] <- r_com_sort[i - 1] + (Hinv_sort_fill[2*i - 3] + 4.*Hinv_sort_fill[2*i - 2] + Hinv_sort_fill[2*i - 1])*(redshifts_sort_fill[2*i - 1] - redshifts_sort_fill[2*i - 3])/6.;
    }

    for (i in 1:n_sne)
      ans[i][1] <-r_com_sort[(unsort_inds[i]-1)*bin+1] * (1+redshifts[i]);
    
    return ans;
  }
*/
}

data{
  int<lower=1> N_sn;
  int<lower=0> N_s_obs;
  int<lower=0> N_s_mis;

  int<lower=1> N_adu_max;

  vector[N_adu_max] adu_s_obs[N_s_obs];
  vector[N_adu_max] adu_s_mis[N_s_mis];

  real<lower=0, upper=2> zs_obs[N_s_obs];
  int<lower=0, upper=1> snIa_obs[N_s_obs];


}

transformed data {

  // data needed for solution for the ode including the zeropoints for z and r
  real x_r[0];
  int x_i[0];

  real r0[1];
  real z0;

  r0[1] <- 0;
  z0 <- 0;

  // variables at which the integral is solved
  real z_int[100];
  z_int[1] = 0;
  for (int i in 2:100)
    z_int[i]=2.5/100*(i-1);
}

parameters{

  // transient parameters
  real <lower=0> alpha_Ia;
  real <lower=0> alpha_nonIa;
  real <lower=0, upper=0.2> sigma;
//  real <lower=0, upper=0.2> sigma_Ia;
//  real <lower=0, upper=0.2> sigma_nonIa;

  // cosmology parameters
  real <lower=0,upper=1> theta[2];

  // galaxy parameters
//  vector<lower=0, upper=2>[N_s_obs] zgal_obs;
//  vector<lower=0, upper=2>[N_s_mis] zgal_mis;

  real<lower=0, upper=1> snIa_rate;

  // missing data spectroscopy
  real<lower=0, upper=1.5> zs_mis[N_s_mis];

 }

transformed parameters{
  real luminosity_distance_obs[N_s_obs,1];
  real luminosity_distance_mis[N_s_mis,1];


// for flux data
  vector[2] lp_s_obs[N_s_obs];
  vector[2] lp_s_mis[N_s_mis];

  // luminosity distance
  int unsort_inds[n_sne];
  real sorted_redshifts[n_sne];
  real luminosity_distance[n_sne,1];

  unsort_inds <- sort_indices_asc(redshifts);
  sorted_redshifts <- sort_asc(redshifts);

  luminosity_distance_obs <- integrate_ode(friedmann, r0, z0, zs_obs, theta, x_r, x_i);

  luminosity_distance_obs <- integrate_ode(friedmann, r0, z0, zs_obs, theta, x_r, x_i);
  for (m in 1:N_s_obs)
    luminosity_distance_obs[m,1] <-luminosity_distance_obs[m,1]*(1+zs_obs[m]);

  luminosity_distance_mis <- luminosity_distance(zs_mis, theta, N_s_mis,10);
  for (m in 1:N_s_mis)
    luminosity_distance_mis[m,1] <-luminosity_distance_mis[m,1]*(1+zs_mis[m]);

// marginalize over type

  for (s in 1:N_s_obs) {
    lp_s_obs[s][1] <- normal_log(adu_s_obs[s], luminosity_distance_obs[s,1]*alpha_Ia , luminosity_distance_obs[s,1]*alpha_Ia*sigma) + bernoulli_log(1, snIa_rate) +  bernoulli_log(snIa_obs[s],1-1e-8);
    lp_s_obs[s][2] <- normal_log(adu_s_obs[s], luminosity_distance_obs[s,1]*alpha_nonIa, luminosity_distance_obs[s,1]*alpha_nonIa*sigma)+bernoulli_log(0, snIa_rate) + bernoulli_log(snIa_obs[s],1e-8) ;
  }

  for (s in 1:N_s_mis){
    //flux
    lp_s_mis[s][1] <- normal_log(adu_s_mis[s], luminosity_distance_mis[s,1]*alpha_Ia , luminosity_distance_mis[s,1]*alpha_Ia*sigma) + bernoulli_log(1, snIa_rate);
    lp_s_mis[s][2] <- normal_log(adu_s_mis[s], luminosity_distance_mis[s,1]*alpha_nonIa, luminosity_distance_mis[s,1]*alpha_nonIa*sigma)+bernoulli_log(0, snIa_rate);
  }
  
}

model{
  
  for (s in 1:N_s_obs){
    increment_log_prob(log_sum_exp(lp_s_obs[s]));
  }
  for (s in 1:N_s_mis)
    increment_log_prob(log_sum_exp(lp_s_mis[s]));

  
}