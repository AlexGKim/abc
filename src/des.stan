functions{

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
 //   real omega_k;
 //   omega_k <- (1-theta[1]-theta[2]);
 //   H2 <- theta[1]*(1+z)^3 + theta[2]+ omega_k*(1+z)^2;
 //   drdz[1] <- sqrt((1-fabs(omega_k)/omega_k*r[1])/H2);

    //do flat for simplicity
    H2 <- theta[1]*(1+z)^3 + (1- theta[1]);
    drdz[1] <- sqrt(1/H2);
    return drdz;
  }
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

  real<lower=0, upper=2> host_zs_obs[N_s_obs];
  real<lower=0, upper=2> host_zs_mis[N_s_mis];
}

transformed data {

  // data needed for doing the conformal distance integral evaluated on a grid
  real x_r[0];
  int x_i[0];

  real r0[1];
  real z0;

  int n_int;
  real z_int[25];

  real galaxy_classification;

  galaxy_classification <- 0.98;

  // the initial condition for solving the cosmology ODE
  r0[1] <- 0;
  z0 <- 0;

  // redshifts at which the integral is solved
  n_int <- 25;

  for (i in 1:n_int)
    z_int[i]<-1.5/n_int*i;
}

parameters{

  // transient parameters
  real <lower=0.2, upper = 4> alpha_Ia;
  real <lower=0.2, upper = 4> alpha_nonIa;
  real <lower=0, upper=0.2> sigma;

  // cosmology parameters
  real <lower=0.1, upper=0.5> Omega_M;
  real <lower=0.5, upper=0.9> Omega_L;

  // relative rate parameter
  real<lower=0, upper=1> snIa_rate;

  // missing data spectroscopy
  real<lower=0, upper=1.5> zs_true_obs[N_s_obs];
  real<lower=0, upper=1.5> zs_true_mis[N_s_mis];

 }

transformed parameters{

  real theta[2];

  // log probability
  vector[2] lp_s_obs[N_s_obs];
  vector[2] lp_s_mis[N_s_mis];

  // log probability
  vector[2] lp_gal_obs[N_s_obs];
  vector[2] lp_gal_mis[N_s_mis];

  //variables for integral
  real <lower=0> luminosity_distance_int[n_int,1];
  
  //variables for spline
  real y2[n_int+1];
  real <lower=0> luminosity_distance_int_s[n_int+1];
  real <lower=0, upper=2.5> z_int_s[n_int+1];

  //predicted adu
  real <lower=0> adu_obs[N_s_obs];
  real <lower=0> adu_mis[N_s_mis];


  theta[1] <- Omega_M;
  theta[2] <- Omega_L;

  // get luminosity distance at nodes
  luminosity_distance_int <- integrate_ode(friedmann, r0, z0, z_int, theta, x_r, x_i);
  luminosity_distance_int_s[1] <-0;
  z_int_s[1] <- 0;
  for (i in 1:n_int){
    luminosity_distance_int_s[i+1] <- luminosity_distance_int[i,1];
    z_int_s[i+1] <- z_int[i];
  }

  // get spline at the desired coordinates
  y2 <- spline(z_int_s, luminosity_distance_int_s,n_int+1, 0., 0., 100);

  for (m in 1:N_s_obs){
    adu_obs[m] <- ((1+zs_obs[m])* splint(z_int_s, luminosity_distance_int_s, y2, n_int+1, zs_true_obs[m]))^(-2);
  }

  for (m in 1:N_s_mis){
    adu_mis[m] <- ((1+zs_true_mis[m])* splint(z_int_s, luminosity_distance_int_s, y2, n_int+1, zs_true_mis[m]))^(-2);
  }

// marginalize over type
  for (s in 1:N_s_obs) {
    lp_s_obs[s][1] <- normal_log(adu_s_obs[s][1], adu_obs[s]*alpha_Ia, adu_obs[s]*alpha_Ia*sigma) + bernoulli_log(1, snIa_rate) +  bernoulli_log(snIa_obs[s],1-1e-6);
    lp_s_obs[s][2] <- normal_log(adu_s_obs[s][1], adu_obs[s]*alpha_nonIa, adu_obs[s]*alpha_nonIa*sigma)+bernoulli_log(0, snIa_rate) + bernoulli_log(snIa_obs[s],1e-6) ;

    lp_gal_obs[s][1] <- normal_log(zs_true_obs[s],host_zs_obs[s],host_zs_obs[s]*0.005) + bernoulli_log(1, galaxy_classification);
    lp_gal_obs[s][2] <- uniform_log(zs_true_obs[s],0.1,2) + bernoulli_log(0, galaxy_classification);
  }

  for (s in 1:N_s_mis){
    //flux
    lp_s_mis[s][1] <- normal_log(adu_s_mis[s][1], adu_mis[s]*alpha_Ia,  adu_mis[s]*alpha_Ia*sigma) + bernoulli_log(1, snIa_rate);
    lp_s_mis[s][2] <- normal_log(adu_s_mis[s][1], adu_mis[s]*alpha_nonIa, adu_mis[s]*alpha_nonIa*sigma)+bernoulli_log(0, snIa_rate);

    lp_gal_mis[s][1] <- normal_log(zs_true_mis[s],host_zs_mis[s],host_zs_mis[s]*0.005) + bernoulli_log(1, galaxy_classification);
    lp_gal_mis[s][2] <- uniform_log(zs_true_mis[s],0.1,2) + bernoulli_log(0, galaxy_classification);
  }
}

model{
  for (s in 1:N_s_obs){
    increment_log_prob(log_sum_exp(lp_s_obs[s]));
    increment_log_prob(log_sum_exp(lp_gal_obs[s]));
    zs_true_obs[s] ~ normal(zs_obs[s],zs_obs[s]*0.005);
  }

  for (s in 1:N_s_mis){
    increment_log_prob(log_sum_exp(lp_s_mis[s]));
    increment_log_prob(log_sum_exp(lp_gal_mis[s]));
  }
}