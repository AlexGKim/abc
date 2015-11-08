functions{

  int testObsIndeces(int N_obs, int N_SNIa,   int[] snIa_obs, int[] index_SNIa,  int[] index_nonIa){

    for (s in 1:N_SNIa){
        if (snIa_obs[index_SNIa[s]] ==0)
          return 0;
    }
    for (s in 1:N_obs-N_SNIa){
        if (snIa_obs[index_nonIa[s]] ==1)
          return 0;
    }
    return 1;
  }

  int testRedshiftIndeces(int N_sn, int N_obs, vector trans_ainv_obs, vector host_zs_mis, vector host2_zs_mis,  real[] ainv_all, int[] ainv_all_ind_obs,  
  int[] ainv_all_ind_mis,  int[] ainv2_all_ind_mis){

    for (s in 1:N_obs){
        if (trans_ainv_obs[s] != ainv_all[ainv_all_ind_obs[s]])
          return 0;
    }
    for (s in 1:N_sn-N_obs){
        if (host_zs_mis[s]+1 != ainv_all[ainv_all_ind_mis[s]])
          return 0;
        if (host2_zs_mis[s]+1 != ainv_all[ainv2_all_ind_mis[s]])
          return 0;
    }
    return 1;
  }

  vector transrate(int type, vector zs, real snIa_rate_0_,real  snIa_rate_1_, real zmax){
    real slope;
    real snIa_rate_1;
    real snIa_rate_0;
    vector[num_elements(zs)] rate;
    if (type ==1){
      snIa_rate_0 <- snIa_rate_0_;
      snIa_rate_1 <- snIa_rate_1_;
    } else {
      snIa_rate_0 <- 1-snIa_rate_0_;
      snIa_rate_1 <- 1-snIa_rate_1_;
    }
    slope <- (snIa_rate_1-snIa_rate_0)/(zmax*1.1);
    rate <- snIa_rate_0+ zs*slope; 
    return rate;
  }


  vector transrate_quad(int type, vector zs, real snIa_rate_0_,real  snIa_rate_1_, real zmax){
    real a2;
    real snIa_rate_1;
    real snIa_rate_0;
    vector[num_elements(zs)] rate;
    if (type ==1){
      snIa_rate_0 <- snIa_rate_0_;
      snIa_rate_1 <- snIa_rate_1_;
    } else {
      snIa_rate_0 <- 1-snIa_rate_0_;
      snIa_rate_1 <- 1-snIa_rate_1_;
    }

    a2 <- (snIa_rate_1-snIa_rate_0)/((zmax*1.1)^2);
    rate <- snIa_rate_0+ a2*(zs .* zs); 
    return rate;
  }


  vector myPhi(real ADU0, vector adu, real alpha, real sigma, real ln10d25){
    vector[num_elements(adu)] ans;
    ans <- adu * alpha;
    ans <- (ADU0 - ans)/sqrt(2.) ./ (ans*sigma*ln10d25);
    # for (s in 1:num_elements(adu)){
    #   ans[s]<-erfc(ans[s]);
    # }
    # ans <- ans/2;

    // use fast approximation
    ans <-(0.07056 * ans .* ans + 1.5976) .* ans;
    ans <- 1.+exp(-ans);  #inverse logit
    return 1.-1. ./ ans;
  }

  vector myRenorm(real ADU0, vector adu, real alpha, real sigma, real ln10d25){
    vector[num_elements(adu)] ans;

    if (ADU0 <=0.){
      for (s in 1:num_elements(adu)){
        ans[s] <- 1.;
      }
    } else {
      ans <- adu * alpha;
      ans <- log(ans);
      for (s in 1:num_elements(adu)){
        ans[s] <- lognormal_ccdf_log(ADU0, ans[s], sigma*ln10d25);
      }
      ans <- exp(ans);
    }
    return ans;
  }

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

  # vector logdifferentialVolume(real Omega_M, real w, vector adu_true, vector ainv, real adu_min, real ainv_min, real adu_max,
  #       real adu_max, real ainv_min){
  #   vector[num_elements(ainv)] logdifferentialVolume;
 
  #   //for already calculated redshifts
  #   logdifferentialVolume <-  Omega_M * ainv .* ainv .* ainv;
  #   for (s in 1:num_elements(ainv)){
  #     logdifferentialVolume[s] <- logdifferentialVolume[s] + (1- Omega_M)*ainv[s]^(3*(1+w));
  #   }
  #   logdifferentialVolume <- -0.5*log(logdifferentialVolume) -log(adu_true) -2*log(ainv);

  #   /*
  #   Calculation of the normalization factor
  #   */
  #   return logdifferentialVolume;
  # }

  real normalization_term3(real ADU0, vector adu, real[] ainv, real adumin, real ainvmin,
      real adumax, real ainvmax, real alpha, real sigma, real alpha2, real sigma2,
      real snIa_rate_0, real snIa_rate_1, real zmax, real ln10d25)
    {
      vector[num_elements(adu)+2] adus;
      vector[num_elements(adu)+2] ainvs;
      vector[num_elements(adu)+2] ans;
      vector[num_elements(adu)+2] rates;
      real sumans;

      adus[1] <- adumin;
      adus[num_elements(adus)] <- adumax;
      ainvs[1] <- ainvmin;
      ainvs[num_elements(adus)] <- ainvmax;
      for (s in 1:num_elements(adu)){
        adus[s+1] <- adu[s];
        ainvs[s+1] <- ainv[s];
      }
      ainvs <- ainvs-1;


      rates <-transrate(1,ainvs,snIa_rate_0,snIa_rate_1,zmax);

      ans <- rates .* myRenorm(ADU0, adus, alpha, sigma, ln10d25);
      ans <- ans + (1-rates) .* myRenorm(ADU0, adus, alpha2, sigma2, ln10d25);
      ans <- ans .* ainvs .* ainvs;

      //Now do trapezoidal rule
      sumans <- 0.;
      for (s in 1:num_elements(adus)-1){
        sumans <- sumans + (ans[s]+ans[s+1])*(ainvs[s+1]-ainvs[s]);
      }
      sumans <- sumans/2.;
      return sumans;
  }

  real normalization_term3_quad(real ADU0, vector adu, real[] ainv, real adumin, real ainvmin,
      real adumax, real ainvmax, real alpha, real sigma, real alpha2, real sigma2,
      real snIa_rate_0, real snIa_rate_1, real zmax, real ln10d25)
    {
      vector[num_elements(adu)+2] adus;
      vector[num_elements(adu)+2] ainvs;
      vector[num_elements(adu)+2] ans;
      vector[num_elements(adu)+2] rates;
      real sumans;

      adus[1] <- adumin;
      adus[num_elements(adus)] <- adumax;
      ainvs[1] <- ainvmin;
      ainvs[num_elements(adus)] <- ainvmax;
      for (s in 1:num_elements(adu)){
        adus[s+1] <- adu[s];
        ainvs[s+1] <- ainv[s];
      }
      ainvs <- ainvs-1;


      rates <-transrate_quad(1,ainvs,snIa_rate_0,snIa_rate_1,zmax);

      ans <- rates .* myRenorm(ADU0, adus, alpha, sigma, ln10d25);
      ans <- ans + (1-rates) .* myRenorm(ADU0, adus, alpha2, sigma2, ln10d25);
      ans <- ans .* ainvs .* ainvs .* ainvs;

      //Now do trapezoidal rule
      sumans <- 0.;
      for (s in 1:num_elements(adus)-1){
        sumans <- sumans + (ans[s]+ans[s+1])*(ainvs[s+1]-ainvs[s]);
      }
      sumans <- sumans/2.;
      return sumans;
    }
  }

data{
  int<lower=1> N_sn;
  int<lower=0> N_obs;

  int<lower=1> N_adu_max;

  int<lower=0> N_SNIa;

  real ADU0;

  real<lower=0> zmin;
  real<lower=zmin> zmax;

  vector<lower = ADU0>[N_obs] adu_obs;
  vector<lower = ADU0>[N_sn-N_obs] adu_mis;

  vector[N_obs] trans_ainv_obs;
  int<lower=0, upper=1> snIa_obs[N_obs];

  vector[N_obs] host_zs_obs;  # not used for now
  vector[N_sn-N_obs] host_zs_mis; 
  vector[N_sn-N_obs] host2_zs_mis; 

  int bias_anal;
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

  // containers for subsets of typed Ia and non-Ia
  vector[N_SNIa] adu_SNIa;
  vector[N_obs-N_SNIa] adu_nonIa;
  vector[N_SNIa] z_SNIa;
  vector[N_obs-N_SNIa] z_nonIa;
  int index_SNIa[N_SNIa];
  int index_nonIa[N_obs-N_SNIa];

  //containers for all redshifts used for calculation of distances.  More efficient if sorted.
  real ainv_all[N_obs+ 2*(N_sn-N_obs)];  
  int ainv_all_ind_obs[N_obs] ;  
  int ainv_all_ind_mis[N_sn-N_obs];  
  int ainv2_all_ind_mis[N_sn-N_obs]; 

  real galaxyProb;

  int N_nonIa;

  real ainvmin[1];
  real ainvmax[1];



  galaxyProb <- 0.98; 

  # if (bias_anal ==1)
  #   galaxyProb <- 0.96;


//  vector[N_sn-N_obs] ainv_zs_mis;  // need to do this way in case N_mis =0, transformed data is used
//  vector[N_sn-N_obs] ainv2_zs_mis;  // need to do this way in case N_mis =0, transformed data is used

//  ainv_zs_mis <- 1+host_zs_mis;
//  ainv2_zs_mis <- 1+host2_zs_mis;

//  real logvolumedensity;

//  vector[N_sn-N_obs] lp_gal_mis_2;

  ainvmin[1] <- 1+zmin;
  ainvmax[1] <- 1+zmax;

  N_mis <- N_sn-N_obs;
  N_nonIa <- N_obs-N_SNIa;



  ln10d25 <- log(10.)/2.5;
  loggalaxyProb <- log(galaxyProb);
  lognotgalaxyProb <- log(1-galaxyProb);

  ainv0 <- 1;

  // vectors for observed SNIa and nonIa used for efficient calculation
  {
    int iaindex;
    int noniaindex;
    iaindex <- 1;
    noniaindex <- 1;
    for (i in 1:N_obs){
      if (snIa_obs[i] == 1){
        adu_SNIa[iaindex] <- adu_obs[i];
        z_SNIa[iaindex] <- trans_ainv_obs[i]-1;
        index_SNIa[iaindex] <- i;
        iaindex <- iaindex+1;
      } else{
        adu_nonIa[noniaindex] <- adu_obs[i];
        z_nonIa[noniaindex] <- trans_ainv_obs[i]-1;
        index_nonIa[noniaindex] <- i;
        noniaindex <- noniaindex+1;
      }
    }
  }
  print("Pass Ia Index Test ",testObsIndeces(N_obs, N_SNIa, snIa_obs,index_SNIa, index_nonIa));

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
    print("Pass Redshift Index Test ",testRedshiftIndeces(N_sn, N_obs, trans_ainv_obs, host_zs_mis, host2_zs_mis, ainv_all, ainv_all_ind_obs,  
      ainv_all_ind_mis, ainv2_all_ind_mis));
  }
}

parameters{

  // transient parameters
  real <lower=0> alpha_Ia;
  real <lower=0> alpha_nonIa;
  real <lower=0, upper=1> sigma_Ia;  //these sigmas are in log space
  real <lower=0, upper=10> sigma_nonIa;


  // cosmology parameters
  real <lower=0.0, upper=1> Omega_M;
//  real <lower=0.0, upper=1> Omega_L;
  real <lower=-2, upper=0> w;

  // relative rate parameter
  # real snIa_rate_0_logit;
  # real snIa_rate_1_logit;
  real <lower=0, upper=1> snIa_rate_0;
  real <lower=0, upper=1> snIa_rate_1;

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

  vector[N_SNIa] adu_true_SNIa;
  vector[N_obs-N_SNIa] adu_true_nonIa;
  vector[N_obs+2*N_mis] adu;

//  real zPDFrenorm;

  # real snIa_rate_0;
  # real snIa_rate_1;

  # snIa_rate_0 <- inv_logit(snIa_rate_0_logit);
  # snIa_rate_1 <- inv_logit(snIa_rate_1_logit);
  {
//    real luminosity_distance_int_s[n_int];
    
//    real y2[n_int];

#    vector[N_mis] ainv_true_mis;
#    vector[N_mis] logainv_true_mis;

    real luminosity_distance[N_obs+2*N_mis,1];
    real adu_min[1,1];
    real adu_max[1,1];

    vector[N_obs] adu_true_obs;
    # vector[N_obs] rate;
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

      and then a spline
      adu_true_SNIa[s] <- splint(ainv_int, luminosity_distance_int_s, y2, n_int, trans_ainv_obs[index_SNIa[s]])^(-2);
      */

      //new direct integration
      luminosity_distance <- integrate_ode(friedmann, r0, ainv0, ainv_all, theta, x_r, x_i);
      for (m in 1:N_obs+2*N_mis){
        luminosity_distance[m,1] <- ainv_all[m]*luminosity_distance[m,1];
      }

      adu_min <- integrate_ode(friedmann, r0, ainv0, ainvmin, theta, x_r, x_i);
      adu_min[1,1] <- ainvmin[1]*adu_min[1,1];
      adu_max <- integrate_ode(friedmann, r0, ainv0, ainvmax, theta, x_r, x_i);
      adu_max[1,1] <- ainvmax[1]*adu_max[1,1];

    }

    /*
      model adu for objects with spectrum used in model section
    */
    for (s in 1:N_SNIa){
      adu_true_SNIa[s] <- luminosity_distance[ainv_all_ind_obs[index_SNIa[s]],1]^(-2);  //powers do not seem to be vectorized in Stan
      adu_true_obs[index_SNIa[s]] <-adu_true_SNIa[s];
    }

    for (s in 1:(N_nonIa)){
      adu_true_nonIa[s] <- luminosity_distance[ainv_all_ind_obs[index_nonIa[s]],1]^(-2);
      adu_true_obs[index_nonIa[s]] <-adu_true_nonIa[s];
    }

    // If there is non Malmquiest bias, this term is independent of parameters and not worth calculating
    for (s in 1:num_elements(ainv_all)){
      adu[s] <- luminosity_distance[s,1]^(-2);
    }

    /*
    if (ADU0 !=0.){
      adu_min[1,1] <- adu_min[1,1]^(-2);
      adu_max[1,1] <- adu_max[1,1]^(-2);
      if (bias_anal ==1){
        zPDFrenorm<- normalization_term3_quad(ADU0, adu, ainv_all, adu_min[1,1], ainvmin[1],
            adu_max[1,1], ainvmax[1], alpha_Ia, sigma_Ia, alpha_nonIa, sigma_nonIa, snIa_rate_0,snIa_rate_1, zmax, ln10d25);
      } else{
        zPDFrenorm<- normalization_term3(ADU0, adu, ainv_all, adu_min[1,1], ainvmin[1],
            adu_max[1,1], ainvmax[1], alpha_Ia, sigma_Ia, alpha_nonIa, sigma_nonIa, snIa_rate_0,snIa_rate_1, zmax, ln10d25);
      }
      zPDFrenorm <- log(zPDFrenorm);
    }
  */
  
        /*

    For the mis set the likelihood is

   p(ADU, g1, g2 | ....) = P_host ( P(ADU| z=g1, zH=g2, T=1, ...) P(T=1| z=g1,...) P(z=g1|...) P(zH=g2|...) 
                          + P(ADU| z=g1, zH=g2, T=0, ...) P(T=0| z=g1,...) P(z=g1|...) P(zH=g2|...) )
                          + (1-P_hosst) (P(ADU| z=g2, zH=g1, T=1, ...) P(T=1| z=g2,...) P(z=g2|...) P(zH=g1|...) 
                           + P(ADU| z=g2, zH=g1, T=0, ...) P(T=0| z=g2,...) P(z=g2|...) P(zH=g1|...) )
  */


  /*
   Do the likelihood for each permuatation of type and redshift
   */
   {
      vector[N_mis] lp_holder[4];
      vector[N_mis] rate;
      vector[N_mis] erfc_Ia;
      vector[N_mis] erfc_nonIa;
      vector[N_mis] renorm;

      vector[N_mis] adu_;

      if (N_mis !=0){
        // do correct redshift first
        if (bias_anal ==1){
          rate <- transrate_quad(1, host_zs_mis, snIa_rate_0, snIa_rate_1, zmax);
        } else {
          rate <- transrate(1, host_zs_mis, snIa_rate_0, snIa_rate_1, zmax);
        }
        for (s in 1:N_mis){
          adu_[s] <- adu[ainv_all_ind_mis[s]];
        }
        for (s in 1:N_mis){
          lp_holder[1][s] <- lognormal_log(adu_mis[s], log(adu_[s]*alpha_Ia), sigma_Ia*ln10d25);
          lp_holder[2][s] <- lognormal_log(adu_mis[s], log(adu_[s]*alpha_nonIa), sigma_nonIa*ln10d25);
        }
        lp_holder[1] <- lp_holder[1] + loggalaxyProb +log(rate);
        lp_holder[2] <- lp_holder[2] + loggalaxyProb +log(1-rate);

        if (ADU0 !=0.){
          erfc_Ia <- myRenorm(ADU0, adu_, alpha_Ia, sigma_Ia, ln10d25);
          erfc_nonIa <- myRenorm(ADU0, adu_, alpha_nonIa, sigma_nonIa, ln10d25);
          renorm <- galaxyProb*( rate .* erfc_Ia + (1-rate) .* erfc_nonIa);

          # erfc_Ia <- log(erfc_Ia);
          # lp_holder[1] <- lp_holder[1] - erfc_Ia;
          # erfc_nonIa <- log(erfc_nonIa);
          # lp_holder[2] <- lp_holder[2] - erfc_nonIa;

#          rate <- rate ./ renorm;
        }

        # if (ADU0 !=0){
        #   for (ind in 1:2)
        #     lp_holder[ind] <- lp_holder[ind]+log(renorm) - zPDFrenorm; // z^2 is independent of parameters
        # }


        #    increment_log_prob(lognormal_log(adu_SNIa, log(adu_true_SNIa*alpha_Ia), sigma_Ia*ln10d25));
        # if (ADU0 !=0){
        #   erfc_Ia <- log(erfc_Ia);
        #   lp_holder[1] <- lp_holder[1] - erfc_Ia;
        #             erfc_nonIa <- log(erfc_nonIa);
        #   lp_holder[2] <- lp_holder[2] - erfc_nonIa;
        # }

        # increment_log_prob(lognormal_log(adu_nonIa, log(adu_true_nonIa*alpha_nonIa), sigma_nonIa*ln10d25));
        # if (ADU0 !=0){
        #   erfc_nonIa <- log(erfc_nonIa);
        #   lp_holder[2] <- lp_holder[2] - erfc_nonIa;
        # }



        // do incorrect redshift next
        if (bias_anal==1){
          rate <- transrate_quad(1, host2_zs_mis, snIa_rate_0, snIa_rate_1, zmax);
        } else{
          rate <- transrate(1, host2_zs_mis, snIa_rate_0, snIa_rate_1, zmax);
        }
        for (s in 1:N_mis){
          adu_[s] <- adu[ainv2_all_ind_mis[s]];
        }
        for (s in 1:N_mis){
          lp_holder[3][s] <- lognormal_log(adu_mis[s], log(adu_[s]*alpha_Ia), sigma_Ia*ln10d25);
          lp_holder[4][s] <- lognormal_log(adu_mis[s], log(adu_[s]*alpha_nonIa), sigma_nonIa*ln10d25);
        }
        lp_holder[3] <- lp_holder[3] + lognotgalaxyProb +log(rate);
        lp_holder[4] <- lp_holder[4] + lognotgalaxyProb +log(1-rate);

        if (ADU0 !=0.){
          erfc_Ia <- myRenorm(ADU0, adu_, alpha_Ia, sigma_Ia, ln10d25);
          erfc_nonIa <- myRenorm(ADU0, adu_, alpha_nonIa, sigma_nonIa, ln10d25);
          renorm <- renorm + (1-galaxyProb)*(rate .* erfc_Ia+ (1-rate) .* erfc_nonIa) ;

          # erfc_Ia <- log(erfc_Ia);
          # lp_holder[3] <- lp_holder[3] - erfc_Ia;
          # erfc_nonIa <- log(erfc_nonIa);
          # lp_holder[4] <- lp_holder[4] - erfc_nonIa;

          renorm <- log(renorm);
          for (s in 1:4){
            lp_holder[s] <- lp_holder[s] - renorm;
          }  
        }        
          # rate <- rate ./ renorm;
        
        # lp_holder[3] <- log(rate);
        # lp_holder[4] <- log(1-rate);
        # if (ADU0 !=0){
        #   for (ind in 3:4)
        #     lp_holder[ind] <- lp_holder[ind]+log(renorm) - zPDFrenorm; // z^2 is independent of parameters
        # }



        #    increment_log_prob(lognormal_log(adu_SNIa, log(adu_true_SNIa*alpha_Ia), sigma_Ia*ln10d25));
        # if (ADU0 !=0){
        #   erfc_Ia <- log(erfc_Ia);
        #   lp_holder[3] <- lp_holder[3] - erfc_Ia;
        # }

        # increment_log_prob(lognormal_log(adu_nonIa, log(adu_true_nonIa*alpha_nonIa), sigma_nonIa*ln10d25));
        # if (ADU0 !=0){
        #   erfc_nonIa <- log(erfc_nonIa);
        #   lp_holder[4] <- lp_holder[4] - erfc_nonIa;
        # }


        for (s in 1:N_mis){
          for (t in 1:4){
            lp_mis[s][t] <-  lp_holder[t][s];
          }
        }
      }
    }
  }
}

model{
  // magnitude zeropoint and intrinsic dispersion constrained by a prior of nearby SNe
  alpha_Ia ~ lognormal(log(2.),0.02*ln10d25);
  sigma_Ia ~ lognormal(log(.1),0.1);

  /*

    For the observed set the likelihood is

   p(ADU, z_S, T_z| ....) = P(ADU| zs, Tz, ...) P(Ts| zs,...) P(zs|...) 
  
  */


  {
    vector[N_obs] rate;
    vector[N_obs] erfc_Ia;
    vector[N_obs] erfc_nonIa;
    vector[N_obs] adu_true;
    vector[N_obs] renorm;

    //second term P(Ts| zs,...)
    if (bias_anal==1){
      rate <- transrate_quad(1,trans_ainv_obs-1, snIa_rate_0, snIa_rate_1, zmax);
    } else{
      rate <- transrate(1,trans_ainv_obs-1, snIa_rate_0, snIa_rate_1, zmax);
    }
    if (ADU0 !=0.){
      for (s in 1:N_SNIa){
        adu_true[index_SNIa[s]] <- adu_true_SNIa[s];
      }
      for (s in 1:N_nonIa){
        adu_true[index_nonIa[s]] <-adu_true_nonIa[s];
      }
      erfc_Ia <- myRenorm(ADU0, adu_true, alpha_Ia, sigma_Ia, ln10d25);
      erfc_nonIa <- myRenorm(ADU0, adu_true, alpha_nonIa, sigma_nonIa, ln10d25);
      renorm <- rate .* erfc_Ia+ (1-rate) .* erfc_nonIa;
      rate <- rate .* erfc_Ia ./ renorm;    
    }

    snIa_obs ~ bernoulli(rate);

    //third term P(zs|...) 
    // if (ADU0 !=0.){
    //   increment_log_prob(log(renorm)); // z^2 is independent of parameters
    //   increment_log_prob(-N_obs*zPDFrenorm);
    // }

    // first term  P(ADU| zs, Tz, ...)

    increment_log_prob(lognormal_log(adu_SNIa, log(adu_true_SNIa*alpha_Ia), sigma_Ia*ln10d25));

    # print(lognormal_log(adu_SNIa[1], log(adu_true_SNIa[1]*alpha_Ia), sigma_Ia*ln10d25));

    if (ADU0 !=0.){
      erfc_Ia <- log(erfc_Ia);
      for (s in 1:N_SNIa)
        increment_log_prob(-erfc_Ia[index_SNIa[s]]);
    }

    if (N_nonIa > 0){
      increment_log_prob(lognormal_log(adu_nonIa, log(adu_true_nonIa*alpha_nonIa), sigma_nonIa*ln10d25));
      if (ADU0 !=0.){
        erfc_nonIa <- log(erfc_nonIa);
        for (s in 1:N_nonIa)
          increment_log_prob(-erfc_nonIa[index_nonIa[s]]);
      }
    }
    
  }

  {
    for (s in 1:N_mis){
      increment_log_prob(log_sum_exp(lp_mis[s]));
    }
  }
}
