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

parameters{

  // transient parameters
  real <lower=-2, upper=2> alpha_Ia;
  real <lower=-2, upper=2> alpha_nonIa;
  real <lower=0, upper=0.2> sigma_Ia;
  real <lower=0, upper=0.2> sigma_nonIa;

  // cosmology parameters
  real <lower=-2, upper=2> beta;

  // galaxy parameters
//  vector<lower=0, upper=2>[N_s_obs] zgal_obs;
//  vector<lower=0, upper=2>[N_s_mis] zgal_mis;

  real<lower=0, upper=1> snIa_rate;

  // missing data spectroscopy
  vector<lower=0, upper=1.5>[N_s_mis] zs_mis;
}

transformed parameters{

// marginalize over type

// for flux data
  vector[2] lp_adu_s_obs[N_s_obs];
  vector[2] lp_adu_s_mis[N_s_mis];

// for spectral data
  vector[2] lp_type_s_obs[N_s_obs];

  for (s in 1:N_s_obs) {
    lp_adu_s_obs[s][1] <- normal_log(adu_s_obs[s], alpha_Ia + beta * zs_obs[s], sigma_Ia) + log(snIa_rate);
    lp_adu_s_obs[s][2] <- normal_log(adu_s_obs[s], alpha_nonIa + beta * zs_obs[s], sigma_nonIa)+log(1-snIa_rate);

    lp_type_s_obs[s][1] <- bernoulli_log(snIa_obs[s],1)+log(snIa_rate);
    lp_type_s_obs[s][2] <- bernoulli_log(snIa_obs[s],0)+log(1-snIa_rate);
  }

  for (s in 1:N_s_mis){
    lp_adu_s_mis[s][1] <- normal_log(adu_s_mis[s], alpha_Ia + beta * zs_mis[s], sigma_Ia) + log(snIa_rate);
    lp_adu_s_mis[s][2] <- normal_log(adu_s_mis[s], alpha_nonIa + beta * zs_mis[s], sigma_nonIa)+log(1-snIa_rate);
  }

}

model{

  for (m in 1:N_s_obs){
    increment_log_prob(log_sum_exp(lp_adu_s_obs[m]));
    increment_log_prob(log_sum_exp(lp_type_s_obs[m]));
  }
  for (m in 1:N_s_mis)
    increment_log_prob(log_sum_exp(lp_adu_s_mis[m]));

}