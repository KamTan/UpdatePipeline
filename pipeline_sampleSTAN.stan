// STAN model for data with no missingness

//
// Data Block
//
data {
  int<lower=1> N;
  int<lower=1> N_uncensored;
  int<lower=1> N_censored;
  int<lower=0> NC;
  //The observed portion is "data"
  // We observe all survival times, and all X 
  matrix[N_censored, NC] X_censored;
  matrix[N_uncensored, NC] X_uncensored;
  vector<lower=0>[N_censored] times_censored;
  vector<lower=0>[N_uncensored] times_uncensored;
  //For the survival model
  vector[NC] prior_mean_betas;
  vector[NC] prior_sd_betas;
  real prior_mean_mu;
  real prior_sd_mu;
}

transformed data{

}


// Betas are log hazard ratios and intercept is gamma
// Parameters block
//
parameters {
  vector[NC] betas;
  real intercept;
  real<lower=0> alpha; //this is the Weibull shape parameter

}

transformed parameters {


}

//
// Model block
//
model {

  //main model
  betas ~ normal(prior_mean_betas,prior_sd_betas);
  intercept ~ normal(prior_mean_mu,prior_sd_mu);
  alpha ~ normal(1,2);
  target += weibull_lpdf(times_uncensored | alpha,
       exp(-(intercept+X_uncensored*betas )/alpha));
  target += weibull_lccdf(times_censored | alpha,
       exp(-(intercept+X_censored*betas )/alpha));
}


//
// Generated Quantities block -- not used
//
generated quantities {
}










