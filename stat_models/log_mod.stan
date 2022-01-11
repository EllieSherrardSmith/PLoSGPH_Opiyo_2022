// bernoulli_logistic transformed data function
data {
  
  int<lower=1> N;                  // rows of data
  
  int<lower=0> n_t[N];             // Total number of mosquitoes counted
  int<lower=0> d_t[N];             // Number mosquites killed during the test
  
  vector<lower=0>[N] time;       // time predictor e.g. months
  
  int<lower=1> N_eff;           // a random effect eg wall type / location / mosquito species etc
  int<lower=1, upper = N_eff> eff[N];
  
}

parameters {
  //Consider death. This is the proportion of mosquitoes dying (d_t) of all tested (n_t)
  real alpha1[N_eff];
  real alpha2[N_eff];
  
}

model {
  real sp[N];
  
  alpha1 ~ normal(0,10);
  alpha2 ~ normal(0,10);
  
  for (n in 1:N) {
    sp[n] = alpha1[eff[n]] + alpha2[eff[n]] * time[n];
  }
  
  d_t ~ binomial_logit(n_t, sp);
}

generated quantities{
  real sp_ppc[N_eff, 365];// this is to predict for 365 time points so adjust time accordingly!
    
    for(v in 1:N_eff){
      for(t in 1:365){
        sp_ppc[v, t] = binomial_rng(365, inv_logit(alpha1[v] + alpha2[v] * t)) / 365.0;
      }
    }
}
