setwd("C:/Users/rober/Desktop/Bayesian Statistics/progetto")
source('mask_creation.R')


### MODEL WITH 3 COVARIATES (ELEDIFF + SLOPE + COVER) ###
model3cov <- '
data{
  int<lower = 1> n;
  vector[n] x1;
  vector[n] x2;
  vector[n] x3;
  int<lower = 0> y[n];
  real b0;
  real b1;
  real b2;
  real b3;
}
parameters{
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  real<lower = 0, upper = 5> sigma_noise;
  vector[n] noise;
}
transformed parameters{
}
model{
  //priors
  target += normal_lpdf(beta0 | b0,5);
  target += normal_lpdf(beta1 | b1,1);
  target += normal_lpdf(beta2 | b2,1);
  target += normal_lpdf(beta3 | b3,1);
  target += uniform_lpdf(sigma_noise | 0,1);
  
  // Prior for the noise
  target += normal_lpdf(noise | 0, sigma_noise);
  
  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3 + noise );
}
generated quantities{
  vector[n] lambda_rep;
  real log_lik[n];
  
  //for (i in 1:n){
  //   sigma_noise_pred = uniform_rng(0,1);
  //  noise_pred[i] = normal_rng(0, sigma_noise_pred);
  // }
  
  lambda_rep = exp(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3 + noise);
  for (i in 1:n){
    log_lik[i] = poisson_log_lpmf(y | beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3 + noise);
  }
}'

real_lambda = q.extracted

b0 <- -1
b_elediff <- 0.1
b_slope <- 0.05
b_nocover <- -0.01

stan_data3 = list(n = length(real_lambda), x1 = as.vector(t(ele.diff.norm)), 
                  x2 = as.vector(t(slope.norm)), x3 = as.vector(t(nocover.vote)), y = real_lambda, 
                  b0 = b0, b1 = b_elediff, b2 = b_slope, b3 = b_nocover)

fit3 <- stan(model_code = model3cov, data = stan_data3, 
             warmup = 4000, iter = 8000, chains = 4)

# saveRDS(fit3, 'model3covariates')

library(loo)

log_lik_3 <- rstan::extract(fit3)['log_lik']
waic(log_lik_3$log_lik)
# Computed from 16000 by 365 log-likelihood matrix
# 
#            Estimate  SE
# elpd_waic -360436.4 0.0
# p_waic      76109.8 0.0
# waic       720872.8 0.0

# 365 (100.0%) p_waic estimates greater than 0.4. We recommend trying loo instead.

loo_3cov <- loo(fit3)
loo_3cov

# Computed from 16000 by 365 log-likelihood matrix
# 
#           Estimate  SE
# elpd_loo -316287.0 0.0
# p_loo      31960.3 0.0
# looic     632574.0 0.0
# ------
#   Monte Carlo SE of elpd_loo is NA.
# Pareto k diagnostic values:
#                      Count Pct. Min. n_eff
# (-Inf, 0.5]   (good)       0     0.0%  <NA>      
#   (0.5, 0.7]   (ok)         0     0.0%  <NA>      
#   (0.7, 1]   (bad)        0     0.0%  <NA>      
#   (1, Inf)   (very bad) 365   100.0%  1  


### MODEL WITH 2 COVARIATES (ELEDIFF + SLOPE) ###
model2cov <- '
data{
  int<lower = 1> n;
  vector[n] x1;
  vector[n] x2;
  int<lower = 0> y[n];
  real b0;
  real b1;
  real b2;
}
parameters{
  real beta0;
  real beta1;
  real beta2;
  real<lower = 0, upper = 5> sigma_noise;
  vector[n] noise;
}
transformed parameters{
}
model{
  //priors
  target += normal_lpdf(beta0 | b0,5);
  target += normal_lpdf(beta1 | b1,1);
  target += normal_lpdf(beta2 | b2,1);
  target += uniform_lpdf(sigma_noise | 0,1);
  
  // Prior for the noise
  target += normal_lpdf(noise | 0, sigma_noise);
  
  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x1 + beta2 * x2 + noise );
}
generated quantities{
  vector[n] lambda_rep;
  real log_lik[n];
  
  //for (i in 1:n){
  //  sigma_noise_pred = uniform_rng(0,1);
  //  noise_pred[i] = normal_rng(0, sigma_noise_pred);
  //}
  
  lambda_rep = exp(beta0 + beta1 * x1 + beta2 * x2 + noise);
  for (i in 1:n){
    log_lik[i] = poisson_log_lpmf(y | beta0 + beta1 * x1 + beta2 * x2 + noise);
  }
}'

stan_data2 = list(n = length(real_lambda), x1 = as.vector(t(ele.diff.norm)), 
                  x2 = as.vector(t(slope.norm)), y = real_lambda, 
                  b0 = b0, b1 = b_elediff, b2 = b_slope)

fit2 <- stan(model_code = model2cov, data = stan_data2, 
             warmup = 4000, iter = 8000, chains = 4)

log_lik_2 <- rstan::extract(fit2)['log_lik']
waic(log_lik_2$log_lik)
# Computed from 16000 by 365 log-likelihood matrix
# 
# Estimate  SE
# elpd_waic -356444.0 0.0
# p_waic      75205.4 0.0
# waic       712888.0 0.0
# 
# 365 (100.0%) p_waic estimates greater than 0.4. We recommend trying loo instead. 

loo_2cov <- loo(fit2)
loo_2cov
# Computed from 16000 by 365 log-likelihood matrix
# 
#          Estimate  SE
# elpd_loo -316755.5 0.0
# p_loo      35516.9 0.0
# looic     633510.9 0.0
# 
# ------
#   Monte Carlo SE of elpd_loo is NA.
# 
# Pareto k diagnostic values:
#                     Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)       0     0.0%  <NA>      
#   (0.5, 0.7]   (ok)         0     0.0%  <NA>      
#   (0.7, 1]   (bad)        0     0.0%  <NA>      
#   (1, Inf)   (very bad) 365   100.0%  2  

library(loo)
loo_compare(loo_3cov, loo_2cov)

#         elpd_diff se_diff
# model1    0.0       0.0     # The best model is displayed on the top (so loo_3cov is better!)
# model2 -468.5       0.0 

library(bridgesampling)
bayes_factor(
  bridge_sampler(fit3, silent = TRUE), # x1
  bridge_sampler(fit2, silent = TRUE)  # x2
)
# Estimated Bayes factor in favor of x1 over x2: 9.26128
# Warning messages:
# 1: 1211 of the 8000 log_prob() evaluations on the proposal draws produced -Inf/Inf. 
# 2: 1198 of the 8000 log_prob() evaluations on the proposal draws produced -Inf/Inf.

# It seems better to add the dummy covariate about the vegatation cover!
