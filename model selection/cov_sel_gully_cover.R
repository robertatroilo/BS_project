setwd("C:/Users/rober/Desktop/Bayesian Statistics/progetto")

source('mask_creation.R')

# STAN POISSON MODEL WITH ONE COVARIATE (ELEDIFF)
mod1cov <- '
data{
  int<lower = 1> n;
  vector[n] x1;
  int<lower = 0> y[n];
  real b0;
  real b1;
}
parameters{
  real beta0;
  real beta1;
}
transformed parameters{
}
model{
  //priors
  target += normal_lpdf(beta0 | b0,5);
  target += normal_lpdf(beta1 | b1,10);
  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x1);
}
generated quantities{
  vector[n] lambda_rep;
  real log_lik[n];
  lambda_rep = exp(beta0 + beta1 * x1);
  for (i in 1:n){
    log_lik[i] = poisson_log_lpmf(y | beta0 + beta1 * x1);
  }
}'

real_lambda = q.extracted
b0 <- -1


# covariate = elediff

b_elediff <- 0.1

stan_data_elediff = list(n = length(real_lambda), x1 = ele.diff.norm, y = real_lambda, 
                  b0 = b0, b1 = b_elediff)

fit <- stan(model_code = mod1cov, data = stan_data_elediff, 
             warmup = 4000, iter = 8000, chains = 4)

beta0 = rstan::extract(fit)['beta0']
mean_beta0 = mean(beta0$beta0)

beta1 = rstan::extract(fit)['beta1']
mean_beta1 = mean(beta1$beta1)

x11()
par(mfrow=c(1,2))
plot(density(beta0$beta0), main="", xlab="Beta0")
abline(v=c(quantile(beta0$beta0, 0.025), mean_beta0, quantile(beta0$beta0, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))
plot(density(beta1$beta1), main="", xlab="Beta1")
abline(v=c(quantile(beta1$beta1, 0.025), mean_beta1, quantile(beta1$beta1, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))

#some predictions 
pred = ceiling(exp(ele.diff.norm*mean_beta1 + mean_beta0))
par(mfrow = c(1,2))
plot(log(pred))
plot(log(q.extracted))

# SLOPE

b_slope <- 0.05
stan_data_slope = list(n = length(real_lambda), x1 = slope.norm, y = real_lambda, 
                         b0 = b0, b1 = b_slope)

fit_slope <- stan(model_code = mod1cov, data = stan_data_slope, 
            warmup = 4000, iter = 8000, chains = 4)

beta0 = rstan::extract(fit_slope)['beta0']
mean_beta0 = mean(beta0$beta0)

beta1 = rstan::extract(fit_slope)['beta1']
mean_beta1 = mean(beta1$beta1)

x11()
par(mfrow=c(2,2))
plot(density(beta0$beta0), main="", xlab="Beta0")
abline(v=c(quantile(beta0$beta0, 0.025), mean_beta0, quantile(beta0$beta0, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))
plot(density(beta1$beta1), main="", xlab="Beta1-Slope")
abline(v=c(quantile(beta1$beta1, 0.025), mean_beta1, quantile(beta1$beta1, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))

plot(seq(1,length(beta0$beta0)), beta0$beta0, type='l', xlab="iter", ylab='val')
plot(seq(1,length(beta1$beta1)), beta1$beta1, type='l', xlab="iter", ylab='val')

#some predictions 
pred = ceiling(exp(slope.norm*mean_beta1 + mean_beta0))
par(mfrow = c(1,2))
plot(log(pred))
plot(log(q.extracted))


library(loo)

log_lik_slope <- rstan::extract(fit_slope)['log_lik']
waic(log_lik_slope$log_lik)

loo_slope <- loo(fit_slope)
loo_slope

# STAN MODEL WITH TWO COVARIATES (SLOPE + NOCOVER)
mod2cov <- '
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
}
transformed parameters{
}
model{
  //priors
  target += normal_lpdf(beta0 | b0,5);
  target += normal_lpdf(beta1 | b1,10);
  target += normal_lpdf(beta2 | b2,10);
  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x1 + beta2 * x2);
}
generated quantities{
  vector[n] lambda_rep;
  real log_lik[n];
  lambda_rep = exp(beta0 + beta1 * x1 + beta2 * x2);
  for (i in 1:n){
    log_lik[i] = poisson_log_lpmf(y | beta0 + beta1 * x1 + beta2 * x2);
  }
}'


b_nocover <- -0.01

stan_data_slope_nocover = list(n = length(real_lambda), x1 = slope.norm, x2 = nocover.vote, y = real_lambda, 
                         b0 = b0, b1 = b_slope, b2 = b_nocover)

fit_slope_nocover <- stan(model_code = mod2cov, data = stan_data_slope_nocover, 
            warmup = 4000, iter = 8000, chains = 4)

beta0 = rstan::extract(fit_slope_nocover)['beta0']
mean_beta0 = mean(beta0$beta0)

beta1 = rstan::extract(fit_slope_nocover)['beta1']
mean_beta1 = mean(beta1$beta1)

beta2 = rstan::extract(fit_slope_nocover)['beta2']
mean_beta2 = mean(beta2$beta2)


x11()
par(mfrow=c(2,3))
plot(density(beta0$beta0), main="", xlab="Beta0")
abline(v=c(quantile(beta0$beta0, 0.025), mean_beta0, quantile(beta0$beta0, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))
plot(density(beta1$beta1), main="", xlab="Beta1-Slope")
abline(v=c(quantile(beta1$beta1, 0.025), mean_beta1, quantile(beta1$beta1, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))
plot(density(beta2$beta2), main="", xlab="Beta2-NoCover")
abline(v=c(quantile(beta2$beta2, 0.025), mean_beta2, quantile(beta2$beta2, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))

plot(seq(1,length(beta0$beta0)), beta0$beta0, type='l', xlab="iter", ylab='val')
plot(seq(1,length(beta1$beta1)), beta1$beta1, type='l', xlab="iter", ylab='val')
plot(seq(1,length(beta2$beta2)), beta2$beta2, type='l', xlab="iter", ylab='val')

#some predictions 
pred = ceiling(exp(slope.norm*mean_beta1 + mean_beta0))
par(mfrow = c(1,2))
plot(log(pred))
plot(log(q.extracted))


library(loo)

log_lik_slope_nocover <- rstan::extract(fit_slope_nocover)['log_lik']
waic(log_lik_slope_nocover$log_lik)

loo_slope_nocover <- loo(fit_slope_nocover)
loo_slope_nocover

loo_compare(loo_slope, loo_slope_nocover)

#          elpd_diff se_diff
# model2     0.0       0.0
# model1 -5376.6       0.0

# The model with nocover seems better

# STAN MODEL WITH THREE COVARIATES (SLOPE + NOCOVER + GULLY)
mod3cov <- '
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
}
transformed parameters{
}
model{
  //priors
  target += normal_lpdf(beta0 | b0,5);
  target += normal_lpdf(beta1 | b1,10);
  target += normal_lpdf(beta2 | b2,10);
  target += normal_lpdf(beta3 | b3,10);
  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3);
}
generated quantities{
  vector[n] lambda_rep;
  real log_lik[n];
  lambda_rep = exp(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3);
  for (i in 1:n){
    log_lik[i] = poisson_log_lpmf(y | beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3);
  }
}'

b_gully = 1

stan_data_slope_nocover_gully = list(n = length(real_lambda), x1 = slope.norm, x2 = nocover.vote, x3 = gully.vote, y = real_lambda, 
                               b0 = b0, b1 = b_slope, b2 = b_nocover, b3 = b_gully)

fit_slope_nocover_gully <- stan(model_code = mod3cov, data = stan_data_slope_nocover_gully, 
                          warmup = 4000, iter = 8000, chains = 4)

beta0 = rstan::extract(fit_slope_nocover_gully)['beta0']
mean_beta0 = mean(beta0$beta0)

beta1 = rstan::extract(fit_slope_nocover_gully)['beta1']
mean_beta1 = mean(beta1$beta1)

beta2 = rstan::extract(fit_slope_nocover_gully)['beta2']
mean_beta2 = mean(beta2$beta2)

beta3 = rstan::extract(fit_slope_nocover_gully)['beta3']
mean_beta3 = mean(beta3$beta3)

x11()
par(mfrow=c(2,4))
plot(density(beta0$beta0), main="", xlab="Beta0")
abline(v=c(quantile(beta0$beta0, 0.025), mean_beta0, quantile(beta0$beta0, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))
plot(density(beta1$beta1), main="", xlab="Beta1-Slope")
abline(v=c(quantile(beta1$beta1, 0.025), mean_beta1, quantile(beta1$beta1, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))
plot(density(beta2$beta2), main="", xlab="Beta2-NoCover")
abline(v=c(quantile(beta2$beta2, 0.025), mean_beta2, quantile(beta2$beta2, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))
plot(density(beta3$beta3), main="", xlab="Beta3-Gully")
abline(v=c(quantile(beta3$beta3, 0.025), mean_beta3, quantile(beta3$beta3, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))

plot(seq(1,length(beta0$beta0)), beta0$beta0, type='l', xlab="iter", ylab='val')
plot(seq(1,length(beta1$beta1)), beta1$beta1, type='l', xlab="iter", ylab='val')
plot(seq(1,length(beta2$beta2)), beta2$beta2, type='l', xlab="iter", ylab='val')
plot(seq(1,length(beta3$beta3)), beta3$beta3, type='l', xlab="iter", ylab='val')

#some predictions 
pred = ceiling(exp(slope.norm*mean_beta1 + nocover.vote*mean_beta2 + gully.vote*mean_beta3 + mean_beta0))
par(mfrow = c(1,2))
plot(log(pred))
plot(log(q.extracted))

library(loo)

log_lik_slope_nocover_gully <- rstan::extract(fit_slope_nocover_gully)['log_lik']
waic(log_lik_slope_nocover_gully$log_lik)

loo_slope_nocover_gully <- loo(fit_slope_nocover_gully)
loo_slope_nocover_gully

loo_compare(loo_slope_nocover, loo_slope_nocover_gully)

#         elpd_diff se_diff
# model2     0.0       0.0
# model1 -6141.6       0.0

# Also the gully categorical variable is significant

