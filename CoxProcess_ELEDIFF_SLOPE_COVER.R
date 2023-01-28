getwd()
setwd("C:/Users/rober/Desktop/Bayesian Statistics/progetto")
list.files()
#Load libraries
#----
library(sf)
library(rgdal)
library(sp)
library(spatstat)
library(rstan)
library(maptools)

#Load data
#----
fname <- "ENTLI.gdb.zip"
fname
st_layers(fname)

landslides = st_read(fname, "ENTLI_Crown")

# consider landslides since 2000
lds_from_2000_with_NA = landslides[which(landslides$YEAR_1 > "2000"),]
# 6629 obs

# data cleaning
u_NA = union(which(lds_from_2000_with_NA$M_WIDTH == 9999),which(lds_from_2000_with_NA$HEADELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$TAILELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$ELE_DIFF == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$SLOPE == 9999))
lds_from_2000=lds_from_2000_with_NA[-u_NA,]
lds = lds_from_2000
# remove less important kinds of cover (C and D) and wrong records of A and B
lds=lds[-cbind(which(lds$COVER=="C"), which(lds$COVER=="D"), which(lds$COVER==" A"), which(lds$COVER==" B")),]

dim(lds)
# 6556 points

# plot M_WIDTH
plot(lds[2])
# plot HEADELEV
plot(lds[7])
# plot TAILELEV
plot(lds[8])
# plot Elev_diff
plot(lds[9])
plot(density(lds$ELE_DIFF))
# plot SLOPE
plot(lds[4])
# plot COVER
plot(lds[5])
#----

#number of quadrants
n_quadr = 30

spatstat.options(npixel=c(n_quadr,n_quadr))#c(dim[1,2]-dim[1,1],dim[2,2]-dim[2,1]))
dim <- rbind(c(801587.1, 859623.4), c(802346, 846054.5))
win <- owin(c(dim[1,1],dim[1,2]), c(dim[2,1],dim[2,2]))

y0 <- seq(win$yrange[1], win$yrange[2],
          length=spatstat.options()$npixel[2])
x0 <- seq(win$xrange[1], win$xrange[2],
          length=spatstat.options()$npixel[1])


define_covariate <- function(northing, easting, npixel = n_quadr, feature){
  covar <- matrix(data = 0, nrow = npixel, ncol = npixel)
  count <- matrix(data = 0, nrow = npixel, ncol = npixel)
  n <- dim(lds)[1]
  for (i in c(1:n)){
    for(j in c(1:(npixel-1))){
      for(z in c(1:(npixel-1))) {
        if(northing[i]>y0[j] && northing[i]<y0[j+1]){
          if(easting[i]>x0[z] && easting[i]<x0[z+1]){
            # update the mean
            covar[j,z] = (covar[j,z]*count[j,z] + feature[i])/(count[j,z]+1)
            count[j,z] = count[j,z] + 1
          }
        }
      }
    }
  }
  return(covar)
}

covar1 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$M_WIDTH)
covar2 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$HEADELEV)
covar3 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$TAILELEV)
covar4 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$ELE_DIFF)
covar5 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$SLOPE)

# create dummy variable
# A: Totally bare of vegetation
# B: Partially bare of vegetation
A <- which(lds$COVER=="A")   
B <- which(lds$COVER=="B")
nA <- length(A)
nB <- length(B)
cover <- rep(0,nA+nB)
cover[A] = 1
length(cover)

covar6 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = cover)

p.sp  <- as(lds, "Spatial")  # Create Spatial* object
lds.ppp <- as(p.sp, "ppp")      # Create ppp object



ppm_stan3 <- '
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
  target += normal_lpdf(beta1 | b1,10);
  target += normal_lpdf(beta2 | b2,10);
  target += normal_lpdf(beta3 | b3,10);
  target += uniform_lpdf(sigma_noise | 0,1);
  
  // Prior for the noise
  target += normal_lpdf(noise | 0, sigma_noise);
  
  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3 + noise );
}
generated quantities{
  vector[n] lambda_rep;
  lambda_rep = exp(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3 + noise);
}'


#Stan Prediction on the real lds point pattern (not the simulated one pp3$pp)
qcounts_real <- quadratcount(lds.ppp, ny=n_quadr, nx=n_quadr)
real_lambda = as.vector(t(qcounts_real))

b0 <- -1
b4 <- 0.1
b5 <- 0.05
b6 <- -0.01

stan_data3 = list(n = length(real_lambda), x1 = as.vector(t(covar4)), 
                  x2 = as.vector(t(covar5)), x3 = as.vector(t(covar6)), y = real_lambda, 
                  b0 = b0, b1 = b4, b2 = b5, b3 = b6)

fit_stan3 <- stan(model_code = ppm_stan3, data = stan_data3, 
                  warmup = 4000, iter = 8000, chains = 4)


print(fit_stan3, pars = c('beta0', 'beta1', 'beta2', 'beta3','sigma_noise'))
plot(fit_stan3, pars = c('beta1', 'beta2'))
plot(fit_stan3, pars = c('beta0', 'beta3'))

lambda_rep3 <- as.data.frame(rstan::extract(fit_stan3)['lambda_rep'])
mean_lambda_rep3 <- apply(lambda_rep3, 2, 'mean')

############# plotting the chains ###############3
beta0_rep3 = rstan::extract(fit_stan3)['beta0']
mean_beta0_rep3 = mean(beta0_rep3$beta0)

beta1_rep3 = rstan::extract(fit_stan3)['beta1']
mean_beta1_rep3 = mean(beta1_rep3$beta1)

beta2_rep3 = rstan::extract(fit_stan3)['beta2']
mean_beta2_rep3 = mean(beta2_rep3$beta2)

beta3_rep3 = rstan::extract(fit_stan3)['beta3']
mean_beta3_rep3 = mean(beta3_rep3$beta3)

x11()
par(mfrow=c(2,4))
plot(density(beta0_rep3$beta0), main="", xlab="Beta0")
abline(v=c(quantile(beta0_rep3$beta0, 0.025), mean_beta0_rep3, quantile(beta0_rep3$beta0, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))
plot(density(beta1_rep3$beta1), main="", xlab="Beta1")
abline(v=c(quantile(beta1_rep3$beta1, 0.025), mean_beta1_rep3, quantile(beta1_rep3$beta1, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))
plot(density(beta2_rep3$beta2), main="", xlab="Beta2")
abline(v=c(quantile(beta2_rep3$beta2, 0.025), mean_beta2_rep3, quantile(beta2_rep3$beta2, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))
plot(density(beta3_rep3$beta3), main="", xlab="Beta3")
abline(v=c(quantile(beta3_rep3$beta3, 0.025), mean_beta3_rep3, quantile(beta3_rep3$beta3, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))

plot(seq(1,length(beta0_rep3$beta0)), beta0_rep3$beta0, type='l', xlab="iter", ylab='val')
plot(seq(1,length(beta1_rep3$beta1)), beta1_rep3$beta1, type='l', xlab="iter", ylab='val')
plot(seq(1,length(beta2_rep3$beta2)), beta2_rep3$beta2, type='l', xlab="iter", ylab='val')
plot(seq(1,length(beta3_rep3$beta3)), beta3_rep3$beta3, type='l', xlab="iter", ylab='val')

#plot for the noise 
x11()
par(mfrow=c(2,1))
sigma_rep3 = rstan::extract(fit_stan3)['sigma_noise']
mean_sigma_rep3 = mean(sigma_rep3$sigma_noise)
plot(density(beta2_rep3$beta2), main="", xlab="Beta2")
abline(v=c(quantile(sigma_rep3$sigma_noise, 0.025), mean_sigma_rep3, quantile(sigma_rep3$sigma_noise, 0.975)), 
       col = c('blue', 'red', 'blue'), lty=c(2,1,2))
plot(seq(1,length(sigma_rep3$sigma_noise)), sigma_rep3$sigma_noise, type='l', xlab="iter", ylab='val')

#############################


library(OpenImageR)
par(mfrow=c(2,2),mar=c(1,1,1,1))
plot(im(flipImage(as.matrix(im(mean_lambda_rep3,xcol=x0,yrow=y0)))),main="Pred Lambda")
plot(im(flipImage(as.matrix(im(real_lambda,xcol=x0,yrow=y0)))),main="Real Lambda")
plot(im(as.matrix(density(lds.ppp))),main="Real Intensity")

rmse <- function(true, observed){sqrt((true - observed)^2)}
sum(rmse(real_lambda, mean_lambda_rep3))
#944.0025

