#Load libraries
#----
library(sf)
library(rgdal)
library(sp)
library(spatstat)
library(rstan)
library(maptools)
library(OpenImageR)

#Load data
#----
fname <- "../ENTLI.gdb.zip"
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
dim(lds)
# 6588 points

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

#would rename this like changeResolution or something 
#just the mean of the feature over the discretized grid 
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


covar4_std = (covar4-min(covar4))/(max(covar4)-min(covar4))
covar5_std = (covar5-min(covar5))/(max(covar5)-min(covar5))

p.sp  <- as(lds, "Spatial")  # Create Spatial* object
lds.ppp <- as(p.sp, "ppp")      # Create ppp object

p.sp1  <- as(lds[2], "Spatial")  # Create Spatial* object
cov1.ppp <- as(p.sp1, "ppp")      # Create ppp object for M_WIDTH
p.sp2  <- as(lds[7], "Spatial")  # Create Spatial* object
cov2.ppp <- as(p.sp2, "ppp")      # Create ppp object for HEADELEV
p.sp3  <- as(lds[8], "Spatial")  # Create Spatial* object
cov3.ppp <- as(p.sp3, "ppp")      # Create ppp object for TAILELEV
p.sp4  <- as(lds[9], "Spatial")  # Create Spatial* object
cov4.ppp <- as(p.sp4, "ppp")      # Create ppp object for ELE_DIFF
p.sp5  <- as(lds[4], "Spatial")  # Create Spatial* object
cov5.ppp <- as(p.sp5, "ppp")      # Create ppp object for SLOPE

#MODEL with ELE_DIFF + SLOPE
genDat_pp3 <- function(covar4, covar5, b0, b4, b5, dim, plotdat = TRUE){
  
  # Define the window of interest
  win <- owin(c(dim[1,1],dim[1,2]), c(dim[2,1],dim[2,2]))
  
  #spatstat.options(npixel=c(n_quadr,n_quadr))
  
  y0 <- seq(win$yrange[1], win$yrange[2],
            length=spatstat.options()$npixel[2])
  x0 <- seq(win$xrange[1], win$xrange[2],
            length=spatstat.options()$npixel[1])

  # Set the coefficients
  beta0 <- b0
  beta4 <- b4
  beta5 <- b5
  
  #gridcov4 <- as.matrix(density(ppp(x=x0,y=y0,window=win,marks=covar4)))
  #gridcov5 <- as.matrix(density(ppp(x=x0,y=y0,window=win,marks=covar5)))

  
  # Simulate the point pattern
  intensity = im(beta0+beta4*covar4+beta5*covar5, xcol = x0, yrow = y0)
  pp <- rpoispp(exp(intensity))#, xcol=x0, yrow=y0)
  
  # We count the number of points in each grid cell, plot it and make a vector out of it
  qcounts <- quadratcount(pp, ny=n_quadr, nx=n_quadr)
  Lambda <- as.vector(t(qcounts))
  lambda_to_plot = im(flipImage(as.matrix(im(Lambda,xcol=x0,yrow=y0)))) #nice 
  
  #Plot the point pattern simulated and the covariates pattern
  if(plotdat == TRUE){
    par(mfrow=c(2,2), mar=c(0.5,0.5,0.5,0.5), mgp=c(1,0.5,0))
    plot(im(covar4), main = 'ELE_DIFF')
    plot(im(covar5), main = 'SLOPE')
    plot(density(lds.ppp), main = 'Real Pattern')
    plot(density(pp), main = 'Generated point pattern')
    #plot(im(flipImage(as.matrix(im(real_lambda,xcol=x0,yrow=y0)))),main="Real Lambda")
    #plot(lambda_to_plot, main = 'Generated Lambda')
    
  }
  
  return(list(Lambda = Lambda, pp = pp, covar4 = covar4, covar5 = covar5)) 
}

set.seed(123)
b0 <- -1.9
b4 <- 0.08
b5 <- 0.06

dim <- rbind(c(0,n_quadr), c(0,n_quadr))

pp3 <- genDat_pp3(covar4, covar5, b0, b4, b5, dim)

ppm_stan3 <- '
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
  target += poisson_log_lpmf(y | beta0 + beta1 * x1 + beta2 * x2 );
}
generated quantities{
  vector[n] lambda_rep;
  lambda_rep = exp(beta0 + beta1 * x1 + beta2 * x2);
}'

#why not start these at 0 like normal?
stan_data3 = list(n = length(pp3$Lambda), x1 = as.vector(t(pp3$covar4)), 
                  x2 = as.vector(t(pp3$covar5)), y = pp3$Lambda, 
                  b0 = b0, b1 = b4, b2 = b5)

fit_stan3 <- stan(model_code = ppm_stan3, data = stan_data3, 
                  warmup = 1000, iter = 8000, chains = 3)



print(fit_stan3, pars = c('beta0', 'beta1', 'beta2'))
plot(fit_stan3, pars = c('beta0', 'beta1', 'beta2'))

######## NATE SANITY CHECK ###################
#manually reconstruct the grid lambda in a different way and compare to chain mean
#just using our model def... 

beta0_rep3 = rstan::extract(fit_stan3)['beta0']
mean_beta0_rep3 = mean(beta0_rep3$beta0)

beta1_rep3 = rstan::extract(fit_stan3)['beta1']
mean_beta1_rep3 = mean(beta1_rep3$beta1)

beta2_rep3 = rstan::extract(fit_stan3)['beta2']
mean_beta2_rep3 = mean(beta2_rep3$beta2)

x1 = as.vector(t(pp3$covar4))
x2 = as.vector(t(pp3$covar5))

mu = exp(mean_beta0_rep3 + mean_beta1_rep3*x1 + mean_beta2_rep3*x2)
#add this plot to chris'
#################################################

lambda_rep3 <- as.data.frame(rstan::extract(fit_stan3)['lambda_rep'])
mean_lambda_rep3 <- apply(lambda_rep3, 2, 'mean')

par(mfrow=c(2,2),mar=c(1,1,1,1))
plot(im(flipImage(as.matrix(im(real_lambda,xcol=x0,yrow=y0)))),main="Real Lambda")
plot(density(lds.ppp),main="Real Intensity")
plot(im((as.matrix(im(mean_lambda_rep3,xcol=x0,yrow=y0)))),main="Pred Lambda")
plot(im(as.matrix(im(mu,xcol=x0,yrow=y0))), main="Nate Lambda") #nate

rmse <- function(true, observed){sqrt((true - observed)^2)}
sum(rmse(grid$lambda_gendata, grid$pred))
#7320.865



#------------------------------------------------
#Stan Prediction on the real lds point pattern (not the simulated one pp3$pp)
qcounts_real <- quadratcount(lds.ppp, ny=n_quadr, nx=n_quadr) 
real_lambda = as.vector(t(qcounts_real))

stan_data3 = list(n = length(real_lambda), x1 = as.vector(t(covar4)), 
                  x2 = as.vector(t(covar5)), y = real_lambda, 
                  b0 = b0, b1 = b4, b2 = b5)

fit_stan3 <- stan(model_code = ppm_stan3, data = stan_data3, 
                  warmup = 800, iter = 5000, chains = 3)

print(fit_stan3, pars = c('beta0', 'beta1', 'beta2'))
plot(fit_stan3, pars = c('beta0', 'beta1', 'beta2'))

lambda_rep3 <- as.data.frame(rstan::extract(fit_stan3)['lambda_rep'])
mean_lambda_rep3 <- apply(lambda_rep3, 2, 'mean')

library(OpenImageR)
par(mfrow=c(2,2),mar=c(1,1,1,1))
plot(im(flipImage(as.matrix(im(mean_lambda_rep3,xcol=x0,yrow=y0)))),main="Pred Lambda")
plot(im(flipImage(as.matrix(im(real_lambda,xcol=x0,yrow=y0)))),main="Real Lambda")
plot(density(lds.ppp),main="Real Intensity")

rmse <- function(true, observed){sqrt((true - observed)^2)}
sum(rmse(mean_lambda_rep3, real_lambda))
