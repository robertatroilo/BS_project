getwd()
setwd("C:/Users/rober/Desktop/Bayesian Statistics/progetto")
list.files()

library(sf)
library(rgdal)
fname <- "ENTLI.gdb.zip"
fname

# new <- st_read(fname)

st_layers(fname)

library(sp)
#methods(st_as_sf)
#methods(st_as_sfc)

landslides = st_read(fname, "ENTLI_Crown")

# consider landslides from 2000 for now
lds_from_2000_with_NA = landslides[which(landslides$YEAR_1 > "2000"),]
# 6629 obs

# data cleaning
u_NA = union(which(lds_from_2000_with_NA$M_WIDTH == 9999),which(lds_from_2000_with_NA$HEADELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$TAILELEV == 9999))
lds_from_2000=lds_from_2000_with_NA[-u_NA,]
lds = lds_from_2000
dim(lds)
# 6591 points

# plot M_WIDTH
plot(lds[2])
# plot HEADELEV
plot(lds[7])
# plot TAILELEV
plot(lds[8])

library(spatstat)

dim <- rbind(c(801587.1, 859623.4), c(802346, 846054.5))
win <- owin(c(dim[1,1],dim[1,2]), c(dim[2,1],dim[2,2]))

spatstat.options(npixel=c(20,20))#c(dim[1,2]-dim[1,1],dim[2,2]-dim[2,1]))

y0 <- seq(win$yrange[1], win$yrange[2],
          length=spatstat.options()$npixel[2])
x0 <- seq(win$xrange[1], win$xrange[2],
          length=spatstat.options()$npixel[1])

define_covariate <- function(northing, easting, npixel = 20, feature){
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

# MODEL with THREE covariates
#p.sp  <- as(lds[2,3,4,7,8,9], "Spatial")  # Create Spatial* object
p.sp1  <- as(lds[2], "Spatial")  # Create Spatial* object
cov1.ppp <- as(p.sp1, "ppp")      # Create ppp object for M_WIDTH
p.sp2  <- as(lds[7], "Spatial")  # Create Spatial* object
cov2.ppp <- as(p.sp2, "ppp")      # Create ppp object for HEADELEV
p.sp3  <- as(lds[8], "Spatial")  # Create Spatial* object
cov3.ppp <- as(p.sp3, "ppp")      # Create ppp object for TAILELEV

genDat_pp3 <- function(cov1.ppp, covar1, cov2.ppp, covar2, cov3.ppp, covar3, b0, b1, b2, b3, dim, plotdat = TRUE){
  
  # Define the window of interest
  win <- owin(c(dim[1,1],dim[1,2]), c(dim[2,1],dim[2,2]))
  
  # set number of pixels to simulate an environmental covariate 
  # (i.e. I decide the number of "small areas" to consider, in each area i will count the number of landslides)
  # dim[1]=dim[2] so that the small areas are squares
  
  spatstat.options(npixel=c(20,20))#c(dim[1,2]-dim[1,1],dim[2,2]-dim[2,1]))
  
  y0 <- seq(win$yrange[1], win$yrange[2],
            length=spatstat.options()$npixel[2])
  x0 <- seq(win$xrange[1], win$xrange[2],
            length=spatstat.options()$npixel[1])
  
  gridcov1 <- cov1.ppp
  gridcov2 <- cov2.ppp
  gridcov3 <- cov3.ppp
  # Set the coefficients
  beta0 <- b0
  beta1 <- b1
  beta2 <- b2
  beta3 <- b3
  
  gridcov1 <- as.matrix(density(gridcov1))
  gridcov2 <- as.matrix(density(gridcov2)) 
  gridcov3 <- as.matrix(density(gridcov3))
  # practically, gridcov1 = gridcov2 = gridcov3 at this point (I could have considered only one but ok)
  
  # Simulate the point pattern
  pp <- rpoispp(im(exp(beta0 + beta1*gridcov1 + beta2*gridcov2 + beta3*gridcov3), xcol=x0, yrow=y0))
  
  # We count the number of points in each grid cell, plot it and make a vector out of it
  qcounts <- quadratcount(pp, ny=20, nx=20)
  Lambda <- as.vector(t(qcounts)) # We have to use t() as we need to construct the vector with the column first
  #Lambda is long as dim[1]*dim[2]
  
  if(plotdat == TRUE){
    par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(1,0.5,0))
    plot(im(covar1), main = 'M_WIDTH')
    plot(im(covar2), main = 'HEADELEV')
    plot(im(covar3), main = 'TAILELEV')
    plot(density(pp), main = 'Intensity')
  }
  # Return a list with which I can play with
  return(list(Lambda = Lambda, pp = pp, gridcov1 = gridcov1, gridcov2 = gridcov2, gridcov3 = gridcov3, 
              covar1 = covar1, covar2 = covar2, covar3 = covar3)) 
}

set.seed(123)
b0 <- 2
b1 <- 5
b2 <- 5
b3 <- 5
dim <- rbind(c(0,20), c(0,20))
#dim <- rbind(c(801587.1, 859623.4), c(802346, 846054.5)) # dimension of the territory we are observing

pp3 <- genDat_pp3(cov1.ppp, covar1, cov2.ppp, covar2, cov3.ppp, covar3, b0, b1, b2, b3, dim)

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
  lambda_rep = exp(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3);
}'

library(rstan)

stan_data3 = list(n = length(pp3$Lambda), x1 = as.vector(t(pp3$covar1)), x2 = as.vector(t(pp3$covar2)), x3 = as.vector(t(pp3$covar3)), y = pp3$Lambda, b0 = b0, b1 = b1, b2 = b2, b3 = b3)
fit_stan3 <- stan(model_code = ppm_stan3, data = stan_data3, 
                  warmup = 500, iter = 2000, chains = 3, control = list(adapt_delta = .99, max_treedepth = 20, stepsize = .8))

# To reduce the diverging transitions I added:
# adapt_delta = 0.99 (0.8 default)
# max_treedepth = 20 (10 default)
# stepsize = 0.8 (1 default)

print(fit_stan3, pars = c('beta0', 'beta1', 'beta2', 'beta3'))
plot(fit_stan3, pars = c('beta0', 'beta1', 'beta2', 'beta3'))

lambda_rep3 <- as.data.frame(rstan::extract(fit_stan3)['lambda_rep'])
mean_lambda_rep3 <- apply(lambda_rep3, 2, 'mean')

pointp <- pp3$pp
pp_sp <- as.SpatialPoints.ppp(pointp)
pp_sf <- st_as_sf(pp_sp)
grid <- st_make_grid(pp_sf, n = c(20,20), what = 'polygons') %>% st_as_sf()
grid$pred_stan <- mean_lambda_rep3
grid$intensity <- as.vector(t(pp3$gridcov1))

plot(grid)

rmse <- function(true, observed){sqrt((true - observed)^2)}
sum(rmse(grid$intensity, grid$pred_stan))
# 3233.047


### TO DO to try to improve the model:

# - change pixel size
# - change prior on lambda
# - change hyperparameters (beta0, beta1, ...)
# - change the covariates
# - create a model with all the non categorical covariates
# - consider landslides before 2000
# - do covariate selection
# - compute AIC to compare models
# - compare models on prediction
# - move to Cox process and Log Gaussian Cox Process
