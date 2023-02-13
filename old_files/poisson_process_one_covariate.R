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
lds_from_2000 = landslides[which(landslides$YEAR_1 > "2000"),]

# data cleaning
lds_from_2000=lds_from_2000[-which(lds_from_2000$HEADELEV == 9999),]

lds = lds_from_2000
dim(lds)
# 6615 points

# plot HEADELEV
plot(lds[7])

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

covar <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$HEADELEV)
covar
#### MODELLO ####

#POINT PATTERN
library(maptools)

# model for ONE covariate
#p.sp  <- as(lds[2,3,4,7,8,9], "Spatial")  # Create Spatial* object
p.sp  <- as(lds[7], "Spatial")  # Create Spatial* object
cov.ppp <- as(p.sp, "ppp")      # Create ppp object for the covariate
cov.ppp

genDat_pp <- function(cov.ppp, covar, b0, b1, dim, plotdat = TRUE){
  
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
  
  # Set the coefficients
  beta0 <- b0
  beta1 <- b1
  
  gridcov <- cov.ppp
  gridcov <- as.matrix(density(gridcov))
  
  # Simulate the point pattern
  pp <- rpoispp(im(exp(beta0 + beta1*gridcov), xcol=x0, yrow=y0))
  
  # We count the number of points in each grid cell, plot it and make a vector out of it
  qcounts <- quadratcount(pp, ny=20, nx=20)
  Lambda <- as.vector(t(qcounts)) # We have to use t() as we need to construct the vector with the column first
  #Lambda is long as dim[1]*dim[2]
  
  if(plotdat == TRUE){
    par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1,0.5,0))
    plot(im(covar), main = 'Covariate')
    plot(density(pp), main = 'Intensity')
  }
  # Return a list with which I can play with
  return(list(Lambda = Lambda, pp = pp, covar = covar, gridcov = gridcov)) 
}

set.seed(123)
b0 <- 2
b1 <- 5
dimen <- rbind(c(0,20), c(0,20))
#dim <- rbind(c(801587.1, 859623.4), c(802346, 846054.5)) # dimension of the territory we are observing

pp <- genDat_pp(cov.ppp, covar, b0, b1, dimen)

ppm_stan <- '
data{
  int<lower = 1> n;
  vector[n] x;
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
  target += poisson_log_lpmf(y | beta0 + beta1 * x);
}
generated quantities{
  vector[n] lambda_rep;
  lambda_rep = exp(beta0 + beta1 * x);
}'

library(rstan)

stan_data = list(n = length(pp$Lambda), x = as.vector(t(pp$covar)), y = pp$Lambda, b0 = b0, b1 = b1)
fit_stan <- stan(model_code = ppm_stan, data = stan_data, 
                 warmup = 500, iter = 2000, chains = 3)

print(fit_stan, pars = c('beta0', 'beta1'))

lambda_rep <- as.data.frame(rstan::extract(fit_stan)['lambda_rep'])
mean_lambda_rep <- apply(lambda_rep, 2, 'mean')

pointp <- pp$pp
pp_sp <- as.SpatialPoints.ppp(pointp)
pp_sf <- st_as_sf(pp_sp)
grid <- st_make_grid(pp_sf, n = 20, what = 'polygons') %>% st_as_sf()
grid$pred_stan <- mean_lambda_rep
grid$intensity <- as.vector(t(pp$gridcov))
plot(grid)

rmse <- function(true, observed){sqrt((true - observed)^2)}
sum(rmse(grid$intensity, grid$pred_stan))
# 3243.583

# it seems that considering only HEADELEV as a covariate we don't get an explanatory model
