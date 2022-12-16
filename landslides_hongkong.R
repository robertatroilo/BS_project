getwd()
setwd("C:/Users/rober/Downloads")
list.files()

library(sf)
library(rgdal)
fname <- "ENTLI.gdb.zip"
fname

# new <- st_read(fname)

st_layers(fname)

library(sp)
methods(st_as_sf)
methods(st_as_sfc)

landslides = st_read(fname, "ENTLI_Crown")
# Crown of a landslide: The practically undisplaced material still in place and 
#                       adjacent to the highest parts of the main scarp.
#plot(landslides, max.plot = 14) # to plot all the 14 features (fields)
#st_write(st_as_sf(landslides), fname, "ENTLI_Crown", 
#         layer_options = "OVERWRITE=true", driver = "OpenFileGDB")
# Error, so the driver-specific options are documented in the driver manual of gdal

trail = st_read(fname, "ENTLI_Trail")

library(ggplot2)
ggplot(landslides) + geom_sf()


lds_from_1980 <- landslides[which(landslides$YEAR_1 > "1980"),]
ggplot(lds_from_1980) + geom_sf()
dim(lds_from_1980)


#--------------------------------------------------------------------------------#
#samuele
lds_from_1980=lds_from_1980[-which(lds_from_1980$ELE_DIFF == 9999),]
lds_from_1980=lds_from_1980[-which(lds_from_1980$SLOPE== 9999),]

# first thing we need to deal with the 9999
lds_from_1980=lds_from_1980[-which(lds_from_1980$ELE_DIFF == 9999),]
lds_from_1980=lds_from_1980[-which(lds_from_1980$SLOPE== 9999),]


#what i want to do is to implement a model using as a covariate, at first the sum of the headelev of each landslide
# in a certain square , the try to use the mean elev


library(spatstat)
library(sf)
library(sp)
library(maptools)
library(raster)
library(rstan)
library(tidyverse)
library(cowplot)



lds = lds_from_1980
#POINT PATTERN
library(maptools)#
library(spatstat)
p.sp  <- as(lds, "Spatial")  # Create Spatial* object
lds.ppp <- as(p.sp, "ppp")      # Create ppp object


#p.sf.utm <- st_transform(lds, 32619) # project from geographic to UTM
#p.sp  <- as(p.sf.utm, "Spatial")      # Create Spatial* object
#p.ppp <- as(p.sp, "ppp")              # Create ppp object
class(lds.ppp)


#prima uso come covariata quella dell'esempio di chris

GenDat_pp <- function(df.ppp, dim, plotdat = TRUE){
  
  # Define the window of interest
  win <- owin(c(0,dim[1]), c(0,dim[2]))
  
  # set number of pixels to simulate an environmental covariate
  spatstat.options(npixel=c(dim[1],dim[2]))
  
  y0 <- seq(win$yrange[1], win$yrange[2],
            length=spatstat.options()$npixel[2])
  x0 <- seq(win$xrange[1], win$xrange[2],
            length=spatstat.options()$npixel[1])
  multiplier <- 1/dim[2]
  
  # Make the environmental covariate
  gridcov <- outer(x0,y0, function (x,y) multiplier*y + 0*x)
  
  # We count the number of points in each grid cell, plot it and make a vector out of it
  qcounts <- quadratcount(df.ppp, ny=dim[1], nx=dim[2])
  dens <- density(df.ppp)
  Lambda <- as.vector(t(qcounts)) # We have to use t() as we need to construct the vector with the column first
  
  if(plotdat == TRUE){
    par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1,0.5,0))
    plot(im(gridcov), main = 'Covariate')
    plot(dens, main = 'Intensity')
  }
  # Return a list with which I can play with
  return(list(Lambda = Lambda, pp = df.ppp, gridcov = gridcov))
}

pp1 = GenDat_pp(lds.ppp,c(20,20),TRUE)

#stan model
ppm_stan <- '
data{
  int<lower = 1> n;
  vector[n] x;
  int<lower = 0> y[n];
}
parameters{
  real beta0;
  real beta1;
}
transformed parameters{
}
model{
  //priors
  target += normal_lpdf(beta0 | 0,5);
  target += normal_lpdf(beta1 | 0,10);

  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x);
}
generated quantities{
  vector[n] lambda_rep;
  lambda_rep = exp(beta0 + beta1 * x);
}'


stan_data = list(n = length(pp1$Lambda), x = as.vector(t(pp1$gridcov)) , y = pp1$Lambda)
fit_stan1 <- stan(model_code = ppm_stan, data = stan_data, 
                 warmup = 1000, iter = 5000, chains = 3)

print(fit_stan1, pars = c('beta0', 'beta1'))

lambda_rep1 <- as.data.frame(rstan::extract(fit_stan1)['lambda_rep'])
mean_lambda_rep1 <- apply(lambda_rep1, 2, 'mean')


grid <- st_make_grid(pp_sf, n = 20, what = 'polygons') %>% st_as_sf()
grid$pred_stan <- mean_lambda_rep1

plot(grid)

#capisce che il posto con intensità maggiore è a sinistra ma per il resto molto male


q=quadrats(lds.ppp,nx=20,ny=20)# here i am trying to figure out which are the extremes of each of the 400 squares 
                               # so that i can compute the total height in the square

som=matrix(0,20,20)


#this code calculates the sum of the heights of each square
for (i in c(1:15580)){
  for(j in c(1:20)){
    for(z in c(1:20)) {
      if(lds_from_1980$NORTHING[i]>q$ygrid[j] && lds_from_1980$NORTHING[i]<q$ygrid[j+1]){
        if(lds_from_1980$EASTING[i]>q$xgrid[z] && lds_from_1980$EASTING[i]<q$xgrid[z+1]){
          som[j,z] = som[j,z]+ lds_from_1980$HEADELEV[i]
        }
      }
    }
  }
}
som

GenDat_pp2 <- function(df.ppp, dim, plotdat = TRUE,cova){
  
  # Define the window of interest
  win <- owin(c(0,dim[1]), c(0,dim[2]))
  
  # set number of pixels to simulate an environmental covariate
  spatstat.options(npixel=c(dim[1],dim[2]))
  
  gridcov <- cova
  
  # We count the number of points in each grid cell, plot it and make a vector out of it
  qcounts <- quadratcount(df.ppp, ny=dim[1], nx=dim[2])
 
  dens <- density(df.ppp)
  Lambda <- as.vector(t(qcounts)) # We have to use t() as we need to construct the vector with the column first
  
  if(plotdat == TRUE){
    par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1,0.5,0))
    plot(im(gridcov), main = 'Covariate')
    plot(dens, main = 'Intensity')
  }
  # Return a list with which I can play with
  return(list(Lambda = Lambda, pp = df.ppp, gridcov = gridcov))
}

pp2 = GenDat_pp2(lds.ppp,c(20,20),TRUE,som)

#stan model
ppm_stan <- '
data{
  int<lower = 1> n;
  vector[n] x;
  int<lower = 0> y[n];
}
parameters{
  real beta0;
  real beta1;
}
transformed parameters{
}
model{
  //priors
  target += normal_lpdf(beta0 | 0,5);
  target += normal_lpdf(beta1 | 0,10);

  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x);
}
generated quantities{
  vector[n] lambda_rep;
  lambda_rep = exp(beta0 + beta1 * x);
}'


stan_data = list(n = length(pp2$Lambda), x = as.vector(t(pp2$gridcov)) , y = pp2$Lambda)
fit_stan2 <- stan(model_code = ppm_stan, data = stan_data, 
                 warmup = 1000, iter = 5000, chains = 3)

print(fit_stan2, pars = c('beta0', 'beta1'))

lambda_rep2 <- as.data.frame(rstan::extract(fit_stan2)['lambda_rep'])
mean_lambda_rep2 <- apply(lambda_rep2, 2, 'mean')


grid <- st_make_grid(pp_sf, n = 20, what = 'polygons') %>% st_as_sf()
grid$pred_stan <- mean_lambda_rep2


plot(grid)

#this one gives a lot of warnings and beta1=0 even if i try to increase the number of iterations 

# Warning messages:
#   1: There were 1164 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
# https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 2: There were 3 chains where the estimated Bayesian Fraction of Missing Information was low. See
# https://mc-stan.org/misc/warnings.html#bfmi-low 
# 3: Examine the pairs() plot to diagnose sampling problems
# 
# 4: The largest R-hat is 3.83, indicating chains have not mixed.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#r-hat 
# 5: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#bulk-ess 
# 6: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#tail-ess 

#now i try to compute the new covariate matrix with the mean height for cell


mean_height=rep(0,400)

for(i in c(1:400)){
    if(pp$Lambda[i] != 0){
      mean_height[i]= as.vector(t(som))[i]/pp$Lambda[i]
    }  
}
mean_height


plot(im(mh_mat))

mh_mat=t(matrix(mean_height,nrow=20,ncol=20))

GenDat_pp2 <- function(df.ppp, dim, plotdat = TRUE,cova){
  
  # Define the window of interest
  win <- owin(c(0,dim[1]), c(0,dim[2]))
  
  # set number of pixels to simulate an environmental covariate
  spatstat.options(npixel=c(dim[1],dim[2]))
  
  gridcov <- cova
  
  # We count the number of points in each grid cell, plot it and make a vector out of it
  qcounts <- quadratcount(df.ppp, ny=dim[1], nx=dim[2])
  
  dens <- density(df.ppp)
  Lambda <- as.vector(t(qcounts)) # We have to use t() as we need to construct the vector with the column first
  
  if(plotdat == TRUE){
    par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1,0.5,0))
    plot(im(gridcov), main = 'Covariate')
    plot(dens, main = 'Intensity')
  }
  # Return a list with which I can play with
  return(list(Lambda = Lambda, pp = df.ppp, gridcov = gridcov))
}

pp3 = GenDat_pp2(lds.ppp,c(20,20),TRUE,mh_mat)

#stan model
ppm_stan <- '
data{
  int<lower = 1> n;
  vector[n] x;
  int<lower = 0> y[n];
}
parameters{
  real beta0;
  real beta1;
}
transformed parameters{
}
model{
  //priors
  target += normal_lpdf(beta0 | 0,5);
  target += normal_lpdf(beta1 | 0,10);

  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x);
}
generated quantities{
  vector[n] lambda_rep;
  lambda_rep = exp(beta0 + beta1 * x);
}'


stan_data = list(n = length(pp3$Lambda), x = as.vector(t(pp3$gridcov)) , y = pp3$Lambda)
fit_stan3 <- stan(model_code = ppm_stan, data = stan_data, 
                 warmup = 500, iter = 2500, chains = 3)

print(fit_stan3, pars = c('beta0', 'beta1'))

lambda_rep3 <- as.data.frame(rstan::extract(fit_stan3)['lambda_rep'])
mean_lambda_rep3 <- apply(lambda_rep3, 2, 'mean')


grid <- st_make_grid(pp_sf, n = 20, what = 'polygons') %>% st_as_sf()
grid$pred_stan <- mean_lambda_rep3


plot(grid)

#still gives warnings and beta1 = 0 and grid is the oppsite of what we expect



