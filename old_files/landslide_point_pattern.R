library(sf)
library(rgdal)
fname <- "ENTLI.gdb.zip"
fname

#fname <- "HLC (Up to Year 2016).gdb" #CATHCMENT
#fname


st_layers(fname)

library(sp)
methods(st_as_sf)
methods(st_as_sfc)

landslides = st_read(fname, "ENTLI_Crown")
# Crown of a landslide: The practically undisplaced material still in place and 
#                       adjacent to the highest parts of the main scarp.
#plot(landslides, max.plot = 14) # to plot all the 14 features (fields)

trail = st_read(fname, "ENTLI_Trail")

lds_from_2000 = landslides[which(landslides$YEAR_1 > "2000"),]
  
library(ggplot2)
#x11()
#ggplot(landslides) + geom_sf()

ggplot(lds_from_2000) + geom_sf()


lds = lds_from_2000

lds_relict <- lds[which(lds$SLIDE_TYPE == "R"),]
dim(lds_relict) # 0

dim(landslides)[1]/dim(landslides)[1] # --> 80% of the recorded landslides is a relict landslide
# indeed
dim(landslides[which(landslides$YEAR_1 > "1975"),]) # == 0 
# relict lanslides are before 1975


lds_flow <- lds[which(lds$SLIDE_TYPE == "C"),]
dim(lds_flow) # 
x11()
ggplot(lds_flow) + geom_sf()

lds_hill <- lds[which(lds$SLIDE_TYPE == "O"),]
dim(lds_hill) # 
x11()
ggplot(lds_hill) + geom_sf()

lds_coastal <- lds[which(lds$SLIDE_TYPE == "S"),]
dim(lds_coastal) # 
x11()
ggplot(lds_coastal) + geom_sf()



ggplot(lds,aes(fill=HEADELEV)) + geom_sf() 


plot(lds$EASTING,lds$NORTHING, col = lds$CLASS)


lds_from_2000 = landslides[which(landslides$YEAR_1 > "2000"),]

lds = lds_from_2000
#POINT PATTERN
library(maptools)#
library(spatstat)
p.sp  <- as(lds, "Spatial")  # Create Spatial* object
lds.ppp <- as(p.sp, "ppp")      # Create ppp object


#p.sf.utm <- st_transform(lds, 32619) # project from geographic to UTM
#p.sp  <- as(p.sf.utm, "Spatial")      # Create Spatial* object
#p.ppp <- as(p.sp, "ppp")              # Create ppp object
class(lds.ppp)


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

pp = GenDat_pp(lds.ppp,c(20,20),TRUE)




#Stan
library(rstan)

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
  target += normal_lpdf(beta0 | 0,10);
  target += normal_lpdf(beta1 | 0,20);

  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x);
}
generated quantities{
  vector[n] lambda_rep;
  lambda_rep = exp(beta0 + beta1 * x);
}'

stan_data = list(n = length(pp$Lambda), x = as.vector(t(pp$gridcov)), y = pp$Lambda)
fit_stan <- stan(model_code = ppm_stan, data = stan_data, 
                 warmup = 10, iter = 100, chains = 3)
#fit_stan <- stan(file = "ppm_stan.stan", data = stan_data, 
#                warmup = 500, iter = 2000, chains = 3)
