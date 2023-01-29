#Load libraries
#----
library(sf)
library(rgdal)
library(sp)
library(spatstat)
library(rstan)
library(maptools)

fname <- "ENTLI.gdb.zip"
fname
st_layers(fname)

landslides = st_read(fname, "ENTLI_Crown")

# we are considering all the landslide to create a new dataset with latitude
lds_from_2000_with_NA = landslides
# data cleaning
u_NA = union(which(lds_from_2000_with_NA$M_WIDTH == 9999),which(lds_from_2000_with_NA$HEADELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$TAILELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$ELE_DIFF == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$SLOPE == 9999))
lds_from_2000=lds_from_2000_with_NA[-u_NA,]
lds = lds_from_2000
all_lds = lds

#map boundaries
hk <- readShapePoly("sf/HKG_adm0.shp")
hkw = as.owin(hk)

#from UTM to lat long
lat.long.df <- data.frame(lds$EASTING, lds$NORTHING) 
coordinates(lat.long.df) <-  ~lds.EASTING + lds.NORTHING
proj4string(lat.long.df)
proj4string(lat.long.df) <- CRS("+init=epsg:2326")
head(lat.long.df)
dist.location <- spTransform(lat.long.df, CRS("+init=epsg:4326"))
dist.location

lds_map <- 
  data.frame(all_lds[,-c(14,15)],
             lat = dist.location$lds.EASTING,
             long = dist.location$lds.NORTHING)


lds_map$Shape = NULL
lds_map <- st_as_sf(lds_map, coords = c(14:15))
lds_map <- lds_map[,c(-13)]

xx <- unlist(lds_map$geometry)

long <- xx[seq(1,length(xx),2)]

lat <- xx[-seq(1,length(xx),2)]

lds_map = cbind(lds_map,cbind(long,lat))

save(lds_map,file="Landslides.Rda")


#Load data
#----
landslides = load("Landslides.Rda")
landslides = lds_map

#they were saved as " A" instead of "A"
landslides$COVER[4616] = 'B'
landslides$COVER[5406] ='A'

dim(landslides)
# 106804 points

# consider landslides from 2000
lds_from_2000 = landslides[which(landslides$YEAR_1 > "2000"),]
lds = lds_from_2000
dim(lds)
# 6588 points

#map boundaries
hk <- readShapePoly("sf/HKG_adm0.shp")
hkw = as.owin(hk)


#this is the default option and is a good approximation
#then we can consider lower dimension for the coovariates approximation
spatstat.options(npixel=c(128,128))
#PPP object with map boundaries
X  <- as(lds, "Spatial")  
#X <- as(p.sp, "ppp") [hkw]
X <- as(X, "ppp")
marks(X) <- NULL
Window(X) = hkw

#Plot point pattern density and tesselation
plot(density(X))
plot(X,add=TRUE,pch='.')
plot(hkw, add=TRUE)

v <- tess(image = X)#,window=hkw)
plot(v)
plot(hkw, add=TRUE)

#PPP object with MARKS (covariates)
X2  <- as(lds, "Spatial")  
#X <- as(p.sp, "ppp") [hkw]
X2 <- as(X2, "ppp")
Window(X2) = hkw
plot(X2,which.marks="ELE_DIFF", cols = c("gray","red","blue","green"))
plot(X2,which.marks="COVER", cols = c("green","red","gray","blue"))

#PP object with 1 Mark: elev_diff
X3 = X
marks(X3) = X2$marks$ELE_DIFF
plot(Smooth(X3,sigma=0.1),main="Elevation Difference")
plot(X3,which.marks=NULL,add=TRUE,pch='.')
plot(hkw, add=TRUE)
#PP object with 1 Mark: HEADELEV
X4 = X
marks(X4) = X2$marks$HEADELEV
plot(Smooth(X4,sigma=0.05),main="HEADELEV")
plot(X4,which.marks=NULL,add=TRUE,pch='.')
plot(hkw, add=TRUE)

#PP object with 1 Mark: SLOPE
X4 = X
marks(X4) = X2$marks$SLOPE
plot(Smooth(X4,sigma=0.015),main="SLOPE")
plot(X4,which.marks=NULL,add=TRUE,pch='.')
plot(hkw, add=TRUE)

#PP object with 1 categorical Mark: GULLY
gully = split.ppp(X2,"GULLY")
gully$Y$n #few samples
gully$N$n
marks(gully$N) = NULL
marks(gully$Y) = NULL
x11()
par(mfrow = c(1,2))
plot(density(gully$N, main = "GULLY N"))
plot(hkw, add=TRUE)
plot(gully$N, add=TRUE,pch='.')
plot(density(gully$Y, main = "GULLY Y"))
plot(hkw, add=TRUE)
plot(gully$Y, add=TRUE,pch='.')

#PP object with 1 categorical Mark: Cover
Cover = split.ppp(X2,"COVER")
plot(density(Cover),main="COVER")
Cover$A$n
Cover$B$n
Cover$C$n 
Cover$D$n 


#number of quadrants
n_quadr = 30

spatstat.options(npixel=c(128,128))
print((X$window))
dim <- rbind(c(113.83459, 114.44097), c(22.153194, 22.562094))
win <- owin(c(dim[1,1],dim[1,2]), c(dim[2,1],dim[2,2])) #RECTANGLE
#??win = X$window #POLYGON

#inverted x0 and y0 because of correspondence between long lat and easting northing
y0 <- seq(win$yrange[1], win$yrange[2],
          length=n_quadr)
x0 <- seq(win$xrange[1], win$xrange[2],
          length=n_quadr)

qcounts_real2 <- quadratcount(X, ny=n_quadr, nx=n_quadr)
real_lambda_polygon = as.vector(t(qcounts_real2))

#Quadratcount visualization example
Q <- quadratcount(X, nx= 5, ny=5)
plot(Q)
plot(intensity(Q, image=TRUE))
plot(Q,add=TRUE)

######
#Here we should create the covariates but how can we create them in
#a irregular polygon ? I tried without success, i procede as always
#
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
          } } } } }
  return(covar)
}

covar1 <- define_covariate(northing = lds$lat, easting = lds$long, feature = lds$M_WIDTH)
covar2 <- define_covariate(northing = lds$lat, easting = lds$long, feature = lds$HEADELEV)
covar3 <- define_covariate(northing = lds$lat, easting = lds$long, feature = lds$TAILELEV)
covar4 <- define_covariate(northing = lds$lat, easting = lds$long, feature = lds$ELE_DIFF)
covar5 <- define_covariate(northing = lds$lat, easting = lds$long, feature = lds$SLOPE)
######

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
  target += uniform_lpdf(sigma_noise | 0,1);
  
  // Prior for the noise
  target += normal_lpdf(noise | 0, sigma_noise);
  
  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x1 + beta2 * x2 + noise );
}
generated quantities{
  vector[n] lambda_rep;
  lambda_rep = exp(beta0 + beta1 * x1 + beta2 * x2  + noise);
}'

#######
#Stan simulation
#Can't manage the polygon window instead of the rectangular we have usually 
#used in the stan simulation
b0 <- 0.1
b2 <- 0.2
b5 <- 0.05

stan_data3 = list(n = length(real_lambda), x1 = as.vector(t(covar4)), 
                  x2 = as.vector(t(covar5)), y = real_lambda, 
                  b0 = b0, b1 = b2, b2 = b5)

fit_stan3 <- stan(model_code = ppm_stan3, data = stan_data3, 
                  warmup = 3000, iter = 10000, chains = 4)


print(fit_stan3, pars = c('beta0', 'beta1', 'beta2','sigma_noise'))
plot(fit_stan3, pars = c('beta0', 'beta1', 'beta2','sigma_noise'))

lambda_rep3 <- as.data.frame(rstan::extract(fit_stan3)['lambda_rep'])
mean_lambda_rep3 <- apply(lambda_rep3, 2, 'mean')
m = matrix(mean_lambda_rep3,n_quadr,n_quadr)

library(OpenImageR)
par(mfrow=c(1,2),mar=c(1,1,1,1))
plot(im(flipImage(as.matrix(im(mean_lambda_rep3,xcol=x0,yrow=y0)))),main="Pred Lambda")
plot(im(flipImage(as.matrix(im(real_lambda,xcol=x0,yrow=y0)))),main="Real Lambda")

rmse <- function(true, observed){sqrt((true - observed)^2)}
sum(rmse(real_lambda, mean_lambda_rep3))

#######
#Noise is n_quadr x n_quadr
#Better to compute mean between the iteration but considering each cell separately
noise <- as.data.frame(rstan::extract(fit_stan3)['noise'])
mean_noise <- apply(noise, 2, 'mean')

#These are the coefficients of the best simulation (with ele_diff and slope)
#The covariates were defined on the rectangular grid 30x30
beta0   =     0.55    
beta1    =   -0.02   
beta2     =   0.01  
mean_noise = mean_noise



#How can we create the new covariates grids that should stay in the polygon 
#region of hong kong
#Even if we create rectangular grid, we don't have the value of the covariates
#over all the map but just in the landslide points recorded

#Function to return the Lambda intensity fitted with stan to predict new
#pattern 
#I use the training rectangular covariates)
lfun.pred <- function(x,y){ 
  exp(beta1 * x + beta0 )}
m_lam = matrix(exp(beta1 * covar4+beta2*covar5 + beta0 + mean_noise), n_quadr)

#In this way we can recrate the point pattern from the lambda (but always 
#in a rectangular grid)
pred_pp <- rpoispp(im(flipImage(as.matrix(im(m_lam)))), win = hkw)
plot(density(pred_pp))
#Window(ppdat) = hkw #polygon window doesn't work


#Trying to fit the fake predictions of the training in a polygon window
#but without success
a = covar4
b = covar5
z = mean_noise
gen.sample.fitted <- function() {
  lfun <- function(a, b) exp(beta0 + beta1 * a + beta2 * b + z)
  pp <- rpoispp(lfun, win = hkw)
  data.frame(x = pp$x, y = pp$y)
}
lfun <- function(a, b) exp(beta0 + beta1 * a + beta2 * b + 12)
pp <- rpoispp(lfun, win = hkw)

q = gen.sample.fitted()
