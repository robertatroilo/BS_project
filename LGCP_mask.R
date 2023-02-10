library(sf)
library(rgdal)
library(sp)
library(spatstat)
library(rstan)
library(maptools)

fname <- "ENTLI.gdb.zip"

######## 1. LOADING IN THE DATA  ##########
# consider landslides since 2000
landslides = st_read(fname, "ENTLI_Crown")
lds_from_2000_with_NA = landslides[which(landslides$YEAR_1 > "2000"),]


######## 2. DATA  CLEANING ##########
#should add some filters here for duplicates or any other stuff 
u_NA = union(which(lds_from_2000_with_NA$M_WIDTH == 9999),which(lds_from_2000_with_NA$HEADELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$TAILELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$ELE_DIFF == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$SLOPE == 9999))
lds_from_2000=lds_from_2000_with_NA[-u_NA,]
lds = lds_from_2000

########## 2.5 chris code for utm to lon/lat ##########
#chris is insanely productive holy 
#from UTM to lat long

lat.long.df <- data.frame(lds$EASTING, lds$NORTHING) 
coordinates(lat.long.df) <-  ~lds.EASTING + lds.NORTHING
proj4string(lat.long.df)
proj4string(lat.long.df) <- CRS("+init=epsg:2326")

head(lat.long.df)
dist.location <- spTransform(lat.long.df, CRS("+init=epsg:4326"))
head(dist.location)

head(dist.location@coords)



######## 3. PRELIM SPATIAL OBJECTS  ##########
#super cool (thanks chris)
hk <- readShapePoly("./hk_boundaries/hh544xg3454.shp")
hkw = as.owin(hk)
hkw = owin(hkw$xrange, hkw$yrange) #the old window preserves boundary info but we want grid 
hkw

features.drop = c("SLIDE_TYPE", "YEAR_1", "ENTLI_NO") #read data dict 
feature.mask = !names(lds)%in%features.drop
lds.ppp = ppp(dist.location@coords[,1], dist.location@coords[,2], hkw, marks=lds[names(lds)[feature.mask]]) #use hk window
lds.ppp = unique(lds.ppp)


#choose a window size preserving the relative dimensions of hk 
win = hkw
n_quadr = 15
r = diff(win$yrange)/diff(win$xrange)

spatstat.options(npixel=c(n_quadr,ceiling(n_quadr*r)))
spatstat.options()$npixel

ref = ppp(lds.ppp$x, lds.ppp$y, win) 
q <- quadratcount(lds.ppp, ny=spatstat.options()$npixel[2], nx=spatstat.options()$npixel[1])


########## 4. generate mask from hk shapefile  #########
#chris again
grid <- makegrid(x=hk)
grid = spsample(hk, 2500,type="regular")
coords = grid@coords

hk_coords.ppp = ppp(coords[,1], coords[,2], hkw)
hk.q = quadratcount(hk_coords.ppp, ny=spatstat.options()$npixel[2], nx=spatstat.options()$npixel[1])

#define mask as points which are non-zero in the actual shapefile 
mask = which(as.vector(hk.q)!=0)
mask

q.real = as.vector(q)
q.extracted = as.vector(q)[mask]

loss = sum(q.real) - sum(q.extracted)
loss/sum(q.real)

length(q.real)
length(q.extracted)


#visualize 
x11()
par(mfrow=c(1,2))
plot(hk_coords.ppp, cex = 0.5, pch = "+", main="real map")
plot(hk.q, add=T)

plot(ref, cex = 0.5, pch = "+", main="landslide data")
plot(q, add=T)


########### 5. STAN implementation & modelling #########

#modified this a bit cause it wasn't working originally (lower time complexity tho)
reduce.covariates <- function(coords, data, q.ref, win.ref) {
  npixel = dim(q.ref)
  reduced = matrix(data = 0, nrow=npixel[1], ncol = npixel[2])
  counts = matrix(data = 0, nrow=npixel[1], ncol = npixel[2])
  
  for(i in c(1:length(data))) {
    col = ceiling(((coords[i,1]-win.ref$xrange[1])/diff(win.ref$xrange))*npixel[2])
    row = ceiling(((coords[i,2]-win.ref$yrange[1])/diff(win.ref$yrange))*npixel[1])
    
    reduced[row, col] = (reduced[row, col]*counts[row,col] + data[i])/(counts[row,col] +1)
    counts[row, col] = counts[row, col] + 1
  }
  
  return(reduced);
}

#covariates 
slope <- reduce.covariates(data.frame(lds.ppp$x, lds.ppp$y), lds.ppp$marks$SLOPE, q, hkw)
ele.diff <- reduce.covariates(data.frame(lds.ppp$x, lds.ppp$y), lds.ppp$marks$ELE_DIFF, q, hkw)
width <- reduce.covariates(data.frame(lds.ppp$x, lds.ppp$y), lds.ppp$marks$M_WIDTH, q, hkw)
length <- reduce.covariates(data.frame(lds.ppp$x, lds.ppp$y), lds.ppp$marks$S_LENGTH, q, hkw)
head.elevation <- reduce.covariates(data.frame(lds.ppp$x, lds.ppp$y), lds.ppp$marks$HEADELEV, q, hkw)

#we can either treat aggregated binary features as continuous through averaging or encode over area
#through voting
nocover = ifelse(lds.ppp$marks$COVER=="A", 1, 0)
nocover = reduce.covariates(data.frame(lds.ppp$x, lds.ppp$y), nocover, q, hkw)
nocover.vote = ifelse(nocover==0, 0, 1)

#apparently na's here 
gully = reduce.covariates(data.frame(lds.ppp$x, lds.ppp$y), ifelse(lds.ppp$marks$GULLY=="N", 0, 1), q, hkw)
gully.vote = ifelse(gully==0, 0, 1)

########### 5.5 Distance Matrix derivation and such ###########

#distance matrix from centroids of the grid cells (will need to convert coords to km so just use relative)
coords.naive = expand.grid(seq(1:dim(q)[1]), seq(1:dim(q)[2]))

#using the mask here - expands column wise not row by row 
DMat = as.matrix(dist(coords.naive[mask,], method="euclidean", diag=T, upper=T))



ppm_stan3 <- '
// Fit an accurate LGCP in Stan with the Exponential covariance structure
functions{
  
    matrix GP(matrix D, real sigma_sq, real scale, real delta) {
        int N = dims(D)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sigma_sq + delta;
          for (j in (i + 1):N) {
            K[i, j] = sigma_sq * exp(- D[i,j] / scale );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sigma_sq + delta;
        return K;
    }
}
data{
  int<lower = 1> N;
  vector[N] x;
  int<lower = 0> y[N];
  matrix[N, N] DMat; // Distance matrix
}
parameters{
  real beta0;
  real beta1;
  vector[N] k;
  
  // GP standard deviation parameters
  real<lower=0> sigma_sq;
  // GP length-scale parameters
  real<lower=0> scale;
}
model{
  matrix[N,N] SIGMA;
  vector[N] mu;

  SIGMA = GP(DMat, sigma_sq, scale, 0.01);
  k ~ multi_normal(rep_vector(0,N), SIGMA);
  
  //priors for the coefficients
  target += normal_lpdf(beta0 | 0,5);
  target += normal_lpdf(beta1 | 0,10);
  
  // Prior for the noise
  target += cauchy_lpdf(sigma_sq | 0, 1);
  target += inv_gamma_lpdf(scale | 3.884416, 0.77454);

  // likelihood
    for(i in 1:N){
    mu[i] = beta0 + beta1 * x[i] + k[i];
  }
  
  target += poisson_log_lpmf(y | mu);
}
generated quantities{
}
'


#extracting relevant cells
slope = as.vector(slope)[mask]
ele.diff = as.vector(ele.diff)[mask]
nocover.vote = as.vector(nocover.vote)[mask]
width = as.vector(width)[mask]
length = as.vector(length)[mask]
head.elevation = as.vector(head.elevation)[mask]
gully.vote = as.vector(gully.vote)[mask]

stan_data3 = list(N = N, x = slope, y = q.extracted, 
                  DMat=DMat)


stan_fit0 <- stan(model_code = ppm_stan3,
                  data = stan_data3,
                  chains = 1, warmup = 1000, iter = 8000)#,
# control = list(adapt_delta = 0.999, max_treedepth=13))

print(stan_fit0, pars = c('beta0', 'beta1', 'sigma_sq', 'scale'))

beta0 = rstan::extract(stan_fit0)['beta0']
beta0 = mean(beta0_rep3$beta0)

beta1 = rstan::extract(stan_fit0)['beta1']
beta1 = mean(beta1_rep3$beta1)

# Get the posterior of the parameters
draws <- rstan::extract(stan_fit0, pars = c('beta0', 'beta1', 'sigma_sq', 'scale'))

# We make a sequence of distance
dist_seq <- seq(from = min(DMat), to = max(DMat), length.out = 100)

# Compute the mean and the standard deviation of the posterior correlation
post_cov <- sapply(dist_seq,function(x)draws$sigma_sq*exp(-draws$scale*x^2))
post_cov_mu <-apply(post_cov,2,mean)
post_cov_sd <-apply(post_cov,2,sd)

# Make a dataframe and plot
library(tibble)
post_df <- tibble(dist = dist_seq,
                  mu = post_cov_mu,
                  sd = post_cov_sd)

library(ggplot2)
ggplot(post_df, aes(x = dist)) +
  geom_line(aes(y = mu), color = "#CD295A", size = 1) +
  geom_ribbon(aes(ymin = mu - sd, ymax = mu + sd), fill = "#38ADAE", alpha = .3) +
  theme_classic() +
  ylab("Covariance") +
  xlab("Distance")

#We can see that the correlation becomes null at a distance close to 2 units

#Extract the predicted point pattern

k = as.data.frame(rstan::extract(stan_fit0)['k'])
k <- apply(k, 2, 'mean')
pred = exp(beta0 + beta1 * slope + k)
#pred = exp(beta0 + beta1 * slope )
#####
df_quadr = data.frame(as.matrix(q))[mask,]
df_pred =  data.frame(pred,x = df_quadr$x, y = df_quadr$y) 
library(tidyr)

#code to remove "[" and ")"
c = df_pred$x[1]
c = gsub("\\[|\\]", "", c)
c = gsub("[()]", "", c)
scan(textConnection(c), sep=",")
###

x_pred = c(1:length(df_pred$pred))
y_pred = c(1:length(df_pred$pred))

for (i in c(1:length(df_pred$pred))){
  x_pred[i] = mean(scan(textConnection(gsub("[()]", "", gsub("\\[|\\]", "", df_pred$x[i]))), sep=","))
  y_pred[i] = mean(scan(textConnection(gsub("[()]", "", gsub("\\[|\\]", "", df_pred$y[i]))), sep=","))
}

df = data.frame(value = pred,x = x_pred, y = y_pred)

npix <- nrow(df)

#each pixel has a probability to have points, weighted with its lambda value
lpix <- df$value

#number of total points drawn from a poisson
mu = sum(df$value) #this number generates too much points equal to the one used in the training
mu = mu
ni <- rpois(1, mu)

#sample indexes from df, each index is sampled a number of times proportional
#to the lambda value associated to its pixel
ii <- sample.int(npix, size=ni, replace=TRUE, prob=lpix/mu)

#Radius of the cells, to give "uncertainty" to the pattern
#points are in a random position within each cell
dx <- abs(114.28-114.32)/2
dy <- (abs(22.3-22.34)-0.005)/2

#generate coordinates of the points in the point pattern
#the second term gives a random position in the cell, around its centroid
xx <- df$x[ii] + runif(ni, -dx, dx)
yy <- df$y[ii] + runif(ni, -dy, dy)


hk <- readShapePoly("sf/HKG_adm0.shp")#file found in "hk_boundaries.zip"
hkwindow = as.owin(hk)
pred_pp <- ppp(xx, yy, window=hkwindow)#, check=FALSE)
Window(pred_pp) = hkwindow

#compare to real pattern
#Real pattern
X  <- ref  
X <- as(X, "ppp")
marks(X) <- NULL
Window(X) = hkwindow

x11()
par(mfrow=c(2,2))
plot(density(pred_pp),main="Pred Density")
plot(density(X),main="Real Density")
plot(pred_pp,main="Pred Pattern",pch='.')
plot(X,main="Real Pattern",pch='.')

plot(pred_pp,main="Pred Pattern",pch='+')
plot(X,main="Real Pattern",pch='+')







