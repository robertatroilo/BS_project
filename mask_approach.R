library(sf)
library(rgdal)
library(sp)
library(spatstat)
library(rstan)
library(maptools)

fname <- "../ENTLI.gdb.zip"

######## 1. LOADING IN THE DATA  ##########
# consider landslides since 2000
landslides = st_read(fname, "ENTLI_Crown")

#c("2000", "2001", "2002")
lds_from_2000_with_NA = landslides[which(landslides$YEAR_1 %in% c("2000", "2001", "2002")),]


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
n_quadr = 20
r = diff(win$yrange)/diff(win$xrange)

spatstat.options(npixel=c(ceiling(n_quadr*r), n_quadr))  #n_quadr is ref for number of cols/x range 
spatstat.options()$npixel

ref = ppp(lds.ppp$x, lds.ppp$y, win) 
q <- quadratcount(lds.ppp, ny=spatstat.options()$npixel[1], nx=spatstat.options()$npixel[2])
dim(q)

plot(ref)
plot(q, add=T)

########## 4. generate mask from hk shapefile  #########
#chris again
grid <- makegrid(x=hk)
grid = spsample(hk, 2500,type="regular")
coords = grid@coords

hk_coords.ppp = ppp(coords[,1], coords[,2], hkw)
hk.q = quadratcount(hk_coords.ppp, ny=spatstat.options()$npixel[1], nx=spatstat.options()$npixel[2])

#define mask as points which are non-zero in the actual shapefile 
mask = which(as.vector(t(hk.q), )!=0)
mask

q.real = as.vector(t(q))
q.extracted = as.vector(t(q))[mask]

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
    #convert coordinate to (row,col) indices
    row = ceiling(((coords[i,2]-win.ref$yrange[1])/diff(win.ref$yrange))*npixel[1])
    col = ceiling(((coords[i,1]-win.ref$xrange[1])/diff(win.ref$xrange))*npixel[2])
    
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

###############################################################
########### 5.5 Distance Matrix derivation and such ###########
###############################################################

#let's actually make the coordinates 
x0 <- seq(hkw$xrange[1], hkw$xrange[2], length=spatstat.options()$npixel[2])
y0 <- seq(hkw$yrange[1], hkw$yrange[2], length=spatstat.options()$npixel[1])

print(length(x0))

#the result of this is expanded row by row for the coords -- need consistency for the mask
coords.actual = expand.grid(x = x0, y = y0)

#with centroids?
x0_extended = seq(hkw$xrange[1], hkw$xrange[2], length=spatstat.options()$npixel[2]+1)
y0_extended <- seq(hkw$yrange[1], hkw$yrange[2], length=spatstat.options()$npixel[1]+1)

x0_centered = x0_extended[1:20]+ 0.5*diff(x0_extended)
y0_centered = y0_extended[1:20]+ 0.5*diff(y0_extended)

coords.centered = expand.grid(x = x0_centered, y = y0_centered)

#distance matrix from centroids of the grid cells (will need to convert coords to km so just use relative)
coords.naive = expand.grid(seq(1:dim(q)[2]), seq(1:dim(q)[1]))

#using the mask here - expands column wise not row by row 
DMat = as.matrix(dist(coords.naive[mask,], method="euclidean", diag=T, upper=T))

#now derive the CAR matrix we employ (from notes on geostats)
dmat.quantile = quantile(as.vector(DMat), 0.1)
A = ifelse(DMat<=dmat.quantile, 1, 0)
diag(A) = 0

Dw = diag(rowSums(A))

rho = 0.5
W.inv = Dw - rho*A
W = solve(W.inv)

det(W%*%W.inv) #verify 


#stan model 
ppm_stan <- '
data{
  int<lower = 1> N;
  int<lower = 1> p;
  matrix[N, p] X; //covariates
  matrix[N,N] W; //CAR model 
  int<lower = 0> y[N];
}
parameters{
  vector[p] beta;
  vector<lower=0>[p] sigma_beta;
  
  // W standard deviation parameter
  vector[N] e;
  real<lower=0> sigma_sq;
}

transformed parameters {
    vector[N] mu;
    for(i in 1:N) {
      mu[i] = exp(row(X, i) * beta + e[i]); //link 
    }
}

model{
  for (s in 1:N) {
        y[s] ~ poisson(mu[s]);  
  } 
  
  //e ~ multi_normal(rep_vector(0,N), sigma_sq*W); 
  e ~ normal(0, sigma_sq);
    
  //prior for the noise on CAR 
  sigma_sq ~ inv_gamma(2,2);
  
  //priors for the coefficients
  for (j in 1:p) {
        beta[j] ~ normal(0.0, sigma_beta);
        sigma_beta[j] ~ inv_gamma(2, 2);
    }
}
generated quantities {
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

#normalize 
slope.norm = (slope - mean(slope))/sqrt(var(slope))
width.norm = (width - mean(width))/sqrt(var(width))
head.elevation.norm = (head.elevation - mean(head.elevation))/sqrt(var(head.elevation))
ele.diff.norm = (ele.diff - mean(ele.diff))/sqrt(var(ele.diff))
length.norm = (length - mean(length))/sqrt(var(length))


#fitting (and praying)
#normalized covariates
X = cbind(rep(1,length(slope.norm)), slope.norm, width.norm, head.elevation.norm,
          ele.diff.norm, length.norm, nocover.vote, gully.vote)

X = cbind(rep(1,length(slope)), slope, width, head.elevation, ele.diff, length, nocover.vote, gully.vote)

N = dim(X)[1]
p = dim(X)[2]

stan_data = list(N = N, p = p, X = X, y = q.extracted, W=W)

stan_fit0 <- stan(model_code = ppm_stan,
                  data = stan_data,
                  chains = 1, warmup = 1000, iter = 16000,
                  control = list(adapt_delta = 0.999, max_treedepth=13))




############# plotting the chains ###############3
beta.chain = rstan::extract(stan_fit0)['beta']
beta.hat = colMeans(beta.chain$beta)

sigma = rstan::extract(stan_fit0)['sigma_sq']
mean(sigma$sigma_sq)

#saving these since they took forever
#save(beta.chain, file="../Rdata/CARModel_beta.Rda")
#save(sigma, file="../Rdata/CARModel_sigma.Rda")

#recall this is a tiled surface which should in theory predict the realization in our data 
#no need to be exactly the same, but makes things hard to analyze 
pred = exp(X%*%beta.hat) n

x11()
par(mfrow=c(1,2))
plot(seq(1,N), log(q.extracted), pch = 19, col="red")
plot(seq(1,N), log(pred), pch = 19, col="blue")


x11()
par(mfrow=c(4,2))
plot(density(beta.chain$beta[,1]), main="", xlab="Beta0")
plot(density(beta.chain$beta[,2]), main="", xlab="Beta1")
plot(density(beta.chain$beta[,3]), main="", xlab="Beta2")
plot(density(beta.chain$beta[,4]), main="", xlab="Beta3")
plot(density(beta.chain$beta[,5]), main="", xlab="Beta4")
plot(density(beta.chain$beta[,6]), main="", xlab="Beta5")
plot(density(beta.chain$beta[,7]), main="", xlab="Beta6")
plot(density(beta.chain$beta[,8]), main="", xlab="Beta7")


plot(density(sigma$sigma), main="", xlab="sigma")

i=4
plot(seq(1,length(beta.chain$beta[,i])), beta.chain$beta[,i], type='l', xlab="iter", ylab='val')


