source("./mask_generation.R")

fname <- "../../ENTLI.gdb.zip"

######## 1. LOADING IN THE DATA  ##########
# consider landslides since 2000
landslides = st_read(fname, "ENTLI_Crown")


#YEARS = c("2016", "2017", "2018")
YEARS= c("2019")
lds_from_2000_with_NA = landslides[which(landslides$YEAR_1 %in% YEARS),]
#lds_from_2000_with_NA = landslides[which(landslides$YEAR_1 > "2000"),]


######## 2. DATA  CLEANING ##########
u_NA = union(which(lds_from_2000_with_NA$M_WIDTH == 9999),which(lds_from_2000_with_NA$HEADELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$TAILELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$ELE_DIFF == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$SLOPE == 9999))
u_NA = union(u_NA,is.na(lds_from_2000_with_NA$GULLY))

lds_from_2000=lds_from_2000_with_NA[-u_NA,]
lds = lds_from_2000

#quick annual contribution plot 
#years = as.factor(lds$YEAR_1)
#x11()
#table(years)
#barplot(table(years), xlab="Year", ylab="recorded landslides")


########## 2.5 chris code for utm to lon/lat ##########
#chris is insanely productive holy 
#from UTM to lat long
lat.long.df <- data.frame(lds$EASTING, lds$NORTHING) 
coordinates(lat.long.df) <-  ~lds.EASTING + lds.NORTHING
proj4string(lat.long.df)
proj4string(lat.long.df) <- CRS("+init=epsg:2326")

head(lat.long.df)
dist.location <- spTransform(lat.long.df, CRS("+init=epsg:4326"))

head(dist.location@coords)


############ create ppp object using extracted coords & our data ################
features.drop = c("YEAR_1", "ENTLI_NO") #read data dict 
feature.mask = !names(lds)%in%features.drop
lds.ppp = ppp(dist.location@coords[,1], dist.location@coords[,2], hkw, marks=lds[names(lds)[feature.mask]]) #use hk window
lds.ppp = unique(lds.ppp)


ref = ppp(lds.ppp$x, lds.ppp$y, hkw) #using the hk window from mask_generation 
q <- quadratcount(lds.ppp, ny=spatstat.options()$npixel[1], nx=spatstat.options()$npixel[2])
dim(q)

plot(ref)
plot(q, add=T)

#complete spatial randomness? -- not really valid since some cells aren't even hong kong 
quadrat.test(q)

########## sample reduction; we can either use hk mask OR only keep cells which are in the sample \
#before we kept all the samples which were inside HK but instead we should really only consider our samples
#idk if this is a little biased or what will happen 

mask.data = which(as.vector(t(q))!=0)

q.real = as.vector(t(q))
q.extracted = as.vector(t(q))[mask.data]

loss = sum(q.real) - sum(q.extracted)
loss/sum(q.real)

length(q.real)
length(q.extracted)


########### 5. STAN implementation & modelling #########

reduce.covariates <- function(coords, data, q.ref, win.ref) {
  npixel = dim(q.ref)
  reduced = matrix(data = 0, nrow=npixel[1], ncol = npixel[2])
  counts = matrix(data = 0, nrow=npixel[1], ncol = npixel[2])
  
  for(i in c(1:length(data))) {
    #convert coordinate to (row,col) indices
    row = npixel[1] - ceiling(((coords[i,2]-win.ref$yrange[1])/diff(win.ref$yrange))*npixel[1]) + 1
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
vegetation = ifelse(lds.ppp$marks$COVER!="A", 1, 0)
vegetation = reduce.covariates(data.frame(lds.ppp$x, lds.ppp$y), vegetation, q, hkw)
vegetation.vote = ifelse(vegetation==0, 0, 1)

#apparently na's here - treat this in the data cleaning when deciding to use this 
erosion = reduce.covariates(data.frame(lds.ppp$x, lds.ppp$y), ifelse(lds.ppp$marks$GULLY=="N", 0, 1), q, hkw)
erosion.vote = ifelse(erosion==0, 0, 1)

###############################################################
########### 5.5 Distance Matrix derivation and such ###########
###############################################################

#extract the grid coords 
tess = as.tess(q) #creates tessellation object from quadratcount
x0 = tess$xgrid #this is the n+1 points in which the grid cells lie
y0 = tess$ygrid

x = x0[1:length(x0)-1]+diff(x0)
y = y0[1:length(y0)-1]+diff(y0)

coords.tess = expand.grid(x,y)

#using the mask here - expands column wise not row by row 
DMat = as.matrix(dist(coords.tess[mask.data,], method="euclidean", diag=T, upper=T))

#k = 5+1 #knn approach 
#for(i in 1:dim(DMat)[1]){
#  where.condition = (DMat[i,]%in%sort(DMat[i,])[1:4])
#}

########## CAR ##################
#now derive the CAR matrix we employ (from notes on geostats) -> this isnt true adjacency but its fine 
dmat.quantile = quantile(as.vector(DMat), 0.1)
A = ifelse(DMat<=dmat.quantile, 1, 0)
diag(A) = 0

Dw = diag(rowSums(A))

rho = 0.5
W.inv = Dw - rho*A
W = solve(W.inv)

det(W%*%W.inv) #verify 
#################################

#stan model 
ppm_stan <- '
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
  int<lower = 1> p;
  matrix[N,p] X; //covariates
  matrix[N,N] DMat;
  int<lower = 0> y[N];
}
parameters{
  //regression coefficients 
  vector[p] beta;
  real<lower=0> nu_sq; //changed it to one variance for the regression coefs 
  
  //independent noise
  vector[N] eps;
  real<lower=0> tau_sq;
  
  //spatial noise 
  vector[N] w;
  real<lower=0> sigma_sq;
  real<lower=0> scale;
}

transformed parameters {
    vector[N] mu;
    for(i in 1:N) {
      mu[i] = exp(row(X, i) * beta + w[i] + eps[i]); //link 
    }
}

model{
  //spatial noise 
  matrix[N,N] SIGMA;
  SIGMA = GP(DMat, sigma_sq, scale, 0.01);
  w ~ multi_normal(rep_vector(0,N), SIGMA);
  
  for (s in 1:N) {
        y[s] ~ poisson(mu[s]);  
  } 
  
  //independent noise 
  eps ~ normal(0, tau_sq); 
  
  //priors for the coefficients
  beta ~ normal(0, nu_sq);
    
  //lpdfs for the tau^2, sigma^2, scale, nu^2 (4)
  target += cauchy_lpdf(sigma_sq | 0, 1); //tutorial - HC noise 
  target += inv_gamma_lpdf(scale | 3.884416, 0.77454); //tutorial
  target += inv_gamma_lpdf(tau_sq | 1, 1); //me
  target += inv_gamma_lpdf(nu_sq | 1, 1); //me
}
generated quantities {
  }
'

#extracting relevant cells
slope = as.vector(t(slope))[mask.data]
ele.diff = as.vector(t(ele.diff))[mask.data]
vegetation.vote = as.vector(t(vegetation.vote))[mask.data]
width = as.vector(t(width))[mask.data]
length = as.vector(t(length))[mask.data]
head.elevation = as.vector(t(head.elevation))[mask.data]
erosion.vote = as.vector(t(erosion.vote))[mask.data]

#normalize 
slope.norm = (slope - mean(slope))/sqrt(var(slope))
width.norm = (width - mean(width))/sqrt(var(width))
head.elevation.norm = (head.elevation - mean(head.elevation))/sqrt(var(head.elevation))
ele.diff.norm = (ele.diff - mean(ele.diff))/sqrt(var(ele.diff))
length.norm = (length - mean(length))/sqrt(var(length))


#fitting (and praying)
#normalized covariates
X = cbind(rep(1,length(q.extracted)), slope.norm, head.elevation.norm, vegetation.vote, erosion.vote)

#non-normalized covariates
#X = cbind(rep(1,length(q.extracted)), slope, width, head.elevation, ele.diff, length, vegetation.vote, erosion.vote)

N = dim(X)[1]
p = dim(X)[2]

stan_data = list(N = N, p = p, X = X, y = q.extracted, DMat = DMat)

stan_fit0 <- stan(model_code = ppm_stan,
                  data = stan_data,
                  chains = 1, warmup = 1000, iter = 8000,
                  control = list(adapt_delta = 0.999, max_treedepth=13))


#saving the model fit so we can relaod later and extract everything
#saveRDS(stan_fit0, file = "./expkernel30.RData")
stan_fit0 = readRDS("./expkernel30.RData")


###################################################################
################# POSTERIOR STUFF #################################
###################################################################

beta.chain = rstan::extract(stan_fit0)['beta']
beta.hat = colMeans(beta.chain$beta)

#sigma
sigma_sq = rstan::extract(stan_fit0)['sigma_sq']
sigma_sq.mean = mean(sigma_sq$sigma_sq)

#tau
tau_sq = rstan::extract(stan_fit0)['tau_sq']
tau_sq.mean = mean(tau_sq$tau_sq)

#phi
phi = rstan::extract(stan_fit0)['scale']
phi.mean = mean(phi$scale)

##################### BASIC KRIGING #######################
#1. save old data 
q.ref = q.extracted
mu.ref = X%*%beta.hat
X.ref = X
coords.ref = coords.tess[mask.data,]
mask.ref = mask.data

#2. generate new data by rerunning above code with only years 2019
q.new = q.extracted
mu.new = X%*%beta.hat
X.new = X
coords.new = coords.tess[mask.data,]
mask.new = mask.data


#2.5 quantify number of new locations we're predicting at -> 25% 
sum(!mask.new%in%mask.ref)/length(mask.new)

#generate the covariance matrix 
coords.augmented = rbind(coords.ref, coords.new)
DMat.augmented = as.matrix(dist(coords.augmented, method="euclidean", diag=T, upper=T))
dim(DMat.augmented)

#sum(DMat.augmented[1:166,1:166] != DMat) #sanity check that upper left block equals old distance matrix

#convert DMat to kernel distances using phi and sigma_sq
kernel = function(d){return (sigma_sq.mean*exp(-d/phi.mean))}
DMat.augmented = kernel(DMat.augmented)

#add the noise 
C = DMat.augmented + diag(rep(tau_sq.mean+0.01), dim(coords.augmented)[1])

#extract block matrices 
C11 = C[1:length(q.ref), 1:length(q.ref)]
C12 = C[1:length(q.ref), (length(q.ref)+1):dim(coords.augmented)[1]]
C21 = C[(length(q.ref)+1):dim(coords.augmented)[1], 1:length(q.ref)]
C22 = C[(length(q.ref)+1):dim(coords.augmented)[1], (length(q.ref)+1):dim(coords.augmented)[1]]
  
#prediction?
pred = mu.new + C21%*%solve(C11)%*%(log(q.ref)-mu.ref)
pred = ceiling(exp(pred)/3)

#compare to real one 
x11()
par(mfrow=c(1,2))
plot(log(pred))
plot(log(q.new))



x11()
par(mfrow=c(2,1))
plot(X%*%beta.hat, ylim=c(0,max(log(q.extracted))))
plot(log(q.extracted))








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


plot(density(sigma_sq$sigma_sq), main="", xlab="sigma")
plot(density(tau_sq$tau_sq), main="", xlab="tau")

