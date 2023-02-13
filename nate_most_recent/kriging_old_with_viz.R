source("./mask_generation.R")
library(raster)
library(tidyverse)

fname <- "../ENTLI.gdb.zip"

######## 1. LOADING IN THE DATA  ##########
# consider landslides since 2000
landslides = st_read(fname, "ENTLI_Crown")


YEARS = c("2016", "2017", "2018")
lds_from_2000_with_NA = landslides[which(landslides$YEAR_1 %in% YEARS),]

# HK outline data
HK_map <- getData("GADM",country='HKG', level = 1)
HK_outline = st_as_sf(HK_map)$geometry




######## 2. DATA  CLEANING ##########
u_NA = union(which(lds_from_2000_with_NA$M_WIDTH == 9999),which(lds_from_2000_with_NA$HEADELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$TAILELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$ELE_DIFF == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$SLOPE == 9999))
u_NA = union(u_NA,is.na(lds_from_2000_with_NA$GULLY))

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

head(dist.location@coords)


############ create ppp object using extracted coords & our data ################
features.drop = c("YEAR_1", "ENTLI_NO") #read data dict 
feature.mask = !names(lds)%in%features.drop
lds.ppp = ppp(dist.location@coords[,1], dist.location@coords[,2], hkw, marks=lds[names(lds)[feature.mask]]) #use hk window
lds.ppp = unique(lds.ppp)


ref = ppp(lds.ppp$x, lds.ppp$y, hkw) #using the hk window from mask_generation 
q <- quadratcount(lds.ppp, ny=spatstat.options()$npixel[1], nx=spatstat.options()$npixel[2])
dim(q)

# initial plot where we can visualize the ppp, the quadrant, and the outline of HK
plot(ref)
plot(q, add=T)
plot(HK_outline, add = T)

########## sample reduction; we can either use hk mask OR only keep cells which are in the sample \
#before we kept all the samples which were inside HK but instead we should really only consider our samples
#idk if this is a little biased or what will happen 

mask.data = which(as.vector(t(q))!=0)



# intensity + pts
plot(intensity(q, image=TRUE), main=NULL, las=1)  # Plot density raster
points(lds.ppp)
plot(HK_outline, add = T)


# intensity + mask + outline
plot(intensity(q, image=TRUE), main=NULL, las=1, zlim=c(0.1,max(intensity(q))))  # Plot density raster
plot(HK_outline, add = T)

# same but messy
plot(intensity(q, image=TRUE), main=NULL, las=1, zlim=c(0.1,max(intensity(q))))  # Plot density raster
plot(q, add = T)
plot(HK_outline, add = T)
points(lds.ppp)




q.real = as.vector(t(q))
q.extracted = as.vector(t(q))[mask.data]


# TEST FOR THE RECONSTRUCTION OF THE INTENSITY MATRIX
reconstruct_q = rep(0, 13*18)
i = 1
for(index in which(q.real > 0)){
  reconstruct_q[index] = q.extracted[i]
  i=i+1
}
reconstruct_q == t(q) # Everything seems ok
reconstruct_q = matrix(reconstruct_q, nrow = 13, byrow = T)
dim(reconstruct_q)
dim(q)

# WE GET THE SAME IMAGE
image(reconstruct_q, useRaster=TRUE, axes=FALSE)
image(q, useRaster=TRUE, axes=FALSE) 
plot(HK_outline, add = T) # WE CANNOT ADD THE OUTLINE ONE TOP OF THIS



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

#now derive the CAR matrix we employ (from notes on geostats) -> this isnt true adjacency but its fine 
dmat.quantile = quantile(as.vector(DMat), 0.07) # 0.05 was not working, had to increase to 0.7
A = ifelse(DMat<=dmat.quantile, 1, 0)
diag(A) = 0

Dw = diag(rowSums(A))


rho = 0.9
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
  
  // W standard deviation parameter -- CAR
  vector[N] e;
  real<lower=0> sigma_sqCAR;
  
  //random noise independent
  vector[N] w;
  real<lower=0> sigma_sqI;
}

transformed parameters {
    vector[N] mu;
    for(i in 1:N) {
      mu[i] = exp(row(X, i) * beta + e[i] + w[i]); //link 
    }
}

model{
  for (s in 1:N) {
        y[s] ~ poisson(mu[s]);  
  } 
  
  e ~ multi_normal(rep_vector(0,N), sigma_sqCAR*W); 
  w ~ normal(0, sigma_sqI); 
    
  //prior for the noises 
  sigma_sqCAR ~ inv_gamma(0.01,0.01); //kinda from lecture notes
  sigma_sqI ~ inv_gamma(1,1);
  
  //priors for the coefficients
  for (j in 1:p) {
        beta[j] ~ normal(0.0, sigma_beta);
        sigma_beta[j] ~ inv_gamma(1, 1);
    }
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
X = cbind(rep(1,length(q.extracted)), slope.norm, width.norm, head.elevation.norm,
          ele.diff.norm, length.norm, vegetation.vote, erosion.vote)

#non-normalized covariates
#X = cbind(rep(1,length(q.extracted)), slope, width, head.elevation, ele.diff, length, vegetation.vote, erosion.vote)

N = dim(X)[1]
p = dim(X)[2]

stan_data = list(N = N, p = p, X = X, y = q.extracted, W=W)

stan_fit0 <- stan(model_code = ppm_stan,
                  data = stan_data,
                  chains = 1, warmup = 1000, iter = 6000,
                  control = list(adapt_delta = 0.999, max_treedepth=13))




############# plotting the chains ###############3
beta.chain = rstan::extract(stan_fit0)['beta']
beta.hat = colMeans(beta.chain$beta)

#sigma = rstan::extract(stan_fit0)['sigma_sq']
#sigma.mean = mean(sigma$sigma_sq)


sigma = rstan::extract(stan_fit0)['sigma_sqI']
#sigma = rstan::extract(stan_fit0)['sigma_sqCAR']
sigma.mean = mean(sigma$sigma_sqI)


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


#plot(density(sigma$sigma), main="", xlab="sigma")
plot(density(sigma$sigma_sqI), main="", xlab="sigma")

#some predictions 
pred = ceiling(exp(X%*%beta.hat + sigma.mean)) #our predictoins
par(mfrow = c(1,2))
plot(log(pred))
plot(log(q.extracted))







# We create a vector that we only populate on the indexes different from 0
pred_matrix = rep(0, 13*18)

i = 1
for(index in which(q.real > 0)){
  pred_matrix[index] = pred[i]
  i=i+1
}

# Turn vector into a matrix
pred_matrix = matrix(pred_matrix, nrow = 13, byrow = T)

# Check that the dimensions match
dim(pred_matrix)
dim(q)

# plot everything
image(pred_matrix, useRaster=TRUE, axes=FALSE)
image(q, useRaster=TRUE, axes=FALSE)

# Can't add the outline 
# plot(HK_outline)

# plot(intensity(q, image=TRUE), main=NULL, las=1, zlim=c(0.1,max(intensity(q))))  # Plot density raster
# plot(HK_outline, add = T)



colnames(pred_matrix) <- colnames(q)
rownames(pred_matrix) <- rownames(q)


image(pred_matrix, useRaster=TRUE, axes=TRUE)

pred_matrix


# WE RECREATE THE PPP


prediction_vector = data.frame(x = NULL, y = NULL)

# go through all line, all column, all observation in each cell
for(y in 1:length(colnames(pred_matrix))){
  for(x in 1:length(rownames(pred_matrix))){
    if(pred_matrix[x, y])
      for(i in 1:pred_matrix[x, y]){
        pt = data.frame(x = rownames(pred_matrix)[x], y = colnames(pred_matrix)[y])
        prediction_vector = rbind(prediction_vector, pt)
      }
        
  }
}

prediction_vector

regex_ = '\\)'

test = prediction_vector
test$x = sapply(prediction_vector$x,  gsub, pattern = regex_, replacement = "\\]" )
test$y = sapply(prediction_vector$y,  gsub, pattern = regex_, replacement = "\\]" )
test

test$x = sapply(test$x,  gsub, pattern = "\\]", replacement = "" )
test$y = sapply(test$y,  gsub, pattern = "\\]", replacement = "" )
test$x = sapply(test$x,  gsub, pattern = "\\[", replacement = "" )
test$y = sapply(test$y,  gsub, pattern = "\\[", replacement = "" )
test



# put everything to the centroids of the tiles

# X's
upper_x = sapply(test$x,  str_split_i, pattern = ",", simplify = T, i = 1) 
upper_x = sapply(upper_x,  as.double) 
upper_x

lower_x = sapply(test$x,  str_split_i, pattern = ",", simplify = T, i = 2) 
lower_x = sapply(lower_x,  as.double) 
lower_x

centroids_x = (upper_x + lower_x)/2
centroids_x

# Y's
upper_y = sapply(test$y,  str_split_i, pattern = ",", simplify = T, i = 1) 
upper_y = sapply(upper_y,  as.double) 
upper_y

lower_y = sapply(test$y,  str_split_i, pattern = ",", simplify = T, i = 2) 
lower_y = sapply(lower_y,  as.double) 
lower_y

centroids_y = (upper_y + lower_y)/2
centroids_y




# I don't know why i had to invert x & y
predicted_coord = data.frame(x = centroids_y, y = centroids_x)

# every predicted landslide has value 1
predicted_dummy_ones = data.frame(rep(1, dim(predicted_coord)[1])) 

predicted_spatial_points_df = SpatialPointsDataFrame(coords = predicted_coord, data = predicted_dummy_ones)
predicted.ppp = as.ppp(predicted_spatial_points_df)

# Check if we have the same window 
Window(ref)
Window(predicted.ppp)

plot(predicted.ppp)

# Now we can vizualize quadrants 
predicted.q <- quadratcount(predicted.ppp, ny=spatstat.options()$npixel[1], nx=spatstat.options()$npixel[2])
plot(predicted.q)
plot(HK_outline, add = T)


# predictions but messy
plot(intensity(predicted.q, image=TRUE), main=NULL, las=1, zlim=c(0.1,max(intensity(q))))  # Plot density raster
plot(predicted.q, add = T)
plot(HK_outline, add = T)

# original plos but messy
plot(intensity(q, image=TRUE), main=NULL, las=1, zlim=c(0.1,max(intensity(q))))  # Plot density raster
plot(q, add = T)
plot(HK_outline, add = T)

# see the difference in scaling
plot(predicted.q)
plot(q, add = T)



