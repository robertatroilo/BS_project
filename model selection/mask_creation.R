#getwd()
#setwd("C:/Users/rober/Desktop/Bayesian Statistics/progetto")
#list.files()
#Load libraries
#----
library(sf)
library(rgdal)
library(sp)
library(spatstat)
library(rstan)
library(maptools)

#Load data
#----
fname <- "ENTLI.gdb.zip"

landslides = st_read(fname, "ENTLI_Crown")

# consider landslides since 2000
lds_from_2000_with_NA = landslides[which(landslides$YEAR_1 > "2000"),]


# data cleaning
u_NA = union(which(lds_from_2000_with_NA$M_WIDTH == 9999),which(lds_from_2000_with_NA$HEADELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$TAILELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$ELE_DIFF == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$SLOPE == 9999))
lds_from_2000=lds_from_2000_with_NA[-u_NA,]
lds = lds_from_2000
# remove less important kinds of cover (C and D) and wrong records of A and B
lds=lds[-cbind(which(lds$COVER=="C"), which(lds$COVER=="D"), which(lds$COVER==" A"), which(lds$COVER==" B")),]

dim(lds)
# 6556 points

# #number of quadrants
# n_quadr = 30
# 
# spatstat.options(npixel=c(n_quadr,n_quadr))#c(dim[1,2]-dim[1,1],dim[2,2]-dim[2,1]))
# dim <- rbind(c(801587.1, 859623.4), c(802346, 846054.5))
# win <- owin(c(dim[1,1],dim[1,2]), c(dim[2,1],dim[2,2]))
# 
# y0 <- seq(win$yrange[1], win$yrange[2],
#           length=spatstat.options()$npixel[2])
# x0 <- seq(win$xrange[1], win$xrange[2],
#           length=spatstat.options()$npixel[1])

# p.sp  <- as(lds, "Spatial")  # Create Spatial* object
# lds.ppp <- as(p.sp, "ppp")      # Create ppp object

lat.long.df <- data.frame(lds$EASTING, lds$NORTHING) 
coordinates(lat.long.df) <-  ~lds.EASTING + lds.NORTHING
proj4string(lat.long.df)
proj4string(lat.long.df) <- CRS("+init=epsg:2326")

head(lat.long.df)
dist.location <- spTransform(lat.long.df, CRS("+init=epsg:4326"))
head(dist.location)

head(dist.location@coords)

# Hong Kong shape

hk <- readShapePoly("sf/hh544xg3454.shp")
hkw = as.owin(hk)
hkw = owin(hkw$xrange, hkw$yrange) #the old window preserves boundary info but we want grid 
hkw

features.drop = c("SLIDE_TYPE", "YEAR_1", "ENTLI_NO") #read data dict 
feature.mask = !names(lds)%in%features.drop
lds.ppp = ppp(dist.location@coords[,1], dist.location@coords[,2], hkw, marks=lds[names(lds)[feature.mask]]) #use hk window
lds.ppp = unique(lds.ppp)


#choose a window size preserving the relative dimensions of hk 
win = hkw
n_quadr = 50
r = diff(win$yrange)/diff(win$xrange)

spatstat.options(npixel=c(ceiling(n_quadr*r), n_quadr))  #n_quadr is ref for number of cols/x range 
spatstat.options()$npixel

ref = ppp(lds.ppp$x, lds.ppp$y, win) 
q <- quadratcount(lds.ppp, ny=spatstat.options()$npixel[1], nx=spatstat.options()$npixel[2])
dim(q)

# plot the actual number of landslides in each quadrant
plot(ref)
plot(q, add=T)

# apply the HK mask on the data
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

# extract covariates from data and apply the mask to them before stan

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

# 
# define_covariate <- function(northing, easting, npixel = n_quadr, feature){
#   covar <- matrix(data = 0, nrow = npixel, ncol = npixel)
#   count <- matrix(data = 0, nrow = npixel, ncol = npixel)
#   n <- dim(lds)[1]
#   for (i in c(1:n)){
#     for(j in c(1:(npixel-1))){
#       for(z in c(1:(npixel-1))) {
#         if(northing[i]>y0[j] && northing[i]<y0[j+1]){
#           if(easting[i]>x0[z] && easting[i]<x0[z+1]){
#             # update the mean
#             covar[j,z] = (covar[j,z]*count[j,z] + feature[i])/(count[j,z]+1)
#             count[j,z] = count[j,z] + 1
#           }
#         }
#       }
#     }
#   }
#   return(covar)
# }
# 
# covar1 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$M_WIDTH)
# covar2 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$HEADELEV)
# covar3 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$TAILELEV)
# covar4 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$ELE_DIFF)
# covar5 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$SLOPE)

# # create dummy variable
# # A: Totally bare of vegetation
# # B: Partially bare of vegetation
# A <- which(lds$COVER=="A")   
# B <- which(lds$COVER=="B")
# nA <- length(A)
# nB <- length(B)
# cover <- rep(0,nA+nB)
# cover[A] = 1
# length(cover)

# covar6 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = cover)

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

#let's actually make the coordinates 
x0 <- seq(hkw$xrange[1], hkw$xrange[2], length=spatstat.options()$npixel[2])
y0 <- seq(hkw$yrange[1], hkw$yrange[2], length=spatstat.options()$npixel[1])

print(length(x0))
print(length(y0))

#the result of this is expanded row by row for the coords -- need consistency for the mask
coords.actual = expand.grid(x = x0, y = y0)

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
