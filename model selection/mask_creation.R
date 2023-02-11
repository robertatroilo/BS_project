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
u_NA = union(u_NA,is.na(lds_from_2000_with_NA$GULLY))
lds_from_2000=lds_from_2000_with_NA[-u_NA,]
lds = lds_from_2000
# remove less important kinds of cover (C and D) and wrong records of A and B
lds=lds[-cbind(which(lds$COVER=="C"), which(lds$COVER=="D"), which(lds$COVER==" A"), which(lds$COVER==" B")),]

dim(lds)
# 6556 points

n_quadr = 18

# Hong Kong shape

hk <- readShapePoly("sf/hh544xg3454.shp")
hkw = as.owin(hk)
hkw = owin(hkw$xrange, hkw$yrange) #the old window preserves boundary info but we want grid 
hkw

#choose a hkwdow size preserving the relative dimensions of hk 
r = diff(hkw$yrange)/diff(hkw$xrange)

spatstat.options(npixel=c(ceiling(n_quadr*r), n_quadr))  #n_quadr is ref for number of cols/x range 
spatstat.options()$npixel

grid <- makegrid(x=hk)
grid = spsample(hk, 2500,type="regular")
coords = grid@coords

hk_coords.ppp = ppp(coords[,1], coords[,2], hkw)
hk.q = quadratcount(hk_coords.ppp, ny=spatstat.options()$npixel[1], nx=spatstat.options()$npixel[2])

#define mask as points which are non-zero in the actual shapefile 
mask.hk = which(as.vector(t(hk.q), )!=0)


lat.long.df <- data.frame(lds$EASTING, lds$NORTHING) 
coordinates(lat.long.df) <-  ~lds.EASTING + lds.NORTHING
proj4string(lat.long.df)
proj4string(lat.long.df) <- CRS("+init=epsg:2326")

head(lat.long.df)
dist.location <- spTransform(lat.long.df, CRS("+init=epsg:4326"))

head(dist.location@coords)

features.drop = c("YEAR_1", "ENTLI_NO") #read data dict 
feature.mask = !names(lds)%in%features.drop
lds.ppp = ppp(dist.location@coords[,1], dist.location@coords[,2], hkw, marks=lds[names(lds)[feature.mask]]) #use hk window
lds.ppp = unique(lds.ppp)

ref = ppp(lds.ppp$x, lds.ppp$y, hkw) #using the hk window from mask_generation 
q <- quadratcount(lds.ppp, ny=spatstat.options()$npixel[1], nx=spatstat.options()$npixel[2])
dim(q)


#visualize 
x11()
par(mfrow=c(1,2))
plot(hk_coords.ppp, cex = 0.5, pch = "+", main="real map")
plot(hk.q, add=T)

plot(ref, cex = 0.5, pch = "+", main="landslide data")
plot(q, add=T)


mask.data = which(as.vector(t(q))!=0)

q.real = as.vector(t(q))
q.extracted = as.vector(t(q))[mask.data]

loss = sum(q.real) - sum(q.extracted)
loss/sum(q.real)

length(q.real)
length(q.extracted)


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
slope = as.vector(slope)[mask.data]
ele.diff = as.vector(ele.diff)[mask.data]
nocover.vote = as.vector(nocover.vote)[mask.data]
width = as.vector(width)[mask.data]
length = as.vector(length)[mask.data]
head.elevation = as.vector(head.elevation)[mask.data]
gully.vote = as.vector(gully.vote)[mask.data]

#normalize 
slope.norm = (slope - mean(slope))/sqrt(var(slope))
width.norm = (width - mean(width))/sqrt(var(width))
head.elevation.norm = (head.elevation - mean(head.elevation))/sqrt(var(head.elevation))
ele.diff.norm = (ele.diff - mean(ele.diff))/sqrt(var(ele.diff))
length.norm = (length - mean(length))/sqrt(var(length))
