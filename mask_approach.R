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
lds_from_2000_with_NA = landslides[which(landslides$YEAR_1 > "2000"),]


######## 2. DATA  CLEANING ##########
#should add some filters here for duplicates or any other stuff 
u_NA = union(which(lds_from_2000_with_NA$M_WIDTH == 9999),which(lds_from_2000_with_NA$HEADELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$TAILELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$ELE_DIFF == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$SLOPE == 9999))
lds_from_2000=lds_from_2000_with_NA[-u_NA,]
lds = lds_from_2000


######## 3. PRELIM SPATIAL OBJECTS  ##########
p.sp  <- as(lds, "Spatial")  # Create Spatial* object
lds.ppp <- as(p.sp, "ppp")      # Create ppp object

#plot(p.sp)
#plot(lds.ppp, which.marks = "ELE_DIFF")

#super cool (thanks chris)
hk <- readShapePoly("./hk_boundaries/hh544xg3454.shp")
hkw = as.owin(hk)
hkw = owin(hkw$xrange, hkw$yrange) #the old window preserves boundary info but we want grid 
hkw

win = as.owin(lds.ppp)


#choose a window size preserving the relative dimensions of hk 
n_quadr = 13
r = diff(win$yrange)/diff(win$xrange)

spatstat.options(npixel=c(n_quadr,ceiling(n_quadr*r)))
spatstat.options()$npixel

ref = ppp(lds$EASTING, lds$NORTHING, win) 
q <- quadratcount(lds.ppp, ny=spatstat.options()$npixel[2], nx=spatstat.options()$npixel[1])

#x11()
#plot(ref, cex = 0.5, pch = "+")
#plot(q, add=T) 


########## 4. generate mask from hk shapefile  #########
grid <- makegrid(x=hk)
grid = spsample(hk, 2500,type="regular")
coords = grid@coords

hk_coords.ppp = ppp(coords[,1], coords[,2], hkw)
hk.q = quadratcount(hk_coords.ppp, ny=spatstat.options()$npixel[2], nx=spatstat.options()$npixel[1])

#define mask as points which are non-zero in the actual shapefile 
mask = which(as.vector(t(hk.q))!=0)
mask

q.real = as.vector(t(q))
q.extracted = as.vector(t(q))[mask]

#we lose data since the shapes we're comparing are slightly different (see plot)
#at high resolutions especially this is bad since there's areas in the middle of hk
#which are 0 in the mask but contain a lot of points in our data 

loss = sum(q.real) - sum(q.extracted)
loss/sum(q.real)

length(q.real)
length(q.extracted)


x11()
par(mfrow=c(1,2))
plot(hk_coords.ppp, cex = 0.5, pch = "+", main="real map")
plot(hk.q, add=T)

plot(ref, cex = 0.5, pch = "+", main="landslide data")
plot(q, add=T)


#using this approach we will incur data losses. this is unavoidable since i can't
#reverse search the utm -> long lat coordinates as our utm convention in the data 
#is currently unknown.  for a refinement of 13x10 our losses are 2% 

#MAIN IDEA: define covariates as before, but apply mask to the covariate grid to 
#extract only those which are relevant 








