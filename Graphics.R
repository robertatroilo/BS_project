library(sf)
library(rgdal)
library(sp)
library(spatstat)
library(rstan)
library(maptools)
library(raster)
#Load data
#----
fname <- "ENTLI.gdb.zip"
fname
st_layers(fname)

landslides = st_read(fname, "ENTLI_Crown")

# consider landslides since 2000
lds_from_2000_with_NA = landslides[which(landslides$YEAR_1 > "2000"),]
#lds_from_2000_with_NA = landslides
# 6629 obs

# data cleaning
u_NA = union(which(lds_from_2000_with_NA$M_WIDTH == 9999),which(lds_from_2000_with_NA$HEADELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$TAILELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$ELE_DIFF == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$SLOPE == 9999))
lds_from_2000=lds_from_2000_with_NA[-u_NA,]
lds = lds_from_2000
#all_lds = lds_from_2000
dim(lds)
# 6588 points

# plot M_WIDTH
plot(lds[2])
# plot HEADELEV
plot(lds[7])
# plot TAILELEV
plot(lds[8])
# plot Elev_diff
plot(lds[9])
plot(density(lds$ELE_DIFF))
# plot SLOPE
plot(lds[4])
#----


#number of quadrants
n_quadr = 50

spatstat.options(npixel=c(n_quadr,n_quadr))#c(dim[1,2]-dim[1,1],dim[2,2]-dim[2,1]))
dim <- rbind(c(801050.7, 859623.4), c(801788.9, 846548.8)) #ALL LDS, LARGER WINDOW
#dim <- rbind(c(801587.1, 859623.4), c(802346, 846054.5))
win <- owin(c(dim[1,1],dim[1,2]), c(dim[2,1],dim[2,2]))

y0 <- seq(win$yrange[1], win$yrange[2],
          length=spatstat.options()$npixel[2])
x0 <- seq(win$xrange[1], win$xrange[2],
          length=spatstat.options()$npixel[1])


define_covariate <- function(northing, easting, npixel = n_quadr, feature){
  covar <- matrix(data = 0, nrow = npixel, ncol = npixel)
  count <- matrix(data = 0, nrow = npixel, ncol = npixel)
  n <- dim(lds)[1]#
  n <- dim(all_lds)[1] # TO CONSIDER ALL THE LANDSLIDE: MORE POINTS, more zeros
  for (i in c(1:n)){
    for(j in c(1:(npixel-1))){
      if(northing[i]>y0[j] && northing[i]<y0[j+1]){
        for(z in c(1:(npixel-1))) {
          if(easting[i]>x0[z] && easting[i]<x0[z+1]){
            # update the mean
            covar[j,z] = (covar[j,z]*count[j,z] + feature[i])/(count[j,z]+1)
            count[j,z] = count[j,z] + 1   
           } }}}}
  return(covar)
}

covar1 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$M_WIDTH)
covar2 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$HEADELEV)
covar3 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$TAILELEV)
covar4 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$ELE_DIFF)
covar5 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$SLOPE)

covar4_all_lds =  define_covariate(northing = all_lds$NORTHING, easting = all_lds$EASTING, feature = all_lds$ELE_DIFF)

plot(density(covar4))
plot(im(covar4_all_lds))

p.sp  <- as(lds, "Spatial")  # Create Spatial* object
lds.ppp <- as(p.sp, "ppp")      # Create ppp object


#Plot tess
b <- quantile(im(covar5), probs = (0:4)/4)
b = b + seq_along(b) * .Machine$double.eps
slope_break <- cut(im(covar5), breaks = b, labels = 1:4)
v <- tess(image = slope_break)

plot(v)


b <- quantile(lds_from_2000$SLOPE, probs = (0:4)/4)
slope_break <- cut(lds_from_2000$SLOPE, breaks = b, labels = 1:4)
v <- tess(image = slope_break)
plot(v)


##########################################################
#PLOT LANDSLIDES OVER HONG KONG BOUNDARIES
library(ggplot2)
shape <- read_sf(dsn = "hk_boundaries", layer = "hh544xg3454")
plot(shape$geometry)
ggplot() + geom_sf(data=shape$geometry) + geom_sf(data=lds_from_2000)

x11()
ggplot() + geom_sf(data=shape$geometry) + 
  geom_sf(data=lds_from_2000,aes(colour = lds_from_2000$SLIDE_TYPE),pch='+')


x11()
ggplot() + geom_sf(data=shape$geometry) + 
  geom_sf(data=landslides[which(landslides$YEAR_1 < "2000"),],aes(colour = landslides[which(landslides$YEAR_1 < "2000"),]$SLIDE_TYPE),pch='+')


#############################################
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

#PPP object with map boundaries
X  <- as(lds_map, "Spatial")  
#X <- as(p.sp, "ppp") [hkw]
X <- as(X, "ppp")
marks(X) <- NULL
Window(X) = hkw
plot(X, pch = '.')

