#copying roby's idea to just include this in subsequent files 

library(sf)
library(rgdal)
library(sp)
library(spatstat)
library(rstan)
library(maptools)

#HYPERPARAMS 
n_quadr = 18

#read from shape file 
hk <- readShapePoly(".././hk_boundaries/hh544xg3454.shp")
hkw = as.owin(hk)
hkw = owin(hkw$xrange, hkw$yrange) #the old hkwdow preserves boundary info but we want grid 
hkw

#choose a hkwdow size preserving the relative dimensions of hk 
r = diff(hkw$yrange)/diff(hkw$xrange)

spatstat.options(npixel=c(ceiling(n_quadr*r), n_quadr))  #n_quadr is ref for number of cols/x range 
spatstat.options()$npixel

#chris again
grid <- makegrid(x=hk)
grid = spsample(hk, 2500,type="regular")
coords = grid@coords

hk_coords.ppp = ppp(coords[,1], coords[,2], hkw)
hk.q = quadratcount(hk_coords.ppp, ny=spatstat.options()$npixel[1], nx=spatstat.options()$npixel[2])

#plot(hk_coords.ppp, cex = 0.5, pch = "+", main="real map")
#plot(hk.q, add=T)

#define mask as points which are non-zero in the actual shapefile 
mask.hk = which(as.vector(t(hk.q), )!=0)

