library(sf)
library(rgdal)
library(sp)
library(ggplot2)
#POINT PATTERN
library(maptools)
library(spatstat)

fname <- "../ENTLI.gdb.zip"
landslides = st_read(fname, "ENTLI_Crown")

trail = st_read(fname, "ENTLI_Trail")
lds_from_2000 = landslides[which(landslides$YEAR_1 > "2000"),]
lds = lds_from_2000


p.sp  <- as(lds, "Spatial")  # Create Spatial* object
lds.ppp <- as(p.sp, "ppp")      # Create ppp object

class(lds.ppp)


#creating `box` for our coords to fit inside based on easting and northing coords
#this box is a crude estimate for encapsulating hong kong - USED LATER MAYBE?
dim = as.integer(c(max(lds$NORTHING)-min(lds$NORTHING), max(lds$EASTING)-min(lds$EASTING))/1000)
dim 

ggplot(lds) + geom_sf() #plot data 

#plot intensity 
dens <- density(lds.ppp)
x11()
plot(dens, main = 'Intensity')



