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
rf = 1000
dim = as.integer(c(max(lds$NORTHING)-min(lds$NORTHING), max(lds$EASTING)-min(lds$EASTING))/rf) + 1
dim 

ggplot(lds) + geom_sf() #plot data 

#plot intensity 
dens <- density(lds.ppp)
x11()
plot(dens, main = 'Intensity')


#can i make a crude little image for the locations?  -- BAD since we loose resolution (also missing vals??)
m = matrix(0, nrow=dim[1], ncol=dim[2])
for(i in 1:dim(lds)[1]){
  j = as.integer((lds$NORTHING[i] - min(lds$NORTHING))/rf)+1
  k = as.integer((lds$EASTING[i] - min(lds$EASTING))/rf)+1
  m[j,k]=m[j,k]+1
}
sum(m)
dim(lds)[1]

x11()
image(t(m))
