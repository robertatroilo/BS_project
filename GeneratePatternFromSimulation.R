
#Matrix of Lambda: values are given from the coefficient predicted with stan
m_lam = matrix(exp(beta1 * covar4+beta2*covar5 + beta0 + mean_noise), n_quadr)
lambda = im(flipImage(m_lam))
plot(lambda)

#Total area
mu <- integral(lambda)

#create dataframe from image with coordinates
df <- as.data.frame(lambda)
#assign to each pixel its coordinates
for (i in c(1:n_quadr)){
  for (j in c(1:n_quadr)){
    df[,1][j+n_quadr*(i-1)] = x0[i]
    df[,2][j+n_quadr*(i-1)] = y0[j]
  }
}


npix <- nrow(df)

#each pixel has a probability to have points, weighted with its lambda value
lpix <- df$value

#number of total points drawn from a poisson
ni <- rpois(1, mu)

#sample indexes from df, each index is sampled a number of times proportional
#to the lambda value associated to its pixel
ii <- sample.int(npix, size=ni, replace=TRUE, prob=lpix/mu)

#Radius of the cells, to give "uncertainty" to the pattern
#points are in a random position within each cell
dx <- abs(x0[1]-x0[2])/2
dy <- abs(y0[1]-y0[2])/2

#generate coordinates of the points in the point pattern
#the second term gives a random position in the cell, around its centroid
xx <- df$x[ii] + runif(ni, -dx, dx)
yy <- df$y[ii] + runif(ni, -dy, dy)

#generate point pattern with rectangular window
dim <- rbind(c(113.83459, 114.44097), c(22.153194, 22.562094))
win <- owin(c(dim[1,1],dim[1,2]), c(dim[2,1],dim[2,2]))
pred_pp <- ppp(xx, yy, window=win, check=FALSE)
#change the window to irregular polygon (hong kong shape)
hk <- readShapePoly("sf/HKG_adm0.shp")#file found in "hk_boundaries.zip"
hkw = as.owin(hk)
Window(pred_pp) = hkw

#compare to real pattern
#Real pattern
X  <- as(lds, "Spatial")  
X <- as(X, "ppp")
marks(X) <- NULL
Window(X) = hkw

par(mfrow=c(2,2))
plot(density(pred_pp),main="Pred Density")
plot(density(X),main="Real Density")
plot(pred_pp,main="Pred Pattern",pch='+')
plot(X,main="Real Pattern",pch='+')




#Alternative Route (?)
hk <- readShapePoly("sf/HKG_adm0.shp")
hkw = as.owin(hk)
grid <- makegrid(x=hk)

grid = spsample(hk, 2500,type="regular")
grid@coords
#CREATE COOVARIATE IN THE CELL OF THE GRID ABOVE
grid@coords[1,2]
#compare if the distance between 
#abs(lds$long-grid@coords[1,1]) + abs(lds$lat-grid@coords[1,2]) < epsilon
