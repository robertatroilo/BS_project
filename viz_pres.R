library(sf)
library(rgdal)
library(sp)
library(ggplot2)
library("RColorBrewer")
library(spatstat)

#devtools::install_github("paleolimbot/ggspatial")
library(ggspatial)
library(raster)
library(ggpubr)



fname <- "ENTLI.gdb.zip"
fname
st_layers(fname)
methods(st_as_sf)
methods(st_as_sfc)
landslides = st_read(fname, "ENTLI_Crown")

# getting rid of the "NA"
lds=landslides[-which(landslides$HEADELEV == 9999),]

# preparing the map
HK<- getData("GADM",country='HKG', level = 1)
HK_test = st_as_sf(HK)


# ploting the head elevation over the map of HK
p1 <- ggplot() + scale_color_gradient(low="blue", high="red") + 
  geom_sf(data = HK_test, linewidth = 2) + 
  geom_sf(data = lds, alpha = 0.7, size = 0.5, aes(color = HEADELEV)) +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    #panel.grid.major = element_blank(), #remove major gridlines
    #panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    #legend.box.background = element_blank(), #transparent legend panel
    axis.text.x = element_text(color="#ffffff",
                               face="bold",
                               angle=0,
                               size = 30),
    axis.text.y = element_text(color="#ffffff",
                               face="bold",
                               angle=0,
                               size = 30),
    legend.title = element_text(color="#ffffff",
                                face="bold",
                                angle=0,
                                size = 30),
    legend.text = element_text(color="#ffffff",
                                face="bold",
                                angle=0,
                                size = 30),
  #  #legend.box.background = element_rect(color = NA)
  )
p1
ggsave('map_headelev.png', p1, bg='transparent', scale = 3)

dev.off()
###################################################
###################################################
###################################################


# consider landslides since 2000
lds_from_2000_with_NA = landslides[which(landslides$YEAR_1 > "2000"),]
# 6629 obs

# data cleaning
u_NA = union(which(lds_from_2000_with_NA$M_WIDTH == 9999),which(lds_from_2000_with_NA$HEADELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$TAILELEV == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$ELE_DIFF == 9999))
u_NA = union(u_NA,which(lds_from_2000_with_NA$SLOPE == 9999))
lds_from_2000=lds_from_2000_with_NA[-u_NA,]
lds = lds_from_2000





# ploting all the covariates KDE



p_kde_1 <- ggplot(data = lds) +  geom_density(aes(x = M_WIDTH), color = "blue", fill = "blue", alpha = 0.1)
p_kde_2 <- ggplot(data = lds) +  geom_density(aes(x = HEADELEV), color = "blue", fill = "blue", alpha = 0.1)
p_kde_3 <- ggplot(data = lds) +  geom_density(aes(x = TAILELEV), color = "blue", fill = "blue", alpha = 0.1)
p_kde_4 <- ggplot(data = lds) +  geom_density(aes(x = ELE_DIFF), color = "blue", fill = "blue", alpha = 0.1)
p_kde_5 <- ggplot(data = lds) +  geom_density(aes(x = SLOPE), color = "blue", fill = "blue", alpha = 0.1)

ggarrange(p_kde_1, p_kde_2, p_kde_3, p_kde_4, p_kde_5)



###################################################
###################################################
###################################################
# other covar plots
# From Samu's code
n_quadr = 30

spatstat.options(npixel=c(n_quadr,n_quadr))#c(dim[1,2]-dim[1,1],dim[2,2]-dim[2,1]))
dim <- rbind(c(801587.1, 859623.4), c(802346, 846054.5))
win <- owin(c(dim[1,1],dim[1,2]), c(dim[2,1],dim[2,2]))

y0 <- seq(win$yrange[1], win$yrange[2],
          length=spatstat.options()$npixel[2])
x0 <- seq(win$xrange[1], win$xrange[2],
          length=spatstat.options()$npixel[1])


define_covariate <- function(northing, easting, npixel = n_quadr, feature){
  covar <- matrix(data = 0, nrow = npixel, ncol = npixel)
  count <- matrix(data = 0, nrow = npixel, ncol = npixel)
  n <- dim(lds)[1]
  for (i in c(1:n)){
    for(j in c(1:(npixel-1))){
      for(z in c(1:(npixel-1))) {
        if(northing[i]>y0[j] && northing[i]<y0[j+1]){
          if(easting[i]>x0[z] && easting[i]<x0[z+1]){
            # update the mean
            covar[j,z] = (covar[j,z]*count[j,z] + feature[i])/(count[j,z]+1)
            count[j,z] = count[j,z] + 1
          }
        }
      }
    }
  }
  return(covar)
}

covar1 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$M_WIDTH)
covar2 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$HEADELEV)
covar3 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$TAILELEV)
covar4 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$ELE_DIFF)
covar5 <- define_covariate(northing = lds$NORTHING, easting = lds$EASTING, feature = lds$SLOPE)

p.sp  <- as(lds, "Spatial")  # Create Spatial* object
lds.ppp <- as(p.sp, "ppp")      # Create ppp object


x11()
par(mfrow=c(1,2))
plot(im(covar4), main = 'ELE_DIFF')
plot(im(covar5), main = 'SLOPE')



library(tidyr)
library(tibble)
library(hrbrthemes)
library(dplyr)

min_easting = min(landslides$EASTING)
max_easting = max(landslides$EASTING)

min_northing = min(landslides$NORTHING)
max_northing = max(landslides$NORTHING)



min_easting = 113.83
max_easting = 114.45

min_northing = 22.19
max_northing = 22.55

eastings_seq = seq(from = min_easting, to = max_easting, length.out = n_quadr)
northings_seq = seq(from = min_northing, to = max_northing, length.out = n_quadr)

test = covar4
test = as_tibble(test)
test = rowid_to_column(test, var="X")
test = gather(test, key="Y", value="Z", -1)
test = mutate(test, Y=as.numeric(gsub("V","",Y)))
test$X = rep(eastings_seq, n_quadr)
test$Y = rep(northings_seq, each = n_quadr)


p <- ggplot() + 
  geom_tile(data = test, aes(X, Y, fill= Z), alpha = 0.9) + 
  geom_sf(data = HK_test, linewidth = 2, fill = "transparent") +
  scale_fill_gradient(low="purple", high="yellow")
p  







p1 <- ggplot() + scale_color_gradient(low="blue", high="red") + 
  geom_sf(data = HK_test, linewidth = 2) + 
  geom_sf(data = covar4, alpha = 0.7, size = 0.5, aes(color = HEADELEV)) +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    #panel.grid.major = element_blank(), #remove major gridlines
    #panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    #legend.box.background = element_blank(), #transparent legend panel
    axis.text.x = element_text(color="#ffffff",
                               face="bold",
                               angle=0,
                               size = 30),
    axis.text.y = element_text(color="#ffffff",
                               face="bold",
                               angle=0,
                               size = 30),
    legend.title = element_text(color="#ffffff",
                                face="bold",
                                angle=0,
                                size = 30),
    legend.text = element_text(color="#ffffff",
                               face="bold",
                               angle=0,
                               size = 30),
    #  #legend.box.background = element_rect(color = NA)
  )
p1



###################################################
###################################################
###################################################




