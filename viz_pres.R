library(sf)
library(rgdal)
fname <- "ENTLI.gdb.zip"
fname

# new <- st_read(fname)

st_layers(fname)

library(sp)
methods(st_as_sf)
methods(st_as_sfc)

landslides = st_read(fname, "ENTLI_Crown")


library(ggplot2)
library("RColorBrewer")

#devtools::install_github("paleolimbot/ggspatial")
library(ggspatial)
library(raster)


lds=landslides[-which(landslides$HEADELEV == 9999),]

HK<- getData("GADM",country='HKG', level = 1)

HK_test = st_as_sf(HK)

p <- ggplot() + scale_color_gradient(low="blue", high="red") + 
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
p
ggsave('map_headelev.png', p, bg='transparent', scale = 3)

dev.off()


