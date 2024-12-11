#install.packages("ggspatial")
#install.packages('ggsn')
sf::sf_use_s2(FALSE)

library(ggplot2)
library(cowplot)
#library(ggspatial)
#library(ggsn)
library(sf)


####Create object for basemap
#Bring ESRI world shapefile into R
setwd("C:/Users/Francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Thesis Draft")


usa <- sf::st_read("C:/Users/francis van oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/gadm36_USA_shp/gadm36_USA_0.shp")

world <-sf::st_read("C:/Users/Francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/gadm36_CAN_shp/gadm36_CAN_0.shp")
#conti.shelf<- sf::st_read("C:/Users/francis van oordt/OneDrive - McGill University/Documents/McGill/Field data/Global Margin shape/ContinentalMargins.shp")

#load map file

locsSIA<-read.csv("locations.csv", header = T, na.strings=c("", "NA"))
locsSIA<-read.csv("C:/Users/Francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/locations.csv", header = T, na.strings=c("", "NA"))

#lod raw df with 63 individuals for analysis
MMraw<-read.csv("C:/Users/Francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/multiDimNiche/raw_data/AASIA_raw_data.csv")
final.locations<-unique(MMraw$Location.x)


#subset for the locations from the clean 63 samples dataset
locsSIA<-locsSIA[locsSIA$Location %in% final.locations,]


#Use geom_sf to add basemap to your plot
ggplot2::theme_set(ggplot2::theme_light())
####Plot points on graph
plot <- ggplot2::ggplot() +
  ggplot2::geom_point(data = locsSIA, ggplot2::aes(x = Longitude, y = Latitude, color = Location), cex= 4) +
  ggplot2::xlab("Longitude")+
  ggplot2::ylab("Latitude")+
  ggplot2::geom_sf(data = world, ggplot2::aes()) + #add basemap
  ggplot2::geom_text(ggplot2::aes(x=LongLab, y=LatLab, label=Location), data=locsSIA, size=3, hjust=-0.15, col = "black")+
  ggplot2::coord_sf(xlim = c(-135, -123), ylim = c(48, 56)) +#Choose coordinate limits
  ggplot2::scale_color_manual(values = c("#999999", "#E69F00", 
                                         "#56B4E9", "#009E73", 
                                         "#F0E442", "#0072B2", 
                                         "#D55E00", "#CC79A7"))+
  ggplot2::theme(legend.position = "none")

plot

ggplot2::ggsave("final_locations_map_SIAA.png", plot, units = 'in', width = 7, height = 7)
