####################################################################
# import packages and upload data-set

packages<-c("ggmap","leaflet","ggplot2", "adehabitatHR","data.table","ggfortify","grid","move","moveVis","OpenStreetMap","pbapply","plotly","rgdal","sp","tidyverse","viridis","ggspatial")
sapply(packages, require, character.only=T)

library(readr)
American_mink_Neovison_vison_space_use_in_Illinois_data_from_Ahlers_et_al_2015_ <- read_csv("American mink (Neovison vison) space use in Illinois (data from Ahlers et al. 2015) .csv")
mink <- American_mink_Neovison_vison_space_use_in_Illinois_data_from_Ahlers_et_al_2015_
colnames(mink) <- c("event.id","visible", "timestamp","location.long","location.lat","location.error.numerical","sensor.type","individual.taxon","tag.local.identifier","individual.local.identifier","study.name","utm.easting","utm.northing","utm.zone")

Mink <- mink %>% slice(1:196)

######################################################################
# make a basic plot of the importata data that is plotted on an 
# interactive map that allows examination of data within each data point

qaqc_plot <- ggplot() + geom_point(data=Mink, aes(Mink$location.long,Mink$location.lat, color=Mink$individual.local.identifier))  +
  labs(x="Longitude", y="Latitude") +
  guides(color=guide_legend("individual.local.identifier"))

ggplotly(qaqc_plot)

#####################################################################
# Create sepearte csv files for each specified individual within your data set

lapply(split(Mink, Mink$individual.local.identifier), 
       function(x)write.csv(x, file = paste(x$individual.local.identifier[1],".csv"), row.names = FALSE))

#files <- list.files(path = ".", pattern = "[1-2]", full.names = TRUE)
files <- c("./1 .csv", "./2 .csv")

#####################################################################
# create a background map using leaflet

colors <- c("brown","yellow","darkgreen")

leaflet() %>% 
  addTiles() %>%
  addCircleMarkers(Mink$location.long,
                   Mink$location.lat,
                   weight = 1,
                   color = "blue",
                   fillColor = "blue",
                   fillOpacity = 0.7)

#######################################################################
# Create a function that can easily be applied to all of the files created
# to generate a minimum convex polygon 

mcp_raster <- function(filename){
  data <- read.csv(file = filename)
  x <- as.data.frame(data$utm.easting)
  y <- as.data.frame(data$utm.northing)
  xy <- c(x,y)
  data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS("+proj=utm +zone=16 +ellps=WGS84 +units=m +no_defs"))
  xy <- SpatialPoints(data.proj@coords)
  mcp.out <- mcp(xy, percent=100, unout="ha")
  mcp.points <- cbind((data.frame(xy)),data$individual.local.identifier)
  colnames(mcp.points) <- c("x","y", "identifier")
  mcp.poly <- fortify(mcp.out, region = "id")
  mcp.plot <- ggplot() + theme_bw() + theme(legend.position="none") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    geom_polygon(data=mcp.poly, aes(x=mcp.poly$long, y=mcp.poly$lat), alpha=0.8) +
    geom_point(data=mcp.points, aes(x=x, y=y)) + 
    labs(x="Easting (m)", y="Northing (m)", title=mcp.points$identifier) +
    theme(legend.position="none", plot.title = element_text(face = "bold"))
  mcp.plot
}

pblapply(files, mcp_raster)

################################################################################
# Create a function that can easily be applied to all of the files created
# to generate a kernal density map 

kde_raster <- function(filename){
  data <- read.csv(file = filename)
  x <- as.data.frame(data$utm.easting)
  y <- as.data.frame(data$utm.northing)
  xy <- c(x,y)
  data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS("+proj=utm +zone=16 +ellps=WGS84 +units=m +no_defs"))
  xy <- SpatialPoints(data.proj@coords)
  kde<-kernelUD(xy, h="href", kern="bivnorm", grid=100)
  ver95 <- getverticeshr(kde, 95)
  ver75 <- getverticeshr(kde, 75)
  ver50 <- getverticeshr(kde, 50)
  kde.points <- cbind((data.frame(data.proj@coords)),data$individual.local.identifier)
  colnames(kde.points) <- c("x","y","identifier")
  kde.poly95 <- fortify(ver95, region = "id")
  kde.poly75 <- fortify(ver75, region = "id")
  kde.poly50 <- fortify(ver50, region = "id")
  # use get_stamenmap to create an object prior to this function without using ggmap. 
  # then follow the instructions on the github page to convert the coordinated from lat/long 
  # to UTM. You will need to figure out what class the "raster" is, it may just be a raster 
  # class though so try that first.
  # or. . . you can also always just use the orignal code and run it in the computer lab so 
  # that open street map will work.
  kde.plot <- ggmap(get_stamenmap(bbox = c(left = -95.80204, bottom = 29.38048, right = -94.92313, top = 30.14344), zoom = 10, maptype = c("terrain"))) +
    theme_bw() + theme(legend.position="none") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    geom_polygon(data=kde.poly95, aes(x=kde.poly95$long, y=kde.poly95$lat, group=group, fill="green"), alpha = 0.8) +
    geom_polygon(data=kde.poly75, aes(x=kde.poly75$long, y=kde.poly75$lat, group=group, fill="blue"), alpha = 0.8) +
    geom_polygon(data=kde.poly50, aes(x=kde.poly50$long, y=kde.poly50$lat, group=group, fill="red"), alpha = 0.8) +
    geom_point(data=kde.points, aes(x=x, y=y)) +
    labs(x="Easting (m)", y="Northing (m)", title=kde.points$identifier) +
    theme(legend.position="none", plot.title = element_text(face = "bold"))
  kde.plot
}

pblapply(files, kde_raster)

##################################################################################
# Creation of a Brownian Bridge movement map (Wont work with this dataset. requires 
# very, very large number of data point. example in calss used around 40,000 for 2 
# individuals)

m1 <- read.csv("./1 .csv")

date <- as.POSIXct(strptime(as.character(m1$timestamp),"%Y-%m-%d %H:%M:%S", tz="America/Chicago"))

m1$date <- date

m1.reloc <- cbind.data.frame(m1$utm.easting, m1$utm.northing,
                                as.vector(m1$individual.local.identifier),
                                as.POSIXct(date))

colnames(m1.reloc) <- c("x","y","id","date")

trajectory <- as.ltraj(m1.reloc, date=date, id="m1")

sig1 <- liker(trajectory, sig2 = 58, rangesig1 = c(0, 5), plotit = FALSE)

m1.traj <- kernelbb(trajectory, sig1 = .7908, sig2 = 58, grid = 100)

bb_ver <- getverticeshr(opha.traj, 95)

bb_poly <- fortify(bb_ver, region = "id", 
                   proj4string = CRS("+proj=utm +zone=16+
                                     ellps=WGS84 +units=m +no_defs"))

colnames(bb_poly) <- c("x","y","order","hole","piece","id","group")

bb_image <- crop(m1.traj, bb_ver, 
                 proj4string = CRS("+proj=utm +zone=16 +
                                   ellps=WGS84 +units=m +no_defs"))

bb_units <- grid.text(paste(round(bb_ver$area,2)," ha"), x=0.85,  y=0.95,
                      gp=gpar(fontface=4, col="white", cex=0.9), draw = FALSE)

bb.plot <- ggplot() + theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_tile(data=bb_image, 
            aes(x=bb_image@coords[,1], y=bb_image@coords[,2],
                fill = bb_image@data$ud)) +
  geom_polygon(data=bb_poly, aes(x=x, y=y, group = group), color = "black", fill = NA) +
  scale_fill_viridis_c(option = "inferno") + annotation_custom(bb_units) +
  labs(x="Easting (m)", y="Northing (m)", title="OPHA1") +
  theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5))

bb.plot

################################################################################
# Creation of a movement map

one <- read_csv("./1 .csv")
m1.move <- move(x=one$location.long, 
                  y=one$location.lat, 
                  time=as.POSIXct(one$timestamp, 
                                  format="%Y-%m-%d %H:%M:%S", tz="America/Chicago"), 
                  proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                  data=one, animal=one$individual.local.identifier, 
                  sensor=one$sensor.type)

movement <- align_move(m1.move, res = "max", digit = 0, unit = "secs")

frames <- frames_spatial(movement, path_colours = "red",
                         map_service = "osm",
                         alpha = 0.5) %>% 
  add_labels(x = "Longitude", y = "Latitude") %>%
  add_northarrow() %>% 
  add_scalebar() %>% 
  add_timestamps(movement, type = "label") %>% 
  add_progress()

animate_frames(frames, fps = 5, overwrite = TRUE,
               out_file = "./moveVis-5fps.gif")

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#ggmap
#data type: stamen
# getstamenmap
# ggmap instead of ggplot

















