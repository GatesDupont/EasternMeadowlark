# Gates Dupont                 #
# Eastern Meadowlark           #
# MassWildlife & Mass Audubon  #
# March 2019                   #
# # # # # # # # # # # # # # # #

library(FedData)
library(cdlTools)
library(dismo)
library(dplyr)
library(tidyr)
library(raster)
library(velox)
library(rgeos)
library(latticeExtra)
library(randomForest)
library(caTools)
library(leaflet)
library(maps)
library(viridis)
library(sp)

# Territory radius = sqrt(24281.1/pi)

#------------------------------------0. Extent------------------------------------

#----*State Contour----
states.full = c("Massachusetts")
us = getData('GADM', country = 'US', level = 1)
st.contour = us[us$NAME_1 %in% states.full,]
st.contour = spTransform(st.contour, CRS("+init=epsg:4326"))
st.contour = aggregate(st.contour) # remove boundaries between states

#----*Extent points---
ur = c(bbox(st.contour)[2,2], bbox(st.contour)[1,2])
ll = c(bbox(st.contour)[2,1], bbox(st.contour)[1,1])
ul = c(ur[1],ll[2])
lr = c(ll[1],ur[2])

extent.long = c(ur[2],ll[2],ul[2],lr[2])
extent.lat = c(ur[1],ll[1],ul[1],lr[1])
extent.df = data.frame(cbind(extent.long, extent.lat))

#----*Spatial points object with embedded extent----
coordinates(extent.df) = ~extent.long+extent.lat
study.extent.corners = SpatialPoints(extent.df, CRS("+init=epsg:4326"))

#------------------------------------1. GET Ocurrence------------------------------------

#----*Pulling GBIF data----
gbif.pull = gbif("Sturnella", "magna", geo = TRUE, removeZeros = 0, ext = extent(st.contour)) 

#----*Cleaning up data----
gbif = gbif.pull %>%
  filter(adm1 == "Massachusetts") %>%
  mutate(date = substr(eventDate,1,10)) %>%
  filter(year > 2015) %>%
  filter(month > 4 & month < 8) %>%
  select(lat=69, long=74, month=76, year=119) %>%
  filter(complete.cases(.)) %>%
  select(lat=1, long=2) %>%
  distinct(.keep_all = T)

#----*Generating sampling grid----
r = raster(st.contour)
res(r) = 0.075 * 2.5
r = extend(r, extent(r)+1)

#----*Grid Sampling----
coordinates(gbif) = ~long+lat
pres.samp = gridSample(gbif, r, n=1)
pres = as.data.frame(pres.samp)

#----*Plotting points----
plot(st.contour)
points(lat~long, data=pres, cex=0.4, col="blue", pch=20)

#----*Plotting on map to remove suspicious locations----
leaflet(pres) %>% addProviderTiles("Esri.WorldImagery", group = "ESRI") %>%
  addMarkers()

#----*Pseudo-absence points----
set.seed(1)
bg = as.data.frame(spsample(st.contour,n=length(pres$lat)*10,"random"))
colnames(bg) = c("long", "lat")
coordinates(bg) = ~long+lat
bg = data.frame(bg@coords)

#----*Combining presence/absence----
species = data.frame(rbind(pres, bg))
species$pa = c(rep(1, length(pres$lat)), rep(0, length(bg$lat)))
species.spatial = SpatialPointsDataFrame(cbind(species$lon, species$lat), data=species, proj4string = CRS("+init=epsg:4326"))

#------------------------------------2. GET CropScape------------------------------------

#----*Pulling cropscape data----
CropScape = getCDL("MA", 2017)
crops = crop(crops.raw, extent(spTransform(study.extent.corners, crs(CropScape$MA2017)))*1.1)

#------------------------------------3. GET NLCD------------------------------------

#----*Puilling nlcd data----
label = "MA"
nlcd = get_nlcd(template=crops, label=label, year = 2011, force.redo = T) # auto-crops

#------------------------------------4. GET Elevation------------------------------------

#----*Pulling elevation data----
lat.increment = (ur[1]-ll[1])/2
long.increment = (ur[2]-ll[2])/3
e.lats = c(ul[1], ul[1]-lat.increment, ul[1]-(2*lat.increment))
e.longs = c(ul[2], ul[2]+long.increment, ul[2]+(2*long.increment), ul[2]+(3*long.increment))

srtm1 = getData("SRTM", lat=e.lats[1], lon=e.longs[1])
srtm2 = getData("SRTM", lat=e.lats[2], lon=e.longs[4])

#----*Stitching elevation tiles----
elevation = mosaic(srtm1, srtm2, fun=max)
elevation = crop(elevation, extent(spTransform(study.extent.corners, crs(elevation)))*1.1)

#------------------------------------5. GET WorldClim------------------------------------
climate = getData('worldclim', var='bio', res=0.5, lat=ll[1], lon=ll[2])
climate = crop(climate, extent(spTransform(study.extent.corners, crs(climate)))*1.1)

#------------------------------------6. EXTRACT CropScape------------------------------------

#----Converting data to CropScape crs----
species.spatial.crops = spTransform(species.spatial, crs(crops))

#----Extracting CropScape values----
crops.vx = velox(stack(crops))
spol = gBuffer(species.spatial.crops, width=350, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = crops.vx$extract(spdf, df=T)
ex.mat$ID_sp = as.numeric(ex.mat$ID_sp)

#----Calculating proportional land cover----
p.crops = ex.mat %>%
  setNames(c("ID", "class")) %>%
  group_by(ID, class) %>%
  summarise(n = n()) %>%
  mutate(pland = n / sum(n)) %>%
  ungroup() %>%
  select(ID, class, pland) %>%
  tidyr::complete(ID, nesting(class), fill = list(pland = 0)) %>% # Fill in implicit landcover 0s
  spread(class, pland)

colnames(p.crops) = c("ID", paste0("cdl", colnames(p.crops[,2:length(colnames(p.crops))])))

#------------------------------------7. EXTRACT nlcd------------------------------------

#----Converting data to NLCD crs----
species.spatial.nlcd = spTransform(species.spatial, crs(nlcd))

#----Extracting NLCD values----
nlcd.vx = velox(stack(nlcd))
spol = gBuffer(species.spatial.nlcd, width=350, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = nlcd.vx$extract(spdf, df=T)
ex.mat$ID_sp = as.numeric(ex.mat$ID_sp)

#----Calculating proportional land cover----
p.nlcd = ex.mat %>%
  setNames(c("ID", "class")) %>%
  group_by(ID, class) %>%
  summarise(n = n()) %>%
  mutate(pland = n / sum(n)) %>%
  ungroup() %>%
  select(ID, class, pland) %>%
  tidyr::complete(ID, nesting(class), fill = list(pland = 0)) %>% # Fill in implicit landcover 0s
  spread(class, pland)

colnames(p.nlcd) = c("ID", paste0("nlcd", colnames(p.nlcd[,2:length(colnames(p.nlcd))])))
