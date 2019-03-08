# Gates Dupont                 #
# Eastern Meadowlark           #
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
  filter(year > 2013) %>%
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
#n.p.pts = round(0.7*(length(pres$lat))) # 70% training 30% testing
#bg = as.data.frame(spsample(st.contour, n = n.p.pts*10 + length(pres$lat)-n.p.pts,"random"))
bg = as.data.frame(spsample(st.contour, n = length(pres$lat)*13,"random"))
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
crops = crop(CropScape$MA2017, extent(spTransform(study.extent.corners, crs(CropScape$MA2017)))*1.1)

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
climate = getData('worldclim', var='bio', res=2.5, lat=ll[1], lon=ll[2])
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
p.crops = p.crops[,-1]

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
p.nlcd = p.nlcd[,-1]

#------------------------------------8. Extract Elevation------------------------------------

#----Converting data to elevation crs----
species.spatial.elevation = spTransform(species.spatial, crs(elevation))
elevation.extVals = raster::extract(elevation, species.spatial.elevation)

#------------------------------------9. Extract Climate------------------------------------

#----Converting data to climate crs----
species.spatial.climate = spTransform(species.spatial, crs(climate))
climate.extVals = data.frame(raster::extract(climate, species.spatial.climate))

#----Bringing everything together----
species.df = as.data.frame(species.spatial@coords)
colnames(species.df) = c("long", "lat")
species.df = cbind(pa = species.spatial@data[,3], species.df, p.crops, p.nlcd, elevation.extVals, climate.extVals)
species.df = species.df[complete.cases(species.df),]

#------------------------------------10. Training and Testing------------------------------------

#----Presence Data----
pres = species.df[species.df$pa==1,] # Subsetting for presence data
n.pres = length(pres[,1]) # Number of pres points
n.pres.train = round(length(pres[,1])*0.7) # Number of pres points for training
n.pres.test = n.pres - n.pres.train # Number of pres points for tetsing

# Split up pres for train and test
train.pres = pres[1:n.pres.train,]
test.pres = pres[(n.pres.train+1):(n.pres.train+n.pres.test),]

# Checking the spatial distribution of pres
plot(st.contour)
points(lat~long, data=train.pres, col="blue", pch=20)
points(lat~long, data=test.pres, col="red", pch=20)

#----Absence Data----
abs = species.df[species.df$pa==0,] # Absence data
n.abs = length(abs[,1]) # Number of absence points
n.rows.samp = 1:n.abs # Vector for sampling absence data

#----Looping to create 10 datasets----
for(k in 1){
  n.rows.samp = 1:n.abs # Vector for sampling absence data
  
  # Train Datasets
  rf.train.dfs = vector("list", 10)
  for(i in 1:10){
    j = sample(n.rows.samp, n.pres.train, replace = F)
    abs.samp = abs[j,]
    rf.train.dfs[[i]] = rbind(train.pres, abs.samp)
    n.rows.samp = n.rows.samp[!n.rows.samp %in% j]
  }
  
  # Test Datasets
  rf.test.dfs = vector("list", 10)
  for(i in 1:10){
    j = sample(n.rows.samp, n.pres.test, replace = F)
    abs.samp = abs[j,]
    rf.test.dfs[[i]] = rbind(test.pres, abs.samp)
    n.rows.samp = n.rows.samp[!n.rows.samp %in% j]
  }
}

#------------------------------------11. Model------------------------------------

#----10 Random Forest models----
rf = vector("list", 10)
for(i in 1:10){
  rf[[i]] = rf(pa ~ ., rf.train.dfs[[i]])
}

#----Plotting variable importance----
par(mfrow=c(2,5))
for(i in 1:10){
  varImpPlot(rf[[i]], type=2)
}
par(mfrow=c(1,1))

#----Calculating AUC score for each model----
rf.evals = vector("list", 10)
rf.AUC.scores = vector("numeric", 10)
for(i in 1:10){
  test.df = rf.test.dfs[[i]]
  rf.evals[[i]] = evaluate(test.df[test.df$pa==1,], test.df[test.df$pa==0,], rf[[i]])
  rf.AUC.scores[i] = rf.evals[[i]]@auc
}
boxplot(rf.AUC.scores, ylim=c(0,1), main="Random Forest Model", ylab="AUC")

