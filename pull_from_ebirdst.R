# Gates Dupont        #
# dupont97@gmail.com  #
# Feb. 4, 2020        #
# # # # # # # # # # # # 

# install.packages("ebirdst")
library(ebirdst)
library(raster)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(ggpubr)
library(viridisLite)
library(dplyr)
# handle namespace conflicts
extract <- raster::extract


"STATE DATA"

# Get all data
us <- getData("GADM", country="USA", level=1)

# Subset to MA only
ma = us[match(toupper("Massachusetts"),toupper(us$NAME_1)),]


"GETTING EBIRD DATA"

# Download data (this takes time (~20 mins for me), and space)
sp_path <- ebirdst_download(species =  "easmea", tifs_only = FALSE)


"ABUNDANCE"

# load trimmed mean abundances
abunds <- load_raster("abundance_umean", path = sp_path)

# Cropping to an area rougly the size of MA
# (Week 23 = June 6-12)
abunds_23_cr =  crop(
  abunds[[23]], 
  extent(c(-6.2e6, -5.6e6,  4.5e6, 4.85e6)))

# define mollweide projection
mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"

# project single layer from stack to mollweide
week23_moll <- projectRaster(
  abunds_23_cr, 
  crs = mollweide, method = "ngb")

# Mask to MA and crop
ma_moll = spTransform(ma, mollweide)
r = mask(week23_moll, ma_moll) %>%
  crop(., ma_moll) %>%
  projectRaster(., crs = crs(ma), method = "ngb")


"VARIABLE IMPORTANCE"

# Load predictor importance
pis = load_pis(sp_path)

# Select region and season
lp_extent <- ebirdst_extent(
  st_as_sf(ma), 
  t = c("2016-06-06", "2016-06-12")) # Models assumed 2016

# Plot centroids and extent of analysis
par(mfrow = c(1, 1), mar = c(0, 0, 0, 6))
calc_effective_extent(sp_path, ext = lp_extent)


"FINAL PLOT"

# Convert raster to data frame for ggplot
r_spdf <- as(r, "SpatialPixelsDataFrame")
r_df <- as.data.frame(r_spdf)
colnames(r_df) <- c("value", "x", "y")

# RELATIVE ABUNDANCE
ggplot() +
  geom_raster(data = r_df , aes(x = x, y = y, fill = value)) + 
  scale_fill_gradientn(colors = abundance_palette(10, season = "breeding")) +
  coord_quickmap() +
  theme_bw() +
  # theme(legend.position = "right", # For exporting as png of 2000x1000
  #       plot.title = element_text(size = 80),
  #       plot.subtitle = element_text(size = 50),
  #       plot.caption = element_text(size = 50),
  #       axis.text = element_text(size = 40),
  #       axis.title = element_text(size = 50),
  #       legend.title = element_text(size = 50),
  #       legend.text = element_text(size = 40)) +
  # guides(fill = guide_colourbar(barwidth = 4, barheight = 40)) +
  labs(title = "Eastern Meadowlark",
       subtitle = "Relative Abundance, June 6-12",
       caption = "Data source: eBird Status and Trends",
       fill = "RA", y = "Latitude", x = "Longitude")

# VARIABLE IMPORTANCE
plot_pis(pis, ext = lp_extent, by_cover_class = TRUE, n_top_pred = 15)
# I modified this function separately to produce stylistically similar plots
