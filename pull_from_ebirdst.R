# install.packages("ebirdst")
library(ebirdst)
library(raster)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(viridisLite)
# handle namespace conflicts
extract <- raster::extract


"STATE DATA"

# Get all data
us <- getData("GADM", country="USA", level=1)

# Subset to MA only
ma = us[match(toupper("Massachusetts"),toupper(us$NAME_1)),]


"GETTING EBIRD DATA"

# Download data
sp_path <- ebirdst_download(species =  "easmea", tifs_only = FALSE)


"ABUNDANCE"

# load trimmed mean abundances
abunds <- load_raster("abundance_umean", path = sp_path)

# Cropping to an area rougly the size of MA
abunds_18_cr =  crop(abunds[[18]], extent(c(-6.2e6, -5.6e6,  4.5e6, 4.85e6)))

# define mollweide projection
mollweide <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"

# project single layer from stack to mollweide
week18_moll <- projectRaster(abunds_18_cr, crs = mollweide, method = "ngb")

# Mask to MA and crop
ma_moll = spTransform(ma, mollweide)
r = mask(week18_moll, ma_moll)
r = crop(r, ma_moll)

# map single layer with full annual extent
plot(r, 
     col = abundance_palette(10, season = "breeding"), 
     axes = F, box = F, 
     maxpixels = ncell(r),
     main = "Eastern Meadowlark, June 1, MA \n (eBird Status and Trends)")


"VARIABLE IMPORTANCE"

pis <- load_pis(sp_path)
pds <- load_pds(sp_path)

# Select region and season
lp_extent <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
                            t = c(0.425, 0.475))
map_centroids(sp_path, ext = lp_extent)

# aggregating fragstats for cover classes
plot_pis(pis, ext = lp_extent, by_cover_class = TRUE, n_top_pred = 15)


