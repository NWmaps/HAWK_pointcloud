#(HAWK job intervew exercise - LiDAR pointcloud processing)
setwd("P:/datahawk/lidarhawk/")

#Load packages
#install.packages(c("lidR", "sf", "terra", "ggplot2", "dplyr",'RCSF'))
library(lidR)
library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(RCSF)

#load data
las <- readLAS("als_data.laz")
boundary <- st_read("stand_boundary.gpkg")

#Explore data and preprocessing####
# Check summary statistics
summary(las)

# visualizae
plot(las, color = "Z")

# optional- Identify and remove extreme values (outliers) 
las <- filter_poi(las, Z > quantile(las@data$Z, 0.01) & Z < quantile(las@data$Z, 0.99))

#check CRS is the same in both files
# Compare CRS, reproject if needed
las_crs <- st_crs(las)
boundary_crs <- st_crs(boundary)

if (las_crs != boundary_crs) {
  boundary <- st_transform(boundary, las_crs)
}

# Compute polygon area in hectares
polygon_area_ha <- as.numeric(st_area(boundary) / 10000)


#clip las
las_clipped <- clip_roi(las, boundary)

# Compute basic metrics
metrics <- cloud_metrics(las_clipped, .stdmetrics)
#write.csv(metrics, "forest_lidar_metrics.csv") #save csv (optional)
#z values are over 400, which means the cloud is not normalized.

#Normalize to ground
# Classify ground points using a filtering algorithm
las_ground <- classify_ground(las_clipped, csf())
#terrain model from ground points
dtm <- rasterize_terrain(las_ground, res = 1, tin())
# Normalized cloud
las_norm <- normalize_height(las_clipped, dtm)
# Canopy Height Model from the normalized cloud
chm <- rasterize_canopy(las_norm, res = 1, p2r())



#Processing
#1. find tree tops####
ttops <- locate_trees(las_clipped, lmf(ws = 5))
plot(las_clipped, bg = "white", size = 2)
add_treetops3d(ttops)


#2. Determine height of all trees ####
tree_heights <- extract(chm, ttops)
ttops$Height <- tree_heights
ttops$Height <- as.numeric(ttops$Height$Z) #ensure is numeric
#check tree heights
summary(ttops$Height)
plot(ttops$Height)
# In general looks good, most trees around 25 m. 3 trees found to be less than a meter height. Artifacts? depending on the aim could explore the cloud visually and remove values if necessary

#check N/ha
tree_density <- nrow(ttops) / polygon_area_ha # 315 trees /ha. = normal value
tree_density


#3.Estimate DBH####
ttops$BHD <- (ttops$Height / 47)^2
#check DBH values
summary(ttops$BHD)
plot(ttops$BHD)
#most values around 30cm. For 25m height trees it makes sense. 2 near zero values


# 4. Create a histogram of dbh distribution ####
# (not sure which sort of figure "Abbildung" refers to)
ggplot(ttops, aes(x = BHD)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Diameter Distribution",
    x = "BHD (m) (5 cm intervals)",
    y = "Freq"
  )



# 5. calculate basal area####

ttops$BasalArea <- pi * ((ttops$BHD) / 2)^2
# Sum of basal areas per hectare
basal_area_ha <- sum(ttops$BasalArea) / (st_area(boundary) / 10000)
basal_area_ha


#summary table
summary_table <- data.frame(
  Total_Trees = nrow(ttops),
  Mean_Height = mean(ttops$Height, na.rm = TRUE),
  Mean_BHD = mean(ttops$BHD, na.rm = TRUE),
  Total_Basal_Area_ha = basal_area_ha
)

write.csv(summary_table, "lidar_summary.csv") #save csv (optional)
