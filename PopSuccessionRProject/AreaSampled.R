#This script calculates the sampling area where soil inocula
#was sourced from

#loads sf library
library(sf)

#Sets path for input files
Input_path="./InputFiles/"

#reads in gps coordinates of sampling locations
SamplingCoords=read.csv(paste(Input_path, "SamplingCoords.csv", sep=""))

# Convert to sf object (WGS84 CRS)
sampling_sf <- st_as_sf(SamplingCoords, coords = c("Longitude", "Latitude"), crs = 4326)

# Transform to UTM zone 12N (EPSG:32612) for southern Utah
sampling_utm <- st_transform(sampling_sf, crs = 32612)

# Create convex hull
convex_hull <- st_convex_hull(st_union(sampling_utm))

# Calculate area in square meters
area_m2 <- st_area(convex_hull)

# Convert to hectares if desired
area_ha <- as.numeric(area_m2) / 10000

# Print result
cat("Area (m²):", round(area_m2, 2), "\n")
cat("Area (ha):", round(area_ha, 2), "\n")
