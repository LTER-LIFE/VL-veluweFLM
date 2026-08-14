rm(list=ls())
setwd("/Users/mariagutierrezfuentesal/Documents/R")
library(terra)
library(sf)
library(tidyverse)
library(openxlsx)


## VELUWE POLYGON
natura2000_sf <- sf::read_sf("https://service.pdok.nl/rvo/natura2000/atom/downloads/natura2000.gpkg") %>%
  sf::st_transform(28992)

# Create the Veluwe sf by filtering the Natura2000 shapefile 
veluwe_sf <- natura2000_sf %>%
  dplyr::filter(objectid == 4049)


veluwe_v <- vect(veluwe_sf)


## FUNCTION: MAPPING IN SAME SCALE
plot_same_scale <- function(r1, r2, r3, depth_names, title, palette = "YlOrBr", probs = c(0, 1)) {
  
  r_stack <- c(r1, r2, r3)
  names(r_stack) <- depth_names
  
  lims <- quantile(values(r_stack), probs = probs, na.rm = TRUE)
  
  r_plot <- clamp(
    r_stack,
    lower = lims[1],
    upper = lims[2],
    values = TRUE
  )
  
  plot(
    r_plot,
    range = lims,
    col = hcl.colors(30, palette, rev = TRUE),
    nc = 3,
    main = paste(title, depth_names)
  )
}

## FUNCTION: CORRELATION BETWEEN SOIL DEPTHS

layer_correlation <- function(l1, l2, l3, l4, l5,l6, layer_names) {
  
  l_stack <- c(l1, l2, l3, l4, l5, l6)
  
  names(l_stack) <- layer_names
  
  l_values <- as.data.frame(l_stack)
  
  pearson <- cor(l_values,  method = "pearson")
  spearman <- cor(l_values,  method = "spearman")
  
  list(
    pearson = pearson,
    spearman = spearman
  )
  
}
#FUNCTION:PLOTS

plot_depth <- function(l1, l2, l3, l4, l5,l6, layer_names, y_label) {
  
  l_stack <- c(l1, l2, l3, l4,l5,l6)
  
  names(l_stack) <- layer_names
  
  profile_df <- data.frame()
  for (i in 1:6) {
    
    r <- l_stack[[i]]
    mat <- as.matrix(r, wide = TRUE)  
    x_coord <- xFromCol(r, 1:ncol(r))
    mean_value <- apply(mat, 2, mean, na.rm = TRUE)
    
    # make small table for this depth
    layer_df <- data.frame(
      coordinate = x_coord,
      value = mean_value,
      depth = names(l_stack)[i]
    )
    
    # add it to the full table
    profile_df <- rbind(profile_df, layer_df)
    
  }
  
  profile_df$depth <- factor(
    profile_df$depth,
    levels = layer_names
  )
  # make plot
  ggplot(profile_df, aes(x = coordinate, y = value, colour = depth)) +
    geom_line(linewidth = 1) +
    theme_minimal() +
    labs(
      x = "X coordinate",
      y = y_label,
      colour = "Soil depth")
  
}

## FUNCTION: CORRELATION ACCROSS YEARS
year_correlation <- function(l1, l2, l3, l4,l5, layer_names) {
  
  l_stack <- c(l1, l2, l3, l4,l5)
  
  names(l_stack) <- layer_names
  
  l_values <- as.data.frame(l_stack)
  
  pearson <- cor(l_values,  method = "pearson")
  spearman <- cor(l_values,  method = "spearman")
  
  list(
    pearson = pearson,
    spearman = spearman
  )
  
}

### SOIL TEXTURE 
##CLAY
#0-5
clay <- rast("internship_2/data_raw/clay/clay_per_d_0_5_QRF_pred_mean_processed.tif")
clay_veluwe <- crop(clay, veluwe_v)
clay_veluwe <- mask(clay_veluwe, veluwe_v)
#5-15 
clay_5_15 <- rast("internship_2/data_raw/clay/clay_per_d_5_15_QRF_pred_mean_processed.tif")
clay_veluwe_5_15 <- crop(clay_5_15, veluwe_v)
clay_veluwe_5_15 <- mask(clay_veluwe_5_15, veluwe_v)
#15-30
clay_15_30 <- rast("internship_2/data_raw/clay/clay_per_d_15_30_QRF_pred_mean_processed.tif")
clay_veluwe_15_30  <- crop(clay_15_30, veluwe_v)
clay_veluwe_15_30  <- mask(clay_veluwe_15_30, veluwe_v)
#30-60
clay_30_60 <- rast("internship_2/data_raw/clay/clay_per_d_30_60_QRF_pred_mean_processed.tif")
clay_veluwe_30_60 <- crop(clay_30_60, veluwe_v)
clay_veluwe_30_60 <- mask(clay_veluwe_30_60, veluwe_v)
#60-100
clay_60_100 <- rast("internship_2/data_raw/clay/clay_per_d_60_100_QRF_pred_mean_processed.tif")
clay_veluwe_60_100 <- crop(clay_60_100, veluwe_v)
clay_veluwe_60_100 <- mask(clay_veluwe_60_100, veluwe_v)
#100-200
clay_100_200 <- rast("internship_2/data_raw/clay/clay_per_d_100_200_QRF_pred_mean_processed.tif")
clay_veluwe_100_200 <- crop(clay_100_200, veluwe_v)
clay_veluwe_100_200 <- mask(clay_veluwe_100_200, veluwe_v)


#Maps
plot_same_scale(
  clay_veluwe,
  clay_veluwe_15_30,
  clay_veluwe_100_200,
  depth_names = c("0-5 cm", "15-30 cm", "100-200 cm"),
  title = "Clay content (%)")

#Correlation
clay_cor <-layer_correlation(clay_veluwe,
                             clay_veluwe_5_15,
                             clay_veluwe_15_30,
                             clay_veluwe_30_60,
                             clay_veluwe_60_100,  
                             clay_veluwe_100_200,
                             c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"))
clay_output <- clay_cor$pearson

#Soil depth plot

clay_depth_plot <- plot_depth(clay_veluwe,
           clay_veluwe_5_15,
           clay_veluwe_15_30,
           clay_veluwe_30_60,
           clay_veluwe_60_100,  
           clay_veluwe_100_200,
           c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"),y_label = "clay (%)")


##SILT
#0-5
silt <- rast("internship_2/data_raw/silt/silt_per_d_0_5_QRF_pred_mean_processed.tif")
silt_veluwe <- crop(silt,veluwe_v)
silt_veluwe <- mask(silt_veluwe,veluwe_v)
#5-15
silt_5_15 <- rast("internship_2/data_raw/silt/silt_per_d_5_15_QRF_pred_mean_processed.tif")
silt_veluwe_5_15 <- crop(silt_5_15, veluwe_v)
silt_veluwe_5_15 <- mask(silt_veluwe_5_15, veluwe_v)
#15-30
silt_15_30 <- rast("internship_2/data_raw/silt/silt_per_d_15_30_QRF_pred_mean_processed.tif")
silt_veluwe_15_30 <- crop(silt_15_30,veluwe_v)
silt_veluwe_15_30 <- mask(silt_veluwe_15_30,veluwe_v)
#30-60
silt_30_60 <- rast("internship_2/data_raw/silt/silt_per_d_30_60_QRF_pred_mean_processed.tif")
silt_veluwe_30_60 <- crop(silt_30_60, veluwe_v)
silt_veluwe_30_60 <- mask(silt_veluwe_30_60, veluwe_v)
#60-100
silt_60_100 <- rast("internship_2/data_raw/silt/silt_per_d_60_100_QRF_pred_mean_processed.tif")
silt_veluwe_60_100 <- crop(silt_60_100, veluwe_v)
silt_veluwe_60_100 <- mask(silt_veluwe_60_100, veluwe_v)
#100-200
silt_100_200 <- rast("internship_2/data_raw/silt/silt_per_d_100_200_QRF_pred_mean_processed.tif")
silt_veluwe_100_200 <- crop(silt_100_200, veluwe_v)
silt_veluwe_100_200 <- mask(silt_veluwe_100_200, veluwe_v)



#Maps
plot_same_scale(
  silt_veluwe,
  silt_veluwe_15_30,
  silt_veluwe_100_200,
  depth_names = c("0-5 cm", "15-30 cm", "100-200 cm"),
  title = "Silt content (%)")

#Correlation
silt_cor <-layer_correlation(silt_veluwe,
                             silt_veluwe_5_15,
                             silt_veluwe_15_30,
                             silt_veluwe_30_60,
                             silt_veluwe_60_100,  
                             silt_veluwe_100_200,
                             c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"))
silt_output <- silt_cor$pearson

#Soil depth plot
silt_depth_plot <- plot_depth(silt_veluwe,
           silt_veluwe_5_15,
           silt_veluwe_15_30,
           silt_veluwe_30_60,
           silt_veluwe_60_100,  
           silt_veluwe_100_200,
           c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"),
           y_label = "silt (%)")
##SAND
#0-5
sand <- rast("internship_2/data_raw/sand/sand_per_d_0_5_QRF_pred_mean_processed.tif")
sand_veluwe <- crop(sand,veluwe_v)
sand_veluwe <- mask(sand_veluwe,veluwe_v)
#5-15
sand_5_15 <- rast("internship_2/data_raw/sand/sand_per_d_5_15_QRF_pred_mean_processed.tif")
sand_veluwe_5_15 <- crop(sand_5_15, veluwe_v)
sand_veluwe_5_15 <- mask(sand_veluwe_5_15, veluwe_v)
#15-30
sand_15_30 <- rast("internship_2/data_raw/sand/sand_per_d_15_30_QRF_pred_mean_processed.tif")
sand_veluwe_15_30 <- crop(sand_15_30,veluwe_v)
sand_veluwe_15_30 <- mask(sand_veluwe_15_30,veluwe_v)
#30-60
sand_30_60 <- rast("internship_2/data_raw/sand/sand_per_d_30_60_QRF_pred_mean_processed.tif")
sand_veluwe_30_60 <- crop(sand_30_60, veluwe_v)
sand_veluwe_30_60 <- mask(sand_veluwe_30_60, veluwe_v)
#60-100
sand_60_100 <- rast("internship_2/data_raw/sand/sand_per_d_60_100_QRF_pred_mean_processed.tif")
sand_veluwe_60_100 <- crop(sand_60_100, veluwe_v)
sand_veluwe_60_100 <- mask(sand_veluwe_60_100, veluwe_v)
#100-200
sand_100_200 <- rast("internship_2/data_raw/sand/sand_per_d_100_200_QRF_pred_mean_processed.tif")
sand_veluwe_100_200 <- crop(sand_100_200, veluwe_v)
sand_veluwe_100_200 <- mask(sand_veluwe_100_200, veluwe_v)

#Maps
plot_same_scale(
  sand_veluwe,
  sand_veluwe_15_30,
  sand_veluwe_100_200,
  depth_names = c("0-5 cm", "15-30 cm", "100-200 cm"),
  title = "Sand content (%)"
)

#Correlation
sand_cor <-layer_correlation(sand_veluwe,
                             sand_veluwe_5_15,
                             sand_veluwe_15_30,
                             sand_veluwe_30_60,
                             sand_veluwe_60_100,  
                             sand_veluwe_100_200,
                             c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"))
sand_output <- sand_cor$pearson

#Soil depth plot
sand_depth_plot <- plot_depth(sand_veluwe,
           sand_veluwe_5_15,
           sand_veluwe_15_30,
           sand_veluwe_30_60,
           sand_veluwe_60_100,  
           sand_veluwe_100_200,
           c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"),y_label = "sand (%)")


##BULK DENSITY
#0-5
BD <- rast("internship_2/data_raw/bulk_density/BD_gcm3_d_0_5_QRF_pred_mean.tif")
BD_veluwe <- crop(BD,veluwe_v)
BD_veluwe <- mask(BD_veluwe,veluwe_v)
#15-30
BD_15_30 <- rast("internship_2/data_raw/bulk_density/BD_gcm3_d_15_30_QRF_pred_mean.tif")
BD_veluwe_15_30 <- crop(BD_15_30,veluwe_v)
BD_veluwe_15_30 <- mask(BD_veluwe_15_30,veluwe_v)
#100-200
BD_100_200 <- rast("internship_2/data_raw/bulk_density/BD_gcm3_d_100_200_QRF_pred_mean.tif")
BD_veluwe_100_200 <- crop(BD_100_200, veluwe_v)
BD_veluwe_100_200 <- mask(BD_veluwe_100_200, veluwe_v)

#Maps
plot_same_scale(
  BD_veluwe,
  BD_veluwe_15_30,
  BD_veluwe_100_200,
  depth_names = c("0-5 cm", "15-30 cm", "100-200 cm"),
  title = "Bulk density (g/cm3)"
)

#analysis was not continued, not much spatial variation and was not given importance in expert meeting

### SOIL CHEMISTRY 
##pH
#0-5 cm layer
pH <- rast("internship_2/data_raw/pH/pH_KCl_d_0_5_QRF_pred_mean.tif")
pH_veluwe <- crop(pH,veluwe_v)
pH_veluwe <- mask(pH_veluwe,veluwe_v)

#5-15
pH_5_15 <- rast("internship_2/data_raw/pH/pH_KCl_d_5_15_QRF_pred_mean.tif")
pH_veluwe_5_15 <- crop(pH_5_15, veluwe_v)
pH_veluwe_5_15 <- mask(pH_veluwe_5_15, veluwe_v)
#5-15
pH_5_15 <- rast("internship_2/data_raw/pH/pH_KCl_d_5_15_QRF_pred_mean.tif")
pH_veluwe_5_15 <- crop(pH_5_15,veluwe_v)
pH_veluwe_5_15 <- mask(pH_veluwe_5_15,veluwe_v)
# 15-30
pH_15_30 <- rast("internship_2/data_raw/pH/pH_KCl_d_15_30_QRF_pred_mean.tif")
pH_veluwe_15_30 <- crop(pH_15_30,veluwe_v)
pH_veluwe_15_30 <- mask(pH_veluwe_15_30,veluwe_v)
#30-60
pH_30_60 <- rast("internship_2/data_raw/pH/pH_KCl_d_30_60_QRF_pred_mean.tif")
pH_veluwe_30_60 <- crop(pH_30_60, veluwe_v)
pH_veluwe_30_60 <- mask(pH_veluwe_30_60, veluwe_v)
#60-100
pH_60_100 <- rast("internship_2/data_raw/pH/pH_KCl_d_60_100_QRF_pred_mean.tif")
pH_veluwe_60_100 <- crop(pH_60_100, veluwe_v)
pH_veluwe_60_100 <- mask(pH_veluwe_60_100, veluwe_v)
#100-200cm 
pH_100_200 <- rast("internship_2/data_raw/pH/pH_KCl_d_100_200_QRF_pred_mean.tif")
pH_veluwe_100_200 <- crop(pH_100_200,veluwe_v)
pH_veluwe_100_200 <- mask(pH_veluwe_100_200,veluwe_v)


#Mapping
plot_same_scale(
  pH_veluwe,
  pH_veluwe_15_30,
  pH_veluwe_100_200,
  depth_names = c("0-5 cm", "15-30 cm", "100-200 cm"),
  title = "pH",
  palette = "Blue-Red 3"
)

#Correlation
pH_cor <-layer_correlation(pH_veluwe,
                           pH_veluwe_5_15,
                           pH_veluwe_15_30,
                           pH_veluwe_30_60,
                           pH_veluwe_60_100,  
                           pH_veluwe_100_200,
                           c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"))
pH_output <- pH_cor$pearson

# Soil depth plot
pH_depth_plot <- plot_depth(pH_veluwe,
           pH_veluwe_5_15,
           pH_veluwe_15_30,
           pH_veluwe_30_60,
           pH_veluwe_60_100,  
           pH_veluwe_100_200,
           c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"), y_label ="pH")

##TOTAL N
#0-5 depth
N <- rast("internship_2/data_raw/total_N/N_tot_mgkg_d_0_5_QRF_pred_mean.tif")
N_veluwe <- crop(N,veluwe_v)
N_veluwe <- mask(N_veluwe,veluwe_v)
#5-15 depth
N_5_15 <- rast("internship_2/data_raw/total_N/N_tot_mgkg_d_5_15_QRF_pred_mean.tif")
N_veluwe_5_15 <- crop(N_5_15,veluwe_v)
N_veluwe_5_15 <- mask(N_veluwe_5_15,veluwe_v)
#15-30 depth
N_15_30 <- rast("internship_2/data_raw/total_N/N_tot_mgkg_d_15_30_QRF_pred_mean.tif")
N_veluwe_15_30 <- crop(N_15_30,veluwe_v)
N_veluwe_15_30 <- mask(N_veluwe_15_30,veluwe_v)
#30-60
N_30_60 <- rast("internship_2/data_raw/total_N/N_tot_mgkg_d_30_60_QRF_pred_mean.tif")
N_veluwe_30_60 <- crop(N_30_60, veluwe_v)
N_veluwe_30_60 <- mask(N_veluwe_30_60, veluwe_v)
#60-100
N_60_100 <- rast("internship_2/data_raw/total_N/N_tot_mgkg_d_60_100_QRF_pred_mean.tif")
N_veluwe_60_100 <- crop(N_60_100, veluwe_v)
N_veluwe_60_100 <- mask(N_veluwe_60_100, veluwe_v)
#100-200 depth
N_100_200 <- rast("internship_2/data_raw/total_N/N_tot_mgkg_d_100_200_QRF_pred_mean.tif")
N_veluwe_100_200 <- crop(N_100_200,veluwe_v)
N_veluwe_100_200 <- mask(N_veluwe_100_200,veluwe_v)


#Maps
plot_same_scale(
  N_veluwe,
  N_veluwe_15_30,
  N_veluwe_100_200,
  depth_names = c("0-5 cm", "15-30 cm", "100-200 cm"),
  title = "total N (mg/kg)",
  palette = "Blue-Red 3"
)

#Correlation
N_cor <-layer_correlation(N_veluwe,
                          N_veluwe_5_15,
                          N_veluwe_15_30,
                          N_veluwe_30_60,
                          N_veluwe_60_100,  
                          N_veluwe_100_200,
                          c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"))
N_output <- N_cor$pearson

# Soil depth plot
N_depth_plot <- plot_depth(N_veluwe,
           N_veluwe_5_15,
           N_veluwe_15_30,
           N_veluwe_30_60,
           N_veluwe_60_100,  
           N_veluwe_100_200,
           c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"), y_label="Total N")
##CEC
#0-5 depth
CEC <- rast("internship_2/data_raw/cation_exchange_capacity/CEC_d_0_5_QRF_pred_mean.tif")
CEC_veluwe <- crop(CEC,veluwe_v)
CEC_veluwe <- mask(CEC_veluwe,veluwe_v)
#5-15
CEC_5_15 <- rast("internship_2/data_raw/cation_exchange_capacity/CEC_d_5_15_QRF_pred_mean.tif")
CEC_veluwe_5_15 <- crop(CEC_5_15, veluwe_v)
CEC_veluwe_5_15 <- mask(CEC_veluwe_5_15, veluwe_v)
#15-30 depth
CEC_15_30 <- rast("internship_2/data_raw/cation_exchange_capacity/CEC_d_15_30_QRF_pred_mean.tif")
CEC_veluwe_15_30 <- crop(CEC_15_30,veluwe_v)
CEC_veluwe_15_30 <- mask(CEC_veluwe_15_30,veluwe_v)
#30-60
CEC_30_60 <- rast("internship_2/data_raw/cation_exchange_capacity/CEC_d_30_60_QRF_pred_mean.tif")
CEC_veluwe_30_60 <- crop(CEC_30_60, veluwe_v)
CEC_veluwe_30_60 <- mask(CEC_veluwe_30_60, veluwe_v)
#60-100
CEC_60_100 <- rast("internship_2/data_raw/cation_exchange_capacity/CEC_d_60_100_QRF_pred_mean.tif")
CEC_veluwe_60_100 <- crop(CEC_60_100, veluwe_v)
CEC_veluwe_60_100 <- mask(CEC_veluwe_60_100, veluwe_v)
#100-200 depth
CEC_100_200 <- rast("internship_2/data_raw/cation_exchange_capacity/CEC_d_100_200_QRF_pred_mean.tif")
CEC_veluwe_100_200 <- crop(CEC_100_200,veluwe_v)
CEC_veluwe_100_200 <- mask(CEC_veluwe_100_200,veluwe_v)

#Mapping
plot_same_scale(
  CEC_veluwe,
  CEC_veluwe_15_30,
  CEC_veluwe_100_200,
  depth_names = c("0-5 cm", "15-30 cm", "100-200 cm"),
  title = "CEC (mmol(c)/kg)",
  palette = "Blue-Red 3"
)

#Correlation
CEC_cor <-layer_correlation(CEC_veluwe,
                            CEC_veluwe_5_15,
                            CEC_veluwe_15_30,
                            CEC_veluwe_30_60,
                            CEC_veluwe_60_100,  
                            CEC_veluwe_100_200,
                            c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"))
CEC_output <- CEC_cor$pearson

# Soil depth map
CEC_depth_plot <- plot_depth(CEC_veluwe,
           CEC_veluwe_5_15,
           CEC_veluwe_15_30,
           CEC_veluwe_30_60,
           CEC_veluwe_60_100,  
           CEC_veluwe_100_200,
           c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"), y_label="CEC")


##SOM (2023)

#0-5 depth
SOM <- rast("internship_2/data_raw/soil_organic_matter/2023/SOM_per_2023_d_0_5_QRF_pred_mean.tif")
SOM_veluwe <- crop(SOM,veluwe_v)
SOM_veluwe <- mask(SOM_veluwe,veluwe_v)

#15-30 depth
SOM_15_30 <- rast("internship_2/data_raw/soil_organic_matter/2023/SOM_per_2023_d_15_30_QRF_pred_mean.tif")
SOM_veluwe_15_30 <- crop(SOM_15_30,veluwe_v)
SOM_veluwe_15_30 <- mask(SOM_veluwe_15_30,veluwe_v)

#100-200
SOM_100_200 <- rast("internship_2/data_raw/soil_organic_matter/2023/SOM_per_2023_d_100_200_QRF_pred_mean.tif")
SOM_veluwe_100_200 <- crop(SOM_100_200,veluwe_v)
SOM_veluwe_100_200 <- mask(SOM_veluwe_100_200,veluwe_v)

#Mapping
plot_same_scale(
  SOM_veluwe,
  SOM_veluwe_15_30,
  SOM_veluwe_100_200,
  depth_names = c("0-5 cm", "15-30 cm", "100-200 cm"),
  title = "SOM content (%)"
)


#2000
SOM_2000 <- rast("internship_2/data_raw/soil_organic_matter/2000/SOM_per_2000_d_0_5_QRF_pred_mean.tif")
SOM_2000_veluwe <- crop(SOM_2000,veluwe_v)
SOM_2000_veluwe <- terra::mask(SOM_2000_veluwe,veluwe_v)
#1980
SOM_1980 <- rast("internship_2/data_raw/soil_organic_matter/1980/SOM_per_1980_d_0_5_QRF_pred_mean.tif")
SOM_1980_veluwe <- crop(SOM_1980,veluwe_v)
SOM_1980_veluwe <- terra::mask(SOM_1980_veluwe,veluwe_v)
#1960
SOM_1960 <- rast("internship_2/data_raw/soil_organic_matter/1960/SOM_per_1960_d_0_5_QRF_pred_mean.tif")
SOM_1960_veluwe <- crop(SOM_1960,veluwe_v)
SOM_1960_veluwe <- terra::mask(SOM_1960_veluwe,veluwe_v)
#1953
SOM_1953 <- rast("internship_2/data_raw/soil_organic_matter/1953/SOM_per_1953_d_0_5_QRF_pred_mean.tif")
SOM_1953_veluwe <- crop(SOM_1953,veluwe_v)
SOM_1953_veluwe <- terra::mask(SOM_1953_veluwe,veluwe_v)
#2023
SOM_2023 <- rast("internship_2/data_raw/soil_organic_matter/2023/SOM_per_2023_d_0_5_QRF_pred_mean.tif")
SOM_2023_veluwe <- crop(SOM_2023,veluwe_v)
SOM_2023_veluwe <- terra::mask(SOM_2023_veluwe,veluwe_v)

#Correlation
SOM_year_cor <-year_correlation(SOM_2023_veluwe,
                                 SOM_2000_veluwe,
                                 SOM_1980_veluwe,
                                 SOM_1960_veluwe,
                                SOM_1953_veluwe,
                             c("2023","2000", "1980","1960","1953"))
SOM_year_output <- SOM_year_cor$pearson

#picking 2000 as its closest to initializing year
#5-15
SOM_2000_5_15 <- rast("internship_2/data_raw/soil_organic_matter/2000/SOM_per_2000_d_5_15_QRF_pred_mean.tif")
SOM_veluwe_2000_5_15 <- crop(SOM_2000_5_15,veluwe_v)
SOM_veluwe_2000_5_15 <- mask(SOM_veluwe_2000_5_15,veluwe_v)
#15-30 depth
SOM_2000_15_30 <- rast("internship_2/data_raw/soil_organic_matter/2000/SOM_per_2000_d_15_30_QRF_pred_mean.tif")
SOM_veluwe_2000_15_30 <- crop(SOM_2000_15_30,veluwe_v)
SOM_veluwe_2000_15_30 <- mask(SOM_veluwe_2000_15_30,veluwe_v)
#30-60
SOM_2000_30_60 <- rast("internship_2/data_raw/soil_organic_matter/2000/SOM_per_2000_d_30_60_QRF_pred_mean.tif")
SOM_veluwe_2000_30_60 <- crop(SOM_2000_30_60,veluwe_v)
SOM_veluwe_2000_30_60 <- mask(SOM_veluwe_2000_30_60,veluwe_v)
#60-100
SOM_2000_60_100 <- rast("internship_2/data_raw/soil_organic_matter/2000/SOM_per_2000_d_60_100_QRF_pred_mean.tif")
SOM_veluwe_2000_60_100 <- crop(SOM_2000_60_100,veluwe_v)
SOM_veluwe_2000_60_100 <- mask(SOM_veluwe_2000_60_100,veluwe_v)
#100-200
SOM_2000_100_200 <- rast("internship_2/data_raw/soil_organic_matter/2000/SOM_per_2000_d_100_200_QRF_pred_mean.tif")
SOM_veluwe_2000_100_200 <- crop(SOM_2000_100_200,veluwe_v)
SOM_veluwe_2000_100_200 <- mask(SOM_veluwe_2000_100_200,veluwe_v)


#Correlation
SOM_cor <-layer_correlation(SOM_2000_veluwe,
                            SOM_veluwe_2000_5_15,
                            SOM_veluwe_2000_15_30,
                            SOM_veluwe_2000_30_60,
                            SOM_veluwe_2000_60_100,  
                            SOM_veluwe_2000_100_200,
                            c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"))
SOM_output <- SOM_cor$pearson

# Soil depth map
SOM_depth_plot <- plot_depth(SOM_2000_veluwe,
           SOM_veluwe_2000_5_15,
           SOM_veluwe_2000_15_30,
           SOM_veluwe_2000_30_60,
           SOM_veluwe_2000_60_100,  
           SOM_veluwe_2000_100_200,  
           c("0-5 cm","5-15 cm", "15-30 cm","30-60 cm","60-100cm","100-200 cm"), y_label="SOM")


##SAVE OUTPUTS
cor_matrices <- list(
  SOM_year = SOM_year_output,
  SOM_2000 = SOM_output,
  CEC = CEC_output,
  Clay = clay_output,
  N = N_output,
  pH = pH_output,
  Sand = sand_output,
  Silt = silt_output
)

cor_matrices <- lapply(cor_matrices, function(mat) {
  df <- as.data.frame(mat)
  df <- cbind(row_name = rownames(df), df)
  rownames(df) <- NULL
  df
})
write.xlsx(cor_matrices, file = "/Users/mariagutierrezfuentesal/Documents/R/internship_2/outputs/correlation_matrices.xlsx")

#ALL DEPTH PLOTS
library(patchwork)

wrap_plots(SOM_depth_plot, N_depth_plot, pH_depth_plot,CEC_depth_plot,clay_depth_plot, sand_depth_plot, ncol = 2)



  (SOM_depth_plot + N_depth_plot) /
  (pH_depth_plot  + CEC_depth_plot) /
  (clay_depth_plot+ sand_depth_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")




