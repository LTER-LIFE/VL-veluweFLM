# Relationship age, height, dbh for Scots pine
# Author: Stefan Vriend
# Created: 2026-04-14

# Sources:
# 'Groei en productie van grove den in Nederland', Jansen et al. 2018, https://edepot.wur.nl/444090
# 'Opbrengsttabellen Nederland 2018', Jansen & Oosterbaan 2018, https://edepot.wur.nl/460211

# Retrieve NBI-7 data from Probos website
url <- "https://probos.nl/DB/NBI-7_openbare_database.zip"
download.file(url, "NBI.zip", mode = "wb")
unzip("NBI.zip", exdir = "NBI")
nbi <- RODBC::odbcConnectAccess2007("NBI/NBI7_openbare database_definitief_UpdateJuli2024.accdb")

nbi7_proefbomen <- RODBC::sqlFetch(nbi, "data_NBI7_proefbomen")
nbi7_plotmetingen <- RODBC::sqlFetch(nbi, "data_NBI7_plotmetingen")
nbi7_plotdefinitie <- RODBC::sqlFetch(nbi, "data_NBI7_plotdefinitie")

# Retrieve Veluwe Natura2000 polygon
veluwe <- sf::read_sf("https://service.pdok.nl/rvo/natura2000/atom/downloads/natura2000.gpkg") |> 
  sf::st_transform(4326) |> 
  dplyr::filter(objectid == 4049) |> 
  dplyr::select("naam_n2k", "geom")

# Select NBI plots within Veluwe
veluwe_plots <- nbi7_plotdefinitie |> 
  dplyr::rename("plotnummer" = "PLOT") |> 
  dplyr::mutate(
    lat = as.integer(stringr::str_extract(`INSPIRE grid`, "(?<=N)\\d+")) * 1000 + 500,
    lon = as.integer(stringr::str_extract(`INSPIRE grid`, "(?<=E)\\d+")) * 1000 + 500,
  ) |> 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 3035) |> 
  sf::st_transform(crs = 4326) |> 
  {\(.) dplyr::mutate(.,
                      x = sf::st_coordinates(.)[,1], 
                      y = sf::st_coordinates(.)[,2])}() |> 
  sf::st_intersection(veluwe) |> 
  dplyr::left_join(nbi7_plotmetingen |> 
                     dplyr::select("plotnummer", "kiemjaar"), 
                   by = "plotnummer")

# Associate age (based on sprouting year, kiemjaar) with sample trees
veluwe_trees <- nbi7_proefbomen |> 
  dplyr::right_join(veluwe_plots |> 
                      dplyr::select("plotnummer", "kiemjaar"),
                    by = "plotnummer") |> 
  dplyr::filter(boomsoort == "GD", !is.na(kiemjaar)) |> 
  dplyr::rename("h" = "boomhoogte_m") |> 
  dplyr::mutate("t" = 2026 - kiemjaar)

# Function to calculate h_top, equals formula 10 in Jansen et al. 2018 (https://edepot.wur.nl/444090)
# c1, c2, c3 derived from Table 3, page 12
h_top <- function(t, 
                  h70,
                  c1 = 1.5162, 
                  c2 = 2500.5706, 
                  c3 = 22.0801) {
  
  Z <- h70 - c3
  R <- Z + sqrt(Z^2 + (2 * c2 * h70) / (70^c1))
  h <- h70 * ((t^c1) * ((70^c1) * R + c2)) / ((70^c1) * ((t^c1) * R + c2))
  
}

# Estimate h70 for Veluwe
# h70 is a site index and measure for 'boniteit'
mod <- nls(h ~ h70 * ((t^1.5162) * ((70^1.5162) * ((h70 - 22.0801) + sqrt((h70 - 22.0801)^2 + (2 * 2500.5706 * h70)/(70^1.5162))) + 2500.5706)) / ((70^1.5162) * ((t^1.5162) * ((h70 - 22.0801) + sqrt((h70 - 22.0801)^2 + (2 * 2500.5706 * h70)/(70^1.5162))) + 2500.5706)),
           data = veluwe_trees , start = list(h70 = 20),
           control = list(maxiter = 100))

# Plot estimated h_top and actual height measurements for Veluwe trees
veluwe_trees |> 
  dplyr::mutate(htop = h_top(t, h70 = coefficients(mod))) |> 
  ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = t, y = h)) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t, y = htop), color = "red") +
  ggplot2::theme_classic()

# Look-up table for relation between boniteit and h70, based on Jansen et al. 2018 (https://edepot.wur.nl/444090)
scots_pine_boni <-
  tibble::tibble(
    boniteit = c("I", "II", "III", "IV", "V"),
    h70 = c(23.8, 20.8, 17.8, 14.8, 11.8)
  )

# Look-up table for diameter, height, age relation as given by Appendix 1 (Bijlage 1) in Jansen et al. 2018 (https://edepot.wur.nl/444090)
scots_pine <- tibble::tibble(
  boniteit = "I",
  t = seq(5, 150, by = 5),
  d = c(1.3, 5.2, 8.1, 9.6, 12.3, 14.8, 17.1, 19.1, 21, 22.8, 24.4, 26.1, 27.8,
        29.3, 30.8, 32.2, 33.6, 34.8, 36.3, 37.8, 39.3, 40.7, 42.2, 43.6,
        44.9, 46.3, 47.7, 49, 50.4, 51.7),
  h = c(1.6, 4.3, 7.2, 9.9, 12.3, 14.4, 16.1, 17.6, 18.8, 19.9, 20.8,
        21.6, 22.3, 22.9, 23.4, 23.9, 24.3, 24.7, 25, 25.4, 25.7,
        25.9, 26.2, 26.4, 26.6, 26.8, 27, 27.2, 27.4, 27.5)
) |> 
  dplyr::add_row(
    boniteit = "II",
    t = seq(5, 150, by = 5),
    d = c(0.3, 3.5, 6.4, 8.3, 9.8, 12, 14, 15.9, 17.6, 19.2, 20.6, 22.2,
          23.8, 25.3, 26.7, 28, 29.3, 30.5, 31.9, 33.3, 34.6, 36,
          37.3, 38.6, 39.9, 41.2, 42.5, 43.8, 45, 46.3),
    h = c(1.1, 3, 5.3, 7.5, 9.5, 11.4, 13, 14.4, 15.7, 16.7, 17.7, 18.5,
          19.2, 19.9, 20.5, 21, 21.4, 21.9, 22.3, 22.6, 22.9, 23.2,
          23.5, 23.8, 24, 24.2, 24.4, 24.6, 24.8, 25)
  ) |> 
  dplyr::add_row(
    boniteit = "III",
    t = seq(5, 150, by = 5),
    d = c(NA, 2.1, 4.4, 6.5, 8.1, 9.2, 11, 12.7, 14.2, 15.7, 17, 18.5, 
          19.9, 21.3, 22.6, 23.8, 25, 26.2, 27.5, 28.8, 30, 31.3, 32.5,
          33.8, 35, 36.2, 37.4, 38.6, 39.8, 40.9),
    h = c(NA, 2.1, 3.7, 5.4, 7.1, 8.6, 10.1, 11.4, 12.6, 13.6, 14.6,
          15.4, 16.2, 16.9, 17.5, 18.1, 18.6, 19, 19.5, 19.9, 20.2, 20.6,
          20.9, 21.2, 21.4, 21.7, 21.9, 22.2, 22.4, 22.6)
  ) |> 
  dplyr::add_row(
    boniteit = "IV",
    t = seq(5, 150, by = 5),
    d = c(NA, 0.9, 2.8, 4.5, 6, 7.5, 8.5, 9.6, 11, 12.3, 13.6, 14.9,
          16.2, 17.4, 18.6, 19.7, 20.9, 21.9, 23.1, 24.3, 25.5, 26.6,
          27.8, 29, 30.1, 31.2, 32.4, 33.5, 34.6, 35.7),
    h = c(NA, 1.4, 2.5, 3.7, 5, 6.3, 7.5, 8.6, 9.7, 10.6, 11.6, 12.4,
          13.2, 13.9, 14.5, 15.1, 15.7, 16.2, 16.7, 17.1, 17.5, 17.9,
          18.2, 18.6, 18.9, 19.2, 19.5, 19.7, 20, 20.2)
  ) |> 
  dplyr::add_row(
    boniteit = "V",
    t = seq(5, 150, by = 5),
    d = c(NA, NA, 1.3, 2.7, 4, 5.2, 6.3, 7.3, 8.2, 9, 10.1, 11.3, 12.6,
          13.7, 14.8, 15.9, 16.9, 17.8, 18.9, 19.9, 21, 22.1, 23.1,
          24.2, 25.2, 26.3, 27.3, 28.4, 29.4, 30.5),
    h = c(NA, NA, 1.6, 2.5, 3.4, 4.3, 5.3, 6.2, 7.1, 7.9, 8.7, 9.5,
          10.2, 10.9, 11.5, 12.1, 12.7, 13.2, 13.7, 14.2, 14.6, 15.1, 
          15.5, 15.8, 16.2, 16.5, 16.8, 17.1, 17.4, 17.7)
  ) |> 
  dplyr::left_join(scots_pine_boni, by = "boniteit")

ggplot2::ggplot(scots_pine) +
  ggplot2::geom_point(mapping = ggplot2::aes(x = t, y = d, color = boniteit)) +
  ggplot2::geom_smooth(mapping = ggplot2::aes(x = t, y = d, color = boniteit, group = boniteit)) +
  ggplot2::theme_classic()

# Find which boniteit and associated growth curve is closest to estimated h70
boni <- scots_pine_boni |> 
  dplyr::mutate(dif = coefficients(mod) - h70) |> 
  dplyr::filter(dif > 0) |> 
  dplyr::filter(dif == min(dif)) |> 
  dplyr::pull(boniteit)

# Derive age from given diameter at breast height (dbh)
dbh <- 48

scots_pine |> 
  dplyr::filter(boniteit == {boni}) |> 
  dplyr::mutate(dif = {dbh} - d) |> 
  dplyr::filter(abs(dif) == min(abs(dif), na.rm = TRUE)) |> 
  dplyr::pull(t)
