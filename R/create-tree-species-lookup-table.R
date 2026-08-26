# Create tree species lookup table
# Author: Stefan Vriend
# Created: 2026-08-26

# Setup
# Load packages
library(RODBC)
library(dplyr)
library(tidyr)

# Load NBI data
nbi <- RODBC::odbcConnectAccess2007("NBI/NBI7_openbare database_definitief_UpdateJuli2024.accdb")

mfv_boomsoort <- RODBC::sqlFetch(nbi, "lookup_MFV_Boomsoort")
nbi6_boomsoort <- RODBC::sqlFetch(nbi, "lookup_NBI6_Boomsoort")
nbi7_boomsoort <- RODBC::sqlFetch(nbi, "lookup_NBI7_Boomsoort")

RODBC::odbcCloseAll()

# Wrangle tables: rename columns to English names and lowercase
mfv_tree_species <- mfv_boomsoort |>
  dplyr::rename("species_numeric_code" = "CODE", "species_name" = "OMSCHRIJVING", "species_code" = "NBI", "scientific_name" = "Wetenschapppelijke naam") |>
  dplyr::mutate("inventory" = "mfv") |>
  tidyr::drop_na("scientific_name")

nbi6_tree_species <- nbi6_boomsoort |>
  dplyr::select("species_code" = "code", "species_name" = "boomsoort", "nbi_species_group" = "Boomsoortgroep", "scientific_name" = "Wetenschappelijke naam") |>
  dplyr::mutate("inventory" = "nbi6") |> 
  tidyr::drop_na("scientific_name")

nbi7_tree_species <- nbi7_boomsoort |>
  dplyr::select("species_code" = "code", "species_name" = "boomsoort", "nbi_species_group" = "Boomsoortgroep", "scientific_name" = "Wetenschappelijke naam") |>
  dplyr::mutate("inventory" = "nbi7") |> 
  tidyr::drop_na("scientific_name")

# Add species group classification for LANDIS-II
# Our list of species of interest:
# berk (birch), beuk (beech), inlandse eik (pedunculate oak), Amerikaanse eik (red oak), Japanese larch (Japanse lariks), douglasspar (douglas fir), fijnspar (Norway spruce), grove den (Scots pine)
# in addition two rest categories: broadleaf and coniferous

shrubs <- c("BW", "DK", "GW", "HU", "HZ", "KM", "KN", "LB", "LG", "MD", "OW", "SD", "ST", "SW", "TX", "VB", "VK", "VL")

other_broadleaf <- c("AB", "AC", "AV", "ED", "EO", "ES", "GE", "HB", "IE", "IL", "IP", "LI", "NE", "OE", "PK", "PL", "PO", "RP", "SA", "TK", "UL", "WI", "ZE", "ZK", "ZP")

other_coniferous <- c("AA", "AG", "CD", "CH", "DO", "JE", "NO", "OA", "OD", "ON", "OS", "PC", "RD", "SO", "SS", "TH", "TS", "WD", "ZD")

tree_species_lookup <- dplyr::bind_rows(mfv_tree_species, nbi6_tree_species, nbi7_tree_species) |>
  dplyr::mutate("species_group" = dplyr::case_when(
    .data$species_code == "BE" ~ "betula",
    .data$species_code == "BU" ~ "fagus sylvatica",
    .data$species_code == "EI" ~ "quercus",
    .data$species_code == "AE" ~ "quercus rubra",
    .data$species_code %in% c("JL", "EL") ~ "larix kaempferi", # European larch included in Japanese larch
    .data$species_code == "DG" ~ "pseudotsuga menziesii",
    .data$species_code == "FS" ~ "picea abies",
    .data$species_code == "GD" ~ "pinus sylvestris",
    .data$species_code %in% {{shrubs}} ~ "shrub",
    .data$species_code %in% {{other_broadleaf}} ~ "broadleaf",
    .data$species_code %in% {{other_coniferous}} ~ "coniferous"
  ))

readr::write_csv(tree_species_lookup, file = "nbi-tree-species-lookup.csv")