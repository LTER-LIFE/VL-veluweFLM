#this script will return all species trait files that biomass succession requires. Currently only made for biomass succession.

library(dplyr)
library(stringr)
library(tidyr)
library(readr)

#block for landis extension in flowchart for which extensions to include
extension_vector <- c("Biomass_succession", "Output_biomass") #etc

# Location data sources
input_path <- file.path("/home/jovyan/Cloud Storage/naa-vre-user-data/")

#retrieve species trait data
#load and prepare TRY database for biomass succession
TRY_biomass <- read.delim(paste0(input_path, "try-data.txt"), fileEncoding = "latin1", dec = ".", quote = "", stringsAsFactors = FALSE)
TRY_biomass <- TRY_biomass %>%
  select(-`X`) %>%
  mutate(
    SpeciesCode = str_split(AccSpeciesName, " ", simplify = TRUE) %>%
      apply(1, function(x) paste0(str_sub(x[1], 1, 4), str_sub(x[2], 1, 4))),
    
    SpeciesCode = case_when(
      SpeciesCode %in% c("Querrobu", "Querpetr") ~ "Querspec",
      SpeciesCode %in% c("Betupend", "Betupube") ~ "Betuspec",
      TRUE ~ SpeciesCode
    )
  )

#prepare data for Landis-II

#function to easily min, mean or max a value for a tree trait in TRY    
summarize_orig_value_by_species <- function(data, stat = c("min", "max", "mean"), name = "result") {
  stat <- match.arg(stat)  # Validate input
  name_sym <- sym(name)    # Convert name to symbol for tidy evaluation
  
  data %>%
    group_by(SpeciesCode) %>%
    summarize(
      !!name_sym := case_when(
        stat == "min"  ~ min(OrigValueStr, na.rm = TRUE),
        stat == "max"  ~ max(OrigValueStr, na.rm = TRUE),
        stat == "mean" ~ mean(OrigValueStr, na.rm = TRUE)
      ),
      .groups = "drop"
    )
}



# ---- core species data ----
#is needed for every succession extension
    
##prepare TRY database for Landis parameter sexual maturity (traitID 155)
TRY_155 <- TRY_biomass %>%
  filter(TraitID == 155) %>%
  mutate(
    OrigValueStr = na_if(OrigValueStr, ""),  # Remove empty strings
    is_date = str_detect(OrigValueStr, "\\d{4}-\\d{2}-\\d{2}"),  # Flag dates
    number_list = str_extract_all(OrigValueStr, "\\d+\\.?\\d*"),  # Extract numbers
    OrigValueStr = sapply(number_list, function(x) {  # Overwrite with cleaned value
      x_num <- suppressWarnings(as.numeric(x))
      if (length(x_num) == 0 || all(is.na(x_num))) return(NA_real_)
      else return(mean(x_num, na.rm = TRUE))
    })
  ) %>%
  filter(!is.na(OrigValueStr), !is_date) #keep only rows with valid values

#makes a column with mean age of flowering of trees, TRY misses Pseumenz and Larikaem
df_flower_age <- summarize_orig_value_by_species(TRY_155, "mean", "Sexual Maturity")
Pseumenz_flowerage <- 12 #future possible pipeline, for now not possible
Larikaem_flowerage <- 10 #future possible pipeline, for now not possible

# Add missing species manually
df_flower_age_missing <- tibble(
  SpeciesCode = c("Pseumenz", "Larikaem"),
  "Sexual Maturity" = c(Pseumenz_flowerage, Larikaem_flowerage)
)

# Combine both
df_flower_age <- bind_rows(df_flower_age, df_flower_age_missing)

##prepare TRY database for Landis parameter plant longevity (traitID 59)
TRY_59 <- TRY_biomass %>%
  filter(TraitID == 59) %>%
  mutate(
    OrigValueStr = na_if(OrigValueStr, ""),  # Remove empty strings
    is_date = str_detect(OrigValueStr, "\\d{4}-\\d{2}-\\d{2}"),  # Flag dates
    OrigValueStr = str_extract(OrigValueStr, "\\d+\\.?\\d*"),  # Extract single number
    OrigValueStr = as.numeric(OrigValueStr)  # Convert to numeric
  ) %>%
  filter(!is.na(OrigValueStr), !is_date)  # Keep only rows with valid values

#makes a column with max age of a tree
df_plant_age <- summarize_orig_value_by_species(TRY_59, "max", "Longevity")

##prepare TRY database for Landis parameter dispersal distance (traitID 193)
TRY_193 <- TRY_biomass %>%
  filter(TraitID == 193) %>%
  mutate(OrigValueStr = na_if(OrigValueStr, ""),  # Remove empty strings
         OrigValueStr = as.numeric(OrigValueStr)  # Convert to numeric
  ) %>%
  filter(!is.na(OrigValueStr))

#makes a column with effective dispersal distance, TRY misses Betuspec, Fagusylv, Piceabie, Querspec and Larikaem
df_eff_seed_dist <- TRY_193 %>% 
  filter(!str_starts(OriglName, "Max")) %>% # Exclude names starting with "Max"
  summarize_orig_value_by_species(., "mean", "Seed Dispersal Dist Effective")
Betuspec_effseeddist <- 50 #future possible pipeline, for now not possible
Fagusylv_effseeddist <- 30 #future possible pipeline, for now not possible
Piceabie_effseeddist <- 45 #future possible pipeline, for now not possible
Querspec_effseeddist <- 100 #future possible pipeline, for now not possible
Larikaem_effseeddist <- 50 #future possible pipeline, for now not possible

# Add missing species manually
df_eff_seed_dist_missing <- tibble(
  SpeciesCode = c("Betuspec", "Fagusylv", "Piceabie", "Querspec", "Larikaem"),
  "Seed Dispersal Dist Effective" = c(Betuspec_effseeddist, Fagusylv_effseeddist, Piceabie_effseeddist, Querspec_effseeddist, Larikaem_effseeddist)
) 

#complete effective_seed_distance
df_eff_seed_dist <- bind_rows(df_eff_seed_dist, df_eff_seed_dist_missing)


#makes a column with maximum dispersal distance, TRY misses Querrubr, Betuspec, Fagusylv, Piceabie, Querspec and Larikaem
df_max_seed_dist <- TRY_193 %>% 
  filter(str_starts(OriglName, "Max")) %>% # Include names starting with "Max"
  filter(!str_starts(OriglName, "Eff")) %>% #Exclude names starting with "Effective"
  summarize_orig_value_by_species(.,"max", "Seed Dispersal Dist Maximum")
Betuspec_maxseeddist <- 1000 #future possible pipeline, for now not possible
Fagusylv_maxseeddist <- 70 #future possible pipeline, for now not possible
Piceabie_maxseeddist <- 200 #future possible pipeline, for now not possible
Querspec_maxseeddist <- 1500 #future possible pipeline, for now not possible
Larikaem_maxseeddist <- 1000 #future possible pipeline, for now not possible
Querrubr_maxseeddist <- 1500 #future possible pipeline, for now not possible, this still needs a source!!!

# Add missing species manually
df_max_seed_dist_missing <- tibble(
  SpeciesCode = c("Betuspec_maxseeddist", "Fagusylv_maxseeddist", "Piceabie_maxseeddist", "Querspec_maxseeddist", "Larikaem_maxseeddist", "Querrubr_maxseeddist"),
  "Seed Dispersal Dist Maximum" = c(Betuspec_maxseeddist, Fagusylv_maxseeddist, Piceabie_maxseeddist, Querspec_maxseeddist, Larikaem_maxseeddist, Querrubr_maxseeddist)
)

#complete effective_seed_distance
df_max_seed_dist <- bind_rows(df_max_seed_dist, df_max_seed_dist_missing)  %>%
  mutate(SpeciesCode = str_sub(SpeciesCode, 1, 8))

##prepare TRY database for resprouting probability (Traitid 819)
TRY_819 <- TRY_biomass %>%
  filter(TraitID == 819) %>%
  filter(!str_detect(OrigUnitStr, "%")) %>%
  mutate(
    OrigValueStr = na_if(OrigValueStr, ""),  # Remove empty strings
    OrigValueStr = case_when(
      OrigValueStr %in% c("No", "no", "1") ~ "0.0", #0.0 for no chance for resprouting
      OrigValueStr %in% c("Yes", "yes", "0", "moderate") ~ "0.8", #0.8 for all species that are able to resprout
      TRUE ~ OrigValueStr
    ),
    OrigValueStr = as.numeric(OrigValueStr)  # Convert to numeric
  ) %>%
  filter(!is.na(OrigValueStr))  # Keep only rows with valid values

#makes a column with resprouting probability of trees
df_resprout_prob <- summarize_orig_value_by_species(TRY_819, "min", "Vegetative Reprod Prob") #min because Larix kaempferi had an able to resprout data point, however that paper did not mention larix kaempferi
#also pinus sylvestris had an entry that had able to resprout, however, in that database the entry was not able to resprout
#also picea abies had an entry that had able to resprout, however, there was no source in the database, and other sources specified that picea abies was definitely not able to resprout

##makes a column with minimum resprouting age, TRY misses all species
Querrubr_minresprout <- 20 #future possible pipeline, for now not possible
Betuspec_minresprout <- 0 #future possible pipeline, for now not possible
Fagusylv_minresprout <- 20 #future possible pipeline, for now not possible
Pseumenz_minresprout <- 0 #future possible pipeline, for now not possible, species that do not resprout need a 0
Piceabie_minresprout <- 0 #future possible pipeline, for now not possible, species that do not resprout need a 0
Pinusylv_minresprout <- 0 #future possible pipeline, for now not possible, species that do not resprout need a 0
Querspec_minresprout <- 20 #future possible pipeline, for now not possible
Larikaem_minresprout <- 0 #future possible pipeline, for now not possible, species that do not resprout need a 0

#complete df minimum resprouting age
df_min_resprout_age <- tibble(
  SpeciesCode = c("Querrubr", "Betuspec", "Fagusylv", "Pseumenz", "Piceabie", "Pinusylv", "Querspec", "Larikaem"),
  "Sprout Age Min" = c(Querrubr_minresprout, Betuspec_minresprout, Fagusylv_minresprout, Pseumenz_minresprout, Piceabie_minresprout, Pinusylv_minresprout, Querspec_minresprout, Larikaem_minresprout)
)

##makes a column with maximum resprouting age, TRY misses all species
Querrubr_maxresprout <- 200 #future possible pipeline, for now not possible
Betuspec_maxresprout <- 100 #future possible pipeline, for now not possible
Fagusylv_maxresprout <- 150 #future possible pipeline, for now not possible
Pseumenz_maxresprout <- 0 #future possible pipeline, for now not possible, species that do not resprout need a 0
Piceabie_maxresprout <- 0 #future possible pipeline, for now not possible, species that do not resprout need a 0
Pinusylv_maxresprout <- 0 #future possible pipeline, for now not possible, species that do not resprout need a 0
Querspec_maxresprout <- 150 #future possible pipeline, for now not possible
Larikaem_maxresprout <- 0 #future possible pipeline, for now not possible, species that do not resprout need a 0

#complete df maximum resprouting age
df_max_resprout_age <- tibble(
  SpeciesCode = c("Querrubr", "Betuspec", "Fagusylv", "Pseumenz", "Piceabie", "Pinusylv", "Querspec", "Larikaem"),
  "Sprout Age Max" = c(Querrubr_maxresprout, Betuspec_maxresprout, Fagusylv_maxresprout, Pseumenz_maxresprout, Piceabie_maxresprout, Pinusylv_maxresprout, Querspec_maxresprout, Larikaem_maxresprout)
)


##prepare TRY database for method of post fire regenaration (Traitid 318 and 819)
TRY_318.819 <- TRY_biomass %>%
  filter(TraitID %in% c(318, 819)) %>%
  mutate(
    OrigValueStr = case_when( #resprout and serotiny are in different TRY traitID, serotiny takes priority over resprout over none
      TraitID == 318 & OriglName == "SeedlEmerg" & OrigValueStr %in% c("low", "yes") ~ "serotiny",
      TraitID == 819 & OrigValueStr %in% c("No", "no", "1") ~ "none", #TRY paper uses 1 for no
      TraitID == 819 & OrigValueStr %in% c("Yes", "yes", "0", "moderate", "70") ~ "resprout", #TRY paper uses 0 for yes
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(OrigValueStr)) %>%
  group_by(SpeciesCode) %>%
  summarize(
    OrigValueStr = case_when(
      any(OrigValueStr == "serotiny") ~ "serotiny",
      any(OrigValueStr == "none") ~ "none",
      any(OrigValueStr == "resprout") ~ "resprout"),
    .groups = "drop"
  )

df_postfire_regen <- TRY_318.819 #not with function, since value is not numerical


#make final input file for core species data
core_species_data <- df_plant_age %>%
  full_join(df_flower_age, by = "SpeciesCode") %>%
  full_join(df_eff_seed_dist, by = "SpeciesCode") %>%
  full_join(df_max_seed_dist, by = "SpeciesCode") %>%
  full_join(df_resprout_prob, by = "SpeciesCode") %>%
  full_join(df_min_resprout_age, by = "SpeciesCode") %>%
  full_join(df_max_resprout_age, by = "SpeciesCode") %>%
  full_join(df_postfire_regen, by = "SpeciesCode")

core_species_data <- apply(core_species_data, 1, function(row) paste(row, collapse = "\t"))

header <- c(
  "LandisData Species"
)

#part that should make the output to the next NaaVRE block
output_lines <- c(header, core_species_data)

writeLines(output_lines, "Core_species_data.txt")

# ---- SppEcoregionData for biomass succession ----

#make an if loop, eventually this file should contain all the code for all species trait data that is required per 
#succession extension. So the following code will only run in case biomass succession is the succession extension that is chosen
if ("Biomass_succession" %in% extension_vector) {
  
#also read in establishment values from Probos, Probos misses Querrubr
df_prob_establish <- read_csv(paste0(input_path, "establishmentTreeVeluwe.csv")) %>%
  select(-1) %>%              # Remove first column if it's row names
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "SpeciesCode", values_to = "ProbEstablish") %>%
  mutate(
    Year = 0,
    EcoregionName = 101,
    ProbEstablish = pmin(ProbEstablish / 8, 1)
  )

#add missing species manually
Querrubr_prob_establish <- 1

#make dataframe for ProbEstablishment missing
df_prob_establish_missing <- tibble(
  Year = 0,
  EcoregionName = 101,
  SpeciesCode = c("Querrubr"),
  ProbEstablish = c(Querrubr_prob_establish)
)

#make final dataframe
df_prob_establish <- bind_rows(df_prob_establish, df_prob_establish_missing)

#read in vitality scores from NBI 7
mortality <- read.csv(paste0(input_path, "probMortality_biomass.txt"))

#make dataframe for Probmortality
df_prob_mortality <- mortality %>%
    dplyr::mutate(EcoregionName = 101,
                  Year = 0)

#prepare ANPPmax, comes from training algorythm, but for now default values
Pseumenz_anpp_max <- 7.12 #ANPPmax_pseumenz 
Larikaem_anpp_max <- 7.12 #ANPPmax_larikaem 
Piceabie_anpp_max <- 7.12 #ANPPmax_piceabie 
Pinusylv_anpp_max <- 7.12 #ANPPmax_pinusylv 
Fagusylv_anpp_max <- 11.3 #ANPPmax_fagusylv 
Querspec_anpp_max <- 11.3 #ANPPmax_querspec 
Querrubr_anpp_max <- 11.3 #ANPPmax_querrubr 
Betuspec_anpp_max <- 11.3 #ANPPmax_betuspec 

#make dataframe for ANPPmax
df_anpp_max <- tibble(
  Year = 0,
  EcoregionName = 101,
  SpeciesCode = c("Larikaem", "Pseumenz", "Piceabie", "Pinusylv", "Fagusylv", "Querspec", "Querrubr", "Betuspec"),
  ANPPmax = c(Larikaem_anpp_max, Pseumenz_anpp_max, Piceabie_anpp_max, Pinusylv_anpp_max, 
              Fagusylv_anpp_max, Querspec_anpp_max, Querrubr_anpp_max, Betuspec_anpp_max)
)

#prepare BiomassMax, comes from training algorythm, but for now default values
Pseumenz_biomass_max <- 26000 #biomassMAX_pseumenz 
Larikaem_biomass_max <- 26000 #biomassMAX_larikaem
Piceabie_biomass_max <- 26000 #biomassMAX_piceabie
Pinusylv_biomass_max <- 26000 #biomassMAX_pinusylv
Fagusylv_biomass_max <- 26000 #biomassMAX_fagusylv
Querspec_biomass_max <- 26000 #biomassMAX_querspec
Querrubr_biomass_max <- 26000 #biomassMAX_querrubr
Betuspec_biomass_max <- 26000 #biomassMAX_betuspec

#make dataframe for BiomassMax
df_biomass_max <- tibble(
  Year = 0,
  EcoregionName = 101,
  SpeciesCode = c("Larikaem", "Pseumenz", "Piceabie", "Pinusylv", "Fagusylv", "Querspec", "Querrubr", "Betuspec"),
  BiomassMax = c(Larikaem_biomass_max, Pseumenz_biomass_max, Piceabie_biomass_max, Pinusylv_biomass_max, 
              Fagusylv_biomass_max, Querspec_biomass_max, Querrubr_biomass_max, Betuspec_biomass_max)
) 

#make final input file for Species ecoregion Data for biomass
SppEcoregionData.biomass <- df_prob_establish %>%
  full_join(df_prob_mortality, by = c("SpeciesCode", "Year", "EcoregionName")) %>%
  full_join(df_anpp_max, by = c("SpeciesCode", "Year", "EcoregionName")) %>%
  full_join(df_biomass_max, by = c("SpeciesCode", "Year", "EcoregionName"))

write_csv(SppEcoregionData.biomass, "SppEcoregionData.csv")
}


# ---- SpeciesDataBiomass for biomass succession ----
if ("Biomass_succession" %in% extension_vector) {
#prepare TRY database for Landis parameter leaf longevity (traitID 12)
TRY_12 <- TRY_biomass %>%
  filter(TraitID == 12) %>%
  mutate(
    # Remove non-numeric characters from OrigValueStr
    OrigValueStr = str_extract(OrigValueStr, "\\d+\\.?\\d*"),
    OrigValueStr = as.numeric(OrigValueStr),
    
    # Normalize OrigUnitStr to lowercase for easier case_when matching
    OrigUnitStr = str_to_lower(OrigUnitStr),
    
    # Apply unit conversion
    OrigValueStr = case_when(
      OrigUnitStr == "days" ~ OrigValueStr / 365,
      OrigUnitStr %in% c("month", "months") ~ OrigValueStr / 12,
      OrigUnitStr %in% c("year", "yr-1", "year-1") ~ OrigValueStr,
      TRUE ~ NA_real_  # Handle unexpected units
    ),
    
    # Enforce minimum value of 1, because deciduous trees in landis should have a value of 1
    OrigValueStr = if_else(OrigValueStr < 1, 1, OrigValueStr),
    
    # Standardize unit label
    OrigUnitStr = "year"
  ) %>%
  filter(!is.na(OrigValueStr))  

#makes a column with mean longevity of leafs, TRY misses Larikaem
df_leaf_longevity <- summarize_orig_value_by_species(TRY_12, "mean", "LeafLongevity")
Larikaem_leaf_longevity <- 1 #future possible pipeline, for now not possible

# Add missing species manually
df_leaf_longevity_missing <- tibble(
  SpeciesCode = c("Larikaem"),
  LeafLongevity = c(Larikaem_leaf_longevity)
)

# Combine both
df_leaf_longevity <- bind_rows(df_leaf_longevity, df_leaf_longevity_missing)  %>%
  mutate(SpeciesCode = str_sub(SpeciesCode, 1, 8))

#prepare TRY database for Landis parameter wood decay rate (traitID 1158)
TRY_1158 <- TRY_biomass %>%
  filter(TraitID == 1158) %>%
  filter(OriglName == "single exponential model: Yt=Yoe-kt") %>%
  filter(OrigValueStr >= 0) %>%
  mutate(OrigValueStr = as.numeric(OrigValueStr))

#makes a column with mean of wood decay rate, k in Yt=Yoe-kt, TRY misses Larikaem
df_wood_decayrate <- summarize_orig_value_by_species(TRY_1158, "mean", "WoodDecayRate")
Larikaem_wood_decayrate <- 0.021 #future possible pipeline, for now not possible

# Add missing species manually
df_wood_decayrate_missing <- tibble(
  SpeciesCode = c("Larikaem"),
  WoodDecayRate = c(Larikaem_wood_decayrate)
)

# Combine both
df_wood_decayrate <- bind_rows(df_wood_decayrate, df_wood_decayrate_missing)  %>%
  mutate(SpeciesCode = str_sub(SpeciesCode, 1, 8))



#prepare mortalitycurve
Pseumenz_mortality_curve <- 5 #future might be different because of different fit, for now 5, because of high maximum age
Larikaem_mortality_curve <- 5 #future might be different because of different fit, for now 5, because of high maximum age
Piceabie_mortality_curve <- 5 #future might be different because of different fit, for now 5, because of high maximum age
Pinusylv_mortality_curve <- 5 #future might be different because of different fit, for now 5, because of high maximum age
Fagusylv_mortality_curve <- 5 #future might be different because of different fit, for now 5, because of high maximum age
Querspec_mortality_curve <- 5 #future might be different because of different fit, for now 5, because of high maximum age
Querrubr_mortality_curve <- 5 #future might be different because of different fit, for now 5, because of high maximum age
Betuspec_mortality_curve <- 5 #future might be different because of different fit, for now 5, because of high maximum age

#make dataframe for mortalitycurve
df_mortality_curve <- tibble(
  SpeciesCode = c("Larikaem", "Pseumenz", "Piceabie", "Pinusylv", "Fagusylv", "Querspec", "Querrubr", "Betuspec"),
  MortalityCurve = c(Larikaem_mortality_curve, Pseumenz_mortality_curve, Piceabie_mortality_curve, Pinusylv_mortality_curve, 
                     Fagusylv_mortality_curve, Querspec_mortality_curve, Querrubr_mortality_curve, Betuspec_mortality_curve)
)

#prepare growthcurve, should come from training algorythm, but for now default values
Pseumenz_growth_curve <- 0.9 #growthcurve_pseumenz
Larikaem_growth_curve <- 0.9 #growthcurve_larikaem
Piceabie_growth_curve <- 0.9 #growthcurve_piceabie
Pinusylv_growth_curve <- 0.9 #growthcurve_pinusylv
Fagusylv_growth_curve <- 0.9 #growthcurve_fagusylv
Querspec_growth_curve <- 0.9 #growthcurve_querspec
Querrubr_growth_curve <- 0.9 #growthcurve_querrubr
Betuspec_growth_curve <- 0.9 #growthcurve_betuspec

#make dataframe for growthcurve
df_growth_curve <- tibble(
  SpeciesCode = c("Larikaem", "Pseumenz", "Piceabie", "Pinusylv", "Fagusylv", "Querspec", "Querrubr", "Betuspec"),
  GrowthCurve = c(Larikaem_growth_curve, Pseumenz_growth_curve, Piceabie_growth_curve, Pinusylv_growth_curve, 
                  Fagusylv_growth_curve, Querspec_growth_curve, Querrubr_growth_curve, Betuspec_growth_curve)
)

#prepare TRY database for Landis parameter leaf lignin content (traitID 87)
TRY_87 <- TRY_biomass %>%
  filter(TraitID == 87) %>%
  mutate(
    # Remove non-numeric characters from OrigValueStr
    OrigValueStr = as.numeric(OrigValueStr)
    
  ) %>%
  filter(!is.na(OrigValueStr))  

#makes a column with mean of leaf lignin content, TRY misses Larikaem, Pseumenz, Piceabie, Fagusylv, Querspec, Querrubr
df_leaf_lignin <- summarize_orig_value_by_species(TRY_87, "mean", "LeafLignin")
Larikaem_leaf_lignin <- 0.21 #future possible pipeline, for now not possible
Pseumenz_leaf_lignin <- 0.19 #future possible pipeline, for now not possible
Piceabie_leaf_lignin <- 0.18 #future possible pipeline, for now not possible
Fagusylv_leaf_lignin <- 0.15 #future possible pipeline, for now not possible
Querspec_leaf_lignin <- 0.167 #future possible pipeline, for now not possible
Querrubr_leaf_lignin <- 0.145 #future possible pipeline, for now not possible


# Add missing species manually
df_leaf_lignin_missing <- tibble(
  SpeciesCode = c("Larikaem", "Pseumenz", "Piceabie", "Fagusylv", "Querspec", "Querrubr"),
  LeafLignin = c(Larikaem_leaf_lignin, Pseumenz_leaf_lignin, Piceabie_leaf_lignin, Fagusylv_leaf_lignin, Querspec_leaf_lignin, Querrubr_leaf_lignin))

# Combine both
df_leaf_lignin <- bind_rows(df_leaf_lignin, df_leaf_lignin_missing)  %>%
  mutate(SpeciesCode = str_sub(SpeciesCode, 1, 8))

#prepare TRY database for Landis parameter shade tolerance (traitID 603)
TRY_603 <- TRY_biomass %>%
  filter(TraitID == 603) %>%
  # Remove rows with "%" in OrigUnitStr
  filter(!str_detect(OrigUnitStr, "%")) %>%
  
  # Normalize OrigValueStr to lowercase for easier case_when usage
  mutate(OrigValueStr = str_to_lower(OrigValueStr)) %>%
  
  # Map textual values to numeric scale
  mutate(
    OrigValueStr = case_when(
      OrigValueStr %in% c("intolerant") ~ 1,
      OrigValueStr %in% c("indifferent - intolerant") ~ 2,
      OrigValueStr %in% c("intermediate") ~ 3,
      OrigValueStr %in% c("indifferent") ~ 4,
      OrigValueStr %in% c("tolerant") ~ 5,
      OrigValueStr %in% as.character(1:9) ~ case_when(
        OrigValueStr == "1" ~ 5,
        OrigValueStr == "2" ~ 4.5,
        OrigValueStr == "3" ~ 4,
        OrigValueStr == "4" ~ 3.5,
        OrigValueStr == "5" ~ 3,
        OrigValueStr == "6" ~ 2.5,
        OrigValueStr == "7" ~ 2,
        OrigValueStr == "8" ~ 1.5,
        OrigValueStr == "9" ~ 1
      ),
      TRUE ~ NA_real_  # Remove other non-numeric or unknown values
    )
  ) %>%
  
  # Remove rows with NA values after transformation
  filter(!is.na(OrigValueStr)) %>%
  
  # Standardize OrigUnitStr label
  mutate(OrigUnitStr = "tolerance_score")


#makes a column with mean of shade tolerance
df_shade_tolerance <- summarize_orig_value_by_species(TRY_603, "mean", "ShadeTolerance") %>%
  mutate(ShadeTolerance = round(ShadeTolerance)) #Landis requires integer input in this category

#prepare TRY database for Landis parameter fire tolerance (traitID 318)
TRY_318 <- TRY_biomass %>%
  filter(TraitID == 318) %>%
  filter(OriglName != "SeedlEmerg") %>%
  # Step 2: Normalize OrigUnitStr to lowercase
  mutate(OrigValueStr = str_to_lower(OrigValueStr)) %>%
  
  # Map textual values to numeric scale
  mutate(
    OrigValueStr = case_when(
      OrigValueStr %in% c("no", "none") ~ 1,
      OrigValueStr %in% c("low", "surface fire") ~ 2,
      OrigValueStr %in% c("yes") ~ 4,
      TRUE ~ NA_real_  # Remove other non-numeric or unknown values
    )
  ) %>%
  
  # Remove rows with NA values after transformation
  filter(!is.na(OrigValueStr)) %>%
  
  # Standardize OrigUnitStr label
  mutate(OrigUnitStr = "tolerance_score")


#makes a column with mean of fire tolerance, TRY misses Querspec and Fausylv
df_fire_tolerance <- summarize_orig_value_by_species(TRY_318, "mean", "FireTolerance") %>%
  mutate(FireTolerance= round(FireTolerance)) #Landis requires integer input in this category
Querspec_fire_tole <- 3 #future possible pipeline, for now not possible
Fagusylv_fire_tole <- 2 #future possible pipeline, for now not possible

# Add missing species manually
df_fire_tole_missing <- tibble(
  SpeciesCode = c("Querspec", "Fagusylv"),
  FireTolerance = c(Querspec_fire_tole, Fagusylv_fire_tole)
)

# Combine both
df_fire_tolerance <- bind_rows(df_fire_tolerance, df_fire_tole_missing)  %>%
  mutate(SpeciesCode = str_sub(SpeciesCode, 1, 8))


#make final input file for SpeciesData for biomass
SpeciesData.biomass <- df_leaf_longevity %>%
  full_join(df_wood_decayrate, by = "SpeciesCode") %>%
  full_join(df_mortality_curve, by = "SpeciesCode") %>%
  full_join(df_growth_curve, by = "SpeciesCode") %>%
  full_join(df_leaf_lignin, by = "SpeciesCode") %>%
  full_join(df_shade_tolerance, by = "SpeciesCode") %>%
  full_join(df_fire_tolerance, by = "SpeciesCode") 


write_csv(SpeciesData.biomass, "SpeciesData.csv")

}

