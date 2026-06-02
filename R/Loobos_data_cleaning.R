#Download data from online source
url <- "https://zenodo.org/records/15721310/files/NL-Loo_treeinventory_100trees_1996-2025.xlsx?download=1"
download.file(url,"~/Masters/Internship-NIOO_KNAW/resources/data/vlm-git/vl-veluweflm/Loobos1996.xlsx", mode = "wb") #Excel workbook with multiple sheets

#Read excel workbook and list available sheets
sheets <- readxl::excel_sheets("Loobos1996.xlsx")
sheets

# Read in specific sheet
trees <- readxl::read_excel("Loobos1996.xlsx", sheet = "1996")
head(trees)
str(trees)
view(trees)

#Separate data into multiple files
tree_data1 <- readxl::read_excel(
    "Loobos1996.xlsx",
    sheet = "1996",
    range = "A13:Y154")

#start Cleaning the data
tree_data_cleaning1 <- tree_data1 %>%
dplyr::slice(-c(1, 2)) %>%
dplyr::filter(!if_all(4:25, is.na)) %>%
dplyr::select(where(~!all(is.na(.)))) %>%
dplyr::mutate(across(everything(), ~ na_if(., "Invalid Number"))) %>%
dplyr::mutate(across(where(is.numeric), ~round(., 2))) %>%
dplyr::rename(
    oldID = "Tree ID", newNum = "Tree...2", oldNum = "Tree...3", circum_cm = "Circum", sapArea_cm2 = "SapArea", crownBase_m = "Crown base",
    crownLength_m = "Crown lenght", crownProj_90_m = "Kroonprojectie (afstand vanaf stam (m))", crownProj_180_m = "...9", crownProj_270_m = "...10",
    crownProj_360_m = "...11", dbh_cm = "dbh", height_m = "height", volume_sb_m3tree = "Volume", volume_ba_m2tree = "...18", agBiomass_kgtree = "biomass...19", bgBiomass_kgtree = "biomass...20",
    branchBiomass_kgtree = "biomass...21", leafBiomass_kgtree = "biomass...22", rootBiomass_kgtree = "biomass...23", sbBiomass_kgtree = "biomass...24", totalABG_kgtree = "total ABG"
    )

#Convert to csv
readr::write_csv(
    tree_data_cleaning1, file = "Loobos1996_clean.csv"
)

#Read and work with the csv
loobos1 <- readr::read_csv("Loobos1996_clean.csv")
loobos1$newNum <- as.character(loobos1$newNum)
loobos1$oldNum <- as.character(loobos1$oldNum)


#Repeat above process for second half of the sheet
#------------------------------------------------
#Separate data into multiple files
tree_data2 <- readxl::read_excel(
    "Loobos1996.xlsx",
    sheet = "1996",
    range = "A158:Y282")

#start Cleaning the data
tree_data_cleaning2 <- tree_data2 %>%
dplyr::slice(-c(1, 2)) %>%
dplyr::filter(!if_all(4:25, is.na)) %>%
dplyr::select(where(~!all(is.na(.)))) %>%
dplyr::mutate(across(where(is.numeric), ~round(., 2)))%>%
dplyr::rename(
    newNum = "...2", oldNum = "...3", circum_cm = "...4", sapArea_cm2 = "...5", crownBase_m = "...6",
    crownLength_m = "...7", crownProj_90_m = "...8", crownProj_180_m = "...9", crownProj_270_m = "...10",
    crownProj_360_m = "...11", dbh_cm = "...15", height_m = "...16", volume_sb_m3tree = "...17", volume_ba_m2tree = "...18", agBiomass_kgtree = "...19", 
    bgBiomass_kgtree = "...20",branchBiomass_kgtree = "...21", leafBiomass_kgtree = "...22", 
    rootBiomass_kgtree = "...23", sbBiomass_kgtree = "...24", totalABG_kgtree = "...25"
    )

# Convert to csv
readr::write_csv(
    tree_data_cleaning2, file = "Loobos1996_clean2.csv"
)

#Read and work with the csv
loobos2 <- readr::read_csv("Loobos1996_clean2.csv")
loobos2$newNum <- as.character(loobos2$newNum)
loobos2$oldNum <- as.character(loobos2$oldNum)
loobos2


#Combine the two datafiles
#-------------------------
loobos_1996 <- dplyr::bind_rows(loobos1, loobos2)

readr::write_csv(loobos_1996, file = "Loobos1996.csv")

#Repeat for additional information
#----------------------------------

#Separate data into multiple files
tree_data_classes <- readxl::read_excel(
    "Loobos1996.xlsx",
    sheet = "1996",
    range = "C290:R306")
tree_data_classes

#start Cleaning the data
tree_data_classes_clean <- tree_data_classes %>%
dplyr::slice(-c(1)) %>%
dplyr::filter(!if_all(1:16, is.na)) %>%
dplyr::select(where(~!all(is.na(.)))) %>%
dplyr::rename(diameter_class = "diameter", mid_diam_class = "Gemiddelde", n_trees = "...3", n_perc = "...4", meanBA_cm2tree = "Mean BA",
BA_cm2ha = "Basal Area", meanSA_cm2tree = "Mean SA", sapwood_m2ha = "Sapwood", meanVol_dm3tree = "Mean volume", volume_m3ha = "volume", count_treeha = "number")



# BAI data cleaning #


BAI <- read.csv2("BAI_NL_Loo.csv")

BAI <- BAI |>
  dplyr::rename(
    year = anno,
    tree = Tree
  ) |>
  stats::na.omit(BAI)

BAI <- BAI |>
  dplyr::mutate(
    BAI_cm = BAI / 100
  )

unique(BAI$tree)

# 2023 inventory cleaning
trees <- readxl::read_excel("Loobos_TreeInventory_2023.xlsx", sheet = "Tree locations")

trees <- trees |>
  dplyr::group_by(TREE_ID) |>
  dplyr::select(TREE_ID, TREE_DBH, TREE_SPP) |>
  dplyr::rename(
    tree = TREE_ID,
    dbh_cm = TREE_DBH
  ) |>
  dplyr::mutate(
    genus = stringr::word(TREE_SPP, 1),
    species = stringr::word(TREE_SPP, 2),
    SpeciesCode = paste0(tolower(genus), tolower(species))
  )

tree_ids <- c(173, 175, 188, 219, 238, 301, 307, 346, 406, 409, 451,
              505, 514, 678, 705, 708, 921, 935)

trees <- trees |>
  dplyr::filter(
    tree %in% tree_ids
  )

trees <- trees |>
  dplyr::select(tree, dbh_cm, SpeciesCode)
