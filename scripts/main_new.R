library(tidyverse)
library(ggplot2)
library(ape)
library(aplot)
library(geoR)
library(ggtree)
library(brotli)
library(brms)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# convert Kinbank format to Mollica data format
# primitive concepts (upper case for female ego, lower case for female ego):
# S/s: son
# D/d: daughter
# F/f: father
# M/m: mother
# B/b: brother
# Z/z: sister

# concepts not in Mollica data format
# H/h: husband
# W/w: wife
# E/e: spouse
# C/c: child
# P/p: parent
# G/g: sibling
# R/r: grandparent (KinBank also used G for grandparent)
# A/a: ancestor

# suffixes (to be added after the concept, instead of before the concepts, as in KinBank):
# y: younger
# e: elder

# kinship concepts that we don't consider:
# born same day
# husband same group
# wifes same group
# concepts related to "not"
# agnatic / cognatic terms
# twins

kinship_map <- list(
  "C" = c("S", "D"),  # Child -> Son, Daughter
  "E" = c("H", "W"),  # Spouse -> Husband, Wife
  "P" = c("F", "M"), # Parent -> Father, Mother
  "G" = c("B", "Z"), # Sibling -> Brother, Sister
  "R" = c("FF", "FM", "MF", "MM"), # Grandparent -> ...
  # male ego
  "c" = c("s", "d"),  # Child -> Son, Daughter
  "e" = c("h", "w"),  # Spouse -> Husband, Wife
  "p" = c("f", "m"), # Parent -> Father, Mother
  "g" = c("b", "z"), # Sibling -> Brother, Sister
  "r" = c("ff", "fm", "MF", "MM") # Grandparent -> ...
)

# Function to expand kinship acronyms
expand_kinship <- function(acronym) {
  # Split acronym into characters
  chars <- str_split(acronym, "")[[1]]
  
  # Replace each character with its specific alternatives
  expanded_list <- map(chars, ~ kinship_map[[.x]] %||% .x)
  
  # Generate all possible combinations
  expanded_combinations <- expand.grid(expanded_list, stringsAsFactors = FALSE) %>%
    unite("expanded", everything(), sep = "") %>%
    pull(expanded)
  
  # Return as underscore-separated string
  str_c(expanded_combinations, collapse = "_")
}


forms <- read.csv('../kinbank/cldf/forms.csv') %>%
  group_by(Language_ID, Parameter_ID, Value, Form) %>%
  filter(row_number() == 1)
languages <- read.csv('../kinbank/cldf/languages.csv')
MDL <- read.csv('../data/MDL_Hyps.csv')



to_exclude <- c('BornSameDay', 'HusbandsSameGroup', 'WifesSameGroup')

parameters <- read.csv('../kinbank/cldf/parameters.csv') %>% 
  filter(!ID %in% to_exclude,
         !grepl('_.*_', ID),
         !grepl('not', ID),
         !grepl('twin', Name)) %>%
  rowwise() %>%
  mutate(ID_new = str_replace_all(ID, '([ye]+)([SDFMBZHWTCPAEG])', '\\2\\1'),
         ID_new2 = ifelse(grepl('^m.*[SDFMBZHWTCPAEG]', ID_new), 
                          str_sub(str_to_lower(ID_new), 2, -1), str_sub(ID_new, 2, -1)),
         ID_new3 = case_when(
           ID == 'mG' ~ 'r',
           ID == 'fG' ~ 'R',
           TRUE ~ ID_new2
         ),
         revised_ID = ID_new3
  ) %>%
  select(-starts_with('ID_new'))


df <- forms %>% 
  left_join(languages %>% select(-Name), by=c('Language_ID' = 'ID')) %>% 
  left_join(parameters, by=c('Parameter_ID' = 'ID')) %>% 
  filter(Glottocode != 'aust1271',
         !(Source == 'BarthW2017' & Glottocode == 'stan1293'),
         !(Source == 'el-khaissi_2020' & Glottocode == 'stan1323'),
         !(Source == 'Evans2016' & Glottocode == 'russ1263'),
         !(Source == 'BarthW2017' & Glottocode == 'stan1295')) %>% 
  select(-Source) %>%
  distinct()

extensions <- df %>% 
  filter(Form != 'TRUE') %>%
  group_by(Glottolog_Name,  Language_ID, Glottocode, Family, Latitude, Longitude, Form) %>% 
  distinct() %>%
  arrange(revised_ID) %>% 
  summarize(extension = list(revised_ID)) %>% 
  rowwise() %>%
  mutate(extension = paste(extension, collapse='_'),
  ) %>% ungroup()
