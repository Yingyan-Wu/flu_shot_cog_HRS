# LHMS explore

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven", "labelled",
       "lubridate")

#---- function ----
source(here::here("scripts", "functions.R"))
# Directory to Box
# path_to_box <- "C:/Users/yingyan_wu/Box"
# path_to_box <- "C:/Users/Yingyan Wu/Box"
path_to_box <- "/Users/yingyanwu/Library/CloudStorage/Box-Box/"

#---- Read the data ----
lhms_2015_2017 <- read_sas(paste0(path_to_box, "HRS_MRG/HRSRawData/",
                        "Life History/LHMS_Aggregate_2015_2017/lhms1517a_r.sas7bdat")) %>%
  unite("HHIDPN", c("hhid", "pn"), sep = "", remove = TRUE) %>%
  mutate_at("HHIDPN", as.numeric)

lhms_2019 <- read_sas(paste0(path_to_box, "HRS_MRG/HRSRawData/",
                             "Life History/LHMS19SPR/lhms19spr_r.sas7bdat")) %>%
  unite("HHIDPN", c("HHID", "PN"), sep = "", remove = TRUE) %>%
  mutate_at("HHIDPN", as.numeric)

# Sanity check participants
table(unique(lhms_2019$HHIDPN) %in% unique(lhms_2015_2017$HHIDPN), useNA = "ifany")

lhms_selected <- lhms_2015_2017 %>% 
  select(HHIDPN, contains(paste0("LH", 24:31)), -contains("LH30A")) %>%
  rbind(lhms_2019 %>% select(-VERSION) %>%
  rename_all(~str_replace(., "_19", "")) %>%
  select(HHIDPN, contains(paste0("LH", 24:31)), -LH29, -contains("LH30A")))

#---- Check the missingness of the variables of interest ----
summary(lhms_selected)

