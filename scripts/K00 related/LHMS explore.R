# LHMS explore

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven", "labelled",
       "lubridate", "gtsummary")

#---- function ----
source(here::here("scripts", "functions.R"))
# Directory to Box
# path_to_box <- "C:/Users/yingyan_wu/Box"
# path_to_box <- "C:/Users/Yingyan Wu/Box"
path_to_box <- "/Users/yingyanwu/Library/CloudStorage/Box-Box"

#---- Life history ----
lhms_2015_2017 <- read_sas(paste0(path_to_box, "/HRS_MRG/HRSRawData",
                                  "/Life History/LHMS_Aggregate_2015_2017/lhms1517a_r.sas7bdat")) %>%
  unite("HHIDPN", c("hhid", "pn"), sep = "", remove = TRUE) %>%
  mutate_at("HHIDPN", as.numeric)

lhms_2019 <- read_sas(paste0(path_to_box, "/HRS_MRG/HRSRawData",
                             "/Life History/LHMS19SPR/lhms19spr_r.sas7bdat")) %>%
  unite("HHIDPN", c("HHID", "PN"), sep = "", remove = TRUE) %>%
  mutate_at("HHIDPN", as.numeric)

# Sanity check participants
table(unique(lhms_2019$HHIDPN) %in% unique(lhms_2015_2017$HHIDPN), useNA = "ifany")

lhms_selected <- lhms_2015_2017 %>% 
  select(HHIDPN, contains(paste0("LH", 24:31)), LH33, -contains("LH30A")) %>%
  rbind(lhms_2019 %>% select(-VERSION) %>%
          rename_all(~str_replace(., "_19", "")) %>%
          select(HHIDPN, contains(paste0("LH", 24:31)), LH33, -LH29, -contains("LH30A")))

##---- Cognitive activities ----
# LH29A Q29A. HIGH SCHOOL TAKE COLLEGE PREP CLASSES
# LH29B Q29B. HIGH SCHOOL TAKE VOCATIONAL CLASSES
# LH29C Q29C. HIGH SCHOOL MONITOR ATTENDANCE
# LH29D Q29D. HIGH SCHOOL DISRUPT CLASSES
# LH30 Q30. HIGH SCHOOL FOREIGN LANGUAGE
# LH31A Q31A. HIGH SCHOOL MUSICAL INSTRUMENT
# LH31B Q31B. HIGH SCHOOL SING
# LH31C Q31C. HIGH SCHOOL WOODWORK
# LH31D Q31D. HIGH SCHOOL CRAFT
# LH31E Q31E. HIGH SCHOOL DANCE
# LH31F Q31F. HIGH SCHOOL ART
cog_mat_hs <- lhms_selected %>% select(paste0("LH", 29, LETTERS[1:4]), LH30, 
                                       paste0("LH", 31, LETTERS[1:6])) %>%
  mutate_all(function(x) ifelse(x == 5, 0, x)) %>%
  mutate(cog_hs = case_when(
    is.na(LH29A) & is.na(LH29B) & is.na(LH29C) & is.na(LH29D) & is.na(LH30) &
      is.na(LH31A) & is.na(LH31B) & is.na(LH31C) & is.na(LH31D) & is.na(LH31E) 
    & is.na(LH31F) ~ NA_real_,
    TRUE ~ rowSums(., na.rm = T)))
view(cog_mat_hs %>% filter(cog_hs == 0))
view(cog_mat_hs %>% filter(is.na(cog_hs)))
summary(cog_mat_hs$cog_hs)
lhms_selected[colnames(cog_mat_hs)] <- cog_mat_hs

##---- Social activities ----
# LH31G Q31G. HIGH SCHOOL MATH OR SCIENCE CLUB
# LH31H Q31H. HIGH SCHOOL DRAFTING
# LH31I Q31I. HIGH SCHOOL TRADE CLASSES
# LH31J Q31J. HIGH SCHOOL DRAMA OR DEBATE CLUB
# LH33 Q33. HIGH SCHOOL SPORTS
social_mat_hs <- lhms_selected %>% select(paste0("LH", 31, LETTERS[7:10]),
                                          LH33) %>%
  mutate_all(function(x) ifelse(x == 5, 0, x)) %>%
  mutate(social_hs = case_when(
    is.na(LH31G) & is.na(LH31H) & is.na(LH31I) & is.na(LH31J) & is.na(LH33) ~ NA_real_,
    TRUE ~ rowSums(., na.rm = T)))
view(social_mat_hs %>% filter(social_hs == 0))
view(social_mat_hs %>% filter(is.na(social_hs)))
summary(social_mat_hs$social_hs)
lhms_selected[colnames(social_mat_hs)] <- social_mat_hs

#---- Leave behind questionnaires ----
# Check if RAND fat files have the variables
path_to_rand_fat <- paste0(path_to_box, "/HRS_MRG/RANDFatFiles")
years <- str_sub(as.character(seq(1998, 2018, 2)), start = 3)
paste(years, collapse = "|")
files <- list.files(path_to_rand_fat, full.names = TRUE) %>%
  purrr::keep(., ~str_detect(., "sas7bdat") & !str_detect(., "sasfmt") &
                str_detect(., paste(years, collapse = "|")))
list <- lapply(files, haven::read_sas)
names(list) <- gsub(".*/(.*)\\..*", "\\1", files)
# This is too large!
# hrs_fat <- join_all(list)
# only look at 2008 maybe!
hrs_08 <- list[[5]]

##---- Social ----
social_08 <- hrs_08 %>%
  select(LLB001E, LLB001F, LG086) %>%
  # LG195) %>%
  dplyr::mutate(LG086 = case_when(LG086 == 1 ~ 1,
                                  LG086 == 5 ~ 0,
                                  TRUE ~ NA_real_),
                across(LLB001E:LLB001F, function(x) case_when(x == 6 ~ 0,
                                                              x %in% 1:5 ~ 1))) %>%
  mutate(social_sum_08 = case_when(
    is.na(LLB001E) & is.na(LLB001F) & is.na(LG086) ~ NA_real_,
    TRUE ~ rowSums(., na.rm = T)))
summary(social_08)

##---- Cognitive ----
cog_08 <- hrs_08 %>%
  select(LLB001D, LLB001H, LLB001I, LLB001J, LLB001K, LLB001L) %>%
  dplyr::mutate(across(LLB001D:LLB001L, function(x) case_when(x == 6 ~ 0,
                                                              x %in% 1:5 ~ 1))) %>%
  mutate(cog_sum_08 = case_when(
    is.na(LLB001D) & is.na(LLB001H) & is.na(LLB001I) & is.na(LLB001J) 
    & is.na(LLB001K) & is.na(LLB001L) ~ NA_real_,
    TRUE ~ rowSums(., na.rm = T)))
summary(cog_08)

edu_08 <- hrs_08 %>% select(HHIDPN) %>% cbind(., cog_08, social_08)

#----all----
edu_prelim <- lhms_selected %>%
  select(HHIDPN, colnames(cog_mat_hs), colnames(social_mat_hs)) %>%
  left_join(edu_08, by = "HHIDPN")
summary(edu_prelim)

#---- RAND ----
number_waves <- 1:14
rand_variables <- c("hhidpn", "ragender", "raracem", "rahispan", "rabmonth", 
                    "rabyear", "rabdate", "radmonth", "radyear", "raddate",
                    paste0("r", number_waves, "agem_e"), "raedyrs", 
                    "r9wtresp",
                    paste0("r", number_waves, "iwstat"))
RAND <- read_dta(paste0(path_to_box, "/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2018v1.dta"),
                 col_select = all_of(rand_variables)) 
colnames(RAND)[1] <- "HHIDPN" #For merging

##---- Age ----
age_m <- RAND %>% dplyr::select(contains("agem_e")) %>% 
  apply(., 1, impute_ages) %>% t() 

#Exact ages
RAND[, paste0("r", number_waves, "age_y")] <- age_m/12
#Ages rounded down to nearest year
RAND[, paste0("r", number_waves, "age_y_int")] <- floor(age_m/12)

# Drop RAND age variables
RAND %<>% dplyr::select(-paste0("r", number_waves, "agem_e"))

##---- gender ----
RAND %<>% 
  mutate("female" = ifelse(ragender == 2, 1, 0))

#Drop RAND gender variable 
RAND %<>% dplyr::select(-ragender)

##---- race-eth ----
#Code any hispanic as 1, else 0
RAND %<>% 
  mutate(hispanic = ifelse(rahispan == 0 | is.na(rahispan), 0, 1), 
         white = ifelse(raracem == 1 & hispanic == 0, 1, 0),
         black = ifelse(raracem == 2 & hispanic == 0, 1, 0),
         other = ifelse(raracem == 3 & hispanic == 0, 1, 0),
         unknown_race_eth = ifelse(is.na(raracem) & hispanic == 0, 1, 0)) %>%
  mutate(raceeth = case_when(hispanic == 0 & white == 1 ~ "Non-Hispanic White",
                             hispanic == 0 & black == 1 ~ "Non-Hispanic Black",
                             hispanic == 1 ~ "Hispanic",
                             other == 1 ~ "Other"))

# #Sanity check
# table(RAND$hispanic, RAND$rahispan, useNA = "ifany")
# table(RAND$hispanic, RAND$raracem, RAND$black, useNA = "ifany")
# table(RAND$hispanic, RAND$raracem, RAND$other, useNA = "ifany")
# table(RAND$hispanic, RAND$raracem, RAND$unknown_race_eth,
#       useNA = "ifany")
# table(RAND$unknown_race_eth, useNA = "ifany")
# colSums(is.na(RAND[, c("hispanic", "white", "black", "other",
#                            "unknown_race_eth")]))
# table(RAND$white, RAND$raracem, useNA = "ifany")
# table(RAND$rahispan, RAND$hispanic, useNA = "ifany")

RAND %<>% 
  #Drop the RAND rahispan variable (recoded as hispanic) and race variables
  dplyr::select(-c("rahispan", "raracem"))

##---- education ----
#Create education categories
RAND %<>% 
  mutate("ed_cat" = 
           case_when(raedyrs < 12 ~ 1, 
                     raedyrs == 12 ~ 2, 
                     raedyrs > 12 & raedyrs < 16 ~ 3, 
                     raedyrs == 16 ~ 4, 
                     raedyrs > 16 ~ 5)) 



#---- merge with edu prelim ----
hrs_edu <- edu_prelim %>%
  select(-contains("LH")) %>%
  left_join(RAND %>% select(HHIDPN, r9age_y_int, female, raceeth, ed_cat, r9wtresp), 
            by = "HHIDPN")

#---- memory score ----
hrs_mem <- read_sas(paste0(path_to_box, 
                           "/HRS_MRG/HRSCleanData/ImputedMemoryScores/",
                           "dpmemimp_feb2023.sas7bdat")) %>%
                           # "dpmemimp_sep2022.sas7bdat")) %>%
  unite("HHIDPN", c("HHID", "PN"), sep = "", remove = TRUE) %>%
  mutate_at("HHIDPN", as.numeric)
hrs_edu_mem <- hrs_edu %>%
  left_join(hrs_mem, by = "HHIDPN")

#---- Dementia Probabilities ----
hrs_dem <- read_sas(paste0(path_to_box, 
                           "/HRS_MRG/HRSCleanData/“，
                           ”hrsdementia_2021_1109/",
                           "hrsdementia_2021_1109.sas7bdat")) %>%
  unite("HHIDPN", c("HHID", "PN"), sep = "", remove = TRUE) %>%
  mutate_at("HHIDPN", as.numeric)

hrs_dem_cleaned <- hrs_dem %>%
  arrange(HHIDPN, hurd_p) %>%
  group_by(HHIDPN) %>%
  slice_tail()
summary(hrs_dem_cleaned)
length(unique(hrs_dem_cleaned$HHIDPN))

hrs_edu_mem_dem <- hrs_edu_mem %>%
  left_join(hrs_dem_cleaned, by = "HHIDPN")

#---- Table 1 ----
hrs_edu_mem_dem %>%
  filter(r9age_y_int > 50) %>% nrow()

hrs_edu_mem_dem %>%
  filter(r9age_y_int > 50 &
           (!is.na(memimp08) |!is.na(memimp10)|!is.na(memimp12)|
              !is.na(memimp14)|!is.na(memimp16)|!is.na(memimp18))) %>%
  nrow()
hrs_edu_mem_dem %>%
  filter(r9age_y_int > 50 & (!is.na(hurd_dem)|!is.na(expert_dem)|
                               !is.na(lasso_dem))) %>% nrow()

hrs_edu_mem_dem %>%
  filter(r9age_y_int > 50 & !is.na(r9wtresp)) %>% nrow()
hrs_edu_svy <- hrs_edu_mem_dem %>%
  filter(r9age_y_int > 50 & !is.na(r9wtresp)) %>%
  select(r9age_y_int, female, raceeth, ed_cat, r9wtresp, cog_hs, social_hs,
         social_sum_08, cog_sum_08, memimp08,) %>%
  labelled::set_variable_labels(
    r9age_y_int = "Age in years at wave 9 (2008)",
    female = "Female (%)",
    raceeth = "Race/Ethnicity (%)",
    ed_cat = "Education level (%)",
    cog_hs = "Number of cognitive activities in High school",
    social_hs = "Number of social activities in High school",
    social_sum_08 = "Number of social activities in 2008",
    cog_sum_08 = "Number of cognitive activities in 2008",
    memimp08 = "Imputed memory score in 2008") %>%
  labelled::drop_unused_value_labels() %>%
  labelled::set_value_labels(female = c("Yes" = 1, "No" = 0),
                             ed_cat = c("Less than high school" = 1, "High school" = 2, 
                                        "Some college" = 3, "Bachelor's degree" = 4, 
                                        "Graduate studies" = 5)) %>%
  modify_if(is.labelled, to_factor) %>%
  # modify_if(is.factor, 
  #           function(x) {forcats::fct_explicit_na(x, na_level = "Missing")}) %>%
  survey::svydesign(~1, data = ., weights = ~r9wtresp)

table1 <- tbl_svysummary(
  data = hrs_edu_svy,
  statistic = list(all_categorical() ~ "{n_unweighted} ({p}%)",
                   all_continuous() ~ c("{mean} ({sd})",
                                        "({min}, {max})")),
  type = list(all_continuous() ~ "continuous2",
              c("cog_hs", "social_hs", "cog_sum_08", "social_sum_08") ~ "continuous2"),
  digits = list(all_continuous() ~ 1,
                all_categorical() ~ c(0, 1))) %>% 
  bold_labels() %>%
  modify_footnote(
    update = list(
      starts_with("stat_") ~
        "Mean (SD) for continuous variables; n (%) for categorical variables")
  ) 

table1 %>% as_flex_table()

table1 %>% as_flex_table() %>% flextable::save_as_docx(
  path = here::here("output", "tables", "table1_K00.docx"))

# save
save(hrs_edu_mem_dem, file = here::here("data", "analysis_data", "hrs_K.RData"))

#---- OLD ----
# path_to_core <- paste0(path_to_box, "HRS_MRG/HRSRawData/Core/")
# read_hrs_core <- function(year, section){
#   hrsdata <- paste0(path_to_core, "h", year, "core/h", year, "da/H", 
#                       year, section, ".DA")
#   # e.g. h98core/h98da/H98E_R.DA
#   hrsdict <- paste0(path_to_core, "h", year, "core/h", year, "sta/H", 
#                       year, section, ".DCT")
#   # e.g. h98core/h98sta/H98E_R.DCT
#   
#   # Read the dictionary file
#   df.dict <- read.table(hrsdict, skip = 1, fill = TRUE, 
#                           stringsAsFactors = FALSE) %>%
#     set_colnames(c("col.num","col.type","col.name","col.width","col.lbl"))
#   
#   # Remove last row which only contains a closing }
#   df.dict <- df.dict[-nrow(df.dict), ]
#   
#   # Extract numeric value from column width field
#   df.dict$col.width <- as.integer(
#     sapply(df.dict$col.width, gsub, pattern = "[^0-9\\.]", replacement = ""))
#   
#   # Convert column types to format to be used with read_fwf function
#   df.dict$col.type <- 
#     sapply(df.dict$col.type, 
#            function(x) case_when(x %in% c("int","byte","long") ~ "i",
#                                  x == "float" ~ "n",
#                                  x == "double" ~ "d",
#                                  TRUE ~ "c"))
#   # Read the data file into a dataframe
#   hrs <- read_fwf(file = hrsdata, 
#                     fwf_widths(widths = df.dict$col.width, 
#                                col_names = df.dict$col.name), 
#                     col_types = paste(df.dict$col.type, collapse = ""))
#   
#   # Add column labels to headers
#   attributes(hrs)$variable.labels <- df.dict$col.lbl
#   return(hrs)
# }
# 
# hrs98E_R <- read_hrs_core("98", "E_R")
# hrs00E_R <- read_hrs_core("00", "E_R")
# hrs02G_R <- read_hrs_core("02", "G_R")
# hrs18G_R <- read_hrs_core("18", "G_R")


