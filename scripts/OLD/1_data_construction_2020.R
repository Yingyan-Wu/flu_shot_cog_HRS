# HRS_data_construction
# Created by: Yingyan Wu
# Some code adapted from Crystal Shaw's and Mingfei Dong's Data construction scripts
# Dec.25.2022

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "haven", "labelled",
       "lubridate", "ggalluvial")

#---- function ----
source(here::here("scripts", "functions.R"))

#---- Load data ----
##---- RAND ----
number_waves <- 1:15
rand_variables <- c("hhidpn", "ragender", "raracem", "rahispan", "rabmonth", 
                    "rabyear", "rabdate", "radmonth", "radyear", "raddate",
                    "radage_y", "hacohort", # the cohort in which the household was originally sampled. 
                    "racohbyr", # cohort based on actual birth year 
                    # Response indicator
                    paste0("inw", number_waves),
                    paste0("r", number_waves, "iwbeg"),# Interview begin dates
                    paste0("r", number_waves, "nhmliv"), # living in nursing home at interview
                    paste0("r", number_waves, "agem_e"), "raedyrs", 
                    paste0("r", number_waves, "iwstat"),
                    paste0("r", number_waves, "wtresp"),
                    "raedegrm", 
                    paste0("r", 3:15, "flusht"),
                    paste0("r", number_waves, "mstat"),
                    paste0("r", number_waves, "bmi"),
                    paste0("r", number_waves, "hibpe"),
                    paste0("r", number_waves, "diabe"),
                    paste0("r", number_waves, "cancre"),
                    paste0("r", number_waves, "stroke"),
                    paste0("r", number_waves, "hearte"),
                    paste0("r", number_waves, "lunge"),
                    # "r1flonex", paste0("r", 2:15, "flone"),
                    paste0("r", 4:9, "memrye"),
                    paste0("r", 10:15, "demen"),
                    paste0("r", 10:15, "alzhee"),
                    paste0("r", number_waves, "conde"),
                    paste0("r", number_waves, "smoken"),
                    paste0("r", 3:15, "drinkd"),
                    paste0("r", 3:15, "drinkn"),
                    paste0("r", 2:15, "cesd"),
                    paste0("r", number_waves, "shlt"), # self rated health
                    # doctor visits
                    paste0("r", number_waves, "doctor"),
                    paste0("r", number_waves, "doctim"),
                    # Physical activities
                    paste0("r", 1:6, "vigact"), # vigorous physical activities 3+/week
                    paste0("r", 7:15, "vgactx"), # Freq vigorous phys activities
                    
                    # Daily functioning
                    paste0("r", 2:15, "adl6a"),
                    paste0("r", 2:15, "iadl5a"),
                    # Household income
                    paste0("h", number_waves, "itot"),
                    paste0("h", number_waves, "hhres"), # household size
                    # Health insurance
                    paste0("r", number_waves, "higov"), # any government plan
                    # medicare, medicaid, military, VA, other
                    paste0("r", number_waves, "prpcnt")) # private insurance plans

RAND <- read_dta(here::here("data", "raw_data", "randhrs1992_2020v2_STATA",
                            "randhrs1992_2020v2.dta"),
                 col_select = all_of(rand_variables)) 

colnames(RAND)[1] <- "HHIDPN" #For merging

##---- Imputed memory score and dem probabilities ----
memimp <- read_sas(here::here("data", "raw_data", "imputed_memory_scores", 
                              "dpmemimp_jan2024.sas7bdat")) %>%
  unite("HHIDPN", c("HHID", "PN"), sep = "", remove = TRUE) %>%
  mutate_at("HHIDPN", as.numeric) %>%
  select(-contains("95"))

# Mapping years to waves
# memimp start in 1996
HRS_years <- seq(1996, 2020, by = 2) # Years
waves <- (as.numeric(HRS_years) - 2000)/2 + 5  #should be 3:15 
# Function to rename variables
rename_to_rwxxx <- function(var_names, years, waves) {
  # Create lookup for year-to-wave
  year_wave_map <- setNames(waves, substr(years, 3, 4))
  # Pattern match and replace
  str_replace(var_names, "(\\D+)(\\d{2})$", function(x) {
    prefix <- str_match(x, "(\\D+)")[, 1]
    year <- str_match(x, "(\\d{2})$")[, 1]
    paste0("r", year_wave_map[year], prefix)
  })
}

memimp %<>%
  # years from 1996-2016 corresponding to waves 3 - 15
  # Change the name to rwmemimp etc.
  dplyr::rename_at(vars(contains(c("memimp", "dementpimp"))),
                   ~rename_to_rwxxx(.x, HRS_years, waves))

# Sanity check
colnames(memimp)
memimp %>% select(contains("dementpimp")) %>% apply(., 2, mean, na.rm = T)

# demprob <- read_sas(here::here("data", "raw_data", "power_hrsdementia",
#                                "hrsdementia_2024_1217.sas7bdat")) %>%
#   unite("HHIDPN", c("HHID", "PN"), sep = "", remove = TRUE) %>%
#   mutate_at("HHIDPN", as.numeric) %>%
#   # years from 2000-2016 corresponding to waves 5 - 13
#   # (seq(2000, 2016, by = 2) - 2000)/2 + 5
#   mutate(wave = (as.numeric(HRS_year) - 2000)/2 + 5)
# 
# # Sanity check
# with(demprob, table(HRS_year, wave, useNA = "ifany"))
# 
# # Include only Expert algorithm variables
# demprob_wide_expert <- demprob %>%
#   mutate(wave = (as.numeric(HRS_year) - 2000)/2 + 5) %>%
#   select(HHIDPN, wave, expert_p, expert_dem) %>%
#   pivot_wider(names_from = wave, values_from = c(expert_p, expert_dem),
#               names_glue = "r{wave}{.value}")
# 
# with(demprob %>%
#        arrange(HHIDPN, HRS_year) %>%
#        group_by(HHIDPN) %>%
#        slice_tail(), table(hurd_dem))
# 
# with(demprob %>%
#        arrange(HHIDPN, HRS_year) %>%
#        group_by(HHIDPN) %>%
#        slice_tail(), table(hurd_dem))
# 
# demprob_wide_expert %>% select(contains("expert_p")) %>% apply(., 2, mean, na.rm = T)

##---- merging data ----
hrs <- RAND %>% left_join(memimp, by = "HHIDPN")
  # left_join(demprob_wide_expert, by = "HHIDPN")

#---- in interview ----
hrs %<>%
  rename_with(~ str_replace(.x, "^inw(\\d+)$", "r\\1inw"), starts_with("inw"))
# # Sanity check
# colnames(hrs %>% select(contains("inw")))
# # Check response indicator and interview dates
# hrs %>% select(contains("inw")) %>%
#   apply(., 2, table, useNA = "ifany")
# 
# for (wave in 3:15){
#   print(hrs %>% filter(!!paste0("r", wave, "inw") == 0 & 
#                    !is.na(paste0("r", wave, "iwbeg"))) %>% nrow())
# }

# Interview dates
hrs %<>%
  mutate(across(c(contains("iwbeg"), raddate, rabdate), 
                ~as.Date(.x, origin = "1960-01-01")))

# Exact age at waves (to calculate exact follow up time)
for (w in number_waves){
  hrs%<>%
    mutate(!!paste0("r", w, "fu_age") := # compare to birthdate
             round(interval(start = rabdate,
                            end = !!sym(paste0("r", w, "iwbeg"))) /
                     duration(num = 1, units = "years"), 1))
}

#---- IWSTAT ----
# Check cross wave to overturn unknown status (iwstat = 9)
hrs %>% select(contains("iwstat")) %>%
  apply(., 2, table, useNA = "ifany")
# Check these 2 observations with iwstat = 9
hrs %>% select(HHIDPN, contains("iwstat")) %>%
  filter(r8iwstat == 9 | r9iwstat == 9) %>% as.data.frame()
# Since this person attended wave 10-wave 12, change r8 and r9 iwstat to 4!

hrs %<>%
  mutate(r8iwstat = case_when(HHIDPN == "74494020" & r8iwstat == 9 ~ 4,
                              TRUE ~ r8iwstat),
         r9iwstat = case_when(HHIDPN == "74494020" & r9iwstat == 9 ~ 4,
                              TRUE ~ r9iwstat))

#---- Death ----
# Create death indicator by wave
for (w in number_waves){
  hrs %<>%
    mutate(!!paste0("r", w, "dead") := 
             case_when(!!sym(paste0("r", w, "iwstat")) %in% c(5, 6) ~ 1,
                       !!sym(paste0("r", w, "iwstat")) == 9 ~ NA_real_, # NR, DK if alive or died
                       TRUE ~ 0))
}

# Sanity check
hrs %>% select(ends_with("dead")) %>% apply(., 2, table, useNA = "ifany")

for (wave in number_waves){
  print(with(hrs, table(get(paste0("r", wave, "inw")), 
                        get(paste0("r", wave, "iwstat")), useNA = "ifany")))
}

for (wave in number_waves){
  print(with(hrs, table(get(paste0("r", wave, "dead")), 
                        !is.na(raddate), useNA = "ifany")))
}

for (wave in number_waves){
  print(with(hrs, table(get(paste0("r", wave, "dead")), 
                        get(paste0("r", wave, "iwstat")), useNA = "ifany")))
}

# Covariates:
# Time-independent: age at baseline, sex/gender, race/ethnicity, education
# Time-dependent: marital status, depressive symptoms, household income, ADLs 
# and IADLs, drinking, doctor visits, physical activities,
# health conditions: diabetes, hypertension, heart disease, stroke, 
# cancer and lung disease, self-rated health

#---- Age ----
age_m <- hrs %>% dplyr::select(contains("agem_e")) %>% 
  apply(., 1, impute_ages) %>% t() 

#Exact ages
hrs[, paste0("r", number_waves, "age_y")] <- age_m/12

# Drop RAND age variables
hrs %<>% dplyr::select(-paste0("r", number_waves, "agem_e"))
rm(age_m)

# Use fu_age for exact age, unless it's missing, then use computed age
for (w in number_waves){
  hrs %<>%
    mutate(!!paste0("r", w, "age_v2") := 
             case_when(!is.na(!!sym(paste0("r", w, "fu_age"))) ~ 
                         !!sym(paste0("r", w, "fu_age")),
                       TRUE ~ !!sym(paste0("r", w, "age_y"))))
}

# # Sanity check
# hrs %>%
#   select(contains("age_v2")) %>%
#   apply(., 2, summary)
# 
# # Check if participants have the same age at different waves
# hrs %>%
#   select(HHIDPN, contains("age_v2")) %>%
#   pivot_longer(cols = contains("age_v2"),
#                names_to = "wave", values_to = "age", 
#                names_pattern = "r(\\d+)age_v2") %>%
#   mutate(newid = paste0(HHIDPN, age)) %>%
#   group_by(newid) %>%
#   count() %>% view()
# 
# hrs %>% select(HHIDPN, contains("age_v2")) %>% filter(HHIDPN == "32570030")
# hrs %>% select(HHIDPN, contains("age_v2")) %>% filter(HHIDPN == "32731010")

summary(hrs$rabyear)
summary(hrs$r3age_v2)

#---- Sex/gender ----
hrs %<>% 
  mutate("female" = ifelse(ragender == 2, 1, 0))

# Drop RAND gender variable 
hrs %<>% dplyr::select(-ragender)

#---- Race/Ethnicity ----
#Code any hispanic as 1, else 0
hrs %<>% 
  mutate(hispanic = ifelse(rahispan == 0 | is.na(rahispan), 0, 1), 
         white = ifelse(raracem == 1 & hispanic == 0, 1, 0),
         black = ifelse(raracem == 2 & hispanic == 0, 1, 0),
         other = ifelse(raracem == 3 & hispanic == 0, 1, 0),
         unknown_race_eth = ifelse(is.na(raracem) & hispanic == 0, 1, 0)) %>%
  mutate(raceeth = case_when(hispanic == 0 & white == 1 ~ "1", # Non-Hispanic White
                             hispanic == 0 & black == 1 ~ "2", # Non-Hispanic Black
                             hispanic == 1 ~ "3", # Hispanic
                             other == 1 ~ "4")) # Other

# #Sanity check
# table(hrs$hispanic, hrs$rahispan, useNA = "ifany")
# table(hrs$hispanic, hrs$raracem, hrs$black, useNA = "ifany")
# table(hrs$hispanic, hrs$raracem, hrs$other, useNA = "ifany")
# table(hrs$hispanic, hrs$raracem, hrs$unknown_race_eth,
#       useNA = "ifany")
# table(hrs$unknown_race_eth, useNA = "ifany")
# colSums(is.na(hrs[, c("hispanic", "white", "black", "other",
#                            "unknown_race_eth")]))
# table(hrs$white, hrs$raracem, useNA = "ifany")
# table(hrs$rahispan, hrs$hispanic, useNA = "ifany")

hrs %<>% 
  #Drop the RAND rahispan variable (recoded as hispanic) and race variables
  dplyr::select(-c("rahispan", "raracem"))

#---- Education ----
#Create education categories
hrs %<>% 
  mutate("ed_cat" = 
           case_when(raedyrs < 12 ~ 1, 
                     raedyrs == 12 ~ 2, 
                     raedyrs > 12 & raedyrs < 16 ~ 3, 
                     raedyrs == 16 ~ 4, 
                     raedyrs > 16 ~ 5)) 

# #Sanity check
# table(is.na(hrs$raedyrs))
# table(hrs$raedyrs)
# table(hrs$raedyrs, hrs$ed_cat, useNA = "ifany")

#---- Marital status ----
# Create marital status categories
for (wave in number_waves){
  hrs %<>% 
    mutate(
      !!paste0("r", wave, "mstat_2cat") := 
        case_when(
          # Married/partnered
          !!sym(paste0("r", wave, "mstat")) %in% c(1, 2, 3) ~ 1, 
          # Separated/divorced/widowed/not married/Widowed
          !!sym(paste0("r", wave, "mstat")) %in% c(4, 5, 6, 7, 8) ~ 0))
}

# Sanity check
for (wave in number_waves){
  print(with(hrs, table(get(paste0("r", wave, "mstat_2cat")),
                        get(paste0("r", wave, "mstat")), useNA = "ifany")))
}

#---- Household income ----
# Calculate individual level household income
# total income/sqrt(household size)
for (wave in number_waves){
  hrs %<>% 
    mutate(
      !!paste0("r", wave, "income_pp") := 
        !!sym(paste0("h", wave, "itot"))/sqrt(!!sym(paste0("h", wave, "hhres"))))
}

# # Sanity check
# hrs %>% select(contains("hhres")) %>% apply(., 2, summary)
# hrs %>% select(contains("income_pp")) %>% apply(., 2, summary)
# hrs %>% select(contains("itot")) %>% apply(., 2, summary)

#---- Insurance info ----
for (wave in number_waves){
  hrs %<>% 
    mutate(
      !!paste0("r", wave, "insurance_2cat") := 
        case_when(
          # NAed both government and private insurance => NA
          is.na(!!sym(paste0("r", wave, "higov"))) &
            is.na(!!sym(paste0("r", wave, "prpcnt"))) ~ NA_real_,
          # have insurance including government and private
          # medicare, medicaid, military, VA, other and private
          !!sym(paste0("r", wave, "higov")) == 1 |
            !!sym(paste0("r", wave, "prpcnt")) > 0 ~ 1, 
          # No government, private insurance
          !!sym(paste0("r", wave, "higov")) == 0 &  
            !!sym(paste0("r", wave, "prpcnt")) == 0 ~ 0))
}

# # Sanity check
# hrs %>% select(contains("insurance_2cat")) %>% apply(., 2, table, useNA = "ifany")
# hrs %>% select(contains("higov")) %>% apply(., 2, table, useNA = "ifany")
# hrs %>% select(contains("prpcnt")) %>% apply(., 2, table, useNA = "ifany")

#---- Number of doctor visits ----
# Since wave 1 asks about doctor visits in prior 12 months, the others asked about 
# prior 2 years (and not using wave 1 anyway), start with cleaning wave 2 doctors

hrs %>% select(contains("doctor")) %>% apply(., 2, table, useNA = "ifany")
hrs %>% select(contains("doctim")) %>% apply(., 2, summary)
# with(hrs, table(r3doctor, r3doctim, useNA = "ifany"))

# categorize at more than 5 visits (median)
for (wave in 2:15){
  hrs %<>% 
    mutate(
      !!paste0("r", wave, "doct_vis_mt5") := 
        case_when(
          !!sym(paste0("r", wave, "doctim")) >= 5 ~ 1,
          !!sym(paste0("r", wave, "doctim")) < 5 |
            !!sym(paste0("r", wave, "doctor")) == 0 ~ 0,
          is.na(!!sym(paste0("r", wave, "doctim"))) ~ NA_real_))
}

# Sanity check
hrs %>%
  select(contains("doct_vis_mt5")) %>%
  apply(., 2, table, useNA = "ifany")
with(hrs %>% filter(is.na(r15doct_vis_mt5)), table(r15doctor, r15doctim, useNA = "ifany"))

#---- Smoking ----
# Smokes now
hrs %>%
  select(contains("smoken")) %>%
  apply(., 2, table, useNA = "ifany")

#---- Depressive symptoms ----
# CES-D
hrs %>%
  select(contains("cesd")) %>%
  apply(., 2, table, useNA = "ifany")

# #---- Loneliness ----
# hrs %<>%
#   mutate(r1flone = case_when(r1flonex %in% 1:3 ~ 1,
#                              r1flonex == 4 ~ 0))
# # Sanity check
# with(hrs, table(r1flonex, r1flone, useNA = "ifany"))
# hrs %>%
#   select(paste0("r", number_waves, "flone")) %>%
#   apply(., 2, table, useNA = "ifany")

#---- Drinking ----
# Based on 2000 Dietary Guidelines for Americans: 
# https://odphp.health.gov/sites/default/files/2020-01/DGA2000.pdf

# Moderation is defined as no more than one drink
# per day for women and no more than two drinks per day for men.

# drinks/day when drinks
colSums(is.na(hrs %>% select(contains("drinkn"))))

for (wave in 3:15){
  hrs %<>% 
    mutate(
      !!paste0("r", wave, "mod_drinking") := 
        case_when(female == 1 & !!sym(paste0("r", wave, "drinkn")) <= 1 ~ 1,
                  female == 1 & !!sym(paste0("r", wave, "drinkn")) > 1 ~ 0,
                  
                  female == 0 & !!sym(paste0("r", wave, "drinkn")) <= 2 ~ 1,
                  female == 0 & !!sym(paste0("r", wave, "drinkn")) > 2 ~ 0,
                  !!sym(paste0("r", wave, "drinkd")) == 0 ~ 1,
                  is.na(!!sym(paste0("r", wave, "drinkd"))) &
                    is.na(!!sym(paste0("r", wave, "drinkn"))) ~ NA_real_))
}

# Sanity check
hrs %>%
  select(contains("mod_drinking")) %>%
  apply(., 2, table, useNA = "ifany")
# hrs %>% select(HHIDPN, contains("r3drink"), r3mod_drinking) %>% 
#   filter(is.na(r3mod_drinking)) %>% view()
# drinks / day
colSums(is.na(hrs %>% select(contains("drinkn"))))
# days/ week
colSums(is.na(hrs %>% select(contains("drinkd"))))

#---- Self-rated health ----
hrs %>% select(contains("shlt")) %>% apply(., 2, table, useNA = "ifany")
for (wave in number_waves){
  hrs %<>% 
    mutate(
      !!paste0("r", wave, "shlt_fair_poor") := 
        case_when(
          # fair or poor
          !!sym(paste0("r", wave, "shlt")) %in% c(4, 5) ~ 1,
          # good, very good, excellent
          !!sym(paste0("r", wave, "shlt")) %in% c(1, 2, 3) ~ 0))
}
# Sanity check
hrs %>% 
  select(contains("shlt_fair_poor")) %>% 
  apply(., 2, table, useNA = "ifany")

#---- Physical activities ----
hrs %>% select(contains("vigact")) %>% apply(., 2, table, useNA = "ifany")
hrs %>% select(contains("vgactx")) %>% apply(., 2, table, useNA = "ifany")

for (wave in number_waves){
  if (wave <= 6){
    hrs %<>% 
      mutate(!!paste0("r", wave, "vigorous_pa") := !!sym(paste0("r", wave, "vigact")))
  } else if (wave %in% 7:15){
    hrs %<>% 
      mutate(
        !!paste0("r", wave, "vigorous_pa") := 
          case_when(
            !!sym(paste0("r", wave, "vgactx")) %in% c(1, 2) ~ 1, # every day, >1 per week
            !!sym(paste0("r", wave, "vgactx")) %in% c(3, 4, 5) ~ 0, # <= 1 per week
            is.na(!!sym(paste0("r", wave, "vgactx"))) ~ NA_real_))
  }
}
# Sanity check
hrs %>% select(contains("vigorous_pa")) %>% apply(., 2, table, useNA = "ifany")


#---- Health conditions ----
# diabetes, hypertension, heart disease, stroke, cancer and lung disease
hrs %>%
  select(contains(c("hibpe", "diabe", "cancre", "stroke", "hearte", "lunge"))) %>%
  apply(., 2, table, useNA = "ifany")

#---- ADLs ----
# bathing, dressing, eating, getting in/out of bed,
# and walking across a room, toileting
hrs %>%
  select(contains("adl6a")) %>%
  apply(., 2, table, useNA = "ifany")

for (wave in 2:15){
  hrs %<>% 
    mutate(
      !!paste0("r", wave, "adl_2cat") := 
        case_when(
          !!sym(paste0("r", wave, "adl6a")) == 0 ~ 0,
          !!sym(paste0("r", wave, "adl6a")) > 0 ~ 1))
}

# # Sanity check
# for (wave in 2:15){
#   print(with(hrs, table(get(paste0("r", wave, "adl_2cat")),
#                         get(paste0("r", wave, "adl6a")), useNA = "ifany")))
# }

#---- IADLs ----
# difficulty using the phone, managing money, taking medications, 
# shopping for groceries, and preparing hot meals
hrs %>%
  select(contains("iadl5a")) %>%
  apply(., 2, table, useNA = "ifany")

for (wave in 2:15){
  hrs %<>% 
    mutate(
      !!paste0("r", wave, "iadl_2cat") := 
        case_when(
          !!sym(paste0("r", wave, "iadl5a")) == 0 ~ 0,
          !!sym(paste0("r", wave, "iadl5a")) > 0 ~ 1))
}

# # Sanity check
# for (wave in 2:15){
#   print(with(hrs, table(get(paste0("r", wave, "iadl_2cat")),
#                         get(paste0("r", wave, "iadl5a")), useNA = "ifany")))
# }

#---- FLU SHOT ----
##---- Check the flusht ----
flusht_wide <- RAND %>% select(HHIDPN, paste0("r", 3:15, "flusht"))
# Check the distribution
# view(flusht_wide %>% filter(!is.na(r4flusht)))
flusht_wide %>% select(-HHIDPN) %>% apply(., 2, function(x) table(!is.na(x)))

# Check the refreshers/new participants flushot question
# From 1998 (w4) to 2018 (w14)
# Most of them answered this question
for (i in 4:14){
  print(with(RAND %>%
               filter(rowSums(across(paste0("r", 1:(i-1), "iwstat"))) == 0 &
                        get(paste0("r", i, "iwstat")) == 1),
             table(!is.na(get(paste0("r", i, "flusht"))))))
}

# Check if the every 4 year rule is followed in the data
seq_nonNAs <- function(x){
  y <- rle(x) # un length encoding, return values and corresponding lengths
  seq_nonNAs <-  case_when(max(y$lengths[!is.na(y$values)]) > 1 ~ 1,
                           sum(is.na(y$values)) == length(y$values) ~ NA_real_,
                           TRUE ~ 0)
  return(seq_nonNAs)
}
# # Test
# test <- c(0, NA, 1, NA, 1, 1, NA, 0, NA, 0, 1, NA, NA)
# test <- c(0, NA, 0, NA, NA)
# seq_nonNAs(test)

flusht_wide$seq_flusht <- flusht_wide %>%
  select(-HHIDPN, -contains("flusht_")) %>%
  apply(., 1, seq_nonNAs)
table(flusht_wide$seq_flusht)
# 9326 participants have sequential flu shot status.

# Number of flushot taken
# All flu shot responses
flusht_wide %<>%
  mutate(num_flusht = case_when(rowSums(!is.na(across(contains("flusht")))) == 0 ~ NA_real_,
                                TRUE ~ rowSums(across(contains("flusht")), na.rm = TRUE)))
with(flusht_wide, table(num_flusht, useNA = "ifany"))

# Only responses at alternative waves from waves 3 to 15
flusht_wide %<>%
  mutate(num_flusht_alt = case_when(
    rowSums(!is.na(across(paste0("r", seq(3, 15, by = 2), "flusht")))) == 0 ~ NA_real_,
    TRUE ~ rowSums(across(paste0("r", seq(3, 15, by = 2), "flusht")), na.rm = TRUE)))
with(flusht_wide, table(num_flusht_alt, useNA = "ifany"))

# flusht_wide %>% filter(seq_flusht == 1) %>% view()

# Flu shot wide 

#---- From the same household ----
hrs %<>%
  mutate(
    # PN: last 3 digits of HHIDPN
    pn = gsub("^(\\d+)(\\d{3})$", "\\2", HHIDPN),
    # HHID: anything in HHIDPN besides pn
    hhid = gsub("^(\\d+)(\\d{3})$", "\\1", HHIDPN))


# Check for hhidpn and hhid, pn
hrs %>% select(HHIDPN, hhid, pn) %>% view()
length(unique(hrs$hhid))

#---- Dropping people ----
##---- 1. Dropping those who still missing age at e every wave ----
table(hrs$female, useNA = "ifany")
hrs %>%
  select(contains("age_v2")) %>%
  apply(., 2, summary)
with(hrs, table(is.na(rabdate))) # n = 1 missing age

hrs_filtered <- hrs %>%
  filter(!is.na(rabdate)) # n = 42404

##---- 2. Filter on participants eligible for 2000 interview ----
# Filter on participants eligible for 2000 interview
# and participated in 1998 interview to be able to have prebaseline covariate values
r5_resp_id <- hrs_filtered %>% filter(r5iwstat %in% c(1, 4)) %>% pull(HHIDPN) %>% unique()
# filter(wave == 3, iwstat %in% c(1, 4)) %>% pull(HHIDPN) %>% unique()
r4_notinsamp_id <-  hrs_filtered %>% filter(r4iwstat == 0) %>% pull(HHIDPN) %>% unique()
r5_resp_r4insamp_id <- r5_resp_id[!r5_resp_id %in% r4_notinsamp_id] # 21918

hrs_filtered %<>% filter(HHIDPN %in% r5_resp_r4insamp_id)
# n=21918
# with(hrs_filtered, table(racohbyr, useNA = "ifany"))

##---- 3. Filter on participants aged [52, 90) for 2000 interview ----
# Since they need to be eligible for both 1998 and 2000 interview
# people aged < 52 in 2000, these people were in the sample since they are spouses

# check if we can do 50+ in 2000
hrs_filtered %>% filter(r5age_v2 >= 50, r5age_v2 < 52) %>% with(., sum(r5wtresp, na.rm = T))
# n=347, weighted n=55272 (too small to represent [50, 52) people)
hrs_filtered %>% filter(r5age_v2 >= 52, r5age_v2 < 54) %>% with(., sum(r5wtresp, na.rm = T))
# n=953, weighted n=5380803 for [52, 54) people

hrs_filtered %>%
  filter(r5age_v2 < 52) %>% nrow() # n=941 (age younger than 52)
hrs_filtered %>%
  filter(r5age_v2 >= 90) %>% nrow() # n=471 (age older than 90)

hrs_filtered %<>%
  filter(r5age_v2 >= 52 & r5age_v2 < 90) #n=20506

##---- 4. filter on participants not living in nursing home for 2000 ----
# Filter on participants not living in nursing home for 2000 
# (if missing, not living in nursing home in previous waves)
with(hrs_filtered, table(r5nhmliv, r5iwstat, useNA = "ifany"))
with(hrs_filtered %>% filter(r5iwstat == 4), 
     table(r5nhmliv, r4nhmliv, useNA = "ifany"))
with(hrs_filtered %>% filter(r5iwstat == 4, is.na(r5nhmliv)), 
     table(r3nhmliv, r4nhmliv, useNA = "ifany"))
hrs_filtered %<>%
  mutate(
    r5nursing_home_v2 = 
      case_when(
        r5nhmliv == 1 ~ 1, # living in nursing home
        is.na(r5nhmliv) & r4nhmliv == 1 ~ 1,
        is.na(r5nhmliv) & is.na(r4nhmliv) & r3nhmliv == 1 ~ 1,
        is.na(r5nhmliv) & is.na(r4nhmliv) & is.na(r3nhmliv) & r2nhmliv == 1 ~ 1,
        is.na(r5nhmliv) & is.na(r4nhmliv) & is.na(r3nhmliv) & is.na(r2nhmliv) &
          r1nhmliv == 1 ~ 1,
        is.na(r5nhmliv) & is.na(r4nhmliv) & is.na(r3nhmliv) & is.na(r2nhmliv) &
          is.na(r1nhmliv) ~ NA_real_, 
        TRUE ~ 0))# not living in nursing home

with(hrs_filtered, table(r5nursing_home_v2, useNA = "ifany")) # n=346

hrs_filtered %<>%
  filter(r5nursing_home_v2 != 1) # n = 20091

#---- Select the variables ----
hrs_selected <- hrs_filtered %>%
  select(HHIDPN, raddate, hacohort, racohbyr, rabyear,
         contains("inw"),
         contains("iwstat"),
         ends_with("dead"),
         contains("wtresp"),
         contains("age_v2"),
         female, raceeth, raedyrs, ed_cat, 
         contains(c("income_pp", "itot", "hhres")),
         contains(c("mstat", "insurance")),
         contains(c("smoken", "cesd")), # "flone")), -(contains("flonex")),
         contains("mod_drinking"),
         contains("vigorous_pa"),
         contains("shlt_fair_poor"),
         contains(c("doct_vis_mt5", "doctor")),
         contains(c("hibpe", "diabe", "cancre", "stroke", "hearte", "lunge")),
         contains(c("adl", "iadl")),
         contains("flusht"),
         contains(c("memimp", "dementpimp")) # "expert_p", "expert_dem"))
  ) %>%
  # Change the names of household variable (h) to individual variable (r)
  rename_at(vars(contains(c("hhres", "itot"))), ~str_replace(.x, "h", "r")) # This will only replace the first h letter

#---- Save the data ----
saveRDS(hrs, here::here("data", "analysis_data", "hrs_cleaned.rds"))
save(hrs_filtered, hrs_selected, file = here::here("data", "analysis_data", "hrs_selected_2000_ot52l90.RData"))

#---- OLD ----
##---- 2. Filter on participants eligible for 1996 interview ----
# # Filter on participants eligible for 1996 interview
# # and participated in 1994 interview to be able to have covariate values
# r3_resp_id <- hrs_filtered %>% filter(r3iwstat %in% c(1, 4)) %>% pull(HHIDPN) %>% unique()
# # filter(wave == 3, iwstat %in% c(1, 4)) %>% pull(HHIDPN) %>% unique()
# r2_notinsamp_id <-  hrs_filtered %>% filter(r2iwstat == 0) %>% pull(HHIDPN) %>% unique()
# r3_resp_r2insamp_id <- r3_resp_id[!r3_resp_id %in% r2_notinsamp_id] # 19845
# # hrs_id <- hrs_filtered %>% filter(hacohort == 3) %>% pull(HHIDPN) %>% unique()
# # table(r3_resp_r2insamp_id %in% hrs_id) # n=12362
# hrs_filtered <- hrs_filtered %>% filter(HHIDPN %in% r3_resp_r2insamp_id)
# # n=19845
# ##---- Creating sub cohorts ----
# # 4 subcohorts of flushot
# # Sub cohort 1: 1996, 2000, 2004
# # Sub cohort 2: 2000, 2004, 2008
# # Sub cohort 3: 2004, 2008, 2012
# # Sub cohort 4: 2008, 2012, 2016
# subcohorts <- c("960004", "000408", "040812", "081216")
# subcohort_hrswaves <- list(seq(3, 7, by = 2),
#                            seq(5, 9, by = 2),
#                            seq(7, 11, by = 2),
#                            seq(9, 13, by = 2))
# 
# for (i in 1:length(subcohorts)){
#   hrs %<>%
#     mutate(!!paste0("flusht_", subcohorts[[i]]) := 
#              case_when(rowSums(is.na(across(
#                paste0("r", subcohort_hrswaves[[i]], "flusht")))) == 0 ~ 1,
#                TRUE ~ 0))
#   
#   flusht_temp <- hrs %>%
#     filter(!!sym(paste0("flusht_", subcohorts[[i]])) == 1) %>%
#     select(HHIDPN, ends_with("flusht")) %>%
#     mutate(w1_flusht = !!sym(paste0("r", subcohort_hrswaves[[i]][[1]], "flusht")),
#            w3_flusht = !!sym(paste0("r", subcohort_hrswaves[[i]][[2]], "flusht")),
#            w5_flusht = !!sym(paste0("r", subcohort_hrswaves[[i]][[3]], "flusht")),
#            w1_hrswave = subcohort_hrswaves[[i]][[1]],
#            w3_hrswave = subcohort_hrswaves[[i]][[2]],
#            w5_hrswave = subcohort_hrswaves[[i]][[3]],
#            subcohort = subcohorts[[i]])
#   assign(paste0("flusht_c", i), flusht_temp)
# }
# 
# flusht_cohorts <- bind_rows(flusht_c1, flusht_c2, flusht_c3, flusht_c4)
# # Sanity check
# hrs %>% select(starts_with("flusht")) %>% 
#   apply(., 2, table, useNA = "ifany")

# # Sub cohort traj
# flusht_cohorts_fortraj <- flusht_cohorts %>%
#   mutate(across(contains("flusht"), ~as.factor(case_when(is.na(.) ~ 999,
#                                                          TRUE ~ .)))) %>%
#   group_by(across(all_of(paste0("w", seq(1, 5, by = 2), "_flusht")))) %>%
#   summarize(count = n()) %>% ungroup()
# 
# flusht_cohorts_fortraj %>%
#   ggplot(aes(y = count,
#              axis1 = w1_flusht,
#              axis2 = w3_flusht,
#              axis3 = w5_flusht)) +
#   geom_alluvium(aes(fill = w1_flusht)) +
#   geom_stratum(alpha = 0) +
#   geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
#   labs(fill = "Flu shot acceptance") +
#   scale_fill_brewer(type = "qual", palette = "Paired",
#                     labels = c("Not accepted", "Accepted", "Missing")) +
#   xlab("Waves") +
#   ylab("Count") +
#   theme_minimal()
# ##---- Create aux variables ----
# # Long dataset
# flusht_long <- RAND %>%
#   select(HHIDPN, contains("flusht")) %>%
#   pivot_longer(!HHIDPN,
#                names_to = "wave",
#                names_pattern = "r(.*)flusht",
#                values_to = "flusht") %>%
#   mutate(wave = as.numeric(wave)) %>%
#   na.omit()
# 
# # Num of taking the flu shot
# flusht_times <- flusht_long %>%
#   group_by(HHIDPN) %>%
#   dplyr::summarise(n_flusht = sum(flusht))
# 
# summarytools::freq(flusht_times$n_flusht)
# 
# # Num of reported flu shot
# flusht_report <- flusht_long %>%
#   dplyr::count(HHIDPN) %>%
#   dplyr::rename(n_flusht_wave = n)
# 
# summarytools::freq(flusht_report$n_flusht_wave)
# 
# # The trajectory of flu shot status among (100 random participants)
# set.seed("62283")
# random_id <- sample(flusht_times$HHIDPN, 100, replace = T)
# 
# flusht_long %>%
#   # filter(HHIDPN %in% random_id) %>%
#   arrange(wave) %>%
#   ggplot(aes(x = wave, y = flusht, group = as.factor(HHIDPN),
#              color = as.factor(HHIDPN))) +
#   geom_line() +
#   geom_point() +
#   theme_minimal() +
#   scale_x_continuous(limits = c(3, 14), 
#                      breaks = seq(3, 14, by = 1)) +
#   theme(legend.position = "none")
# 
# #---- Filter to people have at least 2 wave flu shot status ----
# flu_2_id <- flusht_report %>% filter(n_flusht_wave >= 2) %>%
#   select(HHIDPN) %>% pull()
# 
# flu_2_3_8_id <- flusht_long %>%
#   filter(wave <= 8) %>%
#   dplyr::count(HHIDPN) %>%
#   dplyr::rename(n_flusht_wave = n) %>% filter(n_flusht_wave >= 2) %>%
#   select(HHIDPN) %>% pull()
# 
# flusht_wide_2 <- RAND %>%
#   filter(HHIDPN %in% flu_2_id) %>%
#   select(HHIDPN, contains("flusht")) %>%
#   left_join(flusht_times, by = "HHIDPN") %>%
#   left_join(flusht_report, by = "HHIDPN") %>%
#   mutate(prop_flusht = n_flusht/n_flusht_wave)
# 
# # People at least 2 waves: how many people didn’t change their status
# summarytools::freq(flusht_wide_2$prop_flusht)
# 
# # Prop of flu shot/all the available flu shot status (people at least two waves)
# hist(flusht_wide_2$prop_flusht)
# 
# flusht_wide_2_3_8 <- RAND %>%
#   filter(HHIDPN %in% flu_2_3_8_id) %>%
#   select(HHIDPN, contains("flusht")) %>%
#   left_join(flusht_times, by = "HHIDPN") %>%
#   left_join(flusht_report, by = "HHIDPN") %>%
#   mutate(prop_flusht = n_flusht/n_flusht_wave)
# 
# summarytools::freq(flusht_wide_2_3_8$prop_flusht)
