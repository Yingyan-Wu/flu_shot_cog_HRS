# Create Rdata for superlearner running on cluster
# local script
# Yingyan Wu

#---- Load packages ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "haven", "labelled",
       "survey", "geepack", "splines", "lmtp")
options(scipen = 999)

#---- Load the data ----
hrs_imp <- readRDS(here::here("data", "analysis_data", "imp_pmm_stacked_outcome_2000_2022_ot52l90.RDS"))

#---- Preparing the data ----
# Sanity check
with(hrs_imp, table(wave, imp, useNA = "ifany"))

bl_vars <- c("r5age_v2", "female", "raceeth", "raedyrs", "rabyear")
tv_vars <- c("age_v2", 
             "income_pp", "mstat_2cat", 
             "insurance_2cat", "doct_vis_mt5",
             # Health behaviors
             "smoken", "mod_drinking", "vigorous_pa",
             "cesd", # "flone",
             # Health conditions
             "shlt_fair_poor",
             "hibpe", "diabe", "cancre", "stroke", "hearte", "lunge",
             # ADL and IADLs
             "adl_2cat", "iadl_2cat",
             # exposure
             "flusht", 
             # Follow up status #"inw", "dead",
             "iwstat", "inw_imp", "dead",
             # # Survey weights
             # "wtresp",
             # Outcome (at previous waves)
             "memimp_stdz")

hrs_imp_r5resp_wide <- hrs_imp %>%
  select(HHIDPN, wave, imp, all_of(c(bl_vars, tv_vars)), -r5age_v2, r5wtresp) %>%
  mutate(
    # for those who died this wave, set the flusht value to be 0 for this wave
    # These people were not used anyway as the outcome at this wave is still missing for these people
    # this is just for the lmtp_tmle function to be able to run
    flusht = case_when(iwstat %in% c(5, 6, 7, 9) ~ 0,
                       TRUE ~ flusht)) %>%
  pivot_wider(names_from = "wave", 
              values_from = c(all_of(tv_vars)),
              names_glue = "r{wave}{.value}") %>%
  select(-paste0("r", seq(4, 14, by = 2), "flusht"), -contains("r2memimp"))

# Descriptives:
colnames(hrs_imp_r5resp_wide)
for (w in 5:16){
  print(with(hrs_imp_r5resp_wide, table(!is.na(get(paste0("r", w, "memimp_stdz"))), 
                                        get(paste0("r", w, "inw_imp")), 
                                        useNA = "ifany")))
}
hrs_imp_r5resp_wide %>%
  select(contains("memimp")) %>%
  apply(., 2, summary)
hrs_imp_r5resp_wide %>%
  select(contains("memimp")) %>%
  apply(., 2, sd, na.rm = T)
hrs_imp_r5resp_wide %>%
  select(contains("flusht")) %>%
  apply(., 2, summary)
with(hrs_imp_r5resp_wide %>%
       mutate(num_flusht = rowSums(across(contains("flusht")))), summary(num_flusht))
with(hrs_imp_r5resp_wide %>%
       mutate(num_flusht = rowSums(across(contains("flusht")))), sd(num_flusht))
# hrs_imp_r5resp_wide %>%
#   filter(is.na(r5flusht)) %>%
#   select(HHIDPN, contains("iwstat"))

# For Flusht at waves 4, 6, 8, ..., 14, make them the same as the flusht in the previous wave
# e.g. Y4|do(X3), do(X4) = Y4|do(X3), and thus, the value of X4 does not make a difference
flusht_waves <- seq(5, 16, by = 2)
for (w in flusht_waves[flusht_waves != 16]){
  hrs_imp_r5resp_wide %<>%
    mutate(!!paste0("r", w + 1, "flusht") := !!sym(paste0("r", w, "flusht")))
}
# Sanity check
hrs_imp_r5resp_wide %>%
  select(paste0("r", 5:16, "flusht")) %>%
  apply(., 2, table, useNA = "ifany")

# Since flushot is asked about past flu shot status
# covariates were from previous waves, treatment were from the same wave with outcome

##---- Covariates ----
tv_cov_memimp <- tv_vars[!tv_vars%in% c("flusht", # exposure
                                        "fu_year", "wtresp",
                                        # "memimp_stdz",
                                        "inw", "iwstat", "dead", "inw_imp")]

# Time varying covariates list for each wave (alternative waves from 4, 5:16, 6 waves in total)
waves <- 4:16
# flusht_waves <- seq(3, 16, by = 2)

# Use lapply to create lists for each wave
L_allwaves <- lapply(waves, function(w) {
  # Consider only waves before the current one
  previous_waves <- waves[waves < w]
  lapply(previous_waves, function(w) 
    case_when(w == 4 ~ 
                list(paste0("r", w, tv_cov_memimp[tv_cov_memimp != "memimp_stdz"])),
              TRUE ~ list(paste0("r", w, tv_cov_memimp)))[[1]])
})
L_allwaves[[1]] <- NULL
# 
# # Assign separate list objects to individual variables 
# # e.g. L_r5 = r5 tv vars, r5 tv vars
# for (i in 1:(length(lists) - 1)) {
#   assign(paste0("L_r", waves[i]), lists[[i + 1]])
# }

##---- Censoring ----
# 1 if in the study at this wave, 0 otherwise
waves_fu <- waves[waves != 4]
C_allwaves <- lapply(waves_fu, function(w) {
  # Consider only waves before the current one
  uptonow_waves <- waves_fu[waves_fu <= w]
  paste0("r", uptonow_waves, "inw_imp")
})

# # Assign separate list objects to individual variables
# for (i in seq_along(lists)) {
#   assign(paste0("C_r", waves_fu[i]), lists[[i]])
# }

##---- Treatment ----
# Since flushot is asked about past flu shot status
# covariates were from previous waves, treatment were from the same wave with outcome
# Treatment variables r5 outcome ~ r5 r5 trt
# r6 outcome ~ r5 r5 trt
# r7 outcome ~ r5 trt, r5 trt, r7 trt etc.
waves_fu <- waves[waves != 4]
trt_allwaves <- lapply(waves_fu, function(w) {
  # Consider waves on and before the current one
  uptonow_waves <- waves_fu[waves_fu <= w]
  paste0("r", uptonow_waves, "flusht")
})

##---- check flu shot distribution ----
hrs_imp_r5resp_wide %>% select(contains("flusht"), imp) %>% 
  apply(., 2, function(x) table(x, .$imp, useNA = "ifany"))
# Filter to those who got flu shot every wave
hrs_imp_r5resp_wide %>% filter(if_all(contains("flusht"), ~ .x == 1)) %>%
  with(., table(imp, useNA = "ifany"))
hrs_imp_r5resp_wide %>% filter(if_all(contains("flusht"), ~ .x == 0)) %>%
  with(., table(imp, useNA = "ifany"))

#---- Sens1: restrict to one person in the household ----
# using lowest number of pn
hrs_imp_r5resp_wide_minpn <- hrs_imp_r5resp_wide %>%
  left_join(hrs_imp %>% select(HHIDPN, household_min_pn) %>% distinct(),
            by = "HHIDPN") %>%
  filter(household_min_pn == 1)
# Sanity check
with(hrs_imp_r5resp_wide_minpn, table(imp, useNA = "ifany"))
# n=12908
sum(hrs_imp_r5resp_wide_minpn$r5wtresp)/10 # 48531932
mean(hrs_imp_r5resp_wide_minpn$r5age_v2) # 68.15707

#---- Sens: restrict baseline age to 80 ----
hrs_imp_r5resp_wide_ager5le80 <- hrs_imp_r5resp_wide %>%
  filter(r5age_v2 <= 80) 

# Sanity check
with(hrs_imp_r5resp_wide_ager5le80, table(imp, useNA = "ifany"))
# Descriptives
#n=17655
sum(hrs_imp_r5resp_wide_ager5le80$r5wtresp)/10 # 63964473
mean(hrs_imp_r5resp_wide_ager5le80$r5age_v2) # 65.39377

#---- Save the elements for TMLE ----
# restrict to age [52, 90) at baseline, add 2022 variables
save(hrs_imp_r5resp_wide, bl_vars, L_allwaves, trt_allwaves, C_allwaves, 
     file = here::here("data", "analysis_data", "data_for_cluster_2000_2022.RData"))

# Save sens r5age <= 80 data with non-specific dataset name
# r5age [52, 80]
hrs_imp_r5resp_wide <- hrs_imp_r5resp_wide_ager5le80
save(hrs_imp_r5resp_wide, bl_vars, L_allwaves, trt_allwaves, C_allwaves, 
     file = here::here("data", "analysis_data", "data_for_cluster_ager5le80_2000_2022.RData"))

# Save sens filter to min person number only with non-specific dataset name
hrs_imp_r5resp_wide <- hrs_imp_r5resp_wide_minpn
save(hrs_imp_r5resp_wide, bl_vars, L_allwaves, trt_allwaves, C_allwaves, 
     file = here::here("data", "analysis_data", "data_for_cluster_minpn_2000_2022.RData"))

