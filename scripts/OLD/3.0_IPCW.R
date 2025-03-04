# Create IPCW outside the lmtp_tmle function
# Created by: Yingyan Wu
# Does not finish, figured it won't work with lmtp package anyway

#---- Load packages ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "haven", "labelled",
       "survey", "geepack", "splines", "lmtp")
options(scipen = 999)

#---- Load the data ----
# memimp not imputed
hrs_imp <- readRDS(here::here("data", "analysis_data", "imp_pmm_stacked_ahead_hrs.RDS"))

#---- Preparing the data ----
# Sanity check
with(hrs_imp, table(wave, imp, useNA = "ifany"))

bl_vars <- c("r3age_v2", "female", "raceeth", "raedyrs", "r3wtresp")
tv_vars <- c("age_v2", "mstat", "income_pp", 
             "insurance_2cat", "smoken", "cesd", "flone",
             # Health conditions
             "hibpe", "diabe", "cancre", "stroke", "hearte", "lunge",
             # ADL and IADLs
             "adl_2cat", "iadl_2cat",
             # exposure
             "flusht",
             # outcome
             "memimp",
             "iwstat")
tv_pred <- tv_vars[!tv_vars %in% c("iwstat")]

hrs_imp %<>%
  mutate( # for those who died this wave, set the flusht value to be 0 for this wave
    # These people were not used anyway as the outcome at this wave is still missing for these people
    # this is just for the lmtp_tmle function to be able to run
    flusht = case_when(iwstat %in% c(5, 6, 7, 9) ~ 0,
                       TRUE ~ flusht),
    # Create a index for censoring for those who has a missingness on memimp
    memimp_missing = case_when(is.na(memimp) ~ 1,
                               TRUE ~ 0),
    memimp_nonmissing = 1 - memimp_missing)

# Sanity check
with(hrs_imp, table(memimp_missing, dead, useNA = "ifany"))

##---- pivot wider ----
hrs_imp_wide <- hrs_imp %>%
  select(-r3age_v2, -r3memimp) %>%
  pivot_wider(id_cols = c(HHIDPN, imp),
              names_from = wave,
              values_from = c(all_of(tv_vars), "memimp_missing", "memimp_nonmissing"),
              names_glue = "r{wave}{.value}")
colnames(hrs_imp_wide)

#---- Models ----
data_wide <- hrs_imp %>% filter(imp == 1)

waves_fu <- 4:15
waves <- 3:15
outcome_censor <- "memimp_nonmissing"
pcnum_formula <- paste0("r", waves_fu, outcome_censor, " ~ 1")
pcden_formula <- paste0("r", waves_fu, outcome_censor, "~ ",  
  paste(bl_vars[!str_detect(bl_vars, "r3")], collapse = " + "), " + ",
  paste(paste0("r", waves, tv_vars), collapse = " + "))

# Warnings for IPCW models
w_pcnum <- w_pcden <- numeric()

data_long_svy <- svydesign(ids = ~HHIDPN, weights = ~r3wtresp, 
                           data = data_long)
# Numerator
tryCatch({
  pcnum_model <<- svyglm(as.formula(pcnum_formula),
                         design = data_long_svy,
                         family = quasibinomial(link = "logit"),
                         control = list(maxit = 100))
  # If no warnings, append 0s
  w_pcnum <<- c(w_pcnum, 0)
}, warning = function(w) {
  # If there's warning, append 1s
  w_pcnum <<- c(w_pcnum, 1)
  pcnum_model <<- svyglm(as.formula(pcnum_formula),
                         design = data_long_svy,
                         family = quasibinomial(link = "logit"),
                         control = list(maxit = 100))
})

# Denominator
tryCatch({
  pcden_model <<- svyglm(as.formula(pcden_formula),
                         design = data_long_svy,
                         family = quasibinomial(link = "logit"),
                         control = list(maxit = 100))
  # If no warnings, append 0s
  w_pcden <<- c(w_pcden, 0)
}, warning = function(w) {
  # If there's warning, append 1s
  w_pcden <<- c(w_pcden, 1)
  pcden_model <<- svyglm(as.formula(pcden_formula),
                         design = data_long_svy,
                         family = quasibinomial(link = "logit"),
                         control = list(maxit = 100))
})

data_long %<>%
  mutate(
    pred_pcnum = predict(pcnum_model, newdata = data_long, type = "response"),
    pred_pcden = predict(pcden_model, newdata = data_long, type = "response")) %>%
  group_by(newid) %>% 
  mutate(
    cum_pred_pcnum = cumprod(pred_pcnum), 
    cum_pred_pcden = cumprod(pred_pcden)
  ) %>% 
  ungroup()

data_long %<>%
  mutate(ipcw = 1/cum_pred_pcden,
         sipcw = cum_pred_pcnum/cum_pred_pcden,
         # truncate weights
         # remember to remove na.rm = TRUE after imputation
         ipcw_trunc = case_when(
           ipcw > quantile(ipcw, 0.99) ~ quantile(ipcw, 0.99),
           TRUE ~ ipcw), 
         sipcw_trunc = case_when(
           sipcw > quantile(sipcw, 0.99) ~ quantile(sipcw, 0.99),
           TRUE ~ sipcw))

data_long %<>%
  mutate(svy_ipcw = sipcw_trunc * get(svy_wt))

data_long %<>%
  filter(!is.na(!!sym(y)))
##---- IPCW ----
# since memimp not imputed and we are creating the ipcw outside, 
# the non censoring index = memimp_nonmissing
generate_IPCW <- function(longdata, formula1) {
  # testing arguments
  # longdata <- hrs_imp1
  # formula1 <- censor_formula
  original_vars <- names(longdata)
  
  ## end of membership censoring ----
  # denominator
  # plrFitC <- glm(as.formula(paste0(formula2, collapse = "")), 
  #                data = longdata, family = binomial())
  
  tryCatch(
    {
      plrFitC <<- glm(as.formula(paste0(formula1, collapse = "")), 
                      data = longdata, family = binomial())
      w_endofMem_d <<- c(w_endofMem_d, 0)
    },
    warning = function(w){
      w_endofMem_d <<- c(w_endofMem_d, 1)
      plrFitC <<- glm(as.formula(paste0(formula1, collapse = "")), 
                      data = longdata, family = binomial())
    })
  
  # numerator
  # plrFitCnum <- glm(as.formula(formula2[1]), 
  #                   data = longdata, family = binomial())
  
  tryCatch(
    {
      plrFitCnum <<- glm(as.formula(formula2[1]), 
                         data = longdata, family = binomial())
      w_endofMem_n <<- c(w_endofMem_n, 0)
    },
    warning = function(w){
      w_endofMem_n <<- c(w_endofMem_n, 1)
      plrFitCnum <<- glm(as.formula(formula2[1]), 
                         data = longdata, family = binomial())
    })
  
  longdata <- longdata %>% 
    mutate(
      predC = 1 - predict(plrFitC, newdata = longdata, type = "response"),
      predC_num = 1 - predict(plrFitCnum, newdata = longdata, type = "response")
    ) %>% 
    group_by(new_id) %>% 
    mutate(
      cumPredC = cumprod(predC), 
      cumPredC_num = cumprod(predC_num)
    ) %>% 
    ungroup()
  
  ## calculate weights ----
  longdata <- longdata %>% 
    mutate(
      # censoring due to end of membership
      wt_endmem = 1 / cumPredC, # unstabilized
      wt_endmem_stab = cumPredC_num / cumPredC, # stabilized
      # censoring due to death
      wt_death = 1 / denom_cum_death,
      wt_death_stab = num_cum_death / denom_cum_death,
      
      # combined censoring weights
      wt_cens = wt_endmem * wt_death,
      wt_cens_stab = wt_endmem_stab * wt_death_stab,
      
      # truncate death weights
      # remember to remove na.rm = TRUE after imputation
      wt_death_trunc = ifelse(
        wt_death > quantile(wt_death, 0.99), 
        quantile(wt_death, 0.99), wt_death
      ), 
      wt_death_stab_trunc = ifelse(
        wt_death_stab > quantile(wt_death_stab, 0.99), 
        quantile(wt_death_stab, 0.99), wt_death_stab
      )
    )
  
  # exclude intermediate variables
  longdata <- longdata %>% select(all_of(original_vars), starts_with("wt_"))
  return(longdata)
}

# These are for r3memimp imputed
hrs_imp_bl <- hrs_imp %>%
  select(HHIDPN, imp, all_of(bl_vars), "r3memimp") %>%
  unique()
hrs_imp_r3resp_wide <- hrs_imp_bl %>%
  pivot_wider(names_from = "wave", 
              values_from = c(all_of(tv_vars), "memimp_missing", "memimp_nonmissing"),
              names_glue = "r{wave}{.value}") %>%
  select(-r3memimp, -r3age_v2), by = c("HHIDPN", "imp")) %>%
  select(-paste0("r", seq(2, 14, by = 2), "flusht"), -contains("r2memimp"))

# Descriptives:
colnames(hrs_imp_r3resp_wide)
hrs_imp_r3resp_wide %>%
  select(contains("memimp_nonmissing")) %>%
  apply(., 2, table, useNA = "ifany")
hrs_imp_r3resp_wide %>%
  select(contains("flusht")) %>%
  apply(., 2, summary)
with(hrs_imp_r3resp_wide %>%
       mutate(num_flusht = rowSums(across(contains("flusht")))), summary(num_flusht))
with(hrs_imp_r3resp_wide %>%
       mutate(num_flusht = rowSums(across(contains("flusht")))), sd(num_flusht))
# hrs_imp_r3resp_wide %>%
#   filter(is.na(r5flusht)) %>%
#   select(HHIDPN, contains("iwstat"))

# For Flusht at waves 4, 6, 8, ..., 14, make them the same as the flusht in the previous wave
# e.g. Y4|do(X3), do(X4) = Y4|do(X3), and thus, the value of X4 does not make a difference
flusht_waves <- seq(3, 15, by = 2)
for (w in flusht_waves){
  hrs_imp_r3resp_wide %<>%
    mutate(!!paste0("r", w + 1, "flusht") := !!sym(paste0("r", w, "flusht")))
}
# Sanity check
hrs_imp_r3resp_wide %>%
  select(contains("flusht")) %>%
  apply(., 2, table, useNA = "ifany")

# # These are for imputed all memimp
# hrs_imp_r3resp_wide <- hrs_imp %>%
#   select(HHIDPN, imp, wave, all_of(bl_vars), all_of(tv_vars), -r3age_v2) %>%
#   # for those who died this wave, set the flusht value to be 0 for this wave
#   # These people were not used anyway as the outcome at this wave is still missing for these people
#   # this is just for the lmtp_tmle function to be able to run
#   mutate(flusht = case_when(iwstat %in% c(5, 6, 7, 9) ~ 0,
#                             TRUE ~ flusht)) %>%
#   pivot_wider(names_from = "wave", 
#               values_from = tv_vars,
#               names_glue = "r{wave}{.value}")

# Since flushot is asked about past flu shot status
# covariates were from previous waves, treatment were from the same wave with outcome

##---- Covariates ----
tv_cov_memimp <- tv_vars[!tv_vars%in% c("flusht", # exposure
                                        "fu_year", "wtresp",
                                        # "memimp",
                                        "inw", "iwstat", "dead", "inw_imp")]

# Time varying covariates list for each wave (alternative waves from 2, 3:15, 8 waves in total)
waves <- 2:15
# flusht_waves <- seq(3, 15, by = 2)

# Use lapply to create lists for each wave
L_allwaves <- lapply(waves, function(w) {
  # Consider only waves before the current one
  previous_waves <- waves[waves < w]
  lapply(previous_waves, function(w) 
    case_when(w == 2 ~ 
                list(paste0("r", w, tv_cov_memimp[tv_cov_memimp != "memimp"])),
              TRUE ~ list(paste0("r", w, tv_cov_memimp)))[[1]])
})
L_allwaves[[1]] <- NULL
# 
# # Assign separate list objects to individual variables 
# # e.g. L_r5 = r3 tv vars, r5 tv vars
# for (i in 1:(length(lists) - 1)) {
#   assign(paste0("L_r", waves[i]), lists[[i + 1]])
# }

##---- Treatment ----
# Since flushot is asked about past flu shot status
# covariates were from previous waves, treatment were from the same wave with outcome
# Treatment variables r5 outcome ~ r3 r5 trt
# r6 outcome ~ r3 r5 trt
# r7 outcome ~ r3 trt, r5 trt, r7 trt etc.
waves_fu <- waves[waves != 2]
trt_allwaves <- lapply(waves_fu, function(w) {
  # Consider waves on and before the current one
  uptonow_waves <- waves_fu[waves_fu <= w]
  paste0("r", uptonow_waves, "flusht")
})

# # Assign separate list objects to individual variables
# for (i in 1:(length(lists))) {
#   assign(paste0("trt_r", waves[i]), lists[[i]])
# }