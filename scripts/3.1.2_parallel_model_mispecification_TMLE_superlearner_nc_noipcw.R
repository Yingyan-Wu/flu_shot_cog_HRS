# TMLE with super learner
# running locally
# Yingyan Wu
# Natural course without IPCW
# running on participants without missingness till current wave

#---- Package loading + option specifications ----
library(tidyverse)
library(lmtp)
library(SuperLearner) # Super Learner
library(future) # Parallel folds
library(future.apply) # Parallel apply
# Super learner related packages
library(nnls)
library(randomForest)
library(earth)

#---- Load the data & function ----
# # Paths local
data_folder <- here::here("data", "analysis_data//")
# load the data
load(paste0(data_folder, "data_for_cluster_20250824.RData"))

#---- Scenario tibble ----
n_imp = 10
waves_fu = 5:15

scenario_tib <- expand_grid(imp = 1:n_imp,
                            waves = waves_fu,
                            interventions = "nc")

#---- Natural course function ----
nc_parallel_func <- function(scenario_num){
  # Parallel set up (similar to cluster setting)
  scenario <- scenario_tib[scenario_num, ]
  wave <- scenario$waves
  i <- which(waves_fu == wave)
  intervention <- scenario$interventions
  num_imp <- scenario$imp
  outcome_var <- paste0("r", wave, "memimp_stdz")
  intervention <- scenario$interventions
  num_imp <- scenario$imp
  outcome_var <- paste0("r", wave, "memimp_stdz")
  
  # Filter on nonmissing outcome at each wave
  dt_imp_wide_imp1 <- hrs_imp_r5resp_wide %>% filter(imp == num_imp,
                                                     !is.na(get(outcome_var)))
  nrow(dt_imp_wide_imp1)
  
  lmtp_nc_temp <- lmtp_tmle(dt_imp_wide_imp1,
                            trt = trt_allwaves[[i]],
                            outcome = paste0("r", waves_fu[i], "memimp_stdz"),
                            baseline = bl_vars,
                            time_vary = L_allwaves[[i]],
                            outcome_type = "continuous",
                            mtp = FALSE, 
                            weights = dt_imp_wide_imp1$r5wtresp,
                            folds = 1)
  
  results <- c(lmtp_nc_temp$theta, lmtp_nc_temp$standard_error) %>%
    t() %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("point_est", "se")) %>%
    mutate(type = "nc_noipcw",
           wave = wave,
           imp = num_imp)
  return(results)
}

#---- Parallel set up ----
(start_nc <- Sys.time())
plan(multisession, workers = 10)
lmtp_nc_noipcw <- future_Map(function(i) {
  nc_parallel_func(scenario_num = i)
}, 
i = 1:nrow(scenario_tib), future.seed = TRUE)

plan(sequential)
(end_nc <- Sys.time())
end_nc - start_nc
# 26.02432 mins parallel on 10 cores

nc_results <- do.call(rbind, lmtp_nc_noipcw)

#---- save the results ----
save(lmtp_nc_noipcw, file = here::here("data", "model_results", 
                                       "lmtp_nc_noipcw_model_output_2000_ot52l90.RData"))
save(nc_results, file = here::here("data", "model_results", 
                                   "lmtp_nc_noipcw_results_2000_ot52l90.RData"))
