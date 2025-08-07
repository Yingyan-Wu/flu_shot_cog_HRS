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
data_folder <- here::here("data", "analysis_data/")

# # Paths cluster
# project_folder <- "/u/home/y/yingyanw/flu_cog_HRS/"
# data_folder <- paste0(project_folder, "data/")
# results_folder <- paste0(project_folder, "TMLE_SuperLearner/")

# load the data
load(paste0(data_folder, "data_for_cluster_20250502.RData"))

#---- Scenario tibble ----
n_imp = 10
waves_fu = 3:15
# Specify base learners for Super Learning
# Review available models. 
# listWrappers()
## GLM, random forest, and MARS (multivariate adaptive regression splines)
sl_lib <- c("SL.glm", "SL.randomForest", "SL.earth")

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
  dt_imp_wide_imp1 <- hrs_imp_r3resp_wide %>% filter(imp == num_imp,
                                                     !is.na(get(outcome_var)))
  nrow(dt_imp_wide_imp1)
  
  lmtp_nc_temp <- lmtp_tmle(dt_imp_wide_imp1,
                            trt = trt_allwaves[[i]],
                            outcome = paste0("r", waves_fu[i], "memimp_stdz"),
                            baseline = bl_vars,
                            time_vary = L_allwaves[[i]],
                            outcome_type = "continuous",
                            mtp = FALSE, 
                            weights = dt_imp_wide_imp1$r3wtresp,
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
test <- future_Map(function(i) {
  nc_parallel_func(scenario_num = i)
}, 
i = 1:nrow(scenario_tib), future.seed = TRUE)

plan(sequential)
(end_nc <- Sys.time())
end_nc - start_nc

nc_results <- do.call(rbind, test)

#---- save the results ----
save(nc_results, file = here::here("data", "model_results", "lmtp_nc_noipcw_results.RData"))
