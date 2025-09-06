# TMLE with super learner
# running on hoffman cluster
# Yingyan Wu

#---- Package loading + option specifications ----
library(tidyverse)
library(lmtp)
library(SuperLearner) # Super Learner
library(future) # Parallel folds
library(future.apply)
# Super learner related packages
library(nnls)
library(randomForest)
library(earth)

#---- Load the data & function ----
# # Paths local
data_folder <- here::here("data", "analysis_data//")
results_folder <- here::here("data", "model_results", "local", "TMLE_SuperLearner//")

# Paths cluster
# project_folder <- "/u/home/y/yingyanw/flu_cog_HRS/"
# data_folder <- paste0(project_folder, "data/")
# results_folder <- paste0(project_folder, "TMLE_SuperLearner/")

# load the data
load(paste0(data_folder, "data_for_cluster_20250824.RData"))

#---- Scenario tibble ----
n_imp = 10
waves_fu = 5:15

# Specify base learners for Super Learning
# Review available models. 
# listWrappers()
## GLM, random forest, and MARS (multivariate adaptive regression splines)
sl_lib <- c("SL.glm", "SL.randomForest", "SL.earth")

scenario_tib <- expand_grid(
  imp = 1:n_imp,
  waves = waves_fu,
  interventions = "1sus") %>%
  filter(waves == 15)

#---- LMTP function ----
lmtp_1sus_func <- function(intervention = "1sus"){
  
  #---- Analysis ----
  if (intervention == "nc"){
    # Natural course # Around 1 h for each imputation
    (start_nc <- Sys.time())
    lmtp_nc <-  function(i) {
      lmtp_tmle(dt_imp_wide_imp1,
                trt = trt_allwaves[[i]],
                outcome = paste0("r", waves_fu[i], "memimp_stdz"),
                baseline = bl_vars,
                time_vary = L_allwaves[[i]],
                cens = C_allwaves[[i]],
                outcome_type = "continuous",
                mtp = FALSE, 
                weights = dt_imp_wide_imp1$r5wtresp,
                folds = 1)
    }
    temp <- lapply(seq_along(waves_fu), lmtp_nc)
    (end_nc <- Sys.time())
    end_nc - start_nc
    assign(paste0("lmtp_nc_results_imp_", num_imp), temp)
  }
  
  # "Sustained" with super learner
  if (str_detect(intervention, "sus")){
    lmtp_sus <-  function(i, treated = 1) {
      lmtp_tmle(dt_imp_wide_imp1,
                trt = trt_allwaves[[i]],
                outcome = paste0("r", waves_fu[i], "memimp_stdz"),
                baseline = bl_vars,
                time_vary = L_allwaves[[i]],
                shift = case_when(treated == 1 ~ list(static_binary_on),# Turn all treatment nodes on
                                  treated == 0 ~ list(static_binary_off))[[1]], # Turn all treatment nodes off
                cens = C_allwaves[[i]],
                outcome_type = "continuous",
                mtp = TRUE, 
                weights = dt_imp_wide_imp1$r5wtresp,
                folds = 5,
                learners_outcome = sl_lib,
                learners_trt = sl_lib)
    }
    (start_sus <- Sys.time())
    plan(list(tweak(multisession, workers = 5)))
    set.seed(num_imp)
    (trt <- case_when(intervention == "1sus" ~ 1,
                      intervention == "0sus" ~ 0))
    temp <- lmtp_sus(i, treated = trt)
    plan(sequential)
    
    assign(paste0("y_w", wave, "_", intervention, "_", num_imp), temp)
    (end_sus <- Sys.time())
    print(end_sus - start_sus)
    
    
  }
  save(list = case_when(intervention == "nc" ~ paste0("lmtp_nc_results_imp_", num_imp),
                        TRUE ~ paste0("y_w", wave, "_", intervention, "_", num_imp)), 
       file = paste0(results_folder, "tmle_SL_scenario_", scenario_num, ".RData"))
  
  
}


#---- Parallel ----
options(future.globals.maxSize = 8000*1024^2)

# Half and Half
scenes <- 4:5
# scenes <- 6:10
for (scenario_num in scenes){
  
  print(scenario_num)
  scenario <- scenario_tib[scenario_num, ]
  wave <- scenario$waves
  i <- which(waves_fu == wave)
  intervention <- scenario$interventions
  num_imp <- scenario$imp
  outcome_var <- paste0("r", wave, "memimp_stdz")
  
  print(wave)
  print(i)
  print(intervention)
  print(num_imp)
  
  dt_imp_wide_imp1 <- hrs_imp_r5resp_wide %>% filter(imp == num_imp)
  nrow(dt_imp_wide_imp1) # 20091
  
  (start_1sus <- Sys.time())
  lmtp_1sus_func(intervention = "1sus")
  (end_1sus <- Sys.time())
  print(end_1sus - start_1sus)
  
}

Sys.time()
# 
# plan(list(tweak(multisession, workers = 5)))
# a <- future_Map(function(x){ rnorm(100, x, 3)}, x=1:5, future.seed = TRUE)
# plan(sequential)
