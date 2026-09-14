# TMLE with super learner (r5age <= 80)
# running on hoffman cluster
# Yingyan Wu
# Same analysis as in Script 3.1 (with different dataset, and drop 0sus intervention)

#---- Argument ----
# Read in the arguments listed in the:
# R CMD BATCH --no-save --no-restore "--args scenario_num=$SGE_TASK_ID"  
## expression:
args = (commandArgs(TRUE))

# Check to see if arguments are passed and set default values if not.
# If so, parse the arguments. (I only have one argument here.)
if (length(args) == 0) {
  print("No arguments supplied.")
  ##supply default values
  scenario_num <- 1
} else {
  for (i in 1:length(args)) {
    eval(parse(text = args[[i]])) # parse each argument (only have 1 here)
  }
}

# print values just to make sure:
print(scenario_num)
args

#---- Package loading + option specifications ----
library(tidyverse)
library(lmtp)
library(SuperLearner) # Super Learner
library(future) # Parallel folds
# Super learner related packages
library(nnls)
library(randomForest)
library(earth)

#---- Load the data & function ----
# # Paths local
# data_folder <- here::here("data", "analysis_data/")

# Paths cluster
project_folder <- "/u/home/y/yingyanw/flu_cog_HRS/"
data_folder <- paste0(project_folder, "data/")
results_folder <- paste0(project_folder, "TMLE_SuperLearner_ager5le80/")

# load the data
load(paste0(data_folder, "data_for_cluster_ager5le80_2000_2022.RData"))

#---- Scenario tibble ----
n_imp = 10
waves_fu = 5:16
# Specify base learners for Super Learning
# Review available models. 
# listWrappers()
## GLM, random forest, and MARS (multivariate adaptive regression splines)
sl_lib <- c("SL.glm", "SL.randomForest", "SL.earth")

scenario_tib <- expand_grid(
  imp = 1:n_imp,
  waves = waves_fu,
  interventions = "1sus") %>% # "0sus")) %>%
  add_row(imp = 1:n_imp,
          waves = NA,
          interventions = "nc") %>%
# faster things run first
arrange(desc(interventions), waves)

#---- Scenario variables ----
# Make sure that scenario_num is not greater than all the scenario index
if (scenario_num > nrow(scenario_tib)) {
  print(paste0("out of bound senario number", scenario_num))
  # supply default value
  scenario_num <- 1
}

# print values just to make sure:
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

dt_imp_wide_imp1 <- hrs_imp_r5resp_wide %>% filter(imp == num_imp) # 17655 obs
nrow(dt_imp_wide_imp1)

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
  plan(multisession, workers = 5)
  set.seed(num_imp)
  (trt <- case_when(intervention == "1sus" ~ 1,
                    intervention == "0sus" ~ 0))
  temp <- lmtp_sus(i, treated = trt)
  plan(sequential)
  
  assign(paste0("y_w", wave, "_", intervention, "_", num_imp), temp)
  (end_sus <- Sys.time())
  end_sus - start_sus
}

#---- Save the lmtp object ----
save(list = case_when(intervention == "nc" ~ paste0("lmtp_nc_results_imp_", num_imp),
                      TRUE ~ paste0("y_w", wave, "_", intervention, "_", num_imp)), 
     file = paste0(results_folder, "tmle_SL_scenario_", scenario_num, ".RData"))
