# TMLE with super learner -- positive control
# Death at the last wave.
# running on hoffman cluster
# Yingyan Wu

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
# library(haven) # Read data
library(dplyr)
library(tidyverse)
library(lmtp)
library(SuperLearner) # Super Learner
library(future) # Parallel folds
# Super learner related packages
library(nnls)
library(randomForest)
library(earth)

#---- Load the data & function ----
# Paths local
data_folder <- here::here("data", "analysis_data/")

# # Paths cluster
# project_folder <- "/u/home/y/yingyanw/flu_cog_HRS/"
# data_folder <- paste0(project_folder, "data/")
# results_folder <- paste0(project_folder, "TMLE_SuperLearner_sens_death_pos_ctrl/")

# load the data
load(paste0(data_folder, "data_for_cluster_20250502.RData"))

#---- Data prep -----
# need to trick tmle a bit as e.g. r5dead = 1, time varying covariates are missing
# for wave 5 for people died around wave 5
L_allwaves_death <- L_allwaves[]
trt_allwaves

L <- paste0()
trt <- paste0("r", 3:14, "flusht")

# test
waves_fu <- 3:15
lmtp_test <-  lmtp_tmle(dt_imp_wide_imp1,
                        trt = trt_allwaves[[3]],
                        outcome = paste0("r", waves_fu[[3]], "dead"),
                        baseline = bl_vars,
                        time_vary = L_allwaves[[2]],
                        cens = C_allwaves[[3]],
                        outcome_type = "binomial",
                        mtp = FALSE, 
                        weights = dt_imp_wide_imp1$r3wtresp,
                        folds = 1)

#---- Scenario tibble ----
n_imp = 10
# Specify base learners for Super Learning
# Review available models. 
# listWrappers()
## GLM, random forest, and MARS (multivariate adaptive regression splines)
sl_lib <- c("SL.glm", "SL.randomForest", "SL.earth")

scenario_tib <- expand_grid(
  imp = 1:n_imp,
  interventions = "1sus") %>%
  add_row(imp = 1:n_imp,
          interventions = "nc")

#---- Scenario variables ----
# Make sure that scenario_num is not greater than all the scenario index
if (scenario_num > nrow(scenario_tib)) {
  print(paste0("out of bound senario number", scenario_num))
  # supply default value
  scenario_num <- 1
}
# scenario_num <- 11

# print values just to make sure:
print(scenario_num)
scenario <- scenario_tib[scenario_num, ]
intervention <- scenario$interventions
num_imp <- scenario$imp
outcome_var <- "r15dead"
i = 13 # treatment cumulative to wave 15

print(intervention)
print(num_imp)

dt_imp_wide_imp1 <- hrs_imp_r3resp_wide %>% filter(imp == num_imp)
nrow(dt_imp_wide_imp1) # 19845

dt_imp_wide_imp1_alive_14 <- dt_imp_wide_imp1 %>%
  filter(r14dead == 0, r13dead == 0, r12dead == 0) # 8084 obs

lmtp_result <-  lmtp_tmle(dt_imp_wide_imp1_alive_14,
                          trt = trt_allwaves[[i]],
                          outcome = outcome_var,
                          baseline = bl_vars,
                          time_vary = L_allwaves[[i]],
                          outcome_type = "binomial",
                          mtp = FALSE, 
                          weights = dt_imp_wide_imp1_alive_14$r3wtresp,
                          folds = 1)

#---- Analysis ----
if (intervention == "nc"){
  # Natural course # Around 1 h for each imputation
  (start_nc <- Sys.time())
  lmtp_result <-  lmtp_tmle(dt_imp_wide_imp1,
                        trt = trt_allwaves[[i]],
                        outcome = outcome_var,
                        baseline = bl_vars,
                        time_vary = L_allwaves[[i]],
                        cens = C_allwaves[[i]],
                        outcome_type = "binomial",
                        mtp = FALSE, 
                        weights = dt_imp_wide_imp1$r3wtresp,
                        folds = 1)
  (end_nc <- Sys.time())
  end_nc - start_nc
}

# test
waves_fu <- 3:15
lmtp_test <-  lmtp_tmle(dt_imp_wide_imp1,
                          trt = trt_allwaves[[3]],
                          outcome = paste0("r", waves_fu[[3]], "dead"),
                          baseline = bl_vars,
                          time_vary = L_allwaves[[2]],
                          cens = C_allwaves[[3]],
                          outcome_type = "binomial",
                          mtp = FALSE, 
                          weights = dt_imp_wide_imp1$r3wtresp,
                          folds = 1)



# "Sustained" with super learner
if (str_detect(intervention, "sus")){
  (start_sus <- Sys.time())
  plan(multisession, workers = 5)
  set.seed(num_imp)
  lmtp_result <-  lmtp_tmle(dt_imp_wide_imp1,
                         trt = trt_allwaves[[i]],
                         outcome = outcome_var,
                         baseline = bl_vars,
                         time_vary = L_allwaves[[i]],
                         shift = static_binary_on,# Turn all treatment nodes on
                         cens = C_allwaves[[i]],
                         outcome_type = "binomial",
                         mtp = TRUE, 
                         weights = dt_imp_wide_imp1$r3wtresp,
                         folds = 5,
                         learners_outcome = sl_lib,
                         learners_trt = sl_lib)
  
  plan(sequential)
  (end_sus <- Sys.time())
  end_sus - start_sus
}

assign(paste0("lmtp_death_", intervention, "_results_imp_", num_imp), lmtp_result)

#---- Save the lmtp object ----
save(list = case_when(intervention == "nc" ~ 
                        paste0("lmtp_death_", intervention, "_results_imp_", num_imp),
                      TRUE ~ paste0("lmtp_death_", intervention, "_results_imp_", num_imp)), 
     file = paste0(results_folder, "tmle_SL_death_scenario_", scenario_num, ".RData"))
