# Multiple Imputation
# Created by: Yingyan Wu
# Dec.25.2022

#---- Package loading + preparation ----
rm(list = ls())

# Packages
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}
p_load("haven", "magrittr", "tidyverse", "mice", "miceadds", "DataExplorer")

#---- function ----
source(here::here("scripts", "functions.R"))

#---- Load the dataset ----
hrs_selected <- readRDS(here::here("data", "analysis_data", "hrs_selected.rds")) %>%
  labelled::remove_val_labels()

#---- Pre imp ----
colnames(hrs_selected)
# Baseline: wave 3, yr 1996
baseline_vars <- c("r3age_v2", "female", "raceeth", "raedyrs")
tv_vars <- c("age_v2", "mstat", "flusht",
             # income pp
             "itot", "hhres",
             "insurance_2cat", "smoken", "cesd", "flone",
             # Health conditions
             "hibpe", "diabe", "cancre", "stroke", "hearte", "lunge",
             # ADL and IADLs
             "adl6a", "iadl5a",
             # Survey weights
             "wtresp",
             # Follow up status
             "inw", "dead", "iwstat",
             # Outcome
             "memimp")

flusht_waves <- seq(3, 15, by = 2)

hrs_forimp <- hrs_selected %>%
  select(-contains(c("adl_2cat", "iadl_2cat", "mstat_2cat"))) %>% # drop ADL and IADL derived variables
  select(HHIDPN, all_of(baseline_vars),
         # only select waves for flushot questions for tv_vars
         paste0("r", as.vector(outer(flusht_waves, tv_vars, "paste0"))))

hrs_forimp_long <- hrs_forimp %>%
  pivot_longer(cols = contains(tv_vars),
               names_to = c("wave", ".value"), 
               names_pattern = "r(\\d+)(.*)") %>%
  left_join(hrs_forimp %>% select(HHIDPN, r3age_v2), by = "HHIDPN") %>%
  mutate(wave = as.numeric(wave),
         fu_year = age_v2 - r3age_v2)

# # # Check iwstat
# hrs_forimp_long %>%
#   filter(!iwstat %in% c(1, 4, 5:7)) %>% view()
# hrs_forimp_long %>%
#   filter(iwstat %in% 5:7) %>% view()
# hrs_forimp_long %>%
#   filter(iwstat == 9) %>% view()

#---- Check the missingness ----
impute_data <- hrs_forimp_long %>%
  # only impute for observations when participants are alive and eligible in the sample
  filter(iwstat %in% c(1, 4)) %>%
  select(-dead, -inw, -age_v2, -fu_year)
DataExplorer::plot_missing(impute_data %>% select(-HHIDPN))
missing_summary <- DataExplorer::profile_missing(impute_data %>% 
                                                   select(-HHIDPN)) %>%
  mutate(pct_missing = pct_missing*100) %>%
  arrange(pct_missing) %>%
  print(n = Inf)

writexl::write_xlsx(missing_summary, 
                    here::here("output", "tables", "missing_summary.xlsx"))

ordered_vars <- paste0(missing_summary$feature)

baseline_vars_w_NA <- missing_summary %>% 
  filter(feature %in% baseline_vars, pct_missing != 0) %>% 
  pull(feature) %>% as.character()
tv_vars_w_NA <- missing_summary %>% 
  filter(feature %in% tv_vars, pct_missing != 0) %>% 
  pull(feature) %>% as.character()

#---- Data prep ----
# Check class for all variables
impute_data %>% map_dfr(class) %>% t()

# Variable types
variable_length_mat <- impute_data %>%
  select(-HHIDPN) %>%
  map_dfr(function(x) length(unique(x))) %>%
  pivot_longer(cols = everything(),
               names_to = "Var",
               values_to = "length")
table(variable_length_mat$length, useNA = "ifany")
# view(variable_length_mat %>% filter(length <= 30 & !Var %in% binary_vars))

binary_vars <- 
  variable_length_mat[variable_length_mat$length %in% c(2, 3), "Var"][["Var"]]

ordinal_vars <- c("adl6a", "iadl5a", "cesd")

cat_vars <- c("raceeth", "mstat")

(cont_vars <- variable_length_mat %>%
    filter(!Var %in% c(binary_vars, ordinal_vars, cat_vars)) %>% pull(Var))

# Check
length(c(cont_vars, binary_vars, ordinal_vars, cat_vars))
ncol(impute_data)
sum(duplicated(c(cont_vars, binary_vars, ordinal_vars, cat_vars)))

# Check types
impute_data %>%
  map_dfc(function(x) class(x)[1]) %>%
  t() %>% print()

# Set classes by type
impute_data %<>%
  mutate(across(all_of(binary_vars), \(x) as.factor(x)),
         # # comment next line out to treat ordinal vars as continuous (it seems multilevel MI does not like ordered vars)
         # across(all_of(ordinal_vars), \(x) factor(x, ordered = T)),
         across(all_of(cat_vars), \(x) as.factor(x)),
         across(all_of(cont_vars), \(x) as.numeric(x))) %>%
  #Arrange by order
  select(HHIDPN, all_of(ordered_vars))

#---- Multilevel imputation ----
##---- initial predictor matrix and imputation methods ----
predM <- mice::make.predictorMatrix(data = impute_data) 
impM <- mice::make.method(data = impute_data) 

predM1 <- predM
predM1[, "HHIDPN"] <- 0 # do not predict anything with HHIDPN
predM1["HHIDPN", ] <- 0 # this step is optional
predM1[tv_vars_w_NA, "HHIDPN"] <- -2 # set nesting variable for the tv vars

impM1 <- impM
# 2l methods are for time-varying variables
# 2l.only methods are for subject-level variables
# 2lonly.function is a generic method; we specify the exact methods to use later on
impM1[tv_vars_w_NA] <- "2l.pmm" 
impM1[c(baseline_vars_w_NA)] <- "2lonly.function" 

# define imputation functions
# use pmm for all 2lonly methods
imputationFunction <- rep("pmm", length(baseline_vars_w_NA)) %>% as.list()
names(imputationFunction) <- c(baseline_vars_w_NA)
imputationFunction

# define cluster variable, i.e. nested in HHIDPN
cluster_var <- rep("HHIDPN", length(baseline_vars_w_NA)) %>% as.list()
names(cluster_var) <- c(baseline_vars_w_NA)
cluster_var

##---- Run imputation ----
tictoc::tic()
# pending on change maxit to 10 and m = 10
imp_pmm_10 <- mice::mice(impute_data, m = 10, maxit = 10, 
                         predictorMatrix = predM1,
                         method = impM1,
                         imputationFunction = imputationFunction,
                         cluster_var = cluster_var,
                         seed = 62283)
tictoc::toc()
# m=1, maxit=1: 3.74min
# m = 2, maxit = 1: 5.93 min
# m = 10, maxit = 10: 6.22h
save(imp_pmm_10, file = here::here("data", "analysis_data", "imp_pmm_10.RData"))

load(here::here("data", "analysis_data", "imp_pmm_10.RData"))
imp_pmm_10$loggedEvents %>% view() # dropped iwstat for a bunch of tv vars
plot(imp_pmm_10)
# # Look at density plots by individual variables
# cont_vars
# densityplot(imp_pmm_10, ~raedyrs) # very similar distributions
# densityplot(imp_pmm_10, ~itot) 
# densityplot(imp_pmm_10, ~hhres) 
# densityplot(imp_pmm_10, ~wtresp) 
# densityplot(imp_pmm_10, ~memimp)

#---- Stack imputed dataset ----
load(here::here("data", "analysis_data", "imp_pmm_10.RData"))
imp_num = 10
imp_temp <- list()
for (i in 1:imp_num){
  imp_temp[[i]] <- complete(get(paste0("imp_pmm_", imp_num)), action = i) %>% 
    mutate(iwstat = as.numeric(as.character(iwstat))) %>%
    # get the nonmissing variables (wave, inw, dead and age_v2) from
    # pre imputation long dataset 
    left_join(hrs_forimp_long %>%
                select(HHIDPN, fu_year, wave, iwstat, inw, dead, age_v2),
              by = c("HHIDPN", "fu_year", "iwstat")) %>%
    # also add back iwstat 5:7 and 9 ~ dead, drop sample or unknown
    rbind(., hrs_forimp_long %>% filter(!iwstat %in% c(1, 4)))
  
  imp_temp[[i]][, "imp"] <- i
}

imp_pmm_stacked <- do.call(rbind, imp_temp)

# Sanity check
colnames(imp_pmm_stacked)
with(imp_pmm_stacked, table(imp, useNA = "ifany"))
with(imp_pmm_stacked, table(iwstat, useNA = "ifany"))
# with(imp_pmm_stacked %>% 
#   mutate(test_id = paste0(HHIDPN, "_", imp, "_", fu_year)) %>%
#   group_by(test_id) %>%
#   count(), table(n))

imp_pmm_stacked %>%
  filter(iwstat %in% c(1, 4)) %>%
  select(all_of(binary_vars)) %>%
  apply(., 2, function(x) table(x, useNA = "ifany"))
imp_pmm_stacked %>%
  filter(iwstat %in% c(1, 4)) %>%
  select(all_of(cat_vars)) %>%
  apply(., 2, function(x) table(x, useNA = "ifany"))
imp_pmm_stacked %>%
  select(all_of(ordinal_vars)) %>%
  apply(., 2, function(x) table(x, useNA = "ifany"))
imp_pmm_stacked %>%
  filter(iwstat %in% c(1, 4)) %>%
  select(all_of(cont_vars)) %>%
  apply(., 2, summary)
impute_data %>%
  filter(iwstat %in% c(1, 4)) %>%
  select(all_of(cont_vars)) %>%
  apply(., 2, summary)
imp_pmm_stacked %>%
  map_dfr(class) %>% t()

#---- rederived the analytic variables ----
imp_pmm_stacked %<>%
  mutate(income_pp = itot/sqrt(hhres),
         mstat_2cat = case_when(mstat %in% c(1, 2, 3) ~ 1, 
                                # Separated/divorced/widowed/not married/Widowed
                                mstat %in% c(4, 5, 6, 7, 8) ~ 0),
         adl_2cat = case_when(adl6a == 0 ~ 0,
                              adl6a > 0 ~ 1),
         iadl_2cat = case_when(iadl5a == 0 ~ 0,
                               iadl5a > 0 ~ 1)) %>%
  mutate(across(all_of(c(binary_vars, cat_vars)), \(x) as.numeric(as.character(x))))

#---- Save imputed data ----
saveRDS(imp_pmm_stacked, here::here("data", "analysis_data", "imp_pmm_stacked.RDS"))
