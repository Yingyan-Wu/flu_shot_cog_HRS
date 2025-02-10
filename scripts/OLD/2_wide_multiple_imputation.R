# Multiple Imputation (wide dataset)
# Created by: Yingyan Wu
# Dec.25.2022

# Cannot run

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
tv_cov <- tv_vars[!tv_vars %in% c("flusht", "memimp")]
flusht_waves <- seq(3, 15, by = 2)

hrs_forimp <- hrs_selected %>%
  select(-contains(c("adl_2cat", "iadl_2cat", "mstat_2cat"))) %>% # drop ADL and IADL derived variables
  select(HHIDPN, all_of(baseline_vars),
         # only select waves for flushot questions for tv_vars
         # and 1 wave before flushot questions for the first wave of flushot
         paste0("r", 2, tv_cov),
         paste0("r", as.vector(outer(flusht_waves, tv_vars, "paste0"))))

hrs_forimp  %>%
  select(contains("iwstat")) %>%
  apply(., 2, table, useNA = "ifany")

#---- Check the missingness ----
with(hrs_forimp, table(r3iwstat, r2iwstat, useNA = "ifany"))
impute_data <- hrs_forimp %>%
  # filter on iwstat 0, 1, 4 for r3 # eligible for wave 3 or joined in later waves
  filter(r3iwstat %in% c(0, 1, 4)) %>%
  select(-contains(c("dead", "inw")))
# DataExplorer::plot_missing(impute_data %>% select(-HHIDPN))
missing_summary <- DataExplorer::profile_missing(impute_data %>% 
                                                   select(-HHIDPN)) %>%
  mutate(pct_missing = pct_missing*100) %>%
  arrange(pct_missing) %>%
  print(n = Inf)

ordered_vars <- paste0(missing_summary$feature)

#---- Data prep ----
# Check class for all variables
impute_data <- impute_data[, ordered_vars] %>%
  as.data.frame() %>% print()
ncol(impute_data)
colnames(impute_data)
str(impute_data)

# Variable types
variable_length_mat <- impute_data %>%
  map_dfr(function(x) length(unique(x))) %>%
  pivot_longer(cols = everything(),
               names_to = "Var",
               values_to = "length")
table(variable_length_mat$length, useNA = "ifany")
# view(variable_length_mat %>% filter(length <= 30 & !Var %in% binary_vars))

binary_vars <- 
  variable_length_mat[variable_length_mat$length %in% c(2, 3), "Var"][["Var"]]

ordinal_vars <- variable_length_mat %>% filter(str_detect(Var, "adl6a|iadl5a|cesd")) %>%
  pull(Var)

cat_vars <- variable_length_mat %>% filter(str_detect(Var, "raceeth|mstat")) %>%
  pull(Var)

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
         across(all_of(ordinal_vars), \(x) factor(x, ordered = T)),
         across(all_of(cat_vars), \(x) as.factor(x)),
         across(all_of(cont_vars), \(x) as.numeric(x))) %>%
  #Arrange by order
  select(all_of(ordered_vars))

#---- Multilevel imputation ----
#---- Initiate imputation (PMM) ----
ini_imp <- mice(impute_data, maxit = 0,
                defaultMethod = c("pmm", "pmm", "pmm", "pmm"),
                seed = 62283)
# ini_imp$loggedEvents
(meth <- ini_imp$method)
(pred <- ini_imp$predictorMatrix)

##---- Run imputation ----
tictoc::tic()
# PENDING on change it 10
imp_pmm_10 <- mice::mice(impute_data, m = 2, maxit = 1,
                         pred = pred, meth = meth,
                         defaultMethod = c("pmm", "pmm", "pmm", "pmm"),
                         seed = 62283)
# Can not run. This is huge! 
# Error in solve.default(xtx + diag(pen)) : 
#   system is computationally singular: reciprocal condition number = 1.8764e-16
tictoc::toc()
# m = 2, maxit = 1 min
save(imp_pmm_10, file = here::here("data", "analysis_data", "imp_pmm_wide.RData"))