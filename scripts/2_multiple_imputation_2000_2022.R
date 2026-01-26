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
# AHEAD, CODA, HRS, War babies
load(here::here("data", "analysis_data", "hrs_selected_2000_2022_ot52l90.RData"))
hrs_selected %<>%
  labelled::remove_val_labels()

# Check age dist and missing memimp
hist(hrs_selected$r5age_v2)
table(hrs_selected$r5iwstat)/nrow(hrs_selected)
sum(is.na(hrs_selected$r5memimp))/nrow(hrs_selected)

hrs_selected %>%
  mutate(r5memimp_missing = case_when(is.na(r5memimp) ~ 1,
                                      !is.na(r5memimp) ~ 0)) %>%
  ggplot(aes(x = r5age_v2, color = factor(r5memimp_missing))) +
  geom_density() +
  theme_minimal()

#---- Flu shot ----
# Check flu shot status and racohbyr
with(hrs_selected, table(racohbyr, !is.na(r5flusht), useNA = "ifany"))
# CODA and WBB has flu shot status in 2000 too

# Corresponding to later analysis, impute non flushot waves flu shot data
# using last value carry forward for waves 4, 6, 8....
# since we are only interested in sustained intervention 
# For Flusht at waves 4, 6, 8, ..., 14, make them the same as the flusht in the previous wave
# e.g. Y4|do(X3=1), do(X4=1) = Y4|do(X3=1), and thus, the value of X4 does not make a difference

flusht_waves <- seq(3, 16, by = 2)
for (w in flusht_waves[flusht_waves != 16]){
  hrs_selected %<>%
    mutate(!!paste0("r", w + 1, "flusht") := !!sym(paste0("r", w, "flusht")))
}

# Sanity check
hrs_selected %>%
  select(paste0("r", 3:16, "flusht")) %>%
  apply(., 2, table, useNA = "ifany")

#---- Pre imp ----
colnames(hrs_selected)
# Baseline: wave 5, yr 2000
baseline_vars <- c("female", "raceeth", "raedyrs", # "r5memimp", impute memimp longitudinally
                   "r5wtresp")
tv_vars <- c("age_v2", "flusht",
             # income pp
             "itot", "hhres",
             "mstat", 
             "insurance_2cat", "doct_vis_mt5",
             # Health behaviors
             "smoken", "mod_drinking", "vigorous_pa",
             "cesd",
             # Health conditions
             "shlt_fair_poor",
             "hibpe", "diabe", "cancre", "stroke", "hearte", "lunge",
             # ADL and IADLs
             "adl6a", "iadl5a",
             # Outcome
             "memimp",
             # Survey weights
             "wtresp",
             # Follow up status 
             "inw", "dead", "iwstat")

# Don't impute the outcome for now
tv_vars_forimp <- tv_vars[!tv_vars %in% c("wtresp", "inw", "dead")]
tv_cov_forimp <- tv_vars_forimp[!tv_vars_forimp %in% c("flusht", "memimp")]
flusht_waves <- seq(5, 16, by = 2)
waves <- 4:16

hrs_forimp <- hrs_selected %>%
  select(-contains(c("adl_2cat", "iadl_2cat", "mstat_2cat"))) %>% # drop ADL and IADL derived variables
  select(HHIDPN, all_of(baseline_vars),
         paste0("r", as.vector(outer(waves, tv_vars, "paste0"))))

hrs_forimp_long <- hrs_forimp %>%
  pivot_longer(cols = contains(tv_vars),
               names_to = c("wave", ".value"), 
               names_pattern = "r(\\d+)(.*)") %>%
  mutate(wave = as.numeric(wave)) %>%
  left_join(hrs_selected %>% select(HHIDPN, # r5memimp, 
                                    r5wtresp), by = "HHIDPN")

# # Check iwstat
# hrs_forimp_long %>%
#   filter(!iwstat %in% c(1, 4, 5:7))
# hrs_forimp_long %>%
#   filter(iwstat %in% 5:7) %>% view()
# hrs_forimp_long %>%
#   filter(iwstat == 9)

#---- Check the missingness ----
hrs_preimp <- hrs_forimp_long %>%
  select(HHIDPN, wave, all_of(baseline_vars), all_of(tv_vars_forimp)) %>%
  # only impute for observations when participants are alive/eligible in the sample
  # i.e. did not die, formally drop off
  filter(iwstat %in% c(1, 4))

# Drop wave (use age_v2 at each wave in the multi level multiple imputation)
impute_data <- hrs_preimp %>%
  select(-wave)

# # Sanity check
# length(unique(impute_data$HHIDPN))
# with(impute_data, table(wave))

# Look at missingness for participants imputed for each wave
hrs_preimp %>% 
  select(wave, all_of(tv_vars_forimp)) %>%
  group_by(wave) %>%
  dplyr::summarise_all(function(x) sum(is.na(x))/length(x)*100)

DataExplorer::plot_missing(impute_data %>% select(-HHIDPN))
missing_summary <- DataExplorer::profile_missing(impute_data %>% 
                                                   select(-HHIDPN)) %>%
  mutate(pct_missing = pct_missing*100) %>%
  arrange(pct_missing) %>%
  print(n = Inf)

writexl::write_xlsx(missing_summary,
                    here::here("output", "tables", "missing_summary_2000_2022.xlsx"))

ordered_vars <- paste0(missing_summary$feature)
baseline_vars_w_NA <- missing_summary %>%
  filter(feature %in% baseline_vars, pct_missing != 0) %>%
  pull(feature) %>% as.character()
tv_vars_w_NA <- missing_summary %>%
  filter(feature %in% tv_vars_forimp, pct_missing != 0) %>%
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
ncol(impute_data) # 1 more var: HHIDPN
sum(duplicated(c(cont_vars, binary_vars, ordinal_vars, cat_vars)))

# Check types
impute_data %>%
  map_dfc(function(x) class(x)[1]) %>%
  t() %>% print()

# Set classes by type
impute_data %<>%
  # mutate_all(as.numeric) %>%
  mutate(across(all_of(binary_vars), \(x) as.factor(x)),
         # # comment next line out to treat ordinal vars as continuous (it seems multilevel MI does not like ordered vars)
         # across(all_of(ordinal_vars), \(x) factor(x, ordered = T)),
         across(all_of(cat_vars), \(x) as.factor(x)),
         across(all_of(cont_vars), \(x) as.numeric(x))) %>%
  # arrange by order
  select(HHIDPN, all_of(ordered_vars))

# # Sanity check
impute_data %>%
  map_dfc(function(x) class(x)[1]) %>%
  t() %>% print()

##---- Rescale variables to reduce warning ----
# total household income is very large comparing to other variable values
# Previous mice imputation had warnings about this
# Warning messages:
# Some predictor variables are on very different scales: consider rescaling

# Rescale income
impute_data %<>%
  mutate(itot = itot/1000000) # per million

# Rescale weight
impute_data %<>%
  mutate(r5wtresp = r5wtresp/1000) # per thousand people

impute_data %>%
  select(all_of(cont_vars)) %>%
  apply(., 2, summary)

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

# # use cart for all 2lonly methods (attempt for using where argument)
# imputationFunction_cart <- rep("cart", length(baseline_vars_w_NA)) %>% as.list()
# names(imputationFunction_cart) <- c(baseline_vars_w_NA)
# imputationFunction_cart

# define cluster variable, i.e. nested in HHIDPN
cluster_var <- rep("HHIDPN", length(baseline_vars_w_NA)) %>% as.list()
names(cluster_var) <- c(baseline_vars_w_NA)
cluster_var

##---- Run imputation ----
# (start <- Sys.time())
# imp_pmm_1 <- mice::mice(impute_data, m = 1, maxit = 1,
#                                predictorMatrix = predM1,
#                                method = impM1,
#                                imputationFunction = imputationFunction,
#                                cluster_var = cluster_var,
#                                seed = 62283)
# (end <- Sys.time())
# end - start
# # m = 1, maxit = 1: 2.15 min in total

# futuremice parallel on imputation not iterations.
(start <- Sys.time())
imp_pmm_10 <- mice::futuremice(impute_data, m = 10, maxit = 10,
                               predictorMatrix = predM1,
                               method = impM1,
                               imputationFunction = imputationFunction,
                               cluster_var = cluster_var,
                               parallelseed = 62283,
                               packages = "miceadds",
                               n.core = 10)
(end <- Sys.time())
end - start
# # m = 10, maxit = 10: 31.12985 mins parallel using 10 cores

save(imp_pmm_10, file = here::here("data", "analysis_data", "imp_pmm_10_outcome_2000_2022_ot52l90.RData"))

imp_pmm_10$loggedEvents # dropping iwstat for a bunch of tv vars, can be ignored
complete(imp_pmm_10) %>%
  apply(., 2, function(x) sum(is.na(x))/length(x)*100)
plot(imp_pmm_10)
# # Look at density plots by individual variables
# cont_vars
# densityplot(imp_pmm_10, ~raedyrs) # very similar distributions
# densityplot(imp_pmm_10, ~itot)
# densityplot(imp_pmm_10, ~hhres)
# densityplot(imp_pmm_10, ~r5wtresp)
# densityplot(imp_pmm_10, ~memimp)

#---- Stack imputed dataset ----
load(here::here("data", "analysis_data", "imp_pmm_10_outcome_2000_2022_ot52l90.RData"))

imp_num = 10
imp_temp <- list()
for (i in 1:10){
  temp_bl <- complete(get(paste0("imp_pmm_", imp_num)), action = i) %>%
    select(HHIDPN, all_of(baseline_vars)) %>% unique() %>%
    left_join(hrs_forimp %>% select(HHIDPN, r5age_v2), by = "HHIDPN") %>%
    left_join(hrs_selected %>% select(HHIDPN, rabyear, household_min_pn), by = "HHIDPN")
  temp_tv_long <- complete(get(paste0("imp_pmm_", imp_num)), action = i) %>%
    select(-all_of(baseline_vars)) %>%
    mutate(iwstat = as.numeric(as.character(iwstat))) %>%
    # get iwstat 1:4 (we impute memimp now)
    left_join(hrs_forimp_long %>%
                select(HHIDPN, wave, iwstat, age_v2, dead, inw),
              by = c("HHIDPN", "age_v2", "iwstat")) %>%
    # merge back iwstat 5:7
    rbind(., hrs_forimp_long %>% 
            select(HHIDPN, wave, all_of(tv_vars_forimp), memimp, dead, inw) %>%
            filter(!iwstat %in% c(1, 4)))
  imp_temp[[i]] <- temp_tv_long %>%
    left_join(temp_bl, by = "HHIDPN")
  
  imp_temp[[i]][, "imp"] <- i
}

imp_pmm_stacked <- do.call(rbind, imp_temp)

# Sanity check
colnames(imp_pmm_stacked)
with(imp_pmm_stacked, table(imp, useNA = "ifany"))
with(imp_pmm_stacked, table(iwstat, useNA = "ifany"))
with(imp_pmm_stacked, table(iwstat, inw, useNA = "ifany"))
with(imp_pmm_stacked, table(iwstat, dead, useNA = "ifany"))
with(imp_pmm_stacked %>% filter(imp == 1), length(unique(HHIDPN)))
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
  filter(iwstat %in% c(1, 4)) %>%
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
  mutate(
    r5wtresp = r5wtresp*1000, # scale it back from per thousand
    income_pp = itot/sqrt(hhres)*1000000, # scale it back from per million
    mstat_2cat = case_when(mstat %in% c(1, 2, 3) ~ 1, # Married/partnered
                           # Separated/divorced/widowed/not married/Widowed
                           mstat %in% c(4, 5, 6, 7, 8) ~ 0),
    adl_2cat = case_when(adl6a == 0 ~ 0,
                         adl6a > 0 ~ 1),
    iadl_2cat = case_when(iadl5a == 0 ~ 0,
                          iadl5a > 0 ~ 1),
    inw_imp = case_when(iwstat %in% c(1, 4) ~ 1,
                        iwstat %in% c(0, 5, 6, 7, 9) ~ 0)) %>%
  mutate(across(all_of(c(binary_vars, cat_vars)), \(x) as.numeric(as.character(x))))

# Restandardize memory score
# using weighted mean and sd from wave 5
r5memimp_mean_sd_imp <- imp_pmm_stacked %>% filter(wave == 5) %>%
  group_by(imp) %>%
  summarise(r5memimp_mean = Hmisc::wtd.mean(memimp, weights = r5wtresp),
            r5memimp_sd = sqrt(Hmisc::wtd.var(memimp, weights = r5wtresp)))

imp_pmm_stacked %<>%
  left_join(r5memimp_mean_sd_imp, by = "imp") %>%
  mutate(memimp_stdz = (memimp - r5memimp_mean)/r5memimp_sd)

# # Sanity check
# with(imp_pmm_stacked %>% filter(wave == 5, imp == 1),
#      Hmisc::wtd.mean(memimp_stdz, weights = r5wtresp))
# with(imp_pmm_stacked %>% filter(wave == 5, imp == 1),
#      Hmisc::wtd.var(memimp_stdz, weights = r5wtresp))

# Remove values for flusht for waves 6, ..., 14 (wave 4 represents imputed wave 3 flu shot)
imp_pmm_stacked %<>%
  mutate(flusht = case_when(wave %% 2 == 0 & wave != 4 ~ NA_real_,
                            TRUE ~ flusht))

# Sanity check
with(imp_pmm_stacked, table(flusht, wave, useNA = "ifany"))
imp_pmm_stacked %>% filter(wave == 5) %>% 
  with(., Hmisc::wtd.mean(memimp_stdz, weights = r5wtresp))

#---- Check Memory score imputed analytic sample by cohort ----
hrs_imp_r5resp_memimp_long <- imp_pmm_stacked %>%
  left_join(hrs_selected %>% select(HHIDPN, racohbyr), by = "HHIDPN") %>%
  filter(imp == 1) %>%
  mutate(wave = as.numeric(wave))

# Imputed
hrs_imp_r5resp_memimp_mean_cohort <- hrs_imp_r5resp_memimp_long %>%
  group_by(wave, racohbyr) %>%
  summarise(mean_memimp = mean(memimp, na.rm = TRUE),
            se_memimp = sd(memimp, na.rm = TRUE)/sqrt(sum(!is.na(memimp))),
            mean_memimp_stdz = mean(memimp_stdz, na.rm = TRUE),
            se_memimp_stdz = sd(memimp_stdz, na.rm = TRUE)/sqrt(sum(!is.na(memimp_stdz))))

# Not imputed
hrs_memimp_long <- hrs_selected %>%
  select(HHIDPN, contains("memimp")) %>%
  pivot_longer(cols = contains("memimp"),
               names_to = "wave",
               values_to = "memimp",
               names_pattern = "r(.*)memimp") %>%
  mutate(wave = as.numeric(wave)) %>%
  filter(wave >= 4)

r5memimp_mean <- Hmisc::wtd.mean(hrs_selected$r5memimp, weights = hrs_selected$r5wtresp)
r5memimp_sd <- sqrt(Hmisc::wtd.var(hrs_selected$r5memimp, weights = hrs_selected$r5wtresp))
hrs_memimp_mean_cohort <- hrs_memimp_long %>%
  mutate(memimp_stdz = (memimp - r5memimp_mean)/r5memimp_sd) %>%
  left_join(hrs_selected %>% select(HHIDPN, racohbyr), by = "HHIDPN") %>%
  group_by(wave, racohbyr) %>%
  summarise(mean_memimp = mean(memimp, na.rm = TRUE),
            se_memimp = sd(memimp, na.rm = TRUE)/sqrt(sum(!is.na(memimp))),
            mean_memimp_stdz = mean(memimp_stdz, na.rm = TRUE),
            se_memimp_stdz = sd(memimp_stdz, na.rm = TRUE)/sqrt(sum(!is.na(memimp_stdz))))

hrs_memimp_preandpost_imp_cohort <- hrs_imp_r5resp_memimp_mean_cohort %>%
  mutate(sample = "imputed") %>%
  bind_rows(hrs_memimp_mean_cohort %>%
              mutate(sample = "pre imputation")) %>%
  mutate(birth_cohort = case_when(racohbyr == 1 ~ "AHEAD",
                                  racohbyr == 2 ~ "CODA",
                                  racohbyr == 3 ~ "HRS",
                                  racohbyr == 4 ~ "War babies",
                                  racohbyr == 5 ~ "Early Baby Boomers",
                                  racohbyr == 6 ~ "Mid Baby Boomers",
                                  racohbyr == 7 ~ "Late Baby Boomers"))

hrs_memimp_preandpost_imp_cohort %>%
  filter(racohbyr != 0) %>%
  mutate(year = (wave - 5)*2 + 2000) %>%
  ggplot(aes(x = year, y = mean_memimp, color = birth_cohort)) +
  geom_point() +
  geom_line(aes(linetype = as.factor(sample))) +
  geom_errorbar(aes(ymin = mean_memimp - 1.96*se_memimp,
                    ymax = mean_memimp + 1.96*se_memimp),
                width = 0.2) +
  # scale_y_continuous(limits = c(0.2, 1.3)) +
  scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
  # scale_linetype_manual(values = c("solid", "dashed"), name = "cohort", labels = c("AHEAD", "HRS")) +
  labs(title = "Mean memory score over time",
       x = "year",
       y = "Mean memory score") +
  theme_bw()

#---- Save imputed data ----
saveRDS(imp_pmm_stacked, here::here("data", "analysis_data", "imp_pmm_stacked_outcome_2000_2022_ot52l90.RDS"))

#---- OLD ----
# ##---- Check the distribution of imputed memimp ----
# hrs_imp_r5resp_memimp_long <- imp_pmm_stacked %>%
#   filter(imp == 1) %>%
#   mutate(wave = as.numeric(wave))
# hrs_imp_r5resp_memimp_mean <- hrs_imp_r5resp_memimp_long %>%
#   group_by(wave) %>%
#   summarise(mean_memimp = mean(memimp, na.rm = TRUE),
#             se_memimp = sd(memimp, na.rm = TRUE)/sqrt(sum(!is.na(memimp)))) %>%
# filter(wave != 2)
# 
# # Only imputed values
# hrs_imponly_r5resp_memimp_long <- imp_pmm_stacked %>%
#   filter(imp == 1, iwstat %in% c(1, 4)) %>%
#   select(HHIDPN, wave, memimp, memimp_stdz) %>%
#   left_join(hrs_preimp %>% select(HHIDPN, wave, memimp) %>% rename(memimp_pre = memimp), 
#             by = c("HHIDPN", "wave")) %>%
#   mutate(impute_mem = case_when(is.na(memimp_pre) & !is.na(memimp) ~ 1,
#                                 !is.na(memimp_pre) ~ 0))
# # Sanity check
# with(hrs_imponly_r5resp_memimp_long %>% filter(wave != 2), table(impute_mem, wave, useNA = "ifany"))
# with(hrs_imponly_r5resp_memimp_long %>% filter(impute_mem == 1), summary(memimp))
# with(hrs_imponly_r5resp_memimp_long %>% filter(impute_mem == 0), summary(memimp))
# with(hrs_imponly_r5resp_memimp_long %>% filter(impute_mem == 1, wave == 3), summary(memimp))
# with(hrs_imponly_r5resp_memimp_long %>% filter(impute_mem == 0, wave == 3), summary(memimp))
# with(hrs_imponly_r5resp_memimp_long %>% filter(wave == 3), summary(memimp))
# 
# hrs_imp_r5resp_memimp_mean %>%
#   mutate(year = (wave - 5)*2 + 2000) %>%
#   ggplot(aes(x = year, y = mean_memimp)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = mean_memimp - 1.96*se_memimp,
#                     ymax = mean_memimp + 1.96*se_memimp),
#                 width = 0.2) +
#   # scale_y_continuous(limits = c(0.45, 1.25)) +
#   scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
#   labs(title = "Mean memory score over time (analytic sample with imputation: imp = 1)",
#        x = "year",
#        y = "Mean memory score") +
#   theme_bw()
# 
# ##---- Analytic sample not imputed ----
# hrs_memimp_mean <- hrs_memimp_long %>%
#   group_by(wave) %>%
#   summarise(mean_memimp = mean(memimp, na.rm = TRUE),
#             se_memimp = sd(memimp, na.rm = TRUE)/sqrt(sum(!is.na(memimp))))
# 
# hrs_memimp_mean %>%
#   mutate(year = (wave - 5)*2 + 2000) %>%
#   ggplot(aes(x = year, y = mean_memimp)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = mean_memimp - 1.96*se_memimp,
#                     ymax = mean_memimp + 1.96*se_memimp),
#                 width = 0.2) +
#   scale_y_continuous(limits = c(0.45, 1.25)) +
#   scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
#   labs(title = "Mean memory score over time (analytic sample pre imputation)",
#        x = "year",
#        y = "Mean memory score") +
#   theme_bw()
# 
# hrs_memimp_preandpost_imp <- hrs_imp_r5resp_memimp_mean %>%
#   mutate(sample = "imputed") %>%
#   bind_rows(hrs_memimp_mean %>%
#               mutate(sample = "pre imputation"))
# 
# hrs_memimp_preandpost_imp %>%
#   mutate(year = (wave - 5)*2 + 2000) %>%
#   ggplot(aes(x = year, y = mean_memimp, color = sample)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = mean_memimp - 1.96*se_memimp,
#                     ymax = mean_memimp + 1.96*se_memimp),
#                 width = 0.2) +
#   scale_y_continuous(limits = c(0.45, 1.25)) +
#   scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
#   labs(title = "Mean memory score over time",
#        x = "year",
#        y = "Mean memory score") +
#   theme_bw()