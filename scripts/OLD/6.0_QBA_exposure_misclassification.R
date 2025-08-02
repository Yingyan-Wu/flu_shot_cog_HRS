#  QBA for flu shot coverage
# Yingyan Wu
# Not finished
# On hold for now

#---- Package loading + option specifications ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "here", "episensr", "survival", "gtsummary")

#---- Load data ----
hrs_mci_imp1 <- 
  readRDS(here::here("data", "model_results", "lmtp_results_mci_risks.RDS")) %>%
  filter(imp == 1, type != "0sus")# only 1sus and natural course

load(here::here("data", "analysis_data", "data_for_cluster_20250502.RData"))
load(here::here("data", "model_results", "wtd_mci_threshold_imp_hrs_ahead.RData"))

#---- Prepare the data ----
hrs_imp1_r3resp_wide <- hrs_imp_r3resp_wide %>%
  # Merge MCI to the dataset
  # Observed MCI
  left_join(wtd_mci_threshold %>% filter(imp == 1), by = "imp") %>%
  mutate(across(paste0("r", 3:15, "memimp_stdz"), 
                function(x) case_when(x <= mci_threshold ~ 1,
                                      x > mci_threshold ~ 0),
                .names = "{col}_mci_obs")) %>%
  rename_at(vars(contains("memimp_stdz_mci_obs")), 
            ~str_remove(., "memimp_stdz_")) %>%
  # Filter to imp 1 and select related variables
  filter(imp == 1) %>%
  # TMLE MCI: 1sus
  left_join(hrs_mci_imp1 %>% filter(type == "1sus") %>%
              select(HHIDPN, paste0("po_w", 3:15, "_mci_pe")) %>%
              set_colnames(c("HHIDPN", paste0("r", 3:15, "mci_1sus"))),
            by = "HHIDPN") %>%
  # TMLE MCI: NC
  left_join(hrs_mci_imp1 %>% filter(type == "nc") %>%
              select(HHIDPN, paste0("po_w", 3:15, "_mci_pe")) %>%
              set_colnames(c("HHIDPN", paste0("r", 3:15, "mci_nc"))),
            by = "HHIDPN")


# Sustained flu shot variables
# Create sustained flushot variables without a loop
# Build list of expressions for sustained flu shot indicators
waves <- 3:15
hrs_imp1_r3resp_wide <- hrs_imp1_r3resp_wide %>%
  mutate(across(paste0("r", waves, flusht),
                function(x) case_when(
                  rowSums(across(paste0("r", waves, flusht))) == 7 ~ 1,
                  TRUE ~ 0),
    !!!setNames(
      lapply(waves, function(w) {
        needed_vars <- paste0("r", waves[waves <= w], "flusht")
        expr(if_else(rowSums(across(all_of(!!needed_vars))) == 
                       length(needed_vars), 1, 0))
      }),
      paste0("r", waves, "flusht_1sus")
    )
  )

# Categorize flu shot: sustained 1 till each wave
mutate(
  
  r3flusht_1sus = case_when(
    rowSums(across(paste0("r", 3, "flusht"))) == 7 ~ 1,
    TRUE ~ 0),
  r4flusht_1sus = r3flusht_1sus,
  r5flusht_1sus = case_when(
    rowSums(across(paste0("r", seq(3, 5, by = 2), "flusht"))) == 7 ~ 1,
    TRUE ~ 0),
  r6flusht_1sus = r5flusht_1sus,
  r7flusht_1sus = case_when(
    rowSums(across(paste0("r", seq(3, 7, by = 2), "flusht"))) == 7 ~ 1,
    TRUE ~ 0),
  r8flusht_1sus = r7flusht_1sus,
  
  
  r3flusht_1sus = case_when(
    rowSums(across(paste0("r", seq(3, 15, by = 2), "flusht"))) == 7 ~ 1,
    TRUE ~ 0))

# Sanity check
with(hrs_imp1_r3resp_wide, table(flusht_1sus, useNA = "ifany"))
with(hrs_imp1_r3resp_wide, table(r15mci_1sus, useNA = "ifany"))
with(hrs_imp1_r3resp_wide, table(r15mci_nc, useNA = "ifany"))
with(hrs_imp1_r3resp_wide, table(r15mci_obs, useNA = "ifany"))

# QBA dataset using observed flu shot and observed MCI at survivors at each wave 
# (long for years)
hrs_imp1_qba <- hrs_imp1_r3resp_wide %>%
  select(HHIDPN, paste0("r", 3:15, "flusht"), paste0("r", 3:15, "mci_obs"),
         paste0("r", 3:15, "inw_imp")) %>%
  pivot_longer(cols = -c(HHIDPN),
               names_to = c("wave", ".value"),
               names_pattern = "r(\\d+)(.*)") %>%
  mutate(wave = as.numeric(wave)) %>%
  filter(inw_imp != 0)

#---- FLu shot Proportions ----
# McLean et al., 2014 community cohort in Wiscosin
# Table 2 Current-Season Vaccination and Frequent Vaccinee (>= 4 out of 5) 50 and above
(317+285+318+190) # 1110
# 50 and above study population from Table 1
(175+524+61+657+96+374+46+429+84+269+21+325+46+161+4+186) # 3458

1110/3458 # 0.32

# In a older cohorts in Sweden, it's 38%. (Örtqvist et al., 2018)
# The upper end can not be larger than 0.32

with(hrs_imp1_r3resp_wide %>% filter(rowSums(across(ends_with("inw_imp"))) == 13), 
     table(flusht_1sus, useNA = "ifany"))
# 1sus proportion
# 1159/(1159+3705) 0.24
# Assuming 0.32 is the actual proprotion, sens of exp among controls was 
# calculated by:   # Specificity: P(Test = no disease | Truth = no disease)
# Sensitivity: P(Test = exp | Truth = exp) = P(Test = exp, Truth = exp)/P(Truth = exp)
# 0.24/0.32
# Current sens for flu shot 1sus around = 0.75

#---- QBA ----
flusht_mis <- function(data, exposure, outcome, w) {
  # # Test
  # data = hrs_imp1_qba
  # w = 8
  # exposure = "flusht"
  # outcome = "mci_obs"
  
  crosstab <- data %>% 
    filter(wave == w) %>% 
    # count cases and exposures
    with(table(get(outcome) == 0, get(exposure) == 0, useNA = "ifany"))
  colnames(crosstab) <- c("exposed", "unexposed")
  rownames(crosstab) <- c("mci", "no mci")
  
  sens_exp <- seq(0.75, 1, by = 0.05)
  
  # The specificity should be high for the flu shot status. Specificity should
  # be high for both the flu shot sus arm and natural course arm. 
  
  for (sens_exp_1 in sens_exp) {
    
    res <- misclassification(
      crosstab,
      type = "exposure", # flu shot
      # sens of exp among cases, sens of exp among controls
      # spec of exp among cases, spec of exp among controls
      bias_parms = c(sens_exp_1, 0.9, 0.99, 0.99),
      alpha = 0.05
    )
    
    if (sens_exp_1 == sens_exp[1]) {
      out <- rbind(res$obs.measures[, 1], res$adj.measures[, 1]) %>% 
        as_tibble() %>% 
        mutate(type = c("obs", "adj"), 
               sens_exp = sens_exp_1,
               wave = w) %>% 
        setNames(c("RR", "OR", "Type", "Sens_exp_control", "wave"))
    } else {
      out <- t(res$adj.measures[, 1]) %>% as_tibble() %>% 
        mutate(type = c("adj"), sens_exp = sens_exp_1,
               wave = w) %>% 
        setNames(c("RR", "OR", "Type", "Sens_exp_control", "wave")) %>% 
        add_row(out, .)
    }
  }
  
  return(out)
}

# Starting from wave 8 (after 10 years)
qba_flusht_8_15 <- lapply(8:15, function(x) flusht_mis(hrs_imp1_qba, 
                                                       exposure = "flusht",
                                                       outcome = "mci_obs",
                                                       w = x)) %>% 
  bind_rows() 

qba_flusht_8_15 %>% 
  ggplot(
    aes(x = wave, y = RR, group = interaction(Type, Sens_exp_control), 
        color = Sens_exp_control, linetype = Type)) + 
  geom_line() + 
  theme_bw()



#---- Flu shot coverage observed proportion ----
# in HRS non censoring till the end of follow up
noncens_n <- hrs_imp1_r3resp_wide %>%
  filter(rowSums(across(ends_with("inw_imp"))) == 13) %>% nrow()
# # Sanity check
# with(hrs_imp1_r3resp_wide %>%
#        filter(r4inw_imp == 1, r5inw_imp == 1, r6inw_imp == 1,
#               r7inw_imp == 1, r8inw_imp == 1, r9inw_imp == 1,
#               r10inw_imp == 1, r11inw_imp == 1, r12inw_imp == 1,
#               r13inw_imp == 1, r14inw_imp == 1), table(r15inw_imp, useNA = "ifany"))

hrs_imp1_r3resp_wide %>%
  filter(rowSums(across(ends_with("inw_imp"))) == 13) %>%
  select(paste0("r", seq(3, 15, 2), "flusht")) %>%
  apply(., 2, table, useNA = "ifany")/noncens_n
# r3flusht  r5flusht  r7flusht  r9flusht r11flusht r13flusht r15flusht
# 0 0.6535773 0.4841694 0.3865132 0.3289474 0.2713816 0.2436266  0.208676
# 1 0.3464227 0.5158306 0.6134868 0.6710526 0.7286184 0.7563734  0.791324

# create non censoring version flushot for each wave
for (w in seq(3, 15, 2)){
  hrs_imp1_r3resp_wide %<>%
    mutate(!!paste0("r", w, "flusht_noncens") := case_when(
      get(paste0("r", w, "inw_imp")) == 1 ~ get(paste0("r", w, "flusht")),
      get(paste0("r", w, "inw_imp")) == 0 ~ NA_real_))
}

hrs_imp1_r3resp_wide %>%
  select(paste0("r", seq(3, 15, 2), "flusht_noncens")) %>%
  apply(., 2, table, useNA = "ifany")
hrs_imp1_r3resp_wide %>%
  select(paste0("r", seq(3, 15, 2), "flusht")) %>%
  apply(., 2, table, useNA = "ifany")

hrs_imp1_r3resp_wide %>%
  select(paste0("r", seq(3, 15, 2), "flusht_noncens")) %>%
  tbl_summary(statistic = list(all_categorical() ~ "{n} ({p})"),
              digits = list(all_categorical() ~ c(0, 1)),
              missing = "no")

#---- OLD ----
# TMLE MCI QBA
# QBA dataset (long for years)
hrs_imp1_qba <- hrs_imp1_r3resp_wide %>%
  # Exposed: sustained flu shot from wave 3 to wave 15
  filter(flusht_1sus == 1) %>%
  select(-contains("mci_nc")) %>%
  rename_with(~str_replace(., "mci_1sus", "mci"),
              cols = vars(paste0("r", 3:15, "mci_1sus"))) %>%
  # Unexposed: natural course
  rbind(hrs_imp1_r3resp_wide %>%
          mutate(flusht_1sus = 0) %>%
          select(-contains("mci_1sus")) %>%
          rename_with(~str_replace(., "mci_nc", "mci"),
                      cols = vars(paste0("r", 3:15, "mci_nc")))) %>%
  pivot_longer(cols = -c(HHIDPN, flusht_1sus),
               names_to = c("wave", ".value"),
               names_pattern = "r(\\d+)(.*)") %>%
  mutate(wave = as.numeric(wave)) %>%
  filter(inw_imp != 0)

# crosstab_obs <- # count cases and exposures
#   with(hrs_imp1_qba, table(r15mci == 0, flusht_1sus == 0, useNA = "ifany"))
# colnames(crosstab_obs) <- c("1sus", "nc")
# rownames(crosstab_obs) <- c("mci", "no mci")

#---- FLu shot Proportions ----
# McLean et al., 2014 community cohort in Wiscosin
# Table 2 Current-Season Vaccination and Frequent Vaccinee (>= 4 out of 5) 50 and above
(317+285+318+190) # 1110
# 50 and above study population from Table 1
(175+524+61+657+96+374+46+429+84+269+21+325+46+161+4+186) # 3458

1110/3458 # 0.32

# In a older cohorts in Sweden, it's 38%. (Örtqvist et al., 2018)
# The upper end can not be larger than 0.32

with(hrs_imp1_r3resp_wide %>% filter(rowSums(across(ends_with("inw_imp"))) == 13), 
     table(flusht_1sus, useNA = "ifany"))
# 1sus proportion
# 1159/(1159+3705) 0.24
# Assuming 0.32 is the actual proprotion, sens of exp among controls was 
# calculated by:   # Specificity: P(Test = no disease | Truth = no disease)
# Sensitivity: P(Test = exp | Truth = exp) = P(Test = exp, Truth = exp)/P(Truth = exp)
# 0.24/0.32
# Current sens for flu shot 1sus around = 0.75

#---- QBA ----
flusht_mis <- function(data, w) {
  # data <- hrs_imp1_qba
  # w = 8
  
  crosstab <- data %>% 
    filter(wave == w) %>% 
    # count cases and exposures
    with(table(mci == 0, flusht_1sus == 0, useNA = "ifany"))
  colnames(crosstab) <- c("1sus", "nc")
  rownames(crosstab) <- c("mci", "no mci")
  
  # crosstab <- crosstab_obs
  # crosstab_perf
  sens_exp <- seq(0.70, 1, by = 0.05)
  
  # The specificity should be high for the flu shot status. Specificity should
  # be high for both the flu shot sus arm and natural course arm. 
  
  for (sens_exp_1 in sens_exp) {
    
    res <- misclassification(
      crosstab,
      type = "exposure", # flu shot
      # sens of exp among cases, sens of exp among controls
      # spec of exp among cases, spec of exp among controls
      bias_parms = c(sens_exp_1, 0.9, 0.99, 0.99),
      alpha = 0.05
    )
    
    if (sens_exp_1 == sens_exp[1]) {
      out <- rbind(res$obs.measures[, 1], res$adj.measures[, 1]) %>% 
        as_tibble() %>% 
        mutate(type = c("obs", "adj"), 
               sens_exp = sens_exp_1,
               wave = w) %>% 
        setNames(c("RR", "OR", "Type", "Sens_exp_control", "wave"))
    } else {
      out <- t(res$adj.measures[, 1]) %>% as_tibble() %>% 
        mutate(type = c("adj"), sens_exp = sens_exp_1,
               wave = w) %>% 
        setNames(c("RR", "OR", "Type", "Sens_exp_control", "wave")) %>% 
        add_row(out, .)
    }
  }
  
  return(out)
}

# Starting from wave 8 (after 10 years)

qba_flusht_8_15 <- lapply(8:15, function(x) flusht_mis(hrs_imp1_qba, w = x)) %>% 
  bind_rows() 

qba_flusht_8_15 %>% 
  ggplot(
    aes(x = wave, y = RR, group = interaction(Type, Sens_exp_control), 
        color = Sens_exp_control, linetype = Type)) + 
  geom_line() + 
  theme_bw()

#---- Attempt: manual bias analysis ----
# # To be able to work around episensr, create a crosstab that has the true sens.
# crosstab_obs <- # count cases and exposures
#   with(hrs_imp1_qba, table(r15mci == 0, flusht_1sus == 0, useNA = "ifany"))
# colnames(crosstab_obs) <- c("1sus", "nc")
# rownames(crosstab_obs) <- c("mci", "no mci")
# 
# #          1sus   nc
# # mci     702 3743
# # no mci  457 1121
# 
# # mci n=4445
# # no mci n=1578
# # 1sus: 1159
# # nc: 4864
# # Observed RR: 0.79
# 
# true_p <- function(sens, spec, observed_p){
#   # Sensitivity: P(Test+|Truth+)
#   # Specificity: P(Test-|Truth-)
#   # p: P(Truth+)
#   # observed_p: P(Test+)
#   # observed_p = p*sens + (1 - p)*(1 - spec)
#   # Thus P can be back calculated as follows:
#   p <- (observed_p + spec - 1)/(sens + spec - 1)
#   return(p)
# }
# 
# # Assume no differential misclassification
# # If sensitivity is 0.75, specificity is 0.99 for the exposure, recover the perfect crosstab
# n_mci_1sus <- round(true_p(0.75, 0.99, 702/1159) * 1159) # n=933
# n_no_mci_1sus <- round(true_p(0.75, 0.99, 457/1159) * 1159) # n=602
# # Recovered 2*2
# #          1sus   nc
# # mci     933 4445-933=3512
# # no mci  602 1578-602=976
# # RR: (933/(933+602))/(3512/(3512+976)) = 0.777
# 
# crosstab_perf <- matrix(c(933, 3512, 602, 976),
#                         dimnames = list(c("mci", "no mci"), # rownames: disease
#                                         c("1sus", "nc")), # colnames: exposure
#                         nrow = 2, byrow = TRUE)




