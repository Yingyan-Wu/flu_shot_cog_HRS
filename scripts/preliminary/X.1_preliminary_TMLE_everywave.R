# TMLE (super learner) IPCW every wave covariates
# Created by: Yingyan Wu
# Date: Aug.19.2024

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
hrs_imp <- readRDS(here::here("data", "analysis_data", "imp_pmm_stacked_ahead_hrs_outcome.RDS"))

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
             # Follow up status #"inw", "dead",
             "iwstat", "inw_imp",
             # # Survey weights
             # "wtresp",
             # Outcome (at previous waves)
             "memimp")

hrs_imp_r3resp_wide <- hrs_imp %>%
  select(HHIDPN, wave, imp, all_of(c(bl_vars, tv_vars)), -r3age_v2) %>%
  mutate(
    # for those who died this wave, set the flusht value to be 0 for this wave
    # These people were not used anyway as the outcome at this wave is still missing for these people
    # this is just for the lmtp_tmle function to be able to run
    flusht = case_when(iwstat %in% c(5, 6, 7, 9) ~ 0,
                       TRUE ~ flusht)) %>%
  pivot_wider(names_from = "wave", 
              values_from = c(all_of(tv_vars)),
              names_glue = "r{wave}{.value}") %>%
  select(-paste0("r", seq(2, 14, by = 2), "flusht"), -contains("r2memimp"))

# Descriptives:
colnames(hrs_imp_r3resp_wide)
for (w in 4:15){
  print(with(hrs_imp_r3resp_wide, table(!is.na(get(paste0("r", w, "memimp"))), 
                                        get(paste0("r", w, "inw_imp")), 
                                        useNA = "ifany")))
}
hrs_imp_r3resp_wide %>%
  select(contains("memimp")) %>%
  apply(., 2, summary)
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
for (w in flusht_waves[flusht_waves != 15]){
  hrs_imp_r3resp_wide %<>%
    mutate(!!paste0("r", w + 1, "flusht") := !!sym(paste0("r", w, "flusht")))
}
# Sanity check
hrs_imp_r3resp_wide %>%
  select(paste0("r", 3:15, "flusht")) %>%
  apply(., 2, table, useNA = "ifany")

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

##---- Censoring ----
# 1 if in the study at this wave, 0 otherwise
waves_fu <- waves[waves != 2]
C_allwaves <- lapply(waves_fu, function(w) {
  # Consider only waves before the current one
  uptonow_waves <- waves_fu[waves_fu <= w]
  paste0("r", uptonow_waves, "inw_imp")
})

# # Assign separate list objects to individual variables
# for (i in seq_along(lists)) {
#   assign(paste0("C_r", waves_fu[i]), lists[[i]])
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

#---- Models ----
hrs_imp_r3resp_wide_imp1 <- hrs_imp_r3resp_wide %>% filter(imp == 1)

##---- Natural course ----
# r3 outcome ~ r3 treatment, r2 tv vars, cens r3
# r5 outcome ~ r3, r5 trt, r2, r3 tv vars, cens r3, r5etc.
(start_nc <- Sys.time())
lmtp_nc <-  function(i) {
  lmtp_tmle(hrs_imp_r3resp_wide_imp1,
            trt = trt_allwaves[[i]],
            outcome = paste0("r", waves_fu[i], "memimp"),
            baseline = bl_vars,
            time_vary = L_allwaves[[i]],
            cens = C_allwaves[[i]],
            outcome_type = "continuous",
            mtp = FALSE, 
            weights = hrs_imp_r3resp_wide_imp1$r3wtresp,
            folds = 1)
}
# lmtp_nc(4)
# 
# colSums(is.na(hrs_imp_r3resp_wide_imp1 %>%
#                 select(all_of(c(trt_allwaves[[4]], paste0("r", waves_fu[4], "memimp"),
#                                 bl_vars, unlist(L_allwaves[[4]]), C_allwaves[[4]], "r3wtresp")))))

lmtp_nc_results <- lapply(seq_along(waves_fu), lmtp_nc)
(end_nc <- Sys.time())
end_nc - start_nc # 51.03662 mins

##---- Sustained intervention ----
(start_sus <- Sys.time())
lmtp_sus <-  function(i, treated = 1) {
  lmtp_tmle(hrs_imp_r3resp_wide_imp1,
            trt = trt_allwaves[[i]],
            outcome = paste0("r", waves_fu[i], "memimp"),
            baseline = bl_vars,
            time_vary = L_allwaves[[i]],
            shift = case_when(treated == 1 ~ list(static_binary_on),# Turn all treatment nodes on
                              treated == 0 ~ list(static_binary_off))[[1]], # Turn all treatment nodes off
            cens = C_allwaves[[i]],
            outcome_type = "continuous",
            mtp = TRUE, 
            weights = hrs_imp_r3resp_wide_imp1$r3wtresp,
            folds = 1)
}
# lmtp_sus(4, treated = 1)
lmtp_sus1_results <- lapply(seq_along(waves_fu), function(i) lmtp_sus(i, treated = 1))
lmtp_sus0_results <- lapply(seq_along(waves_fu), function(i) lmtp_sus(i, treated = 0))
(end_sus <- Sys.time())
end_sus - start_sus # 1.993609 hours

##---- Baseline intervention ----
# Policy intervene function for lmtp
# Policy to intervene only on baseline flusht
policy_0 <- function(data, trt){
  if(trt == "r3flusht") {
    return(rep(0, each = length(data[[trt]])))
  } else {
    # this should not be used as the treatment is not the baseline treatment
    data[[trt]]
  }
}

policy_1 <- function(data, trt){
  if(trt == "r3flusht") {
    return(rep(1, each = length(data[[trt]])))
  } else {
    # this should not be used as the treatment is not the baseline treatment
    data[[trt]]
  }
}

(start_base <- Sys.time())
lmtp_base <-  function(i, treated = 1) {
  lmtp_tmle(hrs_imp_r3resp_wide_imp1,
            trt = trt_allwaves[[i]],
            outcome = paste0("r", waves_fu[i], "memimp"),
            baseline = bl_vars,
            time_vary = L_allwaves[[i]],
            shift = case_when(treated == 1 ~ list(policy_1),# Baseline treatment on
                              treated == 0 ~ list(policy_0))[[1]], # Baseline treatment off
            cens = C_allwaves[[i]],
            outcome_type = "continuous",
            mtp = TRUE,
            weights = hrs_imp_r3resp_wide_imp1$r3wtresp,
            folds = 1)
}

# with IPCW, it has to be the same length for C, for L and for trt
lmtp_base1_results <- lapply(seq_along(waves_fu), function(i) lmtp_base(i, treated = 1))
lmtp_base0_results <- lapply(seq_along(waves_fu), function(i) lmtp_base(i, treated = 0))
(end_base <- Sys.time())

#---- Clean results ----
save(lmtp_sus1_results, lmtp_sus0_results, lmtp_nc_results, 
     file = here::here("data", "model_results", 
                       "lmtp_results_imp1_ahead_hrs_memimpimpute.RData"))

load(here::here("data", "model_results", 
                "lmtp_results_imp1_ahead_hrs_memimpimpute.RData"))

pe_se_to_tib <- function(lmtp_list){
  sapply(lmtp_list, function(x) c(x$theta, x$standard_error)) %>%
    t() %>%
    as.data.frame() %>%
    set_colnames(c("point_est", "se")) %>%
    mutate(wave = waves_fu)
}
results_tib <- rbind(
  pe_se_to_tib(lmtp_nc_results) %>% mutate(Intervention = "Natural course"), # NC
  pe_se_to_tib(lmtp_sus1_results) %>% 
    mutate(Intervention = "Continuous flu shot"), # always 1
  pe_se_to_tib(lmtp_sus0_results) %>%
    mutate(Intervention = "Continuous no flu shot")) %>% # always 0
  mutate(survey_year = (wave - 5)*2 + 2000,
         fu_year = survey_year - 1996)
# wave 5 corresponding to 2000, 4 years of follow up

# Trajectory plot
results_tib %>%
  # filter(Intervention != "Continuous no flu shot") %>%
  mutate(Intervention = factor(Intervention, 
                               levels = c("Continuous no flu shot", "Natural course", "Continuous flu shot"))) %>%
  # ggplot(aes(x = survey_year, y = point_est, color = Intervention)) +
  ggplot(aes(x = fu_year, y = point_est, color = Intervention)) +
  geom_point(position = position_dodge((width = 1))) +
  geom_line(position = position_dodge((width = 1))) +
  geom_errorbar(aes(ymin = point_est - 1.96*se, ymax = point_est + 1.96*se),
                position = position_dodge((width = 1)),
                width = 0.8) +
  # scale_y_continuous(breaks = seq(0, 1.25, by = 0.2)) +
  scale_x_continuous(breaks = seq(0, 24, by = 2)) +
  # scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
  scale_color_manual(name = NULL, values = c("#F7C173", "#4B4B4B", "#3271AE"),
                     labels = c("Continuous no flu shot", "Natural course", "Hypothetical flu shot intervention")) +# "#E5A84B",
  # geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Memory score trajectories under the hypothetical flu shot intervention and natural course",
       x = "Follow up years",
       y = "Memory z-score") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(here::here("output", "figures", "tmle_results_ahead_hrs_memimp_impute.jpg"), width = 9, height = 7)

##---- Contrasts ----
waves_fu
lmtp_contrast(lmtp_nc_results[[7]], ref = lmtp_nc_results[[1]])
lmtp_contrast(lmtp_sus1_results[[7]], ref = lmtp_sus1_results[[1]])
# lmtp_contrast(lmtp_sus0_results[[7]], ref = lmtp_sus0_results[[1]])

for (i in 1:length(waves_fu)){
  print(lmtp_contrast(lmtp_sus1_results[[i]], ref = lmtp_nc_results[[i]]))
}

# #---- Check the trend between 2018 and 2020 ----
# ##---- Analytic sample not imputed ----
# hrs_selected <- readRDS(here::here("data", "analysis_data", "hrs_selected.rds"))
# hrs_analytic <- hrs_selected %>%
#   filter(HHIDPN %in% r3_resp_id,
#          !HHIDPN %in% r2_notinsamp_id)
# hrs_memimp_long <- hrs_analytic %>%
#   select(HHIDPN, contains("memimp")) %>%
#   pivot_longer(cols = contains("memimp"),
#                names_to = "wave",
#                values_to = "memimp",
#                names_pattern = "r(.*)memimp") %>%
#   mutate(wave = as.numeric(wave))
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
#   scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
#   labs(title = "Mean memory score over time (analytic sample pre imputation)",
#        x = "year",
#        y = "Mean memory score") +
#   theme_bw()
# 
# ##---- Analytic sample in this script ----
# hrs_imp_r3resp_memimp_long <- hrs_imp_r3resp_wide_imp1 %>%
#   select(HHIDPN, contains("memimp")) %>%
#   pivot_longer(cols = contains("memimp"),
#                names_to = "wave",
#                values_to = "memimp",
#                names_pattern = "r(.*)memimp") %>%
#   mutate(wave = as.numeric(wave)) 
# hrs_imp_r3resp_memimp_mean <- hrs_imp_r3resp_memimp_long %>%
#   group_by(wave) %>%
#   summarise(mean_memimp = mean(memimp, na.rm = TRUE),
#             se_memimp = sd(memimp, na.rm = TRUE)/sqrt(sum(!is.na(memimp))))
# 
# hrs_imp_r3resp_memimp_mean %>%
#   mutate(year = (wave - 5)*2 + 2000) %>%
#   ggplot(aes(x = year, y = mean_memimp)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = mean_memimp - 1.96*se_memimp, 
#                     ymax = mean_memimp + 1.96*se_memimp),
#                 width = 0.2) +
#   scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
#   labs(title = "Mean memory score over time (analytic sample with imputation)",
#        x = "year",
#        y = "Mean memory score") +
#   theme_bw()
# 
# ##---- All imputed data ----
# hrs_imp1 <- hrs_imp %>%
#   filter(imp == 1) %>%
#   filter(wave != 2)
# 
# hrs_imp1_memimp_mean <- hrs_imp1 %>%
#   group_by(wave) %>%
#   summarise(mean_memimp = mean(memimp, na.rm = TRUE),
#             se_memimp = sd(memimp, na.rm = TRUE)/sqrt(sum(!is.na(memimp))))
# 
# length(unique(hrs_imp1$HHIDPN))
# 
# hrs_imp1_memimp_mean %>%
#   mutate(year = (wave - 5)*2 + 2000) %>%
#   ggplot(aes(x = year, y = mean_memimp)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = mean_memimp - 1.96*se_memimp, 
#                     ymax = mean_memimp + 1.96*se_memimp),
#                 width = 0.2) +
#   scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
#   labs(title = "Mean memory score over time (all imputed data)",
#        x = "Wave",
#        y = "Mean memory score") +
#   theme_bw()
# 
# # wave 3 outcome ~ wave 3 treatment, time invarying var, time varying var wave 2
# # wave 5 outcome ~ wave 3 and wave 5 treatment, time invarying, time varying var wave 2 wave 3
# # 
# # model2020 <- lm(r15memimp ~ ,
# #                 data = hrs_imp_r3resp_wide_imp1)
# 
# ##---- wave 3 and wave 15 participants ----
# r15_resp_id <- hrs_imp %>% select(HHIDPN, wave, iwstat) %>%
#   filter(wave == 15, iwstat == 1) %>% pull(HHIDPN) %>% unique()
# hrs_imp_r3resp_wide_imp1 %>%
#   filter(HHIDPN %in% r15_resp_id) %>% nrow()
# hrs_imp_r3resp_wide_imp1 %>%
#   filter(HHIDPN %in% r15_resp_id) %>%
#   select(contains("iwstat")) %>%
#   apply(., 2, table, useNA = "ifany")
# hrs_imp_r3resp_memimp_long_r15 <- hrs_imp_r3resp_wide_imp1 %>%
#   filter(HHIDPN %in% r15_resp_id) %>%
#   select(HHIDPN, contains("memimp")) %>%
#   pivot_longer(cols = contains("memimp"),
#                names_to = "wave",
#                values_to = "memimp",
#                names_pattern = "r(.*)memimp") %>%
#   mutate(wave = as.numeric(wave)) 
# hrs_imp_r3resp_memimp_mean_r15 <- hrs_imp_r3resp_memimp_long_r15 %>%
#   group_by(wave) %>%
#   summarise(mean_memimp = mean(memimp, na.rm = TRUE),
#             se_memimp = sd(memimp, na.rm = TRUE)/sqrt(sum(!is.na(memimp))))
# 
# hrs_imp_r3resp_memimp_mean_r15 %>%
#   mutate(year = (wave - 5)*2 + 2000) %>%
#   ggplot(aes(x = year, y = mean_memimp)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = mean_memimp - 1.96*se_memimp, 
#                     ymax = mean_memimp + 1.96*se_memimp),
#                 width = 0.2) +
#   scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
#   labs(title = "Mean memory score over time imputed (alive and in the sample throughout)",
#        x = "Wave",
#        y = "Mean memory score") +
#   theme_bw()
# 
# ##---- wave 3 and wave 9 participants ----
# r9_resp_id <- hrs_imp %>% select(HHIDPN, wave, iwstat) %>%
#   filter(wave == 9, iwstat == 1) %>% pull(HHIDPN) %>% unique()
# hrs_imp_r3resp_wide_imp1 %>%
#   filter(HHIDPN %in% r9_resp_id) %>% nrow()
# hrs_imp_r3resp_memimp_long_r9 <- hrs_imp_r3resp_wide_imp1 %>%
#   filter(HHIDPN %in% r9_resp_id) %>%
#   select(HHIDPN, contains("memimp")) %>%
#   pivot_longer(cols = contains("memimp"),
#                names_to = "wave",
#                values_to = "memimp",
#                names_pattern = "r(.*)memimp") %>%
#   mutate(wave = as.numeric(wave)) 
# hrs_imp_r3resp_memimp_mean_r9 <- hrs_imp_r3resp_memimp_long_r9 %>%
#   group_by(wave) %>%
#   summarise(mean_memimp = mean(memimp, na.rm = TRUE),
#             se_memimp = sd(memimp, na.rm = TRUE)/sqrt(sum(!is.na(memimp))))
# 
# hrs_imp_r3resp_memimp_mean_r9 %>%
#   ggplot(aes(x = wave, y = mean_memimp)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = mean_memimp - 1.96*se_memimp, ymax = mean_memimp + 1.96*se_memimp)) +
#   labs(title = "Mean memory score over time (wave 2, 3, 5 participants with imputation)",
#        x = "Wave",
#        y = "Mean memory score") +
#   theme_bw()
# 
# #---- OLD ----

# colSums(is.na(hrs_imp_r3resp_wide_imp1 %>%
#                 select(all_of(c(trt_allwaves[[1]], paste0("r", waves_fu[1], "memimp"),
#                                 bl_vars, unlist(L_allwaves[[1]]), C_allwaves[[1]], "r3wtresp")))))
# colSums(is.na(hrs_imp_r3resp_wide_imp1 %>%
#                 select(all_of(c(trt_allwaves[[5]], paste0("r", waves_fu[5], "memimp"),
#                                 bl_vars, unlist(L_allwaves[[5]]), C_allwaves[[5]], "r3wtresp")))))
# 
# with(hrs_imp_r3resp_wide_imp1 %>% filter(is.na(r5flusht)), table(r5iwstat, useNA = "ifany"))
# 
# 
# test_missing <- hrs_imp_r3resp_wide_imp1 %>%
#   select(HHIDPN, contains(c("flusht", "memimp", "inw_imp", "iwstat")))
# colSums(is.na(test_missing))
# with(test_missing, table(r3inw_imp, r5inw_imp))
# with(test_missing, table(r5inw_imp, r7inw_imp))
# with(test_missing, table(r7inw_imp, r9inw_imp))
# with(test_missing, table(r9inw_imp, r11inw_imp))
# 
# # Exclude that person attended wave 11 but did not attend wave 9
# test <- hrs_imp_r3resp_wide_imp1 %>% filter(!(r9inw_imp == 0 & r11inw_imp == 1))
# lmtp_tmle(test,
#           trt = c("r3flusht", "r5flusht", "r7flusht", "r9flusht", "r11flusht"),
#           outcome = "r11memimp",
#           baseline = bl_vars,
#           time_vary = L_allwaves[[5]],
#           cens = C_allwaves[[5]],
#           outcome_type = "continuous",
#           mtp = FALSE,
#           weights = test$r3wtresp,
#           folds = 1)
# # Turns out lmtp cannot deal with participants missing at one wave but come back at another wave
# test_missing %>% filter((r9inw_imp == 0 & r11inw_imp == 1)) %>%
#   view()
# 
# 
# # 
# # lmtp_tmle(hrs_imp_r3resp_wide_imp1,
# #           trt = trt_allwaves[[4]],
# #           outcome = paste0("r", waves_fu[4], "memimp"),
# #           baseline = bl_vars,
# #           time_vary = L_allwaves[[4]],
# #           cens = C_allwaves[[4]],
# #           outcome_type = "continuous",
# #           mtp = FALSE, 
# #           weights = hrs_imp_r3resp_wide_imp1$r3wtresp,
# #           folds = 1)
# # 
# # 
# # colSums(is.na(hrs_imp_r3resp_wide_imp1 %>% 
# #                 select(all_of(c(trt_allwaves[[3]], paste0("r", waves_fu[3], "memimp"),
# #                                 bl_vars, unlist(L_allwaves[[3]]), C_allwaves[[3]], "r3wtresp")))))
# 
