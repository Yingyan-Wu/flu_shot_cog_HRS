# Estimation of RR and RD for memory impairment (Sens HRS only)
# Yingyan Wu
# Same calculation and process as in Script 3.3 (using HRS only Sens dataset and cluster results)

#---- Package loading + option specifications ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "haven", "labelled",
       "splines", "glue", "ggplot2", "Hmisc")

options(sci = 999)

#---- Load data ----
load(here::here("data", "analysis_data", "data_for_cluster_hrsonly_20250502.RData"))
load(here::here("data", "model_results", "lmtp_results_SL_po_sens_HRSonly.RData"))
waves_fu = 3:15

#---- Categorize impairment ----
##---- Prep ----
# Table e-3 Petersen RC, Lopez O, Armstrong MJ, et al. 
# Practice guideline update summary: Mild cognitive impairment. Neurology. 2018;90:126–135. 
# 70-74: 10.1% (95% CI 7.5%–13.5%) I2 = 5.2
# n = 2173
hrs_age_70_74_id <- hrs_imp_r3resp_wide %>%
  filter(r3age_v2 >= 70 & r3age_v2 < 75) %>% select(HHIDPN, imp)
with(hrs_age_70_74_id , table(imp, useNA = "ifany"))
mci_p = 0.101
mci_p_ll = 0.075
mci_p_ul = 0.135

# # 75–79: 14.8% (95% CI 10.1%–21.1%) I2 = 60.7
# hrs_age_75_79_id <- hrs_imp_r3resp_wide %>%
#   filter(r3age_v2 >= 75 & r3age_v2 < 80) %>% select(HHIDPN, imp)
# with(hrs_age_75_79_id , table(imp, useNA = "ifany"))

results_po_70_74 <- results_po %>%
  right_join(hrs_age_70_74_id, by = c("HHIDPN", "imp"))

# Sanity check
with(results_po_70_74, table(imp, type, useNA = "ifany"))

##---- Threshold ----
wtd_mci_threshold <-  results_po_70_74 %>%
  filter(type == "nc") %>%
  group_by(imp) %>%
  summarise(
    mci_threshold = wtd.quantile(po_w3, weights_w3, probs = mci_p, normwt = T),
    mci_threshold_ll = wtd.quantile(po_w3, weights_w3, probs = mci_p_ll, normwt = T),
    mci_threshold_ul = wtd.quantile(po_w3, weights_w3, probs = mci_p_ul, normwt = T)
  )

results_po %<>%
  left_join(wtd_mci_threshold, by = "imp") %>%
  mutate(
    across(paste0("po_w", waves_fu), function(x) case_when(x <= mci_threshold ~ 1,
                                                           x > mci_threshold ~ 0),
           .names = "{col}_mci_pe"),
    across(paste0("po_w", waves_fu), function(x) case_when(x <= mci_threshold_ll ~ 1,
                                                           x > mci_threshold_ll ~ 0),
           .names = "{col}_mci_ll"),
    across(paste0("po_w", waves_fu), function(x) case_when(x <= mci_threshold_ul ~ 1,
                                                           x > mci_threshold_ul ~ 0),
           .names = "{col}_mci_ul"))

# # Sanity check
# results_po %>%
#   select(contains("_mci"), type) %>%
#   apply(., 2, function(x) table(x, useNA = "ifany"))/nrow(results_po)

# Save the results
# saveRDS(results_po, here::here("data", "model_results", "lmtp_results_mci_risks.RDS"))

RD_RR_risk_mci <- results_po %>%
  select(type, imp, contains("_mci")) %>%
  # Use Rubin's rule, mean of mean = mean(all observations) in balanced dataset
  group_by(type) %>%
  summarise(across(contains("_mci"), 
                   ~ mean(.x))) %>%
  pivot_longer(cols = -type,
               names_to = c("wave", "estimate"),
               values_to = "mci_risk",
               names_pattern = "po_w(.*)_mci_(.*)") %>%
  mutate(wave = as.numeric(wave)) %>%
  pivot_wider(names_from = "type",
              values_from = "mci_risk") %>%
  mutate(RD_1sus_nc = `1sus` - nc,
         RR_1sus_nc = `1sus`/nc)

#---- Table and Figures ----
##---- Table ----
RD_RR_fmt_table <- RD_RR_risk_mci %>%
  pivot_wider(names_from = estimate,
              values_from = c(`1sus`, nc, RD_1sus_nc, RR_1sus_nc),
              names_glue = "{estimate}_{.value}") %>%
  mutate_if(is.numeric, ~sprintf(., fmt = '%#.2f')) %>%
  mutate(`1sus` = paste0(pe_1sus, " (", 
                         ll_1sus, ", ",
                         ul_1sus, ")"),
         nc = paste0(pe_nc, " (", 
                     ll_nc, ", ",
                     ul_nc, ")"),
         RD_1sus_nc = pe_RD_1sus_nc,
         # paste0(pe_RD_1sus_nc, " (", 
         #                   ll_RD_1sus_nc, ", ",
         #                   ul_RD_1sus_nc, ")"),
         RR_1sus_nc = pe_RR_1sus_nc
         # paste0(pe_RR_1sus_nc, " (",
         #                   ll_RR_1sus_nc, ", ",
         #                   ul_RR_1sus_nc, ")")
  ) %>%
  select(-starts_with(c("pe", "ll", "ul")))

writexl::write_xlsx(RD_RR_fmt_table, 
                    here::here("output", "tables", 
                               "sens_HRSonly_lmtp_RD_RR_SL_ahead_hrs_memimp_impute.xlsx"))
