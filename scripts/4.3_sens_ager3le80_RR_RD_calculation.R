# Estimation of RR and RD for memory impairment (Sens r3age <= 80)
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
# load(here::here("data", "analysis_data", "data_for_cluster_hrsonly_20250502.RData"))
# load(here::here("data", "model_results", "lmtp_results_SL_po_sens_ager3le80.RData"))
# waves_fu = 3:15
# 
# load(here::here("data", "analysis_data", "data_for_cluster_20250809.RData"))

load(here::here("data", "model_results", "wtd_mci_threshold_imp_2016_manly.RData"))
load(here::here("data", "model_results", "lmtp_results_SL_po_sens_ager5le80_2000_ot52l90.RData"))

waves_fu = 5:15
#---- Merge Threshold to PO ----
results_po %<>%
  left_join(wtd_mci_threshold, by = "imp") %>%
  mutate(
    across(paste0("po_w", waves_fu), function(x) case_when(x <= mci_threshold ~ 1,
                                                           x > mci_threshold ~ 0),
           .names = "{col}_mci_pe"))

# # Sanity check
# results_po %>%
#   select(contains("_mci"), type) %>%
#   apply(., 2, function(x) table(x, useNA = "ifany"))/nrow(results_po)

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
  mutate_if(is.numeric, ~sprintf(., fmt = '%#.2f'))
# mutate(`1sus` = paste0(pe_1sus, " (", 
#                        ll_1sus, ", ",
#                        ul_1sus, ")"),
#        nc = paste0(pe_nc, " (", 
#                    ll_nc, ", ",
#                    ul_nc, ")"),
#        RD_1sus_nc = pe_RD_1sus_nc,
#        # paste0(pe_RD_1sus_nc, " (", 
#        #                   ll_RD_1sus_nc, ", ",
#        #                   ul_RD_1sus_nc, ")"),
#        RR_1sus_nc = pe_RR_1sus_nc
#        # paste0(pe_RR_1sus_nc, " (",
#        #                   ll_RR_1sus_nc, ", ",
#        #                   ul_RR_1sus_nc, ")")
# ) %>%
# select(-starts_with(c("pe", "ll", "ul")))
writexl::write_xlsx(RD_RR_fmt_table, 
                    here::here("output", "tables", 
                               "sens_ager5le80_lmtp_RD_RR_SL_2000_memimp_impute_ot52l90.xlsx"))
