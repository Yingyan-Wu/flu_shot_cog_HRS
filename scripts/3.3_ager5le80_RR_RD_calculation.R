# Estimation of RR and RD for memory impairment (r5age <= 80)
# Yingyan Wu

#---- Package loading + option specifications ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "haven", "labelled",
       "splines", "glue", "ggplot2", "Hmisc")

options(sci = 999)

#---- Load data ----
load(here::here("data", "analysis_data", "data_for_cluster_ager5le80_2000_2022.RData"))
load(here::here("data", "model_results", "lmtp_results_SL_po_ager5le80_2000_2022.RData"))
waves_fu = 5:16

#---- Load HCAP data for weighted age distribution at 2016 ----
# respondent interview data
hcap_r <- read_dta(here::here("data", "raw_data", "hcap", "hc16hp_r.dta")) %>%
  rename_all(tolower) %>%
  unite("hhidpn", c("hhid", "pn"), sep = "", remove = TRUE) %>%
  mutate_at("hhidpn", as.numeric)
summary(hcap_r$r1rage) # missing for age
with(hcap_r%>% filter(is.na(r1rage)), summary(hcap16wgtr)) # no missing

hcap_i <- read_dta(here::here("data", "raw_data", "hcap", "hc16hp_i.dta")) %>%
  rename_all(tolower) %>%
  unite("hhidpn", c("hhid", "pn"), sep = "", remove = TRUE) %>%
  mutate_at("hhidpn", as.numeric)
summary(hcap_i$inf1riwage) # missing for age

# Merge hcap_i and hcap_r for age
hcap_age <- hcap_r %>%
  select(hhidpn, r1rage, hcap16wgtr) %>%
  left_join(hcap_i %>% select(hhidpn, inf1riwage), by = "hhidpn") %>%
  mutate(r1rage = ifelse(is.na(r1rage), inf1riwage, r1rage),
         r1age_cat = case_when(
           r1rage >= 65 & r1rage < 70 ~ "[65, 70)",
           r1rage >= 70 & r1rage < 75 ~ "[70, 75)",
           r1rage >= 75 & r1rage < 80 ~ "[75, 80)",
           r1rage >= 80 & r1rage < 85 ~ "[80, 85)",
           r1rage >= 85 & r1rage < 90 ~ "[85, 90)",
           r1rage >= 90 ~ "[90, +infinity)"))

# Weighted proportion of r1age_cat
hcap_age %>%
  group_by(r1age_cat) %>%
  summarise(n = n(),
            prop = sum(hcap16wgtr)/sum(hcap_age$hcap16wgtr)) %>%
  arrange(desc(prop)) %>%
  print(n = Inf)

# r1age_cat           n   prop
# <chr>           <int>  <dbl>
#   1 [65, 70)          826 0.316 
# 2 [70, 75)          669 0.251 
# 3 [75, 80)          841 0.173 
# 4 [80, 85)          613 0.125 
# 5 [85, 90)          340 0.0752
# 6 [90, +infinity)   207 0.0598

#---- Check age distribution in year 2016 ----
hrs_imp_r5resp_wide %>% filter(imp == 1) %>%
  with(., hist(r13age_v2))

# create cat var for age groups: [65, 70), [70, 75), [75, 80), [80, 85), [85, 90), [90, +infinity)
age_mat <- hrs_imp_r5resp_wide %>% filter(imp == 1) %>%
  mutate(age_cat = case_when(
    r13age_v2 >= 65 & r13age_v2 < 70 ~ "[65, 70)",
    r13age_v2 >= 70 & r13age_v2 < 75 ~ "[70, 75)",
    r13age_v2 >= 75 & r13age_v2 < 80 ~ "[75, 80)",
    r13age_v2 >= 80 & r13age_v2 < 85 ~ "[80, 85)",
    r13age_v2 >= 85 & r13age_v2 < 90 ~ "[85, 90)",
    r13age_v2 >= 90 ~ "[90, +infinity)"))

age_mat %>%
  group_by(age_cat) %>%
  summarise(n = n(),
            prop = sum(r5wtresp)/sum(age_mat$r5wtresp)) %>%
  # prop = n()/nrow(age_mat)) %>%
  print(n = Inf)

# age_cat             n   prop
# <chr>           <int>  <dbl>
#   1 [65, 70)          822 0.0760
# 2 [70, 75)         2846 0.235 
# 3 [75, 80)         4376 0.209 
# 4 [80, 85)         3777 0.175 
# 5 [85, 90)         2862 0.152 
# 6 [90, +infinity)  2973 0.153 

# age 75-80 has the proportion most similar, use age 75-80

#---- Categorize impairment ----
# # Table 2 
# Manly JJ, Jones RN, Langa KM, et al. 
# Estimating the Prevalence of Dementia and Mild Cognitive Impairment in the US: 
# The 2016 Health and Retirement Study Harmonized Cognitive Assessment Protocol Project. 
# JAMA Neurol. 2022;79(12):1242-1249. doi:10.1001/jamaneurol.2022.3543

# # 75-79: 9% (6-11) for dementia and 21% (18-24) for MCI
hrs_age_75_79_id <- hrs_imp_r5resp_wide %>%
  filter(r13age_v2 >= 75 & r13age_v2 <= 79) %>% select(HHIDPN, imp)
with(hrs_age_75_79_id, table(imp, useNA = "ifany"))
mci_p = 0.3

##---- Prep ----
results_po_75_79 <- results_po %>%
  right_join(hrs_age_75_79_id, by = c("HHIDPN", "imp"))

# Sanity check
with(results_po_75_79, table(imp, type, useNA = "ifany"))

##---- Threshold ----
wtd_mci_threshold <-  results_po_75_79 %>%
  filter(type == "nc") %>%
  group_by(imp) %>%
  summarise(
    mci_threshold = wtd.quantile(po_w13, weights_w13, probs = mci_p, normwt = T))

# Save the threshold
save(wtd_mci_threshold,
     file = here::here("data", "model_results", "wtd_mci_threshold_imp_2016_manly_mainle80.RData"))

results_po %<>%
  left_join(wtd_mci_threshold, by = "imp") %>%
  mutate(across(paste0("po_w", waves_fu), function(x) case_when(x <= mci_threshold ~ 1,
                                                                x > mci_threshold ~ 0),
                .names = "{col}_mci_pe"))

# # Sanity check
# results_po %>%
#   select(contains("_mci"), type) %>%
#   apply(., 2, function(x) table(x, useNA = "ifany"))/nrow(results_po)

table(results_po$po_w13_mci_pe)

# Save the results
saveRDS(results_po, here::here("data", "model_results", "lmtp_results_mci_risks_2000_2020_ot52le80.RDS"))

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
                               "ager5le80_lmtp_RD_RR_SL_2000_2022_memimp_impute.xlsx"))
