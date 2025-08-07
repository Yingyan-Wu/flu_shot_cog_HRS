# Estimation of RR and RD for memory impairment
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
load(here::here("data", "analysis_data", "data_for_cluster_20250731.RData"))
load(here::here("data", "model_results", "lmtp_results_SL_po_2000.RData"))
waves_fu = 5:15

#---- Categorize impairment ----
##---- Prep ----
# Table e-3 Petersen RC, Lopez O, Armstrong MJ, et al. 
# Practice guideline update summary: Mild cognitive impairment. Neurology. 2018;90:126–135. 
# 70-74: 10.1% (95% CI 7.5%–13.5%) I2 = 5.2
# n = 2173
hrs_age_70_74_id <- hrs_imp_r5resp_wide %>%
  filter(r5age_v2 >= 70 & r5age_v2 < 75) %>% select(HHIDPN, imp)
with(hrs_age_70_74_id , table(imp, useNA = "ifany"))
mci_p = 0.101
mci_p_ll = 0.075
mci_p_ul = 0.135

# # 75–79: 14.8% (95% CI 10.1%–21.1%) I2 = 60.7
# hrs_age_75_79_id <- hrs_imp_r5resp_wide %>%
#   filter(r5age_v2 >= 75 & r5age_v2 < 80) %>% select(HHIDPN, imp)
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
    mci_threshold = wtd.quantile(po_w5, weights_w5, probs = mci_p, normwt = T),
    mci_threshold_ll = wtd.quantile(po_w5, weights_w5, probs = mci_p_ll, normwt = T),
    mci_threshold_ul = wtd.quantile(po_w5, weights_w5, probs = mci_p_ul, normwt = T)
  )

# # Save the threshold
# save(wtd_mci_threshold, 
#      file = here::here("data", "model_results", "wtd_mci_threshold_imp_hrs_ahead.RData"))

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
saveRDS(results_po, here::here("data", "model_results", "lmtp_results_mci_risks_2000.RDS"))

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
                               "lmtp_RD_RR_SL_2000_memimp_impute.xlsx"))

# #---- Figure ----
# RD_RR_forplot <- RD_RR_risk_mci %>%
#   mutate(wave = as.numeric(wave)) %>%
#   select(wave, estimate, contains(c("RD", "RR"))) %>%
#   pivot_longer(cols = -c(wave, estimate),
#                names_to = "estimand",
#                names_pattern = "(.*)_1sus_nc") %>%
#   pivot_wider(names_from = estimate,
#               values_from = value)
# 
# RD_RR_forplot %>%
#   # ggplot(aes(x = survey_year, y = point_est, color = type)) +
#   ggplot(aes(x = wave, y = pe)) +
#   geom_point() +
#   geom_line() +
#   # geom_errorbar(aes(ymin = ll, ymax = ul),
#   #               width = 0.3) +
#   # scale_y_continuous(breaks = seq(-0.70, 1, by = 0.2)) +
#   scale_x_continuous(breaks = 3:15) +
#   # scale_x_continuous(breaks = seq(0, 24, by = 2)) +
#   # scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
#   facet_grid(rows = vars(estimand), scales="free") +
#   scale_color_manual(name = NULL, values = c("#F7C173", "#4B4B4B", "#3271AE"),
#                      labels = c("Continuous no flu shot", "Natural course", "Hypothetical flu shot intervention")) +# "#E5A84B",
#   labs(# title = "Memory score trajectories under the hypothetical flu shot intervention and natural course",
#     x = "Study waves",
#     y = element_blank()) +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         text = element_text(size = 12, color = "black"), 
#         strip.text.x = element_text(size = 12, color = "black"),
#         strip.text.y = element_text(size = 12, color = "black"),
#         axis.text.x = element_text(size = 12, color = "black"), 
#         axis.text.y = element_text(size = 12, color = "black"),
#         legend.text = element_text(size = 12, color = "black"))
# 
# ggsave(here::here("output", "figures", "lmtp_RD_RR_traj_SL_ahead_hrs_memimp_impute.png"), 
#        device = "png", width = 9, height = 7, units = "in", dpi = 300)

# #---- OLD ----
# Attempt using normal distribution
# #---- Read in data ----
# results_pe_se <- readRDS(here::here("data", "model_results", 
#                                     "lmtp_results_SL_rubin_forplot.RDS"))
# 
# load(here::here("data", "analysis_data", "data_for_cluster_20250304.RData"))
# 
# #---- Memory impairment prev ----
# sd(hrs_imp_r5resp_wide$r15memimp, na.rm = T)
# 
# n_samp <- 19845
# results_pe_sd_wide <- results_pe_se %>%
#   mutate(sd_rubin = se_rubin*sqrt(n_samp)) %>%
#   filter(type != "Continuous no flu shot") %>%
#   mutate(type = case_when(type == "Natural course" ~ "nc",
#                           type == "Continuous flu shot" ~ "sus1")) %>%
#   select(type, wave, survey_year, fu_year, point_est_rubin, sd_rubin) %>%
#   pivot_wider(names_from = type, 
#               values_from = c(point_est_rubin, sd_rubin),
#               names_glue = "{type}_{.value}")
# 
# 
# hrs_memimpair_prev_all <- tibble(
#   wave = 3:15,
#   # Table e-3 Petersen RC, Lopez O, Armstrong MJ, et al. 
#   # Practice guideline update summary: Mild cognitive impairment. Neurology. 2018;90:126–135. 
#   memimpair_p_nc = 0.158,
#   z_q_nc = qnorm(memimpair_p_nc)) %>%
#   left_join(results_pe_se_wide, by = "wave") %>%
#   mutate(
#     nc_sd_rubin = 1.02,
#     nc_z_mem = nc_point_est_rubin + z_q_nc*nc_sd_rubin,
#     z_q_sus1 = (nc_z_mem - sus1_point_est_rubin)/nc_sd_rubin,
#     memimpair_p_sus1 = pnorm(z_q_sus1),
#     RR = memimpair_p_sus1/memimpair_p_nc) %>%
#   select(wave, survey_year, fu_year, z_q_nc, z_q_sus1, 
#          memimpair_p_nc, memimpair_p_sus1, RR,
#          everything())
# 
# 
# #---- example normal distribution plot ----
# results_pe_se %>%
#   filter(fu_year == 12)
# # 12 years of follow up: 
# # Natural course: mem scores = 0.22 (0.02)
# # Continuous flu shot: mem scores = 0.25 (0.01)
# set.seed(62283)
# exp_norm_dist <- tibble(
#   nc = rnorm(1000, mean = 0.22, sd = 0.02),
#   sus1 = rnorm(1000, mean = 0.25, sd = 0.02))
# 
# exp_norm_dist %>%
#   ggplot(aes(x = nc)) +
#   stat_function(fun = dnorm, args = list(mean = 0.22, 
#                                          sd = 0.02)) +
#   # geom_vline(xintercept = 0.22, linetype = "dashed") +
#   stat_function(fun = dnorm, args = list(mean = 0.25,
#                                          sd = 0.02), 
#                 color = "red") +
#   scale_x_continuous(limits = c(0.1, 0.35)) +
#   # geom_vline(xintercept = 0.25, linetype = "dashed", color = "red") +
#   theme_classic()
# 
# for (w in 3:15){
#   memimpair_prev <- hrs_memimpair_prev %>%
#     filter(wave == w) %>%
#     pull(mem_impair_prev)
#   hrs_memimpair_prev[hrs_memimpair_prev$wave == w, "memimp_cutoff"] <- 
#     with(hrs_imp_r5resp_wide, quantile(get(paste0("r", w, "memimp")), 
#                                        probs = memimpair_prev, na.rm = T))
# }
# 
# hrs_memimpair_prev %>%
#   select(wave, memimp_cutoff) %>%
#   pivot_wider(names_from = wave, 
#               values_from = memimp_cutoff,
#               names_glue = "r{wave}{.value}")
# 
# hrs_imp_r5resp_wide %>%
#   
#   mutate(memimp_cutoff = quantile(hrs_imp_r5resp_wide$age_v2, 
#                                   probs = prevalence, na.rm = T)) %>%
#   
#   hrs_imp_r5resp_wide %>%
#   select(contains("memimp")) %>%
#   apply(., 2, function(x) quantile(x, probs = memimp_prev, na.rm = T))
# 
# 
# #---- Calculate RR and RD ----
# # Using memory impairment prevalence: 15.8%
# memimp_prev <- 15.8/100
# 
# 
# hrs_imp_r5resp_wide %>%
#   select(contains("memimp")) %>%
#   apply(., 2,summary)
# 
# 
# 
# 
# tibble(wave = 3:15,
#        survey_year = (wave - 5)*2 + 2000) %>% t()
# 
# results_pe_se_maxresults_pe_se_rubin <- results_pe_se %>%
#   group_by(type, wave) %>%
#   reframe(
#     point_est_rubin = mean(point_est),
#     vw = sum((se)^2)/num_imp,
#     vb = sum((point_est - point_est_rubin)^2)/(num_imp-1),
#     se_rubin = sqrt(vw + vb + vb/num_imp)) %>%
#   select(-vb, -vw) %>%
#   mutate(num_imp = num_imp,
#          ll_rubin = point_est_rubin - qnorm(0.975)*se_rubin,
#          ul_rubin = point_est_rubin + qnorm(0.975)*se_rubin,
#          survey_year = (wave - 5)*2 + 2000,
#          fu_year = survey_year - 1996)
# 
# #---- OLD ----
# ##---- RR calculation by age ----
# # Table e-3 Petersen RC, Lopez O, Armstrong MJ, et al.
# # Practice guideline update summary: Mild cognitive impairment. Neurology. 2018;90:126–135.
# memimpair_prev_tib <- tibble(
#   age_group = c("60+", "65+", "75+", "80+", "85+"),
#   mem_impair_prev = c(0.158, 0.166, 0.194, 0.286, 0.376)) %>%
#   mutate(z_nc_q = qnorm(mem_impair_prev))
# 
# hrs_imp_r5resp_wide %>%
#   select(contains("age_v2")) %>%
#   apply(., 2, summary)
# # First quantile of age
# # waves 3-6: 60-<65
# # waves 7-11: 65-<75
# # waves 12-13: 75-<80
# # waves 14-15: 80-<85
# 
# results_pe_sd_wide <- results_pe_se %>%
#   mutate(sd_rubin = se_rubin*sqrt(n_samp)) %>%
#   filter(type != "Continuous no flu shot") %>%
#   mutate(type = case_when(type == "Natural course" ~ "nc",
#                           type == "Continuous flu shot" ~ "sus1")) %>%
#   select(type, wave, survey_year, fu_year, point_est_rubin, sd_rubin) %>%
#   pivot_wider(names_from = type,
#               values_from = c(point_est_rubin, sd_rubin),
#               names_glue = "{type}_{.value}")
# 
# hrs_memimpair_prev <- tibble(
#   wave = 3:15,
#   age_group = case_when(wave %in% 3:6 ~ "60+",
#                         wave %in% 7:11 ~ "65+",
#                         wave %in% 12:13 ~ "75+",
#                         wave %in% 14:15 ~ "80+")) %>%
#   left_join(memimpair_prev_tib, by = "age_group") %>%
#   left_join(results_pe_sd_wide, by = "wave") %>%
#   mutate(nc_z_mem = nc_point_est_rubin + z_nc_q*nc_sd_rubin,
#          z_sus1_q = (nc_z_mem - sus1_point_est_rubin)/nc_sd_rubin,
#          mem_impair_prev_sus1 = pnorm(z_sus1_q),
#          RR = mem_impair_prev_sus1/mem_impair_prev) %>%
#   select(wave, age_group, survey_year, fu_year, z_nc_q, z_sus1_q,
#          mem_impair_prev, mem_impair_prev_sus1, RR,
#          everything())
# 
