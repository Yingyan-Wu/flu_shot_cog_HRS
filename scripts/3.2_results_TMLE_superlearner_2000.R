# Obtain estimates from the cluster results
# local script
# Yingyan Wu

#---- Package loading + option specifications ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "haven", "labelled",
       "splines", "lmtp", "glue")

options(sci = 99999)

#---- Read in the lmtp objects ----
total_scenarios <- 1:120
num_imp <- 10
waves_fu <- 5:15
results_folder <- here::here("data", "model_results", "TMLE_SuperLearner/")
results_files <- list.files(path = results_folder, full.names = TRUE)
scenarios_res <- gsub(".*_(.*)\\..*", "\\1", results_files) %>% as.numeric() %>%
  sort()

total_scenarios[which(!total_scenarios %in% scenarios_res)]
# All scenarios have results now

for (s in scenarios_res){
  load(paste0(results_folder, "tmle_SL_scenario_", s, ".RData"))
}

# PENDING: Work on imputations with full results
full_imps <- c(2, 5:10)
full_waves <- 5:15

# lmtp_tmle function returned these values for every participant
# 1. eif
# 2. id
# 3. weights: normalized weights (the weights fed into the weights argument, sampling weights
# in our case)
# 4. outcome_reg: An n x Tau + 1 matrix of outcome regression predictions. 
# The (weighted) mean of the first column is used for calculating theta (PE) using
# `weights` (normalized weights)
# 5. density_ratios: An n x Tau matrix of the estimated, non-cumulative, density ratios.

# #---- Fits ----
# for (wave in waves_fu){
#   for (imp in full_imps){
#     for (intervention in c("1sus", "0sus")){
#       tryCatch({
#         print(get(paste0("y_w", wave, "_", intervention, "_", imp))$fits_m)
#         print(get(paste0("y_w", wave, "_", intervention, "_", imp))$fits_r)
#       }, error = function(e) print(e))
#     }
#   }
# }

#---- Estimates ----
lmtp_pe_se <- function(wave = 5, intervention, imp){
  if (intervention == "nc"){
    # All waves were in one list
    lmtp_list <- get(paste0("lmtp_nc_results_imp_", imp))
    return(sapply(lmtp_list, function(x) c(x$theta, x$standard_error)) %>%
             t() %>%
             as.data.frame() %>%
             set_colnames(c("point_est", "se")) %>%
             mutate(type = intervention,
                    wave = waves_fu,
                    imp = imp))
  } else {
    # waves were separated from each object
    lmtp_object <- get(paste0("y_w", wave, "_", intervention, "_", imp))
    return(
      c(lmtp_object$theta, lmtp_object$standard_error) %>%
        t() %>%
        as.data.frame() %>%
        magrittr::set_colnames(c("point_est", "se")) %>%
        mutate(type = intervention,
               wave = w,
               imp = imp))
  }
}

results_pe_se <- 
  lapply(seq_along(full_imps), 
         function(x) lmtp_pe_se(intervention = "nc", imp = x)) %>%
  do.call(rbind,. )

for (w in full_waves){
  for (imp in full_imps){
    tryCatch({
      results_pe_se %<>% rbind(., lmtp_pe_se(w, "1sus", imp))
    }, error = function(e) print(e))
  }
}

results_pe_se %<>%
  select(type, wave, point_est, se, imp) %>%
  mutate(type = case_when(type == "nc" ~ "Natural course",
                          type == "1sus" ~ "Continuous flu shot"))

#---- Contrasts ----
# lmtp_contrast(y_w10_0sus_1, ref = y_w10_1sus_1)

lmtp_res_contrast <- function(w_index, w_ref, intv_index, intv_ref, imp){
  w_num_index <-  which(waves_fu == w_index)
  w_num_ref <- which(waves_fu == w_ref)
  
  if(intv_index == "nc"){
    lmtp_index_object <- get(paste0("lmtp_nc_results_imp_", imp))[[w_num_index]]
  } else {
    lmtp_index_object <- get(paste0("y_w", w_index, "_", intv_index, "_", imp))
  }
  
  if (intv_ref == "nc"){
    lmtp_ref_object <- get(paste0("lmtp_nc_results_imp_", imp))[[w_num_ref]]
  } else {
    lmtp_ref_object <- get(paste0("y_w", w_ref, "_", intv_ref, "_", imp))
  }
  
  contrast <- lmtp_contrast(lmtp_index_object, ref = lmtp_ref_object)
  return(contrast$vals %>% mutate(wave_index = w_index, wave_ref = w_ref, 
                                  intervention_index = intv_index, 
                                  intervention_ref = intv_ref, imp = imp))
}

contrast_scenarios <- expand_grid(w_index = full_waves,
                                  w_ref = full_waves,
                                  intv_index = c("1sus", "nc"),
                                  intv_ref = c("1sus", "nc"),
                                  imp = full_imps) %>%
  filter(((intv_index != intv_ref & w_index == w_ref) | 
            (intv_index == intv_ref & w_index > w_ref)),
         ((intv_index == "1sus" & intv_ref == "nc") |
            (intv_index == intv_ref)))

# Pending
# # Due to scale differences after digit 10 (mbp had slight difference for scaling the weights)
# table(round(y_w13_1sus_1$weights, 9) == round(y_w5_1sus_1$weights, 9))
# table(round(y_w14_1sus_1$weights, 9) == round(y_w5_1sus_1$weights, 9))
# table(round(y_w15_1sus_1$weights, 9) == round(y_w5_1sus_1$weights, 9))
# y_w13_1sus_1$weights <- y_w5_1sus_1$weights
# y_w14_1sus_1$weights <- y_w5_1sus_1$weights
# y_w15_1sus_1$weights <- y_w5_1sus_1$weights

raw_results_contrasts <- mapply(lmtp_res_contrast, 
                                contrast_scenarios$w_index, 
                                contrast_scenarios$w_ref, 
                                contrast_scenarios$intv_index, 
                                contrast_scenarios$intv_ref, 
                                contrast_scenarios$imp) 

results_contrasts <- raw_results_contrasts %>% 
  t() %>% as_tibble() %>%
  mutate_all(unlist) %>%
  select(intervention_index, intervention_ref, wave_index, wave_ref, 
         theta, std.error, imp) %>%
  rename(point_est = theta, se = std.error) %>%
  mutate(across(c(intervention_index, intervention_ref),
                ~case_when(. == "nc" ~ "Natural course",
                           . == "1sus" ~ "Continuous flu shot")))

results_contrasts_memory_change <- results_contrasts %>%
  mutate(across(c(intervention_index, intervention_ref),
                ~case_when(. == "Natural course" ~ "nc",
                           . == "Continuous flu shot" ~ "1sus"))) %>%
  filter(intervention_index == intervention_ref, wave_ref == 5) %>%
  select(imp, wave_index, intervention_index, pe = point_est, se) %>%
  pivot_wider(names_from = intervention_index, values_from = c(pe, se), 
              names_glue = "{.value}_{intervention_index}") %>%
  mutate(pe_1susnc = pe_1sus - pe_nc,
         se_1susnc = sqrt(se_1sus^2 + se_nc^2))

#---- Potential outcomes ----
# Obtain predicted potential outcomes and normalized sampling weights
# normalized sampling weights were slightly different across imputation
# Theta was calculated by weighted mean of potential outcomes from the first column
# of `outcome_reg`

# # Test
# table(lmtp_nc_results_imp_1[[1]]$weights == lmtp_nc_results_imp_1[[11]]$weights)
# table(round(y_w3_1sus_1$weights, 1) == round(y_w3_1sus_3$weights, 1))
# hist(y_w3_1sus_1$weights)
# hist(y_w3_1sus_3$weights)

load(here::here("data", "analysis_data", "data_for_cluster_20250809.RData"))

# The outcome_reg is organized by the order of participants in the input data
# tested in later chunk
HHIDPN_vec <- hrs_imp_r5resp_wide %>% filter(imp == 1) %>%
  pull(HHIDPN)

get_po <- function(wave = 5, intervention, imp){
  if (intervention == "nc"){
    lmtp_list <- get(paste0("lmtp_nc_results_imp_", imp))
    po_dt <- lapply(lmtp_list, function(x) tibble(x$outcome_reg[, 1], x$weights)) %>%
      do.call(cbind, .) %>%
      set_colnames(as.vector(outer(c("po", "weights"), 
                                   paste0("w", waves_fu), paste, sep = "_"))) %>%
      mutate(HHIDPN = HHIDPN_vec,
             type = intervention,
             imp = imp)
    return(po_dt)
  } else {
    lmtp_object <- get(paste0("y_w", wave, "_", intervention, "_", imp))
    po_dt_w <- tibble(lmtp_object$outcome_reg[, 1], lmtp_object$weights) %>%
      set_colnames(c("po", "weights")) %>%
      mutate(HHIDPN = HHIDPN_vec,
             type = intervention,
             wave = wave,
             imp = imp)
    return(po_dt_w)
  }
}

results_po <- lapply(full_imps, 
                     function(x) get_po(intervention = "nc", imp = x)) %>%
  do.call(rbind,. )

results_po_sus <- tibble() 
for (w in full_waves){
  for (imp in full_imps){
    tryCatch({
      results_po_sus  %<>% rbind(., get_po(w, "1sus", imp))
    }, error = function(e) print(e))
  }
}

results_po_sus_wide <- results_po_sus %>%
  pivot_wider(id_cols = c("HHIDPN", "type", "imp"), 
              names_from = "wave", values_from = c("po", "weights"),
              names_glue = "{.value}_w{wave}")

# results_po %<>% rbind(results_po_sus_wide)

results_po %<>% 
  # select(-contains(paste0("w", waves_fu[!waves_fu %in% full_waves]))) %>%
  rbind(results_po_sus_wide)

# Sanity check
with(results_po, table(type, imp))

# # Test: PO is sorted by the order of input data
# hrs_imp1 <- hrs_imp_r5resp_wide %>% filter(imp == 1) %>%
#   select(HHIDPN, r5age_v2, r5memimp, r5wtresp) %>%
#   mutate(r5wtresp_norm = r5wtresp/mean(r5wtresp))
# results_po_imp1 <- results_po %>%
#   filter(imp == 1, type == "nc")
# weights_w3_mat <- hrs_imp1 %>%
#   select(HHIDPN, r5wtresp_norm) %>%
#   cbind(results_po_imp1 %>%
#           select(weights_w3)) %>%
#   mutate(check = case_when(round(r5wtresp_norm, 6) == round(weights_w3, 6) ~ 1,
#                            TRUE ~ 0))
# with(weights_w3_mat, table(check, useNA = "ifany"))

# #---- Save the cleaned data results ----
save(results_pe_se, results_contrasts, results_contrasts_memory_change,
     file = here::here("data", "model_results", "lmtp_results_SL_pe_se_contrast_2000_ot52.RData"))
save(results_po,
     file = here::here("data", "model_results", "lmtp_results_SL_po_2000_ot52.RData"))

#---- Rubin's rule ----
##---- PE CI ----
load(here::here("data", "model_results", "lmtp_results_SL_pe_se_contrast_2000_ot52.RData"))
num_imp <- length(full_imps)
results_pe_se_rubin <- results_pe_se %>%
  group_by(type, wave) %>%
  reframe(
    point_est_rubin = mean(point_est),
    vw = sum((se)^2)/num_imp,
    vb = sum((point_est - point_est_rubin)^2)/(num_imp-1),
    se_rubin = sqrt(vw + vb + vb/num_imp)) %>%
  select(-vb, -vw) %>%
  mutate(num_imp = num_imp,
         ll_rubin = point_est_rubin - qnorm(0.975)*se_rubin,
         ul_rubin = point_est_rubin + qnorm(0.975)*se_rubin,
         survey_year = (wave - 5)*2 + 2000,
         fu_year = survey_year - 2000)
saveRDS(results_pe_se_rubin,
        here::here("data", "model_results", "lmtp_results_SL_rubin_forplot_2000_ot52.RDS"))

results_contrasts_rubin <- results_contrasts %>%
  group_by(intervention_index, intervention_ref, wave_index, wave_ref) %>%
  reframe(
    point_est_rubin = mean(point_est),
    vw = sum((se)^2)/num_imp,
    vb = sum((point_est - point_est_rubin)^2)/(num_imp-1),
    se_rubin = sqrt(vw + vb + vb/num_imp)) %>%
  select(-vb, -vw) %>%
  mutate(num_imp = num_imp,
         ll_rubin = point_est_rubin - qnorm(0.975)*se_rubin,
         ul_rubin = point_est_rubin + qnorm(0.975)*se_rubin,
         survey_year_index = (wave_index - 5)*2 + 2000,
         fu_year_index = survey_year_index - 2000)

results_contrasts_memory_change_rubin <- results_contrasts_memory_change %>%
  pivot_longer(cols = -c(imp, wave_index),
               names_to = c(".value", "intervention"),
               names_pattern = "(.*)_(.*)") %>%
  group_by(wave_index, intervention) %>%
  reframe(
    pe_rubin = mean(pe),
    vw = sum((se)^2)/num_imp,
    vb = sum((pe - pe_rubin)^2)/(num_imp-1),
    se_rubin = sqrt(vw + vb + vb/num_imp)) %>%
  select(-vb, -vw) %>%
  mutate(num_imp = num_imp,
         ll_rubin = pe_rubin - qnorm(0.975)*se_rubin,
         ul_rubin = pe_rubin + qnorm(0.975)*se_rubin,
         survey_year_index = (wave_index - 5)*2 + 2000,
         fu_year_index = survey_year_index - 2000)


#---- Table and plots ----
Hmisc::wtd.mean(hrs_imp_r5resp_wide$r5memimp_stdz, 
                weights = hrs_imp_r5resp_wide$r5wtresp, na.rm= T)
# 0.0000000000000004949333
mean(hrs_imp_r5resp_wide$r5memimp_stdz) #-0.1235858

##---- Trajectory plot ----
results_pe_se_rubin <-
  readRDS(here::here("data", "model_results", "lmtp_results_SL_rubin_forplot_2000_ot52.RDS"))
results_pe_se_rubin %>%
  filter(type != "Continuous no flu shot") %>%
  mutate(type = factor(type, levels = c("Natural course", "Continuous flu shot"))) %>%
  # mutate(type = factor(type, levels = c("Continuous no flu shot", "Natural course", "Continuous flu shot"))) %>%
  # ggplot(aes(x = survey_year, y = point_est, color = type)) +
  ggplot(aes(x = fu_year, y = point_est_rubin, color = type)) +
  geom_point(position = position_dodge((width = 1))) +
  geom_line(position = position_dodge((width = 1))) +
  geom_errorbar(aes(ymin = ll_rubin, ymax = ul_rubin),
                position = position_dodge((width = 1)),
                width = 0.8) +
  # # SER version
  # geom_point(size = 2, position = position_dodge((width = 1))) +
  # geom_line(size = 1, position = position_dodge((width = 1))) +
  # geom_errorbar(aes(ymin = ll_rubin, ymax = ul_rubin),
  #               size = 1,
  #               position = position_dodge((width = 1)),
  #               width = 1) +
  scale_y_continuous(breaks = seq(-3, 0, by = 0.2)) +
  scale_x_continuous(breaks = seq(0, 24, by = 2)) +
  # scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
  scale_color_manual(name = NULL, values = c("#4B4B4B", "#3271AE"),
                     labels = c("Natural course", "Hypothetical flu shot intervention")) +# "#E5A84B",
  # scale_color_manual(name = NULL, values = c("#F7C173", "#4B4B4B", "#3271AE"),
  #                    labels = c("Continuous no flu shot", "Natural course", "Hypothetical flu shot intervention")) +# "#E5A84B",
  labs(# title = "Memory score trajectories under the hypothetical flu shot intervention and natural course",
    x = "Follow-up years",
    y = "Memory z-score") +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12, color = "black"),
        strip.text.x = element_text(size = 12, color = "black"),
        strip.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"))
# # SER version
# theme(legend.position = "bottom",
#       text = element_text(size = 16, color = "black"),
#       strip.text.x = element_text(size = 16, color = "black"),
#       strip.text.y = element_text(size = 16, color = "black"),
#       axis.text.x = element_text(size = 16, color = "black"),
#       axis.text.y = element_text(size = 16, color = "black"),
#       legend.text = element_text(size = 16, color = "black"),
#       strip.text = element_text(size = 16, color = "black"))
# ggsave(here::here("output", "figures", "lmtp_traj_SL_ahead_hrs_memimp_impute_SER.png"),
#        device = "png", width = 9, height = 7, units = "in", dpi = 600)
ggsave(here::here("output", "figures", "lmtp_traj_SL_2000_memimp_impute.png"),
       device = "png", width = 9, height = 7, units = "in", dpi = 300)

##---- Trajectory table ----
results_pe_se_year_table <- results_pe_se_rubin %>%
  filter(type != "Continuous no flu shot") %>%
  bind_rows(results_contrasts_rubin %>% 
              filter(intervention_index == "Continuous flu shot",
                     intervention_ref == "Natural course",
                     wave_index == wave_ref) %>%
              mutate(type = "Mean difference") %>%
              select(survey_year = survey_year_index, 
                     fu_year = fu_year_index,
                     type,
                     point_est_rubin, ll_rubin, ul_rubin)) %>%
  mutate(
    type = factor(type, levels = 
                    c("Natural course", "Continuous flu shot", "Mean difference")),
    across(point_est_rubin:ul_rubin, 
           ~sprintf(., fmt = '%#.2f')),
    pe_ci = paste0(point_est_rubin, " (", ll_rubin, ", ", ul_rubin, ")")) %>%
  arrange(type) %>%
  select(survey_year, fu_year, type, pe_ci) %>%
  pivot_wider(names_from = type, values_from = pe_ci, 
              names_glue = "{type} (95% CI)")

results_pe_se_year_table %>%
  writexl::write_xlsx(here::here("output", "tables",
                                 "lmtp_pe_SL_2000_memimp_impute.xlsx"))

##---- Contrast table ----
results_contrasts_year_memory_change_table <- results_contrasts_memory_change_rubin %>%
  mutate(across(c("pe_rubin", "ll_rubin", "ul_rubin"), 
                ~sprintf(., fmt = '%#.2f'))) %>%
  mutate(pe_ci = paste0(pe_rubin, " (", ll_rubin, ", ", ul_rubin, ")")) %>%
  select(survey_year_index, fu_year_index, intervention, pe_ci) %>%
  pivot_wider(id_cols = c(survey_year_index, fu_year_index), 
              names_from = intervention, values_from = "pe_ci",
              names_glue = "{.value}_{intervention}") %>%
  select(fu_year_index, pe_ci_nc, pe_ci_1sus, pe_ci_1susnc)

results_contrasts_year_memory_change_table %>%
  set_colnames(c("Follow up year", 
                 "Natural course (95% CI)", "Sustained flu shot (95% CI)", 
                 "Sustained - natural course (95% CI)")) %>%
  writexl::write_xlsx(
    here::here("output", "tables",
               "lmtp_contrasts_SL_2000_memimp_impute.xlsx"))

#---- Comparing observed memory trajectories ----
# Compare based on pooled results
hrs_imp_r5resp_memimp_long <- hrs_imp_r5resp_wide %>%
  select(HHIDPN, imp, r5wtresp, contains("memimp_stdz")) %>%
  pivot_longer(cols = -c(HHIDPN, imp, r5wtresp), names_to = "wave", 
               values_to = "memimp_stdz",
               names_pattern = "r(.*)memimp_stdz") %>%
  mutate(wave = as.numeric(wave)) %>%
  na.omit()

with(hrs_imp_r5resp_memimp_long, table(wave, imp, useNA = "ifany"))
# Size equals across imputations at different waves

nc_observed_comp_tib <- hrs_imp_r5resp_memimp_long %>%
  filter(wave > 4) %>%
  select(wave, memimp_stdz, r5wtresp) %>%
  group_by(wave) %>%
  dplyr::summarize(n = n(),
                   mean_memimpstdz = Hmisc::wtd.mean(memimp_stdz, weights = r5wtresp,
                                                     na.rm = TRUE),
                   se_memimpstdz = sqrt(Hmisc::wtd.var(memimp_stdz, weights = r5wtresp,
                                                       na.rm = TRUE))/sqrt(n)) %>%
  ungroup() %>%
  mutate(fu_yr = (wave - 5) * 2)

results_pe_se_rubin <- results_pe_se %>%
  group_by(type, wave) %>%
  reframe(
    point_est_rubin = mean(point_est),
    vw = sum((se)^2)/num_imp,
    vb = sum((point_est - point_est_rubin)^2)/(num_imp-1),
    se_rubin = sqrt(vw + vb + vb/num_imp)) %>%
  select(-vb, -vw) %>%
  mutate(num_imp = num_imp,
         ll_rubin = point_est_rubin - qnorm(0.975)*se_rubin,
         ul_rubin = point_est_rubin + qnorm(0.975)*se_rubin,
         survey_year = (wave - 5)*2 + 2000,
         fu_year = survey_year - 2000)

# From the lmtp no ipcw nc results
load(here::here("data", "model_results", "lmtp_nc_noipcw_results_2000_ot52.RData"))
nc_noipcw_imp1 <- nc_results %>%
  filter(imp == 1) %>%
  rename(mean_ponc = point_est,
         se_ponc = se)

nc_observed_comp_tib %<>%
  left_join(nc_noipcw_imp1, by = "wave")
# No IPCW NC PO is the same as observed mean!!!

nc_observed_comp_tib %>%
  pivot_longer(cols = -c(wave, n, fu_yr),
               names_to = c(".value", "type"),
               names_pattern = "(.*)_(.*)") %>%
  filter(type == "memimpstdz") %>%
  rbind(., nc_observed_comp_tib %>%
          pivot_longer(cols = -c(wave, n, fu_yr),
                       names_to = c(".value", "type"),
                       names_pattern = "(.*)_(.*)") %>%
          filter(type == "memimpstdz") %>% mutate(type = "po")) %>%
  ggplot(aes(x = fu_yr, y = mean, group = type, color = type)) +
  geom_point(position = position_dodge((width = 0.5))) +
  geom_line(aes(linetype = type), position = position_dodge((width = 0.5))) +
  geom_errorbar(aes(ymin = mean - qnorm(0.975) * se, 
                    ymax = mean + qnorm(0.975) * se),
                position = position_dodge((width = 0.5)),
                width = 0.8) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  # scale_y_continuous(limits = c(-0.8, 0.2), breaks = seq(-0.8, 0.2, by = 0.2)) +
  scale_linetype_manual(name = NULL, values = c("dashed", "solid"),
                        labels = c("Observed memory z-score", "Natural course Potential outcomes")) +
  scale_color_manual(name = NULL, values = c("#74B72E", "#3A5311"),
                     labels = c("Observed memory z-score", "Natural course Potential outcomes")) +
  theme_bw() +
  labs(x = "Follow-up year",
       y = "Memory z-score",
       color = "Type") +
  theme(legend.position = "bottom")

ggsave(here::here("output", "figures", "lmtp_nc_observed_comp_SL_2000_memimp_impute.png"),
       device = "png", width = 9, height = 7, units = "in", dpi = 300)
