# Obtain estimates from the cluster results (Sens r3age <= 80)
# local script
# Yingyan Wu
# Same calculation and process as in Script 3.2 (using HRS only Sens dataset and cluster results)

#---- Package loading + option specifications ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "haven", "labelled",
       "splines", "lmtp", "glue")

options(sci = 999)

#---- Read in the lmtp objects ----
total_scenarios <- 1:120
num_imp <- 10
waves_fu <- 5:15
results_folder <- here::here("data", "model_results", "TMLE_SuperLearner_sens_ager5le80/")
results_files <- list.files(path = results_folder, full.names = TRUE)
scenarios_res <- gsub(".*_(.*)\\..*", "\\1", results_files) %>% as.numeric() %>%
  sort()

total_scenarios[which(!total_scenarios %in% scenarios_res)] # All scenarios have results now

for (s in scenarios_res){
  load(paste0(results_folder, "tmle_SL_scenario_", s, ".RData"))
}

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
#   for (imp in 1:num_imp){
#     for (intervention in c("1sus")){
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
  lapply(seq_along(1:num_imp), 
         function(x) lmtp_pe_se(wave = 5, intervention = "nc", imp = x)) %>%
  do.call(rbind,. )

for (w in waves_fu){
  for (imp in 1:num_imp){
    tryCatch({
      results_pe_se %<>% rbind(., lmtp_pe_se(w, "1sus", imp))
    }, error = function(e) print(e))
  }
}

results_pe_se %<>%
  select(type, wave, point_est, se, imp) %>%
  mutate(type = case_when(type == "nc" ~ "Natural course",
                          type == "1sus" ~ "Continuous flu shot"))
# type == "0sus" ~ "Continuous no flu shot"))

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

contrast_scenarios <- expand_grid(w_index = waves_fu,
                                  w_ref = waves_fu,
                                  intv_index = c("1sus", "nc"), # "0sus"
                                  intv_ref = c("1sus", "nc"), # "0sus"
                                  imp = 1:num_imp) %>%
  filter(((intv_index != intv_ref & w_index == w_ref) | 
            (intv_index == intv_ref & w_index > w_ref)),
         ((intv_index == "1sus" & intv_ref == "nc") |
            # (intv_index == "1sus" & intv_ref == "0sus") |
            (intv_index == intv_ref)))

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
# . == "0sus" ~ "Continuous no flu shot")))

#---- Potential outcomes ----
# Obtain predicted potential outcomes and normalized sampling weights
# normalized sampling weights were slightly different across imputation
# Theta was calculated by weighted mean of potential outcomes from the first column
# of `outcome_reg`

load(here::here("data", "analysis_data", "data_for_cluster_ager5le80_20250731.RData"))
# PENDING
# The outcome_reg is organized by the order of participants in the input data
# tested in later chunk
HHIDPN_vec <- hrs_imp_r3resp_wide %>% filter(imp == 1) %>%
  pull(HHIDPN)

get_po <- function(wave = 3, intervention, imp){
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

results_po <- lapply(seq_along(1:num_imp), 
                     function(x) get_po(intervention = "nc", imp = x)) %>%
  do.call(rbind,. )

results_po_sus <- tibble() 
for (w in waves_fu){
  for (imp in 1:num_imp){
    tryCatch({
      results_po_sus  %<>% rbind(., get_po(w, "1sus", imp))
    }, error = function(e) print(e))
  }
}

results_po_sus_wide <- results_po_sus %>%
  pivot_wider(id_cols = c("HHIDPN", "type", "imp"), 
              names_from = "wave", values_from = c("po", "weights"),
              names_glue = "{.value}_w{wave}")

results_po %<>% rbind(results_po_sus_wide)
# Sanity check
with(results_po, table(type, imp))

#---- Save the cleaned data results ----
save(results_pe_se, results_contrasts, 
     file = here::here("data", "model_results", "lmtp_results_SL_pe_se_contrast_sens_ager3le80.RData"))
save(results_po, 
     file = here::here("data", "model_results", "lmtp_results_SL_po_sens_ager3le80.RData"))

#---- Rubin's rule ----
##---- PE CI ----
load(here::here("data", "model_results", "lmtp_results_SL_pe_se_contrast_sens_ager3le80.RData"))
results_pe_se_sens <- results_pe_se
results_contrasts_sens <- results_contrasts
num_imp = 10

results_pe_se_rubin_sens <- results_pe_se_sens %>%
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
         fu_year = survey_year - 1996)
saveRDS(results_pe_se_rubin_sens, 
        here::here("data", "model_results", "lmtp_results_SL_rubin_forplot_sens_ager3le80.RDS"))

results_contrasts_rubin_sens <- results_contrasts_sens %>%
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
         fu_year_index = survey_year_index - 1996)

#---- Table and plots ----
##---- Trajectory plot ----
results_pe_se_rubin_sens %>%
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
  scale_y_continuous(breaks = seq(-3, 0, by = 0.2)) +
  scale_x_continuous(breaks = seq(0, 24, by = 2)) +
  # scale_x_continuous(breaks = seq(1996, 2020, by = 2)) +
  scale_color_manual(name = NULL, values = c("#4B4B4B", "#3271AE"),
                     labels = c("Natural course", "Hypothetical flu shot intervention")) +# "#E5A84B",
  # scale_color_manual(name = NULL, values = c("#F7C173", "#4B4B4B", "#3271AE"),
  #                    labels = c("Continuous no flu shot", "Natural course", "Hypothetical flu shot intervention")) +# "#E5A84B",
  labs(# title = "Memory score trajectories under the hypothetical flu shot intervention and natural course",
    x = "Follow up years",
    y = "Memory z-score") +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12, color = "black"), 
        strip.text.x = element_text(size = 12, color = "black"),
        strip.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"))
ggsave(here::here("output", "figures", "sens_ager3le80_lmtp_traj_SL_ahead_hrs_memimp_impute.png"), 
       device = "png", width = 9, height = 7, units = "in", dpi = 300)

##---- Trajectory table ----
results_pe_se_year_table_sens <- results_pe_se_rubin_sens %>%
  filter(type != "Continuous no flu shot") %>%
  bind_rows(results_contrasts_rubin_sens %>% 
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

results_pe_se_year_table_sens %>%
  writexl::write_xlsx(here::here("output", "tables", 
                                 "sens_ager3le80_lmtp_pe_SL_ahead_hrs_memimp_impute.xlsx"))

##---- Contrast table ----
results_contrasts_year_table_sens <- results_contrasts_rubin_sens %>%
  mutate(across(c(intervention_index, intervention_ref),
                ~case_when(. == "Natural course" ~ "nc",
                           . == "Continuous flu shot" ~ "1sus"))) %>%
  filter(intervention_index == intervention_ref, wave_ref == 3) %>%
  select(fu_year_index, intervention_index, point_est_rubin, se_rubin) %>%
  rename(pe = point_est_rubin, se = se_rubin) %>%
  pivot_wider(names_from = intervention_index, values_from = c(pe, se), 
              names_glue = "{.value}_{intervention_index}") %>%
  mutate(pe_1sus_nc = pe_1sus - pe_nc,
         se_1sus_nc = sqrt(se_1sus^2 + se_nc^2))

for (type in c("nc", "1sus", "1sus_nc")){
  results_contrasts_year_table_sens %<>%
    mutate(
      !!paste0("ll_", type) := !!sym(paste0("pe_", type)) - qnorm(0.975)*!!sym(paste0("se_", type)),
      !!paste0("ul_", type) := !!sym(paste0("pe_", type)) + qnorm(0.975)*!!sym(paste0("se_", type))) %>%
    mutate(across(c(paste0("pe_", type), paste0("ll_", type), paste0("ul_", type)), 
                  ~sprintf(., fmt = '%#.2f'))) %>%
    mutate(!!paste0("pe_ci_", type) := 
             paste0(!!sym(paste0("pe_", type)), " (", !!sym(paste0("ll_", type)), ", ", 
                    !!sym(paste0("ul_", type)), ")"))
}

results_contrasts_year_table_sens %>%
  select(fu_year_index, contains("pe_ci_")) %>%
  writexl::write_xlsx(
    here::here("output", "tables", 
               "sens_ager3le80_lmtp_contrasts_SL_ahead_hrs_memimp_impute.xlsx"))


