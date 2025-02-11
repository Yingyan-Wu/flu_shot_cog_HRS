# Preliminary analysis: no time varying. baseline confounders and overall flusht acceptance with repeated measure of outcome
# Created by: Yingyan Wu
# Date: Aug.19.2024

#---- Load packages ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "haven", "labelled",
       "survey", "geepack", "splines")
options(scipen = 999)

#---- function ----
ns_1knot <- function(x){
  temp <- quantile(x, c(.05, .5, .95))
  splines::ns(x, Boundary.knots = temp[c(1, 3)], knots= temp[2])
}

#---- Multiple imputed ----
hrs_imp <- readRDS(here::here("data", "analysis_data", "imp_pmm_stacked.RDS"))

#---- Flu shot ----
# Num of flushot
hrs_imp %<>%
  group_by(HHIDPN, imp) %>%
  mutate(num_flusht = sum(flusht, na.rm = T)) %>%
  ungroup()  
# Sanity check
with(hrs_imp, table(num_flusht, useNA = "ifany"))
with(hrs_imp %>% filter(imp == 1), table(num_flusht, useNA = "ifany"))

# Flusht wide
hrs_imp_flusht_wide <- hrs_imp %>%
  select(HHIDPN, imp, wave, flusht) %>%
  pivot_wider(names_from = wave, values_from = flusht, 
              names_glue = "r{wave}{.value}")
hrs_imp %<>%
  left_join(hrs_imp_flusht_wide, by = c("HHIDPN", "imp"))

# # Sanity check
# with(hrs_imp, table(r3flusht, useNA = "ifany"))
# with(hrs_imp_r3flusht, table(r3flusht, useNA = "ifany"))

# First 3 flu sht
hrs_imp %<>%
  mutate(num_flusht_3_7 = 
           case_when(
             rowSums(!is.na(across(paste0("r", seq(3, 7, by = 2), "flusht")))) == 0 ~ NA_real_,
             TRUE ~ rowSums(across(paste0("r", seq(3, 7, by = 2), "flusht")), na.rm = TRUE)))

#---- Filter data ----
# Filter on participants eligible for 1996 interview and not died, ltfu throughout
r3_resp_id <- hrs_imp %>% select(HHIDPN, wave, iwstat) %>%
  filter(wave == 3, iwstat %in% c(1, 4)) %>% pull(HHIDPN) %>% unique()
hrs_imp_r3resp <- hrs_imp %>% filter(HHIDPN %in% r3_resp_id, iwstat %in% c(1, 4))

hrs_imp_r3resp_imp1 <- hrs_imp_r3resp %>% filter(imp == 1)

#---- Empirical plot ----
# imputation 1 overall flusht
unadj_mat_mean_imp_balanced <- hrs_imp_r3resp_imp1 %>%
  select(num_flusht, wave, memimp) %>%
  na.omit() %>%
  group_by(num_flusht, wave) %>%
  summarise_at(vars(memimp), 
               list(n = length, mean = mean, sd = sd)) %>%
  mutate(se = sd/sqrt(n), 
         num_flusht = as.factor(num_flusht),
         lower = mean - 1.96*se,
         upper = mean + 1.96*se)
n_balanced <- hrs_imp_r3resp_imp1 %>%
  select(HHIDPN) %>% unique() %>% nrow()

unadj_mat_mean_imp_balanced %>%
  ggplot(aes(x = wave, y = mean, color = num_flusht, group = num_flusht)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = 3:15) +
  labs(title = paste0("Balanced model participants (wave 3 - 15) memory scores (observed), n = ", n_balanced))

# imputation 1 baseline flusht
unadj_mat_mean_imp_balanced_bl <- hrs_imp_r3resp_imp1 %>%
  select(r3flusht, wave, memimp) %>%
  na.omit() %>%
  group_by(r3flusht, wave) %>%
  summarise_at(vars(memimp), 
               list(n = length, mean = mean, sd = sd)) %>%
  mutate(se = sd/sqrt(n), 
         r3flusht = as.factor(r3flusht),
         lower = mean - 1.96*se,
         upper = mean + 1.96*se)
n_balanced <- hrs_imp_r3resp_imp1 %>%
  select(HHIDPN) %>% unique() %>% nrow()

unadj_mat_mean_imp_balanced_bl %>%
  ggplot(aes(x = wave, y = mean, color = r3flusht, group = r3flusht)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = 3:15) +
  labs(title = paste0("Balanced model participants (wave 3 - 15) memory scores (observed), n = ", n_balanced))


# imputation 1 first 3 flusht
unadj_mat_mean_imp_balanced_3_7 <- hrs_imp_r3resp_imp1 %>%
  select(num_flusht_3_7, wave, memimp) %>%
  na.omit() %>%
  group_by(num_flusht_3_7, wave) %>%
  summarise_at(vars(memimp), 
               list(n = length, mean = mean, sd = sd)) %>%
  mutate(se = sd/sqrt(n), 
         num_flusht_3_7 = as.factor(num_flusht_3_7),
         lower = mean - 1.96*se,
         upper = mean + 1.96*se)
n_balanced <- hrs_imp_r3resp_imp1 %>%
  select(HHIDPN) %>% unique() %>% nrow()

unadj_mat_mean_imp_balanced_3_7 %>%
  ggplot(aes(x = wave, y = mean, color = num_flusht_3_7, group = num_flusht_3_7)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = 3:15) +
  labs(title = paste0("Balanced model participants (wave 3 - 15) memory scores (observed), n = ", n_balanced))

#---- Preliminary GEE ----
long_gcomp_gee <- geepack::geeglm(
  memimp ~ as.factor(flusht) + as.factor(flusht) * ns_1knot(fu_year) + 
    (r3age_v2 + female + raceeth + raedyrs +
       mstat_2cat + income_pp + insurance_2cat + smoken + cesd + flone +
       hibpe + diabe + cancre + stroke + hearte + lunge + adl_2cat + iadl_2cat) * ns_1knot(fu_year),
  data = hrs_imp_r3resp_imp1,
  id = HHIDPN, waves = wave,
  weights = wtresp, family = gaussian)

broom::tidy(long_gcomp_gee) # %>% view()

hrs_imp_r3resp_imp1 %<>%
  mutate(memimp_pred = predict(long_gcomp_gee, newdata = hrs_imp_r3resp_imp1, type = "response"))

pred_mat_mean <- hrs_imp_r3resp_imp1 %>%
  select(flusht, wave, memimp_pred) %>%
  group_by(flusht, wave) %>%
  summarise_at(vars(memimp_pred), 
               list(n = length, mean = mean, sd = sd)) %>%
  mutate(se = sd/sqrt(n), 
         flusht_sts = as.factor(flusht),
         lower = mean - 1.96*se,
         upper = mean + 1.96*se)

pred_mat_mean %>%
  ggplot(aes(x = wave, y = mean, color = flusht_sts, group = flusht_sts)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8) +
  geom_line() +
  theme_bw()

#---- Preliminary linear regression (wave 3 memory score) ----
qmodel_lr <- glm(
  memimp ~ as.factor(flusht) + 
    (r3age_v2 + female + raceeth + raedyrs +
       mstat_2cat + income_pp + insurance_2cat + smoken + cesd + flone +
       hibpe + diabe + cancre + stroke + hearte + lunge + adl_2cat + iadl_2cat),
  data = hrs_imp_r3resp_imp1 %>% filter(wave == 3), family = gaussian)
summary(qmodel_lr)


#---- Original data -----
#---- Load the data ----
hrs_selected <- readRDS(here::here("data", "analysis_data", "hrs_selected.rds")) %>%
  dplyr::mutate(outcome_balanced_3_15 = case_when(rowSums(is.na(across(contains("memimp")))) == 0 ~ 1,
                                                  TRUE ~ 0),
                outcome_balanced_5_15 = case_when(rowSums(is.na(across(paste0("r", 5:15, "memimp")))) == 0 ~ 1,
                                                  TRUE ~ 0),
                outcome_balanced_7_15 = case_when(rowSums(is.na(across(paste0("r", 7:15, "memimp")))) == 0 ~ 1,
                                                  TRUE ~ 0),
                outcome_balanced_10_15 = case_when(rowSums(is.na(across(paste0("r", 10:15, "memimp")))) == 0 ~ 1,
                                                   TRUE ~ 0))

# Sanity check
with(hrs_selected, table(outcome_balanced_3_15))
with(hrs_selected, table(outcome_balanced_5_15))



#---- preliminary data ----
hrs_flusht_cohort1 <- hrs_selected %>%
  filter(flusht_960004 == 1) %>%
  mutate(flusht_sts = case_when(r3flusht == 1 & r5flusht == 1 & r7flusht == 1 ~ 2,
                                r3flusht == 1 | r5flusht == 1 | r7flusht == 1 ~ 1,
                                r3flusht == 0 & r5flusht == 0 & r7flusht == 0 ~ 0))

hrs_flusht_cohort1_long <- hrs_flusht_cohort1 %>%
  pivot_longer(cols = contains(c("memimp", "dead", "inw")),
               names_to = c("wave", ".value"), 
               names_pattern = "r(\\d+)(.*)") %>%
  mutate(wave = as.numeric(wave),
         fu_year = (wave - 3)*2)

hrs_flusht_cohort1_long_cc <- hrs_flusht_cohort1 %>%
  # Conditioning on survival after the end of the exposure period (r8, 2006)
  filter(r8dead == 0) %>%
  pivot_longer(cols = contains(c("memimp", "dead", "inw")),
               names_to = c("wave", ".value"), 
               names_pattern = "r(\\d+)(.*)") %>%
  mutate(wave = as.numeric(wave),
         fu_year = (wave - 3)*2) %>%
  select(HHIDPN, wave, inw, dead, fu_year, flusht_sts, memimp, r3wtresp, r3age_v2, female,
         raceeth, raedyrs, r3mstat_2cat, r3income_pp, r3insurance_2cat, r3smoken,
         r3cesd, r3flone, r3hibpe, r3diabe, r3cancre, r3stroke, r3hearte, r3lunge,
         r3adl_2cat, r3iadl_2cat) %>%
  na.omit()

# Sanity check
with(hrs_flusht_cohort1_long_cc, table(wave, fu_year, useNA = "ifany"))
summary(hrs_flusht_cohort1_long_cc$memimp)

onetime_gcomp_gee <- geepack::geeglm(
  memimp ~ as.factor(flusht_sts) + as.factor(flusht_sts) * ns_1knot(fu_year) + 
    # first_memimp +
    (r3age_v2 + female + raceeth + raedyrs +
       r3mstat_2cat + r3income_pp + r3insurance_2cat + r3smoken + r3cesd + r3flone +
       r3hibpe + r3diabe + r3cancre + r3stroke + r3hearte + r3lunge + r3adl_2cat + r3iadl_2cat) * ns_1knot(fu_year),
  data = hrs_flusht_cohort1_long_cc,
  id = HHIDPN, waves = wave,
  weights = r3wtresp, family = gaussian)

broom::tidy(onetime_gcomp_gee) # %>% view()

# onetime_gcomp_lmer <- lme4::lmer(
#   memimp ~ as.factor(flusht_sts) + as.factor(flusht_sts) * ns_1knot(fu_year) + 
#     # first_memimp +
#     r3age_v2 + female + raceeth + raedyrs +
#     r3mstat_2cat + r3income_pp + r3insurance_2cat + r3smoken + r3cesd + r3flone +
#     r3hibpe + r3diabe + r3cancre + r3stroke + r3hearte + r3lunge + r3adl_2cat + r3iadl_2cat +
#     (1 + fu_year|HHIDPN),
#   data = hrs_flusht_cohort1_long_cc)


hrs_flusht_cohort1_long_cc %<>%
  mutate(memimp_pred = predict(onetime_gcomp_gee, newdata = hrs_flusht_cohort1_long_cc, type = "response"))

pred_mat_mean <- hrs_flusht_cohort1_long_cc %>%
  select(flusht_sts, fu_year, memimp_pred) %>%
  group_by(flusht_sts, fu_year) %>%
  summarise_at(vars(memimp_pred), 
               list(n = length, mean = mean, sd = sd)) %>%
  mutate(se = sd/sqrt(n), 
         flusht_sts = as.factor(flusht_sts),
         lower = mean - 1.96*se,
         upper = mean + 1.96*se)

pred_mat_mean %>%
  ggplot(aes(x = fu_year, y = mean, color = flusht_sts, group = flusht_sts)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8) +
  geom_line() +
  theme_bw()



#---- unadjusted plot cohort 1----
unadj_pred_mat_mean <- hrs_flusht_cohort1_long_cc %>%
  select(flusht_sts, fu_year, memimp) %>%
  group_by(flusht_sts, fu_year) %>%
  summarise_at(vars(memimp), 
               list(n = length, mean = mean, sd = sd)) %>%
  mutate(se = sd/sqrt(n), 
         flusht_sts = as.factor(flusht_sts),
         lower = mean - 1.96*se,
         upper = mean + 1.96*se)

unadj_pred_mat_mean %>%
  ggplot(aes(x = fu_year, y = mean, color = flusht_sts, group = flusht_sts)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8) +
  geom_line() +
  theme_bw()

#---- Check drop outs----
with(hrs_flusht_cohort1, table(r8dead, flusht_sts, useNA = "ifany"))
with(hrs_flusht_cohort1, table(r10dead, flusht_sts, useNA = "ifany"))
with(hrs_flusht_cohort1, table(r11dead, flusht_sts, useNA = "ifany"))


with(hrs_flusht_cohort1 %>% filter(r8dead == 0), table(r8inw, flusht_sts, useNA = "ifany"))
with(hrs_flusht_cohort1 %>% filter(r9dead == 0), table(r9inw, flusht_sts, useNA = "ifany"))


#---- Overall flu shot exposure ----
hrs_selected %<>%
  mutate(num_flusht = case_when(rowSums(!is.na(across(paste0("r", 3:13, "flusht")))) == 0 ~ NA_real_,
                                TRUE ~ rowSums(across(paste0("r", 3:13, "flusht")), na.rm = TRUE)),
         overall_flusht_3cat = case_when(num_flusht == 0 ~ 0,
                                         num_flusht %in% 1:3 ~ 1,
                                         num_flusht %in% 4:6 ~ 2))

# Sanity check
with(hrs_selected, table(num_flusht, useNA = "ifany"))
with(hrs_selected, table(num_flusht, overall_flusht_3cat, useNA = "ifany"))

hrs_selected_long <- hrs_selected %>%
  pivot_longer(cols = contains(c("memimp", "dead", "inw")),
               names_to = c("wave", ".value"), 
               names_pattern = "r(\\d+)(.*)") %>%
  mutate(wave = as.numeric(wave),
         fu_year = (wave - 3)*2)

hrs_selected_long_cc <- hrs_selected_long %>%
  select(HHIDPN, wave, inw, dead, fu_year, num_flusht, overall_flusht_3cat,
         memimp, r3wtresp, r3age_v2, female,
         raceeth, raedyrs, r3mstat_2cat, r3income_pp, r3insurance_2cat, r3smoken,
         r3cesd, r3flone, r3hibpe, r3diabe, r3cancre, r3stroke, r3hearte, r3lunge,
         r3adl_2cat, r3iadl_2cat) %>%
  na.omit()

onetime_gcomp_gee <- geepack::geeglm(
  memimp ~ as.factor(overall_flusht_3cat) * ns_1knot(fu_year) + 
    # first_memimp +
    (r3age_v2 + female + raceeth + raedyrs +
       r3mstat_2cat + r3income_pp + r3insurance_2cat + r3smoken + r3cesd + r3flone +
       r3hibpe + r3diabe + r3cancre + r3stroke + r3hearte + r3lunge + r3adl_2cat + r3iadl_2cat) * ns_1knot(fu_year),
  data = hrs_selected_long_cc,
  id = HHIDPN, waves = wave,
  weights = r3wtresp, family = gaussian)

broom::tidy(onetime_gcomp_gee) %>% view()

hrs_selected_long_cc %<>%
  mutate(memimp_pred = predict(onetime_gcomp_gee, newdata = hrs_selected_long_cc, type = "response"))

pred_mat_mean <- hrs_selected_long_cc %>%
  select(overall_flusht_3cat, fu_year, memimp_pred) %>%
  group_by(overall_flusht_3cat, fu_year) %>%
  summarise_at(vars(memimp_pred), 
               list(n = length, mean = mean, sd = sd)) %>%
  mutate(se = sd/sqrt(n), 
         overall_flusht_3cat = as.factor(overall_flusht_3cat),
         lower = mean - 1.96*se,
         upper = mean + 1.96*se)

pred_mat_mean %>%
  ggplot(aes(x = fu_year, y = mean, color = overall_flusht_3cat, group = overall_flusht_3cat)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8) +
  geom_line() +
  theme_bw()

#---- unadjusted plot overall ----
# CC
unadj_overall_mat_mean_CC <- hrs_selected_long_cc %>%
  select(num_flusht, wave, fu_year, memimp) %>%
  na.omit() %>%
  group_by(num_flusht, wave, fu_year) %>%
  summarise_at(vars(memimp), 
               list(n = length, mean = mean, sd = sd)) %>%
  mutate(se = sd/sqrt(n), 
         num_flusht = as.factor(num_flusht),
         lower = mean - 1.96*se,
         upper = mean + 1.96*se)

unadj_overall_mat_mean_CC %>%
  filter(num_flusht != 7) %>%
  ggplot(aes(x = wave, y = mean, color = num_flusht, group = num_flusht)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = 3:15) +
  labs(title = "CC model participants (wave 3 - 15) memory scores (observed), n = 13937")

with(hrs_selected_long_cc, length(unique(HHIDPN)))

unadj_overall_mat_mean_balanced_3_15 <- hrs_selected_long %>%
  filter(outcome_balanced_3_15 == 1) %>%
  select(num_flusht, wave, fu_year, memimp) %>%
  na.omit() %>%
  group_by(num_flusht, wave, fu_year) %>%
  summarise_at(vars(memimp), 
               list(n = length, mean = mean, sd = sd)) %>%
  mutate(se = sd/sqrt(n), 
         num_flusht = as.factor(num_flusht),
         lower = mean - 1.96*se,
         upper = mean + 1.96*se)

unadj_overall_mat_mean_balanced_3_15 %>%
  filter(num_flusht != 7) %>%
  ggplot(aes(x = wave, y = mean, color = num_flusht, group = num_flusht)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = 3:15) +
  labs(title = "Balanced (wave 3 - 15) memory scores (observed), n = 1946")

unadj_overall_mat_mean_balanced_3_15_flu_cat <- hrs_selected_long %>%
  filter(outcome_balanced_3_15 == 1) %>%
  select(overall_flusht_3cat, wave, fu_year, memimp) %>%
  na.omit() %>%
  group_by(overall_flusht_3cat, wave, fu_year) %>%
  summarise_at(vars(memimp), 
               list(n = length, mean = mean, sd = sd)) %>%
  mutate(se = sd/sqrt(n), 
         overall_flusht_3cat = as.factor(overall_flusht_3cat),
         lower = mean - 1.96*se,
         upper = mean + 1.96*se)

unadj_overall_mat_mean_balanced_3_15_flu_cat %>%
  ggplot(aes(x = wave, y = mean, color = overall_flusht_3cat, group = overall_flusht_3cat)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.6) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = 3:15) +
  labs(title = "Balanced (wave 3 - 15) memory scores (observed), n = 1946")

with(hrs_selected_long %>%
       filter(outcome_balanced_3_15 == 1), length(unique(HHIDPN))) # n = 1946

# #---- OLD ----
# #---- Binary exposure ----
# hrs_flusht_cohort1_long_cc %<>%
#   mutate(flusht_sts_2cat = case_when(flusht_sts >= 1 ~ 1, TRUE ~ 0))
# onetime_gcomp_gee <- geepack::geeglm(
#   memimp ~ as.factor(flusht_sts_2cat) + as.factor(flusht_sts_2cat) * ns_1knot(fu_year) + 
#     # first_memimp +
#     r3age_v2 + female + raceeth + raedyrs +
#     r3mstat_2cat + r3income_pp + r3insurance_2cat + r3smoken + r3cesd + r3flone +
#     r3hibpe + r3diabe + r3cancre + r3stroke + r3hearte + r3lunge + r3adl_2cat + r3iadl_2cat,
#   data = hrs_flusht_cohort1_long_cc,
#   id = HHIDPN, waves = wave,
#   weights = r3wtresp, family = gaussian)
# 
# broom::tidy(onetime_gcomp_gee)
# 
# 
# hrs_flusht_cohort1_long_cc %<>%
#   mutate(memimp_pred = predict(onetime_gcomp_gee, 
#                                newdata = hrs_flusht_cohort1_long_cc, type = "response"))
# 
# pred_mat_mean <- hrs_flusht_cohort1_long_cc %>%
#   select(flusht_sts_2cat, fu_year, memimp_pred) %>%
#   group_by(flusht_sts_2cat, fu_year) %>%
#   summarise_at(vars(memimp_pred), 
#                list(n = length, mean = mean, sd = sd)) %>%
#   mutate(se = sd/sqrt(n), 
#          flusht_sts_2cat = as.factor(flusht_sts_2cat),
#          lower = mean - 1.96*se,
#          upper = mean + 1.96*se)
# 
# pred_mat_mean %>%
#   ggplot(aes(x = fu_year, y = mean, color = flusht_sts_2cat, group = flusht_sts_2cat)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8) +
#   geom_line() +
#   theme_bw()
