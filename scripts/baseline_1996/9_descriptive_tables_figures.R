# Tables
# Created by: Yingyan Wu
# Sep.2.2024

#---- Load packages ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "haven", "labelled", 
       "gtsummary", "survey", "ggalluvial")
options(scipen = 999)

#---- Load the data ----
hrs_imp_r3resp <- readRDS(here::here("data", "analysis_data", "imp_pmm_stacked_ahead_hrs_outcome.RDS"))
n_imp = 10

#---- Preparing the data ----
# Sanity check
with(hrs_imp_r3resp, table(wave, imp, useNA = "ifany"))

bl_vars <- c("r3age_v2", "female", "raceeth", "raedyrs", "r3wtresp")
tv_vars <- c("age_v2", "mstat_2cat", "income_pp", 
             "insurance_2cat", "smoken", "cesd", "flone",
             # Health conditions
             "hibpe", "diabe", "cancre", "stroke", "hearte", "lunge",
             # ADL and IADLs
             "adl_2cat", "iadl_2cat",
             # exposure
             "flusht", 
             # Follow up status #"inw", "dead",
             "iwstat", "inw_imp",
             # Outcome (at previous waves)
             "memimp")

hrs_imp_r3resp_wide <- hrs_imp_r3resp %>%
  mutate(fu_time = case_when(iwstat %in% c(1, 4) ~ age_v2 - r3age_v2,
                             iwstat %in% c(5, 6, 7, 9) ~ NA_real_),
         # for those who died this wave, set the flusht value to be 0 for this wave
         # These people were not used anyway as the outcome at this wave is still missing for these people
         # this is just for the lmtp_tmle function to be able to run
         flusht = case_when(iwstat %in% c(5, 6, 7, 9) ~ 0,
                            TRUE ~ flusht)) %>%
  select(HHIDPN, wave, imp, all_of(c(bl_vars, tv_vars)), fu_time, -r3age_v2) %>%
  pivot_wider(names_from = "wave", 
              values_from = c(all_of(tv_vars), fu_time),
              names_glue = "r{wave}{.value}") %>%
  select(-paste0("r", seq(2, 14, by = 2), "flusht"), -contains("r2memimp")) %>%
  mutate(max_fu_time = 
           exec(pmax, !!!syms(paste0("r", 3:15, "fu_time")), na.rm = TRUE),
         num_flusht = rowSums(across(paste0("r", seq(3, 15, by = 2), "flusht"))),
         flusht_taker = case_when(num_flusht == 0 ~ 0,
                                  num_flusht == 7 ~ 2,
                                  num_flusht %in% 1:6 ~ 1),
         flusht_cat = case_when(num_flusht == 7 ~ 1,
                                TRUE ~ 0))

# Sanity check f/u time
with(hrs_imp_r3resp_wide, summary(max_fu_time))
with(hrs_imp_r3resp_wide, summary(num_flusht))
with(hrs_imp_r3resp_wide, table(num_flusht, flusht_taker, useNA = "ifany"))/10
hrs_imp_r3resp_wide %>% group_by(flusht_cat) %>% dplyr::summarize(
  n = n(),
  mean_age = mean(r3age_v2, na.rm = T))
with(hrs_imp_r3resp %>%
       filter(imp == 1, wave == 3),
     summary(r3age_v2))

#---- Baseline age per cohort ----
load(here::here("data", "analysis_data", "hrs_selected.RData"))
hrs_imp_r3resp_wide_imp1 <- hrs_imp_r3resp_wide %>% 
  filter(imp == 1) %>%
  left_join(hrs_selected %>% select(HHIDPN, hacohort), by = "HHIDPN")


hrs_imp_r3resp_wide_imp1 %>%
  mutate(hacohort = case_when(hacohort == 0 ~ "HRS/AHEAD overlap",
                              hacohort == 1 ~ "AHEAD",
                              hacohort == 3 ~ "HRS") %>%
           factor(., levels = c("HRS", "HRS/AHEAD overlap", "AHEAD"))) %>%
  ggplot() +
  geom_histogram(aes(x = r3age_v2, fill = hacohort), position = "identity", 
                 alpha = 0.5, 
                 binwidth = 1) +
  labs(x = "Baseline age in years", y = "Count", fill = "Sample cohort") +
  theme_bw()

ggsave(here::here("output", "figures", "baseline_age_per_cohort.png"),
       width = 8, height = 5, dpi = 300)



#---- Table 1. imputed description ----
##---- Labelled data for Table 1 ----
labelled_data <- hrs_imp_r3resp_wide %>%
  select(HHIDPN, imp, "r3age_v2", "female", "raceeth", "raedyrs", "r3wtresp",
         "flusht_cat", "num_flusht", "flusht_taker",
         starts_with("r3"),  
         -contains(c("iwstat", "fu_time", "inw_imp")), "max_fu_time") %>%
  dplyr::mutate(across(c("female", "r3insurance_2cat", "r3smoken",
                         "r3flone", "r3hibpe", "r3diabe", "r3cancre",
                         "r3stroke", "r3hearte", "r3lunge", "r3adl_2cat", 
                         "r3iadl_2cat"), 
                       function(x) factor(x, levels = c(1, 0),
                                          labels = c("Yes", "No")))) %>%
  dplyr::mutate(
    flusht_cat = factor(flusht_cat, levels = c(1, 0),
                        labels = c("Always had flu shot", "Other")),
    raceeth = factor(raceeth, levels = c(1, 2, 3, 4),
                     labels = c("Non-Hispanic White", "Non-Hispanic Black", 
                                "Hispanic", "Other")),
    r3mstat_2cat = factor(r3mstat_2cat, levels = c(1, 0),
                          labels = c("Married/partnered", 
                                     "Separated/divorced/widowed/not married")),
    flusht_taker = factor(flusht_taker, levels = c(2, 1, 0),
                          labels = c("Always taker", "Some flu shot", "No flu shot"))) %>%
  labelled::set_variable_labels(
    r3age_v2 = "Baseline age in years",
    female = "Women (%)",
    raceeth = "Race and ethnicity (%)",
    raedyrs = "Years of education (%)",
    r3mstat_2cat = "Marital status (%)",
    r3income_pp = "Household income per person (%)",
    r3insurance_2cat = "Having health insurance (%)",
    r3smoken = "Current smoker (%)",
    r3cesd = "Center for Epidemiological Studies-Depression score",
    r3flone = "Feeling lonely (%)",
    r3hibpe = "Self-reported hypertension (%)",
    r3diabe = "Self-reported diabetes (%)",
    r3cancre = "Self-reported cancer (%)",
    r3stroke = "Self-reported stroke (%)",
    r3hearte = "Self-reported heart disease (%)",
    r3lunge = "Self-reported lung disease (%)",
    r3adl_2cat = "Activities of daily living: having difficulty (%)",
    r3iadl_2cat = "Instrumental activities of daily living: having difficulty (%)",
    num_flusht = "Number of Flu shot at eligibility time reported during period",
    max_fu_time = "Follow-up time in years")

##---- Stratified table 1 ----
###---- categorical ----
cat_imp_svy <- labelled_data %>%
  select(female, raceeth, r3mstat_2cat, r3insurance_2cat, r3smoken, 
         r3flone, r3hibpe, r3diabe, r3cancre, r3stroke, r3hearte, r3lunge,
         r3adl_2cat, r3iadl_2cat, flusht_cat, flusht_taker,
         # baseline weights
         r3wtresp) %>%
  survey::svydesign(~ 1, data = ., weights = ~r3wtresp)

temp_cat_weighted <- tbl_svysummary(
  data = cat_imp_svy,
  statistic = list(all_categorical() ~ c("{n_unweighted} {p}")), # For latter splitting
  digits = all_categorical() ~ 1,
  by = "flusht_cat") %>%
  add_overall() %>%
  bold_labels()

temp_cat_weighted %>% as_flex_table()

unwtd_n <- temp_cat_weighted$table_styling$header$modify_stat_n_unweighted[6:8]
wtd_n <- temp_cat_weighted$table_styling$header$modify_stat_n[6:8]

table1_headers <- c(
  "Characteristic", 
  paste0("Overall \nn=", prettyNum(round(unwtd_n[1]/n_imp), big.mark = ","),
         "\nWeighted n=",  prettyNum(round(wtd_n[1]/n_imp),
                                     big.mark = ",")),
  paste0("Always had flu shot \nn=", prettyNum(round(unwtd_n[2]/n_imp), big.mark = ","),
         "\nWeighted n=", prettyNum(round(wtd_n[2]/n_imp), big.mark = ",")),
  paste0("Other \nn=", prettyNum(round(unwtd_n[3]/n_imp), big.mark = ","),
         "\nWeighted n=", prettyNum(round(wtd_n[3]/n_imp), big.mark = ",")))

table1_cat_weighted <- temp_cat_weighted %>% as_tibble() %>%
  set_colnames(c("Characteristic", "overall", "flushtcat", "no_flushtcat")) %>%
  filter(row_number() <= n() - 1) %>% # remove the last row in the table of median IQR for weights.
  separate(overall, c("overall_count", "overall_prop"), sep = " ") %>%
  separate(flushtcat, c("flushtcat_count", "flushtcat_prop"), sep = " ") %>%
  separate(no_flushtcat, c("no_flushtcat_count", "no_flushtcat_prop"), sep = " ") %>%
  dplyr::mutate(across(contains("count"),
                       function(x) prettyNum(round(as.numeric(gsub(",", "", x))/n_imp), 
                                             big.mark = ","))) %>%
  mutate(
    `Overall` = paste0(overall_count, " (", overall_prop, ")"),
    `Always had flu shot` = paste0(flushtcat_count, " (", flushtcat_prop, ")"),
    `Other` = paste0(no_flushtcat_count, " (", no_flushtcat_prop, ")")) %>%
  dplyr::select(-contains("_")) %>%
  set_colnames(table1_headers)
table1_cat_weighted[table1_cat_weighted == "NA (NA)"] <- NA

###---- continuous ----
# Every continuous variable was imputed.
cont_vars <- c("r3age_v2", "raedyrs", "r3income_pp", "r3cesd", "num_flusht",
               "max_fu_time")

hrs_imp_svy <- labelled_data %>%
  survey::svydesign(~ 1, data = ., weights = ~r3wtresp)
hrs_imp_svy_flushtcat <- subset(hrs_imp_svy, 
                                flusht_cat == "Always had flu shot")
hrs_imp_svy_no_flushtcat <-subset(hrs_imp_svy, 
                                  flusht_cat == "Other")
temp_cont_mean <- svymean(reformulate(cont_vars), hrs_imp_svy, na.rm = T) %>% 
  as.data.frame() %>% select(-SE) %>% set_colnames("overall_mean") %>%
  cbind(svymean(reformulate(cont_vars), hrs_imp_svy_flushtcat, na.rm = T) %>% 
          as.data.frame() %>% select(-SE) %>% set_colnames("flusht_mean"),
        svymean(reformulate(cont_vars), hrs_imp_svy_no_flushtcat, na.rm = T) %>% 
          as.data.frame() %>% select(-SE) %>% set_colnames("no_flusht_mean"))

temp_cont_tib <- data.frame()
for (i in 1:n_imp){
  temp_cont_svy <- labelled_data %>%
    filter(imp == i) %>%
    select(all_of(cont_vars), flusht_cat, r3wtresp) %>%
    survey::svydesign(~1, data = ., weights = ~r3wtresp)
  temp_cont_svy_flushtcat <- subset(temp_cont_svy, 
                                    flusht_cat == "Always had flu shot")
  temp_cont_svy_no_flushtcat <-subset(temp_cont_svy, 
                                      flusht_cat == "Other")
  temp_svyvar <- svyvar(reformulate(cont_vars), temp_cont_svy, na.rm = T) %>% as.data.frame() %>%
    select(-SE) %>% unique() %>% as.matrix() %>% diag() %>% t() %>% as_tibble() %>%
    set_colnames(cont_vars) %>% mutate(flusht = "Overall") %>%
    rbind(svyvar(reformulate(cont_vars), temp_cont_svy_flushtcat, na.rm = T) %>% 
            as.data.frame() %>%
            select(-SE) %>% unique() %>% as.matrix() %>% diag() %>% t() %>% as_tibble() %>%
            set_colnames(cont_vars) %>% mutate(flusht = "Always had flu shot"),
          svyvar(reformulate(cont_vars), temp_cont_svy_no_flushtcat, na.rm = T) %>% 
            as.data.frame() %>%
            select(-SE) %>% unique() %>% as.matrix() %>% diag() %>% t() %>% as_tibble() %>%
            set_colnames(cont_vars) %>% mutate(flusht = "Other"))
  
  temp_cont_tib %<>% rbind(temp_svyvar %>%
                             mutate(imp = i))
}

temp_cont_var <- temp_cont_tib %>%
  mutate(flusht = factor(flusht, 
                         levels = c("Overall", "Always had flu shot", "Other"))) %>%
  group_by(flusht) %>%
  dplyr::summarize_all(mean) %>% as.data.frame()

table1_cont_weighted <- temp_cont_mean %>%
  cbind(t(temp_cont_var %>% select(-imp, -flusht )) %>%
          set_colnames(paste0(c("overall", "flusht", "no_flusht"), "_var"))) %>%
  mutate(Characteristic = var_label(labelled_data[, cont_vars]) %>% unlist(), 
         across(contains("var"), ~ sqrt(.x), 
                .names = "{str_replace(.col, '_var', '_sd')}"),
         overall_mean_sd = paste0(sprintf(overall_mean, fmt = '%#.1f'), " (", 
                                  sprintf(overall_sd, fmt = '%#.1f'), ")"),
         flushtcat_mean_sd = paste0(sprintf(flusht_mean, fmt = '%#.1f'), " (", 
                                    sprintf(flusht_sd, fmt = '%#.1f'), ")"),
         no_flushtcat_mean_sd = paste0(sprintf(no_flusht_mean, fmt = '%#.1f'), " (", 
                                       sprintf(no_flusht_sd, fmt = '%#.1f'), ")")) %>%
  select(Characteristic, contains("mean_sd")) %>%
  as_tibble() %>%
  set_colnames(table1_headers)

##---- Full table 1 ----
table1_weighted <- rbind(table1_cont_weighted[1, ], # Baseline age
                         table1_cat_weighted[1:6, ], # female, raceeth
                         table1_cont_weighted[2, ], # years of edu
                         table1_cat_weighted[7:9, ], # marital status
                         table1_cont_weighted[3, ], # income
                         table1_cat_weighted[10:20, ],# health related
                         table1_cont_weighted[c(4, 6), ] # depression & follow up
) 
writexl::write_xlsx(table1_weighted, 
                    here::here("output", "tables", "table1_imp_weighted_stratified.xlsx"))

##---- Flu shot descriptive table ----
# Always taker, some flu shot, no flu shot
table_flusht_weighted <- rbind(table1_cat_weighted[21:24, ], # flu shot
                               table1_cont_weighted[5, ]) # flu shot numbers
writexl::write_xlsx(table_flusht_weighted, 
                    here::here("output", "tables", "table_flusht.xlsx"))

#---- Flu shot per wave ----
# Flu shot acceptance prop at each eligible time (give a sense of natural course?)
# All intervened proportion at each eligible time

# Current flu shot proportion per wave
flusht_per_wave <- hrs_imp_r3resp %>%
  filter(wave %in% seq(3, 15, by = 2)) %>%
  group_by(wave) %>%
  dplyr::summarize(
    n_sample_size = sum(inw_imp)/10,
    n_flusht = round(sum(flusht, na.rm = T)/10),
    prop = mean(flusht, na.rm = T) * 100,
    wtd.prop = Hmisc::wtd.mean(flusht, r3wtresp, na.rm = T) * 100) %>%
  mutate(`Flu shot reported at that wave` = 
           paste0(n_flusht, " (", sprintf(wtd.prop, fmt = '%#.1f'), ")"))

# # Sanity check
# with(hrs_imp_r3resp, table(wave, flusht, useNA = "ifany"))/10
# with(hrs_imp_r3resp, table(wave, inw_imp, useNA = "ifany"))/10

# Sustained flu shot per wave
flusht_sus_long <- hrs_imp_r3resp %>%
  filter(wave %in% seq(3, 15, by = 2), inw_imp == 1) %>%
  arrange(imp, HHIDPN, wave) %>%
  group_by(imp, HHIDPN) %>%
  mutate(flu_cumsum = cumsum(ifelse(flusht == 1, 1, 0)),
         flu_count = row_number(), 
         sustained = ifelse(flu_cumsum == flu_count, 1, 0)) %>%
  ungroup()

with(flusht_sus_long, table(wave, flu_cumsum, useNA = "ifany"))


flusht_sus_per_wave <- flusht_sus_long %>%
  filter(!is.na(flusht)) %>%
  group_by(wave) %>%
  summarise(n_sus_flusht = round(sum(sustained, na.rm = T)/10),
            n_sample_size = sum(inw_imp)/10,
            wtd.prop = Hmisc::wtd.mean(sustained, r3wtresp, na.rm = T) * 100,
            prop = mean(sustained, na.rm = T) * 100) %>%
  mutate(`Sustained flu shot reported at that wave` = 
           paste0(n_sus_flusht, " (", sprintf(wtd.prop, fmt = '%#.1f'), ")"))

# # Sanity check
# hrs_imp_r3resp_wide %>%
#   filter(r3flusht == 1, r5flusht == 1) %>%
#   with(., table(r5inw_imp, useNA = "ifany"))/10
# hrs_imp_r3resp_wide %>%
#   filter(r3flusht == 1, r5flusht == 1, r7flusht == 1, r9flusht == 1, 
#          r11flusht == 1, r13flusht == 1, r15flusht == 1) %>%
#   with(., table(r15inw_imp, useNA = "ifany"))/10


flusht_perwave_sus_tbl <- flusht_per_wave %>% 
  select(wave, n_sample_size, starts_with("Flu shot")) %>%
  left_join(flusht_sus_per_wave %>% select(-n_sus_flusht, -prop, -wtd.prop), 
            by = join_by(wave, n_sample_size))

writexl::write_xlsx(flusht_perwave_sus_tbl,
                    here::here("output", "tables", "table_flusht_per_wave.xlsx"))



#---- Fig 1. flu shot descriptive ----
#---- Trajectory plot ----
with(hrs_imp_r3resp %>% filter(!wave %in% seq(2, 14, by = 2)), 
     table(flusht, iwstat, useNA = "ifany"))

flusht_fortraj <- hrs_imp_r3resp %>%
  select(HHIDPN, imp, wave, flusht) %>%
  pivot_wider(names_from = "wave", 
              values_from = "flusht",
              names_glue = "r{wave}flusht") %>%
  select(paste0("r", seq(3, 15, by = 2), "flusht")) %>%
  group_by(across(all_of(paste0("r", seq(3, 15, by = 2), "flusht")))) %>%
  summarize(count = round(n()/n_imp)) %>% ungroup() %>%
  mutate(across(paste0("r", seq(3, 15, by = 2), "flusht"), 
                ~case_when(. == 0 ~ "0",
                           . == 1 ~ "1",
                           is.na(.) ~ "NA"))) %>%
  mutate(across(paste0("r", seq(3, 15, by = 2), "flusht"), 
                ~factor(., levels = c("0", "1", "NA"), 
                        labels = c("No flu shot", "Flu shot", "Drop out/death"))))

flusht_fortraj %>%
  ggplot(aes(y = count,
             axis1 = r3flusht,
             axis2 = r5flusht,
             axis3 = r7flusht,
             axis4 = r9flusht,
             axis5 = r11flusht,
             axis6 = r13flusht,
             axis7 = r15flusht)) +
  geom_alluvium(aes(fill = r3flusht)) +
  geom_stratum(alpha = 0) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(name = "Flu shot acceptance", values = c("#94CEF0", "#3271AE")) +
  scale_x_continuous(name = "Survey years",
                     breaks = 1:7,
                     minor_breaks = seq(1, 7, by = 0.5),
                     labels = seq(1996, 2020, by = 4)) +
  scale_y_continuous(name = "Count",
                     breaks = seq(0, 20000, by = 2500)) +
  theme_bw() +
  theme(legend.position = "None",
        text = element_text(size = 12, color = "black"),
        strip.text.x = element_text(size = 12, color = "black"),
        strip.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"))
#   # SER version
#   theme(legend.position = "None",
#         text = element_text(size = 14, color = "black"),
#         strip.text.x = element_text(size = 14, color = "black"),
#         strip.text.y = element_text(size = 14, color = "black"),
#         axis.text.x = element_text(size = 14, color = "black"),
#         axis.text.y = element_text(size = 14, color = "black"),
#         legend.text = element_text(size = 14, color = "black"),
#         strip.text = element_text(size = 14, color = "black"))
# ggsave(here::here("output", "figures", "fig_flusht_traj_SER.png"),
#        device = "png", width = 9, height = 7, units = "in", dpi = 600)
ggsave(here::here("output", "figures", "fig_flusht_traj.png"), 
       device = "png", width = 9, height = 7, units = "in", dpi = 300)

#---- Table S1: missingness ----
# Missing summary
missing_summary <- readxl::read_xlsx(here::here("output", "tables", 
                                                "missing_summary.xlsx"))

#---- Table S: death/ loss to follow up over time ----
# Sanity check on iwstat not change for each imp
with(hrs_imp_r3resp, table(imp, iwstat, useNA = "ifany"))

# Death and loss to follow up over time
hrs_imp_r3resp %>%
  filter(imp == 1) %>%
  select(HHIDPN, wave, iwstat) %>%
  mutate(iwstat = case_when(iwstat %in% c(1, 4) ~ "In Sample",
                            iwstat %in% c(5, 6) ~ "Death",
                            iwstat == 7 ~ "Formally drop out from the sample")) %>%
  ggplot() +
  geom_bar(aes(x = wave, fill = iwstat), position = "fill") +
  scale_fill_manual(name = "Follow up status", 
                    values = c("#94CEF0", "#3271AE", "#B0B0B0")) +
  geom_text(stat = "count", aes(label = after_stat(count)), 
            position = position_fill(vjust = 0.5), size = 3)
  
  
pivot_wider(names_from = "wave", 
              values_from = "iwstat",
              names_glue = "r{wave}iwstat") %>%
  select(-HHIDPN, -r2iwstat) %>%
  tbl_summary() %>%
  bold_labels()



cat_imp_svy <- labelled_data %>%
  select(female, raceeth, r3mstat_2cat, r3insurance_2cat, r3smoken, 
         r3flone, r3hibpe, r3diabe, r3cancre, r3stroke, r3hearte, r3lunge,
         r3adl_2cat, r3iadl_2cat, flusht_cat, flusht_taker,
         # baseline weights
         r3wtresp) %>%
  survey::svydesign(~ 1, data = ., weights = ~r3wtresp)

temp_cat_weighted <- tbl_svysummary(
  data = cat_imp_svy,
  statistic = list(all_categorical() ~ c("{n_unweighted} {p}")), # For latter splitting
  digits = all_categorical() ~ 1,
  by = "flusht_cat") %>%
  add_overall() %>%
  bold_labels()

temp_cat_weighted %>% as_flex_table()
