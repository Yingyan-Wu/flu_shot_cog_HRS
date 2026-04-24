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
hrs_imp_r5resp <- readRDS(here::here("data", "analysis_data", 
                                     "imp_pmm_stacked_outcome_2000_2022_ot52l90.RDS"))
n_imp = 10

#---- Preparing the data ----
# Sanity check
with(hrs_imp_r5resp, table(wave, imp, useNA = "ifany"))

bl_vars <- c("r5age_v2", "female", "raceeth", "raedyrs", "r5wtresp")
tv_vars <- c("age_v2", 
             "income_pp", "mstat_2cat", 
             "insurance_2cat", "doct_vis_mt5",
             # Health behaviors
             "smoken", "mod_drinking", "vigorous_pa",
             "cesd", # "flone",
             # Health conditions
             "shlt_fair_poor",
             "hibpe", "diabe", "cancre", "stroke", "hearte", "lunge",
             # ADL and IADLs
             "adl_2cat", "iadl_2cat",
             # exposure
             "flusht", 
             # Follow up status #"inw", "dead",
             "iwstat", "inw_imp", "dead",
             # # Survey weights
             # "wtresp",
             # Outcome (at previous waves)
             "memimp_stdz")

hrs_imp_r5resp_wide <- hrs_imp_r5resp %>%
  mutate(fu_time = case_when(iwstat %in% c(1, 4) ~ age_v2 - r5age_v2,
                             iwstat %in% c(5, 6, 7, 9) ~ NA_real_),
         # for those who died this wave, set the flusht value to be 0 for this wave
         # These people were not used anyway as the outcome at this wave is still missing for these people
         # this is just for the lmtp_tmle function to be able to run
         flusht = case_when(iwstat %in% c(5, 6, 7, 9) ~ 0,
                            TRUE ~ flusht)) %>%
  select(HHIDPN, wave, imp, all_of(c(bl_vars, tv_vars)), fu_time, -r5age_v2) %>%
  pivot_wider(names_from = "wave", 
              values_from = c(all_of(tv_vars), fu_time),
              names_glue = "r{wave}{.value}") %>%
  select(-paste0("r", seq(4, 16, by = 2), "flusht"), -contains("r2memimp")) %>%
  mutate(max_fu_time = 
           exec(pmax, !!!syms(paste0("r", 5:16, "fu_time")), na.rm = TRUE),
         num_flusht = rowSums(across(paste0("r", seq(5, 16, by = 2), "flusht"))),
         flusht_taker = case_when(num_flusht == 0 ~ 0,
                                  num_flusht == 6 ~ 2,
                                  num_flusht %in% 1:5 ~ 1),
         flusht_cat = case_when(num_flusht == 6 ~ 1,
                                TRUE ~ 0))

# Sanity check f/u time
with(hrs_imp_r5resp_wide, summary(max_fu_time))
with(hrs_imp_r5resp_wide, summary(num_flusht))
with(hrs_imp_r5resp_wide, table(num_flusht, flusht_taker, useNA = "ifany"))/10
hrs_imp_r5resp_wide %>% group_by(flusht_cat) %>% dplyr::summarize(
  n = n(),
  mean_age = mean(r5age_v2, na.rm = T))
with(hrs_imp_r5resp %>%
       filter(imp == 1, wave == 5),
     summary(r5age_v2))
# death
with(hrs_imp_r5resp_wide %>%
       filter(imp == 1),
     table(r16dead, useNA = "ifany"))
hrs_imp_r5resp_wide %>% filter(r16dead == 1) %>% with(., sum(r5wtresp))/
  sum(hrs_imp_r5resp_wide$r5wtresp) # 0.5694435

# INW
hrs_imp_r5resp_wide %>% filter(r16inw_imp == 1) %>% with(., sum(r5wtresp))/
  sum(hrs_imp_r5resp_wide$r5wtresp) #  0.2869848
hrs_imp_r5resp_wide %>% filter(r15inw_imp == 1) %>% with(., sum(r5wtresp))/
  sum(hrs_imp_r5resp_wide$r5wtresp) # 0.3423993

# drop out
hrs_imp_r5resp_wide %>% filter(r16iwstat == 7) %>% with(., sum(r5wtresp))/
  sum(hrs_imp_r5resp_wide$r5wtresp) #  0.1435717

#---- Age distribution at baseline and end of follow-up ----
hrs_imp_r5resp_wide %>%
  filter(imp == 1) %>%
  select(paste0("r", c(5, 16), "age_v2"),
         paste0("r", c(5, 16), "iwstat"),
         r5wtresp) %>%
  pivot_longer(cols = c(-r5wtresp), 
               names_to = c("wave", ".value"), 
               names_pattern = "r(\\d+)(.*)") %>%
  filter(iwstat %in% c(1, 4)) %>%
  mutate(wave = factor(wave, levels = c(5, 16))) %>%
  ggplot() +
  # Histogram (show count in y axis)
  geom_histogram(aes(x = age_v2, fill = wave), 
                 position = "identity", alpha = 0.5, binwidth = 1) +
  labs(x = "Age in years", y = "Counts", fill = "wave") +
  scale_fill_manual(values = c("steelblue", "darkgrey"), 
                    labels = c("Baseline (year 2000)", 
                               "End of follow-up (year 2022)")) +
  scale_x_continuous(breaks = seq(50, 115, by = 5)) +
  scale_y_continuous(breaks = seq(0, 1000, by = 100)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave(here::here("output", "figures", "age_distribution_baseline_endfu.png"),
       width = 9, height = 7, dpi = 300)

# # ---- Baseline age per cohort ----
# load(here::here("data", "analysis_data", "hrs_selected_2000_2022_ot52l90.RData"))
# hrs_imp_r5resp_wide_imp1 <- hrs_imp_r5resp_wide %>%
#   filter(imp == 1) %>%
#   left_join(hrs_selected %>% select(HHIDPN, hacohort, racohbyr), by = "HHIDPN")
# 
# 
# hrs_imp_r5resp_wide_imp1 %>%
#   # mutate(hacohort = case_when(hacohort == 0 ~ "HRS/AHEAD overlap",
#   #                             hacohort == 1 ~ "AHEAD",
#   #                             hacohort == 3 ~ "HRS") %>%
#   #          factor(., levels = c("HRS", "HRS/AHEAD overlap", "AHEAD"))) %>%
#   ggplot() +
#   geom_density(aes(x = r5age_v2, color = as.factor(racohbyr)), position = "identity",
#                  alpha = 0.5,
#                  binwidth = 1) +
#   labs(x = "Baseline age in years", y = "Count", color = "Sample cohort") +
#   theme_bw()
# 
# # ggsave(here::here("output", "figures", "baseline_age_per_cohort.png"),
# #        width = 8, height = 5, dpi = 300)



#---- Table 1. imputed description ----
##---- Labelled data for Table 1 ----
labelled_data <- hrs_imp_r5resp_wide %>%
  select(HHIDPN, imp, "r5age_v2", "female", "raceeth", "raedyrs", "r5wtresp",
         "flusht_cat", "num_flusht", "flusht_taker",
         starts_with("r5"),  
         -contains(c("iwstat", "fu_time", "inw_imp")), "max_fu_time") %>%
  dplyr::mutate(across(c("female", "r5insurance_2cat", "r5smoken","r5mstat_2cat",
                         # "r5flone", 
                         "r5hibpe", "r5diabe", "r5cancre",
                         "r5stroke", "r5hearte", "r5lunge", "r5adl_2cat", 
                         "r5iadl_2cat"), 
                       function(x) factor(x, levels = c(1, 0),
                                          labels = c("Yes", "No")))) %>%
  dplyr::mutate(
    flusht_cat = factor(flusht_cat, levels = c(1, 0),
                        labels = c("Always had flu shot", "Other")),
    raceeth = factor(raceeth, levels = c(1, 2, 3, 4),
                     labels = c("Non-Hispanic White", "Non-Hispanic Black", 
                                "Hispanic", "Other")),
    # r5mstat_2cat = factor(r5mstat_2cat, levels = c(1, 0),
    #                       labels = c("Married/partnered", 
    #                                  "Separated/divorced/widowed/not married")),
    flusht_taker = factor(flusht_taker, levels = c(2, 1, 0),
                          labels = c("Always taker", "Some flu shot", "No flu shot"))) %>%
  labelled::set_variable_labels(
    r5age_v2 = "Baseline age in years",
    female = "Women",
    raceeth = "Race and ethnicity",
    raedyrs = "Years of education",
    r5mstat_2cat = "Married/partnered",
    r5income_pp = "Household income per person",
    r5insurance_2cat = "Having health insurance",
    r5doct_vis_mt5 = "≥ 5 doctor visits in the past 2 years",
    r5mod_drinking = "Moderate drinking", 
    r5vigorous_pa = "≥ 3 times of vigorous physical activities per week",
    r5shlt_fair_poor = "Self-rated health: Fair/poor",
    r5smoken = "Current smoker",
    r5cesd = "Center for Epidemiological Studies-Depression score",
    # r5flone = "Feeling lonely",
    r5hibpe = "Self-reported hypertension",
    r5diabe = "Self-reported diabetes",
    r5cancre = "Self-reported cancer",
    r5stroke = "Self-reported stroke",
    r5hearte = "Self-reported heart disease",
    r5lunge = "Self-reported lung disease",
    r5adl_2cat = "Activities of daily living: having difficulty",
    r5iadl_2cat = "Instrumental activities of daily living: having difficulty",
    num_flusht = "Number of Flu shot at eligibility time reported during period",
    max_fu_time = "Follow-up time in years")

##---- Stratified table 1 ----
###---- categorical ----
cat_imp_svy <- labelled_data %>%
  select(female, raceeth, r5mstat_2cat, r5insurance_2cat, r5doct_vis_mt5,
         r5smoken, r5mod_drinking, r5vigorous_pa,
         # r5flone, 
         r5shlt_fair_poor,
         r5hibpe, r5diabe, r5cancre, r5stroke, r5hearte, r5lunge,
         r5adl_2cat, r5iadl_2cat, flusht_cat, flusht_taker,
         # baseline weights
         r5wtresp) %>%
  survey::svydesign(~ 1, data = ., weights = ~r5wtresp)

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
cont_vars <- c("r5age_v2", "raedyrs", "r5income_pp", "r5cesd", "num_flusht",
               "max_fu_time")

hrs_imp_svy <- labelled_data %>%
  survey::svydesign(~ 1, data = ., weights = ~r5wtresp)
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
    select(all_of(cont_vars), flusht_cat, r5wtresp) %>%
    survey::svydesign(~1, data = ., weights = ~r5wtresp)
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
                         table1_cat_weighted[7, ], # marital status
                         table1_cont_weighted[3, ], # income
                         table1_cat_weighted[8:21, ],# health related
                         table1_cont_weighted[c(4, 6), ] # depression & follow up
) 
writexl::write_xlsx(table1_weighted, 
                    here::here("output", "tables", "table1_imp_weighted_stratified_2000_2022_ot52l90.xlsx"))

##---- Flu shot descriptive table ----
# Always taker, some flu shot, no flu shot
table1_cont_weighted[5, ]
# table_flusht_weighted <- rbind(table1_cat_weighted[25:28, ], # flu shot
#                                table1_cont_weighted[5, ]) # flu shot numbers
# writexl::write_xlsx(table_flusht_weighted,
#                     here::here("output", "tables", "table_flusht.xlsx"))

#---- eTable. pre imputation table 1 ----
load(here::here("data", "analysis_data", "hrs_selected_2000_2022_ot52l90.RData"))

hrs_selected_fu_time <- hrs_selected %>%
  select(HHIDPN, contains(c("iwstat", "age_v2"))) %>%
  pivot_longer(cols = -HHIDPN,
               names_to = c("wave", ".value"),
               values_to = c("iwstat", "age_v2"),
               names_pattern = "r(\\d+)(.*)") %>%
  left_join(hrs_selected %>% select(HHIDPN, r5age_v2), by = "HHIDPN") %>%
  mutate(fu_time = case_when(iwstat %in% c(1, 4) ~ age_v2 - r5age_v2,
                             iwstat %in% c(5, 6, 7, 9) ~ NA_real_)) %>%
  select(-r5age_v2) %>%
  pivot_wider(names_from = "wave",
              values_from = c("iwstat", "age_v2", "fu_time"),
              names_glue = "r{wave}{.value}") %>%
  mutate(max_fu_time = 
           exec(pmax, !!!syms(paste0("r", 5:16, "fu_time")), na.rm = TRUE)) %>%
  select(HHIDPN, max_fu_time)
# Sanity check
with(hrs_selected_fu_time, summary(max_fu_time))
hrs_selected <- hrs_selected %>%
  left_join(hrs_selected_fu_time, by = "HHIDPN")
hrs_selected %>% with(table(is.na(r5wtresp), r5iwstat, useNA = "ifany"))

labelled_data_before_imp <- hrs_selected %>%
  filter(!is.na(r5wtresp)) %>%
  select(r5age_v2, female, raceeth, raedyrs, 
         r5mstat_2cat, r5income_pp, r5insurance_2cat, r5doct_vis_mt5,
         # Health behaviors
         r5smoken, r5mod_drinking, r5vigorous_pa,
         # r5flone,
         # Health conditions
         r5shlt_fair_poor,
         r5hibpe, r5diabe, r5cancre, r5stroke, r5hearte, r5lunge,
         # ADL and IADLs
         r5adl_2cat, r5iadl_2cat, r5cesd, max_fu_time, r5wtresp) %>%
  dplyr::mutate(across(c("female", "r5insurance_2cat", "r5smoken",
                         "r5mstat_2cat", # "r5flone", 
                         "r5hibpe", "r5diabe", "r5cancre",
                         "r5stroke", "r5hearte", "r5lunge", "r5adl_2cat", 
                         "r5iadl_2cat", "r5doct_vis_mt5", "r5vigorous_pa"), 
                       function(x) factor(x, levels = c(1, 0),
                                          labels = c("Yes", "No")))) %>%
  dplyr::mutate(
    raceeth = factor(raceeth, levels = c(1, 2, 3, 4),
                     labels = c("Non-Hispanic White", "Non-Hispanic Black", 
                                "Hispanic", "Other")),
    # r5mstat_2cat = factor(r5mstat_2cat, levels = c(1, 0),
    #                       labels = c("Married/partnered", 
    #                                  "Separated/divorced/widowed/not married")),
    raedyrs = as.numeric(raedyrs)) %>%
  labelled::set_variable_labels(
    r5age_v2 = "Baseline age in years",
    female = "Women",
    raceeth = "Race and ethnicity",
    raedyrs = "Years of education",
    r5mstat_2cat = "Married/partnered",
    r5income_pp = "Household income per person",
    r5insurance_2cat = "Having health insurance",
    r5doct_vis_mt5 = "≥ 5 doctor visits in the past 2 years",
    r5mod_drinking = "No to moderate alcohol use", 
    r5vigorous_pa = "≥ 3 times of vigorous physical activities per week",
    r5shlt_fair_poor = "Self-rated health: Fair/poor",
    r5smoken = "Current smoker",
    r5cesd = "Center for Epidemiological Studies-Depression score",
    # r5flone = "Feeling lonely",
    r5hibpe = "Self-reported hypertension",
    r5diabe = "Self-reported diabetes",
    r5cancre = "Self-reported cancer",
    r5stroke = "Self-reported stroke",
    r5hearte = "Self-reported heart disease",
    r5lunge = "Self-reported lung disease",
    r5adl_2cat = "Activities of daily living: having difficulty",
    r5iadl_2cat = "Instrumental activities of daily living: having difficulty",
    max_fu_time = "Follow-up time in years") %>%
  filter(!is.na(r5wtresp))

data_before_imp_svy <- labelled_data_before_imp %>%
  survey::svydesign(~ 1, data = ., weights = ~r5wtresp)

table1_before_imp_weighted <- tbl_svysummary(
  data = data_before_imp_svy,
  type = list(r5cesd ~ "continuous"),
  statistic = list(all_continuous() ~ "{mean} ({sd})",
                   all_categorical() ~ "{n_unweighted} ({p})"),
  missing = "no",
  digits = list(all_continuous() ~ 1,
                all_categorical() ~ c(0, 1))) %>%
  bold_labels()

table1_before_imp_weighted %>% as_flex_table()

table1_before_imp_weighted %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = here::here("output", "tables",
                                            "table1_before_imp_2000_2022.docx"))

#---- eTable. Flu shot per wave ----
# Flu shot acceptance prop at each eligible time (give a sense of natural course?)
# All intervened proportion at each eligible time

# Current flu shot proportion per wave
flusht_per_wave <- hrs_imp_r5resp %>%
  filter(wave %in% seq(5, 16, by = 2)) %>%
  group_by(wave) %>%
  dplyr::summarize(
    n_sample_size = sum(inw_imp)/10,
    n_flusht = round(sum(flusht, na.rm = T)/10),
    prop = mean(flusht, na.rm = T) * 100,
    wtd.prop = Hmisc::wtd.mean(flusht, r5wtresp, na.rm = T) * 100) %>%
  mutate(`Flu shot reported at that wave` = 
           paste0(n_flusht, " (", sprintf(wtd.prop, fmt = '%#.1f'), ")"))

# # Sanity check
# with(hrs_imp_r5resp, table(wave, flusht, useNA = "ifany"))/10
# with(hrs_imp_r5resp, table(wave, inw_imp, useNA = "ifany"))/10

# Sustained flu shot per wave
flusht_sus_long <- hrs_imp_r5resp %>%
  filter(wave %in% seq(5, 16, by = 2), inw_imp == 1) %>%
  arrange(imp, HHIDPN, wave) %>%
  group_by(imp, HHIDPN) %>%
  mutate(flu_cumsum = cumsum(ifelse(flusht == 1, 1, 0)),
         flu_count = row_number(), 
         sustained = ifelse(flu_cumsum == flu_count, 1, 0),
         never = ifelse(flu_cumsum == 0, 1, 0)) %>%
  ungroup()

with(flusht_sus_long, table(wave, flu_cumsum, useNA = "ifany"))
with(flusht_sus_long, table(wave, never, useNA = "ifany"))


flusht_sus_per_wave <- flusht_sus_long %>%
  filter(!is.na(flusht)) %>%
  group_by(wave) %>%
  summarise(n_sus_flusht = round(sum(sustained, na.rm = T)/10),
            n_never_flusht = round(sum(never, na.rm = T)/10),
            n_sample_size = sum(inw_imp)/10,
            sus_wtd_prop = Hmisc::wtd.mean(sustained, r5wtresp, na.rm = T) * 100,
            sus_prop = mean(sustained, na.rm = T) * 100,
            never_wtd_prop = Hmisc::wtd.mean(never, r5wtresp, na.rm = T) * 100,
            never_prop = mean(never, na.rm = T) * 100) %>%
  mutate(`Sustained flu shot reported at that wave` = 
           paste0(n_sus_flusht, " (", sprintf(sus_wtd_prop, fmt = '%#.1f'), ")"),
         `Never reported` = 
           paste0(n_never_flusht, " (", sprintf(never_wtd_prop, fmt = '%#.1f'), ")"))

# # Sanity check
# hrs_imp_r5resp_wide %>%
#   filter(r5flusht == 1, r5flusht == 1) %>%
#   with(., table(r5inw_imp, useNA = "ifany"))/10
# hrs_imp_r5resp_wide %>%
#   filter(r5flusht == 1, r5flusht == 1, r7flusht == 1, r9flusht == 1, 
#          r11flusht == 1, r13flusht == 1, r16flusht == 1) %>%
#   with(., table(r16inw_imp, useNA = "ifany"))/10


flusht_perwave_sus_tbl <- flusht_per_wave %>% 
  select(wave, n_sample_size, starts_with("Flu shot")) %>%
  left_join(flusht_sus_per_wave %>% select(-contains("_"), n_sample_size), 
            by = join_by(wave, n_sample_size))

writexl::write_xlsx(flusht_perwave_sus_tbl,
                    here::here("output", "tables", "table_flusht_per_wave_2000_2022_ot52l90.xlsx"))



#---- Fig 1. flu shot descriptive ----
with(hrs_imp_r5resp %>% filter(!wave %in% seq(2, 16, by = 2)), 
     table(flusht, iwstat, useNA = "ifany"))

flusht_fortraj <- hrs_imp_r5resp %>%
  select(HHIDPN, imp, wave, flusht) %>%
  pivot_wider(names_from = "wave", 
              values_from = "flusht",
              names_glue = "r{wave}flusht") %>%
  select(paste0("r", seq(5, 16, by = 2), "flusht")) %>%
  group_by(across(all_of(paste0("r", seq(5, 16, by = 2), "flusht")))) %>%
  dplyr::summarize(count = round(n()/n_imp)) %>% ungroup() %>%
  mutate(across(paste0("r", seq(5, 16, by = 2), "flusht"), 
                ~case_when(. == 0 ~ "0",
                           . == 1 ~ "1",
                           is.na(.) ~ "NA"))) %>%
  mutate(across(paste0("r", seq(5, 16, by = 2), "flusht"), 
                ~factor(., levels = c("0", "1", "NA"), 
                        labels = c("Not vaccinated", "Vaccinated", "Drop out/death"))))

flusht_fortraj %>%
  ggplot(aes(y = count,
             axis1 = r5flusht,
             axis3 = r7flusht,
             axis4 = r9flusht,
             axis5 = r11flusht,
             axis6 = r13flusht,
             axis7 = r15flusht)) +
  geom_alluvium(aes(fill = r5flusht)) +
  geom_stratum(alpha = 0) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(name = "Flu shot acceptance", values = c("#94CEF0", "#3271AE")) +
  scale_x_continuous(name = "Follow-up years",
                     breaks = 1:6,
                     minor_breaks = seq(1, 6, by = 0.5),
                     labels = seq(0, 20, by = 4)) +
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
ggsave(here::here("output", "figures", "fig_flusht_traj_2000_2022.png"), 
       device = "png", width = 9, height = 7, units = "in", dpi = 300)

#---- Table S1: missingness ----
# Missing summary
missing_summary <- readxl::read_xlsx(here::here("output", "tables", 
                                                "missing_summary_2000_2022.xlsx"))

missing_summary_fmt <- missing_summary %>%
  mutate(labels = case_when(feature == "iwstat" ~ "Interview status",
                            feature == "r5wtresp" ~ "Baseline survey weight",
                            feature == "female" ~ "Sex/gender",
                            feature == "age_v2" ~ "Age",
                            feature == "raceeth" ~ "Race and ethnicity",
                            feature == "raedyrs" ~ "Years of education",
                            feature == "mstat_2cat" ~ "Marital status",
                            feature == "income_pp" ~ "Household income per person",
                            feature == "insurance_2cat" ~ "Health insurance",
                            feature == "doct_vis_mt5" ~ "Number of doctor visits",
                            feature == "mod_drinking" ~ "Alcohol drinking",
                            feature == "vigorous_pa" ~ "Physical activity",
                            feature == "shlt_fair_poor" ~ "Self-rated health",
                            feature == "smoken" ~ "Smoking status",
                            feature == "cesd" ~ "Center for Epidemiological Studies-Depression score",
                            # feature == "flone" ~ "Loneliness status",
                            feature == "hibpe" ~ "Self-reported hypertension",
                            feature == "diabe" ~ "Self-reported diabetes",
                            feature == "cancre" ~ "Self-reported cancer",
                            feature == "stroke" ~ "Self-reported stroke",
                            feature == "hearte" ~ "Self-reported heart disease",
                            feature == "lunge" ~ "Self-reported lung disease",
                            feature == "adl_2cat" ~ "Activities of daily living",
                            feature == "flusht" ~ "Influenza vaccination",
                            feature == "iadl_2cat" ~ "Instrumental activities of daily living",
                            feature == "memimp" ~ "Memory score",
                            feature == "hhres" ~ "Household size",
                            feature == "itot" ~ "Household income"),
         pct_missing = round(pct_missing, 2)) %>%
  filter(!is.na(labels)) %>%
  select(labels, pct_missing)

writexl::write_xlsx(missing_summary_fmt,
                    here::here("output", "tables", "missing_summary_2000_2022_fmt.xlsx"))

#---- death/ loss to follow up over time ----
# Sanity check on iwstat not change for each imp
with(hrs_imp_r5resp, table(imp, iwstat, useNA = "ifany"))
with(hrs_imp_r5resp, table(wave, iwstat, useNA = "ifany"))

# Death and loss to follow up over time
hrs_death_ltfu_prop <- hrs_imp_r5resp %>%
  filter(wave >= 5) %>%
  select(HHIDPN, wave, iwstat) %>%
  mutate(iwstat = case_when(
    iwstat %in% c(1, 4) ~ "In sample",
    iwstat %in% c(5, 6) ~ "Death",
    iwstat == 7 ~ "Formally drop out from the sample")) %>%
  count(wave, iwstat) %>%
  group_by(wave) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(fu_yr = (wave - 5) * 2)

hrs_death_ltfu_prop %>%
  ggplot(aes(x = fu_yr, y = prop, fill = iwstat)) +
  geom_col(position = "fill") +
  geom_text(aes(label = scales::percent(prop)),
            position = position_fill(vjust = 0.5),
            size = 3) +
  scale_fill_manual(
    name = "Follow up status",
    values = c("Dark Grey", "light Grey", "light blue")) +
  labs(x = "Follow-up years", y = "Proportion %") +
  scale_x_continuous(breaks = seq(0, 22, by = 2)) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12, color = "black"),
        strip.text.x = element_text(size = 12, color = "black"),
        strip.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"))
ggsave(here::here("output", "figures", "fig_death_ltfu_2000_2022.png"), 
       device = "png", width = 9, height = 7, units = "in", dpi = 300)
