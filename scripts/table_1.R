# HRS_data_table_1_proposal
# Created by: Yingyan Wu
# Nov.05.2022

#---- Package loading + options ----
rm(list = ls())
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven", "labelled",
       "lubridate", "gtsummary")

# Directory to Box
# path_to_box <- "C:/Users/yingyan_wu/Box"
path_to_box <- "C:/Users/Yingyan Wu/Box"

#---- Load RAND data ----
load(here::here("data", "analysis_data", "rand_flusht.RData"))

#---- Table 1 ----
table_1 <- RAND %>% 
  filter(!is.na(n_flusht_wave)) %>% # at least report one wave for flusht
  filter(r3age_y_int >= 50) %>%
  select(r3age_y_int, female, raceeth, ed_cat, r3mstat,
         r3hibpe, r3diabe, r3hearte, r3stroke, r3cancre, r3lunge, 
         smoker, r3cesd, n_flusht_wave, n_flusht) %>%
  labelled::set_variable_labels(
    r3age_y_int = "Age in years at wave 3 (1996)",
    female = "Female (%)",
    raceeth = "Race/Ethnicity (%)",
    ed_cat = "Education level (%)",
    r3mstat = "Marital status (%)",
    r3hibpe = "Hypertension(Diagnosed) (%)",
    r3diabe = "Diabetes (%)",
    r3hearte = "Heart disease (%)",
    r3stroke = "Stroke (%)",
    r3cancre = "Cancer (%)",
    r3lunge = "Lung disease (%)",
    smoker = "Smoking status (%)",
    r3cesd = "Wave 3 (1996) CES-D",
    n_flusht_wave = "# of self-reported flushot status", 
    n_flusht = "# of self-reported flushot injections") %>%
  labelled::drop_unused_value_labels() %>%
  labelled::set_value_labels(female = c("Yes" = 1, "No" = 0),
                             ed_cat = c("Less than high school" = 1, "High school" = 2, 
                                        "Some college" = 3, "Bachelor's degree" = 4, 
                                        "Graduate studies" = 5),
                             r3mstat = c("Married/Partnered" = 1,
                                         "Not married/partnered" = 2,
                                         "Widowed" = 3),
                             r3hibpe = c("Yes" = 1, "No" = 0),
                             r3diabe = c("Yes" = 1, "No" = 0),
                             r3hearte = c("Yes" = 1, "No" = 0),
                             r3stroke = c("Yes" = 1, "No" = 0),
                             r3cancre = c("Yes" = 1, "No" = 0),
                             r3lunge = c("Yes" = 1, "No" = 0),
                             smoker = c("Ever smoked" = 1, "Never smoked" = 0)) %>%
  modify_if(is.labelled, to_factor) %>%
  modify_if(is.factor, 
            function(x) {forcats::fct_explicit_na(x, na_level = "Missing")}) %>%
  tbl_summary(
    statistic = list(all_categorical() ~ "{n} ({p}%)",
                     all_continuous() ~ c("{mean} ({sd})",
                                          "({min}, {max})")),
    type = list(c("r3cesd", "n_flusht_wave", "n_flusht") ~ "continuous2",
      all_continuous() ~ "continuous2"),
    digits = list(all_continuous() ~ 1,
                  all_categorical() ~ c(0, 1)),
    missing = "ifany",
    missing_text = "Missing") %>% 
  # modify_spanning_header(starts_with("stat_") ~ "**Baseline CES-D**") %>%
  bold_labels() %>%
  modify_footnote(
    update = list(
      starts_with("stat_") ~
        "Mean (SD) for continuous variables; n (%) for categorical variables")
  )

table_1 %>% as_flex_table() %>% flextable::save_as_docx(
  path = here::here("output", "tables", "table1_proposal.docx"))

#---- **save the table ----
# The table will be formatted in a linked excel
writexl::write_xlsx(as_tibble(table_1), 
                    here::here("output", "tables", "table1_proposal.xlsx"))
