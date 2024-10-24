# Flu -> Cog HRS
# 7.13.2022

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "plyr", "haven", "labelled",
       "lubridate")

#---- function ----
source(here::here("scripts", "functions.R"))
# Directory to Box
# path_to_box <- "C:/Users/yingyan_wu/Box"
# path_to_box <- "C:/Users/Yingyan Wu/Box"
path_to_box <- "/Users/yingyanwu/Library/CloudStorage/Box-Box/"

#---- Load RAND data ----
number_waves <- 1:14
rand_variables <- c("hhidpn", "ragender", "raracem", "rahispan", "rabmonth", 
                    "rabyear", "rabdate", "radmonth", "radyear", "raddate",
                    paste0("r", number_waves, "agem_e"), "raedyrs", 
                    paste0("r", number_waves, "iwstat"),
                    "raedegrm", 
                    paste0("r", 3:14, "flusht"),
                    paste0("r", number_waves, "mstat"),
                    # paste0("r", seq(8, 14, by = 1), "pmwght"), 
                    # paste0("r", number_waves, "weight"),
                    # paste0("r", seq(8, 14, by = 1), "pmhght"), 
                    # paste0("r", number_waves, "height"),
                    # paste0("r", number_waves, "bmi"), 
                    # paste0("r", seq(8, 14, by = 1), "pmbmi"),
                    # paste0("r", seq(8, 14, by = 1), "pmwaist"),
                    paste0("r", number_waves, "hibpe"),
                    paste0("r", number_waves, "diabe"),
                    paste0("r", number_waves, "cancre"),
                    paste0("r", number_waves, "stroke"),
                    paste0("r", number_waves, "hearte"),
                    paste0("r", number_waves, "lunge"),
                    paste0("r", seq(4 ,9, by = 1), "memrye"),
                    paste0("r", number_waves, "conde"),
                    paste0("r", number_waves, "smoken"),
                    paste0("r", seq(3, 14, by = 1), "drinkd"),
                    paste0("r", seq(3, 14, by = 1), "drinkn"),
                    paste0("r", seq(2, 14, by = 1), "cesd"),
                    paste0("r", number_waves, "shlt"))
RAND <- read_dta(paste0(path_to_box, "/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2018v1.dta"),
                 col_select = all_of(rand_variables)) 

colnames(RAND)[1] <- "HHIDPN" #For merging

#---- Age ----
age_m <- RAND %>% dplyr::select(contains("agem_e")) %>% 
  apply(., 1, impute_ages) %>% t() 

#Exact ages
RAND[, paste0("r", number_waves, "age_y")] <- age_m/12
#Ages rounded down to nearest year
RAND[, paste0("r", number_waves, "age_y_int")] <- floor(age_m/12)

# Drop RAND age variables
RAND %<>% dplyr::select(-paste0("r", number_waves, "agem_e"))

#---- gender ----
RAND %<>% 
  mutate("female" = ifelse(ragender == 2, 1, 0))

#Drop RAND gender variable 
RAND %<>% dplyr::select(-ragender)

#---- race-eth ----
#Code any hispanic as 1, else 0
RAND %<>% 
  mutate(hispanic = ifelse(rahispan == 0 | is.na(rahispan), 0, 1), 
         white = ifelse(raracem == 1 & hispanic == 0, 1, 0),
         black = ifelse(raracem == 2 & hispanic == 0, 1, 0),
         other = ifelse(raracem == 3 & hispanic == 0, 1, 0),
         unknown_race_eth = ifelse(is.na(raracem) & hispanic == 0, 1, 0)) %>%
  mutate(raceeth = case_when(hispanic == 0 & white == 1 ~ "Non-Hispanic White",
                             hispanic == 0 & black == 1 ~ "Non-Hispanic Black",
                             hispanic == 1 ~ "Hispanic",
                             other == 1 ~ "Other"))

# #Sanity check
# table(RAND$hispanic, RAND$rahispan, useNA = "ifany")
# table(RAND$hispanic, RAND$raracem, RAND$black, useNA = "ifany")
# table(RAND$hispanic, RAND$raracem, RAND$other, useNA = "ifany")
# table(RAND$hispanic, RAND$raracem, RAND$unknown_race_eth,
#       useNA = "ifany")
# table(RAND$unknown_race_eth, useNA = "ifany")
# colSums(is.na(RAND[, c("hispanic", "white", "black", "other",
#                            "unknown_race_eth")]))
# table(RAND$white, RAND$raracem, useNA = "ifany")
# table(RAND$rahispan, RAND$hispanic, useNA = "ifany")

RAND %<>% 
  #Drop the RAND rahispan variable (recoded as hispanic) and race variables
  dplyr::select(-c("rahispan", "raracem"))

#---- education ----
#Create education categories
RAND %<>% 
  mutate("ed_cat" = 
           case_when(raedyrs < 12 ~ 1, 
                     raedyrs == 12 ~ 2, 
                     raedyrs > 12 & raedyrs < 16 ~ 3, 
                     raedyrs == 16 ~ 4, 
                     raedyrs > 16 ~ 5)) 

# #Sanity check
# table(is.na(RAND$raedyrs))
# table(RAND$raedyrs)
# table(RAND$raedyrs, RAND$ed_cat, useNA = "ifany")

#---- smoking ----
RAND %<>% 
  mutate(smoker = RAND %>% 
           dplyr::select(paste0("r", number_waves, "smoken")) %>%
           apply(1, function(x) 
           {if (sum(is.na(x)) == length(x)){
             return(NA)} else {x[min(which(!is.na(x)))]}}))
# #Sanity check
# View(RAND[, c(paste0("r", number_waves, "smoken"), "smoker")])
# table(RAND$smoker, useNA = "ifany")

#Drop RAND's smoking variables
RAND %<>% dplyr::select(-paste0("r", number_waves, "smoken"))

#---- marital status ----

# Create marital status categories
mstat_mat <- RAND %>% select(contains("mstat"))

#Consolidate categories
mstat_mat[mstat_mat == 1 | mstat_mat == 2 | mstat_mat == 3] <- 1
mstat_mat[mstat_mat == 4 | mstat_mat == 5 | mstat_mat == 6 | 
            mstat_mat == 8] <- 2
mstat_mat[mstat_mat == 7] <- 3

RAND[, colnames(mstat_mat)] <- mstat_mat

for(wave in seq(1, 14)){
  level = 1
  for(cat in c("married_partnered", "not_married_partnered", "widowed")){
    RAND[, paste0("r", wave, cat)] <- 
      (RAND[, paste0("r", wave, "mstat")] == level)*1
    level = level + 1
  }
}

# #Sanity check
# View(RAND %>%
#        dplyr::select("r4mstat", "r4married_partnered",
#                      "r4not_married_partnered", "r4widowed"))
# View(RAND %>%
#        dplyr::select("r9mstat_impute", "r9married_partnered",
#                      "r9not_married_partnered", "r9widowed"))

#---- FLU SHOT ----
# Wide dataset
flusht_wide <- RAND %>%
  select(HHIDPN, contains("flusht")) %>%
  mutate(flusht_960004 = case_when(!is.na(r3flusht) & !is.na(r5flusht) & 
                                     !is.na(r7flusht) ~ 1,
                                   TRUE ~ 0),
         flusht_980206 = case_when(!is.na(r4flusht) & !is.na(r6flusht) & 
                                     !is.na(r8flusht) ~ 1,
                                   TRUE ~ 0),
         flusht_000408 = case_when(!is.na(r5flusht) & !is.na(r7flusht) & 
                                     !is.na(r9flusht) ~ 1,
                                   TRUE ~ 0),
         flusht_020610 = case_when(!is.na(r6flusht) & !is.na(r8flusht) & 
                                     !is.na(r10flusht) ~ 1,
                                   TRUE ~ 0),
         flusht_040812 = case_when(!is.na(r7flusht) & !is.na(r9flusht) & 
                                     !is.na(r11flusht) ~ 1,
                                   TRUE ~ 0),
         flusht_061014 = case_when(!is.na(r8flusht) & !is.na(r10flusht) & 
                                     !is.na(r12flusht) ~ 1,
                                   TRUE ~ 0),
         flusht_081216 = case_when(!is.na(r9flusht) & !is.na(r11flusht) & 
                                     !is.na(r13flusht) ~ 1,
                                   TRUE ~ 0))

# Check the distribution
# view(flusht_wide %>% filter(!is.na(r4flusht)))

flusht_wide %>% select(-HHIDPN) %>% apply(., 2, function(x) table(!is.na(x)))
flusht_wide %>% select(contains("flusht_")) %>%
  apply(., 2, function(x) table(x, useNA = "ifany"))

with(flusht_wide, table(flusht_960004, flusht_000408))
with(flusht_wide, table(flusht_000408, flusht_081216))
flusht_wide %>% filter(flusht_960004 == 1 & flusht_000408 == 1 &
                         flusht_081216 == 1) %>% 
  nrow()

# # Check the refreshers/new participants
# # 1998 (w4)
# with(RAND %>%
#        filter(r1iwstat == 0 & r2iwstat == 0 & r3iwstat == 0 & r4iwstat == 1),
#      table(!is.na(r4flusht)))
# # 2000
# with(RAND %>%
#        filter(r1iwstat == 0 & r2iwstat == 0 & r3iwstat == 0 & r4iwstat == 0 &
#                 r5iwstat == 1),
#      table(!is.na(r5flusht)))
# # 2002
# with(RAND %>%
#        filter(r1iwstat == 0 & r2iwstat == 0 & r3iwstat == 0 & r4iwstat == 0 &
#                 r5iwstat == 0 & r6iwstat == 1),
#      table(!is.na(r6flusht)))
# # 2004 (w7)
# with(RAND %>%
#        filter(r1iwstat == 0 & r2iwstat == 0 & r3iwstat == 0 & r4iwstat == 0 &
#                 r5iwstat == 0 & r6iwstat == 0 & r7iwstat == 1),
#      table(!is.na(r7flusht)))
# # 2006
# with(RAND %>%
#        filter(r1iwstat == 0 & r2iwstat == 0 & r3iwstat == 0 & r4iwstat == 0 &
#                 r5iwstat == 0 & r6iwstat == 0  & r7iwstat == 0  & r8iwstat == 1),
#      table(!is.na(r8flusht)))
# # 2008
# with(RAND %>%
#        filter(r1iwstat == 0 & r2iwstat == 0 & r3iwstat == 0 & r4iwstat == 0 &
#                 r5iwstat == 0 & r6iwstat == 0  & r7iwstat == 0  & 
#                 r8iwstat == 0 & r9iwstat == 1),
#      table(!is.na(r9flusht)))
# # 2010 (w10)
# with(RAND %>%
#        filter(r1iwstat == 0 & r2iwstat == 0 & r3iwstat == 0 & r4iwstat == 0 &
#                 r5iwstat == 0 & r6iwstat == 0  & r7iwstat == 0  & 
#                 r8iwstat == 0 & r9iwstat == 0 & r10iwstat == 1),
#      table(!is.na(r10flusht)))
# # 2012
# with(RAND %>%
#        filter(r1iwstat == 0 & r2iwstat == 0 & r3iwstat == 0 & r4iwstat == 0 &
#                 r5iwstat == 0 & r6iwstat == 0  & r7iwstat == 0  & 
#                 r8iwstat == 0 & r9iwstat == 0 & r10iwstat == 0 &
#                 r11iwstat == 1),
#      table(!is.na(r11flusht)))
# # 2014
# with(RAND %>%
#        filter(r1iwstat == 0 & r2iwstat == 0 & r3iwstat == 0 & r4iwstat == 0 &
#                 r5iwstat == 0 & r6iwstat == 0  & r7iwstat == 0  & 
#                 r8iwstat == 0 & r9iwstat == 0 & r10iwstat == 0 &
#                 r11iwstat == 0 & r12iwstat == 1),
#      table(!is.na(r12flusht)))
# # 2016 (w13)
# with(RAND %>%
#        filter(r1iwstat == 0 & r2iwstat == 0 & r3iwstat == 0 & r4iwstat == 0 &
#                 r5iwstat == 0 & r6iwstat == 0  & r7iwstat == 0  & 
#                 r8iwstat == 0 & r9iwstat == 0 & r10iwstat == 0 &
#                 r11iwstat == 0 & r12iwstat == 0 & r13iwstat == 1),
#      table(!is.na(r13flusht)))
# # 2018
# with(RAND %>%
#        filter(r1iwstat == 0 & r2iwstat == 0 & r3iwstat == 0 & r4iwstat == 0 &
#                 r5iwstat == 0 & r6iwstat == 0  & r7iwstat == 0  & 
#                 r8iwstat == 0 & r9iwstat == 0 & r10iwstat == 0 &
#                 r11iwstat == 0 & r12iwstat == 0 & r13iwstat == 0 & 
#                 r14iwstat == 1),
#      table(!is.na(r14flusht)))

seq_nonNAs <- function(x){
  y <- rle(x) # un length encoding, return values and corresponding lengths
  seq_nonNAs <-  case_when(max(y$lengths[!is.na(y$values)]) > 1 ~ 1,
                           sum(is.na(y$values)) == length(y$values) ~ NA_real_,
                           TRUE ~ 0)
  return(seq_nonNAs)
}
# # Test
# test <- c(0, NA, 1, NA, 1, 1, NA, 0, NA, 0, 1, NA, NA)
# test <- c(0, NA, 0, NA, NA)
# seq_nonNAs(test)

# (y = rle(test))
# max(y$lengths[y$values == 1 & !is.na(y$values)])

# Check the every 4 year rule
flusht_wide %>%
  select(-HHIDPN, -contains("flusht_")) %>%
  apply(., 1, seq_nonNAs) %>%
  table(., useNA = "ifany")
# 9202 participants have sequential flu shot status.

# Long dataset
flusht_long <- RAND %>%
  select(HHIDPN, contains("flusht")) %>%
  pivot_longer(!HHIDPN,
               names_to = "wave",
               names_pattern = "r(.*)flusht",
               values_to = "flusht") %>%
  mutate(wave = as.numeric(wave)) %>%
  na.omit()

# Num of taking the flu shot
flusht_times <- flusht_long %>%
  group_by(HHIDPN) %>%
  dplyr::summarise(n_flusht = sum(flusht))

summarytools::freq(flusht_times$n_flusht)

# Num of reported flu shot
flusht_report <- flusht_long %>%
  dplyr::count(HHIDPN) %>%
  dplyr::rename(n_flusht_wave = n)

summarytools::freq(flusht_report$n_flusht_wave)

# The trajectory of flu shot status among (100 random participants)
set.seed("62283")
random_id <- sample(flusht_times$HHIDPN, 100, replace = T)

flusht_long %>%
  # filter(HHIDPN %in% random_id) %>%
  arrange(wave) %>%
  ggplot(aes(x = wave, y = flusht, group = as.factor(HHIDPN),
             color = as.factor(HHIDPN))) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  scale_x_continuous(limits = c(3, 14), 
                     breaks = seq(3, 14, by = 1)) +
  theme(legend.position = "none")

#---- Filter to people have at least 2 wave flu shot status ----
flu_2_id <- flusht_report %>% filter(n_flusht_wave >= 2) %>%
  select(HHIDPN) %>% pull()

flu_2_3_8_id <- flusht_long %>%
  filter(wave <= 8) %>%
  dplyr::count(HHIDPN) %>%
  dplyr::rename(n_flusht_wave = n) %>% filter(n_flusht_wave >= 2) %>%
  select(HHIDPN) %>% pull()

flusht_wide_2 <- RAND %>%
  filter(HHIDPN %in% flu_2_id) %>%
  select(HHIDPN, contains("flusht")) %>%
  left_join(flusht_times, by = "HHIDPN") %>%
  left_join(flusht_report, by = "HHIDPN") %>%
  mutate(prop_flusht = n_flusht/n_flusht_wave)

# People at least 2 waves: how many people didn’t change their status
summarytools::freq(flusht_wide_2$prop_flusht)

# Prop of flu shot/all the available flu shot status (people at least two waves)
hist(flusht_wide_2$prop_flusht)

flusht_wide_2_3_8 <- RAND %>%
  filter(HHIDPN %in% flu_2_3_8_id) %>%
  select(HHIDPN, contains("flusht")) %>%
  left_join(flusht_times, by = "HHIDPN") %>%
  left_join(flusht_report, by = "HHIDPN") %>%
  mutate(prop_flusht = n_flusht/n_flusht_wave)

summarytools::freq(flusht_wide_2_3_8$prop_flusht)

#---- Save the data ----
RAND %<>% left_join(flusht_report, by = "HHIDPN") %>% 
  left_join(flusht_times, by = "HHIDPN")

save(RAND, file = here::here("data", "analysis_data", "rand_flusht.RData"))
