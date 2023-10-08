
p_load("lme4")
# detectable effect analysis

options(scipen = 999)
#---- Aim 2 cognition----
load(here::here("data", "analysis_data", "hrs_K.RData"))
years <- str_sub(as.character(seq(2008, 2018, 2)), start = 3)
hrs_long <- hrs_edu_mem %>%
  filter(r9age_y_int > 50) %>%
  select(HHIDPN, cog_hs, cog_sum_08, social_hs, social_sum_08, paste0("memimp", years)) %>%
  pivot_longer(cols = !HHIDPN:social_sum_08,
               names_to = c(".value", "year"),
               names_pattern = "(memimp)(.*)") %>%
  mutate(year = as.numeric(year))

model0 <- lmer(memimp ~ year + (1 + year|HHIDPN), data = hrs_long)
summary(model0)
# SE of slope = 0.004
# varratio
# 0.042653/0.0004144^2

model1_cog <- lmer(memimp ~ cog_hs*year + cog_sum_08 + (1 + year|HHIDPN),
                   hrs_long,
                   control = lmerControl(optimizer = "Nelder_Mead"))
summary(model1_cog)
car::Anova(model1_cog)

model1_social <- lmer(memimp ~ social_hs*year + social_sum_08 + (1 + year|HHIDPN),
                      hrs_long,
                      control = lmerControl(optimizer = "Nelder_Mead"))
summary(model1_social)
car::Anova(model1_social)

p_load(sjPlot)
tab_model(list(model0, model1_cog, model1_social), show.ci = F)

power.long <- function(n, d, C, vare, varu, alpha){
  ncp <- d/sqrt(4*(vare/(varu*C)+1)/n)
  crit <- qt(1-alpha/2, n-2)
  power <- 1-pt(crit, n-2, ncp)
  return(power)
}
# variance of slope 
seu <- 0.004 # From above
varu <- seu^2
x <- 0:10 # Time range 2008 - 2018: 10 years
ct <- poly(x, degree=1)*sqrt(sum((x-mean(x))^2))
C <- sum(ct[1:length(x),]^2)
C

# for cog
hrs_edu_mem %>% filter(r9age_y_int > 50 & 
                         (!is.na(cog_hs) | !is.na(cog_sum_08)) &
                         (!is.na(memimp08) |!is.na(memimp10)|!is.na(memimp12)|
                            !is.na(memimp14)|!is.na(memimp16)|!is.na(memimp18))) %>% nrow()
# 9271
hrs_edu_mem %>% filter(r9age_y_int > 50 & 
                         (!is.na(social_hs) | !is.na(social_sum_08)) &
                         (!is.na(memimp08) |!is.na(memimp10)|!is.na(memimp12)|
                            !is.na(memimp14)|!is.na(memimp16)|!is.na(memimp18))) %>% nrow()
# 9470

# Cog activities
power.long(9271, 0.00024/0.004, C, seu^2*6, seu^2, 0.05)
power.long(9271, 0.000242/0.004, C, seu^2*7, seu^2, 0.05)
power.long(9271, 0.000243/0.004, C, seu^2*8, seu^2, 0.05)
power.long(9271, 0.000243/0.004, C, seu^2*9, seu^2, 0.05)
power.long(9271, 0.000244/0.004, C, seu^2*10, seu^2, 0.05)
power.long(9271, 0.000245/0.004, C, seu^2*11, seu^2, 0.05)
power.long(9271, 0.000246/0.004, C, seu^2*12, seu^2, 0.05)
# The values for beta comes from various round of running the power.long
aim2 <- tibble(
  beta = c(0.00024, 0.000242, 0.000243, 0.000243, 0.000244, 0.000245, 0.000246),
  # beta = seq(0.0002, 0.00025, by = 0.00001),
  d = beta/seu, # from a study
  varratio = 6:12,
  n = 9271,
  vare = varu*varratio,
  C = C)

aim2 %<>%
  mutate(power = power.long(n, d, C, vare, varu, 0.05))
aim2 %>% print(n = Inf)

# Social activities
power.long(9470, 0.000237/0.004, C, seu^2*6, seu^2, 0.05)
power.long(9470, 0.000237/0.004, C, seu^2*7, seu^2, 0.05)
power.long(9470, 0.000239/0.004, C, seu^2*8, seu^2, 0.05)
power.long(9470, 0.000240/0.004, C, seu^2*9, seu^2, 0.05)
power.long(9470, 0.000241/0.004, C, seu^2*10, seu^2, 0.05)
power.long(9470, 0.000242/0.004, C, seu^2*11, seu^2, 0.05)
power.long(9470, 0.000243/0.004, C, seu^2*12, seu^2, 0.05)
# The values for beta comes from various round of running the power.long
aim2 <- tibble(
  beta = c(0.000237, 0.000237, 0.000239, 0.000240, 0.000241, 0.000242, 0.000243),
  # beta = seq(0.0002, 0.00025, by = 0.00001),
  d = beta/seu, # from a study
  varratio = 6:12,
  n = 9470,
  vare = varu*varratio,
  C = C)

aim2 %<>%
  mutate(power = power.long(n, d, C, vare, varu, 0.05))
aim2 %>% print(n = Inf)

#----- Aim 2: dementia ----
# Using SAS proc power
summary(hrs_edu_mem_dem %>% select(contains(c("cog_", "social_"))))
hrs_edu_mem_dem %>% filter(cog_hs >= 2 | cog_sum_08 >= 4) %>% nrow() # 5135
hrs_edu_mem_dem %>% filter(!(cog_hs >= 2 | cog_sum_08 >= 4)) %>% nrow() # 1558

with(hrs_edu_mem_dem, table(cog_hs, useNA = "ifany"))
with(hrs_edu_mem_dem, table(cog_sum_08, useNA = "ifany"))
#Ref proportion 0.15
with(hrs_edu_mem_dem %>% filter(cog_hs < 2), table(expert_dem, useNA = "ifany"))/(3093+552)
# Index proportion: 0.10
with(hrs_edu_mem_dem %>% filter(cog_hs >= 2), table(expert_dem, useNA = "ifany"))/(1294+150)

# RR = 0.667 from the crude model
0.10/0.15
