# TMLE test for lmtp package
library(haven)
library(lmtp)
policy_0 <- function(data, trt){
  if(trt == "x0") {
    return(rep(0, each = length(data[[trt]])))
  }else{
    # this should not be used as the treatment is not the baseline treatment
    data[[trt]]
  }
}

policy_1 <- function(data, trt){
  if(trt == "x0") {
    return(rep(1, each = length(data[[trt]])))
  }else{
    # this should not be used as the treatment is not the baseline treatment
    data[[trt]]
  }
}

#---- load the data ----
ucla211epi <- read_sas('/Users/yingyanwu/Library/CloudStorage/OneDrive-Personal/Epidemiology Studies/6th year/EPI M211/programming/R/ucla211epi.sas7bdat')


#---- Natural course ----
nc <- lmtp_tmle(ucla211epi,
                trt = c("x0", "x1"),
                outcome = "y",
                baseline = "z0",
                time_vary = list("z0", "z1"),
                outcome_type = "binomial",
                mtp = FALSE,
                folds = 1)

#---- Sustained treatment effect ----
binary_on <- lmtp_tmle(ucla211epi,
                       trt = c("x0", "x1"),
                       outcome = "y",
                       baseline = "z0",
                       time_vary = list("z0", "z1"),
                       shift = static_binary_on,
                       outcome_type = "binomial",
                       mtp = FALSE, 
                       folds = 1)
# Estimate: 0.4158
# Std. error: 0.0034
# 95% CI: (0.4091, 0.4224)

binary_off <- lmtp_tmle(ucla211epi,
                        trt = c("x0", "x1"),
                        outcome = "y",
                        baseline = "z0",
                        time_vary = list("z0", "z1"),
                        # cens = ,
                        shift = static_binary_off,
                        outcome_type = "binomial",
                        mtp = FALSE, 
                        folds = 1)
# Population intervention estimate
# Estimate: 0.2208
# Std. error: 0.0031
# 95% CI: (0.2147, 0.2268)

lmtp_contrast(binary_on, ref = binary_off, type = c("additive"))
# Match with M211 Y11 - Y00: 0.069+0.0745+0.051=0.195

#---- Baseline treatment effect ----
# Truth
glm(y ~ x0 + z0, data = ucla211epi, family = binomial(link = "identity"))
# (Intercept)           x0           z0  
# 0.22680      0.13900      0.02384 

# LMTP
baseline_on <- lmtp_tmle(ucla211epi,
                         trt = c("x0", "x1"),
                         outcome = "y",
                         baseline = "z0",
                         time_vary = list("z0", "z1"),
                         # cens = ,
                         shift = policy_1,
                         outcome_type = "binomial",
                         mtp = FALSE, 
                         folds = 1)

baseline_off <- lmtp_tmle(ucla211epi,
                          trt = c("x0", "x1"),
                          outcome = "y",
                          baseline = "z0",
                          time_vary = list("z0", "z1"),
                          # cens = ,
                          shift = policy_0,
                          outcome_type = "binomial",
                          mtp = FALSE, 
                          folds = 1)
# Population intervention estimate
# Estimate: 0.3134
# Std. error: 0.0019
# 95% CI: (0.3097, 0.317)

lmtp_contrast(baseline_on, ref = baseline_off, type = c("additive"))
# X0's effect on Y: 0.139

#---- Test on two time points outcome ----
policy_0 <- function(data, trt){
  if(trt == "X3") {
    return(rep(0, each = length(data[[trt]])))
  }else{
    # this should not be used as the treatment is not the baseline treatment
    data[[trt]]
  }
}

policy_1 <- function(data, trt){
  if(trt == "X3") {
    return(rep(1, each = length(data[[trt]])))
  }else{
    # this should not be used as the treatment is not the baseline treatment
    data[[trt]]
  }
}

##---- Modified 211 data ----
ucla211epi_2 <- ucla211epi %>%
  rename(X3 = x0, L3 = z1, Y3 = x1, Y4 = y) %>%
  mutate(X4 = X3)

##---- Y4|X3 ----
Y4X3_1 <- lmtp_tmle(ucla211epi_2,
                    trt = c("X3", "X4"),
                    outcome = "Y4",
                    baseline = "z0",
                    time_vary = list("z0", "L3"),
                    # cens = ,
                    shift = policy_1,
                    outcome_type = "binomial",
                    mtp = FALSE, 
                    folds = 1)
Y4X3_0 <- lmtp_tmle(ucla211epi_2,
                    trt = c("X3", "X4"),
                    outcome = "Y4",
                    baseline = "z0",
                    time_vary = list("z0", "L3"),
                    # cens = ,
                    shift = policy_0,
                    outcome_type = "binomial",
                    mtp = FALSE, 
                    folds = 1)
ucla211epi_11 <- ucla211epi_2 %>%
  filter(X3 == 1)
Y4X3_1_2 <- lmtp_tmle(ucla211epi_11,
                    trt = c("X3", "X4"),
                    outcome = "Y4",
                    baseline = "z0",
                    time_vary = list("z0", "L3"),
                    # cens = ,
                    shift = policy_1,
                    outcome_type = "binomial",
                    mtp = FALSE, 
                    folds = 1)

Y4X3_nc <- lmtp_tmle(ucla211epi_2,
                     trt = "X3",
                     outcome = "Y4",
                     baseline = "z0",
                     time_vary = list("L3"),
                     outcome_type = "binomial",
                     mtp = FALSE, 
                     folds = 1)

lmtp_contrast(Y4X3_1, ref = Y4X3_0, type = c("additive"))
# X3's effect on Y4 (X0's effect on Y): is and should be 0.139

##---- Y4|X3, X4 ----
Y4_X3X4_1 <- lmtp_tmle(ucla211epi_2,
                       trt = c("X3", "X4"),
                       outcome = "Y4",
                       baseline = "z0",
                       time_vary = list("z0", "L3"),
                       shift = static_binary_on,
                       outcome_type = "binomial",
                       mtp = FALSE, 
                       folds = 1)
Y4_X3X4_1_2 <- lmtp_tmle(ucla211epi_11,
                       trt = c("X3", "X4"),
                       outcome = "Y4",
                       baseline = "z0",
                       time_vary = list("z0", "L3"),
                       shift = static_binary_on,
                       outcome_type = "binomial",
                       mtp = FALSE, 
                       folds = 1)

Y4_X3X4_0 <- lmtp_tmle(ucla211epi_2,
                       trt = c("X3", "X4"),
                       outcome = "Y4",
                       baseline = "z0",
                       time_vary = list("z0", "L3"),
                       shift = static_binary_off,
                       outcome_type = "binomial",
                       mtp = FALSE, 
                       folds = 1)
ucla211epi_00 <- ucla211epi_2 %>%
  filter(X3 == 0, X4 == 0)
Y4_X3X4_0_2 <- lmtp_tmle(ucla211epi_00,
                         trt = c("X3", "X4"),
                         outcome = "Y4",
                         baseline = "z0",
                         time_vary = list("z0", "L3"),
                         shift = static_binary_off,
                         outcome_type = "binomial",
                         mtp = FALSE, 
                         folds = 1)
Y4_X3X4_nc <- lmtp_tmle(ucla211epi_2,
                     trt = c("X3", "X4"),
                     outcome = "Y4",
                     baseline = "z0",
                     time_vary = list("z0", "L3"),
                     outcome_type = "binomial",
                     mtp = FALSE, 
                     folds = 1)

lmtp_contrast(Y4_X3X4_1, ref = Y4_X3X4_0, type = c("additive"))
# joint effect of X3 X4's effect on Y4 is 0.139

##---- Self Simulated data ----
set.seed(62283)
n = 25000
data_test <- tibble(
  Z = rbinom(n, size = 1, prob = 0.6),
  X3 = rbinom(n, 1, prob = 0.3 + Z*0.1),
  Y3 = rbinom(n, 1, prob = (0.4 + Z*0.01 + X3*0.3)),
  L3 = rbinom(n, 1, prob = (0.1 + Z*0.01 + X3*0.1)),
  Y4 = rbinom(n, 1, prob = (0.1 + Z*0.01 + X3*0.4 + Y3*0.075 + L3*0.1)),
  X4 = X3)

###---- Y4|X3 ----
Y4X3_1 <- lmtp_tmle(data_test,
                    trt = c("X3", "X4"),
                    outcome = "Y4",
                    baseline = "Z",
                    time_vary = list("Z", "L3"),
                    # cens = ,
                    shift = policy_1,
                    outcome_type = "binomial",
                    mtp = FALSE, 
                    folds = 1)

Y4X3_0 <- lmtp_tmle(data_test,
                    trt = c("X3", "X4"),
                    outcome = "Y4",
                    baseline = "Z",
                    time_vary = list("Z", "L3"),
                    # cens = ,
                    shift = policy_0,
                    outcome_type = "binomial",
                    mtp = FALSE, 
                    folds = 1)

lmtp_contrast(Y4X3_1, ref = Y4X3_0, type = c("additive"))
# X3's effect on Y4: 0.431

###---- Y4|X3, X4 ----
Y4_X3X4_1 <- lmtp_tmle(data_test,
                       trt = c("X3", "X4"),
                       outcome = "Y4",
                       baseline = "Z",
                       time_vary = list("Z", "L3"),
                       shift = static_binary_on,
                       outcome_type = "binomial",
                       mtp = FALSE, 
                       folds = 1)

Y4_X3X4_0 <- lmtp_tmle(data_test,
                        trt = c("X3", "X4"),
                        outcome = "Y4",
                        baseline = "Z",
                        time_vary = list("Z", "L3"),
                        shift = static_binary_off,
                        outcome_type = "binomial",
                        mtp = FALSE, 
                        folds = 1)
lmtp_contrast(Y4_X3X4_1, ref = Y4_X3X4_0, type = c("additive"))
# X3's effect on Y4: 0.431