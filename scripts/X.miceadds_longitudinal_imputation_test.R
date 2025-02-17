set.seed(123)
plouf <- data.frame(ID = rep(LETTERS[1:10],each = 10), x = sample(10,100,replace = T), y = sample(10,100,replace = T))
plouf[sample(100,30), "x"] <- NA
set.seed(999)
plouf[sample(100,30), "y"] <- NA

# Baseline variables
plouf %<>%
  mutate(age = sample(18:65, 100, replace = TRUE),
    baseline = case_when(ID %in% c("A", "C", "H", "J") ~ rnorm(1, mean = 30, sd = 49),
                              ID %in% c("B", "I") ~ NA_real_,
                              TRUE ~ 0))
DataExplorer::plot_missing(plouf)

plouf %<>%
  select(ID, age, baseline, x, y)

#missing value imputation
where_test <- make.where(plouf)
whr <- data.frame(ID = rep(FALSE,100), age = is.na(plouf$age), 
                  baseline = is.na(plouf$baseline), x = rep(FALSE,100), y = is.na(plouf$y))
plouf.imp <- mice(plouf, m = 1, method="pmm", maxit=5, where = whr)
test <- complete(plouf.imp)
sum(is.na(test$y))


# longitudinal imputation
predM <- mice::make.predictorMatrix(data = plouf) 
impM <- mice::make.method(data = plouf) 

predM1 <- predM
predM1[, "ID"] <- 0 # do not predict anything with HHIDPN
predM1["ID", ] <- 0 # this step is optional
predM1[c("x", "y"), "ID"] <- -2 # set nesting variable for the tv vars

imputationFunction <- rep("pmm", length(c("baseline"))) %>% as.list()
names(imputationFunction) <- c("baseline")
imputationFunction

# define cluster variable, i.e. nested in HHIDPN
cluster_var <- rep("ID", length("baseline")) %>% as.list()
names(cluster_var) <- c("baseline")
cluster_var

test_imp_long <- mice(plouf, m = 1, maxit = 1,
                   predictorMatrix = predM1,
                   method = impM1,
                   where = whr,
                   imputationFunction = imputationFunction,
                   cluster_var = cluster_var)
