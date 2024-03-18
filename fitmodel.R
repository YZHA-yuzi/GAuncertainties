############################################################################
# This R script is used to illustrate how to fit the proposed model
# using self-defined R functions based on a toy dataset.
# NOTEs: Please make sure the working directory is set to the folder where 
# the R script "FUNs.R" is located.
# Author: yuzi.zhang@emory.edu
############################################################################

library(Rcpp)
library(RcppArmadillo)
library(dplyr) # left_join
library(tidyr) # gather
library(BayesLogit) # rpg
library(Matrix) # Diagonal
library(MASS) # mvrnorm
library(splines)

## load self-defined functions 
sourceCpp("FUNs_rcpp.cpp")
source("FUNs.R")

## read in the toy dataset ###
load("./data/Data_toy.rda") # dat.toy
# This toy dataset contains 5000 observations, 
# and 6 variables (ID, CENSUS, gest_error_1, gest_error_2, BIRTHDATE, SMOKE)
# ID: the id of the observation
# CENSUS: 11 digits FIPS code for the census tract 
# (sampled from the our real data and will be used to match with the exposure data)
# gest_error_1, gest_error_2: l-th error-prone gestational age, l = 1, 2
# BIRTHDATE: birth date
# SMOKE: smoking status 

dim(dat.toy)
head(dat.toy) 

### read in exposure data ###
load("./data/Data_exp_toy.rda") # dat.exp.toy
dim(dat.exp.toy)
head(dat.exp.toy)
## convert the the exposure data to the wide format
geoid.vec <- unique(dat.exp.toy$GEOID10)
dat.exp.wide.list <- list()
for(i in 1:length(geoid.vec)){
  dat.exp.sub <- subset(dat.exp.toy, GEOID10 == geoid.vec[i])
  dat.exp.sub$MV7 <- zoo::rollmean(dat.exp.sub$ozone, 
                                   k = 7, na.pad = T, align = "left")
  # Week 25 to Week 44
  index <- c(mapply(seq, 1:nrow(dat.exp.sub), 20*7 + (1:nrow(dat.exp.sub)) - 1,
                    MoreArgs = list(by = 7))) 
  dat.exp.wide <- matrix(dat.exp.sub$MV7[index], ncol = 20, byrow = T)
  dat.exp.wide <- data.frame(dat.exp.sub$GEOID10, dat.exp.sub$date,
                             dat.exp.wide)
  colnames(dat.exp.wide) <- c("GEOID", "date", paste0("Wk", 25:44))
  dat.exp.wide.list[[i]] <- dat.exp.wide
  if(i%%100 == 0){cat(i, ",")}
}
# "date" in dat.exp.wide: first date of week 25
dat.exp <- do.call(rbind.data.frame, dat.exp.wide.list)
dat.exp$GEOID <- as.character(dat.exp$GEOID)


### fit the proposed model based on the toy dataset ###
### 25 MCMC iterations take about 20s on a Mac running M2 chip
# nMCMC = 25000
# burn_in = 20000
nMCMC = 10
burn_in = 5
start = proc.time()[3]
print(paste0("The number of MCMC iterations is ", nMCMC))
cov.sim.0 <- c("factor(SMOKE)")
# print(paste(cov.sim.0, collapse = ":"))
fit <- fit.sim(dat = dat.toy, dat.exp = dat.exp,
               niter = nMCMC, burn_in =  burn_in, 
               seed = 12345,
               exp.name = "trim3.cum",
               get.Xdesign = get.Xdesign.tvarying,
               cov.name.vec = list(c("factor(Week)", cov.sim.0),
                                   c("factor(Week)", cov.sim.0)),
               date.start = as.Date("2010-01-01"), 
               date.end = as.Date("2010-12-31"),
               season = c(F, F), 
               concep.splines = c(F, F), df = c(NULL, NULL),
               exposure = c(T, T),
               dist.name = "1st", 
               ntol = 3)
proc.time()[3] - start

root.save = "/Users/zhangyuzi/Documents/Howard/ErrorOutcome/Manuscript/Rcodes/results/"
save(fit, file = paste0(root.save, "Res_toy.rda"))

### Outputs from the self-defined R function fit.sim
names(fit)
## fit = a list of length 8:
# (1) beta: a list of length 2 contains samples (before burn-in samples is discarded)
# of regression coefficients in the two discrete hazard models; 
# the 1st element is for the hazard model for preterm birth
# the 2nd element is for the hazard model for term birth
# (2) parm: a matrix contains samples (before burn-in samples is discarded) of
# parameters included in exponential covariance functions used in the outcome 
# misclassification model; 
# theta1_l = sigma_l^2, theta2_l = phi_l, for l=1, 2, in eqn.(1) in the manuscript  
# (3) gestwk: a matrix contains samples (before burn-in samples is discarded) of
# true gestational ages
# (4) alpha: a list of length 2 with each element further contains a list of 
# length 18. The list nested in the list "alpha" 
# contains random effects alpha_{l,kj} given in eqn. (1) in the manuscript
# (5) accept.mat: a matrix contains binary indicators of whether the proposed 
# value was accepted in M-H algorithms used for updating phi_l, l = 1, 2
# (6) nobs.vec: a vector contains the 
# number of observations included into the study population of interest 
# at each MCMC iteration
# (7) mis.parm: contains individual-specific parameters summarizing 
# misclassifications of gestational ages in term of identifying preterm birth.
# a list of length 2 with each element further contains a 
# list of length M, where M is the number of posterior samples after 
# the burn-in period is discarded.
# specifically, for each observation, we compute
# SE = sensitivity = Pr(Y_li < 37 | T_i < 37)
# SP = specificity = Pr(Y_li >= 37 | T_i >= 37)
# PPV = positive predictive value = Pr(T_i < 37 | Y_li < 37)
# NPV = negative predictive value = Pr(T_i >= 37 | Y_li >= 37)
# (8) hazard: a list of length M, where M is the number of posterior samples after 
# the burn-in period is discarded; m-th element in the list "hazard" contains
# a matrix of dimension n by 18, each row represents 
# individual-specific probability Pr(Ti = t), for t = 27, ..., 44, 
# i = 1, ..., n (n = the number of observations)

# We also describe how to produce simulation results included in the paper in this document. 

