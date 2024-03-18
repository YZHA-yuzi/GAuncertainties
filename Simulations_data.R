############################################################################
# This R script is used to generate simulated datasets used in 
# simulation studies included in Section 4
# NOTEs: (1) Please make sure the working directory is set to the folder where 
# the R script "FUNs.R" and "FUNs_rcpp.cpp" are located.
# (2) Two exposures randomly sampled from the real data were provided in 
# the folder "data": Xdesign_smk.rda and Xdesign_ozone.rda
# (3) Parameters used for generating simulation datasets were also provided in
# the foler "data": "Tab_Supp_simssetup.xlsx"
# (4) Detailed steps of generating simulated datasets can be found in 
# Section S2.1 of Supplmentary Materials
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
library(readxl) 

if(!file.exists("./data")){
  dir.create("data")
}

## load self-defined functions 
source("FUNs.R")
sourceCpp("FUNs_rcpp.cpp")

## read in exposures used for simulating data
load("./data/Xdesign_smk.rda") # X.design.fix.all
load("./data/Xdesign_ozone.rda") # dat.sim.exp.all

## BEGIN GENERATING SIMULATED DATASETS ##
dat.sim.all <- list()
hazard.sim.all <- list()
parm.sim.all <- list()
## loop over different scenarios s = 1, 2, 3 
for(s in 1:3){
  
  ## prepare time-invariant exposures for two hzard models
  X.design.fix.0 <- X.design.fix.all[[s]]
  X.design.fix <- list(X.design.fix.0, X.design.fix.0)
  ## prepare time-varying exposures for two hazard models 
  dat.sim.exp <- dat.sim.exp.all[[s]]
  
  dat.sim.0 <- data.frame("SMOKE" = as.numeric(X.design.fix.0[,1]),
                          "ID" = 1:nrow(X.design.fix.0),
                          "fake" = 1)
  
  if(s %in% c(1, 3)){
    set.seed(1345)
  }else if(s == 2){
    set.seed(12345)
  }
  
  beta.week.true.list <- list()
  beta.fix.true.list <- list()
  beta.o3.true <- rep(NA, 2)
  ## read in beta coefficients in the hazard model for preterm birth
  inter.i <- readxl::read_xlsx(path = "./data/Tab_Supp_simssetup.xlsx",
                               sheet = "hazard_preterm")
  beta.week.true.list[[1]] <- 
    as.numeric(inter.i$value[inter.i$parm %in% paste0("gamma_0,",27:36)])
  beta.fix.true.list[[1]] <- 
    as.numeric(inter.i$value[inter.i$parm %in% "smk"])
  beta.o3.true[1] <- as.numeric(inter.i$value[inter.i$parm %in% "ozone"])
  
  ## read in beta coefficients in the hazard model for preterm birth
  inter.i <- readxl::read_xlsx(path = "./data/Tab_Supp_simssetup.xlsx",
                               sheet = "hazard_fullterm")
  beta.week.true.list[[2]] <- 
    as.numeric(inter.i$value[inter.i$parm %in% paste0("psi_0,",37:44)])
  beta.fix.true.list[[2]] <- 
    as.numeric(inter.i$value[inter.i$parm %in% "smk"])
  beta.o3.true[2] <- as.numeric(inter.i$value[inter.i$parm %in% "ozone"])
  
  multi.prob.forsim <- list()
  ## read in outcome mislcassification probabilities ##
  for(ll in 1:2){
    inter.i <- readxl::read_xlsx(path = "./data/Tab_Supp_simssetup.xlsx",
                                 sheet = paste0("S", s, "_pi_l", ll))
    inter.i <- t(as.matrix(inter.i[,-1]))
    rownames(inter.i) <- paste0("w", 27:44)
    colnames(inter.i) <- paste0("t", 27:44)
    multi.prob.forsim[[ll]] <- inter.i
  }

  re.sim.list <- gen.dat.sim(dat.sim.0 = dat.sim.0,
                             X.design.fix = X.design.fix,
                             dat.sim.exp = dat.sim.exp, 
                             nsims = 100, exp.name = "trim3.cum", 
                             beta.week.true.list = beta.week.true.list, 
                             scale.factor = 1,
                             beta.true = beta.o3.true, 
                             beta.fix.true.list = beta.fix.true.list, 
                             beta.season.true = list(NULL, NULL),
                             beta.splines.true = list(NULL, NULL),
                             season = c(F, F), concep.splines = c(F, F), 
                             df = c(NULL, NULL), 
                             prob.list = multi.prob.forsim,
                             val.name = "o3")
  ## save simulated data
  dat.sim.all[[s]] <- re.sim.list[["data"]]
  ## save individual-specific Pr(Ti = t), for t = 27,...,44
  hazard.sim.all[[s]] <- re.sim.list[["hazard"]][[1]]
  
  ## compute individual-specific SE, SP, PPV, NPV in terms of being identified
  ## as preterm birth
  get.mis.parm <- function(den.w, hazard.vec){
    prob.joint.ind <- do.call(rbind, 
                              lapply(1:18, function(y) den.w[y, ]*hazard.vec))
    se.ind <- sum(prob.joint.ind[1:10, 1:10])/sum(prob.joint.ind[, 1:10])
    sp.ind <- sum(prob.joint.ind[11:18, 11:18])/sum(prob.joint.ind[, 11:18])
    ppv.ind <- sum(prob.joint.ind[1:10, 1:10])/sum(prob.joint.ind[1:10, ])
    npv.ind <- sum(prob.joint.ind[11:18, 11:18])/sum(prob.joint.ind[11:18, ])
    return(c(se.ind, sp.ind, ppv.ind, npv.ind))
  }
  den.w1 <- multi.prob.forsim[[1]]
  den.w2 <- multi.prob.forsim[[2]]
  parm.ind.w1 <- 
    do.call(rbind, lapply(1:nrow(hazard.sim.all[[s]]), function(z) 
      get.mis.parm(den.w = den.w1, 
                   hazard.vec = hazard.sim.all[[s]][z, ])))
  parm.ind.w2 <- 
    do.call(rbind, lapply(1:nrow(hazard.sim.all[[s]]), function(z) 
      get.mis.parm(den.w = den.w2, 
                   hazard.vec = hazard.sim.all[[s]][z, ])))
  
  parm.ind.w2 <- as.data.frame(parm.ind.w2)
  parm.ind.w1 <- as.data.frame(parm.ind.w1)
  colnames(parm.ind.w1) <- colnames(parm.ind.w2) <- c("SE", "SP", "PPV", "NPV")
  parm.list <- list("error1" = parm.ind.w1, "error2" = parm.ind.w2)
  parm.sim.all[[s]] <- parm.list
  
}
## under each simulation scenario, 100 replicates were generated
## save simulated datasets across S1-S3 simulation scenarios
save(dat.sim.all, file = paste0("./data/Data_sims.rda"))
## save true individual-specific Pr(Ti = t), t = 27,...,44
save(hazard.sim.all, file = paste0("./data/Data_hazard_sims.rda"))
## save true individual-specific SE, SP, PPV, NPV
save(parm.sim.all, file = paste0("./data/Data_misparm_sims.rda"))

