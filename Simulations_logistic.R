##############################################################################
# This R script is used to 
# fit regular logistic regression model based on simulated datasets 
# NOTEs: (1) Please make sure the working directory is set to the folder where 
# the R script "FUNs.R" and "FUNs_rcpp.cpp" are located.
# (2) Conventional analyses of preterm birth were conducted solely 
# relying on error-prone GAs for the purpose of comparison.
# Author: yuzi.zhang@emory.edu
###############################################################################

args <-  commandArgs(trailingOnly = TRUE)
sceindex = eval( parse(text=args[1]) )
## NOTEs: If you run this R script via RStudio, please specify 
## sceindex = index of simulation scenario, s = 1, 2, 3

print(paste0("S = ", sceindex))

library(Rcpp)
library(RcppArmadillo)
library(dplyr) # left_join
library(tidyr) # gather
library(BayesLogit) # rpg
library(Matrix) # Diagonal
library(MASS) # mvrnorm
library(splines)
library(readxl) 


if(!file.exists("./interres")){
  dir.create("interres")
}


## load self-defined functions 
sourceCpp("FUNs_rcpp.cpp")
source("FUNs.R")


## read in exposure data 
load("./data/Data_exp_toy.rda") # dat.exp.toy
# dim(dat.exp.toy)
# head(dat.exp.toy)
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


## read in simulated data 
load("./data/Data_sims.rda") # dat.sim.all
dat.sim.list <- dat.sim.all[[sceindex]]
# nsim = length(dat.sim.list)
nsim = 5
# A total of 100 replicates were included under each simulation scenario
# Each simulated dataset contains 5000 observations, 
# and 11 variables 
# The following 6 varaibes will be used when fitting the proposed model
# (ID, CENSUS, gest_error_1, gest_error_2, BIRTHDATE, SMOKE)
# ID: the id of the observation
# CENSUS: 11 digits FIPS code for the census tract 
# (sampled from the our real data and will be used to match with the exposure data)
# gest_error_1, gest_error_2: l-th error-prone gestational age, l = 1, 2
# BIRTHDATE: birth date
# SMOKE: smoking status 

week.risk.ptb <- 27:36; week.risk.birth <- 37:44
week.risk.list <- list(week.risk.ptb, week.risk.birth)
date.start = as.Date("2010-01-01")
date.end = as.Date("2010-12-31")

re.ptb.list <- re.birth.list <- list()
re.glm.all <- list()
ptb.glm.all <- list()
dat.exp.sim <- dat.exp
for(i in 1:nsim){
  dat.use <- dat.sim.list[[i]]
  
  concep.obs <- 
    data.frame("true" = get.concep(db = dat.use$BIRTHDATE, 
                                   ga_wk = dat.use$gest_true),
               "error1" = get.concep(db = dat.use$BIRTHDATE, 
                                     ga_wk = dat.use$gest_error_1),
               "error2" = get.concep(db = dat.use$BIRTHDATE, 
                                     ga_wk = dat.use$gest_error_2))
  keep.true <- date.start <= concep.obs$true & concep.obs$true <= date.end
  keep.error1 <- date.start <= concep.obs$error1 & 
    concep.obs$error1 <= date.end
  keep.error2 <- date.start <= concep.obs$error2 & 
    concep.obs$error2 <= date.end
  keep.mat <- cbind(keep.true, keep.error1, keep.error2)
  
  X.fix.full.ptb <- get.Xdesign.fix(dat = dat.use, 
                                    week.risk = week.risk.ptb, 
                                    cov.name.vec = c("factor(Week)",
                                                     "factor(SMOKE)"))
  X.fix.full.birth <- get.Xdesign.fix(dat = dat.use, 
                                      week.risk = week.risk.birth, 
                                      cov.name.vec = c("factor(Week)",
                                                       "factor(SMOKE)"))
  X.fix.full.list <- list(X.fix.full.ptb, X.fix.full.birth)
  
  ptb.glm.list <- list("overall" = list(), "ind" = list())
  fit.glm.all <- list()
  event.name <- c("PTB", "FULL")
  out.name <- c("gest_true", "gest_error_1", "gest_error_2")
  re.glm.i <- list()
  name.i <- c(paste0(event.name[1], "_", out.name),
              paste0(event.name[2], "_", out.name))
  count = 1
  for(l in 1:2){
    fit.glm.l <- list()
    for(ll in 1:3){
      fit.glm.l[[ll]] <- fit.initi(dat = dat.use, out.init = out.name[ll],
                                   exp.name = "trim3.cum",
                                   dat.exp = dat.exp.sim,
                                   cov.name.vec = c("factor(Week)",
                                                    "factor(SMOKE)"),
                                   season = F, concep.splines = F, df = NULL,
                                   date.start = as.Date("2010-01-01"),
                                   date.end = as.Date("2010-12-31"),
                                   exposure = T, event_name = event.name[l])
      re.glm.i[[count]] <- summary(fit.glm.l[[ll]])$coefficients
      re.glm.i[[count]] <- 
        data.frame(re.glm.i[[count]][,1:2],
                   "lci" = re.glm.i[[count]][,1] - 1.96*re.glm.i[[count]][,2],
                   "uci" = re.glm.i[[count]][,1] + 1.96*re.glm.i[[count]][,2])
      count = count + 1
      
      ## compute individual-specific probability of being identified as 
      ## preterm birth, i.e., Pr(T_i < 37) following steps outlined in 
      ## Section S3 in Supplementary Materials
      if(l == 1){
        
        inter.ll <- summary(fit.glm.l[[ll]])
        beta.glm.ll <- 
          as.data.frame(matrix(as.numeric(inter.ll$coefficients[,1]), nrow = 1))
        colnames(beta.glm.ll) <- rownames(inter.ll$coefficients)
        beta.cov.glm.ll <- inter.ll$cov.unscaled
        
        ## draw samples of MVN 
        beta.glm.post <- mvrnorm(1000, mu = as.numeric(beta.glm.ll),
                                 Sigma = beta.cov.glm.ll)
        beta.glm.post <- rbind(beta.glm.ll, beta.glm.post)
        ## compute individual-specific Pr(Ti <= 36)
        re.X.tvarying <- 
          get.Xdesign.tvarying(dat = dat.use,
                               dat.exp = dat.exp.sim,
                               t.curr = dat.use[[out.name[ll]]],
                               exp.name = "trim3.cum",
                               season = F, concep.splines = F, df = NULL,
                               date.start = as.Date("2010-01-01"),
                               date.end = as.Date("2010-12-31"),
                               week.risk = week.risk.list[[1]])
        ## design matrix of cumulative average of O3
        X.tvarying.full <- re.X.tvarying$exp.tvarying.long
        index.sel <- colnames(beta.glm.ll)%in%colnames(X.fix.full.list[[1]])
        
        ptb.ind.l <- matrix(NA, nrow = nrow(beta.glm.post), 
                            ncol = nrow(dat.use))
        for(kk in 1:nrow(beta.glm.post)){
          pred.comp1 <- as.numeric(X.tvarying.full[,3]*(beta.glm.post[kk, 1]))
          pred.comp2 <- as.numeric(
            X.fix.full.list[[1]]%*%
              matrix(as.numeric(beta.glm.post[kk, index.sel]), ncol = 1))
          pred.val.kk <- pred.comp1 + pred.comp2
          ### Compute Pr(ti = t) for each observation t \in 27:36
          h0.mat.1 <- matrix(expit(pred.val.kk), 
                             ncol = length(week.risk.ptb), 
                             byrow = FALSE)
          colnames(h0.mat.1) <- paste0("wk", week.risk.ptb)
          hazard.curr.mat.1 <- comp_hazard(h0.mat.1)
          ptb.ind.l[kk, ] <- as.numeric(rowSums(hazard.curr.mat.1))
        } # loop over samples of bets
        
        ptb.mat.l <- as.data.frame(matrix(NA, ncol = 4, nrow = nrow(dat.use)))
        colnames(ptb.mat.l) <- c("Est", "SE", "lci", "uci")
        ptb.mat.l$Est <- ptb.ind.l[1, ]
        ptb.mat.l$SE <- apply(ptb.ind.l[-1, ], 2, sd)
        ptb.mat.l$lci <- apply(ptb.ind.l[-1, ], 2, quantile, 0.025)
        ptb.mat.l$uci <- apply(ptb.ind.l[-1, ], 2, quantile, 0.975)
        ptb.glm.list[["ind"]][[ll]] <- ptb.mat.l
        
        ptb.vec.l <- rowMeans(ptb.ind.l[-1, ])
        ptb.glm.list[["overall"]][[ll]] <- 
          data.frame("Outcome" = out.name[ll],
                     "Est" = mean(ptb.ind.l[1, ]),
                     "SE" = sd(ptb.vec.l),
                     "lci" = quantile(ptb.vec.l, 0.025),
                     "uci" = quantile(ptb.vec.l, 0.975))
        
        ptb.bin.l <- mean(dat.use[[out.name[ll]]] <= 36)
        se.bin.l <- sqrt(ptb.bin.l*(1-ptb.bin.l)/nrow(dat.use))
        ptb.glm.list[["overall"]][[ll]][["Est.bin"]] <- ptb.bin.l
        ptb.glm.list[["overall"]][[ll]][["SE.bin"]] <- se.bin.l
        ptb.glm.list[["overall"]][[ll]][["lci.bin"]] <- ptb.bin.l - 1.96*se.bin.l
        ptb.glm.list[["overall"]][[ll]][["uci.bin"]] <- ptb.bin.l + 1.96*se.bin.l
      }
      
    } # END loop over different GA estimates
    names(fit.glm.l) <- out.name
    fit.glm.all[[l]] <- fit.glm.l
  }
  names(re.glm.i) <- name.i
  get.est.i <- function(x){
    re = summary(x)$coefficients
    sel = rownames(re) %in% c("trim3.cum", "factor(SMOKE)1")
    return(re[sel, ])
  }
  names(fit.glm.all) <- event.name
  re.ptb <- lapply(fit.glm.all[[1]], get.est.i)
  re.birth <- lapply(fit.glm.all[[2]], get.est.i)
  re.ptb.list[[i]] <- re.ptb
  re.birth.list[[i]] <- re.birth
  re.glm.i.1 <- lapply(re.glm.i[1:3],
                       function(x) x[rownames(x) %in%
                                       paste0("factor(Week)", 27:36), ])
  re.glm.i.2 <- lapply(re.glm.i[4:6],
                       function(x) x[rownames(x) %in%
                                       paste0("factor(Week)", 37:44), ])
  re.glm.all[[i]] <- c(re.glm.i.1, re.glm.i.2)
  
  ptb.glm.all[[i]] <- ptb.glm.list
  cat(i, ",")
  
} # END loop over simulations
names(re.glm.all) <- paste0("sim", 1:nsim)
names(ptb.glm.all) <- paste0("sim", 1:nsim)

beta1.ptb <- lapply(re.ptb.list,
                    function(y) do.call(rbind, lapply(y, function(x) x[1, ])))
names(beta1.ptb) <- paste0(paste0("sim", 1:nsim), "_ptb_o3")

beta1.birth <- lapply(re.birth.list,
                      function(y) do.call(rbind, lapply(y, function(x) x[1, ])))
names(beta1.birth) <- paste0(paste0("sim", 1:nsim), "_birth_o3")

beta2.ptb <- lapply(re.ptb.list,
                    function(y) do.call(rbind, lapply(y, function(x) x[2, ])))
names(beta2.ptb) <- paste0(paste0("sim", 1:nsim), "_ptb_smoke")

beta2.birth <- lapply(re.birth.list,
                      function(y) do.call(rbind, lapply(y, function(x) x[2, ])))
names(beta2.birth) <- paste0(paste0("sim", 1:nsim), "_birth_smoke")

beta.all <- c(beta1.ptb, beta1.birth, beta2.ptb, beta2.birth)
for(l in 1:length(beta.all)){
  beta.all[[l]] <- as.data.frame(beta.all[[l]])
  beta.all[[l]] <- data.frame(outcome = c("true", "error_1", "error_2"),
                              beta.all[[l]])
}
beta.int.all <- re.glm.all

res.glm.all <- list("beta.all" = beta.all, 
                    "beta.int.all" = re.glm.all,
                    "ptb" = ptb.glm.all)
save(res.glm.all,
     file = paste0("./interres/Res_glm_S", sceindex, ".rda"))
