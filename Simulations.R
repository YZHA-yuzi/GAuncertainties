##############################################################################
# This R script is used to fit the proposed model based on simulated datasets 
# NOTEs: (1) Please make sure the working directory is set to the folder where 
# the R script "FUNs.R" and "FUNs_rcpp.cpp" are located.
# Author: yuzi.zhang@emory.edu
###############################################################################

## NOTEs: If you run this R script via RStudio, please specificy 
## sceindex = index of simulation scenario, s = 1, 2, 3
## simindex = index of simulation within each simulation scenario,
## simindex = 1, ..., 100.
args <-  commandArgs(trailingOnly = TRUE)
sceindex = eval( parse(text=args[1]) )
simindex = eval( parse(text=args[2]) )

print(paste0("S = ", sceindex, "; sim = ", simindex))

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
dat.sim <- dat.sim.all[[sceindex]][[simindex]]
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


###### MODEL FITTING #####
## 1. fit the proposed model
### 25 MCMC iterations take about 20s on a Mac running M2 chip
nMCMC = 25000
burn_in = 20000
# nMCMC = 10
# burn_in = 5
start = proc.time()[3]
print(paste0("The number of MCMC iterations is ", nMCMC))
cov.sim.0 <- c("factor(SMOKE)")
fit <- fit.sim(dat = dat.sim, dat.exp = dat.exp,
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


## get posterior samples of true GA 
gest.post <- fit$gestwk[-c(1:burn_in), ]
############## BEGIN COMPUTE PTB RATE #######################
### compute point-wise Pr(Ti < 37)
ptb.hazard.mat <- do.call(rbind, lapply(fit$hazard, 
                                        function(x) rowSums(x[, 1:10])))
ptb.hazard <- sapply(fit$hazard, 
                     function(x) mean(rowSums(x[, 1:10])))

load(paste0("./data/Data_hazard_sims.rda"))
ptb.hazard.true <- mean(rowSums(hazard.sim.all[[sceindex]][,1:10]))
ptb.hazard.true.vec <- rowSums(hazard.sim.all[[sceindex]][,1:10])
ptb.mat <- data.frame("ptb_error_1" = mean(dat.sim$gest_error_1 < 37),
                      "ptb_error_2" = mean(dat.sim$gest_error_2 < 37),
                      "true" = ptb.hazard.true,
                      "est" = mean(ptb.hazard),
                      "sd" = sd(ptb.hazard),
                      "lci" = quantile(ptb.hazard, 0.025),
                      "uci" = quantile(ptb.hazard, 0.975))
# individual-specific 
ptb.mat.pw <- 
  data.frame("true_hazard" = ptb.hazard.true.vec,
             "est_hazard" = colMeans(ptb.hazard.mat),
             "sd_hazard" = apply(ptb.hazard.mat, 2, sd),
             "lci_hazard" = apply(ptb.hazard.mat, 2, quantile, 0.025),
             "uci_hazard" = apply(ptb.hazard.mat, 2, quantile, 0.975))

ptb.mat.pw[["include"]] <- 
  (ptb.mat.pw$lci_hazard <= ptb.mat.pw$true_hazard) & 
  (ptb.mat.pw$uci_hazard >= ptb.mat.pw$true_hazard)
ptb.list <- list("ptb.mat" = ptb.mat, "ptb.mat.pw" = ptb.mat.pw)

## save estimated overall PTB rate and individual-specific probability 
## of being identified as preterm birth (i.e., Pr(T_i < 37)). The later
## quantity will be used to compute the coverage for the overall PTB rate
## (i.e., the empirical coverage was computed by 
## averaging over simulations and individuals)
save(ptb.list, file = 
       paste0("./interres/Res_PTBrate_S", sceindex, "_sim", simindex, ".rda"))
############## END COMPUTE PTB RATE #######################

############## BEGIN COMPUTE SE, SP, PPV, NPV #######################
load(paste0("./data/Data_misparm_sims.rda")) # parm.sim.all
parm.list <- parm.sim.all[[sceindex]]
mis.parm.true.1 <- parm.list$error1
mis.parm.true.2 <- parm.list$error2
mis.parm.true.list <- list(mis.parm.true.1, mis.parm.true.2)
mis.parm.sum.list <- list("error1" = list(), "error2" = list())
for(ksim in 1:2){
  kparm.vec <- colnames(fit$mis.parm[[ksim]][[1]])
  for(kparm in 1:4){
    re.kk <- do.call(cbind, lapply(fit$mis.parm[[ksim]], 
                                   function(x) x[[kparm]]))
    re.mse.kk <- lapply(1:length(mis.parm.true.list[[ksim]][,kparm]),
                        function(x) 
                          (mean(re.kk[x, ]) - 
                             mis.parm.true.list[[ksim]][x,kparm])^2)
    re.bias.kk <- lapply(1:length(mis.parm.true.list[[ksim]][,kparm]),
                         function(x) 
                           mean((re.kk[x, ] - 
                                   mis.parm.true.list[[ksim]][x,kparm])))
    
    re.sum.kk <- data.frame(rowMeans(re.kk), apply(re.kk, 1, sd),
                            apply(re.kk, 1, quantile, 0.025),
                            apply(re.kk, 1, quantile, 0.975),
                            do.call(c, re.mse.kk), 
                            do.call(c, re.bias.kk))
    colnames(re.sum.kk) <- c(paste0(kparm.vec[kparm],".est"),
                             paste0(kparm.vec[kparm],".sd"),
                             paste0(kparm.vec[kparm],".lci"),
                             paste0(kparm.vec[kparm],".uci"),
                             paste0(kparm.vec[kparm],".mse"),
                             paste0(kparm.vec[kparm],".bias"))
    re.sum.kk <- data.frame(re.sum.kk, 
                            re.sum.kk[,paste0(kparm.vec[kparm],".lci")] <= 
                              mis.parm.true.list[[ksim]][,kparm] & 
                              re.sum.kk[,paste0(kparm.vec[kparm],".uci")] >= 
                              mis.parm.true.list[[ksim]][,kparm])
    colnames(re.sum.kk)[7] <- paste0(kparm.vec[kparm],".include")
    mis.parm.sum.list[[ksim]][[kparm]] <- re.sum.kk
    
  }
  names(mis.parm.sum.list[[ksim]]) <- kparm.vec
}
save(mis.parm.sum.list, file = 
       paste0("./interres/Res_misparm_S", sceindex, "_sim", simindex, ".rda"))
############## END COMPUTE SE, SP, PPV, NPV #######################


###### BEING COMPUTE WAIC for the outcome misclassification model ######
t.post.all <- fit$gestwk
ntol = 3
nweeks = 18; t.unique = 27:44
w.range.list <- list()
for(ll in 1:nweeks){
  t.curr <- t.unique[ll]
  w.range <- max(27, t.curr-ntol):min(44, t.curr+ntol)
  w.range.list[[ll]] <- w.range
}
npost.range = (1:nrow(fit$beta[[1]]))[-c(1:burn_in)]
npost = length(npost.range)
index.use.final <- rep(T, nrow(dat.sim.all[[sceindex]][[simindex]]))
n.in.tol <- n.in.1 <- n.in.2 <-
  matrix(NA, ncol = length(t.unique), nrow = npost)
ll.w1 <- ll.w2 <- matrix(NA, ncol = sum(index.use.final), nrow = npost)
for(j in 1:npost){
  jindex <- npost.range[j]
  t.curr <- t.post.all[jindex, ]
  den.w1.0 <- get.den.mismodel.cons(t.unique = t.unique,
                                    alph.list = fit$alpha$error1,
                                    index = jindex,
                                    w.range.list = w.range.list)
  den.w2.0 <- get.den.mismodel.cons(t.unique = t.unique,
                                    alph.list = fit$alpha$error2,
                                    index = jindex,
                                    w.range.list = w.range.list)
  
  w1.j = dat.sim$gest_error_1[index.use.final]
  w2.j = dat.sim$gest_error_2[index.use.final]
  t.j = t.curr[index.use.final]
  for(k in 1:length(27:44)){
    ##### Observations with two error-prone GAs
    t.interest = t.unique[k]
    indexforj <- which(t.j == t.interest)
    if(length(indexforj) == 0){
      n.in.tol[j, k] <- length(indexforj)
      n.in.1[j, k] <- 0
      n.in.2[j, k] <- 0
    }else{
      
      n.in.tol[j, k] <- length(indexforj)
      y.kj.1 <- lapply(indexforj, function(x)
        get.yind(w = w1.j[x], t = t.j[x], ntol = ntol))
      y.kj.2 <- lapply(indexforj, function(x)
        get.yind(w = w2.j[x], t = t.j[x], ntol = ntol))
      
      t.index = t.interest - 26
      index = which(w.range.list[[t.index]] == t.interest)
      y.bin.1 <- as.numeric(sapply(y.kj.1, function(x) x[index] == 1))
      y.bin.2 <- as.numeric(sapply(y.kj.2, function(x) x[index] == 1))
      
      row.index = paste0("w",w.range.list[[t.index]])
      prob.w1.j <- den.w1.0[row.index,t.index]
      prob.w2.j <- den.w2.0[row.index,t.index]
      
      for(ii in 1:length(indexforj)){
        if(sum(y.kj.1[[ii]]) == 0){
          ll.w1[j, indexforj[ii]] <- NA
        }else{
          ll.w1[j, indexforj[ii]] <- dmultinom(x = y.kj.1[[ii]], size = 1,
                                               prob = prob.w1.j, log = T)
        }
        if(sum(y.kj.2[[ii]]) == 0){
          ll.w2[j, indexforj[ii]] <- NA
        }else{
          ll.w2[j, indexforj[ii]] <- dmultinom(x = y.kj.2[[ii]], size = 1,
                                               prob = prob.w2.j, log = T)
        }
      }
      n.in.1[j, k] <- sum(!is.na(ll.w1[j, indexforj]))
      n.in.2[j, k] <- sum(!is.na(ll.w2[j, indexforj]))
    }
    ### END ll computation for observations only having one error-prone GA
  } # END loop over different weeks
  if(j < 10){cat(j, ",")}
  if(j %% 1000 == 0){ cat(j, ",") }
}
ll.list <- list(n.in.tol = n.in.tol, n.in.1 = n.in.1, n.in.2 = n.in.2,
                ll.w1 = ll.w1, ll.w2 = ll.w2,
                nobs.include = sum(index.use.final))

## log pointwise predictive density (col: observations, row: post samples)
WAIC.mat <- matrix(NA, nrow = 1, ncol = 3)
WAIC.mat <- as.data.frame(WAIC.mat)
colnames(WAIC.mat) <- c("lppd", "pWAIC", "WAIC")
rownames(WAIC.mat) <- c("multi")

log.pd = ll.list$ll.w1 + ll.list$ll.w2
inter.lppd = apply(exp(log.pd), 2, mean, na.rm = T)
lppd = sum(log(inter.lppd[!is.na(inter.lppd)]))
pWAIC = sum(apply(log.pd, 2, var, na.rm = T), na.rm = T)
WAIC = -2*(lppd - pWAIC)
WAIC.mat[1, ] <- c(lppd, pWAIC, WAIC)

re.WAIC.list <- list(n.in.tol = colMeans(n.in.tol), 
                     n.in.1 = colMeans(n.in.1), n.in.2 = colMeans(n.in.2),
                     WAIC = WAIC.mat,
                     nobs.include = sum(index.use.final))
save(re.WAIC.list, file = 
       paste0("./interres/Res_WAIC_outmis_S", sceindex, "_sim", simindex, ".rda"))
############## END COMPUTE WAIC #######################


##### COMPUTE PMF of GAs (true, error_1, error_2, estimated) #####
##### (1) based on posterior samples of the true GA
##### (2) based on two hazard models Pr(Ti = t) 
#### (1) get PMF based on posterior samples of the true GA
prop.post <- list()
count = 1
for(k in (1:nrow(t.post.all))[-c(1:burn_in)]){
  prop.post[[count]] <-
    as.data.frame(prop.table(table(factor(t.post.all[k,], 
                                          levels = 27:44))))
  count = count + 1
  if(count%%1000 == 0){  cat(count, ",") }
}
prop.post.mat <- do.call(cbind, lapply(prop.post, function(x) x[,2]))
prop.est <- rowMeans(prop.post.mat)
prop.true <- as.data.frame(prop.table(table(factor(dat.sim$gest_true,
                                                   levels = 27:44))))
prop.error.1 <- as.data.frame(prop.table(table(factor(dat.sim$gest_error_1,
                                                      levels = 27:44))))
prop.error.2 <- as.data.frame(prop.table(table(factor(dat.sim$gest_error_2,
                                                      levels = 27:44))))
#### (2) get PMF based on two hazard models 
hazard.sum.list <- list()
for(kparm in 1:18){
  re.kk <- do.call(cbind, lapply(fit$hazard, 
                                 function(x) x[,kparm]))
  re.mse.kk <- lapply(1:length(hazard.sim.all[[sceindex]][,kparm]),
                      function(x) 
                        (mean(re.kk[x, ]) - 
                           hazard.sim.all[[sceindex]][x,kparm])^2)
  re.bias.kk <- lapply(1:length(hazard.sim.all[[sceindex]][,kparm]),
                       function(x) 
                         mean((re.kk[x, ] - 
                                 hazard.sim.all[[sceindex]][x,kparm])))
  
  re.sum.kk <- data.frame(rowMeans(re.kk), apply(re.kk, 1, sd),
                          apply(re.kk, 1, quantile, 0.025),
                          apply(re.kk, 1, quantile, 0.975),
                          do.call(c, re.mse.kk), 
                          do.call(c, re.bias.kk))
  colnames(re.sum.kk) <- c(paste0("t", kparm + 26, ".est"),
                           paste0("t", kparm + 26,".sd"),
                           paste0("t", kparm + 26,".lci"),
                           paste0("t", kparm + 26,".uci"),
                           paste0("t", kparm + 26,".mse"),
                           paste0("t", kparm + 26,".bias"))
  re.sum.kk <- data.frame(re.sum.kk, 
                          re.sum.kk[,paste0("t", kparm + 26,".lci")] <= 
                            hazard.sim.all[[sceindex]][,kparm] & 
                            re.sum.kk[,paste0("t", kparm + 26,".uci")] >= 
                            hazard.sim.all[[sceindex]][,kparm])
  colnames(re.sum.kk)[7] <- paste0("t", kparm + 26,".include")
  hazard.sum.list[[kparm]] <- re.sum.kk
}
names(hazard.sum.list) <- paste0("t", 27:44)
prop.mat.all <- data.frame("Week" = 27:44, 
                           "true" = prop.true$Freq, 
                           "error1" = prop.error.1$Freq, 
                           "error2" = prop.error.2$Freq,
                           "est" = prop.est,
                           "true_num" = colMeans(hazard.sim.all[[sceindex]]),
                           "est_num" = 
                             sapply(hazard.sum.list, 
                                    function(x) mean(x[,1])),
                           "bias_num" = 
                             sapply(hazard.sum.list, 
                                    function(x) mean(x[,6])),
                           "relbias_num" = 
                             sapply(1:18, function(x) 
                               mean(hazard.sum.list[[x]][,6]/
                                      hazard.sim.all[[sceindex]][,x])),
                           "mse_num" = sapply(hazard.sum.list, 
                                              function(x) mean(x[,5])),
                           "rmse_num" = sapply(hazard.sum.list, 
                                               function(x) mean(sqrt(x[,5]))))
prop.sum.list <- list("tab_sum" = prop.mat.all,
                      "include" = 
                        do.call(cbind, lapply(hazard.sum.list, 
                                              function(x) x[,7])))
save(prop.sum.list, file = 
       paste0("./interres/Res_PMF_S", sceindex, "_sim", simindex, ".rda"))
############## END COMPUTE PMF #######################


#### drop three elements in the list "fit" to save storage space 
fit <- fit[!(names(fit) %in%
               c("gestwk", "mis.parm", "hazard"))]
save(fit, file = 
       paste0("./interres/Res_S", sceindex, "_sim", simindex, ".rda"))


