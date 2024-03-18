##############################################################################
# This R script is used to summarize simulation results
# NOTEs: (1) Please make sure the working directory is set to the folder where 
# the R script "FUNs.R" and "FUNs_rcpp.cpp" are located.
# Author: yuzi.zhang@emory.edu
###############################################################################

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

## A function to compute posterior mean, sd, 2.5 and 97.5 percentiles
get.sum.multi <- function(post.mat, burn_in){
  beta.post.i <- post.mat[-c(1:burn_in), ]
  f.beta <- function(x){
    re = data.frame(est = colMeans(x), sd = apply(x, 2, sd),
                    lci = apply(x, 2, quantile, 0.025),
                    uci = apply(x, 2, quantile, 0.975))
  }
  re.ptb <- data.frame(cov = colnames(post.mat),
                       f.beta(beta.post.i))
  return(re.ptb)
}

### read in true values
load(paste0("./data/Data_hazard_sims.rda"))
## read in results obtained from conducting the 
## conventional preterm bith analysis using the regular logistic regression 
## based on error-prone GA estimates
# nobs = number of observations within each simulation dataset, 
nobs = 5000
ptb.true.ind <- matrix(NA, ncol = 3, nrow = nobs)
for(s in 1:3){
  ptb.true.ind[,s] <- rowSums(hazard.sim.all[[s]][,1:10])
}

## I. Summarize results for estimated health effects on preterm birth 
## I.(1) read in estimates from the proposed model and 
## the conventional analysis 
## !!!!!! here I set burn-in = 5 just for the purpose demonstration,
## it should be 20,000 !!!!! 
# burn_in = 20000
burn_in = 5
tab.bets.comp.all <- tab.est.multi.all <- list()
beta.glm.all <- list()
for(sceindex in 1:3){
  load(paste0("./interres/Res_glm_S", sceindex, ".rda")) # res.glm.all 
  nsim = length(res.glm.all$beta.int.all)
  # only extract effects of exposures on preterm birth 
  beta.all <- res.glm.all$beta.all 
  beta.glm.all[[sceindex]] <- beta.all
  beta.all.1 <- lapply(beta.all, function(x) x[,-1])
  beta1.ptb <- beta.all.1[paste0("sim", 1:nsim, "_ptb_o3")]
  beta2.ptb <- beta.all.1[paste0("sim", 1:nsim, "_ptb_smoke")]
  
  # read in true coefs used in simulation studies for the hazard model
  # for preterm birth
  inter.i <- readxl::read_xlsx(path = "./data/Tab_Supp_simssetup.xlsx",
                               sheet = "hazard_preterm")
  beta.true.ptb <- c(as.numeric(inter.i$value[inter.i$parm %in% "ozone"]),
                     as.numeric(inter.i$value[inter.i$parm %in% "smk"]))
  names(beta.true.ptb) <- c("ozone", "smk")

  cover.index <- matrix(NA, ncol = 4, nrow = nsim)
  colnames(cover.index) <- 
    c(paste0("glm_", c("true", "error1", "error2")), "proposed")
  cover.index.list <- list("ptb_O3" = cover.index,
                           "ptb_Smoke" = cover.index)
  tab.est.multi.s <- tab.bets.comp.s <- list()
  sim.range = 1:nsim
  # sim.range = 1:5
  for(i in sim.range){
    # read in results from the proposed model
    load(paste0("./interres/Res_S", sceindex, "_sim", i, ".rda"))
    
    # read in estimated betas for hazard model for preterm birth
    # obtained from proposed model and 
    # conventional analyses
    tab.est.multi <- tab.bets.comp <- list()
    l = 1
    tab.est.multi[[l]] <- get.sum.multi(post.mat = fit$beta[[l]], 
                                        burn_in = burn_in)
    sel.index = rownames(tab.est.multi[[l]]) %in% 
      c("trim3.cum", "factor(SMOKE)1")
    tab.bets.comp[[l]] <- data.frame(true = beta.true.ptb,
                                     rbind(beta1.ptb[[i]][,1],
                                           beta2.ptb[[i]][,1]),
                                     tab.est.multi[[l]][sel.index, 2])
    lci.l = cbind(rbind(c(beta1.ptb[[i]][,1] - 1.96*beta1.ptb[[i]][,2]),
                        c(beta2.ptb[[i]][,1] - 1.96*beta2.ptb[[i]][,2])),
                  tab.est.multi[[l]][sel.index, "lci"])
    uci.l = cbind(rbind(c(beta1.ptb[[i]][,1] + 1.96*beta1.ptb[[i]][,2]),
                        c(beta2.ptb[[i]][,1] + 1.96*beta2.ptb[[i]][,2])),
                  tab.est.multi[[l]][sel.index, "uci"])
    cover.index.list$ptb_O3[i, ] <- 
      c(as.numeric(sapply(1:4, function(x) 
        lci.l[1, x] <= beta.true.ptb[1] & beta.true.ptb[1] <= uci.l[1, x])))
    cover.index.list$ptb_Smoke[i, ] <- 
      c(as.numeric(sapply(1:4, function(x) 
        lci.l[2, x] <= beta.true.ptb[2] & beta.true.ptb[2] <= uci.l[2, x])))
    colnames(tab.bets.comp[[l]]) <- 
      c("True", paste0("glm_", c("true", "error1", "error2")), "proposed")
    rownames(tab.bets.comp[[l]]) <- paste0("ptb_", c("O3", "Smoke"))
    names(tab.bets.comp) <- names(tab.est.multi) <- c("ptb")
    tab.est.multi.s[[i]] <- tab.est.multi
    tab.bets.comp.s[[i]] <- tab.bets.comp
    cat(i, ",")
  } # end loop over different simulations
  tab.est.multi.all[[sceindex]] <- tab.est.multi.s
  tab.bets.comp.all[[sceindex]] <- tab.bets.comp.s
} # end loop over different simulation scenarios

## I. (2) get table summarize estimated health effects on PTB rate
tab.sum.beta.all <- list()
for(sceindex in 1:3){
  tab.sum.list <- list()
  cov.name.vec <- c("trim3.cum", "factor(SMOKE)1")
  cov.name.vec.1 <- c("o3", "smoke")
  cov.name.vec.2 <- c("O3", "Smoke")
  count = 1
  for(out.i in c("ptb")){
    for(cov.i in 1:2){
      tab0.comp <- do.call(rbind.data.frame, 
                           lapply(tab.bets.comp.all[[sceindex]], function(x) 
                             x[[out.i]][cov.i, ]))
      tab0 <- do.call(rbind.data.frame,
                      lapply(tab.est.multi.all[[sceindex]],
                             function(x) x[[out.i]][rownames(x[[out.i]]) %in%
                                                      cov.name.vec[cov.i], ]))
      index = paste0("sim", 1:nsim, "_", out.i, "_", cov.name.vec.1[cov.i])
      tab0.glm <- do.call(rbind, lapply(beta.glm.all[[sceindex]][index], 
                                        function(x) x[,3]))
      tab0.glm.lci <- do.call(rbind, lapply(beta.glm.all[[sceindex]][index], 
                                            function(x) x[,2] - 1.96*x[,3]))
      tab0.glm.uci <- do.call(rbind, lapply(beta.glm.all[[sceindex]][index], 
                                            function(x) x[,2] + 1.96*x[,3]))
      lci.mat <- cbind(tab0.glm.lci, tab0$lci)
      uci.mat <- cbind(tab0.glm.uci, tab0$uci)
      
      # read in true coefs used in simulation studies for the hazard model
      # for preterm birth
      inter.i <- readxl::read_xlsx(path = "./data/Tab_Supp_simssetup.xlsx",
                                   sheet = "hazard_preterm")
      beta.true.ptb <- c(as.numeric(inter.i$value[inter.i$parm %in% "ozone"]),
                         as.numeric(inter.i$value[inter.i$parm %in% "smk"]))
      names(beta.true.ptb) <- c("ozone", "smk")
      
      get.mse <- function(x, true.val){
        return(mean((x-true.val)^2))
      }
      tab.sum.i <- 
        data.frame("True" = as.numeric(beta.true.ptb[cov.i]), 
                   "Est" = colMeans(tab0.comp[,-1]),
                   "emp.SD" = apply(tab0.comp[,-1], 2, sd),
                   "avg.SE" = c(colMeans(tab0.glm), mean(tab0$sd)),
                   "cover" = 
                     colMeans(cover.index.list[[paste0(out.i, "_", 
                                                       cov.name.vec.2[cov.i])]], 
                              na.rm = T),
                   "avg.width" = sapply(1:ncol(lci.mat), 
                                        function(x) mean(uci.mat[,x] - lci.mat[,x])),
                   "MSE" = apply(tab0.comp[,-1], 2, 
                                 get.mse, true.val = 
                                   beta.true.ptb[cov.i]))
      tab.sum.i$bias <- tab.sum.i$Est - tab.sum.i$True
      tab.sum.i$rel.bias <-
        (tab.sum.i$Est - tab.sum.i$True)/tab.sum.i$True
      tab.sum.i <- data.frame("Scenario" = paste0("S", sceindex),
                              "Model" = rownames(tab.sum.i), 
                              "Cov" = paste0(out.i, "_", cov.name.vec.2[cov.i]),
                              tab.sum.i)
      tab.sum.list[[count]] <- tab.sum.i
      count = count + 1
    }
  }
  names(tab.sum.list) <- c(paste0("ptb_", c("O3", "Smoke")))
  tab.sum.beta.all[[sceindex]] <- tab.sum.list
} # end loop over different simulation scenarios
names(tab.sum.beta.all) <- paste0("S", 1:3)

tab.sum.beta.ozone <- do.call(rbind.data.frame,
                              lapply(tab.sum.beta.all, function(x) x$ptb_O3))
tab.sum.beta.smk <- do.call(rbind.data.frame,
                            lapply(tab.sum.beta.all, function(x) x$ptb_Smoke))

tab.print.final.1 <-
  data.frame(tab.sum.beta.ozone[, c("Scenario", "Model", "Cov")],
             "rel.bias" = format(round(tab.sum.beta.ozone$rel.bias, 3),
                                 nsmall = 3),
             "avg.SE" = format(round(tab.sum.beta.ozone$avg.SE, 3),
                               nsmall = 3),
             "cover" = format(round(tab.sum.beta.ozone$cover*100),
                              nsmall = 0))
tab.print.final.1$Model <- 
  factor(tab.print.final.1$Model,
         levels = c("proposed", paste0("glm_", c("true", "error1", "error2"))))

tab.print.final.2 <-
  data.frame(tab.sum.beta.smk[, c("Scenario", "Model", "Cov")],
             "rel.bias" = format(round(tab.sum.beta.smk$rel.bias, 3),
                                 nsmall = 3),
             "avg.SE" = format(round(tab.sum.beta.smk$avg.SE, 3),
                               nsmall = 3),
             "cover" = format(round(tab.sum.beta.smk$cover*100),
                              nsmall = 0))
tab.print.final.2$Model <- 
  factor(tab.print.final.2$Model,
         levels = c("proposed", paste0("glm_", c("true", "error1", "error2"))))

tab.bet.final <- dplyr::left_join(tab.print.final.1,
                                  tab.print.final.2, 
                                  by = c("Scenario" = "Scenario",
                                         "Model" = "Model"))
tab.bet.final <- tab.bet.final[order(tab.bet.final$Scenario, 
                                     tab.bet.final$Model), ]
rownames(tab.bet.final) <- NULL




## II. Summarize estimated overall preterm birth rate
# A function to summarize PTB rate from the conventional analysis
f.glm.ptb <- function(est.mat, true.val){
  re = data.frame("True" = true.val,
                  "Est" = mean(est.mat[,"Est"]),
                  "avg.bias" = mean(est.mat[,"Est"]) - true.val,
                  "rel.bias" = (mean(est.mat[,"Est"]) - true.val)/true.val,
                  "emp.SD" = sd(est.mat[,"Est"]),
                  "avg.SE" = mean(est.mat[,"SE"]),
                  "cover" = mean(est.mat[,"lci"] <= true.val &
                                   true.val <= est.mat[,"uci"]),
                  "mis.low" = mean(est.mat[,"uci"] < true.val),
                  "mis.high" = mean(est.mat[,"lci"] > true.val))
  return(re)
}
# A function to summarize PTB rate from the proposed model
f.prop.ptb <- function(est.mat){
  true.val = mean(est.mat$true_hazard)
  re = data.frame("True" = true.val,
                  "Est" = mean(est.mat[,"est_hazard"]),
                  "avg.bias" = mean(est.mat[,"est_hazard"]) - true.val,
                  "rel.bias" = (mean(est.mat[,"est_hazard"]) - true.val)/true.val,
                  "emp.SD" = sd(est.mat[,"est_hazard"]),
                  "avg.SE" = mean(est.mat[,"sd_hazard"]),
                  "cover" = mean(est.mat[,"lci_hazard"] <= true.val &
                                   true.val <= est.mat[,"uci_hazard"]),
                  "mis.low" = mean(est.mat[,"uci_hazard"] < true.val),
                  "mis.high" = mean(est.mat[,"lci_hazard"] > true.val))
  return(re)
}

### II. (1) read in estimated overall PTB rate from the conventional analysis
tab.sum.ind.all <- list()
out.vec <- c("gest_true", paste0("gest_error_", 1:2))
for(sceindex in 1:3){
  load(paste0("./interres/Res_glm_S", sceindex, ".rda")) # res.glm.all
  res.inde.glm <- 
    lapply(1:3, function(z) lapply(1:nobs, function(y)
      do.call(rbind.data.frame,
              lapply(res.glm.all$ptb, function(x) x[["ind"]][[z]][y, ]))))
  names(res.inde.glm) <- c("gest_true", paste0("gest_error_", 1:2))
  tab.sum.ind.all.i <- list()
  for(l in 1:3){
    tab.sum.ind.il <- lapply(1:nobs, function(x)
      data.frame("Scenario" = paste0("S", sceindex),
                 "Outcome" = out.vec[l],
                 "obs" = x,
                 f.glm.ptb(res.inde.glm[[l]][[x]], ptb.true.ind[x, sceindex])))
    tab.sum.ind.all.i[[l]] <- do.call(rbind.data.frame, tab.sum.ind.il)
    cat(paste0(sceindex, ", ", l), ":::")
  }
  tab.sum.ind.all[[sceindex]] <- tab.sum.ind.all.i
}
tab.print.ind <- list()
c = 1
for(sceindex in 1:3){
  for(l in 1:3){
    inter.il <- tab.sum.ind.all[[sceindex]][[l]]
    inter.il.1 <-
      as.data.frame(matrix(colMeans(
        inter.il[,c("True", "rel.bias", "cover")]), nrow = 1))
    colnames(inter.il.1) <- c("True", "rel.bias", "cover")
    tab.print.ind[[c]] <-
      cbind.data.frame(data.frame("Scenario" = inter.il$Scenario[1],
                                  "Outcome" = inter.il$Outcome[1]),
                       inter.il.1)
    c = c + 1
  }
}
tab.sum.ptb.glm <- do.call(rbind.data.frame, tab.print.ind)

### II. (2) summarize estimated PTB rate obtained using the proposed model
sum.ptb.prop.list <- list()
for(sceindex in 1:3){
  est.ptb.s.list <- list()
  for(simindex in 1:nsim){
    load(paste0("./interres/Res_PTBrate_S", sceindex, "_sim", 
                simindex, ".rda"))
    est.ptb.s.list[[simindex]] <-
      lapply(1:nobs, function(x)
        data.frame("Scenario" = paste0("S", sceindex),
                   "obs" = paste0("obs", x),
                   "sim" = paste0("sim", simindex),
                   ptb.list$ptb.mat.pw[x, ]))
  }
  est.ptb.s.mat <- lapply(1:nobs, function(y)
    do.call(rbind.data.frame,
            lapply(est.ptb.s.list, function(x) x[[y]])))
  sum.ptb.s.mat <- do.call(rbind.data.frame, lapply(est.ptb.s.mat, f.prop.ptb))
  sum.ptb.prop.list[[sceindex]] <- 
    sum.ptb.s.mat[,c("True", "rel.bias", "cover")]
}
tab.sum.ptb.prop <-
  data.frame("Scenario" = paste0("S", 1:3),
             "Outcome" = "proposed",
             do.call(rbind, lapply(sum.ptb.prop.list, colMeans)))

## II. (3) combine results from the proposed model and 
## conventional analyses
tab.ptb.comp <- rbind.data.frame(tab.sum.ptb.glm, tab.sum.ptb.prop)
tab.ptb.final <- data.frame(tab.ptb.comp[, c("Scenario", "Outcome")],
                            "rel.bias" = format(round(tab.ptb.comp$rel.bias, 3),
                                                nsmall = 3),
                            "cover" = format(round(tab.ptb.comp$cover*100),
                                             nsmall = 0))
colnames(tab.ptb.final) <- c("Scenario", "Model", "rel.bias", "cover")
tab.ptb.final$Model[tab.ptb.final$Model %in% "gest_true"] <- "glm_true"
tab.ptb.final$Model[tab.ptb.final$Model %in% "gest_error_1"] <- "glm_error1"
tab.ptb.final$Model[tab.ptb.final$Model %in% "gest_error_2"] <- "glm_error2"
tab.ptb.final$Model <- 
  factor(tab.ptb.final$Model,
         levels = c("proposed", paste0("glm_", c("true", "error1", "error2"))))
tab.ptb.final <- tab.ptb.final[order(tab.ptb.final$Scenario, 
                                     tab.ptb.final$Model), ]

## COMBINE results for estimated health effects and PTB rate
## from different methods across all three simulation scenarios
tab.simres.final.0 <- dplyr::left_join(tab.bet.final, tab.ptb.final,
                                       by = c("Scenario" = "Scenario",
                                              "Model" = "Model"))
tab.simres.final <- 
  data.frame("Scenario" = tab.simres.final.0$Scenario,
             "Method" = rep(c("Proposed", "Logistic", 
                              "Logistic", "Logistic"), 3),
             "Outcome" = rep(c("Two error-prone", "GA without errors", 
                               "1st error-prone GA", "2nd error-prone GA"), 3),
             tab.simres.final.0[,-c(1:3,7)])
colnames(tab.simres.final)[-c(1:3)] <- 
  c(paste0(c("rel.bias", "SE", "CP"), c("_ozone")),
    paste0(c("rel.bias", "SE", "CP"), c("_smk")),
    paste0(c("rel.bias", "CP"), c("_PTB_rate")))

## write simulation results
write.csv(tab.simres.final,
          file = paste0("./interres/Tab_simres.csv"), row.names = FALSE)











