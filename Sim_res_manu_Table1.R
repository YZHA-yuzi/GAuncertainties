##############################################################################
# This R script is used to generate Table 1 included manuscript 
# with intermediate results provided.
# NOTEs: (1) Please make sure the working directory is set to the folder where 
# the R script "FUNs.R" and "FUNs_rcpp.cpp" are located.
# Author: yuzi.zhang@emory.edu
###############################################################################

#### Simulation Set up for Table 1 in Manu ######
### 1. read in true SE, SP 
load("./interres/for_manu/Simulation_trueSESP.rda")
misparm.true.1 <- do.call(rbind, lapply(misparm.true.manu, 
                                        function(x) colMeans(x$error1)))
misparm.true.1 <- data.frame("Scenario" = paste0("S", 1:3),
                             misparm.true.1)

misparm.true.2 <- do.call(rbind, lapply(misparm.true.manu, 
                                        function(x) colMeans(x$error2)))
misparm.true.2 <- data.frame("Scenario" = paste0("S", 1:3),
                             misparm.true.2)
misparm.true.print <- 
  left_join(misparm.true.1[,-c(4:5)], misparm.true.2[,-c(4:5)],
            by = c("Scenario" = "Scenario"))
misparm.true.print[,-1] <- format(round(misparm.true.print[,-1], 2), 
                                  nsmall = 2)
colnames(misparm.true.print)[-1] <- c(paste0(c("SE", "SP"), "_error1"),
                                      paste0(c("SE", "SP"), "_error2"))



#### 2. read in estimated health effects across simulation scenarios
load("./interres/for_manu/Res_bets_formanu.rda")
tab.print.o3 <- tab.print.smk <- list()
for(sceindex in 1:3){
  tab.sum.beta <- do.call(rbind.data.frame, tab.bet.all[[sceindex]])
  tab.sum.beta$bias <- tab.sum.beta$Est - tab.sum.beta$True
  tab.sum.beta$rel.bias <-
    (tab.sum.beta$Est - tab.sum.beta$True)/tab.sum.beta$True
  tab.sum.beta$Model[tab.sum.beta$Model == "multi"] <- "proposed"
  tab.sum.beta$Scenario <- paste0("S", sceindex)
  tab.sum.beta.1 <- tab.sum.beta[,c("Scenario", "Model", "Cov", 
                                    "rel.bias", "avg.SE", "cover")]
  tab.sum.beta.1$Model <- factor(tab.sum.beta.1$Model,
                                 levels = c("proposed", "glm_true",
                                            "glm_error1", "glm_error2"))
  tab.sum.beta.1 <- tab.sum.beta.1[order(tab.sum.beta.1$Model), ]
  tab.print.o3[[sceindex]] <- subset(tab.sum.beta.1, Cov == "ptb_O3")
  tab.print.smk[[sceindex]] <- subset(tab.sum.beta.1, Cov == "ptb_Smoke")
}

tab.sum.beta.ozone <- do.call(rbind.data.frame, tab.print.o3)
tab.sum.beta.smk <- do.call(rbind.data.frame, tab.print.smk)

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


#### 3. read in estimated PTB rate across simulation scenarios
load("./interres/for_manu/Res_PTBrate_formanu.rda") # tab.ptb
tab.ptb.final <- data.frame(tab.ptb[, c("Scenario", "Outcome")],
                            "rel.bias" = format(round(tab.ptb$rel.bias, 3),
                                                nsmall = 3),
                            "cover" = format(round(tab.ptb$cover*100, 2),
                                             nsmall = 2))
tab.ptb.final$Outcome <- as.character(tab.ptb.final$Outcome)
colnames(tab.ptb.final) <- c("Scenario", "Model", "rel.bias", "cover")
tab.ptb.final$Model[tab.ptb.final$Model %in% "gest_true"] <- "glm_true"
tab.ptb.final$Model[tab.ptb.final$Model %in% "gest_error_1"] <- "glm_error1"
tab.ptb.final$Model[tab.ptb.final$Model %in% "gest_error_2"] <- "glm_error2"
tab.ptb.final$Model[tab.ptb.final$Model %in% "adjusted"] <- "proposed"
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

misparm.true.print$comb_name <- 
  paste0(paste0("(", misparm.true.print$SE_error1, ", ",
                misparm.true.print$SP_error1, ")"),
         "; ",
         paste0("(", misparm.true.print$SE_error2, ", ",
                misparm.true.print$SP_error2, ")"))

tab.simres.final <- left_join(tab.simres.final,
                              misparm.true.print[,c("Scenario", "comb_name")], 
                              by = c("Scenario" = "Scenario"))
colnames(tab.simres.final)[colnames(tab.simres.final) == "comb_name"] <-
  "(SE1,SP1);(SE2,SP2)"

tab.simres.final <- 
  tab.simres.final[,c("Scenario", "(SE1,SP1);(SE2,SP2)", "Outcome", 
                      "rel.bias_ozone", "SE_ozone", "CP_ozone",
                      "rel.bias_smk", "SE_smk", "CP_smk",
                      "rel.bias_PTB_rate", "CP_PTB_rate")]
## write simulation results
write.csv(tab.simres.final,
          file = paste0("./interres/for_manu/Table1.csv"), row.names = FALSE)







