############################################################################
#                      SELF-DEFINED FUNCTIONS 
############################################################################

# A function to compute conception date
get.concep <- function(db, ga_wk){
  return(db - ga_wk*7 + 1 + 11)
}

# A function to get binary indicators for I(Wi=k,Ti=j) 
# w \in { min{27, j-3}, j, max{44, j+3} }
get.yind <- function(w, t, ntol = 3){
  # ntol is the parameter d introduced in Section 3
  w.pos.val <- max(27, t-ntol):min(44, t+ntol) 
  npos.val = length(w.pos.val)
  re0 <- rep(0, npos.val)
  re0[which(w.pos.val == w)] <- 1
  return(re0)
}

# A function to compute exp(-|k-k'|/theta2) 
# (this is the component of the covariance matrix)
get.expcov0 <- function(psi.val.all, theta2, dist.name){
  # Input: pis.val.all: all possible combinations
  if(dist.name == "1st"){
    cov.psi.vec <- exp(-abs(psi.val.all[,1] - psi.val.all[,2])/theta2)
  }else if(dist.name == "2nd"){
    cov.psi.vec <- exp(-(psi.val.all[,1] - psi.val.all[,2])^2/theta2)
  }
  npost.0 <- length(unique(psi.val.all[,1]))
  cov.psi.star <- matrix(cov.psi.vec, ncol = npost.0, 
                         nrow = npost.0, byrow = F)
  return(cov.psi.star)
}
# A function to compute exponential covariance function
get.expcov <- function(psi.val.all, theta1, theta2, dist.name){
  cov.psi.star <- get.expcov0(psi.val.all = psi.val.all,
                              theta2 = theta2, dist.name = dist.name)
  cov.psi <- theta1*cov.psi.star
  return(cov.psi)
}


# Functions to compute conditional normal based on MVN 
get.premulti <- function(x, y){
  # x: covariance function
  # y: index of R.V. of interest
  x12 = matrix(x[y, -y], nrow = 1)
  x22 = x[-y, -y]
  return(x12 %*% solve(x22))
}
get.varcond <- function(x, y, premulti){
  # x: covariance function
  # y: index of R.V. of interest
  # premulti: sigma12 %in% sigma22^-1
  x11 = x[y, y]
  x12 = matrix(x[y, -y], nrow = 1)
  re = x11 - premulti %*% t(x12)
  return(as.numeric(re))
}
get.mucond <- function(y, premulti){
  # y: values of conditional R.V.
  # premulti: sigma12 %in% sigma22^-1
  y0 <- matrix(y, ncol = 1)
  return(as.numeric(premulti %*% y0))
}




# A function to generate error-prone gestational weeks
# based on Pr(w | t) computed from parameters estimated from the real data  
get.error.gest <- function(gest.true, den.mat){
  index <- sapply(gest.true, function(x) which(27:44 %in% x))
  gest.error.index <- t(do.call(cbind, lapply(index, function(x) 
    rmultinom(1, 1, as.numeric(den.mat[,x])))))
  gest.error <- apply(gest.error.index, 1, function(x) (27:44)[which(x == 1)])
  return(gest.error)
}

expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

# A function to create follow-up and event indicator
event = function(t, risk_period = c(27, 36)){ 
  p1 = risk_period[1]
  p2 = risk_period[2]
  nintervals = length(p1:p2)
  o = rep(0, nintervals)
  if (t < p1){ o[1:nintervals] = NA }
  if (t >= p1 & t < p2 ){ o[t-p1+1] = 1; o[(t-p1+2):nintervals]=NA}
  if (t == p2){ o[nintervals] = 1}
  return(o)
}

get.season <- function(x){
  index1 = x %in% c(12, 1, 2)
  index2 = x %in% 3:5
  index3 = x %in% 6:8
  index4 = x %in% 9:11
  y = rep(NA, length(x))
  y[index1] = "Winter"
  y[index2] = "Spring"
  y[index3] = "Summer"
  y[index4] = "Autumn"
  y <- factor(y, levels = c("Autumn", "Winter", "Spring", "Summer"))
  return(y)
}


# A function to compute Pr(Wi = k | Ti = j)
# based on alpha_kj, where the reference group is j 
# Note: rowsum is 1; K = 18 (27:44)
get.mismodel <- function(alph.vec, nweeks = 18){
  # alph.vec = (alpha_1j, ..., alpha_Kj)
  c.vec <- sapply(1:nweeks, function(x) sum(exp(alph.vec)[-x]))
  nu.vec <- exp(alph.vec)/c.vec
  re <- nu.vec/(1 + nu.vec)
  re[is.infinite(nu.vec)] <- 1
  return(re)
}

# A function to compute log-likelihood of outcome mis-measured model
# row is w_i and column is t
### constrainted version
### given Ti = j, Wi \in {max(27, j-3), j, min(44, j+3)}
get.den.mismodel.cons <- function(t.unique, alph.list, index, w.range.list){
  f.full <- function(x, vec.full){
    vec.re <- vec.full
    vec.re[names(vec.full) %in% names(x)] <- x
    return(vec.re)
  }
  nintervals.list = lapply(w.range.list, length)
  nweeks = length(t.unique)
  den.mismodel.list <- lapply(1:nweeks,
                              function(x) 
                                get.mismodel(alph.vec = alph.list[[x]][index, ],
                                             nweeks = nintervals.list[[x]]))
  vec.full <- rep(0, nweeks)
  names(vec.full) <- paste0("alph", t.unique)
  den.mismodel <- do.call(cbind, lapply(1:nweeks, function(x) 
    f.full(den.mismodel.list[[x]], vec.full)))
  rownames(den.mismodel) <- paste0("w", as.character(t.unique))
  colnames(den.mismodel) <- paste0("t", t.unique)
  return(den.mismodel)
}


get.Xdesign.tvarying <- function(dat, t.curr, dat.exp, exp.name, 
                                 season = T, concep.splines =T, df = 5,
                                 date.start, date.end, week.risk){
  
  # compute the first day of week 25
  exp.date.wk25 <- dat$BIRTHDATE - t.curr*7 + 1 + (25-1)*7
  # compute the first day of week 1
  exp.date.wk1 <- dat$BIRTHDATE - t.curr*7 + 1
  # identify those obs with birth date greater than 2003-12-31 
  exp.date <- exp.date.wk25
  dat$exp.date <- exp.date
  dat$concep.date <- exp.date.wk1
  month.vec <- lubridate::month(dat$concep.date)
  dat$concep.season <- get.season(month.vec)
  
  dat.comb <- left_join(dat, dat.exp, by = c("CENSUS" = "GEOID", 
                                             "exp.date" = "date"))
  week.risk.exp <- 27:36
  trim2.mat <- dat.comb[,c("ID",paste0("Wk", week.risk.exp))]
  colnames(trim2.mat) <- c("ID", as.character(week.risk.exp))
  
  if(identical(week.risk, 27:36)){
    if(exp.name == "trim3.week"){
      trim2.mat.long <- 
        gather(trim2.mat, Week, trim3.week, as.character(week.risk.exp))
    }else if(exp.name == "trim3.cum"){
      trim3.mat <- t(apply(trim2.mat[,-1], 1, dplyr::cummean))
      trim3.mat <- data.frame(ID = trim2.mat$ID, trim3.mat)
      colnames(trim3.mat)[-1] <- week.risk.exp
      trim2.mat.long <- 
        gather(trim3.mat, Week, trim3.cum, as.character(week.risk.exp))
    }
    D <- data.frame(ID = rep(unique(dat.comb$ID), each = length(week.risk)),
                    Week = rep(week.risk, nrow(dat.comb)))
    # Merge in time-invariant covariate by ID
    D <- left_join(D, dat.comb, by = c("ID" = "ID"))
    D <- D[order(D$Week), ]
    D$fake_event <- 1
  }else if(identical(week.risk, 37:44)){
    exp.name = "trim3.cum"
    D <- data.frame(ID = rep(unique(dat.comb$ID), each = length(week.risk)),
                    Week = rep(week.risk, nrow(dat.comb)))
    trim3.mat <- data.frame(ID = trim2.mat$ID, 
                            trim3.cum = rowMeans(trim2.mat[,-1]))
    trim2.mat.long <- left_join(D, trim3.mat, by = c("ID" = "ID"))
    trim2.mat.long <- trim2.mat.long[order(trim2.mat.long$Week), ]
    # Merge in time-invariant covariate by ID
    D <- left_join(D, dat.comb, by = c("ID" = "ID"))
    D <- D[order(D$Week), ]
    D$fake_event <- 1
  }
  concep.form <- paste0("ns(concep.date, df = ", df, 
                        ", Boundary.knots = c(", as.numeric(date.start),
                        ", ", as.numeric(date.end), ")" ,")")
  
  if(season & concep.splines){
    X.concep <- model.matrix(as.formula(
      paste0("fake_event ~ ",
             paste(c("factor(concep.season)", concep.form), 
                   collapse = " + "))), 
      data = D)
  }else if(season & (!concep.splines)){
    X.concep <- model.matrix(as.formula(
      paste0("fake_event ~ ",
             paste(c("factor(concep.season)"), 
                   collapse = " + "))), 
      data = D)
  }else if((!season) & concep.splines){
    X.concep <- model.matrix(as.formula(
      paste0("fake_event ~ ",
             paste(c(concep.form), 
                   collapse = " + "))), 
      data = D)
  }else if((!season) & (!concep.splines)){
    X.concep <- NULL
  }
  return(list(exp.tvarying.long = trim2.mat.long,
              concep.long = X.concep[,-1]))
  
}


# A function to obtain design matrix for time-invariant covariates
get.Xdesign.fix <- function(dat, week.risk, cov.name.vec){
  # INPUT
  # dat: datset of interest
  # week.risk: at-risk window, PTB 25:36, full term birth 37:44
  # cov.name.vec: covariates
  D.1 <- data.frame(ID = rep(unique(dat$ID), each = length(week.risk)),
                    Week = rep(as.character(week.risk), nrow(dat)))
  D.1$fake_event <- 1
  D.1 <- left_join(D.1, dat, by = c("ID" = "ID"))
  ## sort by week 
  D.1 <- D.1[order(D.1$Week), ]
  X.design.fix.full <- model.matrix(as.formula(paste0("fake_event ~ -1 + ",
                                                      paste(c(cov.name.vec), 
                                                            collapse = " + "))), 
                                    data = D.1)
  return(X.design.fix.full)
}

fit.initi <- function(dat, out.init,
                      exp.name, dat.exp, 
                      cov.name.vec = c("factor(Week)"),
                      season = T, concep.splines = T, df = 5,
                      date.start, date.end,
                      exposure = T, event_name = "PTB"){
  
  # event_name = "PTB" (pre-term birth) or "FULL" (full-term birth)
  date.start <- as.Date(date.start)
  date.end <- as.Date(date.end)
  
  concep.use <- (dat$BIRTHDATE - dat[[out.init]]*7 + 1) + 11
  dat <- subset(dat, concep.use <= date.end & date.start <= concep.use)
  
  # compute the first day of week 25
  exp.date.wk25 <- dat$BIRTHDATE - dat[[out.init]]*7 + 1 + (25-1)*7
  # compute the first day of week 1
  exp.date.wk1 <- dat$BIRTHDATE - dat[[out.init]]*7 + 1
  # identify those obs with birth date greater than the 
  # maximum date in exposure data
  exp.date <- exp.date.wk25
  dat$exp.date <- exp.date
  # control for the first day of gestation 
  dat$concep.date <- exp.date.wk1
  month.vec <- lubridate::month(dat$concep.date)
  dat$concep.season <- get.season(month.vec)
  
  if(exposure){
    exp.form = "trim3.cum"
    if(event_name == "PTB"){
      week.risk <- 27:36
      week.risk.exp <- 27:36
      dat.comb <- left_join(dat, dat.exp, by = c("CENSUS" = "GEOID", 
                                                 "exp.date" = "date"))
      trim2.mat <- dat.comb[,c("ID", paste0("Wk", week.risk.exp))]
      colnames(trim2.mat) <- c("ID", as.character(week.risk.exp))
      
      trim3.mat <- t(apply(trim2.mat[,-1], 1, dplyr::cummean))
      trim3.mat <- data.frame(ID = trim2.mat$ID, trim3.mat)
      colnames(trim3.mat)[-1] <- week.risk.exp
      keep.index = apply(!is.na(trim3.mat[,-1]), 1, all)
      
      trim3.mat <- trim3.mat[keep.index, ]
      trim3.mat.long <- gather(trim3.mat, Week, trim3.cum, 
                               as.character(week.risk.exp))
      trim3.mat.long$Week <- as.integer(trim3.mat.long$Week)
      dat.comb <- dat.comb[keep.index, ]
      
      D <- data.frame(ID = rep(unique(dat.comb$ID), each = length(week.risk)),
                      Week = rep(c(week.risk), nrow(dat.comb)))
      # Merge in time-invariante covariate by ID
      D <- left_join(D, dat.comb, by = c("ID" = "ID"))
      D <- left_join(D, trim3.mat.long, by = c("ID" = "ID", "Week" = "Week"))
      D$event = unlist(mapply(event, dat.comb[[out.init]], SIMPLIFY = FALSE, 
                              MoreArgs = list(risk_period = 
                                                c(week.risk[1], 
                                                  week.risk[length(week.risk)]))))
    }else if(event_name == "FULL"){
      week.risk <- 37:44
      week.risk.exp <- 27:36
      dat.comb <- left_join(dat, dat.exp, by = c("CENSUS" = "GEOID", 
                                                 "exp.date" = "date"))
      trim2.mat <- dat.comb[,c("ID", paste0("Wk", week.risk.exp))]
      colnames(trim2.mat) <- c("ID", as.character(week.risk.exp))
      trim2.mat <- data.frame(ID = trim2.mat$ID,
                              trim3.cum = rowMeans(trim2.mat[,-1]))
      D <- data.frame(ID = rep(unique(dat.comb$ID), each = length(week.risk)),
                      Week = rep(c(week.risk), nrow(dat.comb)))
      # Merge in time-invariante covariate by ID
      D <- left_join(D, dat.comb, by = c("ID" = "ID"))
      D <- left_join(D, trim2.mat, by = c("ID" = "ID"))
      D$event = unlist(mapply(event, dat.comb[[out.init]], SIMPLIFY = FALSE, 
                              MoreArgs = list(risk_period = 
                                                c(week.risk[1], 
                                                  week.risk[length(week.risk)]))))
    }
  }else{
    exp.form = NULL
    if(event_name == "PTB"){
      week.risk <- 27:36
    }else if(event_name == "FULL"){
      week.risk <- 37:44
    }
    dat.comb <- dat
    D <- data.frame(ID = rep(unique(dat.comb$ID), each = length(week.risk)),
                    Week = rep(week.risk, nrow(dat.comb)))
    # Merge in time-invariant covariate by ID
    D <- left_join(D, dat.comb, by = c("ID" = "ID"))
    D$event = unlist(mapply(event, dat.comb[[out.init]], SIMPLIFY = FALSE, 
                            MoreArgs = list(risk_period = 
                                              c(week.risk[1], 
                                                week.risk[length(week.risk)]))))
  }
  concep.form <- paste0("ns(concep.date, df = ", df, 
                        ", Boundary.knots = c(", as.numeric(date.start),
                        ", ", as.numeric(date.end), ")" ,")")
  if(season & concep.splines){
    model.form <- paste0("event ~ -1 + ",
                         paste(c(exp.form, cov.name.vec, 
                                 "factor(concep.season)", concep.form), 
                               collapse = " + "))
  }else if(season & (!concep.splines)){
    model.form <- paste0("event ~ -1 + ",
                         paste(c(exp.form, cov.name.vec, 
                                 "factor(concep.season)"), 
                               collapse = " + "))
  }else if((!season) & concep.splines){
    model.form <- paste0("event ~ -1 + ",
                         paste(c(exp.form, cov.name.vec, concep.form), 
                               collapse = " + "))
  }else if((!season) & (!concep.splines)){
    model.form <- paste0("event ~ -1 + ", paste(c(exp.form, cov.name.vec), 
                                                collapse = " + "))
  }
  res.initi <- glm(as.formula(model.form), data = D, family = binomial)
  return(res.initi)
}


### A function to fit the proposed model ###
fit.sim <- function(dat, dat.exp,
                    niter = 25000, burn_in = 20000, seed = 1234,
                    exp.name = "trim3.cum",
                    get.Xdesign,
                    cov.name.vec = list(c("factor(Week)"),
                                        c("factor(Week)")),
                    date.start, date.end,
                    season = c(T, T), concep.splines = c(T, T), df = c(3, 3),
                    exposure = c(T, T),
                    dist.name = "1st", ntol = 3){
  ## INPUTs:
  ## dat: data frame
  ## w1.vec: "Clinical" estimate of gestational age in weeks
  ## w2.vec: Last menstrual period estimate of gestational age in weeks
  ## niter (# iterations)
  ## dist.name = "abs" = "1st"/"Euclidean" = "2nd"
  set.seed(seed)
  t.unique = 27:44
  nweeks = length(t.unique); nobs = nrow(dat)
  index11 <- (!is.na(dat$gest_error_1)) & (!is.na(dat$gest_error_2))
  index11.pos <- which(index11)
  index10 <- (!is.na(dat$gest_error_1)) & (is.na(dat$gest_error_2))
  index10.pos <- which(index10)
  index01 <- (is.na(dat$gest_error_1)) & (!is.na(dat$gest_error_2))
  index01.pos <- which(index01)
  
  index1dot <- !is.na(dat$gest_error_1)
  index2dot <- !is.na(dat$gest_error_2)
  index.all <- (index11 | index10 | index01)
  
  out.init = "gest_error_2"
  coef.init.list <- list()
  event_vec <- c("PTB", "FULL")
  for(ll in 1:2){
    # fit glm to get initial values for betas
    res.init <- fit.initi(dat = dat, out.init = out.init,
                          exp.name = exp.name, 
                          dat.exp = dat.exp,
                          cov.name.vec = cov.name.vec[[ll]],
                          season = season[ll], 
                          concep.splines = concep.splines[ll],
                          df = df[ll], 
                          date.start = date.start,
                          date.end = date.end,
                          exposure = exposure[ll],
                          event_name = event_vec[ll])
    coef.init.list[[ll]] <- coef(res.init)
  }
  nbetas.vec = sapply(coef.init.list, length)
  if(sum(names(coef.init.list[[2]]) %in% paste0("factor(Week)", 37:44)) < 8){
    nbetas.vec[2] <- nbetas.vec[2] + (8 - 
                                        sum(names(coef.init.list[[2]]) %in% paste0("factor(Week)", 37:44)))
  }
  
  ## priors ##
  ## beta (model coefficients)
  nbetas.ptb <- nbetas.vec[1]
  c.beta.1 = matrix(rep(0, nbetas.ptb), ncol = 1)
  C.beta.pre.1 = diag(rep(1/100, nbetas.ptb))
  nbetas.birth <- nbetas.vec[2]
  c.beta.2 = matrix(rep(0, nbetas.birth), ncol = 1)
  C.beta.pre.2 = diag(rep(1/100, nbetas.birth))
  
  C.beta.pre <- list(C.beta.pre.1, C.beta.pre.2)
  c.beta <- list(c.beta.1, c.beta.2)
  
  ## IG for theta1 and theta2 ##
  a1 = b1 = 0.01; a2 = 1; b2 = 1
  ## tuning parameter for theta2_1 and theta2_2
  theta2.tun = c("theta2_1" = 0.5, "theta2_2" = 0.5)
  
  ## initialize data frames to store results 
  # true gestational age 
  t.mat <- matrix(NA, ncol = nobs, nrow = niter + 1)
  # beta in outcome-prediction model (preterm birth and full-term birth)
  beta.list <- list()
  for(ll in 1:2){
    if(ll == 2){
      cond.ll = names(coef.init.list[[ll]]) %in% 
        paste0("factor(Week)", 37:44)
      if(sum(cond.ll) < 8){
        beta.list[[ll]] <- 
          as.data.frame(matrix(NA, ncol = nbetas.vec[ll], nrow = niter + 1))
        colnames(beta.list[[ll]]) <- 
          c("trim3.cum", paste0("factor(Week)", 37:44), "factor(SMOKE)1")
        beta.list[[ll]][1, names(coef.init.list[[ll]])] <- 
          as.numeric(coef.init.list[[ll]])
        not.ll <- ! c("trim3.cum", paste0("factor(Week)", 37:44), 
                      "factor(SMOKE)1") %in% names(coef.init.list[[ll]])
        beta.list[[ll]][1, not.ll] <- 
          max(as.numeric(coef.init.list[[ll]]))
      }else{
        beta.list[[ll]] <- 
          as.data.frame(matrix(NA, ncol = nbetas.vec[ll], nrow = niter + 1))
        colnames(beta.list[[ll]]) <- names(coef.init.list[[ll]])
        beta.list[[ll]][1, ] <- as.numeric(coef.init.list[[ll]])
      }
    }else if(ll == 1){
      beta.list[[ll]] <- 
        as.data.frame(matrix(NA, ncol = nbetas.vec[ll], nrow = niter + 1))
      colnames(beta.list[[ll]]) <- names(coef.init.list[[ll]])
      beta.list[[ll]][1, ] <- as.numeric(coef.init.list[[ll]])
    }
  }
  
  ll.health <- list("ptb" = rep(NA, niter + 1),
                    "birth" = rep(NA, niter + 1))
  
  # parameters in outcome mis-classified model
  parm.mat <- matrix(NA, ncol = 4, nrow = niter + 1)
  colnames(parm.mat) <- c(paste0("theta", 1:2, "_1"),
                          paste0("theta", 1:2, "_2"))
  parm.mat[1, ] <- rep(c(1, 0.5), 2)
  
  accept.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = niter + 1))
  colnames(accept.mat) <- paste0("theta2_", 1:2)
  
  w.curr <- dat$gest_error_1
  w.curr[index01.pos] <- dat$gest_error_2[index01.pos]
  t.mat[1, ] <- w.curr
  t.index <- sample(index11.pos, round(length(index11.pos)/2))
  t.mat[1, t.index] <- dat$gest_error_2[t.index]
  
  # A function to assign names to alpha_kj matrices
  w.range.list <- alph.list <- psi.val.list <- list()
  ref.index.vec <- c()
  for(ll in 1:nweeks){
    t.curr <- t.unique[ll]
    w.range <- max(27, t.curr-ntol):min(44, t.curr+ntol) 
    w.range.list[[ll]] <- w.range
    inter <- matrix(NA, ncol = length(w.range), nrow = niter + 1)
    colnames(inter) <- paste0("alph", w.range)
    # since always pick w = j as the reference group
    inter[,colnames(inter) %in% paste0("alph", ll+26)] <- 0
    alph.list[[ll]] <- inter
    psi.val.list[[ll]] <- expand.grid(w.range[w.range!=t.curr], 
                                      w.range[w.range!=t.curr])
    ref.index.vec <- c(ref.index.vec, which(w.range == t.curr))
  }
  names(alph.list) <- paste0("T", 1:nweeks)
  alph.all.list <- list("error1" = alph.list,
                        "error2" = alph.list)
  
  ## compute the number of alpha_kj, this value will be used to update
  ## theta1
  ig.shape <- sum(sapply(w.range.list, length) - 1)
  
  ### get initial values of alpha_kj using multinomial logistic regression
  ### use function multinom from package nnet 
  for(l in 1:2){
    freq.mat.int <- table(factor(t.mat[1, ], levels = t.unique), 
                          factor(dat[[paste0("gest_error_", l)]], 
                                 levels = t.unique))
    GA.names.vec <- paste0("Wk", t.unique)
    for(j in 1:nweeks){
      w.range.j <- w.range.list[[j]]
      dat.sub.0 <- data.frame(GA = paste0("Wk", w.range.j), 
                              freq = freq.mat.int[j, w.range.j-26])
      dat.sub.0 <- do.call(c, lapply(1:length(w.range.j), function(x) 
        rep(dat.sub.0$GA[x], dat.sub.0$freq[x])))
      dat.sub.0 <- data.frame(GA.c = dat.sub.0)
      dat.sub.0$GA.f <- factor(dat.sub.0$GA.c,
                               levels = paste0("Wk", w.range.j))
      ## !! using the j-th group as the reference group  
      dat.sub.0$GA.f <- relevel(dat.sub.0$GA.f, ref = GA.names.vec[j])
      if(length(unique(dat.sub.0$GA.f)) > 1){
        fit0 <- nnet::multinom(GA.f ~ 1, data = dat.sub.0)
        alph.int <- coef(fit0)
        index.0 <- which(colnames(alph.all.list[[l]][[j]]) %in% 
                           paste0("alph", j+26))
        if(is.null(rownames(alph.int))){
          index.1 <- which(colnames(alph.all.list[[l]][[j]]) %in% 
                             unique(dat.sub.0$GA.c[dat.sub.0$GA.c != GA.names.vec[j]]))
          alph.all.list[[l]][[j]][1, index.1] <- as.numeric(alph.int)
          alph.all.list[[l]][[j]][1, -c(index.0, index.1)] <- -12
        }else{
          index.1 <- which(colnames(alph.all.list[[l]][[j]]) %in% 
                             paste0("alph", substr(rownames(alph.int), 3, 4)))
          alph.all.list[[l]][[j]][1, index.1] <- as.numeric(alph.int)
          alph.all.list[[l]][[j]][1, -c(index.0, index.1)] <- -12
        }
      }else{
        index.0 <- which(colnames(alph.all.list[[l]][[j]]) %in% 
                           paste0("alph", j+26))
        alph.all.list[[l]][[j]][1, -index.0] <- -12
      }
    } # end loop over ti
  } # end loop over two mis-measured outcome
  
  week.risk.ptb <- 27:36; week.risk.birth <- 37:44
  X.fix.full.ptb <- get.Xdesign.fix(dat = dat, week.risk = week.risk.ptb, 
                                    cov.name.vec = cov.name.vec[[1]])
  X.fix.full.birth <- get.Xdesign.fix(dat = dat, week.risk = week.risk.birth, 
                                      cov.name.vec = cov.name.vec[[2]])
  X.fix.full.list <- list(X.fix.full.ptb, X.fix.full.birth)
  
  nobs.vec <- rep(NA, niter)
  hazard.list <- list()
  mis.parm.list <- list("error1" = list(), "error2" = list())
  #### BEGIN MCMC procedure
  # start = proc.time()[3]
  for(i in 1:niter){
    
    t.curr <- t.mat[i, ]
    concep.t.i <- dat$BIRTHDATE - t.curr*7 + 1 + 11
    include.index <- concep.t.i <= as.Date(date.end) & 
      as.Date(date.start) <= concep.t.i
    
    index1dot.i <- include.index & index1dot
    index2dot.i <- include.index & index2dot
    indexdot.i <- list(index1dot.i, index2dot.i)
    
    w1.vec.i <- dat$gest_error_1[index1dot.i]
    w2.vec.i <- dat$gest_error_2[index2dot.i]
    w.error.vec.i <- list(w1.vec.i, w2.vec.i)
    
    index.all.i <- (index11 | index10 | index01) & include.index
    index11.i <- index11 & include.index
    index10.i <- index10 & include.index
    index01.i <- index01 & include.index
    nobs.vec[i] <- sum(index.all.i)
    
    #### I. update parameters in outcome misclassified model
    # 1. update alpha_kj for w1i and w2i, conjugate (Poly-Gamma Approach)
    # alpha_kj has prior MVN, 
    # Cov(alpha_kj, alpha_k'j) = theta1 exp(-|k-k'|/theta2) 
    # start = proc.time()[3]
    cov.psi.star.curr <- list()
    for(l in 1:2){
      w.l = w.error.vec.i[[l]]; nobs.l = length(w.l)
      t.l = t.mat[i, indexdot.i[[l]]]
      ## a list of length nobs, which contains I(Wi = k, Ti = j) 
      # y.kj.i <- lapply(1:nobs.l, function(x) 
      #   get.yind(w = w.l[x], t = t.l[x], nweeks = nweeks, ntol = ntol))
      
      ## for given k and j, update alpha_kj for wl (l = 1, 2)
      cov.psi.star.curr[[l]] <- list()
      for(j in 1:nweeks){
        ## pre-compute conditional normal
        psi.val.j <- psi.val.list[[j]]
        cov.psi.star <- get.expcov0(psi.val.all = psi.val.j, 
                                    theta2 = parm.mat[i, paste0("theta2_", l)], 
                                    dist.name = dist.name)
        cov.psi.star.curr[[l]][[j]] <- cov.psi.star
        cov.psi <- parm.mat[i, paste0("theta1_", l)]*cov.psi.star
        npos.j <- length(w.range.list[[j]]) - 1
        premulti.val <- lapply(1:npos.j, function(z) 
          get.premulti(x = cov.psi, y = z))
        val.cond <- sapply(1:npos.j, function(z) 
          get.varcond(x = cov.psi, y = z, premulti = premulti.val[[z]]))
        names(premulti.val) <- names(val.cond) <- 
          paste0("alph", w.range.list[[j]][w.range.list[[j]] != t.unique[j]])
        
        ## get indices for the reference group and group that will be updated
        ref.index <- which(w.range.list[[j]] == t.unique[j])
        windex.update <- which(w.range.list[[j]] != t.unique[j]) 
        alph.j.curr <- alph.all.list[[l]][[j]][i, ]
        
        indexforj <- which(t.l == t.unique[j])
        
        if(length(indexforj) == 0){
          alph.all.list[[l]][[j]][i+1, ] <- alph.j.curr
        }else{
          y.kj.i <- lapply(indexforj, function(x) 
            get.yind(w = w.l[x], t = t.l[x], ntol = ntol))
          nobs.kj <- length(y.kj.i)
          for(k.i in windex.update){
            c.kj.curr <- sum(exp(alph.j.curr)[-k.i])
            y.kj.vec <- sapply(y.kj.i, function(x) x[k.i])
            eta.kj.curr <- alph.j.curr[k.i] - log(c.kj.curr)
            omega.kj.vec <- rpg(nobs.kj, 1, eta.kj.curr)
            z.kj.vec <- (y.kj.vec - 1/2)/omega.kj.vec
            z.kj.foralph.vec <- z.kj.vec + log(c.kj.curr)
            ## posterior mean and variance of alpha_kj ##
            k.ii <- paste0("alph", w.range.list[[j]][k.i])
            M <- 1/(sum(omega.kj.vec) + 1/val.cond[k.ii])
            mu.cond.curr <- get.mucond(y = alph.j.curr[-c(k.i, ref.index)], 
                                       premulti = premulti.val[[k.ii]])
            m <- M*(sum(omega.kj.vec*z.kj.foralph.vec) + 
                      (1/val.cond[k.ii])*mu.cond.curr)
            alph.j.curr[k.i] <- rnorm(1, m, sqrt(M))
          } # end loop over k
          alph.all.list[[l]][[j]][i+1, ] <- alph.j.curr
        }
      } # end loop over j
    } # end loop over l
    # proc.time()[3] - start
    
    ##### Cov(psi_kj, psi_k'j) = theta_1 exp(-|k-k'|/theta_2), 
    ##### k = 28 - 44 
    ## 3. Update theta1 (reference groups don't be included)
    for(l in 1:2){
      cov.psi.star.inv <- lapply(cov.psi.star.curr[[l]], solve)
      ig.rate.post <- sum(sapply(1:nweeks, 
                                 function(x) 
                                   t(alph.all.list[[l]][[x]][i+1, -ref.index.vec[x]])%*%
                                   cov.psi.star.inv[[x]]%*%
                                   alph.all.list[[l]][[x]][i+1, -ref.index.vec[x]]))
      parm.mat[i+1, paste0("theta1_", l)] <- 1/rgamma(1, ig.shape/2 + a1, 
                                                      b1 + ig.rate.post/2)
    }
    
    ## 4. Update theta2
    # sample theta2 from lognormal(log(theta2_{i}), sdlog = theta2.tun) #
    for(l in 1:2){
      theta2.c <- parm.mat[i, paste0("theta2_", l)]
      theta1.c <- parm.mat[i+1, paste0("theta1_", l)]
      theta2.prop <- rlnorm(1, log(theta2.c), theta2.tun[l])
      
      Sigma.prop <- lapply(1:nweeks, function(x) 
        get.expcov(psi.val.all = psi.val.list[[x]], theta1 = theta1.c,
                   theta2 = theta2.prop, dist.name = dist.name))
      Sigma.c <- lapply(1:nweeks, function(x) 
        get.expcov(psi.val.all = psi.val.list[[x]], theta1 = theta1.c,
                   theta2 = theta2.c, dist.name = dist.name))
      
      ll.prop <- sum(sapply(1:nweeks, function(y)
        mvtnorm::dmvnorm(alph.all.list[[l]][[y]][i+1, -ref.index.vec[y]], 
                         sigma = Sigma.prop[[y]], log = T))) + 
        dgamma(theta2.prop, a2, b2, log = T) + log(theta2.prop)
      ll.c <- sum(sapply(1:nweeks, function(y)
        mvtnorm::dmvnorm(alph.all.list[[l]][[y]][i+1, -ref.index.vec[y]], 
                         sigma = Sigma.c[[y]], log = T))) + 
        dgamma(theta2.c, a2, b2, log = T) + log(theta2.c)
      
      ratio <- min(1, exp(ll.prop - ll.c))
      if(ratio > runif(1)){
        parm.mat[i+1, paste0("theta2_", l)]  <- theta2.prop
        accept.mat[i+1, paste0("theta2_", l)] <- 1
      }else{
        parm.mat[i+1, paste0("theta2_", l)]  <- theta2.c
        accept.mat[i+1, paste0("theta2_", l)] <- 0
      }
    }
    
    # 2. update true gestational weeks t.mat
    ### faster procedure 
    # get design matrix for time-varying covariates
    pred.val.list <- list()
    week.risk.list <- list(week.risk.ptb, week.risk.birth)
    for(ll in 1:2){
      re.X.tvarying <- get.Xdesign(dat = dat[index.all, ],
                                   dat.exp = dat.exp,
                                   t.curr = t.mat[i, index.all],
                                   exp.name = exp.name,
                                   season = season[ll],
                                   concep.splines = concep.splines[ll],
                                   df = df[ll],
                                   date.start = date.start,
                                   date.end = date.end,
                                   week.risk = week.risk.list[[ll]])
      ## design matrix of cumulative average of O3
      X.tvarying.full <- re.X.tvarying$exp.tvarying.long
      ## design matrix of conception date
      X.concep <- re.X.tvarying$concep.long
      # compute predictive values separately 
      if(exposure[ll]){
        pred.comp1 <- as.numeric(X.tvarying.full[,3]*beta.list[[ll]][i, 1])
      }else{pred.comp1 <- 0}
      index.sel <- colnames(beta.list[[ll]])%in%colnames(X.fix.full.list[[ll]])
      pred.comp2 <- as.numeric(
        X.fix.full.list[[ll]]%*%
          matrix(as.numeric(beta.list[[ll]][i, index.sel]), ncol = 1))
      if(concep.splines[ll]|season[ll]){
        index.sel <- colnames(beta.list[[ll]]) %in% colnames(X.concep)
        pred.comp3 <- as.numeric(
          X.concep %*%
            matrix(as.numeric(beta.list[[ll]][i, index.sel]), ncol = 1))
        pred.val.list[[ll]] <- pred.comp1 + pred.comp2 + pred.comp3
      }else{
        pred.val.list[[ll]] <- pred.comp1 + pred.comp2
      }
    } # end up loop over different health models
    ### Compute Pr(ti = t) for each observation t \in 27:36
    h0.mat.1 <- matrix(expit(pred.val.list[[1]]), ncol = length(week.risk.ptb), 
                       byrow = FALSE)
    colnames(h0.mat.1) <- paste0("wk", week.risk.ptb)
    hazard.curr.mat.1 <- comp_hazard(h0.mat.1)
    prob.cond.birth <- 1 - rowSums(hazard.curr.mat.1)
    ### Compute Pr(ti = t) for each observation t \in 37:44
    h0.mat.2 <- matrix(expit(pred.val.list[[2]]), ncol = length(week.risk.birth), 
                       byrow = FALSE)
    colnames(h0.mat.2) <- paste0("wk", week.risk.birth)
    hazard.curr.mat.2 <- comp_hazard(h0.mat.2)
    hazard.curr.mat.2 <- 
      do.call(rbind, 
              lapply(1:nrow(hazard.curr.mat.1), 
                     function(x) prob.cond.birth[x]*hazard.curr.mat.2[x, ]))
    hazard.curr.mat <- cbind(hazard.curr.mat.1, hazard.curr.mat.2)
    colnames(hazard.curr.mat) <- paste0("wk", c(week.risk.ptb, week.risk.birth))
    
    #############################################################
    ## keep tracking hazard 
    ## compute SE, SP, PPV, NPV
    if(i > burn_in){
      hazard.list[[i - burn_in]] <- hazard.curr.mat
    }
    den.w1.0 <- get.den.mismodel.cons(t.unique = t.unique,
                                      alph.list = alph.all.list[[1]], 
                                      index = i,
                                      w.range.list = w.range.list)
    den.w2.0 <- get.den.mismodel.cons(t.unique = t.unique,
                                      alph.list = alph.all.list[[2]], 
                                      index = i,
                                      w.range.list = w.range.list)
    # start = proc.time()[3]
    parm.ind.w1 <- 
      do.call(rbind, lapply(1:nrow(hazard.curr.mat), function(z) 
        get.mis.parm(den.w = den.w1.0, hazard.vec = hazard.curr.mat[z, ])))
    parm.ind.w1 <- as.data.frame(parm.ind.w1)
    parm.ind.w2 <- 
      do.call(rbind, lapply(1:nrow(hazard.curr.mat), function(z) 
        get.mis.parm(den.w = den.w2.0, hazard.vec = hazard.curr.mat[z, ])))
    parm.ind.w2 <- as.data.frame(parm.ind.w2)
    colnames(parm.ind.w1) <- colnames(parm.ind.w2) <- c("SE", "SP", "PPV", "NPV")
    # proc.time()[3] - start
    if(i > burn_in){
      mis.parm.list[["error1"]][[i - burn_in]] <- parm.ind.w1
      mis.parm.list[["error2"]][[i - burn_in]] <- parm.ind.w2
    }
    #############################################################
    
    # row sums of hazard.curr.mat are 1
    den.w1 <- data.frame(w = t.unique, 
                         get.den.mismodel.cons(t.unique = t.unique,
                                               alph.list = alph.all.list[[1]], 
                                               index = i+1,
                                               w.range.list = w.range.list))
    den.w2 <- data.frame(w = t.unique, 
                         get.den.mismodel.cons(t.unique = t.unique,
                                               alph.list = alph.all.list[[2]], 
                                               index = i+1,
                                               w.range.list = w.range.list))
    l.comp1.11 = left_join(data.frame(w = dat$gest_error_1[index11]), 
                           den.w1, by = c("w" = "w"))
    l.comp2.11 = left_join(data.frame(w = dat$gest_error_2[index11]), 
                           den.w2, by = c("w" = "w"))
    l.comp1.10 = left_join(data.frame(w = dat$gest_error_1[index10]), 
                           den.w1, by = c("w" = "w"))
    l.comp2.01 = left_join(data.frame(w = dat$gest_error_2[index01]), 
                           den.w2, by = c("w" = "w"))
    nobs.i <- nrow(l.comp1.11) + nrow(l.comp1.10) + nrow(l.comp2.01)
    ll.comp.all <- matrix(NA, nrow = nobs.i, ncol = nweeks)
    ll.comp.all[index11, ] <- as.matrix(l.comp1.11[,-1]*l.comp2.11[,-1], 
                                        ncol = nweeks)
    ll.comp.all[index10, ] <- as.matrix(l.comp1.10[,-1], ncol = nweeks)
    ll.comp.all[index01, ] <- as.matrix(l.comp2.01[,-1], ncol = nweeks)
    
    post.prob = ll.comp.all*hazard.curr.mat
    ## get indices that the mis-classification model implies that
    ## there is no such true GA can produce such two mis-measured GAs 
    ## (e.g., w1 = 30 and w2 = 37)
    index.zero = apply(ll.comp.all, 1, function(x) all(x == 0))
    post.prob[index.zero, ] <- hazard.curr.mat[index.zero, ]
    
    post.prob = t(apply(post.prob, 1, function(x) x/sum(x)))
    if(sum(is.na(post.prob)) > 0){
      rowindex = apply(post.prob, 1, function(x) all(is.na(x)))
      post.prob[rowindex, ] <- 1/ncol(post.prob)
    }
    t.mat[i+1, ] <- apply(post.prob, 1, 
                          function(y) sample(x = t.unique, size = 1, prob = y))
    
    # prepare data for fitting logistic regressions using the PG method
    ### update index using updated GA 
    t.curr <- t.mat[i+1, ]
    concep.t.i <- dat$BIRTHDATE - t.curr*7 + 1 + 11
    include.index <- concep.t.i <= as.Date(date.end) & 
      as.Date(date.start) <= concep.t.i
    
    index1dot.i <- include.index & index1dot
    index2dot.i <- include.index & index2dot
    indexdot.i <- list(index1dot.i, index2dot.i)
    index.all.i <- (index11 | index10 | index01) & include.index
    
    ### Update parameters in health modes
    ### discrete survival model for preterm birth
    ### discrete survival model for full-term birth
    for(ll in 1:2){
      npositive.mat <- 
        data.frame(ID = rep(unique(dat$ID), each = length(week.risk.list[[ll]])),
                   Week = rep(week.risk.list[[ll]], nrow(dat)))
      npositive.mat$event = 
        unlist(mapply(event, t.mat[i+1, ], SIMPLIFY = FALSE,
                      MoreArgs = 
                        list(risk_period = c(week.risk.list[[ll]][1], 
                                             week.risk.list[[ll]][length(week.risk.list[[ll]])]))))
      npositive.mat <- npositive.mat[order(npositive.mat$Week), ]
      index.long.all.i <- npositive.mat$ID %in% dat$ID[index.all.i]
      npositive.mat <- npositive.mat[index.long.all.i, ]
      npositive.vec <- npositive.mat$event
      keep.index.long = !is.na(npositive.vec)
      npositive.vec = npositive.vec[keep.index.long]
      ntest.vec = rep(1, length(npositive.vec))
      
      ## re-compute design matrix for time-varying covariates,
      ## since t has been updated 
      re.X.tvarying <- get.Xdesign(dat = dat[index.all.i, ],
                                   dat.exp = dat.exp,
                                   t.curr = t.mat[i+1, index.all.i],
                                   exp.name = exp.name,
                                   season = season[ll],
                                   concep.splines = concep.splines[ll],
                                   df = df[ll],
                                   date.start = date.start,
                                   date.end = date.end,
                                   week.risk = week.risk.list[[ll]])
      X.tvarying.full <- re.X.tvarying$exp.tvarying.long
      X.concep <- re.X.tvarying$concep.long
      if(exposure[ll]){
        X1 <- as.matrix(X.tvarying.full[keep.index.long, 3], ncol = 1)
      }
      X2 <- X.fix.full.list[[ll]][index.long.all.i, ][keep.index.long, ]
      if(concep.splines[ll]|season[ll]){
        X3 <- X.concep[keep.index.long, ]
      }
      
      beta.names.ll <- colnames(beta.list[[ll]])
      if(exposure[ll]){
        index.X1 <- beta.names.ll %in% "trim3.cum"
        pred.comp1 <- as.numeric(X1*beta.list[[ll]][i, index.X1])
      }else{pred.comp1 <- 0}
      index.X2 <- beta.names.ll %in% colnames(X.fix.full.list[[ll]])
      pred.comp2 <- as.numeric(X2 %*% 
                                 matrix(as.numeric(beta.list[[ll]][i, index.X2]), 
                                        ncol = 1))
      if(concep.splines[ll]|season[ll]){
        index.X3 <- beta.names.ll %in% colnames(X.concep)
        pred.comp3 <- as.numeric(
          X3 %*%matrix(as.numeric(beta.list[[ll]][i, index.X3]), ncol = 1))
        pg.parm2 <- pred.comp1 + pred.comp2 + pred.comp3
      }else{
        pg.parm2 <- pred.comp1 + pred.comp2
      }
      numerator <- npositive.vec - 0.5*ntest.vec
      ntotal = length(pg.parm2)
      
      ll.health[[ll]][i] <- sum(dbinom(x = npositive.vec, 
                                       size = ntest.vec, 
                                       prob = expit(pg.parm2), log = T))
      # 3. update omega 
      omega.vec <- rpg(ntotal, ntest.vec, pg.parm2)
      sqrt.omega <- sqrt(omega.vec)
      Omega <- Diagonal(x = omega.vec)
      
      # 4. update betas 
      z.vec <- numerator/(omega.vec)
      ## posterior mean and variance of beta ##
      if(exposure[ll]){
        index.X1 <- which(index.X1)
      }
      index.X2 <- which(index.X2)
      if(concep.splines[ll]|season[ll]){
        index.X3 <- which(index.X3)
      }
      
      M0 <- matrix(NA, ncol = length(beta.names.ll), nrow = length(beta.names.ll))
      rownames(M0) <- colnames(M0) <- beta.names.ll
      if(exposure[ll]){
        M0[index.X1, index.X1] <- crossprod(sqrt.omega*X1)
        M0[index.X1, index.X2] <- crossprod(sqrt.omega*X1, sqrt.omega*X2)
        if(concep.splines[ll]|season[ll]){
          M0[index.X1, index.X3] <- crossprod(sqrt.omega*X1, sqrt.omega*X3)
        }
        
        M0[index.X2, index.X1] <- t(M0[index.X1, index.X2])
        M0[index.X2, index.X2] <- crossprod(sqrt.omega*X2)
        if(concep.splines[ll]|season[ll]){
          M0[index.X2, index.X3] <- crossprod(sqrt.omega*X2, sqrt.omega*X3)
        }
      }else{
        M0[index.X2, index.X2] <- crossprod(sqrt.omega*X2)
        if(concep.splines[ll]|season[ll]){
          M0[index.X2, index.X3] <- crossprod(sqrt.omega*X2, sqrt.omega*X3)
        }
      }
      if(concep.splines[ll]|season[ll]){
        M0[index.X3, index.X1] <- t(M0[index.X1, index.X3])
        M0[index.X3, index.X2] <- t(M0[index.X2, index.X3])
        M0[index.X3, index.X3] <- crossprod(sqrt.omega*X3)
      }
      
      M <- solve(M0 + C.beta.pre[[ll]])
      if(exposure[ll]){
        if(concep.splines[ll]|season[ll]){
          preM <- rbind(rbind(t(sqrt.omega*X1),
                              t(sqrt.omega*X2)),
                        t(sqrt.omega*X3))
        }else{
          preM <- rbind(t(sqrt.omega*X1),
                        t(sqrt.omega*X2))
        }
      }else{
        if(concep.splines[ll]|season[ll]){
          preM <- rbind(t(sqrt.omega*X2),
                        t(sqrt.omega*X3))
        }else{
          preM <- t(sqrt.omega*X2)
        }
      }
      
      m <- M%*%(C.beta.pre[[ll]]%*%c.beta[[ll]] +
                  preM%*%(sqrt(omega.vec)*matrix(z.vec, ncol = 1)))
      beta.list[[ll]][i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    } # end loop over different health models 
    
    if(i < 10){ cat(i, ",") }
    if(i%%100 == 0){ cat(i) }
    
    if(i <= 5000 & i%%100 == 0){
      accept.rate <- colMeans(accept.mat[2:i, ])
      for(l in 1:2){
        if(accept.rate[l] > 0.6){theta2.tun[l] = 1.1*theta2.tun[l]}
        if(accept.rate[l] < 0.2){theta2.tun[l] = 0.8*theta2.tun[l]}
      }
    }
    
  } # end MCMC
  re <- list(beta = beta.list,
             parm = parm.mat,
             gestwk = t.mat,
             alpha = alph.all.list,
             accept.mat = accept.mat,
             nobs.vec = nobs.vec,
             mis.parm = mis.parm.list,
             hazard = hazard.list)
  return(re)
}


get.sesp <- function(x.true, x.error){
  Xtab <- as.data.frame(table(x.true, x.error))
  # colnames(Xtab)[1:2] <- c("x.true", "x.error")
  index.se = Xtab$x.true == 1 & Xtab$x.error == 1
  se = Xtab$Freq[index.se]/(sum(Xtab$Freq[Xtab$x.true == 1]))
  index.sp = Xtab$x.true == 0 & Xtab$x.error == 0
  sp = Xtab$Freq[index.sp]/(sum(Xtab$Freq[Xtab$x.true == 0]))
  return(data.frame("se" = se, "sp" = sp))
}


get.mis.parm <- function(den.w, hazard.vec){
  prob.joint.ind <- do.call(rbind, 
                            lapply(1:18, function(y) den.w[y, ]*hazard.vec))
  se.ind <- sum(prob.joint.ind[1:10, 1:10])/sum(prob.joint.ind[, 1:10])
  sp.ind <- sum(prob.joint.ind[11:18, 11:18])/sum(prob.joint.ind[, 11:18])
  ppv.ind <- sum(prob.joint.ind[1:10, 1:10])/sum(prob.joint.ind[1:10, ])
  npv.ind <- sum(prob.joint.ind[11:18, 11:18])/sum(prob.joint.ind[11:18, ])
  return(c(se.ind, sp.ind, ppv.ind, npv.ind))
}




# A function to simulate data with confounders  
gen.dat.sim <- function(dat.sim.0, X.design.fix,
                        dat.sim.exp, nsims, exp.name, 
                        beta.true,
                        beta.week.true.list, scale.factor = 1,
                        beta.fix.true.list,
                        beta.season.true = list(NULL, NULL),
                        beta.splines.true = list(NULL, NULL),
                        season = c(T, T), 
                        concep.splines = c(T, T), df = c(3, 3),
                        prob.list,
                        val.name = "o3"){
  
  dat.sim.list <- list()
  hazard.list <- list()
  nobs = nrow(dat.sim.0)
  
  beta.week.true.ptb <- beta.week.true.list[[1]]
  beta.week.true.full <- beta.week.true.list[[2]]
  beta.fix.true.ptb <- beta.fix.true.list[[1]]
  beta.fix.true.full <- beta.fix.true.list[[2]]
  beta0.true.ptb <- as.numeric(logit(expit(beta.week.true.ptb)*scale.factor))
  beta0.true.birth <- as.numeric(logit(expit(beta.week.true.full)))
  beta.true.ptb <- beta.true[1]
  beta.true.birth <- beta.true[2]
  
  for(i in 1:nsims){
    
    ## prepare exposure data: 3rd trimester cumulative average 
    week.risk.ptb = 27:36
    exp.sim.mat <- t(apply(dat.sim.exp[,paste0("Wk", 27:36)], 1, dplyr::cummean))
    colnames(exp.sim.mat) <- paste0("wk", 27:36)
    pred.val.1 <- exp.sim.mat*beta.true.ptb
    
    week.risk.birth = 37:44
    exp.sim.vec <- rowMeans(dat.sim.exp[,paste0("Wk", 27:36)])
    pred.val.2 <- exp.sim.vec*beta.true.birth
    
    prred.val.list <- list(pred.val.1, pred.val.2)
    beta.fix.true <- list(beta.fix.true.ptb, beta.fix.true.full)
    pred.val.all <- list()
    for(l in 1:2){
      pred.val.fix <- as.numeric(X.design.fix[[l]] %*% 
                                   matrix(beta.fix.true[[l]], ncol = 1))
      if(season[l] & concep.splines[l]){
        # get first date of conception date 
        cday.df <- data.frame(ID = 1:nobs,
                              concep.date = dat.sim.exp$date - (25-1)*7)
        month.vec <- lubridate::month(cday.df$concep.date)
        cday.df$concep.season <- get.season(month.vec)
        cday.df$fake <- 1
        concep.form <- paste0("ns(concep.date, df = ", df[l], 
                              ", Boundary.knots = c(", as.numeric(date.start),
                              ", ", as.numeric(date.end), ")" ,")")
        X.design.concep <- model.matrix(as.formula(
          paste0("fake ~ ",
                 paste(c("factor(concep.season)", concep.form), 
                       collapse = " + "))), 
          data = cday.df)[,-1]
        beta.concep.true <- c(beta.season.true[[l]], beta.splines.true[[l]])
        pred.val.concep <- as.numeric(X.design.concep %*% 
                                        matrix(beta.concep.true, ncol = 1))
        pred.val.fix.all <- pred.val.fix + pred.val.concep
        if(l == 1){
          pred.val.all[[l]] <- apply(prred.val.list[[l]], 
                                     2, function(x) x + pred.val.fix.all)
        }else if(l == 2){
          pred.val.all[[l]] <- prred.val.list[[l]] + pred.val.fix.all
        }
      }else if(season[l] & (!concep.splines[l])){
        # get conception date 
        cday.df <- data.frame(ID = 1:nobs,
                              concep.date = dat.sim.exp$date - (25-1)*7)
        month.vec <- lubridate::month(cday.df$concep.date)
        cday.df$concep.season <- get.season(month.vec)
        cday.df$fake <- 1
        X.design.concep <- model.matrix(as.formula(
          paste0("fake ~ ",
                 paste(c("factor(concep.season)"), 
                       collapse = " + "))), 
          data = cday.df)[,-1]
        beta.concep.true <- c(beta.season.true[[l]], beta.splines.true[[l]])
        pred.val.concep <- as.numeric(X.design.concep %*% 
                                        matrix(beta.concep.true, ncol = 1))
        pred.val.fix.all <- pred.val.fix + pred.val.concep
        if(l == 1){
          pred.val.all[[l]] <- apply(prred.val.list[[l]], 
                                     2, function(x) x + pred.val.fix.all)
        }else if(l == 2){
          pred.val.all[[l]] <- prred.val.list[[l]] + pred.val.fix.all
        }
      }else if((!season[l]) & concep.splines[l]){
        # get conception date 
        cday.df <- data.frame(ID = 1:nobs,
                              concep.date = dat.sim.exp$date - (25-1)*7)
        cday.df$fake <- 1
        concep.form <- paste0("ns(concep.date, df = ", df[l], 
                              ", Boundary.knots = c(", as.numeric(date.start),
                              ", ", as.numeric(date.end), ")" ,")")
        X.design.concep <- model.matrix(as.formula(
          paste0("fake ~ ",
                 paste(c(concep.form), 
                       collapse = " + "))), 
          data = cday.df)[,-1]
        beta.concep.true <- c(beta.season.true[[l]], beta.splines.true[[l]])
        pred.val.concep <- as.numeric(X.design.concep %*% 
                                        matrix(beta.concep.true, ncol = 1))
        pred.val.fix.all <- pred.val.fix + pred.val.concep
        if(l == 1){
          pred.val.all[[l]] <- apply(prred.val.list[[l]], 
                                     2, function(x) x + pred.val.fix.all)
        }else if(l == 2){
          pred.val.all[[l]] <- prred.val.list[[l]] + pred.val.fix.all
        }
      }else if((!season[l]) & (!concep.splines[l])){
        if(l == 1){
          pred.val.all[[l]] <- apply(prred.val.list[[l]], 
                                     2, function(x) x + pred.val.fix)
        }else if(l == 2){
          pred.val.all[[l]] <- prred.val.list[[l]] + pred.val.fix
        }
      }
    } # end loop over different error-prone outcomes
    
    ### Compute Pr(ti = t) for each observation t \in 27:36
    week.risk.ptb = 27:36
    h0.mat.1 <- expit(t(apply(pred.val.all[[1]], 1, 
                              function(x) beta0.true.ptb + x)))
    colnames(h0.mat.1) <- paste0("wk", week.risk.ptb)
    hazard.curr.mat.1 <- comp_hazard(h0.mat.1)
    prob.cond.birth <- 1 - rowSums(hazard.curr.mat.1)
    
    ### Compute Pr(ti = t) for each observation t \in 37:44
    week.risk.birth = 37:44
    h0.mat.2 <- expit(do.call(rbind, lapply(pred.val.all[[2]], 
                                            function(x) beta0.true.birth + x)))
    colnames(h0.mat.2) <- paste0("wk", week.risk.birth)
    hazard.curr.mat.2 <- comp_hazard(h0.mat.2)
    hazard.curr.mat.2 <- 
      do.call(rbind, 
              lapply(1:nrow(hazard.curr.mat.1), 
                     function(x) prob.cond.birth[x]*hazard.curr.mat.2[x, ]))
    hazard.curr.mat <- cbind(hazard.curr.mat.1, hazard.curr.mat.2)
    colnames(hazard.curr.mat) <- paste0("wk", c(week.risk.ptb, week.risk.birth))
    ## check the row sum of the matrix hazard.curr.mat
    # summary(rowSums(hazard.curr.mat))
    
    ### sample TRUE gestational week based on 
    ### the multinomial probabilities given in hazard.curr.mat
    gest_true_index <- t(apply(hazard.curr.mat, 1, 
                               function(x) rmultinom(1, 1, x)))
    gest_true <- apply(gest_true_index, 1, function(x) (27:44)[which(x == 1)])
    
    ### get two error-prone gestational week using parameters estiamted 
    ### from the real data 
    t.unique = 27:44
    den.w1 <- prob.list[[1]]
    den.w2 <- prob.list[[2]]
    
    gest_error_1 <- get.error.gest(gest.true = gest_true, 
                                   den.mat = den.w1)
    gest_error_2 <- get.error.gest(gest.true = gest_true, 
                                   den.mat = den.w2)
    dat.sim <- data.frame(ID = 1:nobs,
                          CENSUS = dat.sim.exp$GEOID,
                          date = dat.sim.exp$date,
                          gest_true = gest_true,
                          gest_error_1 = gest_error_1,
                          gest_error_2 = gest_error_2)
    # use the first day of Week 27 to compute the date of Week 1
    dat.sim$firstday <- dat.sim$date - 7*(25-1) 
    dat.sim$cday <- dat.sim$firstday + 11
    dat.sim$BIRTHDATE <- dat.sim$firstday + dat.sim$gest_true*7 - 1
    dat.sim$preterm.true <- as.numeric(dat.sim$gest_true <= 36 & 
                                         dat.sim$gest_true >= 27)
    dat.sim.fix <- dat.sim.0
    dat.sim.fix$ID <- 1:nobs
    dat.sim <- left_join(dat.sim, dat.sim.fix, by = c("ID" = "ID"))
    dat.sim <- dat.sim[colnames(dat.sim) != "fake"]
    
    dat.sim.list[[i]] <- dat.sim
    hazard.list[[i]] <- hazard.curr.mat
    if(i%%50 == 0){cat(i, ",")}
  }
  return(list("data" = dat.sim.list, "hazard" = hazard.list))
}
