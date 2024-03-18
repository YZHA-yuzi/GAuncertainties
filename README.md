README
================

In this document, we illustrate the implementation of the Bayesian
hierarchical model proposed in the paper titled “Time-to-Event Analysis
of Preterm Birth Accounting for Gestational Age Uncertainties” with a
simulated toy dataset.

We also provided instructions on reproducing simulation results included
in Section 4 in the manuscript in the second half of this document.

# Toy Example

**Note**: Please make sure the working directory is set to the folder
where the R script “FUNs.R” is located.

Load in packages

``` r
library(Rcpp)
library(RcppArmadillo)
library(dplyr) # left_join
library(tidyr) # gather
library(BayesLogit) # rpg
library(Matrix) # Diagonal
library(MASS) # mvrnorm
library(splines)
```

Read in self-defined functions

``` r
sourceCpp("FUNs_rcpp.cpp")
source("FUNs.R")
```

Load the simulated toy dataset

This toy dataset contains 5000 observations and 6 variables.

- `ID` = the ID of the observation
- `CENSUS` = 11 digits FIPS code for the census tract (sampled from the
  our real data and will be used when compute the time-varying exposure
  of interest)
- `gest_error_1` = 1st error-prone gestational age (GA) estimate
- `gest_error_2` = 2nd error-prone GA estimate
- `BIRTHDATE` = birth date
- `SMOKE` = smoking status

``` r
load("./data/Data_toy.rda") # dat.toy
dim(dat.toy)
```

    ## [1] 5000    6

``` r
head(dat.toy) 
```

    ##   ID      CENSUS gest_error_1 gest_error_2  BIRTHDATE SMOKE
    ## 1  1 20125950300           38           39 2010-11-01     0
    ## 2  2 20091052303           37           36 2010-12-27     0
    ## 3  3 20191962500           38           38 2010-12-03     0
    ## 4  4 20169000500           37           36 2011-08-03     0
    ## 5  5 20033967600           38           38 2011-03-17     1
    ## 6  6 20013480600           35           33 2011-07-23     0

Load the exposure data (ozone concentrations from our real data
application)

``` r
load("./data/Data_exp_toy.rda") # dat.exp.toy
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
  # if(i%%100 == 0){cat(i, ",")}
}
# "date" in dat.exp.wide: first date of week 25
dat.exp <- do.call(rbind.data.frame, dat.exp.wide.list)
dat.exp$GEOID <- as.character(dat.exp$GEOID)
```

Fit the proposed model using the self-defined R function `fit.sim`.

Arguments in the function `fit.sim`:

- `dat` = a dataframe contains observed error-prone GA estimates, birth
  dates, and time-invariant covariates. Please make sure the variable
  for the birth date is named as “BIRTHDATE”. In addition, the dataframe
  must contain a variable for observations’ ID. In this toy example, a
  variable named “CENSUS” is required when computing the time-varying
  exposure of interest (i.e., cumulative ozone average during 3rd
  trimester of pregnancy).

- `dat.exp` = a dataframe contains the environmental exposure of
  interest in a wide-format.

- `niter` = the number of MCMC iterations (the default value is 25,000).

- `burn_in` = the number of burn-in samples (the default value is
  20,000).

- `seed` = the random seed.

- `exp.name` = a character used to specify the time-varying exposure of
  interest. In this function, the only value allowed is “trim3.cum”
  representing the cumulative average during 3rd trimester of pregnancy.

- `get.Xdesign` = a function that used to compute time-varying exposures
  of interest. The function `get.Xdesign.tvarying` included in the R
  script `FUNs.R` is defined to compute time-varying functions of
  interest.

- `cov.name.vec` = a list of length 2 contains names of covariates
  included in the two hazard models. The 1-st element is for the hazard
  model for preterm birth, and the 1-nd element is for the hazard model
  for term birth. Note that `factor(Week)` must be included for
  week-specific intercepts.

- `date.start` = the first date of the study time period of interest in
  the format “YYYY-MM-DD”.

- `date.end` = the last date of the study time period of interest in the
  format “YYYY-MM-DD”.

- `season` = a vector of length two contains binary indicators
  indicating whether conception season is controlled for in the hazard
  model. The 1-st element is for the hazard model for preterm birth, and
  the 1-nd element is for the hazard model for term birth. For example,
  `season = c(T, F)` implies conception season is included for the
  hazard model for preterm birth, but is not included for the model for
  term birth. The default value is `c(T, T)`.

- `concep.splines` = a vector of length two contains binary indicators
  indicating whether conception dates are controlled for in the hazard
  model using natural cubic splines with the degrees of freedom
  specified via the argument `df`. The 1-st element is for the hazard
  model for preterm birth, and the 1-nd element is for the hazard model
  for term birth. For example, `concep.splines = c(T, F)` implies
  conception date is included for the hazard model for preterm birth,
  but is not included for the model for term birth. The default value is
  `c(T, T)`.

- `df` = a vector of length two specifies the degrees of freedom of
  natural cubic splines used for modeling conception season (the default
  value is `c(3, 3)`). The 1-st element is for the hazard model for
  preterm birth, and the 1-nd element is for the hazard model for term
  birth. Note that if there is `FALSE` in the argument `concep.splines`,
  then the corresponding element in `df` should be set as `NULL`.

- `exposure` = a vector of length two contains binary indicators
  indicating whether time-varying environmental exposure is controlled
  for in the hazard model. The 1-st element is for the hazard model for
  preterm birth, and the 1-nd element is for the hazard model for term
  birth. For example, `exposure = c(T, F)` implies time-varying
  environmental exposure is included for the hazard model for preterm
  birth, but is not included for the model for term birth. The default
  value is `c(T, T)`. **Note that if the time-varying exposure is
  included, this function currently only allows to consider the
  cumulative average during 3rd trimester of pregnancy.**

- `dist` = a character specifies the value of
  ![d](https://latex.codecogs.com/png.latex?d "d") included in
  exponential covariance functions,
  ![\left(\Sigma\_{l}\right)\_{k,k'} = \sigma\_{l}^2 \exp \left( \frac{-\|k-k'\|^d}{\phi\_{l}} \right)](https://latex.codecogs.com/png.latex?%5Cleft%28%5CSigma_%7Bl%7D%5Cright%29_%7Bk%2Ck%27%7D%20%3D%20%5Csigma_%7Bl%7D%5E2%20%5Cexp%20%5Cleft%28%20%5Cfrac%7B-%7Ck-k%27%7C%5Ed%7D%7B%5Cphi_%7Bl%7D%7D%20%5Cright%29 "\left(\Sigma_{l}\right)_{k,k'} = \sigma_{l}^2 \exp \left( \frac{-|k-k'|^d}{\phi_{l}} \right)"),
  used in the outcome misclassification model, for
  ![l = 1,2](https://latex.codecogs.com/png.latex?l%20%3D%201%2C2 "l = 1,2"),
  ![k, k' \in \\27, \dots, 44\\](https://latex.codecogs.com/png.latex?k%2C%20k%27%20%5Cin%20%5C%7B27%2C%20%5Cdots%2C%2044%5C%7D "k, k' \in \{27, \dots, 44\}"),
  where ![l](https://latex.codecogs.com/png.latex?l "l") is the index
  for the error-prone GA. `dist` = `"1st"` indicates
  ![d = 1](https://latex.codecogs.com/png.latex?d%20%3D%201 "d = 1"),
  and `dist` = `"2nd"` indicates
  ![d = 2](https://latex.codecogs.com/png.latex?d%20%3D%202 "d = 2").

- `ntol` = the number of weeks that the error-prone GA is allowed to
  deviate from the true GA (the default value is 3).

**Notes**: In this toy example, for the purpose of the demonstration, we
only ran 10 MCMC iterations and discarded the first 5 samples as
burn-in, because the computation time is a concern. Specifically, the
run-time for the default 25,000 MCMC iterations is about 5.5 hours which
was approximated using a Mac with Apple M2 Pro and 32 Gb memory.

``` r
# Apply the proposed model to the simulated toy dataset.
# 10 MCMC iterations take about 8s on a Mac with Apple M2 pro and 32 Gb memory.
nMCMC = 10
burn_in = 5

# time-invariant exposure (i.e., smoking) was controlled for both hazard models.
cov.sim.0 <- c("factor(SMOKE)")

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
               dist.name = "1st", ntol = 3)
```

The output `fit` is a list of length 8 which contains:

- `beta`: a list of length 2 contains samples (before the burn-in
  samples are discarded) of regression coefficients in the two discrete
  hazard models. The 1-st element is for the hazard model for preterm
  birth, and the 2-nd element is for the hazard model for term birth.

- `parm`: a matrix contains samples (before the burn-in samples are
  discarded) of parameters included in exponential covariance functions
  used in the outcome misclassification model; `theta1_l` =
  ![\sigma\_{l}^2](https://latex.codecogs.com/png.latex?%5Csigma_%7Bl%7D%5E2 "\sigma_{l}^2"),
  `theta2_l` =
  ![\phi\_{l}](https://latex.codecogs.com/png.latex?%5Cphi_%7Bl%7D "\phi_{l}"),
  for ![l=1,2](https://latex.codecogs.com/png.latex?l%3D1%2C2 "l=1,2"),
  in eqn.(1) in the manuscript.

- `gestwk`: a matrix contains samples (before the burn-in samples are
  discarded) of true gestational ages.

- `alpha`: a list of length 2 with each element further contains a list
  of length 18. The list nested in the list `alpha` contains random
  effects
  ![\alpha\_{l,kj}](https://latex.codecogs.com/png.latex?%5Calpha_%7Bl%2Ckj%7D "\alpha_{l,kj}")
  given in eqn. (1) in the manuscript.

- `accept.mat`: a matrix contains binary indicators of whether the
  proposed value was accepted in M-H algorithms used for updating
  ![\phi\_{l}](https://latex.codecogs.com/png.latex?%5Cphi_%7Bl%7D "\phi_{l}"),
  for
  ![l = 1,2](https://latex.codecogs.com/png.latex?l%20%3D%201%2C2 "l = 1,2").

- `nobs.vec`: a vector contains the number of observations included into
  the study population of interest at each MCMC iteration.

- `mis.parm`: contains individual-specific parameters summarizing
  misclassifications of gestational ages in term of identifying preterm
  birth. It is a list of length 2 with each element further contains a
  list of length ![M](https://latex.codecogs.com/png.latex?M "M"), where
  ![M](https://latex.codecogs.com/png.latex?M "M") is the number of
  posterior samples after the burn-in period is discarded. Specifically,
  for each observation, we compute SE = sensitivity =
  ![Pr(Y\_{li} \< 37 \| T_i \< 37)](https://latex.codecogs.com/png.latex?Pr%28Y_%7Bli%7D%20%3C%2037%20%7C%20T_i%20%3C%2037%29 "Pr(Y_{li} < 37 | T_i < 37)"),
  SP = specificity =
  ![Pr(Y\_{li} \>= 37 \| T_i \>= 37)](https://latex.codecogs.com/png.latex?Pr%28Y_%7Bli%7D%20%3E%3D%2037%20%7C%20T_i%20%3E%3D%2037%29 "Pr(Y_{li} >= 37 | T_i >= 37)"),
  PPV = positive predictive value =
  ![Pr(T_i \< 37 \| Y\_{li} \< 37)](https://latex.codecogs.com/png.latex?Pr%28T_i%20%3C%2037%20%7C%20Y_%7Bli%7D%20%3C%2037%29 "Pr(T_i < 37 | Y_{li} < 37)"),
  and NPV = negative predictive value =
  ![Pr(T_i \>= 37 \| Y\_{li} \>= 37)](https://latex.codecogs.com/png.latex?Pr%28T_i%20%3E%3D%2037%20%7C%20Y_%7Bli%7D%20%3E%3D%2037%29 "Pr(T_i >= 37 | Y_{li} >= 37)"),
  where
  ![T\_{i}](https://latex.codecogs.com/png.latex?T_%7Bi%7D "T_{i}") is
  the true GA, and
  ![Y\_{li}](https://latex.codecogs.com/png.latex?Y_%7Bli%7D "Y_{li}")
  is ![l](https://latex.codecogs.com/png.latex?l "l")-th error-prone GA
  of ![i](https://latex.codecogs.com/png.latex?i "i")-th observation,
  for
  ![l = 1, 2](https://latex.codecogs.com/png.latex?l%20%3D%201%2C%202 "l = 1, 2").

- `hazard`: a list of length
  ![M](https://latex.codecogs.com/png.latex?M "M"), where
  ![M](https://latex.codecogs.com/png.latex?M "M") is the number of
  posterior samples after the burn-in period is discarded;
  ![m](https://latex.codecogs.com/png.latex?m "m")-th element in the
  list `hazard` contains a matrix of dimension
  ![n](https://latex.codecogs.com/png.latex?n "n") by 18, each row
  represents individual-specific probability
  ![Pr(T_i = t)](https://latex.codecogs.com/png.latex?Pr%28T_i%20%3D%20t%29 "Pr(T_i = t)"),
  for
  ![t = 27, \dots, 44](https://latex.codecogs.com/png.latex?t%20%3D%2027%2C%20%5Cdots%2C%2044 "t = 27, \dots, 44"),
  ![i = 1, \dots, n](https://latex.codecogs.com/png.latex?i%20%3D%201%2C%20%5Cdots%2C%20n "i = 1, \dots, n")
  (![n](https://latex.codecogs.com/png.latex?n "n") = the number of
  observations).

Let’s take a look of posterior samples (the first 5 MCMC iterations) of
regression coefficients included in the hazard model for preterm birth.
Specifically, there are 12 parameters in total containing 10
week-specific intercepts (`factor(Week)t`), for
![t = 27, \dots, 36](https://latex.codecogs.com/png.latex?t%20%3D%2027%2C%20%5Cdots%2C%2036 "t = 27, \dots, 36"),
1 coefficient for the time-varying exposure (`trim3.cum`), and 1
coefficient for the time-invariant exposure (`factor(SMOKE)1`).

``` r
head(fit$beta[[1]])
```

    ##    trim3.cum factor(Week)27 factor(Week)28 factor(Week)29 factor(Week)30
    ## 1 0.02580238      -6.372968      -5.557214      -5.545919      -5.111862
    ## 2 0.02568448      -6.328621      -5.554749      -5.647020      -5.130221
    ## 3 0.02146494      -6.049382      -5.404133      -5.445572      -4.924018
    ## 4 0.02143971      -6.060304      -5.465734      -5.428333      -4.918771
    ## 5 0.02390152      -6.156415      -5.496100      -5.600031      -4.942563
    ## 6 0.02211815      -5.970233      -5.452723      -5.530572      -4.945926
    ##   factor(Week)31 factor(Week)32 factor(Week)33 factor(Week)34 factor(Week)35
    ## 1      -4.794967      -4.364012      -4.013955      -3.723886      -2.943456
    ## 2      -4.799917      -4.382178      -4.072129      -3.657333      -2.910938
    ## 3      -4.610733      -4.267317      -3.873896      -3.538430      -2.916713
    ## 4      -4.606527      -4.348453      -3.891108      -3.588422      -2.915954
    ## 5      -4.628889      -4.332790      -3.996844      -3.712439      -2.960796
    ## 6      -4.679006      -4.271380      -3.920701      -3.646105      -2.886307
    ##   factor(Week)36 factor(SMOKE)1
    ## 1      -2.254837      0.3727102
    ## 2      -2.297857      0.3608940
    ## 3      -2.157926      0.4118332
    ## 4      -2.144735      0.3769862
    ## 5      -2.295315      0.4410398
    ## 6      -2.121778      0.4186761

Let’s take a look of posterior samples (the first 5 MCMC iterations) of
true gestational ages for the first 3 observations.

``` r
fit$gestwk[1:5, 1:3]
```

    ##      [,1] [,2] [,3]
    ## [1,]   38   37   38
    ## [2,]   38   37   38
    ## [3,]   38   37   38
    ## [4,]   38   37   38
    ## [5,]   38   36   38

# Reproducibility

## R scripts

- `FUNs.R` (self-defined R functions)
- `FUNs_rcpp.cpp` (self-defined functions)
- `Simulations_data.R` (data generation)
- `Simulations.R` (fit the proposed model)
- `Simulations_logistic.R` (conduct the conventional preterm birth
  analysis)
- `Simulations_sumres.R` (summarize simulation results)
- `Sim_res_manu_Table1.R` (generate Table 1 included in the manuscript
  based on intermediate results provided in the subfolder
  `./interres/for_manu/`)

## Bash shell scripts

- `sim_data.sh`
- `sim_proposed.sh`
- `sim_logistic.sh`
- `sim_sumres.sh`
- `sim_Tab1.sh`

## Intermediate results

Given a simulated dataset including 5,000 observations, the proposed
method with 25,000 MCMC iterations takes about 6 hours to run (the
computation time was approximated using a Mac running M2 chip). Thus, we
have provided intermediate results that can be used to generate the
table 1 included in the manuscript.

## Generate simulated datasets

We provided an R script `Simulations_data.R` to generate simulated
datasets used for the simulation study included in the manuscript. The
time-varying and time-invariant exposures used for the data generation
were also provided in the sub-folder `data`. By running this R script,
simulated datasets will be generated and saved to the sub-folder `data`.
Under each simulation scenario, 100 simulated datasets will be
generated.

## Instructions for use

There are two options to run those provided R scripts: (1) RStudio or
(2) Bash Shell scripts. When running R codes in RStudio, please create a
new R project and make sure R scripts are stored at the same folder
where the created R project `*.Rproj` is located; when running R codes
using provided bash shell scripts, please set the working directory to
the folder where R scripts and bash shell scripts are located.

Run-time was approximated using a Mac with Apple M2 Pro chip and 32 Gb
memory.

Steps to produce results presented and discussed in Section 4 are given
below.

**Step 1**. Generate data

Run `Simulations_data.R` in RStudio; or execute the bash script
`sim_data.sh` using the command `bash sim_data.sh` in terminal.

The simulated data `Data_sims.rda` will be saved to the sub-folder
`./data/`. In addition, true values of some individual-specific
quantities, for example, the probability of being identified as preterm
birth
![Pr(T\_{i} \le 36)](https://latex.codecogs.com/png.latex?Pr%28T_%7Bi%7D%20%5Cle%2036%29 "Pr(T_{i} \le 36)")
are also saved to the sub-folder `./data/`.

**Step 2**. Fit the proposed model based on simulated data

Run `Simulations.R` in RStudio (arguments `sceindex` and `simindex` have
to be specified); or execute the bash script `sim_proposed.sh`. One
simulation takes about 6 hours to run, we have a total of 300
simulations (3 simulation scenarios, and 100 simulations were conducted
within each scenario).

`Res_S*j_sim*i.rda` and `Res_PTBrate_S*j_sim*i.rda` will be saved to the
sub-folder named `./interres/`, where
![j=1,...,3](https://latex.codecogs.com/png.latex?j%3D1%2C...%2C3 "j=1,...,3")
and
![i=1,...,100](https://latex.codecogs.com/png.latex?i%3D1%2C...%2C100 "i=1,...,100").

**NOTE**: 25,000 MCMC iterations is the default for fitting the proposed
model.

**Step 3**. Conduct the conventional preterm birth analysis based on
simulated data

Run `Simulations_logistic.R` in RStudio (arguments `sceindex` have to be
specified); or execute the bash script `sim_logistic.sh`.

`Res_glm_S*j.rda` will be saved to the sub-folder `./interres/`, where
![j=1,...,3](https://latex.codecogs.com/png.latex?j%3D1%2C...%2C3 "j=1,...,3").

**Step 4**. Summarize simulation results based on results generated in
Steps 2 and 3

Run `Simulations_sumres.R` in RStudio; or execute the bash script
`sim_sumres.sh`. `Tab_simres.csv` will be saved to the sub-folder
`./interres/`

**Step 5**. Summarize simulation results based on intermediate results
provided upfront

Run `Sim_res_manu_Table1.R` in RStudio; or execute the bash script
`sim_Tab1.sh`. Intermediate results were provided upfront located at the
sub-folder `./interres/for_manu/`

`Table1.csv` will be saved to the sub-folder `./interres/for_manu/`
