setwd("F:/datapackage_threestepMILC/simulationstudy/10_LC_results_ML2")

require(brew)                                                                   # load brew package
require(plyr)                                                                   # load plyr package

options(scipen=999)                                                             # scientific notation off


source("10_results.R")                                                          # source results functions
source("10_mean_results.R")

nsim  <- 1000                                                                   # number of simulations
nboot <- 5                                                                      # number of bootstrap samples
ncomb <- 32                                                                     # number of possible combinations of scores                
nvar  <- 5                                                                      # number of variables


scens <- list(ss     =c(200, 500, 1000),                                        # all simulation conditions
              c.prob =c(0.70, 0.80, 0.90, 0.95, 0.99),
              z.coef =c(0.01, 0.05, 0.10, 0.20),
              q.coef =c(0.5, 0.6224593, 0.7310586, 0.8807971))

scen.mat <- expand.grid(scens)                                                  

for (i in 1:nrow(scen.mat)){                                                    # for every simulation condition
  
  ssize <- scen.mat[i,1]                                                        # save sample size
  
  cat(i)
  set.seed(123)
  
  create.results(i)                                                             # run create.results function
  mean.results(i)                                                               # and mean results function
}


