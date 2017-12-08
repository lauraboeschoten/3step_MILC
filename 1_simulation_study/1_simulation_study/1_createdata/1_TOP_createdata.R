setwd("F:/datapackage_threestepMILC/1_simulation_study/1_createdata")

options(scipen=999)                                                             # scientific notation off

require(brew)                                                                   # load package brew
require(plyr)                                                                   # load package plyr

source("1_fun_createdata.R")                                                    # load the createdata function
source("1_createdata_template.R")                                               # load the template function

nsim     <- 1000                                                                # number of simulations
nboot    <- 5                                                                   # number of bootstrap samples
ncomb    <- 32                                                                  # number of possible combinations of scores                
nvar     <- 5                                                                   # number of variables
brew     <- "1_create_data.brew"                                                # name of brew file

scens    <- list(ss     =c(200, 500, 1000),                                     # all simulation conditions
                 c.prob =c(0.70, 0.80, 0.90, 0.95, 0.99),
                 z.coef =c(0.01, 0.05, 0.10, 0.20),
                 q.coef =c(0.5, 0.6224593, 0.7310586, 0.8807971))

scen.mat <- expand.grid(scens)                                                  # list with all simulation conditions


createdata(scen.mat, nsim, brew)                                                # run the createdata function
                                                                                # on all specified simulation
                                                                                # conditions