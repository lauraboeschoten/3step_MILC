setwd("F:/datapackage_threestepMILC/simulationstudy/2_bootstrap")

options(scipen=999)                                                             # scientific notation off

require(brew)                                                                   # load package brew
require(plyr)                                                                   # load package plyr

source("2_fun_create_bootstraps.r")                                             # source the function to create
                                                                                # the bootstrap samples
nsim  <- 1000                                                                   # number of simulations
nboot <- 5                                                                      # number of bootstrap samples
ncomb <- 32                                                                     # number of possible combinations of scores                
nvar  <- 5                                                                      # number of variables

scens <- list(ss     =c(200, 500, 1000),                                        # all simulation conditions
              c.prob =c(0.70, 0.80, 0.90, 0.95, 0.99),
              z.coef =c(0.01, 0.05, 0.10, 0.20),
              q.coef =c(0.5, 0.6224593, 0.7310586, 0.8807971))

scen.mat <- expand.grid(scens)                                                  # list with all simulation conditions


for (i in 1:nrow(scen.mat)){                                                    # loop over scenarios 
  
  set.seed(123)                                                                 # every scenario starts with same seed
  cat(i)                                                                        # count scenarios
  
  folder.name <- paste0("F:/datapackage_threestepMILC/simulationstudy/1_createdata/pops_", 
                        scen.mat[i,1],"_",                                      # specify folder where data can be found
                        substring(scen.mat[i,2], 3),"_",
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3))
  
  setwd(folder.name)                                                            # set location to this folder
  
  data.name   <- paste0("dat_",scen.mat[i,1],"_",                               # specify name of the dataset
                         substring(scen.mat[i,2], 3),"_",
                         substring(scen.mat[i,3], 3),"_",
                         substring(scen.mat[i,4], 3),".txt")
  
  
  ssize             <- scen.mat[i,1]                                            # specify sample size of this scenario
  samples           <- read.delim(data.name, header=T, sep="\t")                # load the sample data
  bootstrap.samples <- create.bootstrap.samples(nsim, nboot, ssize, ncomb,      # run the bootstrap.samples function
                                                samples, i)
}
  