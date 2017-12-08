setwd("F:/datapackage_threestepMILC/simulationstudy/4_impute")

require(brew)                                                                   # load brew package
require(psych)                                                                  # load psych package
require(testthat)                                                               # load testthat package
require(plyr)                                                                   # load plyr package

options(scipen=999)                                                             # scientific notation off

source("4_createprofile.R")                                                     # source function to create profile
source("4_posteriors.R")                                                        # source function to create posteriors
source("4_impute.R")                                                            # source function to create imputations

nsim  <- 1000                                                                   # number of simulations
nboot <- 5                                                                      # number of bootstrap samples
ncomb <- 32                                                                     # number of possible combinations of scores                
nvar  <- 5                                                                      # number of variables

scens <- list(ss     =c(200, 500, 1000),                                        # all simulation conditions
              c.prob =c(0.70, 0.80, 0.90, 0.95, 0.99),
              z.coef =c(0.01, 0.05, 0.10, 0.20),
              q.coef =c(0.5, 0.6224593, 0.7310586, 0.8807971))

scen.mat <- expand.grid(scens) 


for (i in 1:nrow(scen.mat)){                                                    # loop over all scenario's
  
  set.seed(123)                                                                 # start every scenario with same seed
  cat(i)                                                                        # count scenario's
  ssize <- scen.mat[i,1]                                                        # specify sample size
  
  folder.name <- paste0("F:/datapackage_threestepMILC/simulationstudy/3_LC_models/LC_", 
                        scen.mat[i,1],"_",                                      # specify folder where the LC
                        substring(scen.mat[i,2], 3),"_",                        # model is located
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3))
  
  setwd(folder.name)
  
  profiles   <- vector("list", nsim)                                            # create list to place profiles
  posteriors <- vector("list", nsim)                                            # create list to place posteriors
  
  for (s in 1:nsim){
    profiles[[s]]        <- vector("list",nboot)                                # create a list for the 5 profiles for each 
                                                                                # for each iteration in the simulation
    for (j in 1:nboot){                                                         # for each bootstrap sample within the iteration
      select.model       <- paste0("model.",s,j,".lst")                         # get the corresponding model output
      profiles[[s]][[j]] <- create.profile(select.model)                        # and use the profile function to create the profile
    }
    
    posteriors[[s]]      <- posterior.probabilities(profiles[[s]])              # we now use the profile to create the posteriors
  }
  
  folder.name <- paste0("F:/datapackage_threestepMILC/simulationstudy/4_impute/files_", 
                        scen.mat[i,1],"_",                                      # specify the folder to put the imputed data in
                        substring(scen.mat[i,2], 3),"_",
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3))
  setwd(folder.name)
  
  for (s in 1:nsim){                                                            # for every iteration in the simulation
    for (j in 1:nboot){                                                         # and every bootstrap sample in the iteration
      write.table(profiles[[s]][[j]], paste0("profile",s,j,".txt"), sep="\t")   # write a file containing the profile   
    }
    write.table(posteriors[[s]], paste0("posteriors",s,".txt"), sep="\t")       # and write a file containing the posteriors
  }
  
  for (s in 1:nsim){                                                            # for each iteration in the simulation
    posterior        <- posteriors[[s]]                                         # put the posteriors in a list
    bootstrap.sample <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/2_bootstrap/samples_", 
                               scen.mat[i,1],"_",                               # open the datafile containing the bootstraps
                               substring(scen.mat[i,2], 3),"_",
                               substring(scen.mat[i,3], 3),"_",
                               substring(scen.mat[i,4], 3),
                               "\\bootstrap_sample_",s,".txt"), header=TRUE)
    imputations      <- impute.samples(posterior, bootstrap.sample, ssize)      # create imputations
    
    write.table(imputations, paste0("imputations_",s,".txt"), sep="\t",         # save a dataset containing the multiple imputations
                row.names = FALSE, quote = FALSE)
  }
}
