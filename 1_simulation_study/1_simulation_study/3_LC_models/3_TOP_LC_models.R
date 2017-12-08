setwd("F:/datapackage_threestepMILC/simulationstudy/3_LC_models")

options(scipen=999)                                                             # scientific notation off

require(brew)                                                                   # load package brew
require(plyr)                                                                   # load package plyr

source("3_template_LCmodel.r")                                                  # source the lG template

nsim  <- 1000                                                                   # number of simulations
nboot <- 5                                                                      # number of bootstrap samples
ncomb <- 32                                                                     # number of possible combinations of scores                
nvar  <- 5                                                                      # number of variables
brew  <- "3_brew_LCmodel.brew"                                                  # name of the brew file

scens <- list(ss     =c(200, 500, 1000),                                        # all simulation conditions
              c.prob =c(0.70, 0.80, 0.90, 0.95, 0.99),
              z.coef =c(0.01, 0.05, 0.10, 0.20),
              q.coef =c(0.5, 0.6224593, 0.7310586, 0.8807971))

scen.mat <- expand.grid(scens)                                                  # list with all simulation conditions

for (i in 1:nrow(scen.mat)){                                                    # loop over all scenarios
  
  set.seed(123)                                                                 # set seed for every scenario
  cat(i)                                                                        # count the scenario's 
  
  y.min  = log((1-scen.mat[i,2])/scen.mat[i,2])                                 # obtain coefficients to put in
  y.plus = log(scen.mat[i,2]/(1-scen.mat[i,2]))                                 # brew file for indicators
  
  folder.name <- paste0("F:/datapackage_threestepMILC/simulationstudy/3_LC_models/LC_", 
                        scen.mat[i,1],"_",                                      # folder to place the LC models
                        substring(scen.mat[i,2], 3),"_",
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3))
  setwd(folder.name)                                                            # set to the folder
  
  bmods       <- list(nsim=1:nboot, nmod=1:nsim)                                # create names for the outfiles
  bmods.mat   <- expand.grid(bmods)
  bmods.files <- vector(length=(nrow(bmods.mat)))
  
  for (b in 1:nrow(bmods.mat)){
    bmods.files[b] <- paste0("model.",bmods.mat[b,2],bmods.mat[b,1])
  }
  
  bcnames <- c("b1","b2","b3","b4","b5")                                        # specify the names of the bootstrap samples                                    
  
  for(s in 1:nsim){                                                             # run over all iterations of the simulation
    outfiles <- bmods.files[(1+(s-1)*5): ((1+(s-1)*5)+4)]                       # names of the outfiles
    for (j in 1:nboot){                                                         # run over the bootstraps
    envir       <- new.env()
    template_LC(template.path      = brew,                                      # run the LG template
                envir              = envir, 
                temp.filename.base = outfiles[j], 
                s                  = s, 
                j                  = j,
                bcname             = bcnames[j], 
                y.min              = y.min, 
                y.plus             = y.plus)
    }
  }
}
