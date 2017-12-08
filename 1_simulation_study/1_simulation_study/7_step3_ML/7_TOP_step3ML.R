setwd("F:/datapackage_threestepMILC/simulationstudy/7_step3_ML")

require(brew)                                                                   # load brew package
require(psych)                                                                  # load psych package
require(testthat)                                                               # load testthat package
require(plyr)                                                                   # load plyr package

options(scipen=999)                                                             # scientific notation off

source("7_template_step3ML.R")                                                  # source the template function

nsim  <- 1000                                                                   # number of simulations
nboot <- 5                                                                      # number of bootstrap samples
ncomb <- 32                                                                     # number of possible combinations of scores                
nvar  <- 5                                                                      # number of variables

brew  <- "7_step3ML_brew.brew"                                                  # brew file

scens <- list(                                                                  # all simulation conditions
  ss     =c(200, 500, 1000),
  c.prob =c(0.70, 0.80, 0.90, 0.95, 0.99),
  z.coef =c(0.01, 0.05, 0.10, 0.20),
  q.coef =c(0.5, 0.6224593, 0.7310586, 0.8807971))

scen.mat <- expand.grid(scens) 

for (i in 1:nrow(scen.mat)){                                                    # for every scenario
  
  cat(i)
  set.seed(123)
  ssize <- scen.mat[i,1]                                                        # save sample size
  
  zcoef = log((2*scen.mat[i,3])/(1-(2*scen.mat[i,3])))                          # obtain coefficients to put in brew 
  qcoef = log(scen.mat[i,4]/(1-scen.mat[i,4]))                                  # file for the covariates
  covar = qcoef/zcoef
  
  folder.name <- paste0("F:/datapackage_threestepMILC/simulationstudy/7_step3_ML/step3ML_", 
                        scen.mat[i,1],"_",                                      # specify folder for output
                        substring(scen.mat[i,2], 3),"_",
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3))
  
  setwd(folder.name)
  
  bmods       <- list(nsim=1:nboot, nmod=1:nsim)                                # create outfile names
  bmods.mat   <- expand.grid(bmods)
  bmods.files <- vector(length=(nrow(bmods.mat)))
  
  for (b in 1:nrow(bmods.mat)){
    bmods.files[b] <- paste0("model.",bmods.mat[b,2],bmods.mat[b,1])
  }
  
  for(s in 1:nsim){
    cat(s)
    outfiles <- bmods.files[(1+(s-1)*5): ((1+(s-1)*5)+4)]                       # define outfiles

    ctables  <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/6_classification_tables\\table_", 
                                  scen.mat[i,1],"_",                            # load the classification tables
                                  substring(scen.mat[i,2], 3),"_",
                                  substring(scen.mat[i,3], 3),"_",
                                  substring(scen.mat[i,4], 3),
                                  "\\class_mat_",s,".txt"), header=T)
    ctables  <- as.matrix(ctables)                                              # put the classification tables
    ctables  <- matrix(ctables, ncol = ncol(ctables), dimnames = NULL)          # in a matrix
    imps <- c("w1","w2","w3","w4","w5")                                         # specify the imputed variables
    
    for (j in 1:nboot){                                                         # apply the three-step ML method
      envir       <- new.env()                                                  # to all bootstrap samples
      template_step3(template.path   = brew,
                  envir              = envir, 
                  temp.filename.base = outfiles[j], 
                  s                  = s, 
                  j                  = j,
                  bcname             = imps[j], 
                  a                  = ctables[j,],
                  qcoef              = qcoef,
                  zcoef              = zcoef,
                  covar              = covar)
    }
  }
}
