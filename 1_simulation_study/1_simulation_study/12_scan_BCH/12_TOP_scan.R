setwd("F:/datapackage_threestepMILC/simulationstudy/12_scan_BCH")

require(brew)                                                                   # load brew package
require(plyr)                                                                   # load plyr package

options(scipen=999)                                                             # scientific notation off

nsim  <- 1000                                                                   # number of simulations
nboot <- 5                                                                      # number of bootstrap samples
ncomb <- 32                                                                     # number of possible combinations of scores                
nvar  <- 5                                                                      # number of variables


scens <- list(ss     =c(200, 500, 1000),                                        # all simulation conditions
              c.prob =c(0.70, 0.80, 0.90, 0.95, 0.99),
              z.coef =c(0.01, 0.05, 0.10, 0.20),
              q.coef =c(0.5, 0.6224593, 0.7310586, 0.8807971))


scen.mat <- expand.grid(scens) 

for (i in 1:nrow(scen.mat)){
  cat(i)
  ssize <- scen.mat[i,1]
  negatives <- matrix(NA, 1000, 5)                                              # create matrix to store results for bootstraps
  negativeboot <- matrix(NA, 1000, 1)                                           # create matrix to store final results
  
  for (s in 1:nsim){                                                            # for every iteration
    BCH <- vector("list", 5)                                                    # create list of 5
    for(c in 1:nboot){                                                          # for every bootstrap sample
      BCH[[c]] <- round(read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/11_BCH/BCH_",
                                        scen.mat[i,1],"_",                      # load amatrices
                                        substring(scen.mat[i,2], 3),"_",
                                        substring(scen.mat[i,3], 3),"_",
                                        substring(scen.mat[i,4], 3),
                                        "/BCH_amatrix_",s,c,".txt"), header=T), 8)
      negatives[s,c] <- any(rowSums(BCH[[c]])==0)                               # see if there are negative values
      nega <- negatives[s,]                                                     
      negativeboot[s] <- any(nega==TRUE)                                       
    }
  }
    write.table(negatives, paste0("F:/datapackage_threestepMILC/simulationstudy/12_scan_BCH/", 
                                     scen.mat[i,1],"_",                         # save which iteration contains negatives
                                     substring(scen.mat[i,2], 3),"_",
                                     substring(scen.mat[i,3], 3),"_",
                                     substring(scen.mat[i,4], 3),
                                     "negatives.txt"), sep="\t", col.names=TRUE, row.names=FALSE)
    write.table(negativeboot, paste0("F:/datapackage_threestepMILC/simulationstudy/12_scan_BCH/", 
                                  scen.mat[i,1],"_",
                                  substring(scen.mat[i,2], 3),"_",
                                  substring(scen.mat[i,3], 3),"_",
                                  substring(scen.mat[i,4], 3),
                                  "negativeboot.txt"), sep="\t", col.names=TRUE, row.names=FALSE)
}

