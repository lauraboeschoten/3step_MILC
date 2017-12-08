setwd("F:/datapackage_threestepMILC/simulationstudy/13_impute_BCH")

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

for (i in 1:nrow(scen.mat)){                                                    # for every simulation condition
  cat(i)
  set.seed(123)
  
  ssize      <- scen.mat[i,1]                                                   # save sample size
  short_imps <- vector("list", nsim)                                            # create list for imputations
  amatrixx   <- vector("list", nsim)                                            # create list for amatrices
  postimps   <- list(NA, nsim)                                                  # create list for posteriors
  
  for (s in 1:nsim){                                                            # for every iteration
    imputations  <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/4_impute/files_", 
                                      scen.mat[i,1],"_",                        # open the original imputations
                                      substring(scen.mat[i,2], 3),"_",
                                      substring(scen.mat[i,3], 3),"_",
                                      substring(scen.mat[i,4], 3),
                                      "/imputations_",s,".txt"), header=TRUE)
    short_imps[[s]] <-     ddply(imputations, .(q,z),nrow, .drop=FALSE)         # create an aggregated dataset
    
    for (b in 1:nboot) {                                                        # for every bootstrap sample
      amatrixx[[s]][[b]] <- round(read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/11_BCH/BCH_", 
                                                    scen.mat[i,1],"_",          # open the amatrices
                                                    substring(scen.mat[i,2], 3),"_",
                                                    substring(scen.mat[i,3], 3),"_",
                                                    substring(scen.mat[i,4], 3),
                                                    "/BCH_amatrix_",s,b,".txt"), header=TRUE), 8)
    }
    postimps[[s]] <- cbind(short_imps[[s]], matrix(NA, 4,5))                    # create space to save the posteriors
    colnames(postimps[[s]]) <- c("q","z","imp","post1","post2","post3","post4","post5")
    for (b in 1:nboot){
      postimps[[s]][1,b+3] <- amatrixx[[s]][[b]][1,1]/sum(amatrixx[[s]][[b]][1,]) # create posteriors
      postimps[[s]][2,b+3] <- amatrixx[[s]][[b]][2,1]/sum(amatrixx[[s]][[b]][2,])
      postimps[[s]][3,b+3] <- amatrixx[[s]][[b]][3,1]/sum(amatrixx[[s]][[b]][3,])
      postimps[[s]][4,b+3] <- amatrixx[[s]][[b]][4,1]/sum(amatrixx[[s]][[b]][4,])
    }
    write.table(postimps[[s]], paste0(paste0("F:/datapackage_threestepMILC/simulationstudy/13_impute_BCH/BCH_", 
                                             scen.mat[i,1],"_",                 # save dataset with posteriors
                                             substring(scen.mat[i,2], 3),"_",
                                             substring(scen.mat[i,3], 3),"_",
                                             substring(scen.mat[i,4], 3),
                                             "/newposteriors_",s,".txt")), sep="\t", 
                 row.names = FALSE, quote = FALSE)
  
    if(any(postimps[[s]]=="NaN")) {                                             # if there are any inadmissable
      new.imp <- matrix(NA, ssize, 13)                                          # solutions, create an empty matrix
    } else {                                                                    # otherwise
    new.matrix <- NULL                                                          # turn the short data into long data
    for (j in 1:4){
      old.matrix <- matrix(rep(postimps[[s]][j,],postimps[[s]][j,3]),          
                           postimps[[s]][j,3],8, byrow=T)
      new.matrix <- rbind(old.matrix, new.matrix)
    }
    new.result <- matrix(unlist(new.matrix), ncol=ncol(new.matrix), 
                         dimnames=list(NULL, colnames(new.matrix)))
    
    to.impute <- matrix(NA, nrow(new.result),5)                                 # create matrix to store imputations
    new.imp <- cbind(new.result, to.impute)
    for (l in 1:nrow(new.imp)){
      for (j in 1:nboot){
        new.imp[l,j+8] <- rbinom(1,1,new.imp[l,j+3])                            # create imputations
      }
    }
  }
    
  colnames(new.imp) <- c("q","z","n","post1","post2","post3","post4","post5",
                           "w1","w2","w3","w4","w5")
  write.table(new.imp, paste0(paste0("F:\\sim_kerst\\13_impute_BCH\\BCH_",      # save imputed dataset
                                       scen.mat[i,1],"_",
                                       substring(scen.mat[i,2], 3),"_",
                                       substring(scen.mat[i,3], 3),"_",
                                       substring(scen.mat[i,4], 3),
                                       "\\newimputations_",s,".txt")), sep="\t", 
              row.names = FALSE, quote = FALSE)
    
  }
}

