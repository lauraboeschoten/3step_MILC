setwd("F:/datapackage_threestepMILC/simulationstudy/8_LC_results_ML1")

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
  ssize <- scen.mat[i,1]                                                        # save sample size
  for (s in 1:nsim){
    imputations  <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/4_impute/files_",
                                      scen.mat[i,1],"_",                        # load imputed files
                                      substring(scen.mat[i,2], 3),"_",
                                      substring(scen.mat[i,3], 3),"_",
                                      substring(scen.mat[i,4], 3),
                                      "/imputations_",s,".txt"), header=T)
    
    ctables  <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/6_classification_tables/table_", 
                                  scen.mat[i,1],"_",                            # load classification tables
                                  substring(scen.mat[i,2], 3),"_",
                                  substring(scen.mat[i,3], 3),"_",
                                  substring(scen.mat[i,4], 3),
                                  "/class_mat_",s,".txt"), header=T)
    
    amatrix <- vector("list", 5)                                                # create list fot 5 amatrices
    
    for(b in 1:nboot){
      ematrix     <- matrix(NA, 4,2)                                            # create ematrices
      ematrix[1,] <- table(imputations$z, imputations[,10+b],imputations$q)[1,,1]
      ematrix[3,] <- table(imputations$z, imputations[,10+b],imputations$q)[1,,2]
      
      if (length(table(imputations$z, imputations[,10+b],imputations$q)[,,1])==4){
        ematrix[2,] <- table(imputations$z, imputations[,10+b],imputations$q)[2,,1]
        ematrix[4,] <- table(imputations$z, imputations[,10+b],imputations$q)[2,,2]
      } else {
        ematrix[2,] <- 0
        ematrix[4,] <- 0
      }

      ematrix <- ematrix/ssize
      dmatrix <- matrix(unlist(ctables[1,]),2, byrow=TRUE)                      # create dmatrices
      
      qpsolve <- function(e,d,iequal){
        nr <- nrow(e)
        nc <- ncol(e)
        ncel <- nr*nc
        evec <- as.vector(e)
        id <- diag(nr)
        p <- kronecker(t(d),id)
        dmat <- kronecker(d %*% t(d),id)
        dvec <- as.vector(evec %*% p)
        im <- diag(ncel)
        i1 <- iequal
        i2 <- setdiff(1:ncel,i1)
        index <- c(i1,i2)
        im2 <- im[index,]
        at <- rbind(rep(1,ncel),im2)
        amat <- t(at)
        bvec <- c(1,rep(0,ncel))
        meq <- 1 + length(iequal)
        res <- solve.QP(dmat,dvec,amat,bvec,meq)
        return(res)
      }
    
      amatrix[[b]] <- matrix(qpsolve(ematrix,dmatrix,c(2,4))$solution, ncol=2)  # create amatrices
      
      write.table(amatrix[[b]], paste0("F:/datapackage_threestepMILC/simulationstudy/11_BCH/BCH_", 
                                       scen.mat[i,1],"_",                       # save amatrices
                                       substring(scen.mat[i,2], 3),"_",
                                       substring(scen.mat[i,3], 3),"_",
                                       substring(scen.mat[i,4], 3),
                                       "/BCH_amatrix_",s,b,".txt"), sep="\t", col.names=TRUE, row.names=FALSE)
    }
  }
}