impute <- function(i){
  
  for(s in 1:nsim){                                                             # for every iteration in the simulation
    for(c in 1:nboot){                                                          # for every bootstrap sample
      data <- read.delim(paste0("F:/datapackage_threestepMILC/simulationstudy/7_step3_ML/step3ML_", 
                                scen.mat[i,1],"_",                              # read in the datasets containing
                                substring(scen.mat[i,2], 3),"_",                # the posterior membership probabilities
                                substring(scen.mat[i,3], 3),"_",                # obtained after applying ML 3 step
                                substring(scen.mat[i,4], 3),"/LC3step",s,c,".dat"), 
                         header=T, sep="\t")
    
    comdata           <- matrix(NA, ssize, 5)                                   # create empty matrix to store results
    colnames(comdata) <- c("q","z","wold","post","wnew")
    comdata[,"q"]     <- data[,"q"]                                             # put the covariates and old results
    comdata[,"z"]     <- data[,"z"]                                             # in the matrix 
    comdata[,"wold"]  <- data[,1]                                               # and the relevant posteriors
    comdata[,"post"]  <- data[,4]

    for(k in 1:ssize){                                                          # sample from the posteriors to create
      comdata[k,5] <- rbinom(1, 1, comdata[k,"post"])                           # new imputations
    }

    write.table(round(comdata, 4), paste0("F:/datapackage_threestepMILC/simulationstudy/9_ML_impute/MLimpute_", 
                                          scen.mat[i,1],"_",                    # write table with new imputations
                                          substring(scen.mat[i,2], 3),"_",
                                          substring(scen.mat[i,3], 3),"_",
                                          substring(scen.mat[i,4], 3),"/imputations_",s,c,".txt"), 
                sep="\t", row.names = FALSE, quote = FALSE)
    }
  }
}