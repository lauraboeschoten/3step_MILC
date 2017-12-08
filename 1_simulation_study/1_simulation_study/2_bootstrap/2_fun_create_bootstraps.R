create.bootstrap.samples <- function(nsim, nboot, ssize, ncomb, samples, i) {   
  
  # make list with separate sample files
  list_samples  <- split(samples, samples$simulation)                           # make a list containing the sample files for all iterations 
  short_samples <- list(NA)                                                     # create an empty list to put the aggregated samples in
  
  for (j in 1:nsim) {                                                           # change every sample from a long format into
    short_samples[[j]] <- ddply(list_samples[[j]],                              # the aggregated short format
                                .(q,z,y1,y2,y3),nrow, 
                                .drop=FALSE)
    
    if(nrow(short_samples[[j]])==16) {                                          # if a combination of scores does not 
      extra <- matrix(NA,16,6)                                                  # occur in the data, the number of combinations
      extra[,1] <- short_samples[[j]][,1]                                       # is different and we need to adjust this
      extra[,2] <- 2
      extra[,3] <- short_samples[[j]][,3]
      extra[,4] <- short_samples[[j]][,4]
      extra[,5] <- short_samples[[j]][,5]
      extra[,6] <- 0

      part1 <- as.matrix(short_samples[[j]][1:8,])
      part2 <- as.matrix(short_samples[[j]][9:16,])
      yes <- rbind(part1,extra[1:8,],part2,extra[9:16,])
      short_samples[[j]] <- as.data.frame(yes, row.names = FALSE)
    }
    
  }

  ncols                  <- length(short_samples[[1]])                          # here we create an array where we place
  boot_samples           <- array(NA, dim=c(ncomb, ncols+nboot, nsim))          # all the bootstrap samples we are going
  colnames(boot_samples) <- c(names(short_samples[[1]])[1:(ncols-1)],"original",# to create
                              paste0("b",seq(1:nboot)))    
  
  # each layer in array gets a sample dataset and bootstrap samples are taken
  for (j in 1:nsim){                                                            # now each layer in the array gets a sample
    boot_samples[,1:ncols,j] <- as.matrix(short_samples[[j]])                   # of an iteration and the corresponding 
    boot_samples[,((ncols+1):(ncols+nboot)),j] <- rmultinom(nboot, ssize,       # bootstrap samples
                                                        short_samples[[j]][,ncols]/
                                                        sum(short_samples[[j]][,ncols]))           
  }                                                                             
  
  
  for(j in 1:nsim){                                                             # we save a set of bootstrap samples for 
    write.table(boot_samples[,,j],                                              # each iteration of the simulation
                file=paste0("F:/datapackage_threestepMILC/simulationstudy/2_bootstrap/samples_", 
                            scen.mat[i,1],"_",                                  # in the folder corresponding to the 
                            substring(scen.mat[i,2], 3),"_",                    # simulation condition
                            substring(scen.mat[i,3], 3),"_",                    # in a separate .txt file
                            substring(scen.mat[i,4], 3),"/bootstrap_sample_",j,".txt"), 
                row.names=FALSE, quote=FALSE)
  }
}
