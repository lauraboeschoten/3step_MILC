impute.samples <- function(posterior, bootstrap.sample, ssize) {
  
  new.matrix <- NULL                                                            # create an empty matrix
  
  for (j in 1:ncomb){                                                           # for every combination of scores
    old.matrix <- matrix(rep(posterior[j,1:10],bootstrap.sample[j,6]),          # create a matrix for each combination
                         bootstrap.sample[j,6],10, byrow=T)                     # of scores and make it the length of 
    new.matrix <- rbind(old.matrix, new.matrix)                                 # the frequency corresponding to this 
  }                                                                             # combination of scores
  
  new.result <- matrix(unlist(new.matrix), ncol=ncol(new.matrix),               # put it in a matrix   
                       dimnames=list(NULL, colnames(new.matrix)))

  to.impute <- matrix(NA, ssize, nboot)                                         # create matrix for the imputations
  posterior.probabilities <- cbind(new.result, to.impute)                       # combine the new dataset with the empty
                                                                                # matrix
  for (k in 1:ssize){                                                           # loop over the length of the data
    for (m in 1:nboot){                                                         # and over the bootstrap samples
      posterior.probabilities[k,m+10] <- rbinom(1, 1, posterior.probabilities[k,m+5])                                                          
    }                                                                           # sample from the posteriors to create
  }                                                                             # imputations
  
  colnames(posterior.probabilities) <- c("q","z","y1","y2","y3",                # give names to this new dataset
                                         "post1","post2","post3","post4","post5",
                                         "w1","w2","w3","w4","w5")
  imputations <- posterior.probabilities                                        # save this new dataset under imputations
  
  return(imputations)
}


