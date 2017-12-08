mean.results <- function(i) {
  
  folder.name <- paste0("F:/datapackage_threestepMILC/simulationstudy/10_LC_results_ML2/ML2_results_", 
                        scen.mat[i,1],"_",                                      # set folder name to store results
                        substring(scen.mat[i,2], 3),"_",
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3))
  setwd(folder.name)
  
  final_q <- array(NA, dim=c(nsim,10,2))                                        # create arrays to store results in
  final_z <- array(NA, dim=c(nsim,10,4))
  
  
  for(s in 1:nsim){
    results_q           <- read.table(paste0("results_q_",s,".txt"), header=T)  # open results for q from the iterations
    final_q[s,c(1:9),1] <- as.matrix(results_q[1,c(11:19)])                     # and store the relevant output
    final_q[s,c(1:9),2] <- as.matrix(results_q[2,c(11:19)])
    final_q[s,10,1]     <- final_q[s,7,1]^2
    final_q[s,10,2]     <- final_q[s,7,2]^2
    results_z           <- read.table(paste0("results_z_",s,".txt"), header=T)  # open results for z and store relevant
    final_z[s,c(1:9),1] <- as.matrix(results_z[1,c(11:19)])                     # output
    final_z[s,c(1:9),2] <- as.matrix(results_z[2,c(11:19)])
    final_z[s,c(1:9),3] <- as.matrix(results_z[3,c(11:19)])
    final_z[s,c(1:9),4] <- as.matrix(results_z[4,c(11:19)])
    final_z[s,10,1]     <- final_z[s,7,1]^2
    final_z[s,10,2]     <- final_z[s,7,1]^2
    final_z[s,10,3]     <- final_z[s,7,1]^2
    final_z[s,10,4]     <- final_z[s,7,1]^2
  }
 
  finalresult_q <- matrix(NA, 2,10)                                             # create matrices for final results
  colnames(finalresult_q) <- c("Qbar","Ubar","B","T","sqT","pop","bias","cov","sesd","RMSE")
  rownames(finalresult_q) <- c("intercept","coefficient")
  
  finalresult_z <- matrix(NA, 4,10)
  colnames(finalresult_z) <- c("Qbar","Ubar","B","T","sqT","pop","bias","cov","sesd","RMSE")
  rownames(finalresult_z) <- c("w1z1","w2z1","w1z2","w2z2")
  
  for(u in 1:2){
    finalresult_q[u,c(1:7)] <- colMeans(final_q[,c(1:7),u],na.rm = TRUE)        # take the mean over nsim results for q
    finalresult_q[u,8]      <- sum(as.vector(na.omit(final_q[,8,u])) < finalresult_q[u,"pop"] & # coverage
                               finalresult_q[u,"pop"] < as.vector(na.omit(final_q[,9,u])))/
                               (nsim-(sum(is.na(final_q[,1,2]))))
    finalresult_q[u,9]      <- finalresult_q[u,"sqT"]/sd(final_q[,1,u], na.rm = TRUE) #se/sd
    finalresult_q[u,10]     <- sqrt(mean(final_q[,10,u], na.rm = TRUE))         # rmse
  }
  
  for(u in 1:4){
    finalresult_z[u,c(1:7)] <- colMeans(final_z[,c(1:7),u],na.rm = TRUE)        # take the mean over nsim results for z
    finalresult_z[u,8]      <- sum(as.vector(na.omit(final_z[,8,u])) < finalresult_z[u,"pop"] & # coverage
                                  finalresult_z[u,"pop"] < as.vector(na.omit(final_z[,9,u])))/
                                  (nsim-(sum(is.na(final_z[,1,2]))))   
    finalresult_z[u,9]      <- finalresult_z[u,"sqT"]/sd(final_z[,1,u], na.rm = TRUE) #se/sd
    finalresult_z[u,10]     <- sqrt(mean(final_z[,10,u], na.rm = TRUE))         # rmse
  }
  
  write.table(finalresult_q, "finalresult_q.txt", sep="\t")                     # write tables with final results
  write.table(finalresult_z, "finalresult_z.txt", sep="\t")
  nrmis <- sum(is.na(final_q[,1,2]))                                            # sum the number of inadmissable solutions
  write.table(nrmis, "nrmis.txt", sep="\t")                                     # and save this 
}
  

  
  
  