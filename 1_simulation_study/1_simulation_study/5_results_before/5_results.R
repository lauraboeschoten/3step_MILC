create.results <- function(i) {
  
  folder.name <- paste0("F:/datapackage_threestepMILC/simulationstudy/5_results_before/results_", 
                        scen.mat[i,1],"_",                                      # set to folder to store results
                        substring(scen.mat[i,2], 3),"_",
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3))
  
  setwd(folder.name)
  
  ssize <- scen.mat[i,1]                                                        # save sample size
  
  for (s in 1:nsim){
    data <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/4_impute/files_", 
                              scen.mat[i,1],"_",                                # open the imputed dataset
                              substring(scen.mat[i,2], 3),"_",
                              substring(scen.mat[i,3], 3),"_",
                              substring(scen.mat[i,4], 3),
                              "/imputations_",s,".txt"), header=TRUE)
    
    septables <- matrix(NA, 5,4)                                                # if statement to check for separation
    septables[1,] <- table(data$q, data$w1)
    septables[2,] <- table(data$q, data$w2)
    septables[3,] <- table(data$q, data$w3)
    septables[4,] <- table(data$q, data$w4)
    septables[5,] <- table(data$q, data$w5)
    
    if(any(septables==0)){
      
      results_q           <- matrix(NA, 2, 19)                                  # results saved if there is separation
      colnames(results_q) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",         # by saving just an empty matrix
                               "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                               "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
      write.table(results_q, paste0("results_q_",s,".txt"), sep="\t")
      
      results_z           <- matrix(NA, 4, 19)
      colnames(results_z) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",
                               "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                               "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
      
      write.table(results_z, paste0("results_z_",s,".txt"), sep="\t")
      
    } else {                                                                    # for all other solutions
      
      results_q           <- matrix(NA, 2, 19)                                  # create matrix to store results
      colnames(results_q) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",
                               "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                               "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
    
      data$q <- ifelse(data$q==2,1,0)                                           # recode the covariate from 1-2 in 1-0
    
      results_q[,1]  <- summary(glm(q ~ w1, "binomial", data))$coefficients[,1] # apply logistic regressions
      results_q[,2]  <- summary(glm(q ~ w2, "binomial", data))$coefficients[,1] # and save the coefficients
      results_q[,3]  <- summary(glm(q ~ w3, "binomial", data))$coefficients[,1] # and standard errors
      results_q[,4]  <- summary(glm(q ~ w4, "binomial", data))$coefficients[,1]
      results_q[,5]  <- summary(glm(q ~ w5, "binomial", data))$coefficients[,1]
    
      results_q[,6]  <- summary(glm(q ~ w1, "binomial", data))$coefficients[,2]^2
      results_q[,7]  <- summary(glm(q ~ w2, "binomial", data))$coefficients[,2]^2
      results_q[,8]  <- summary(glm(q ~ w3, "binomial", data))$coefficients[,2]^2
      results_q[,9]  <- summary(glm(q ~ w4, "binomial", data))$coefficients[,2]^2
      results_q[,10] <- summary(glm(q ~ w5, "binomial", data))$coefficients[,2]^2
    
    
      for (c in 1:2){                                                           # now mean over all obtained results
        
        results_q[c,"Qbar"] <- mean(results_q[c,1:5])
        results_q[c,"Ubar"] <- mean(results_q[c,6:10])
        results_q[c,"B"]    <- var(results_q[c,1:5])                            # and obtain between and  total variance
        results_q[c,"T"]    <- results_q[c,"Ubar"]+(5+1)*(results_q[c,"B"]/5) 
        
      }
      
      results_q[,"sqT"]  <- sqrt(results_q[,"T"])                               # sqrt of total variance
      results_q[,"pop"]  <- c(0,(log(scen.mat[i,4]/(1-scen.mat[i,4]))))         # population value
      results_q[,"bias"] <- results_q[,"pop"] - results_q[,"Qbar"]              # bias and confidence intervals
      results_q[,"ll"]   <- results_q[,"Qbar"]-(qt(.975,ssize-1)*results_q[,"sqT"]) 
      results_q[,"ul"]   <- results_q[,"Qbar"]+(qt(.975,ssize-1)*results_q[,"sqT"])
    
      write.table(results_q, paste0("results_q_",s,".txt"), sep="\t")           # write the results for covariate q
    
      results_z           <- matrix(NA, 4, 19)                                  # create matrix for the crosstable with z
      colnames(results_z) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",
                               "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                               "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
    
      for (c in 1:nboot) {
        results_z[,c]   <- (table(data[,c+10], data[,"z"]))/ssize               # obtain proportions and se's for the 4 cells
        results_z[,c+5] <- (results_z[,c]*(1-results_z[,c]))/ssize
      }
      
      for (c in 1:4) {
        results_z[c,"Qbar"] <- mean(results_z[c,1:5])                           # pool the proportions
        results_z[c,"Ubar"] <- mean(results_z[c,6:10])
        results_z[c,"B"]    <- var(results_z[c,1:5])
        results_z[c,"T"]    <- results_z[c,"Ubar"]+(5+1)*(results_z[c,"B"]/5)  
      }
      
      results_z[,"sqT"]  <- sqrt(results_z[,"T"])                               # obtain bias etc
      results_z[,"pop"]  <- c(0.5, (0.5-scen.mat[i,3]), 0.0, scen.mat[i,3])
      results_z[,"bias"] <- results_z[,"pop"]-results_z[,"Qbar"]
      results_z[,"ll"]   <- results_z[,"Qbar"]-(qt(.975,ssize-1)*results_z[,"sqT"]) 
      results_z[,"ul"]   <- results_z[,"Qbar"]+(qt(.975,ssize-1)*results_z[,"sqT"])
    
      write.table(results_z, paste0("results_z_",s,".txt"), sep="\t")           # write results for crosstable z
      
    }
  }
}
  