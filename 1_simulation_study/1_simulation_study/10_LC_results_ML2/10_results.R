create.results <- function(i) {

  folder.name <- paste0("F:/datapackage_threestepMILC/simulationstudy/10_LC_results_ML2/ML2_results_", 
                        scen.mat[i,1],"_",                                      # specify folder with results
                        substring(scen.mat[i,2], 3),"_",
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3))
  
  setwd(folder.name)
  
  for (s in 1:nsim){
    data1 <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/9_ML_impute/MLimpute_", 
                              scen.mat[i,1],"_",                                # obtain imputed datasets
                              substring(scen.mat[i,2], 3),"_",
                              substring(scen.mat[i,3], 3),"_",
                              substring(scen.mat[i,4], 3),
                              "/imputations_",s,1,".txt"), header=TRUE)
    data2 <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/9_ML_impute/MLimpute_", 
                               scen.mat[i,1],"_",
                               substring(scen.mat[i,2], 3),"_",
                               substring(scen.mat[i,3], 3),"_",
                               substring(scen.mat[i,4], 3),
                               "/imputations_",s,2,".txt"), header=TRUE)
    data3 <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/MLimpute_", 
                               scen.mat[i,1],"_",
                               substring(scen.mat[i,2], 3),"_",
                               substring(scen.mat[i,3], 3),"_",
                               substring(scen.mat[i,4], 3),
                               "/imputations_",s,3,".txt"), header=TRUE)
    data4 <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/MLimpute_", 
                               scen.mat[i,1],"_",
                               substring(scen.mat[i,2], 3),"_",
                               substring(scen.mat[i,3], 3),"_",
                               substring(scen.mat[i,4], 3),
                               "/imputations_",s,4,".txt"), header=TRUE)
    data5 <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/9_ML_impute/MLimpute_", 
                               scen.mat[i,1],"_",
                               substring(scen.mat[i,2], 3),"_",
                               substring(scen.mat[i,3], 3),"_",
                               substring(scen.mat[i,4], 3),
                               "/imputations_",s,5,".txt"), header=TRUE)

    septables <- matrix(NA, 5,4)                                                # check if there is any separation
    septables[1,] <- table(data1$q, data1$wnew)
    septables[2,] <- table(data2$q, data2$wnew)
    septables[3,] <- table(data3$q, data3$wnew)
    septables[4,] <- table(data4$q, data4$wnew)
    septables[5,] <- table(data5$q, data5$wnew)
    
    if(any(septables==0)){                                                      # create empty matriced for the 
                                                                                # iterations that contained separation
      results_q           <- matrix(NA, 2, 19)
      colnames(results_q) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",
                               "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                               "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
      write.table(results_q, paste0("results_q_",s,".txt"), sep="\t")
      
      results_z <- matrix(NA, 4, 19)
      colnames(results_z) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",
                               "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                               "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
      write.table(results_z, paste0("results_z_",s,".txt"), sep="\t")
      
    } else {                                                                    # otherwise create empty matrix
                                                                                # to store results in
      results_q           <- matrix(NA, 2, 19)
      colnames(results_q) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",
                               "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                               "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
      data1$q <- ifelse(data1$q==2,0,1)                                         # recode q
      data2$q <- ifelse(data2$q==2,0,1)
      data3$q <- ifelse(data3$q==2,0,1)
      data4$q <- ifelse(data4$q==2,0,1)
      data5$q <- ifelse(data5$q==2,0,1)
    
      results_q[,1]  <- summary(glm(q ~ wnew, "binomial", data1))$coefficients[,1] # apply logistic regression
      results_q[,2]  <- summary(glm(q ~ wnew, "binomial", data2))$coefficients[,1]
      results_q[,3]  <- summary(glm(q ~ wnew, "binomial", data3))$coefficients[,1]
      results_q[,4]  <- summary(glm(q ~ wnew, "binomial", data4))$coefficients[,1]
      results_q[,5]  <- summary(glm(q ~ wnew, "binomial", data5))$coefficients[,1]
    
      results_q[,6]  <- summary(glm(q ~ wnew, "binomial", data1))$coefficients[,2]^2
      results_q[,7]  <- summary(glm(q ~ wnew, "binomial", data2))$coefficients[,2]^2
      results_q[,8]  <- summary(glm(q ~ wnew, "binomial", data3))$coefficients[,2]^2
      results_q[,9]  <- summary(glm(q ~ wnew, "binomial", data4))$coefficients[,2]^2
      results_q[,10] <- summary(glm(q ~ wnew, "binomial", data5))$coefficients[,2]^2

      for (c in 1:2){                                                           # pool results for q
        results_q[c,"Qbar"] <- mean(results_q[c,1:5])
        results_q[c,"Ubar"] <- mean(results_q[c,6:10])
        results_q[c,"B"]    <- var(results_q[c,1:5])
        results_q[c,"T"]    <- results_q[c,"Ubar"]+(5+1)*(results_q[c,"B"]/5) 
      }
      results_q[,"sqT"]  <- sqrt(results_q[,"T"])                               # add population values, bias, ci
      results_q[,"pop"]  <- c(0,(log(scen.mat[i,4]/(1-scen.mat[i,4]))))
      results_q[,"bias"] <- results_q[,"pop"] - results_q[,"Qbar"]
      results_q[,"ll"]   <- results_q[,"Qbar"]-(qt(.975,ssize-1)*results_q[,"sqT"]) 
      results_q[,"ul"]   <- results_q[,"Qbar"]+(qt(.975,ssize-1)*results_q[,"sqT"])
    
      write.table(results_q, paste0("results_q_",s,".txt"), sep="\t")           # write results for every iteration in file
    
      results_z <- matrix(NA, 4, 19)                                            # create empty matrix for z
      colnames(results_z) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",
                               "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                               "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
    
      results_z[,1]   <- (table(data1[,"wnew"], data1[,"z"]))/ssize             # obtain cell proportions and se's
      results_z[,2]   <- (table(data2[,"wnew"], data2[,"z"]))/ssize
      results_z[,3]   <- (table(data3[,"wnew"], data3[,"z"]))/ssize
      results_z[,4]   <- (table(data4[,"wnew"], data4[,"z"]))/ssize
      results_z[,5]   <- (table(data5[,"wnew"], data5[,"z"]))/ssize
      
      results_z[,6]  <- (results_z[,1]*(1-results_z[,1]))/ssize
      results_z[,7]  <- (results_z[,2]*(1-results_z[,2]))/ssize
      results_z[,8]  <- (results_z[,3]*(1-results_z[,3]))/ssize
      results_z[,9]  <- (results_z[,4]*(1-results_z[,4]))/ssize
      results_z[,10] <- (results_z[,5]*(1-results_z[,5]))/ssize

      for (c in 1:4){                                                           # pool results for z
        results_z[c,"Qbar"] <- mean(results_z[c,1:5])
        results_z[c,"Ubar"] <- mean(results_z[c,6:10])
        results_z[c,"B"]    <- var(results_z[c,1:5])
        results_z[c,"T"]    <- results_z[c,"Ubar"]+(5+1)*(results_z[c,"B"]/5)  
      }
      results_z[,"sqT"]  <- sqrt(results_z[,"T"])                               # obtain population values, bias and ci
      results_z[,"pop"]  <- c((0.5-scen.mat[i,3]), 0.5, scen.mat[i,3], 0.0)
      results_z[,"bias"] <- results_z[,"pop"]-results_z[,"Qbar"]
      results_z[,"ll"]   <- results_z[,"Qbar"]-(qt(.975,ssize-1)*results_z[,"sqT"]) 
      results_z[,"ul"]   <- results_z[,"Qbar"]+(qt(.975,ssize-1)*results_z[,"sqT"])
    
      write.table(results_z, paste0("results_z_",s,".txt"), sep="\t")           # write table with results
    }
  }
}