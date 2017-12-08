create.results <- function(i) {

  for (s in 1:nsim){                                                            # for every iteration in the simulation
    step3.model <- vector("list", nboot)                                        # create an empty list for the 5 results
    results_q           <- matrix(NA, 2, 19)                                    # create empty matrix q for results
    colnames(results_q) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",
                             "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                             "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
  
    results_z           <- matrix(NA, 4, 19)                                    # and matrix z
    colnames(results_z) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",
                             "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                             "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
    
    for (c in 1:nboot){                                                         # for each bootstrap sample
      step3.model[[c]] <- paste0("F:/datapackage_threestepMILC/simulationstudy/7_step3_ML/step3ML_",
                                 scen.mat[i,1],"_",                             # obtain the LC output
                                 substring(scen.mat[i,2], 3),"_",
                                 substring(scen.mat[i,3], 3),"_",
                                 substring(scen.mat[i,4], 3),
                                 "\\model.",s,c,".lst")
      
      scan.mod <- readLines(step3.model[[c]])
      lin.reg  <- grep('Regression',scan.mod)                                   # define on which line 'Regression' is found
      
      results_q[1,c]   <- as.numeric(strsplit(scan.mod[lin.reg+3], "\t")[[1]][4]) # obtain logistic regression output
      results_q[2,c]   <- as.numeric(strsplit(scan.mod[lin.reg+7], "\t")[[1]][4])
      results_q[1,c+5] <- as.numeric(strsplit(scan.mod[lin.reg+3], "\t")[[1]][5])^2
      results_q[2,c+5] <- as.numeric(strsplit(scan.mod[lin.reg+7], "\t")[[1]][5])^2
  
      
      text <- scan.mod[grep("^Profile", scan.mod) + 7]                          # define on which line 'Profile' is found
      text <- substr(text, start = 5, nchar(text))                              # create the marginal table by selecting
      text <- gsub("\t\\.\t", "\t0\t", text)                                    # the relevant numbers from the LG output
      v    <- scan(text = text, quiet = TRUE)
      p    <- v[seq(1, 28, by = 2)]
      marginals <- p[1:6]
      names(marginals) <- paste0(rep(c("x", "lq", "lz"), each=2), 1:2)
      joint <- cbind(Pr = p[-(1:6)], expand.grid(lz = 1:2, lq = 1:2, x = 1:2)[, 3:1])
      tab   <- xtabs(Pr ~ x+lq+lz, data = joint)
      tab   <- aperm(tab, c(1,3,2))
      tabXZ <- apply(tab, 1:2, sum)                                             # put them in a table

      results_z[1,c] <- tabXZ[1,1]
      results_z[2,c] <- tabXZ[2,1]
      results_z[3,c] <- tabXZ[1,2]
      results_z[4,c] <- tabXZ[2,2]
      
    }
    
    if(any(results_q[2,1:5]>(11.51292)|results_q[2,1:5]<(-11.51292))) {         # select inadmissable solutions again
      
      results_q           <- matrix(NA, 2, 19)
      colnames(results_q) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",         # and create empty matrices for them
                               "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                               "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
      write.table(results_q, paste0("F:/datapackage_threestepMILC/simulationstudy/ML1_", 
                                    scen.mat[i,1],"_", substring(scen.mat[i,2], 3),"_",
                                    substring(scen.mat[i,3], 3),"_", substring(scen.mat[i,4], 3),
                                    "/results_q_",s,".txt"), sep="\t")
      
      results_z           <- matrix(NA, 4, 19)
      colnames(results_z) <- c("Qhat1","Qhat2","Qhat3","Qhat4","Qhat5",
                               "Uhat1","Uhat2","Uhat3","Uhat4","Uhat5",
                               "Qbar","Ubar","B","T","sqT","pop","bias","ll","ul")
      write.table(results_z, paste0("F:/datapackage_threestepMILC/simulationstudy/ML1_", 
                                    scen.mat[i,1],"_", substring(scen.mat[i,2], 3),"_",
                                    substring(scen.mat[i,3], 3),"_", substring(scen.mat[i,4], 3),
                                    "/results_z_",s,".txt"), sep="\t")
    } else {
      for (c in 1:2){                                                           # for the rest of the solutions we 
        results_q[c,"Qbar"] <- mean(results_q[c,1:5])                           # save the relevant results
        results_q[c,"Ubar"] <- mean(results_q[c,6:10])                          # and pool them
        results_q[c,"B"]    <- var(results_q[c,1:5])
        results_q[c,"T"]    <- results_q[c,"Ubar"]+(nboot+1)*(results_q[c,"B"]/nboot) 
      }
      results_q[,"sqT"]     <- sqrt(results_q[,"T"])                            # also save bias and confidence intervals
      results_q[,"pop"]     <- c(0,(log(scen.mat[i,4]/(1-scen.mat[i,4]))))      # for the q results
      results_q[,"bias"]    <- results_q[,"pop"] - results_q[,"Qbar"]
      results_q[,"ll"]      <- results_q[,"Qbar"]-(qt(.975,ssize-1)*results_q[,"sqT"])
      results_q[,"ul"]      <- results_q[,"Qbar"]+(qt(.975,ssize-1)*results_q[,"sqT"])
      
      for(k in 1:4){                                                            # also obtain the z results
        for(c in 1:nboot){
          results_z[k,c+5]  <- results_z[k,c]*(1-results_z[k,c])/ssize
        }
        results_z[k,"Qbar"] <- mean(results_z[k,1:5])                           # and pool them
        results_z[k,"Ubar"] <- mean(results_z[k,6:10])
        results_z[k,"B"]    <- var(results_z[k,1:5])
        results_z[k,"T"]    <- results_z[k,"Ubar"]+(5+1)*(results_z[k,"B"]/5)  
      }
      results_z[,"sqT"]  <- sqrt(results_z[,"T"])                               # obtain population values, bias and 
      results_z[,"pop"]  <- c(0.5, (0.5-scen.mat[i,3]), 0.0, scen.mat[i,3])     # confidence intervals
      results_z[,"bias"] <- results_z[,"pop"]-results_z[,"Qbar"]
      results_z[,"ll"]   <- results_z[,"Qbar"]-(qt(.975,ssize-1)*sqrt(results_z[,"T"])) 
      results_z[,"ul"]   <- results_z[,"Qbar"]+(qt(.975,ssize-1)*sqrt(results_z[,"T"]))
      
      write.table(results_q, paste0("F:/datapackage_threestepMILC/simulationstudy/8_LC_results_ML1/ML1_", 
                         scen.mat[i,1],"_", substring(scen.mat[i,2], 3),"_",    # save the results tables
                         substring(scen.mat[i,3], 3),"_", substring(scen.mat[i,4], 3),
                         "/results_q_",s,".txt"), sep="\t")
      
      write.table(results_z, paste0("F:/datapackage_threestepMILC/simulationstudy/8_LC_results_ML1\\ML1_", 
                         scen.mat[i,1],"_", substring(scen.mat[i,2], 3),"_",
                         substring(scen.mat[i,3], 3),"_", substring(scen.mat[i,4], 3),
                         "/results_z_",s,".txt"), sep="\t")
    }
  }
}
    
    
    