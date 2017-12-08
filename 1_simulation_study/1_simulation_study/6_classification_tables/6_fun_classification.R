classification <- function(i){
  
  set.seed(123)
  cat(i)
  class_table <- list(NA)                                                       # create an empty list to store
                                                                                # the classification tables
  for (s in 1:nsim){                                                            # for every iteration in the simulation
    imputed.sample <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/4_impute/files_", 
                                        scen.mat[i,1],"_",                      # obtain the imputed dataset
                                        substring(scen.mat[i,2], 3),"_",
                                        substring(scen.mat[i,3], 3),"_",
                                        substring(scen.mat[i,4], 3),
                                        "/imputations_",s,".txt"), header=TRUE)
    
    short_impsam      <- list(NA)                                               # create a list to store the intermediate
    short_impsam[[1]] <- ddply(imputed.sample, .(y1,y2,y3,post1), summarise,    # results
                               W1 = sum(w1), W2 = sum(1-w1), N  = (W1+W2))      # create an aggregated dataset with the
    short_impsam[[2]] <- ddply(imputed.sample, .(y1,y2,y3,post2), summarise,    # number of assignments to each class for
                               W1 = sum(w2), W2 = sum(1-w2), N  = (W1+W2))      # every combination of scores
    short_impsam[[3]] <- ddply(imputed.sample, .(y1,y2,y3,post3), summarise,    # and do this for every bootstrap sample
                               W1 = sum(w3), W2 = sum(1-w3), N  = (W1+W2))      # also save the total number of assigned 
    short_impsam[[4]] <- ddply(imputed.sample, .(y1,y2,y3,post4), summarise,    # scores per combination of scores
                               W1 = sum(w4), W2 = sum(1-w4), N  = (W1+W2))
    short_impsam[[5]] <- ddply(imputed.sample, .(y1,y2,y3,post5), summarise,
                               W1 = sum(w5), W2 = sum(1-w5), N  = (W1+W2))
    
    # obtain relevant information for classification tables
    for (j in 1:5){                                                             # for every bootstrap
      for (k in 1:nrow(short_impsam[[j]]))                                      # and every combination of scores
      {
        short_impsam[[j]][k,"PY"] <- short_impsam[[j]][k,"N"]/                  # percentage of persons with response pattern Y
          sum(short_impsam[[j]][ ,"N"])
        short_impsam[[j]][k,"PW1"] <- short_impsam[[j]][k,"W1"]/                # probability of assigning W1 given response
          sum(short_impsam[[j]][k,c("W1","W2")])                                # pattern Y and probability of assigning W2
        short_impsam[[j]][k,"PW2"] <- short_impsam[[j]][k,"W2"]/                # given response pattern Y
          sum(short_impsam[[j]][k,c("W1","W2")])
        short_impsam[[j]][k,"W1X1"] <- (short_impsam[[j]][k,"PY"]*              # probability if assigning W given X 
          short_impsam[[j]][k,paste0("post",j)]*short_impsam[[j]][k,"PW1"])
        short_impsam[[j]][k,"W1X2"] <- (short_impsam[[j]][k,"PY"]*
          (1-short_impsam[[j]][k,paste0("post",j)])*short_impsam[[j]][k,"PW1"])
        short_impsam[[j]][k,"W2X1"] <- (short_impsam[[j]][k,"PY"]*
          short_impsam[[j]][k,paste0("post",j)]*short_impsam[[j]][k,"PW2"])
        short_impsam[[j]][k,"W2X2"] <- (short_impsam[[j]][k,"PY"]*
          (1-short_impsam[[j]][k,paste0("post",j)])*short_impsam[[j]][k,"PW2"])
        short_impsam[[j]][k,"PX1"] <- short_impsam[[j]][k,paste0("post",j)]*    # probability of X
          (short_impsam[[j]][k,"N"]/sum(short_impsam[[j]][,"N"]))
        short_impsam[[j]][k,"PX2"] <- (1-short_impsam[[j]][k,paste0("post",j)])*
          (short_impsam[[j]][k,"N"]/sum(short_impsam[[j]][,"N"]))
      }
    }
    
    class_mat <- matrix(NA,5,4)                                                 # create empty matrix for the 5 
                                                                                # classification tables
    for (k in 1:5){
      class_mat[k,] <- c(sum(short_impsam[[k]][,"W1X1"])/                       # fill the matrix with the 5
                           sum(short_impsam[[k]][,"PX1"]),                      # classification tables
                         sum(short_impsam[[k]][,"W2X1"])/
                           sum(short_impsam[[k]][,"PX1"]),
                         sum(short_impsam[[k]][,"W1X2"])/
                           sum(short_impsam[[k]][,"PX2"]),
                         sum(short_impsam[[k]][,"W2X2"])/
                           sum(short_impsam[[k]][,"PX2"]))
    }
  
    write.table(class_mat, paste0("F:/datapackage_threestepMILC/simulationstudy/6_classification_tables/table_", 
                                  scen.mat[i,1],"_",                            # save a classification file for 
                                  substring(scen.mat[i,2], 3),"_",              # every iteration in the simulation
                                  substring(scen.mat[i,3], 3),"_",
                                  substring(scen.mat[i,4], 3),
                                  "/class_mat_",s,".txt"), sep="\t", col.names=TRUE, row.names=FALSE)
  }
}

