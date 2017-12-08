
scens <- list(
  ss     =c(200,500, 1000),
  c.prob =c(0.70, 0.80, 0.90, 0.95, 0.99),
  z.coef =c(0.01, 0.05, 0.10, 0.20),
  q.coef =c(0.5, 0.6224593, 0.7310586, 0.8807971))

scen.mat <- expand.grid(scens) 
list.out <- vector("list",nrow(scen.mat))

# model 1

for(a in 1:nrow(scen.mat)){
  list.out[[a]] <- read.table(paste0("F:\\sim_kerst\\5_results_before\\results_",
                                     scen.mat[a,1],"_",
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"\\finalresult_q.txt"))
}


results <- matrix(NA, nrow(scen.mat),8)

for (i in 1:240){
  results[i,1] <- scen.mat[i,1]
  results[i,2] <- scen.mat[i,2]
  results[i,3] <- scen.mat[i,3]
  results[i,4] <- scen.mat[i,4]
  results[i,5] <- as.matrix(list.out[[i]][2,7])
  results[i,6] <- as.matrix(list.out[[i]][2,8])
  results[i,7] <- as.matrix(list.out[[i]][2,9])
  results[i,8] <- as.matrix(list.out[[i]][2,10])
}

colnames(results) <- c("ss","c.prob","z.coef","q.coef","bias","cov","sesd","rmse")
results <- as.data.frame(results)

library(xtable)
t200_70 <- round(subset(results, ss==200&c.prob==0.70), 4)
t200_80 <- round(subset(results, ss==200&c.prob==0.80), 4)
t200_90 <- round(subset(results, ss==200&c.prob==0.90), 4)
t200_95 <- round(subset(results, ss==200&c.prob==0.95), 4)
t200_99 <- round(subset(results, ss==200&c.prob==0.99), 4)

t200 <- cbind(t200_70, t200_80[,c(5:8)], t200_90[,c(5:8)], t200_95[,c(5:8)], t200_99[,c(5:8)])

t500_70 <- round(subset(results, ss==500&c.prob==0.70), 4)
t500_80 <- round(subset(results, ss==500&c.prob==0.80), 4)
t500_90 <- round(subset(results, ss==500&c.prob==0.90), 4)
t500_95 <- round(subset(results, ss==500&c.prob==0.95), 4)
t500_99 <- round(subset(results, ss==500&c.prob==0.99), 4)

t500 <- cbind(t500_70, t500_80[,c(5:8)], t500_90[,c(5:8)], t500_95[,c(5:8)], t500_99[,c(5:8)])

t1000_70 <- round(subset(results, ss==1000&c.prob==0.70), 4)
t1000_80 <- round(subset(results, ss==1000&c.prob==0.80), 4)
t1000_90 <- round(subset(results, ss==1000&c.prob==0.90), 4)
t1000_95 <- round(subset(results, ss==1000&c.prob==0.95), 4)
t1000_99 <- round(subset(results, ss==1000&c.prob==0.99), 4)

t1000 <- cbind(t1000_70, t1000_80[,c(5:8)], t1000_90[,c(5:8)], t1000_95[,c(5:8)], t1000_99[,c(5:8)])

tableNO <- rbind(t200, t500, t1000)

rownames(tableNO) <- NULL
tableNO <- tableNO[,-2]
tableNO <- tableNO[,c(1,3,2,4:23)]

xtable(tableNO, digits=c(0,0,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))


#####################################################################################
# ML E
for(a in 1:nrow(scen.mat)){
  list.out[[a]] <- read.table(paste0("F:\\sim_kerst\\8_LC_results_ML1\\ML1_",
                                     scen.mat[a,1],"_",
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"\\finalresult_q.txt"))
}



results <- matrix(NA, nrow(scen.mat),8)

for (i in 1:240){
  results[i,1] <- scen.mat[i,1]
  results[i,2] <- scen.mat[i,2]
  results[i,3] <- scen.mat[i,3]
  results[i,4] <- scen.mat[i,4]
  results[i,5] <- as.matrix(list.out[[i]][2,7])
  results[i,6] <- as.matrix(list.out[[i]][2,8])
  results[i,7] <- as.matrix(list.out[[i]][2,9])
  results[i,8] <- as.matrix(list.out[[i]][2,10])
}


colnames(results) <- c("ss","c.prob","z.coef","q.coef","bias","cov","sesd","rmse")
results <- as.data.frame(results)

library(xtable)
t200_70 <- round(subset(results, ss==200&c.prob==0.70), 4)
t200_80 <- round(subset(results, ss==200&c.prob==0.80), 4)
t200_90 <- round(subset(results, ss==200&c.prob==0.90), 4)
t200_95 <- round(subset(results, ss==200&c.prob==0.95), 4)
t200_99 <- round(subset(results, ss==200&c.prob==0.99), 4)

t200 <- cbind(t200_70, t200_80[,c(5:8)], t200_90[,c(5:8)], t200_95[,c(5:8)], t200_99[,c(5:8)])

t500_70 <- round(subset(results, ss==500&c.prob==0.70), 4)
t500_80 <- round(subset(results, ss==500&c.prob==0.80), 4)
t500_90 <- round(subset(results, ss==500&c.prob==0.90), 4)
t500_95 <- round(subset(results, ss==500&c.prob==0.95), 4)
t500_99 <- round(subset(results, ss==500&c.prob==0.99), 4)

t500 <- cbind(t500_70, t500_80[,c(5:8)], t500_90[,c(5:8)], t500_95[,c(5:8)], t500_99[,c(5:8)])

t1000_70 <- round(subset(results, ss==1000&c.prob==0.70), 4)
t1000_80 <- round(subset(results, ss==1000&c.prob==0.80), 4)
t1000_90 <- round(subset(results, ss==1000&c.prob==0.90), 4)
t1000_95 <- round(subset(results, ss==1000&c.prob==0.95), 4)
t1000_99 <- round(subset(results, ss==1000&c.prob==0.99), 4)

t1000 <- cbind(t1000_70, t1000_80[,c(5:8)], t1000_90[,c(5:8)], t1000_95[,c(5:8)], t1000_99[,c(5:8)])

tableMLE <- rbind(t200, t500, t1000)

rownames(tableMLE) <- NULL

xtable(tableMLE, digits=c(0,0,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))


#####################################################################################

#####################################################################################
# ML E
for(a in 1:nrow(scen.mat)){
  list.out[[a]] <- read.table(paste0("F:\\sim_kerst\\10_LC_results_ML2\\ML2_results_",
                                     scen.mat[a,1],"_",
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"\\finalresult_q.txt"))
}



results <- matrix(NA, nrow(scen.mat),8)

for (i in 1:240){
  results[i,1] <- scen.mat[i,1]
  results[i,2] <- scen.mat[i,2]
  results[i,3] <- scen.mat[i,3]
  results[i,4] <- scen.mat[i,4]
  results[i,5] <- as.matrix(list.out[[i]][2,7])
  results[i,6] <- as.matrix(list.out[[i]][2,8])
  results[i,7] <- as.matrix(list.out[[i]][2,9])
  results[i,8] <- as.matrix(list.out[[i]][2,10])
}


colnames(results) <- c("ss","c.prob","z.coef","q.coef","bias","cov","sesd","rmse")
results <- as.data.frame(results)

library(xtable)
t200_70 <- round(subset(results, ss==200&c.prob==0.70), 4)
t200_80 <- round(subset(results, ss==200&c.prob==0.80), 4)
t200_90 <- round(subset(results, ss==200&c.prob==0.90), 4)
t200_95 <- round(subset(results, ss==200&c.prob==0.95), 4)
t200_99 <- round(subset(results, ss==200&c.prob==0.99), 4)

t200 <- cbind(t200_70, t200_80[,c(5:8)], t200_90[,c(5:8)], t200_95[,c(5:8)], t200_99[,c(5:8)])

t500_70 <- round(subset(results, ss==500&c.prob==0.70), 4)
t500_80 <- round(subset(results, ss==500&c.prob==0.80), 4)
t500_90 <- round(subset(results, ss==500&c.prob==0.90), 4)
t500_95 <- round(subset(results, ss==500&c.prob==0.95), 4)
t500_99 <- round(subset(results, ss==500&c.prob==0.99), 4)

t500 <- cbind(t500_70, t500_80[,c(5:8)], t500_90[,c(5:8)], t500_95[,c(5:8)], t500_99[,c(5:8)])

t1000_70 <- round(subset(results, ss==1000&c.prob==0.70), 4)
t1000_80 <- round(subset(results, ss==1000&c.prob==0.80), 4)
t1000_90 <- round(subset(results, ss==1000&c.prob==0.90), 4)
t1000_95 <- round(subset(results, ss==1000&c.prob==0.95), 4)
t1000_99 <- round(subset(results, ss==1000&c.prob==0.99), 4)

t1000 <- cbind(t1000_70, t1000_80[,c(5:8)], t1000_90[,c(5:8)], t1000_95[,c(5:8)], t1000_99[,c(5:8)])

tableMLI <- rbind(t200, t500, t1000)

rownames(tableMLI) <- NULL
tableMLI <- tableMLI[,-2]
tableMLI <- tableMLI[,c(1,3,2,4:23)]

xtable(tableMLI, digits=c(0,0,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))


#####################################################################################

# BCH
for(a in 1:nrow(scen.mat)){
  list.out[[a]] <- read.table(paste0("F:\\sim_kerst\\14_results_BCH\\BCH_",
                                     scen.mat[a,1],"_",
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"\\finalresult_q.txt"))
}



results <- matrix(NA, nrow(scen.mat),8)

for (i in 1:240){
  results[i,1] <- scen.mat[i,1]
  results[i,2] <- scen.mat[i,2]
  results[i,3] <- scen.mat[i,3]
  results[i,4] <- scen.mat[i,4]
  results[i,5] <- as.matrix(list.out[[i]][2,7])
  results[i,6] <- as.matrix(list.out[[i]][2,8])
  results[i,7] <- as.matrix(list.out[[i]][2,9])
  results[i,8] <- as.matrix(list.out[[i]][2,10])
}


colnames(results) <- c("ss","c.prob","z.coef","q.coef","bias","cov","sesd","rmse")
results <- as.data.frame(results)

library(xtable)
t200_70 <- round(subset(results, ss==200&c.prob==0.70), 4)
t200_80 <- round(subset(results, ss==200&c.prob==0.80), 4)
t200_90 <- round(subset(results, ss==200&c.prob==0.90), 4)
t200_95 <- round(subset(results, ss==200&c.prob==0.95), 4)
t200_99 <- round(subset(results, ss==200&c.prob==0.99), 4)

t200 <- cbind(t200_70, t200_80[,c(5:8)], t200_90[,c(5:8)], t200_95[,c(5:8)], t200_99[,c(5:8)])

t500_70 <- round(subset(results, ss==500&c.prob==0.70), 4)
t500_80 <- round(subset(results, ss==500&c.prob==0.80), 4)
t500_90 <- round(subset(results, ss==500&c.prob==0.90), 4)
t500_95 <- round(subset(results, ss==500&c.prob==0.95), 4)
t500_99 <- round(subset(results, ss==500&c.prob==0.99), 4)

t500 <- cbind(t500_70, t500_80[,c(5:8)], t500_90[,c(5:8)], t500_95[,c(5:8)], t500_99[,c(5:8)])

t1000_70 <- round(subset(results, ss==1000&c.prob==0.70), 4)
t1000_80 <- round(subset(results, ss==1000&c.prob==0.80), 4)
t1000_90 <- round(subset(results, ss==1000&c.prob==0.90), 4)
t1000_95 <- round(subset(results, ss==1000&c.prob==0.95), 4)
t1000_99 <- round(subset(results, ss==1000&c.prob==0.99), 4)

t1000 <- cbind(t1000_70, t1000_80[,c(5:8)], t1000_90[,c(5:8)], t1000_95[,c(5:8)], t1000_99[,c(5:8)])

tableBCH <- rbind(t200, t500, t1000)

rownames(tableBCH) <- NULL
tableBCH <- tableBCH[,-2]
tableBCH <- tableBCH[,c(1,3,2,4:23)]

xtable(tableBCH, digits=c(0,0,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))
