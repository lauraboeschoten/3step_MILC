library(ggplot2)                                                                # load ggplot package
options(scipen=999)                                                             # scientific notation off

scens <- list(                                                                  # scenarios we show in the graph
  ss     =c(1000),
  c.prob =c(0.70, 0.80, 0.90, 0.95, 0.99),
  z.coef =c(0.20),
  q.coef =c(0.5, 0.6224593, 0.7310586, 0.8807971))

scen.mat <- expand.grid(scens) 
list.out <- vector("list",nrow(scen.mat))                                       # list for output

# model 1
for(a in 1:nrow(scen.mat)){
  list.out[[a]] <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/5_results_before/results_",
                                     scen.mat[a,1],"_",                         # load the output for 'no correction'
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"/finalresult_q.txt"))
}

coef <- cbind(scen.mat, cov=rep(NA,length(list.out)))                           # save the conditions and the output
for(i in 1:length(list.out)){
  coef[i,5] <- list.out[[i]][2,8]  
}


# model step 3 MLE
for(a in 1:nrow(scen.mat)){
  list.out[[a]] <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/8_LC_results_ML1/ML1_",
                                     scen.mat[a,1],"_",                         # load the output for 'MLE'
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"/finalresult_q.txt"))
}

coef1 <- cbind(scen.mat, cov=rep(NA,length(list.out)))                          # save the conditions and the output
for(i in 1:length(list.out)){
  coef1[i,5] <- list.out[[i]][2,8]  
}

# MODEL MLI
for(a in 1:nrow(scen.mat)){
  list.out[[a]] <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/10_LC_results_ML2/ML2_results_",
                                     scen.mat[a,1],"_",                         # load the output for 'MLI'
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"/finalresult_q.txt"))
}

coef3 <- cbind(scen.mat, cov=rep(NA,length(list.out)))                          # save the conditions and the output
for(i in 1:length(list.out)){
  coef3[i,5] <- list.out[[i]][2,8]  
}

# BCH 
for(a in 1:nrow(scen.mat)){
  list.out[[a]] <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/14_results_BCH/BCH_",
                                     scen.mat[a,1],"_",                         # load the output for 'BCH'
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"/finalresult_q.txt"))
}

coefBCH <- cbind(scen.mat, cov=rep(NA,length(list.out)))                        
for(i in 1:length(list.out)){                                                   # save the conditions and the output
  coefBCH[i,5] <- list.out[[i]][2,8]  
}

######

models        <- rep(c("no correction", "ML-E","ML-I","BCH"), each=length(list.out)) # create a vector with condition names
tabel         <- rbind(coef,coef1,coef3,coefBCH)                                # combine the results
tabel         <- cbind(tabel, models)
tabel$q.coef  <- factor(tabel$q.coef, labels = c("logit = 0","logit = 0.5","logit = 1","logit = 2")) # q conditions
tabel$coml    <- rep(1:5, 16)                                                   # number the outputs

ggplot(tabel, aes(x=coml, y=cov, group=models)) +                               # create a ggplot
  geom_line(aes(linetype=models), size=1) + 
  scale_linetype_manual(values = c("solid","dashed","dotted","dotdash")) +
  geom_point(aes(shape=models), size=4, fill="grey") + 
  scale_shape_manual(values=c(21,22,23,24)) + 
  facet_wrap(~q.coef, ncol=4) +
  theme_bw() +
  geom_hline(yintercept = 0.95) +
  scale_x_continuous(breaks=c(1,2,3,4,5),
                     labels=c("0.70", "0.80", "0.90","0.95","0.99")) +
  guides(colour= guide_legend("method")) +
  ylab("coverage of the 95 percent confidence interval") +
  xlab("conditional probability")
