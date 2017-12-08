library(ggplot2)                                                                # load ggplot package
options(scipen=999)                                                             # scientific notation off

scens <- list(                                                                  # scenarios we show in the graph
  ss     =c(1000),
  c.prob =c(0.70, 0.80, 0.90, 0.95, 0.99),
  z.coef =c(0.20),
  q.coef =c(0.5, 0.6224593, 0.7310586, 0.8807971))

scen.mat <- expand.grid(scens) 
list.out <- vector("list",nrow(scen.mat))                                       # list for output

# step 1 model
for(a in 1:nrow(scen.mat)){
  ssize <- scen.mat[a,1]
  list.out[[a]] <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/5_results_before/results_",
                                     scen.mat[a,1],"_",                         # load the output for 'no correction'
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"/finalresult_z.txt"))
}

prop1 <- matrix(NA, length(list.out), 1)                                        # save the conditions and the output
for(i in 1:length(list.out)){
  prop1[i,1] <- list.out[[i]][3,1]*ssize
}

# ML 1
for(a in 1:nrow(scen.mat)){
  list.out[[a]] <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/8_LC_results_ML1/ML1_",
                                     scen.mat[a,1],"_",                         # load the output for MLE
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"/finalresult_z.txt"))
}

prop2 <- matrix(NA, length(list.out), 1)                                        # save the conditions and the output
for(i in 1:length(list.out)){
  prop2[i,1] <- list.out[[i]][3,1]*ssize
}

# ML 3
for(a in 1:nrow(scen.mat)){
  list.out[[a]] <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/10_LC_results_ML2/ML2_results_",
                                     scen.mat[a,1],"_",                         # load the output for  MLI
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"/finalresult_z.txt"))
}

prop3 <- matrix(NA, length(list.out), 1)                                        # save the conditions and the output
for(i in 1:length(list.out)){
  prop3[i,1] <- list.out[[i]][4,1]*ssize
}

## BCH
for(a in 1:nrow(scen.mat)){
  list.out[[a]] <- read.table(paste0("F:/datapackage_threestepMILC/simulationstudy/14_results_BCH/BCH_",
                                     scen.mat[a,1],"_",                         # load the output for BCH
                                     substring(scen.mat[a,2], 3),"_",
                                     substring(scen.mat[a,3], 3),"_",
                                     substring(scen.mat[a,4], 3),"/finalresult_z.txt"))
}

prop4 <- matrix(NA, length(list.out), 1)                                        # save the conditions and the output
for(i in 1:length(list.out)){
  prop4[i,1] <- list.out[[i]][4,1]*ssize
}

models       <- rep(c("no correction", "ML-E","ML-I","BCH"), each=length(list.out)) # create a vector with condition names
scens4       <- rbind(scen.mat, scen.mat, scen.mat, scen.mat)                   # combine the results
props        <- rbind(prop1, prop2, prop3, prop4)
tabel        <- cbind(models, scens4, props)
tabel$q.coef <- factor(tabel$q.coef, labels = c("logit = 0","logit = 0.5","logit = 1","logit = 2"))
tabel$coml   <- rep(1:5, 16)
tabel$models <- factor(tabel$models, as.character(tabel$models))
tabel$imp    <- factor(tabel$imp, labels = c("P(Z=2)=0.01","P(Z=2)=0.05","P(Z=2)=0.10","P(Z=2)=0.20"))
tabel$ss     <- factor(tabel$ss, labels = c("sample size 200"))
tabel$c.prob <- factor(tabel$c.prob, labels = c("0.70","0.80","0.90","0.95","0.99"))
tabel$q.coef <- factor(tabel$q.coef, labels = c("0","0.5","1","2"))

p <- ggplot(tabel, aes(x=c.prob, y=props))                                      # make a ggplot
p + geom_bar(stat="identity", position = "dodge") +
  facet_wrap(~models, ncol=4) +
  theme_bw() +
  ylab("number of impossible combinations created after imputation") +
  xlab("conditional probability") +
  geom_hline(yintercept = 0)
