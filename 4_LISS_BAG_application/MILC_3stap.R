set.seed(123)                                                                   
                                                                                
library(brew)
library(plyr)
library(quadprog)

setwd("F:/datapackage_threestepMILC/4_LISS_BAG_application")

source("runtemplate.R")                                                         

samples <- read.table("1_bootstrap_samples.txt", comment.char="", header=TRUE)  # data including bootstrap samples

################################################################################
                                                                                
outfiles           <- paste0("3_model.",1:5)                                    # run 5 LC models without covariate
bootstrap.colnames <- tail(colnames(samples),5)

envir <- new.env()                                                             
run.template(template.path="2_LC1.brew", envir=envir,temp.filename.base=outfiles[1],                                    
             bootstrap.colnames=bootstrap.colnames[1])                           
                                                                                
envir <- new.env()
run.template(template.path="2_LC2.brew", envir=envir,temp.filename.base=outfiles[2], 
             bootstrap.colnames=bootstrap.colnames[2])

envir <- new.env()
run.template(template.path="2_LC3.brew", envir=envir,temp.filename.base=outfiles[3], 
             bootstrap.colnames=bootstrap.colnames[3])

envir <- new.env()
run.template(template.path="2_LC4.brew", envir=envir,temp.filename.base=outfiles[4], 
             bootstrap.colnames=bootstrap.colnames[4])

envir <- new.env()
run.template(template.path="2_LC5.brew", envir=envir,temp.filename.base=outfiles[5], 
             bootstrap.colnames=bootstrap.colnames[5])

################################################################################
                                                                                
lc1 <- read.table("2_LC1.txt", comment.char="", header=TRUE)                    # obtain posteriors           
lc2 <- read.table("2_LC2.txt", comment.char="", header=TRUE)                    
lc3 <- read.table("2_LC3.txt", comment.char="", header=TRUE)                  
lc4 <- read.table("2_LC4.txt", comment.char="", header=TRUE)
lc5 <- read.table("2_LC5.txt", comment.char="", header=TRUE)

samples$LC1 <- lc1[,"Cluster.1"]
samples$LC2 <- lc2[,"Cluster.1"]
samples$LC3 <- lc3[,"Cluster.1"]
samples$LC4 <- lc4[,"Cluster.1"]
samples$LC5 <- lc5[,"Cluster.1"]
 
write.table(samples, "LCsamples.txt", na=".")                                   # save data including posteriors                    
samples     <- read.table("LCsamples.txt", na=".", header=TRUE)
samples     <- as.data.frame(samples)
new.samples <- NULL                                                             # turn into long format
sampmat     <- as.matrix(samples)                                               

for (i in 1:nrow(samples)){
  old.samples <- matrix(rep(sampmat[i,1:16],samples[i,6]), ncol=16, 
                        nrow=samples[i,6], byrow=T)
  new.samples <- rbind(new.samples,old.samples)
}

new.samples           <- as.data.frame(new.samples)
colnames(new.samples) <- c("woningREG","woningBACK","woningHOUSE","burgstat",
                           "toeslag","original","b1","b2","b3","b4","b5","p1",
                           "p2","p3","p4","p5")
new.samples$p1        <- as.numeric(as.character(new.samples$p1))
new.samples$p2        <- as.numeric(as.character(new.samples$p2))
new.samples$p3        <- as.numeric(as.character(new.samples$p3))
new.samples$p4        <- as.numeric(as.character(new.samples$p4))
new.samples$p5        <- as.numeric(as.character(new.samples$p5))

for (i in 1:nrow(new.samples)){                                                 
  for (j in 1:5){                                                               # create imputations
    new.samples[i,j+16] <- rbinom(1, 1, new.samples[i,j+11])
  }
}

colnames(new.samples) <- c("woningREG","woningBACK","woningHOUSE","burgstat",
                           "toeslag","original","b1","b2","b3","b4","b5","p1",
                           "p2","p3","p4","p5","i1","i2","i3","i4","i5")

new.samples$i1 <- as.factor(new.samples$i1)
new.samples$i2 <- as.factor(new.samples$i2)
new.samples$i3 <- as.factor(new.samples$i3)
new.samples$i4 <- as.factor(new.samples$i4)
new.samples$i5 <- as.factor(new.samples$i5)

levels(new.samples$i1) <- c("koop","huur")
levels(new.samples$i2) <- c("koop","huur")
levels(new.samples$i3) <- c("koop","huur")
levels(new.samples$i4) <- c("koop","huur")
levels(new.samples$i5) <- c("koop","huur")

gl1 <- glm(formula=i1~burgstat, family="binomial", new.samples)                 # logistic regression
gl2 <- glm(formula=i2~burgstat, family="binomial", new.samples)
gl3 <- glm(formula=i3~burgstat, family="binomial", new.samples)
gl4 <- glm(formula=i4~burgstat, family="binomial", new.samples)
gl5 <- glm(formula=i5~burgstat, family="binomial", new.samples)

Q1hat1 <- summary(gl1)$coefficients[1,1]
Q1hat2 <- summary(gl2)$coefficients[1,1]
Q1hat3 <- summary(gl3)$coefficients[1,1]
Q1hat4 <- summary(gl4)$coefficients[1,1]
Q1hat5 <- summary(gl5)$coefficients[1,1]

Q2hat1 <- summary(gl1)$coefficients[2,1]
Q2hat2 <- summary(gl2)$coefficients[2,1]
Q2hat3 <- summary(gl3)$coefficients[2,1]
Q2hat4 <- summary(gl4)$coefficients[2,1]
Q2hat5 <- summary(gl5)$coefficients[2,1]

U1hat1 <- summary(gl1)$coefficients[1,2]
U1hat2 <- summary(gl2)$coefficients[1,2]
U1hat3 <- summary(gl3)$coefficients[1,2]
U1hat4 <- summary(gl4)$coefficients[1,2]
U1hat5 <- summary(gl5)$coefficients[1,2]

U2hat1 <- summary(gl1)$coefficients[2,2]
U2hat2 <- summary(gl2)$coefficients[2,2]
U2hat3 <- summary(gl3)$coefficients[2,2]
U2hat4 <- summary(gl4)$coefficients[2,2]
U2hat5 <- summary(gl5)$coefficients[2,2]

Qbar1 <- mean(c(Q1hat1, Q1hat2, Q1hat3, Q1hat4, Q1hat5))
Qbar2 <- mean(c(Q2hat1, Q2hat2, Q2hat3, Q2hat4, Q2hat5))

Ubar1 <- (mean(c(U1hat1, U1hat2, U1hat3, U1hat4, U1hat5)))^2
Ubar2 <- (mean(c(U2hat1, U2hat2, U2hat3, U2hat4, U2hat5)))^2

Uhat1 <- c(U1hat1, U1hat2, U1hat3, U1hat4, U1hat5)
Uhat2 <- c(U2hat1, U2hat2, U2hat3, U2hat4, U2hat5)

B1 <- var(Uhat1)
B2 <- var(Uhat2)

T1 <- Ubar1+(5+1)*(B1/5)
T2 <- Ubar2+(5+1)*(B2/5)

sqrtT1  <- sqrt(T1)   
sqrtT2  <- sqrt(T2) 

ll1 <- Qbar1 - qt(.975, sum(table(new.samples$i1))-1)*sqrtT1
ll2 <- Qbar2 - qt(.975, sum(table(new.samples$i1))-1)*sqrtT1

ul1 <- Qbar1 + qt(.975, sum(table(new.samples$i1))-1)*sqrtT1
ul2 <- Qbar2 + qt(.975, sum(table(new.samples$i1))-1)*sqrtT1

ci1 <- ul1 - ll1
ci2 <- ul2 - ll2

model <- matrix(NA, nrow=2, ncol=8)
colnames(model) <- c("Qbar","Ubar","B","T","sqrtT","ll","ul","ci")
model[,"Qbar"]  <- c(Qbar1, Qbar2)
model[,"Ubar"]  <- c(Ubar1, Ubar2)
model[,"B"]     <- c(B1, B2)
model[,"T"]     <- c(T1, T2)
model[,"sqrtT"] <- c(sqrtT1, sqrtT2)
model[,"ll"]    <- c(ll1, ll2)
model[,"ul"]    <- c(ul1, ul2)
model[,"ci"]    <- c(ci1, ci2)

round(model, 4)


################################################################################
                                                                                # apply ML method
new.samples$i1 <- as.numeric(new.samples$i1)                                    # classification table
new.samples$i2 <- as.numeric(new.samples$i2)
new.samples$i3 <- as.numeric(new.samples$i3)
new.samples$i4 <- as.numeric(new.samples$i4)
new.samples$i5 <- as.numeric(new.samples$i5)

new.samples$i1 <- ifelse(new.samples$i1==2,0,1)
new.samples$i2 <- ifelse(new.samples$i2==2,0,1)
new.samples$i3 <- ifelse(new.samples$i3==2,0,1)
new.samples$i4 <- ifelse(new.samples$i4==2,0,1)
new.samples$i5 <- ifelse(new.samples$i5==2,0,1)

write.table(new.samples, "4_imputations.txt", row.names=FALSE, quote=FALSE)

short_impsam      <- list(NA)                                                 
short_impsam[[1]] <- ddply(new.samples, .(woningREG,woningBACK,p1), summarise,  # short datasets
                           W1=sum(1-i1), W2=sum(i1), N=(W1+W2))          
short_impsam[[2]] <- ddply(new.samples, .(woningREG,woningBACK,p2), summarise,  
                           W1=sum(1-i2), W2=sum(i2), N=(W1+W2))
short_impsam[[3]] <- ddply(new.samples, .(woningREG,woningBACK,p3), summarise,
                           W1=sum(1-i3), W2=sum(i3), N=(W1+W2))
short_impsam[[4]] <- ddply(new.samples, .(woningREG,woningBACK,p4), summarise,
                           W1=sum(1-i4), W2=sum(i4), N=(W1+W2))
short_impsam[[5]] <- ddply(new.samples, .(woningREG,woningBACK,p5), summarise,
                           W1=sum(1-i5), W2=sum(i5), N=(W1+W2))

for (j in 1:5){                                                                 
  for (k in 1:nrow(short_impsam[[j]]))
  {
    short_impsam[[j]][k,"PY"]   <- short_impsam[[j]][k,"N"]/sum(short_impsam[[j]][ ,"N"], na.rm=T)   # P(Y)                                  # P(Y)
    short_impsam[[j]][k,"PW1"]  <- short_impsam[[j]][k,"W1"]/sum(short_impsam[[j]][k,c("W1","W2")])  # P(W|Y)                             
    short_impsam[[j]][k,"PW2"]  <- short_impsam[[j]][k,"W2"]/sum(short_impsam[[j]][k,c("W1","W2")])
    short_impsam[[j]][k,"W1X1"] <- (short_impsam[[j]][k,"PY"]*short_impsam[[j]][k,paste0("p",j)]*
                                    short_impsam[[j]][k,"PW1"])                                      # P(W|X)
    short_impsam[[j]][k,"W1X2"] <- (short_impsam[[j]][k,"PY"]*(1-short_impsam[[j]][k,paste0("p",j)])*
                                    short_impsam[[j]][k,"PW1"])
    short_impsam[[j]][k,"W2X1"] <- (short_impsam[[j]][k,"PY"]*short_impsam[[j]][k,paste0("p",j)]*
                                    short_impsam[[j]][k,"PW2"])
    short_impsam[[j]][k,"W2X2"] <- (short_impsam[[j]][k,"PY"]*(1-short_impsam[[j]][k,paste0("p",j)])*
                                    short_impsam[[j]][k,"PW2"])
    short_impsam[[j]][k,"PX1"] <- short_impsam[[j]][k,paste0("p",j)]*(short_impsam[[j]][k,"N"]/      # P(X)
                                    sum(short_impsam[[j]][,"N"], na.rm=T))
    short_impsam[[j]][k,"PX2"] <- (1-short_impsam[[j]][k,paste0("p",j)])*(short_impsam[[j]][k,"N"]/
                                    sum(short_impsam[[j]][,"N"], na.rm=T))
  }
}

class_mat <- matrix(NA,5,4)                                                     
                                                                                
for (k in 1:5){                                                                 # classification errors
  class_mat[k,] <- c(sum(short_impsam[[k]][,"W1X1"], na.rm=T)/sum(short_impsam[[k]][,"PX1"], na.rm=T),
                     sum(short_impsam[[k]][,"W2X1"], na.rm=T)/sum(short_impsam[[k]][,"PX1"], na.rm=T),
                     sum(short_impsam[[k]][,"W1X2"], na.rm=T)/sum(short_impsam[[k]][,"PX2"], na.rm=T),
                     sum(short_impsam[[k]][,"W2X2"], na.rm=T)/sum(short_impsam[[k]][,"PX2"], na.rm=T))
  }

outfiles           <- paste0("6_model.",1:5)                                    
bootstrap.colnames <- tail(colnames(samples),5)

envir <- new.env()                                                              
run.template(template.path="6_ML1.brew", envir=envir,temp.filename.base=outfiles[1], 
             bootstrap.colnames=bootstrap.colnames[1])
envir <- new.env()
run.template(template.path="6_ML2.brew", envir=envir,temp.filename.base=outfiles[2], 
             bootstrap.colnames=bootstrap.colnames[2])
envir <- new.env()
run.template(template.path="6_ML3.brew", envir=envir,temp.filename.base=outfiles[3], 
             bootstrap.colnames=bootstrap.colnames[3])
envir <- new.env()
run.template(template.path="6_ML4.brew", envir=envir,temp.filename.base=outfiles[4], 
             bootstrap.colnames=bootstrap.colnames[4])
envir <- new.env()
run.template(template.path="6_ML5.brew", envir=envir,temp.filename.base=outfiles[5], 
             bootstrap.colnames=bootstrap.colnames[5])


ML1 <- read.delim(paste0("6_ML1.txt"), header=T, sep="\t")                      # posteriors
ML2 <- read.delim(paste0("6_ML2.txt"), header=T, sep="\t")                      
ML3 <- read.delim(paste0("6_ML3.txt"), header=T, sep="\t")
ML4 <- read.delim(paste0("6_ML4.txt"), header=T, sep="\t")
ML5 <- read.delim(paste0("6_ML5.txt"), header=T, sep="\t")

ML_posts <- cbind(ML1[,c(2:3)], ML2[,3], ML3[,3], ML4[,3], ML5[,3])           
imps     <- matrix(NA, nrow(ML_posts),5)                                      
ML_imps  <- cbind(ML_posts, imps)
colnames(ML_imps) <- c("burgstat","p1","p2","p3","p4","p5",                     
                       "i1","i2","i3","i4","i5")

for(k in 1:nrow(ML_posts)){                                                     # impute with posteriors
  ML_imps[k,"i1"] <- rbinom(1, 1, ML_imps[k,"p1"])
  ML_imps[k,"i2"] <- rbinom(1, 1, ML_imps[k,"p2"])
  ML_imps[k,"i3"] <- rbinom(1, 1, ML_imps[k,"p3"])
  ML_imps[k,"i4"] <- rbinom(1, 1, ML_imps[k,"p4"])
  ML_imps[k,"i5"] <- rbinom(1, 1, ML_imps[k,"p5"])
}

gl1 <- glm(formula=i1~burgstat, family="binomial", ML_imps)                     # regressies doen
gl2 <- glm(formula=i2~burgstat, family="binomial", ML_imps)
gl3 <- glm(formula=i3~burgstat, family="binomial", ML_imps)
gl4 <- glm(formula=i4~burgstat, family="binomial", ML_imps)
gl5 <- glm(formula=i5~burgstat, family="binomial", ML_imps)

Q1hat1 <- summary(gl1)$coefficients[1,1]
Q1hat2 <- summary(gl2)$coefficients[1,1]
Q1hat3 <- summary(gl3)$coefficients[1,1]
Q1hat4 <- summary(gl4)$coefficients[1,1]
Q1hat5 <- summary(gl5)$coefficients[1,1]

Q2hat1 <- summary(gl1)$coefficients[2,1]
Q2hat2 <- summary(gl2)$coefficients[2,1]
Q2hat3 <- summary(gl3)$coefficients[2,1]
Q2hat4 <- summary(gl4)$coefficients[2,1]
Q2hat5 <- summary(gl5)$coefficients[2,1]

U1hat1 <- summary(gl1)$coefficients[1,2]
U1hat2 <- summary(gl2)$coefficients[1,2]
U1hat3 <- summary(gl3)$coefficients[1,2]
U1hat4 <- summary(gl4)$coefficients[1,2]
U1hat5 <- summary(gl5)$coefficients[1,2]

U2hat1 <- summary(gl1)$coefficients[2,2]
U2hat2 <- summary(gl2)$coefficients[2,2]
U2hat3 <- summary(gl3)$coefficients[2,2]
U2hat4 <- summary(gl4)$coefficients[2,2]
U2hat5 <- summary(gl5)$coefficients[2,2]

Qbar1 <- mean(c(Q1hat1, Q1hat2, Q1hat3, Q1hat4, Q1hat5))
Qbar2 <- mean(c(Q2hat1, Q2hat2, Q2hat3, Q2hat4, Q2hat5))

Ubar1 <- (mean(c(U1hat1, U1hat2, U1hat3, U1hat4, U1hat5)))^2
Ubar2 <- (mean(c(U2hat1, U2hat2, U2hat3, U2hat4, U2hat5)))^2

Uhat1 <- c(U1hat1, U1hat2, U1hat3, U1hat4, U1hat5)
Uhat2 <- c(U2hat1, U2hat2, U2hat3, U2hat4, U2hat5)

B1 <- var(Uhat1)
B2 <- var(Uhat2)

T1 <- Ubar1+(5+1)*(B1/5)
T2 <- Ubar2+(5+1)*(B2/5)

sqrtT1  <- sqrt(T1)   
sqrtT2  <- sqrt(T2) 

ll1 <- Qbar1 - qt(.975, sum(table(ML_imps$i1))-1)*sqrtT1
ll2 <- Qbar2 - qt(.975, sum(table(ML_imps$i1))-1)*sqrtT1

ul1 <- Qbar1 + qt(.975, sum(table(ML_imps$i1))-1)*sqrtT1
ul2 <- Qbar2 + qt(.975, sum(table(ML_imps$i1))-1)*sqrtT1

ci1 <- ul1 - ll1
ci2 <- ul2 - ll2

model <- matrix(NA, nrow=2, ncol=8)
colnames(model) <- c("Qbar","Ubar","B","T","sqrtT","ll","ul","ci")
model[,"Qbar"]  <- c(Qbar1, Qbar2)
model[,"Ubar"]  <- c(Ubar1, Ubar2)
model[,"B"]     <- c(B1, B2)
model[,"T"]     <- c(T1, T2)
model[,"sqrtT"] <- c(sqrtT1, sqrtT2)
model[,"ll"]    <- c(ll1, ll2)
model[,"ul"]    <- c(ul1, ul2)
model[,"ci"]    <- c(ci1, ci2)

round(model, 4)

################################################################################

qpsolve <- function(e,d,iequal){                                                # function constrained BCH
  nr <- nrow(e)                                                                 
  nc <- ncol(e)
  ncel <- nr*nc
  evec <- as.vector(e)
  id <- diag(nr)
  p <- kronecker(t(d),id)
  dmat <- kronecker(d %*% t(d),id)
  dvec <- as.vector(evec %*% p)
  im <- diag(ncel)
  i1 <- iequal
  i2 <- setdiff(1:ncel,i1)
  index <- c(i1,i2)
  im2 <- im[index,]
  at <- rbind(rep(1,ncel),im2)
  amat <- t(at)
  bvec <- c(1,rep(0,ncel))
  meq <- 1 + length(iequal)
  res <- solve.QP(dmat,dvec,amat,bvec,meq)
  return(res)
}

new.samples <- read.table("4_imputations.txt", header=TRUE)

D1 <- matrix(c(0.9680783, 0.03192170, 0.014753061, 0.9852469), ncol=2, byrow=T) # D matrix = class table
E1 <- matrix(table(new.samples$burgstat, new.samples$i1), ncol=2)               
E1[1,] <- E1[1,]/sum(E1[1,])                                                    # E matrix = P(W,Q)
E1[2,] <- E1[2,]/sum(E1[2,])                                                    
                                                                                
D2 <- matrix(c(0.9694548, 0.03054524, 0.014522061, 0.9854779), ncol=2, byrow=T)
E2 <- matrix(table(new.samples$burgstat, new.samples$i2), ncol=2)
E2[1,] <- E2[1,]/sum(E2[1,])
E2[2,] <- E2[2,]/sum(E2[2,])

D3 <- matrix(c(0.9807727, 0.01922726, 0.007157617, 0.9928424), ncol=2, byrow=T)
E3 <- matrix(table(new.samples$burgstat, new.samples$i3), ncol=2)
E3[1,] <- E3[1,]/sum(E3[1,])
E3[2,] <- E3[2,]/sum(E3[2,])

D4 <- matrix(c(0.9742938, 0.02570619, 0.012251459, 0.9877485), ncol=2, byrow=T)
E4 <- matrix(table(new.samples$burgstat, new.samples$i4), ncol=2)
E4[1,] <- E4[1,]/sum(E4[1,])
E4[2,] <- E4[2,]/sum(E4[2,])

D5 <- matrix(c(0.9653227, 0.03467734, 0.025041735, 0.9749583), ncol=2, byrow=T)
E5 <- matrix(table(new.samples$burgstat, new.samples$i5), ncol=2)
E5[1,] <- E5[1,]/sum(E5[1,])
E5[2,] <- E5[2,]/sum(E5[2,])

iequal <- c()                                                                   # no constraints

A1 <- qpsolve(E1,D1,iequal)                                                     # Calculate A matrix 
A2 <- qpsolve(E2,D2,iequal)                                                     
A3 <- qpsolve(E3,D3,iequal)
A4 <- qpsolve(E4,D4,iequal)
A5 <- qpsolve(E5,D5,iequal)

##########

BCHpost <- matrix(NA, 4, 12)                                                    # korte dataset met 4 rijen
rownames(BCHpost) <- c("b0k0","b1k0","b0k1","b1k1")                             # voor de 4 scores X*Z
colnames(BCHpost) <- c("f1","f2","f3","f4","f5",
                       "p1","p2","p3","p4","p5","b","k")
BCHpost[,1] <- table(new.samples$burgstat, new.samples$i1)                      # frequencies van de
BCHpost[,2] <- table(new.samples$burgstat, new.samples$i2)                      # bootstraps
BCHpost[,3] <- table(new.samples$burgstat, new.samples$i3)
BCHpost[,4] <- table(new.samples$burgstat, new.samples$i4)
BCHpost[,5] <- table(new.samples$burgstat, new.samples$i5)

BCHpost[,6] <- c(A1$solution[1]/sum(A1$solution[c(1,3)]),                       # posteriors voor elke
                 A1$solution[2]/sum(A1$solution[c(2,4)]),                       # boostrap sample
                 A1$solution[3]/sum(A1$solution[c(1,3)]),
                 A1$solution[4]/sum(A1$solution[c(2,4)]))

BCHpost[,7] <- c(A2$solution[1]/sum(A2$solution[c(1,3)]),
                 A2$solution[2]/sum(A2$solution[c(2,4)]),
                 A2$solution[3]/sum(A2$solution[c(1,3)]),
                 A2$solution[4]/sum(A2$solution[c(2,4)]))

BCHpost[,8] <- c(A3$solution[1]/sum(A3$solution[c(1,3)]),
                 A3$solution[2]/sum(A3$solution[c(2,4)]),
                 A3$solution[3]/sum(A3$solution[c(1,3)]),
                 A3$solution[4]/sum(A3$solution[c(2,4)]))

BCHpost[,9] <- c(A4$solution[1]/sum(A4$solution[c(1,3)]),
                 A4$solution[2]/sum(A4$solution[c(2,4)]),
                 A4$solution[3]/sum(A4$solution[c(1,3)]),
                 A4$solution[4]/sum(A4$solution[c(2,4)]))

BCHpost[,10] <- c(A5$solution[1]/sum(A5$solution[c(1,3)]),
                 A5$solution[2]/sum(A5$solution[c(2,4)]),
                 A5$solution[3]/sum(A5$solution[c(1,3)]),
                 A5$solution[4]/sum(A5$solution[c(2,4)]))

BCHpost[,11] <- c(1,2,1,2)
BCHpost[,12] <- c(1,1,2,2)

##########

new.matrix <- NULL                                                              # change data into long format
for (j in 1:4){                                                                 
  old.matrix <- matrix(rep(BCHpost[j,c(6,11:12)],BCHpost[j,1]),                 
                       BCHpost[j,1],3, byrow=T)
  new.matrix <- rbind(old.matrix, new.matrix)
}
data1 <- new.matrix

new.matrix <- NULL
for (j in 1:4){
  old.matrix <- matrix(rep(BCHpost[j,c(7,11:12)],BCHpost[j,2]),
                       BCHpost[j,2],3, byrow=T)
  new.matrix <- rbind(old.matrix, new.matrix)
}
data2 <- new.matrix

new.matrix <- NULL
for (j in 1:4){
  old.matrix <- matrix(rep(BCHpost[j,c(8,11:12)],BCHpost[j,3]),
                       BCHpost[j,3],3, byrow=T)
  new.matrix <- rbind(old.matrix, new.matrix)
}
data3 <- new.matrix

new.matrix <- NULL
for (j in 1:4){
  old.matrix <- matrix(rep(BCHpost[j,c(9,11:12)],BCHpost[j,4]),
                       BCHpost[j,4],3, byrow=T)
  new.matrix <- rbind(old.matrix, new.matrix)
}
data4 <- new.matrix

new.matrix <- NULL
for (j in 1:4){
  old.matrix <- matrix(rep(BCHpost[j,c(10,11:12)],BCHpost[j,4]),
                       BCHpost[j,4],3, byrow=T)
  new.matrix <- rbind(old.matrix, new.matrix)
}
data5 <- new.matrix

##########

imp1 <- rep(NA, nrow(data1))                                                    
imp2 <- rep(NA, nrow(data2))                                                    
imp3 <- rep(NA, nrow(data3))
imp4 <- rep(NA, nrow(data4))
imp5 <- rep(NA, nrow(data5))

data1 <- cbind(data1, imp1)
data2 <- cbind(data2, imp2)
data3 <- cbind(data3, imp3)
data4 <- cbind(data4, imp4)
data5 <- cbind(data5, imp5)

colnames(data1) <- colnames(data2) <- colnames(data3) <- colnames(data4) <-     
  colnames(data5) <- c("post","burg","koop","imp")

for (i in 1:nrow(data1)){                                                       
  data1[i,"imp"] <- rbinom(1, 1, data1[i,"post"])                               # new imputations
}                                                                               

for (i in 1:nrow(data2)){
  data2[i,"imp"] <- rbinom(1, 1, data2[i,"post"])
}

for (i in 1:nrow(data3)){
  data3[i,"imp"] <- rbinom(1, 1, data3[i,"post"])
}

for (i in 1:nrow(data4)){
  data4[i,"imp"] <- rbinom(1, 1, data4[i,"post"])
}

for (i in 1:nrow(data5)){
  data5[i,"imp"] <- rbinom(1, 1, data5[i,"post"])
}

data1 <- as.data.frame(data1)
data2 <- as.data.frame(data2)
data3 <- as.data.frame(data3)
data4 <- as.data.frame(data4)
data5 <- as.data.frame(data5)

data1$imp <- as.factor(data1$imp)
data2$imp <- as.factor(data2$imp)
data3$imp <- as.factor(data3$imp)
data4$imp <- as.factor(data4$imp)
data5$imp <- as.factor(data5$imp)

levels(data1$imp) <- levels(data2$imp) <- levels(data3$imp) <- 
  levels(data4$imp) <- levels(data5$imp) <- c("huur","koop")

##########

gl1 <- glm(formula=imp~burg, family="binomial", data1)                          # logistic regressions
gl2 <- glm(formula=imp~burg, family="binomial", data2)
gl3 <- glm(formula=imp~burg, family="binomial", data3)
gl4 <- glm(formula=imp~burg, family="binomial", data4)
gl5 <- glm(formula=imp~burg, family="binomial", data5)

Q1hat1 <- summary(gl1)$coefficients[1,1]
Q1hat2 <- summary(gl2)$coefficients[1,1]
Q1hat3 <- summary(gl3)$coefficients[1,1]
Q1hat4 <- summary(gl4)$coefficients[1,1]
Q1hat5 <- summary(gl5)$coefficients[1,1]

Q2hat1 <- summary(gl1)$coefficients[2,1]
Q2hat2 <- summary(gl2)$coefficients[2,1]
Q2hat3 <- summary(gl3)$coefficients[2,1]
Q2hat4 <- summary(gl4)$coefficients[2,1]
Q2hat5 <- summary(gl5)$coefficients[2,1]


U1hat1 <- summary(gl1)$coefficients[1,2]
U1hat2 <- summary(gl2)$coefficients[1,2]
U1hat3 <- summary(gl3)$coefficients[1,2]
U1hat4 <- summary(gl4)$coefficients[1,2]
U1hat5 <- summary(gl5)$coefficients[1,2]

U2hat1 <- summary(gl1)$coefficients[2,2]
U2hat2 <- summary(gl2)$coefficients[2,2]
U2hat3 <- summary(gl3)$coefficients[2,2]
U2hat4 <- summary(gl4)$coefficients[2,2]
U2hat5 <- summary(gl5)$coefficients[2,2]

Qbar1 <- mean(c(Q1hat1, Q1hat2, Q1hat3, Q1hat4, Q1hat5))
Qbar2 <- mean(c(Q2hat1, Q2hat2, Q2hat3, Q2hat4, Q2hat5))

Ubar1 <- (mean(c(U1hat1, U1hat2, U1hat3, U1hat4, U1hat5)))^2
Ubar2 <- (mean(c(U2hat1, U2hat2, U2hat3, U2hat4, U2hat5)))^2

Uhat1 <- c(U1hat1, U1hat2, U1hat3, U1hat4, U1hat5)
Uhat2 <- c(U2hat1, U2hat2, U2hat3, U2hat4, U2hat5)

B1 <- var(Uhat1)
B2 <- var(Uhat2)

T1 <- Ubar1+(5+1)*(B1/5)
T2 <- Ubar2+(5+1)*(B2/5)

sqrtT1  <- sqrt(T1)   
sqrtT2  <- sqrt(T2) 

ll1 <- Qbar1 - qt(.975, sum(table(data1$imp))-1)*sqrtT1
ll2 <- Qbar2 - qt(.975, sum(table(data1$imp))-1)*sqrtT1

ul1 <- Qbar1 + qt(.975, sum(table(data1$imp))-1)*sqrtT1
ul2 <- Qbar2 + qt(.975, sum(table(data1$imp))-1)*sqrtT1

ci1 <- ul1 - ll1
ci2 <- ul2 - ll2

model <- matrix(NA, nrow=2, ncol=8)
colnames(model) <- c("Qbar","Ubar","B","T","sqrtT","ll","ul","ci")
model[,"Qbar"]  <- c(Qbar1, Qbar2)
model[,"Ubar"]  <- c(Ubar1, Ubar2)
model[,"B"]     <- c(B1, B2)
model[,"T"]     <- c(T1, T2)
model[,"sqrtT"] <- c(sqrtT1, sqrtT2)
model[,"ll"]    <- c(ll1, ll2)
model[,"ul"]    <- c(ul1, ul2)
model[,"ci"]    <- c(ci1, ci2)

round(model, 4)
