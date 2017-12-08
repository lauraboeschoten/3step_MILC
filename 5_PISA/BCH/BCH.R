setwd("F:datapackage_threestepMILC/5_PISA/BCH")

library(quadprog)

options(scipen=999)
set.seed(123)

qpsolve <- function(e,d,iequal){                                                # function for BCH
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


pdata1 <- read.delim("6_ML1.txt", header=T)                                     # imputed data from model
pdata2 <- read.delim("6_ML2.txt", header=T)                                     # without covariate (copy'd from ML folder)
pdata3 <- read.delim("6_ML3.txt", header=T)
pdata4 <- read.delim("6_ML4.txt", header=T)
pdata5 <- read.delim("6_ML5.txt", header=T)

D1 <- matrix(c(0.95902850, 0.0409715, 0.04918273, 0.9508173), ncol=2, byrow=T)  # matrices copy'd from ML
E1 <- matrix(table(pdata1$ST006Q04TA, pdata1$imp), ncol=2)               
E1 <- E1/sum(E1)

D2 <- matrix(c(0.04555412, 0.95444588, 0.98935465, 0.01064535), ncol=2, byrow=F)
E2 <- matrix(table(pdata2$ST006Q04TA, pdata2$imp), ncol=2)
E2 <- E2/sum(E2)

D3 <- matrix(c(0.98769676, 0.01230324, 0.05727163, 0.94272837), ncol=2, byrow=F)
E3 <- matrix(table(pdata3$ST006Q04TA, pdata3$imp), ncol=2)
E3 <- E3/sum(E3)

D4 <- matrix(c(0.04325209, 0.95674791, 0.97356798, 0.02643202), ncol=2, byrow=F)
E4 <- matrix(table(pdata4$ST006Q04TA, pdata4$imp), ncol=2)
E4 <- E4/sum(E4)

D5 <- matrix(c(0.9543027, 0.0456973, 0.0463535, 0.9536465), ncol=2, byrow=F)
E5 <- matrix(table(pdata5$ST006Q04TA, pdata5$imp), ncol=2)
E5 <- E5/sum(E5)

iequal <- c()                                                                   # no constraints

A1 <- qpsolve(E1,D1,iequal)                                                     # estimate A matrix
A2 <- qpsolve(E2,D2,iequal)                                                     
A3 <- qpsolve(E3,D3,iequal)
A4 <- qpsolve(E4,D4,iequal)
A5 <- qpsolve(E5,D5,iequal)


##########

BCHpost1 <- matrix(NA, 2, 2)                                                    
colnames(BCHpost1) <- c("theta1","theta2")                             
rownames(BCHpost1) <- c("g1","g2")
BCHpost1[1,1] <- A1$solution[1]/sum(A1$solution[c(1,3)])
BCHpost1[1,2] <- A1$solution[3]/sum(A1$solution[c(1,3)])
BCHpost1[2,1] <- A1$solution[2]/sum(A1$solution[c(2,4)])
BCHpost1[2,2] <- A1$solution[4]/sum(A1$solution[c(2,4)])                        # turn into posteriors

freq1 <- table(pdata1$ST006Q04TA, pdata1$imp)   
mom <- c("no","yes")
BCHpost1 <- cbind(BCHpost1, mom)


new.matrix <- NULL                                                              # turn data in long format            
for (j in 1:2){                                                                 
  old.matrix <- matrix(rep(BCHpost1[j,],sum(freq1[j,])),                 
                       sum(freq1[j,]),3, byrow=T)
  new.matrix <- rbind(old.matrix, new.matrix)
}
data1 <- new.matrix
imp1 <- rep(NA, nrow(data1)) 
data1 <- cbind(data1, imp1)
colnames(data1) <- c("theta1","theta2","mom","imp1")

for (i in 1:nrow(data1)){                                                       # create imputations 
  data1[i,"imp1"] <- which(rmultinom(1, 1, as.numeric(as.character(data1[i,c(1:2)])))==1) # check for every bootstrap
}                                                                               # if you have label switching! 

data1 <- as.data.frame(data1)
data1$imp1 <- as.factor(data1$imp1)

lm1 <- glm(formula=mom~imp1, family="binomial", data1)
summary(lm1)

################################################################################

BCHpost2 <- matrix(NA, 2, 2)                                          
colnames(BCHpost2) <- c("theta1","theta2")                            
rownames(BCHpost2) <- c("g1","g2")
BCHpost2[1,1] <- A2$solution[1]/sum(A2$solution[c(1,3)])
BCHpost2[1,2] <- A2$solution[3]/sum(A2$solution[c(1,3)])
BCHpost2[2,1] <- A2$solution[2]/sum(A2$solution[c(2,4)])
BCHpost2[2,2] <- A2$solution[4]/sum(A2$solution[c(2,4)])
#BCHpost2 <- round(BCHpost2, 10)

freq2 <- table(pdata2$ST006Q04TA, pdata2$imp)   
mom <- c("no","yes")
BCHpost2 <- cbind(BCHpost2, mom)

new.matrix <- NULL                                                    
for (j in 1:2){                                                       
  old.matrix <- matrix(rep(BCHpost2[j,],sum(freq2[j,])),              
                       sum(freq2[j,]),3, byrow=T)
  new.matrix <- rbind(old.matrix, new.matrix)
}
data2 <- new.matrix
imp2  <- rep(NA, nrow(data2)) 
data2 <- cbind(data2, imp2)
colnames(data2) <- c("theta1","theta2","mom","imp2")

for (i in 1:nrow(data2)){                                             
  data2[i,"imp2"] <- which(rmultinom(1, 1, as.numeric(as.character(data2[i,c(2,1)])))==1)
} 

data2 <- as.data.frame(data2)
data2$imp2 <- as.factor(data2$imp2)
lm2 <- glm(formula=mom~imp2, family="binomial", data2)  
summary(lm2)

################################################################################

BCHpost3 <- matrix(NA, 2, 2)                                                 
colnames(BCHpost3) <- c("theta1","theta2")                             
rownames(BCHpost3) <- c("g1","g2")
BCHpost3[1,1] <- A3$solution[1]/sum(A3$solution[c(1,3)])
BCHpost3[1,2] <- A3$solution[3]/sum(A3$solution[c(1,3)])
BCHpost3[2,1] <- A3$solution[2]/sum(A3$solution[c(2,4)])
BCHpost3[2,2] <- A3$solution[4]/sum(A3$solution[c(2,4)])
BCHpost3 <- round(BCHpost3, 10)

freq3 <- table(pdata3$ST006Q04TA, pdata3$imp)   
mom <- c("no","yes")
BCHpost3 <- cbind(BCHpost3, mom)

new.matrix <- NULL                                                     
for (j in 1:2){                                                        
  old.matrix <- matrix(rep(BCHpost3[j,],sum(freq3[j,])),               
                       sum(freq3[j,]),3, byrow=T)
  new.matrix <- rbind(old.matrix, new.matrix)
}
data3 <- new.matrix
imp3 <- rep(NA, nrow(data3)) 
data3 <- cbind(data3, imp3)
colnames(data3) <- c("theta1","theta2","mom","imp3")

for (i in 1:nrow(data3)){                                                      
  data3[i,"imp3"] <- which(rmultinom(1, 1, as.numeric(as.character(data3[i,c(1:2)])))==1) 
} 

data3 <- as.data.frame(data3)
data3$imp3 <- as.factor(data3$imp3)
lm3 <- glm(formula=mom~imp3, family="binomial", data3) 
summary(lm3)

################################################################################

BCHpost4 <- matrix(NA, 2, 2)                                                   
colnames(BCHpost4) <- c("theta1","theta2")                             
rownames(BCHpost4) <- c("g1","g2")
BCHpost4[1,1] <- A4$solution[1]/sum(A4$solution[c(1,3)])
BCHpost4[1,2] <- A4$solution[3]/sum(A4$solution[c(1,3)])
BCHpost4[2,1] <- A4$solution[2]/sum(A4$solution[c(2,4)])
BCHpost4[2,2] <- A4$solution[4]/sum(A4$solution[c(2,4)])
BCHpost4 <- round(BCHpost4, 10)

freq4 <- table(pdata4$ST006Q04TA, pdata4$imp)   
mom <- c("no","yes")
BCHpost4 <- cbind(BCHpost4, mom)


new.matrix <- NULL                                                            
for (j in 1:2){                                                         
  old.matrix <- matrix(rep(BCHpost4[j,],sum(freq4[j,])),                
                       sum(freq4[j,]),3, byrow=T)
  new.matrix <- rbind(old.matrix, new.matrix)
}
data4 <- new.matrix
imp4 <- rep(NA, nrow(data4)) 
data4 <- cbind(data4, imp4)
colnames(data4) <- c("theta1","theta2","mom","imp4")

for (i in 1:nrow(data4)){                                               
  data4[i,"imp4"] <- which(rmultinom(1, 1, as.numeric(as.character(data4[i,c("theta2","theta1")])))==1)
} 

data4 <- as.data.frame(data4)

data4$imp4 <- as.factor(data4$imp4)
lm4 <- glm(formula=mom~imp4, family="binomial", data4)  
summary(lm4)

################################################################################

BCHpost5 <- matrix(NA, 2, 2)                                                    
colnames(BCHpost5) <- c("theta1","theta2")                             
rownames(BCHpost5) <- c("g1","g2")
BCHpost5[1,1] <- A5$solution[1]/sum(A5$solution[c(1,3)])
BCHpost5[1,2] <- A5$solution[3]/sum(A5$solution[c(1,3)])
BCHpost5[2,1] <- A5$solution[2]/sum(A5$solution[c(2,4)])
BCHpost5[2,2] <- A5$solution[4]/sum(A5$solution[c(2,4)])
BCHpost5 <- round(BCHpost5, 10)

freq5 <- table(pdata5$ST006Q04TA, pdata5$imp)   
mom <- c("no","yes")
BCHpost5 <- cbind(BCHpost5, mom)

new.matrix <- NULL                                                     
for (j in 1:2){                                                        
  old.matrix <- matrix(rep(BCHpost5[j,],sum(freq5[j,])),               
                       sum(freq5[j,]),3, byrow=T)
  new.matrix <- rbind(old.matrix, new.matrix)
}
data5 <- new.matrix
imp5 <- rep(NA, nrow(data5)) 
data5 <- cbind(data5, imp5)
colnames(data5) <- c("theta1","theta2","mom","imp5")

for (i in 1:nrow(data5)){                                              
  data5[i,"imp5"] <- which(rmultinom(1, 1, as.numeric(as.character(data5[i,c("theta1","theta2")])))==1)
} 

data5 <- as.data.frame(data5)

data5$imp5 <- as.factor(data5$imp5)
lm5 <- glm(formula=mom~imp5, family="binomial", data5)  
summary(lm5)

################################################################################

Q1hat1 <- summary(lm1)$coefficients[1,1]                                        # pool results
Q1hat2 <- summary(lm2)$coefficients[1,1]
Q1hat3 <- summary(lm3)$coefficients[1,1]
Q1hat4 <- summary(lm4)$coefficients[1,1]
Q1hat5 <- summary(lm5)$coefficients[1,1]

Q2hat1 <- summary(lm1)$coefficients[2,1]
Q2hat2 <- summary(lm2)$coefficients[2,1]
Q2hat3 <- summary(lm3)$coefficients[2,1]
Q2hat4 <- summary(lm4)$coefficients[2,1]
Q2hat5 <- summary(lm5)$coefficients[2,1]

U1hat1 <- summary(lm1)$coefficients[1,2]
U1hat2 <- summary(lm2)$coefficients[1,2]
U1hat3 <- summary(lm3)$coefficients[1,2]
U1hat4 <- summary(lm4)$coefficients[1,2]
U1hat5 <- summary(lm5)$coefficients[1,2]

U2hat1 <- summary(lm1)$coefficients[2,2]
U2hat2 <- summary(lm2)$coefficients[2,2]
U2hat3 <- summary(lm3)$coefficients[2,2]
U2hat4 <- summary(lm4)$coefficients[2,2]
U2hat5 <- summary(lm5)$coefficients[2,2]

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

ll1 <- Qbar1 - qt(.975, (nrow(data1)-1)-1)*sqrtT1
ll2 <- Qbar2 - qt(.975, (nrow(data1)-1)-1)*sqrtT2

ul1 <- Qbar1 + qt(.975, (nrow(data1)-1)-1)*sqrtT1
ul2 <- Qbar2 + qt(.975, (nrow(data1)-1)-1)*sqrtT2

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
