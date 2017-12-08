setwd("F:datapackage_threestepMILC/5_PISA/ML")

set.seed(123)
library(plyr)
library(brew)

source("runtemplate.R")

indicatoren <- read.table("pisa.txt", header=TRUE)
imp1        <- read.table("pdata1.txt", header=TRUE)
all1        <- cbind(indicatoren, imp1)
kort1       <- ddply(all1,.(CM033Q01S,CM474Q01S,CM155Q01S,CM155Q04S,
                            CM411Q01S,CM411Q02S,CM803Q01S,CM442Q02S,
                            DM462Q01C,CM034Q01S,CM305Q01S,CM496Q01S,
                            CM496Q02S,CM423Q01S,DM406Q01C,DM406Q02C,
                            CM603Q01S,CM571Q01S,CM564Q01S,CM564Q02S,
                            X1,X2), summarise, W1=sum(imp==1),
                            W2=sum(imp==2), N=(W1+W2))

for (k in 1:nrow(kort1))
{
  kort1[k,"PY"] <- kort1[k,"N"]/sum(kort1[ ,"N"], na.rm=T)                      # P(Y)
  kort1[k,"PW1"] <- kort1[k,"W1"]/sum(kort1[k,c("W1","W2")])                    # P(W|Y)
  kort1[k,"PW2"] <- kort1[k,"W2"]/sum(kort1[k,c("W1","W2")])
  kort1[k,"W1X1"] <- (kort1[k,"PY"]*kort1[k,paste0("X1")]*kort1[k,"PW1"])       # P(W|X)
  kort1[k,"W1X2"] <- (kort1[k,"PY"]*kort1[k,paste0("X2")]*kort1[k,"PW1"])
  kort1[k,"W2X1"] <- (kort1[k,"PY"]*kort1[k,paste0("X1")]*kort1[k,"PW2"])
  kort1[k,"W2X2"] <- (kort1[k,"PY"]*kort1[k,paste0("X2")]*kort1[k,"PW2"])
  kort1[k,"PX1"] <- kort1[k,paste0("X1")]*(kort1[k,"N"]/sum(kort1[,"N"], na.rm=T)) # P(X)
  kort1[k,"PX2"] <- kort1[k,paste0("X2")]*(kort1[k,"N"]/sum(kort1[,"N"], na.rm=T))
}

class_mat1           <- matrix(NA,2,2)                                                    
class_mat1[1,1]      <- sum(kort1[,"W1X1"], na.rm=T)/sum(kort1[,"PX1"], na.rm=T)
class_mat1[1,2]      <- sum(kort1[,"W2X1"], na.rm=T)/sum(kort1[,"PX1"], na.rm=T)
class_mat1[2,1]      <- sum(kort1[,"W1X2"], na.rm=T)/sum(kort1[,"PX2"], na.rm=T)
class_mat1[2,2]      <- sum(kort1[,"W2X2"], na.rm=T)/sum(kort1[,"PX2"], na.rm=T)
colnames(class_mat1) <- c("W1","W2")
rownames(class_mat1) <- c("X1","X2")                                            # classification error
outfiles             <- paste0("6_model.",1:5)                                    
bootstrap.colnames   <- c("b1","b2","b3","b4","b5")

envir <- new.env()                                                              # run 3 step model
run.template(template.path="6_ML1.brew", envir=envir, temp.filename.base=outfiles[1], 
             bootstrap.colnames=bootstrap.colnames[1])

lc1           <- read.delim("6_ML1.txt", header=T)
lc1imp        <- lc1[,c("imp","ST006Q04TA","theta.1","theta.2")]
lc1imp$newimp <- NULL

for (i in 1:nrow(lc1imp)){                                                      # create imputations
  lc1imp[i,"newimp"] <- which(rmultinom(1,1,lc1imp[i,c("theta.1","theta.2")])==1) # again watch out
}                                                                               # for label switching!

lc1imp$newimp <- as.factor(lc1imp$newimp)
gl1 <- glm(formula=ST006Q04TA~newimp, family="binomial", lc1imp)              
summary(gl1)

################################################################################

imp2 <- read.table("pdata2.txt", header=TRUE)

all2 <- cbind(indicatoren, imp2)

kort1 <- ddply(all2,.(CM033Q01S,CM474Q01S,CM155Q01S,CM155Q04S,
                      CM411Q01S,CM411Q02S,CM803Q01S,CM442Q02S,
                      DM462Q01C,CM034Q01S,CM305Q01S,CM496Q01S,
                      CM496Q02S,CM423Q01S,DM406Q01C,DM406Q02C,
                      CM603Q01S,CM571Q01S,CM564Q01S,CM564Q02S,
                      X1,X2), summarise, W1=sum(imp==1),
               W2=sum(imp==2), N=(W1+W2))

for (k in 1:nrow(kort1))
{
  kort1[k,"PY"] <- kort1[k,"N"]/sum(kort1[ ,"N"], na.rm=T)                  
  kort1[k,"PW1"] <- kort1[k,"W1"]/sum(kort1[k,c("W1","W2")])             
  kort1[k,"PW2"] <- kort1[k,"W2"]/sum(kort1[k,c("W1","W2")])
  # kans op W gegeven X
  kort1[k,"W1X1"] <- (kort1[k,"PY"]*kort1[k,paste0("X1")]*kort1[k,"PW1"])
  kort1[k,"W1X2"] <- (kort1[k,"PY"]*kort1[k,paste0("X2")]*kort1[k,"PW1"])
  kort1[k,"W2X1"] <- (kort1[k,"PY"]*kort1[k,paste0("X1")]*kort1[k,"PW2"])
  kort1[k,"W2X2"] <- (kort1[k,"PY"]*kort1[k,paste0("X2")]*kort1[k,"PW2"])
  
  kort1[k,"PX1"] <- kort1[k,paste0("X1")]*(kort1[k,"N"]/sum(kort1[,"N"], na.rm=T))
  kort1[k,"PX2"] <- kort1[k,paste0("X2")]*(kort1[k,"N"]/sum(kort1[,"N"], na.rm=T))
}

class_mat1 <- matrix(NA,2,2)                                                  
class_mat1[1,1] <- sum(kort1[,"W1X1"], na.rm=T)/sum(kort1[,"PX1"], na.rm=T)
class_mat1[1,2] <- sum(kort1[,"W2X1"], na.rm=T)/sum(kort1[,"PX1"], na.rm=T)

class_mat1[2,1] <- sum(kort1[,"W1X2"], na.rm=T)/sum(kort1[,"PX2"], na.rm=T)
class_mat1[2,2] <- sum(kort1[,"W2X2"], na.rm=T)/sum(kort1[,"PX2"], na.rm=T)

colnames(class_mat1) <- c("W1","W2")
rownames(class_mat1) <- c("X1","X2")

class_mat1

outfiles           <- paste0("6_model.",1:5)                            
bootstrap.colnames <- c("b1","b2","b3","b4","b5")

envir <- new.env()                                                          
run.template(template.path="6_ML2.brew", envir=envir, 
             temp.filename.base=outfiles[2], 
             bootstrap.colnames=bootstrap.colnames[1])

lc2 <- read.delim("6_ML2.txt", header=T)
lc2imp <- lc2[,c("imp","ST006Q04TA","theta.1","theta.2")]
lc2imp$newimp <- NULL

for (i in 1:nrow(lc2imp)){
  lc2imp[i,"newimp"] <- which(rmultinom(1,1,lc2imp[i,c("theta.2","theta.1")])==1)
}

lc2imp$newimp <- as.factor(lc2imp$newimp)
gl2 <- glm(formula=ST006Q04TA~newimp, family="binomial", lc2imp)          
summary(gl2)

################################################################################

imp3 <- read.table("pdata3.txt", header=TRUE)

all3 <- cbind(indicatoren, imp3)

kort1 <- ddply(all3,.(CM033Q01S,CM474Q01S,CM155Q01S,CM155Q04S,
                      CM411Q01S,CM411Q02S,CM803Q01S,CM442Q02S,
                      DM462Q01C,CM034Q01S,CM305Q01S,CM496Q01S,
                      CM496Q02S,CM423Q01S,DM406Q01C,DM406Q02C,
                      CM603Q01S,CM571Q01S,CM564Q01S,CM564Q02S,
                      X2,X1), summarise, W1=sum(imp==1),
               W2=sum(imp==2), N=(W1+W2))

for (k in 1:nrow(kort1))
{
  kort1[k,"PY"] <- kort1[k,"N"]/sum(kort1[ ,"N"], na.rm=T)            
  kort1[k,"PW1"] <- kort1[k,"W1"]/sum(kort1[k,c("W1","W2")])             
  kort1[k,"PW2"] <- kort1[k,"W2"]/sum(kort1[k,c("W1","W2")])
  # kans op W gegeven X
  kort1[k,"W1X1"] <- (kort1[k,"PY"]*kort1[k,paste0("X1")]*kort1[k,"PW1"])
  kort1[k,"W1X2"] <- (kort1[k,"PY"]*kort1[k,paste0("X2")]*kort1[k,"PW1"])
  kort1[k,"W2X1"] <- (kort1[k,"PY"]*kort1[k,paste0("X1")]*kort1[k,"PW2"])
  kort1[k,"W2X2"] <- (kort1[k,"PY"]*kort1[k,paste0("X2")]*kort1[k,"PW2"])
  
  kort1[k,"PX1"] <- kort1[k,paste0("X1")]*(kort1[k,"N"]/sum(kort1[,"N"], na.rm=T))
  kort1[k,"PX2"] <- kort1[k,paste0("X2")]*(kort1[k,"N"]/sum(kort1[,"N"], na.rm=T))
}

class_mat1 <- matrix(NA,2,2)                                                     
class_mat1[1,1] <- sum(kort1[,"W1X1"], na.rm=T)/sum(kort1[,"PX1"], na.rm=T)
class_mat1[1,2] <- sum(kort1[,"W2X1"], na.rm=T)/sum(kort1[,"PX1"], na.rm=T)
class_mat1[2,1] <- sum(kort1[,"W1X2"], na.rm=T)/sum(kort1[,"PX2"], na.rm=T)
class_mat1[2,2] <- sum(kort1[,"W2X2"], na.rm=T)/sum(kort1[,"PX2"], na.rm=T)
colnames(class_mat1) <- c("W1","W2")
rownames(class_mat1) <- c("X1","X2")

class_mat1

outfiles           <- paste0("6_model.",1:5)                  
bootstrap.colnames <- c("b1","b2","b3","b4","b5")

envir <- new.env()                                                              
run.template(template.path="6_ML3.brew", envir=envir, 
             temp.filename.base=outfiles[3], 
             bootstrap.colnames=bootstrap.colnames[1])

lc3 <- read.delim("6_ML3.txt", header=T)
lc3imp <- lc3[,c("imp","ST006Q04TA","theta.1","theta.2")]
lc3imp$newimp <- NULL

for (i in 1:nrow(lc3imp)){
  lc3imp[i,"newimp"] <- which(rmultinom(1,1,lc3imp[i,c("theta.1","theta.2")])==1)
}

lc3imp$newimp <- as.factor(lc3imp$newimp)
gl3 <- glm(formula=ST006Q04TA~newimp, family="binomial", lc3imp)              
summary(gl3)

################################################################################

imp4 <- read.table("pdata4.txt", header=TRUE)

all4 <- cbind(indicatoren, imp4)

kort1 <- ddply(all4,.(CM033Q01S,CM474Q01S,CM155Q01S,CM155Q04S,
                      CM411Q01S,CM411Q02S,CM803Q01S,CM442Q02S,
                      DM462Q01C,CM034Q01S,CM305Q01S,CM496Q01S,
                      CM496Q02S,CM423Q01S,DM406Q01C,DM406Q02C,
                      CM603Q01S,CM571Q01S,CM564Q01S,CM564Q02S,
                      X1,X2), summarise, W1=sum(imp==1),
               W2=sum(imp==2), N=(W1+W2))

for (k in 1:nrow(kort1))
{
  kort1[k,"PY"] <- kort1[k,"N"]/sum(kort1[ ,"N"], na.rm=T)                  
  kort1[k,"PW1"] <- kort1[k,"W1"]/sum(kort1[k,c("W1","W2")])             
  kort1[k,"PW2"] <- kort1[k,"W2"]/sum(kort1[k,c("W1","W2")])
  # kans op W gegeven X
  kort1[k,"W1X1"] <- (kort1[k,"PY"]*kort1[k,paste0("X1")]*kort1[k,"PW1"])
  kort1[k,"W1X2"] <- (kort1[k,"PY"]*kort1[k,paste0("X2")]*kort1[k,"PW1"])
  kort1[k,"W2X1"] <- (kort1[k,"PY"]*kort1[k,paste0("X1")]*kort1[k,"PW2"])
  kort1[k,"W2X2"] <- (kort1[k,"PY"]*kort1[k,paste0("X2")]*kort1[k,"PW2"])
  
  kort1[k,"PX1"] <- kort1[k,paste0("X1")]*(kort1[k,"N"]/sum(kort1[,"N"], na.rm=T))
  kort1[k,"PX2"] <- kort1[k,paste0("X2")]*(kort1[k,"N"]/sum(kort1[,"N"], na.rm=T))
}

class_mat1 <- matrix(NA,2,2)                                                     
class_mat1[1,1] <- sum(kort1[,"W1X1"], na.rm=T)/sum(kort1[,"PX1"], na.rm=T)
class_mat1[1,2] <- sum(kort1[,"W2X1"], na.rm=T)/sum(kort1[,"PX1"], na.rm=T)

class_mat1[2,1] <- sum(kort1[,"W1X2"], na.rm=T)/sum(kort1[,"PX2"], na.rm=T)
class_mat1[2,2] <- sum(kort1[,"W2X2"], na.rm=T)/sum(kort1[,"PX2"], na.rm=T)

colnames(class_mat1) <- c("W1","W2")
rownames(class_mat1) <- c("X1","X2")

class_mat1


outfiles           <- paste0("6_model.",1:5)                          
bootstrap.colnames <- c("b1","b2","b3","b4","b5")

envir <- new.env()                                                             
run.template(template.path="6_ML4.brew", envir=envir, 
             temp.filename.base=outfiles[4], 
             bootstrap.colnames=bootstrap.colnames[1])

lc4 <- read.delim("6_ML4.txt", header=T)
lc4imp <- lc4[,c("imp","ST006Q04TA","theta.1","theta.2")]
lc4imp$newimp <- NULL

for (i in 1:nrow(lc4imp)){
  lc4imp[i,"newimp"] <- which(rmultinom(1,1,lc4imp[i,c("theta.2","theta.1")])==1)
}

lc4imp$newimp <- as.factor(lc4imp$newimp)
gl4 <- glm(formula=ST006Q04TA~newimp, family="binomial", lc4imp)              
summary(gl4)

################################################################################

imp5 <- read.table("pdata5.txt", header=TRUE)

all5 <- cbind(indicatoren, imp5)

kort1 <- ddply(all5,.(CM033Q01S,CM474Q01S,CM155Q01S,CM155Q04S,
                      CM411Q01S,CM411Q02S,CM803Q01S,CM442Q02S,
                      DM462Q01C,CM034Q01S,CM305Q01S,CM496Q01S,
                      CM496Q02S,CM423Q01S,DM406Q01C,DM406Q02C,
                      CM603Q01S,CM571Q01S,CM564Q01S,CM564Q02S,
                      X2,X1), summarise, W1=sum(imp==1),
               W2=sum(imp==2), N=(W1+W2))

for (k in 1:nrow(kort1))
{
  kort1[k,"PY"] <- kort1[k,"N"]/sum(kort1[ ,"N"], na.rm=T)                    
  kort1[k,"PW1"] <- kort1[k,"W1"]/sum(kort1[k,c("W1","W2")])          
  kort1[k,"PW2"] <- kort1[k,"W2"]/sum(kort1[k,c("W1","W2")])
  # kans op W gegeven X
  kort1[k,"W1X1"] <- (kort1[k,"PY"]*kort1[k,paste0("X1")]*kort1[k,"PW1"])
  kort1[k,"W1X2"] <- (kort1[k,"PY"]*kort1[k,paste0("X2")]*kort1[k,"PW1"])
  kort1[k,"W2X1"] <- (kort1[k,"PY"]*kort1[k,paste0("X1")]*kort1[k,"PW2"])
  kort1[k,"W2X2"] <- (kort1[k,"PY"]*kort1[k,paste0("X2")]*kort1[k,"PW2"])
  
  kort1[k,"PX1"] <- kort1[k,paste0("X1")]*(kort1[k,"N"]/sum(kort1[,"N"], na.rm=T))
  kort1[k,"PX2"] <- kort1[k,paste0("X2")]*(kort1[k,"N"]/sum(kort1[,"N"], na.rm=T))
}

class_mat1 <- matrix(NA,2,2)                                                     
class_mat1[1,1] <- sum(kort1[,"W1X1"], na.rm=T)/sum(kort1[,"PX1"], na.rm=T)
class_mat1[1,2] <- sum(kort1[,"W2X1"], na.rm=T)/sum(kort1[,"PX1"], na.rm=T)

class_mat1[2,1] <- sum(kort1[,"W1X2"], na.rm=T)/sum(kort1[,"PX2"], na.rm=T)
class_mat1[2,2] <- sum(kort1[,"W2X2"], na.rm=T)/sum(kort1[,"PX2"], na.rm=T)

colnames(class_mat1) <- c("W1","W2")
rownames(class_mat1) <- c("X1","X2")

class_mat1


outfiles           <- paste0("6_model.",1:5)                         
bootstrap.colnames <- c("b1","b2","b3","b4","b5")

envir <- new.env()                                                             
run.template(template.path="6_ML5.brew", envir=envir, 
             temp.filename.base=outfiles[5], 
             bootstrap.colnames=bootstrap.colnames[1])

lc5 <- read.delim("6_ML5.txt", header=T)
lc5imp <- lc5[,c("imp","ST006Q04TA","theta.1","theta.2")]
lc5imp$newimp <- NULL

for (i in 1:nrow(lc5imp)){
  lc5imp[i,"newimp"] <- which(rmultinom(1,1,lc5imp[i,c("theta.1","theta.2")])==1)
}

lc5imp$newimp <- as.factor(lc5imp$newimp)
gl5 <- glm(formula=ST006Q04TA~newimp, family="binomial", lc5imp)                
summary(gl5)

################################################################################
# poolen
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

ll1 <- Qbar1 - qt(.975, (nrow(lc1imp)-1)-1)*sqrtT1
ll2 <- Qbar2 - qt(.975, (nrow(lc1imp)-1)-1)*sqrtT2

ul1 <- Qbar1 + qt(.975, (nrow(lc1imp)-1)-1)*sqrtT1
ul2 <- Qbar2 + qt(.975, (nrow(lc1imp)-1)-1)*sqrtT2

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
