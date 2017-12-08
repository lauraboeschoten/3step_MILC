setwd("F:datapackage_threestepMILC/5_PISA")

library(foreign)
library(plyr)
library(brew)  

options(scipen=999)
set.seed(123)

data1 <- read.spss("NL_kids.sav", to.data.frame = TRUE)
data2 <- read.spss("PISA_NL_book43.sav", to.data.frame = TRUE)

data  <- merge(data1, data2, by="CNTSTUID")                                     # combine data

setwd("F:datapackage_threestepMILC/5_PISA/included")

gdata <- NULL    
for (i in 1:ncol(data)){                                                        # remove spaces 
  gdata <- cbind(gdata,gsub(" ","",data[,i]))
}

gdata        <- data.frame(gdata)
names(gdata) <- names(data)

short_data   <- ddply(gdata,.(CM033Q01S,CM474Q01S,CM155Q01S,CM155Q04S,CM411Q01S,# create short data
                              CM411Q02S,CM803Q01S,CM442Q02S,DM462Q01C,CM034Q01S, 
                              CM305Q01S,CM496Q01S,CM496Q02S,CM423Q01S,DM406Q01C,
                              DM406Q02C,CM603Q01S,CM571Q01S,CM564Q01S,CM564Q02S,
                              ST006Q04TA,ST008Q04TA),nrow)

boots           <- rmultinom(5, nrow(gdata), short_data$V1/sum(short_data$V1))  # bootstrap
colnames(boots) <- c("b1","b2","b3","b4","b5")                                 
bdata           <- cbind(short_data, boots)

write.table(bdata, "pisa.txt", row.names=FALSE, quote=FALSE)

##############################################################
                                                                
source("runtemplate.R")                                         

outfiles           <- paste0("3_model.",1:5)                                  
bootstrap.colnames <- tail(colnames(bdata),5)

##########

envir <- new.env()                                                              # run LC models
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

scan.mod        <- readLines("3_model.1.lst")                                   # read in LC output
lin.size        <- grep('Size',scan.mod)              
overall         <- matrix(NA, 1, 2)
overall[1,1:2]  <- strsplit(scan.mod[lin.size], "\t")[[1]][c(2,4)]
overall         <- as.numeric(overall)

a <- grep('dependent',scan.mod)[[2]]                                            
b <- grep('independent',scan.mod)[[1]]                                          

vars <- unlist(strsplit(scan.mod[a:(b-1)], " "))[c(6,8,10,12,14,16,             
                                                   32,34,36,38,40,42,
                                                   58,60,62,64,66,68,
                                                   84,86)]

pr             <- grep('Profile',scan.mod)                                      # obtain profile                          
matlist        <- vector("list", length(vars))                                   
names(matlist) <- vars                                                          

for (i in 1:length(vars)){
  c <- grep(vars[i],scan.mod)[[which(grep(vars[i],scan.mod) > pr)[1]]]          # for indicators
  if(i!=length(vars)){
    d <- grep(vars[i+1],scan.mod)[[which(grep(vars[i+1],scan.mod) > pr)[1]]]    
  } else {d <- grep('ProbMeans-Posterior',scan.mod)}                            
  matlist[[i]] <- matrix(unlist(strsplit(scan.mod[(c+1):(d-1)], "\t")),         
                         ncol=7, byrow=T)[,c(1,2,4)]
}

varsco <- unlist(strsplit(scan.mod[b+1], " "))[c(14,16)]                        # and covariates
matco        <- vector("list", length(varsco))                                  
names(matco) <- varsco
b2 <- grep('Independent',scan.mod)[[2]]  

for (i in 1:length(varsco)){                                                  
  c <- grep(varsco[i],scan.mod)[[which(grep(varsco[i],scan.mod) > b2)[1]]]      
  if(i!=length(varsco)){                                                      
    d <- grep(varsco[i+1],scan.mod)[[which(grep(varsco[i+1],scan.mod) > b2)[1]]]# put in a list 
  } else {d <- grep('V1',scan.mod)[[3]]}
  matco[[i]] <- matrix(unlist(strsplit(scan.mod[(c+1):(d-1)], "\t")), 
                       ncol=3, byrow=T)
}

for(i in 1:length(matco)){                                                       
  for(j in 2:3){                                                                
    matco[[i]][,j] <- as.numeric(matco[[i]][,j])/sum(as.numeric(matco[[i]][,j]))
  }
}

matlist <- c(matlist, matco)                                                    # combine indicators and covariates

for(i in 1:length(matlist)){                                                    # rename "NA"'s 
  matlist[[i]][matlist[[i]]=="NA"] <- NA                                      
}

ar <- array(NA, dim=c(nrow(bdata),length(matlist),2))                            

for(j in 1:nrow(bdata)){                                                        # fill in all conditionals
  for(i in 1:length(matlist)){                                                  # make numeric
    kkk <- bdata[j,i]                                                           
    m.k <- match(kkk,matlist[[i]][,1])                                          
    if(is.na(m.k)){                                                             # if NA due to bootstrap, 
      ar[j,i,] <- c(1,1)                                                        # fill in 1 
    } else {
      ar[j,i,] <- as.numeric(matlist[[i]][m.k,2:3])
    }}
}

post <- matrix(NA, nrow(bdata), 2)                                              # matrix for posteriors

for(j in 1:nrow(post)){                                                         # create posteriors
  for (i in 1:2){                                                               
    post[j,i] <- prod(overall[i],ar[j,,i], na.rm = FALSE)/
      sum(prod(overall[1],ar[j,,1]),
          prod(overall[2],ar[j,,2]))
  }
}

post[is.na(post)] <- 0.5                                                        # solution for NaN
pdata1 <- cbind(bdata[,c("ST006Q04TA","ST008Q04TA","V1")],post)           

pdata1$imp <- NULL                                                              # sample from posteriors
for (i in 1:nrow(pdata1)){                                                    
  pdata1[i,"imp"] <- which(rmultinom(1,1,pdata1[i,c("2","1")])==1)              # check if there is label switching! 
}                                                                               # this differs for every LC model! 

write.table(pdata1, "pdata1.txt", quote=FALSE, sep=" ", row.names=FALSE)

pdata1 <- read.table("pdata1.txt", header=TRUE)

pdata1$imp <- as.factor(pdata1$imp)
lm1 <- glm(ST006Q04TA~imp, family="binomial", pdata1)
summary(lm1)

################################################################################

scan.mod        <- readLines("3_model.2.lst")                                   
lin.size        <- grep('Size',scan.mod)                                        
overall         <- matrix(NA, 1, 2)
overall[1,1:2]  <- strsplit(scan.mod[lin.size], "\t")[[1]][c(2,4)]
overall         <- as.numeric(overall)

a <- grep('dependent',scan.mod)[[2]]                                            
b <- grep('independent',scan.mod)[[1]]                                          

vars <- unlist(strsplit(scan.mod[a:(b-1)], " "))[c(6,8,10,12,14,16,             
                                                   32,34,36,38,40,42,
                                                   58,60,62,64,66,68,
                                                   84,86)]

pr             <- grep('Profile',scan.mod)                                      
matlist        <- vector("list", length(vars))                                   
names(matlist) <- vars                                                          

for (i in 1:length(vars)){
  c <- grep(vars[i],scan.mod)[[which(grep(vars[i],scan.mod) > pr)[1]]]          
  if(i!=length(vars)){
    d <- grep(vars[i+1],scan.mod)[[which(grep(vars[i+1],scan.mod) > pr)[1]]]    
  } else {d <- grep('ProbMeans-Posterior',scan.mod)}                            
  matlist[[i]] <- matrix(unlist(strsplit(scan.mod[(c+1):(d-1)], "\t")),         
                         ncol=7, byrow=T)[,c(1,2,4)]
}

varsco <- unlist(strsplit(scan.mod[b+1], " "))[c(14,16)]                     
matco        <- vector("list", length(varsco))                               
names(matco) <- varsco
b2 <- grep('Independent',scan.mod)[[2]]  

for (i in 1:length(varsco)){                                                 
  c <- grep(varsco[i],scan.mod)[[which(grep(varsco[i],scan.mod) > b2)[1]]]   
  if(i!=length(varsco)){                                                     
    d <- grep(varsco[i+1],scan.mod)[[which(grep(varsco[i+1],scan.mod) > b2)[1]]]
  } else {d <- grep('V1',scan.mod)[[3]]}
  matco[[i]] <- matrix(unlist(strsplit(scan.mod[(c+1):(d-1)], "\t")), 
                       ncol=3, byrow=T)
}

for(i in 1:length(matco)){                                                       
  for(j in 2:3){                                                                
    matco[[i]][,j] <- as.numeric(matco[[i]][,j])/sum(as.numeric(matco[[i]][,j]))
  }
}

matlist <- c(matlist, matco)                                                    

for(i in 1:length(matlist)){                                                     
  matlist[[i]][matlist[[i]]=="NA"] <- NA                                        
}

ar <- array(NA, dim=c(nrow(bdata),length(matlist),2))                            

for(j in 1:nrow(bdata)){                                                        
  for(i in 1:length(matlist)){                                                  
    kkk <- bdata[j,i]                                                           
    m.k <- match(kkk,matlist[[i]][,1])                                          
    if(is.na(m.k)){
      ar[j,i,] <- c(1,1)
    } else {
      ar[j,i,] <- as.numeric(matlist[[i]][m.k,2:3])
    }}
}

post <- matrix(NA, nrow(bdata), 2)                                              

for(j in 1:nrow(post)){                                                          
  for (i in 1:2){                                                               
    post[j,i] <- prod(overall[i],ar[j,,i], na.rm = FALSE)/
      sum(prod(overall[1],ar[j,,1]),
          prod(overall[2],ar[j,,2]))
  }
}

post[is.na(post)] <- 0.5 
pdata2 <- cbind(bdata[,c("ST006Q04TA","ST008Q04TA","V1")],post)   
pdata2$imp <- NULL                                                 
for (i in 1:nrow(pdata1)){                                         
  pdata2[i,"imp"] <- which(rmultinom(1,1,pdata2[i,c("2","1")])==1) 
}

write.table(pdata2, "pdata2.txt", quote=FALSE, sep=" ", row.names=FALSE)

pdata2 <- read.table("pdata2.txt", header=TRUE)

pdata2$imp <- as.factor(pdata2$imp)
lm2 <- glm(ST006Q04TA~imp, family="binomial", pdata2)
summary(lm2)

################################################################################

scan.mod        <- readLines("3_model.3.lst")                                 
lin.size        <- grep('Size',scan.mod)                                      
overall         <- matrix(NA, 1, 2)
overall[1,1:2]  <- strsplit(scan.mod[lin.size], "\t")[[1]][c(2,4)]
overall         <- as.numeric(overall)

a <- grep('dependent',scan.mod)[[2]]                                          
b <- grep('independent',scan.mod)[[1]]                                

vars <- unlist(strsplit(scan.mod[a:(b-1)], " "))[c(6,8,10,12,14,16,           
                                                   32,34,36,38,40,42,
                                                   58,60,62,64,66,68,
                                                   84,86)]

pr             <- grep('Profile',scan.mod)                                      
matlist        <- vector("list", length(vars))                                  
names(matlist) <- vars                                                      

for (i in 1:length(vars)){
  c <- grep(vars[i],scan.mod)[[which(grep(vars[i],scan.mod) > pr)[1]]]        
  if(i!=length(vars)){
    d <- grep(vars[i+1],scan.mod)[[which(grep(vars[i+1],scan.mod) > pr)[1]]
  } else {d <- grep('ProbMeans-Posterior',scan.mod)}                          
  matlist[[i]] <- matrix(unlist(strsplit(scan.mod[(c+1):(d-1)], "\t")),       
                         ncol=7, byrow=T)[,c(1,2,4)]
}

varsco <- unlist(strsplit(scan.mod[b+1], " "))[c(14,16)]                  
matco        <- vector("list", length(varsco))                            
names(matco) <- varsco
b2 <- grep('Independent',scan.mod)[[2]]  

for (i in 1:length(varsco)){                                                    
  c <- grep(varsco[i],scan.mod)[[which(grep(varsco[i],scan.mod) > b2)[1]]]  
  if(i!=length(varsco)){                                                      
    d <- grep(varsco[i+1],scan.mod)[[which(grep(varsco[i+1],scan.mod) > b2)[1]]]
  } else {d <- grep('V1',scan.mod)[[3]]}
  matco[[i]] <- matrix(unlist(strsplit(scan.mod[(c+1):(d-1)], "\t")), 
                       ncol=3, byrow=T)
}

for(i in 1:length(matco)){                                                      
  for(j in 2:3){                                                          
    matco[[i]][,j] <- as.numeric(matco[[i]][,j])/sum(as.numeric(matco[[i]][,j]))
  }
}

matlist <- c(matlist, matco)                                              

for(i in 1:length(matlist)){                                            
  matlist[[i]][matlist[[i]]=="NA"] <- NA                              
}

ar <- array(NA, dim=c(nrow(bdata),length(matlist),2))                       

for(j in 1:nrow(bdata)){                                                       
  for(i in 1:length(matlist)){                                                 
    kkk <- bdata[j,i]                                                           
    m.k <- match(kkk,matlist[[i]][,1])                                        
    if(is.na(m.k)){
      ar[j,i,] <- c(1,1)
    } else {
      ar[j,i,] <- as.numeric(matlist[[i]][m.k,2:3])
    }}
}

post <- matrix(NA, nrow(bdata), 2)                   

for(j in 1:nrow(post)){                                                         
  for (i in 1:2){                                                               
    post[j,i] <- prod(overall[i],ar[j,,i], na.rm = FALSE)/
      sum(prod(overall[1],ar[j,,1]),
          prod(overall[2],ar[j,,2]))
  }
}

post[is.na(post)] <- 0.5 
pdata3 <- cbind(bdata[,c("ST006Q04TA","ST008Q04TA","V1")],post)     
# scores
pdata3$imp <- NULL                                                              
for (i in 1:nrow(pdata1)){                                                      
  pdata3[i,"imp"] <- which(rmultinom(1,1,pdata3[i,c("2","1")])==1)        
}

write.table(pdata3, "pdata3.txt", quote=FALSE, sep=" ", row.names=FALSE)
pdata3 <- read.table("pdata3.txt", header=TRUE)
pdata3$imp <- as.factor(pdata3$imp)
lm3 <- glm(ST006Q04TA~imp, family="binomial", pdata3)
summary(lm3)
 
################################################################################

scan.mod        <- readLines("3_model.4.lst")                                   
lin.size        <- grep('Size',scan.mod)                                         
overall         <- matrix(NA, 1, 2)
overall[1,1:2]  <- strsplit(scan.mod[lin.size], "\t")[[1]][c(2,4)]
overall         <- as.numeric(overall)

a <- grep('dependent',scan.mod)[[2]]                                            
b <- grep('independent',scan.mod)[[1]]                                          

vars <- unlist(strsplit(scan.mod[a:(b-1)], " "))[c(6,8,10,12,14,16,             
                                                   32,34,36,38,40,42,
                                                   58,60,62,64,66,68,
                                                   84,86)]

pr             <- grep('Profile',scan.mod)                                      
matlist        <- vector("list", length(vars))                                   
names(matlist) <- vars                                                          

for (i in 1:length(vars)){
  c <- grep(vars[i],scan.mod)[[which(grep(vars[i],scan.mod) > pr)[1]]]          
  if(i!=length(vars)){
    d <- grep(vars[i+1],scan.mod)[[which(grep(vars[i+1],scan.mod) > pr)[1]]]    
  } else {d <- grep('ProbMeans-Posterior',scan.mod)}                            
  matlist[[i]] <- matrix(unlist(strsplit(scan.mod[(c+1):(d-1)], "\t")),         
                         ncol=7, byrow=T)[,c(1,2,4)]
}

varsco <- unlist(strsplit(scan.mod[b+1], " "))[c(14,16)]                     
matco        <- vector("list", length(varsco))                               
names(matco) <- varsco
b2 <- grep('Independent',scan.mod)[[2]]  

for (i in 1:length(varsco)){                                                 
  c <- grep(varsco[i],scan.mod)[[which(grep(varsco[i],scan.mod) > b2)[1]]]   
  if(i!=length(varsco)){                                                     
    d <- grep(varsco[i+1],scan.mod)[[which(grep(varsco[i+1],scan.mod) > b2)[1]]]
  } else {d <- grep('V1',scan.mod)[[3]]}
  matco[[i]] <- matrix(unlist(strsplit(scan.mod[(c+1):(d-1)], "\t")), 
                       ncol=3, byrow=T)
}

for(i in 1:length(matco)){                                                      
  for(j in 2:3){                                                                
    matco[[i]][,j] <- as.numeric(matco[[i]][,j])/sum(as.numeric(matco[[i]][,j]))
  }
}

matlist <- c(matlist, matco)                                                    

for(i in 1:length(matlist)){                                              
  matlist[[i]][matlist[[i]]=="NA"] <- NA                                  
}

ar <- array(NA, dim=c(nrow(bdata),length(matlist),2))                      

for(j in 1:nrow(bdata)){                                                  
  for(i in 1:length(matlist)){                                            
    kkk <- bdata[j,i]                                                     
    m.k <- match(kkk,matlist[[i]][,1])                                    
    if(is.na(m.k)){
      ar[j,i,] <- c(1,1)
    } else {
      ar[j,i,] <- as.numeric(matlist[[i]][m.k,2:3])
    }}
}

post <- matrix(NA, nrow(bdata), 2)                                        

for(j in 1:nrow(post)){                                                    
  for (i in 1:2){                                                         
    post[j,i] <- prod(overall[i],ar[j,,i], na.rm = FALSE)/
      sum(prod(overall[1],ar[j,,1]),
          prod(overall[2],ar[j,,2]))
  }
}

post[is.na(post)] <- 0.5 
pdata4 <- cbind(bdata[,c("ST006Q04TA","ST008Q04TA","V1")],post)   
# scores

pdata4$imp <- NULL                                                              
for (i in 1:nrow(pdata4)){                                                    
  pdata4[i,"imp"] <- which(rmultinom(1,1,pdata4[i,c("2","1")])==1)        
}

write.table(pdata4, "pdata4.txt", quote=FALSE, sep=" ", row.names=FALSE)

pdata4 <- read.table("pdata4.txt", header=TRUE)

pdata4$imp <- as.factor(pdata4$imp)
lm4 <- glm(ST006Q04TA~imp, family="binomial", pdata4)
summary(lm4)

################################################################################

scan.mod        <- readLines("3_model.5.lst")                                   
lin.size        <- grep('Size',scan.mod)                                        
overall         <- matrix(NA, 1, 2)
overall[1,1:2]  <- strsplit(scan.mod[lin.size], "\t")[[1]][c(2,4)]
overall         <- as.numeric(overall)

a <- grep('dependent',scan.mod)[[2]]                                            
b <- grep('independent',scan.mod)[[1]]                                          

vars <- unlist(strsplit(scan.mod[a:(b-1)], " "))[c(6,8,10,12,14,16,             
                                                   32,34,36,38,40,42,
                                                   58,60,62,64,66,68,
                                                   84,86)]

pr             <- grep('Profile',scan.mod)                                      
matlist        <- vector("list", length(vars))                                  
names(matlist) <- vars                                                          

for (i in 1:length(vars)){
  c <- grep(vars[i],scan.mod)[[which(grep(vars[i],scan.mod) > pr)[1]]]          
  if(i!=length(vars)){
    d <- grep(vars[i+1],scan.mod)[[which(grep(vars[i+1],scan.mod) > pr)[1]]]    
  } else {d <- grep('ProbMeans-Posterior',scan.mod)}                            
  matlist[[i]] <- matrix(unlist(strsplit(scan.mod[(c+1):(d-1)], "\t")),         
                         ncol=7, byrow=T)[,c(1,2,4)]
}

varsco <- unlist(strsplit(scan.mod[b+1], " "))[c(14,16)]                     
matco        <- vector("list", length(varsco))                               
names(matco) <- varsco
b2 <- grep('Independent',scan.mod)[[2]]  

for (i in 1:length(varsco)){                                                    
  c <- grep(varsco[i],scan.mod)[[which(grep(varsco[i],scan.mod) > b2)[1]]]      
  if(i!=length(varsco)){                                                        
    d <- grep(varsco[i+1],scan.mod)[[which(grep(varsco[i+1],scan.mod) > b2)[1]]]
  } else {d <- grep('V1',scan.mod)[[3]]}
  matco[[i]] <- matrix(unlist(strsplit(scan.mod[(c+1):(d-1)], "\t")), 
                       ncol=3, byrow=T)
}

for(i in 1:length(matco)){                                                      
  for(j in 2:3){                                                                
    matco[[i]][,j] <- as.numeric(matco[[i]][,j])/sum(as.numeric(matco[[i]][,j]))
  }
}

matlist <- c(matlist, matco)                                                  

for(i in 1:length(matlist)){                                                     
  matlist[[i]][matlist[[i]]=="NA"] <- NA                                        
}

ar <- array(NA, dim=c(nrow(bdata),length(matlist),2))                            

for(j in 1:nrow(bdata)){                                                        
  for(i in 1:length(matlist)){                                                  
    kkk <- bdata[j,i]                                                           
    m.k <- match(kkk,matlist[[i]][,1])                                          
    if(is.na(m.k)){
      ar[j,i,] <- c(1,1)
    } else {
      ar[j,i,] <- as.numeric(matlist[[i]][m.k,2:3])
    }}
}

post <- matrix(NA, nrow(bdata), 2)                                              

for(j in 1:nrow(post)){                                                          
  for (i in 1:2){                                                               
    post[j,i] <- prod(overall[i],ar[j,,i], na.rm = FALSE)/
      sum(prod(overall[1],ar[j,,1]),
          prod(overall[2],ar[j,,2]))
  }
}

post[is.na(post)] <- 0.5 
pdata5 <- cbind(bdata[,c("ST006Q04TA","ST008Q04TA","V1")],post)     
# scores
pdata5$imp <- NULL                                                             
for (i in 1:nrow(pdata5)){                                                     
  pdata5[i,"imp"] <- which(rmultinom(1,1,pdata5[i,c("1","2")])==1)         
}

write.table(pdata5, "pdata5.txt", quote=FALSE, sep=" ", row.names=FALSE)
pdata5 <- read.table("pdata5.txt", header=TRUE)

pdata5$imp <- as.factor(pdata5$imp)
lm5 <- glm(ST006Q04TA~imp, family="binomial", pdata5)
summary(lm5)

################################################################################

Q1hat1 <- summary(lm1)$coefficients[1,1]
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

ll1 <- Qbar1 - qt(.975, (nrow(pdata1)-1)-1)*sqrtT1
ll2 <- Qbar2 - qt(.975, (nrow(pdata1)-1)-1)*sqrtT2

ul1 <- Qbar1 + qt(.975, (nrow(pdata1)-1)-1)*sqrtT1
ul2 <- Qbar2 + qt(.975, (nrow(pdata1)-1)-1)*sqrtT2

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
