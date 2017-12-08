posterior.probabilities <- function(conditionals) {
  
  v1 <- v2 <- v3 <- v4 <- v5 <- c(1,2)                                          # create variables
  combinations           <- expand.grid(v1=v1,v2=v2,v3=v3,v4=v4,v5=v5)          # make a grid with all combinations of 
  combinations           <- combinations[, rev(seq_len(ncol(combinations)))]    # scores on the two variables
  colnames(combinations) <- c("q","z","y1","y2","y3")                           # give the variables names
  
  combinations$posterior.5   <- combinations$posterior.4 <-                     # create empty variables for posteriors
    combinations$posterior.3 <- combinations$posterior.2 <- 
    combinations$posterior.1 <- NA
  
  for (a in 1:nboot) {                                                          # for every bootstrap sample
    for (b in 1:ncomb) {                                                        # for every combination of scores
      pr1a  <- conditionals[[a]][1,1]                                           # prior is overall class size
      if(combinations[b,1]==1) {q1 <- 1} else {q1 <- 1}                         # grab relevant conditionals from profile
      if(combinations[b,2]==1) {z1 <- 1} else {z1 <- 1}
      if(combinations[b,3]==1) {y1a <- conditionals[[a]][2,1]} else {y1a <- conditionals[[a]][3,1]}  
      if(combinations[b,4]==1) {y2a <- conditionals[[a]][4,1]} else {y2a <- conditionals[[a]][5,1]}  
      if(combinations[b,5]==1) {y3a <- conditionals[[a]][6,1]} else {y3a <- conditionals[[a]][7,1]}  
      pr1b  <- conditionals[[a]][1,2]
      if(combinations[b,1]==1) {q2 <- 1} else {q2 <- 1}
      if(combinations[b,2]==1) {z2 <- 1} else {z2 <- 1}
      if(combinations[b,3]==1) {y1b <- conditionals[[a]][2,2]} else {y1b <- conditionals[[a]][3,2]}  
      if(combinations[b,4]==1) {y2b <- conditionals[[a]][4,2]} else {y2b <- conditionals[[a]][5,2]}  
      if(combinations[b,5]==1) {y3b <- conditionals[[a]][6,2]} else {y3b <- conditionals[[a]][7,2]} 
      
      combinations[b,a+5] <- (pr1a*q1*z1*y1a*y2a*y3a)/((pr1a*q1*z1*y1a*y2a*y3a)+(pr1b*q2*z2*y1b*y2b*y3b))
    }                                                                           # apply the posterior formula to all
  }                                                                             # combinations
  combinations[is.na(combinations)] <- 0                                        # if in both clusters the probability is 0
  return(combinations)                                                          # NaN's are created, then we replace this
}                                                                               # by 0. 
