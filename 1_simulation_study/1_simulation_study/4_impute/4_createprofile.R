create.profile <- function(model){
  
  cond.probs   <- matrix(NA, 7, 2)                                              # create a matrix to store the results
  scan.mod     <- readLines(model)                                              # read the LG results file
  lin.size     <- grep('Size',scan.mod)                                         # find the row where 'Size' appears 
  cond.probs[1,1:2]  <- strsplit(scan.mod[lin.size], "\t")[[1]][c(2,4)]         # grab the relevant values from the LG file
  cond.probs[2,1:2]  <- strsplit(scan.mod[lin.size+2], "\t")[[1]][c(2,4)]       # and store in matrix
  cond.probs[3,1:2]  <- strsplit(scan.mod[lin.size+3], "\t")[[1]][c(2,4)]
  cond.probs[4,1:2]  <- strsplit(scan.mod[lin.size+5], "\t")[[1]][c(2,4)]
  cond.probs[5,1:2]  <- strsplit(scan.mod[lin.size+6], "\t")[[1]][c(2,4)]
  cond.probs[6,1:2]  <- strsplit(scan.mod[lin.size+8], "\t")[[1]][c(2,4)]
  cond.probs[7,1:2]  <- strsplit(scan.mod[lin.size+9], "\t")[[1]][c(2,4)]
  
  cond.probs <- data.frame(cl1=as.numeric(cond.probs[,1]),                      # make strings numeric
                           cl2=as.numeric(cond.probs[,2]),
                           row.names=c("size","1y1","1y2","2y1",                # and give relevant names
                                       "2y2","3y1","3y2"))
  return(cond.probs)
}