run.template <- function(template.path, envir, temp.filename.base, i, j , bootstrap.colnames,
                         cp.min, cp.plus, imp.comp, co.coef) {
  
  temp.filename <- paste(temp.filename.base, ".lgs", sep="")                          
  template      <- file(template.path, 'r')
  brew(template, output=temp.filename) 
  close(template)
  shell(sprintf('"C:\\Program Files\\LatentGOLD5.0\\lg50.exe" %s', temp.filename))    
  }
