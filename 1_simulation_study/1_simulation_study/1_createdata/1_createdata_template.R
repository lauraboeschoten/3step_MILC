createdata.template <- function(template.path, envir, temp.filename.base,       # this is the createdata.template
                                in.file, out.file, nsim, y.min, y.plus,         # which is used to run Latent Gold 
                                z.coef, q.coef) {                               # from R
  
  temp.filename <- temp.filename.base                         
  template      <- file(template.path, 'r')
  brew(template, output=temp.filename) 
  close(template)
  shell(sprintf('"C:/Users/Documents/LatentGOLD5.0/lg50.exe" %s', temp.filename))
  }                                                                             # specify with shell where the Latent
                                                                                # Gold program file is located