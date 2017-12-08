template_LC <- function(template.path, envir, temp.filename.base,s,j,           # this is the createdata.template
                        bcname, y.min, y.plus) {                                # which is used to run Latent Gold 
                                                                                # from R
  
  temp.filename <- paste(temp.filename.base, ".lgs", sep="")                          
  template      <- file(template.path, 'r')
  brew(template, output=temp.filename) 
  close(template)
  shell(sprintf('"C:/Users/Documents/LatentGOLD5.0/lg50.exe" %s', temp.filename))
}                                                                               # specify with shell where the Latent
                                                                                # Gold program file is located