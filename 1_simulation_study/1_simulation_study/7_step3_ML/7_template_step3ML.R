template_step3<- function(template.path, envir, temp.filename.base, s, j , bcname,
                         a, qcoef, zcoef, covar) {                              # this template is used to run latent
                                                                                # gold from R. 
  temp.filename <- paste(temp.filename.base, ".lgs", sep="")                          
  template      <- file(template.path, 'r')
  brew(template, output=temp.filename) 
  close(template)
  shell(sprintf('"C:/Users/Documents/LatentGOLD5.0/lg50.exe" %s', temp.filename))
}                                                                               # specify with shell where the latent gold
                                                                                # program file is located on your pc