createdata <- function(scen.mat,nsim, brew){
  
  for(i in 1:nrow(scen.mat)){                                                   # loop over all simulation conditions
    set.seed(123)                                                               # start every loop with this seed
    cat(i)                                                                      # count the loops
    
    folder.name <- paste0("F:/datapackage_threestepMILC/1_simulation_study/1_createdata/pops_",                 
                          scen.mat[i,1],"_",                                    # create the folder name to put
                          substring(scen.mat[i,2], 3),"_",                      # the files in we are going to
                          substring(scen.mat[i,3], 3),"_",                      # create
                          substring(scen.mat[i,4], 3))
    setwd(folder.name)                                                          # and set to this folder

    create.name    <- paste0("create_ss",scen.mat$ss[i],".txt")                 # create an imput file for latent gold
    create.content <- c("y1 y2 y3 q z freq", paste("0 0 0 0 0",scen.mat$ss[i])) # with variable names and 
    create.file    <- file(create.name)                                         # sample size
    writeLines(create.content, create.file)
    close(create.file)  
    
    # Create names for lg and txt files 
    lg.name   <- paste0("LG_",scen.mat[i,1],"_",                                # we also create the names for the 
                        substring(scen.mat[i,2], 3),"_",                        # .lgs and .txt files
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3),
                        ".lgs")
    
    data.name <- paste0("dat_",scen.mat[i,1],"_",
                        substring(scen.mat[i,2], 3),"_",
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3),
                        ".txt")
    
    in.file   <- paste0("F:/datapackage_threestepMILC/1_simulation_study/1_createdata/pops_",
                        scen.mat[i,1],"_",                                      # and create the corresponding files
                        substring(scen.mat[i,2], 3),"_",
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3),"/",
                        create.name)
    
    out.file  <- paste0("F:/datapackage_threestepMILC/1_simulation_study/1_createdata/pops_",
                        scen.mat[i,1],"_",
                        substring(scen.mat[i,2], 3),"_",
                        substring(scen.mat[i,3], 3),"_",
                        substring(scen.mat[i,4], 3),"/",
                        data.name)

    y.min  <- log((1-scen.mat[i,2])/scen.mat[i,2])                              # specify the coefficient values
    y.plus <- log(scen.mat[i,2]/(1-scen.mat[i,2]))                              # we need for the LG syntax
    z.coef <- log((2*scen.mat[i,3])/(1-2*scen.mat[i,3]))
    q.coef <- log(scen.mat[i,4]/(1-scen.mat[i,4]))
    
    envir       <- new.env()                                                    # next we run the template which
    createdata.template(template.path      = brew,                              # creates the data samples
                        envir              = envir,                             # for all iterations of the 
                        temp.filename.base = lg.name,                           # simulation
                        in.file            = in.file, 
                        out.file           = out.file, 
                        nsim               = nsim, 
                        y.min              = y.min, 
                        y.plus             = y.plus, 
                        z.coef             = z.coef, 
                        q.coef             = q.coef)
    
    setwd("F:/datapackage_threestepMILC/1_simulation_study/1_createdata")       # and leave the folder of the 
  }                                                                             # simulation condition to go
}                                                                               # to the next one