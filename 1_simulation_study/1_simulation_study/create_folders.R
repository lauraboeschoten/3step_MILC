

scens <- list(                                                                  # all the simulation conditions
  ss     =c(200, 500, 1000),                                                    # sample sizes
  c.prob =c(0.70, 0.80, 0.90, 0.95, 0.99),                                      # conditional probabilities
  z.coef =c(0.01, 0.05, 0.10, 0.20),                                            # proportions for P(Z=2)
  q.coef =c(0.5, 0.6224593, 0.7310586, 0.8807971))                              # proportions for P(Q=2)

scen.mat <- expand.grid(scens)                                                  # Create a list with full factorial
                                                                                # overview of all simulation 
                                                                                # conditions

for (i in 1:nrow(scen.mat)){                                                    # create for all 14 steps:
  # 1. create data
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/'1_createdata"       # the main direction folder
  subDir <- paste0("pops_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",   # and subfolders for every 
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3)) # simulation condition
  dir.create(file.path(mainDir, subDir))
  
  # 2. bootstrap
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/2_bootstrap"
  subDir <- paste0("samples_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))
  
  # 3. LC models
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/3_LC_models"
  subDir <- paste0("LC_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))
  
  # 4. impute
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/4_impute"
  subDir <- paste0("files_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))
  
  # 5. results before
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/5_results_before"
  subDir <- paste0("results_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))
  
  # 6 classification tables
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/6_classification_tables"
  subDir <- paste0("table_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))  

  # 7 step 3
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/7_step3_ML"
  subDir <- paste0("step3ML_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))  
  
  # 8 LC results step 3
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/8_LC_results_ML1"
  subDir <- paste0("ML1_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))
  
  # 9 impute update
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/9_ML_impute"
  subDir <- paste0("MLimpute_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))
  
  # 10 update results
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/10_LC_results_ML2"
  subDir <- paste0("ML2_results_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))
  
  # 11 BCH correction
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/11_BCH"
  subDir <- paste0("BCH_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))
  
  # 12 step 3 BCH 
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/12_step3_BCH"
  subDir <- paste0("BCH_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))
  
  # 13 impute BCH
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/13_impute_BCH"
  subDir <- paste0("BCH_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))
  
  # 14 results BCH
  mainDir <- "F:/datapackage_threestepMILC/simulationstudy/14_results_BCH"
  subDir <- paste0("BCH_",scen.mat[i,1],"_",substring(scen.mat[i,2], 3),"_",
                   substring(scen.mat[i,3], 3),"_",substring(scen.mat[i,4], 3))
  dir.create(file.path(mainDir, subDir))
  
}
