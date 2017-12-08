The function 'create_folders.R' created folders for every condition in the simulation. 

The function in '13_TOP_BCH_impute.R' transforms the amatrices produced in step 11 into
posterior membership probabilities of belonging to a specific latent class given a 
combination of scores on the covariates and the initial imputation. Based on these
posterior membership probabilities, imputations are created and stored in the 
folders belonging to the simulation conditions. Each file containts the 5 new imputations
and  the file names are 'newimputations_<iterationnr1-1000>.txt'.