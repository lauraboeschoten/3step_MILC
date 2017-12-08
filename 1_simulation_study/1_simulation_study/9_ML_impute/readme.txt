The function 'create_folders.R' created folders for every condition in the simulation. 

Next, the function in '9_TOP_impute3' is opened, and the function '9_fun_impute3.R' is 
sourced here. These functions create imputations based on the 3-step latent class model 
that was run in '7_step3_ML'. In '8_LC_results_ML1', the logit coefficients were directly
obtained from the LG output, while we now create imputations based on these latent class
models. 