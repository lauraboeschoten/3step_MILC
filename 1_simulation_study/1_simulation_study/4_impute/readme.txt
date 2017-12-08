The function 'create_folders.R' created folders for every condition in the simulation. 

You run the code in '4_TOP_impute.R', here the code from '4_createprofile.R', '4_posteriors.R'
and '4_impute.R' is sourced. When you are finished, all folders contain 1000 files, with 
imputations, Each file in named 'imputations_<nr1-1000>.txt' and contains the corresponding
sample dataset in long format and 5 imputations. The 1000 files with names 
'posteriors<nr1-1000>.txt' contain for each sample dataset in short format the posterior 
membership probabilities of being in Latent Class 1 obtained using information from the 
5 latent class models. The 5000 files with names 'profile<nr1-1000><nr1-5>.txt' contain
the conditional probabilities for every bootstrap sample in every iteration which are used
to calculate the posterior membership probabilities. 