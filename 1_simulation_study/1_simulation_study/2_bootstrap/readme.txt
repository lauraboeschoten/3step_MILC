The function 'create_folders.R' created folders for every condition in the simulation. 

You run the code in '2_TOP_bootstrap.R', here the code from '2_fun_create_bootstraps.R' is
sourced. After you are finished, all folders will contain 1000 files, one for every iteration
of the simulation. Each file is called 'bootstrap_sample_<nr1-1000>.txt' and contains an
aggregated dataset, with a row for every combination of scores on the variables, and with
the frequency of each combination of scores in the original sample, and with five bootstrap
samples drawn from this. 