The function 'create_folders.R' created folders for every condition in the simulation. 

The next step is to open and run the code from '8_TOP_LC_results_ML1.R'. Here, the 
functions from '8_results.R' and '8_mean_results.R' are sourced. 

When you are finished, all folders contain 1000 files with 
results for the relation between Q and W with names 'results_q_<nr1-1000>.txt'. The 
folders also contain 1000 files with results of the crosstable between Z and Q with names
'results_z_<1-1000>.txt'. All folders also contain the files 'finalresult_q.txt' and
'finalresult_z.txt' with average solutions over the 1000 simulation iterations. At last,
all folders contain a file 'nrmis.txt' which states if there were any inadmissable solutions
produced and is so, how many. Inadmissable solutions are generally obtained due to 
boundary issues. 

After running the whole procedure for every simulation condition, a check by hand was done to 
see if nonconvergence occured for any of the iterations in any of the conditions. To do this,
the results created by '5_results.R' were checked with a piece of code 
'which(results[i][2,2]>20)'. Here, we see which iteration resulted in an exceptionally 
large estimate for the variance, which can be an indicator of nonconvergence. In the folder
'7_step3_ML', the file 'extra_iterations_overview.txt' contains the number of iterations which
possibly dealt with nonconvergence. The latent class models for these conditions were
run again using 5000 instead of 500 iterations. Afterwards, the results were obtained again
using the R-files in this folder.