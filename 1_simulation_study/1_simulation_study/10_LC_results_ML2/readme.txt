The function 'create_folders.R' created folders for every condition in the simulation. 

You run the code in '10_TOP_update_results.R', here the code from '10_results.R' and 
'10_mean_results.R' is sourced. When you are finished, all folders contain 1000 files with 
results for the relation between Q and W with names 'results_q_<nr1-1000>.txt'. The 
folders also contain 1000 files with results of the crosstable between Z and Q with names
'results_z_<1-1000>.txt'. All folders also contain the files 'finalresult_q.txt' and
'finalresult_z.txt' with average solutions over the 1000 simulation iterations. At last,
all folders contain a file 'nrmis.txt' which states if there were any inadmissable solutions
produced and is so, how many. Inadmissable solutions are generally obtained due to 
boundary issues. 