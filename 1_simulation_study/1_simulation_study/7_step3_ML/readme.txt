The function 'create_folders.R' created folders for every condition in the simulation. 
The next step is that the file '7_step3ML_brew.brew' needs to be put in every condition folder
by hand. This is a copy paste procedure that is not automated so should be done by hand. 
You should also make sure that the infile information in the '7_step3ML_brew.brew' file 
corresponds to the simulation condition.

The next step is to open the '7_TOP_step3ML.R' file and run the code inside it. It will
source the function saved under '7_template_step3ML.R'. 

After this procedure is finished, every simulation condition folder contains

- the '7_step3ML_brew.brew' file you put in there by hand and for which you adjusted the 
  'infile' by hand
- For every iteration and every bootstrap sample within that iteration you find a file 
  containing the posterior membership probabilities. They are named:
  LC3step_<simulation iteration><bootstrap>.dat
- the model.<simulation iteration><bootstrap>.lgs file contains the Latent Gold syntax 
  you produced
- the model.<simulation iteration><bootstrap>.lst file contains the corresponding Latent 
  Gold output.

In the next step, it was investigated which iterations in which simulation conditions
did not converge, and an overview of these iterations is given in 
'extra_iterations_overview.txt'. These iterations are rerun while the number of 
emiterations, nriterations and iterations were all specified to 5000 instead of 500.