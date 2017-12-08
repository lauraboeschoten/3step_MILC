The function 'create_folders.R' created folders for every condition in the simulation. 
The next step is that the file '3_brew_LCmodel.brew' needs to be put in every condition folder
by hand. This is a copy paste procedure that is not automated so should be done by hand. 
You should also make sure that the infile information in the '3_brew_LCmodel.brew' file 
corresponds to the simulation condition.

The next step is to open the '3_TOP_LC_models.R' file and run the code inside it. It will
source the functions saved under '3_template_LCmodel.R'. 

After this procedure is finished, every simulation condition folder contains

- the 3_brew_LCmodel.brew file you put in there by hand and for which you adjusted the 'infile' by hand
- For every iteration and every bootstrap sample within that iteration you find a file 
  containing the posterior membership probabilities. They are named:
  LC_<simulation iteration><bootstrap>.dat
- the model.<simulation iteration><bootstrap>.lgs file contains the Latent Gold syntax 
  you produced
- the model.<simulation iteration><bootstrap>.lst file contains the corresponding Latent 
  Gold output