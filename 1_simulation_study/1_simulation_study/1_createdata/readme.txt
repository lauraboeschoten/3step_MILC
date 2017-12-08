The function 'create_folders.R' created folders for every condition in the simulation. 
The next step is that the file '1_create_data.brew' needs to be put in every condition folder
by hand. This is a copy paste procedure that is not automated so should be done by hand. 

The next step is to open the '1_TOP_createdata.R' file and run the code inside it. It will
source the functions saved under '1_fun_createdata.R' and '1_createdata_template.R'. 

When this procedure is finished, every simulation condition folder contains: 
- the 1_create_data.brew file you put in there by hand
- a small create_ss<sample size>.txt file which is used as imput for Latent Gold
- the dat_<sample size>_<classification probability>_<P(Z=2)>_<P(Q=2)>.txt file containing
  all datasets that are simulated for this condition
- the LG_<sample size>_<classification probability>_<P(Z=2)>_<P(Q=2)>.lgs file contains the 
  Latent Gold syntax you produced
- the LG_<sample size>_<classification probability>_<P(Z=2)>_<P(Q=2)>.lst file contains the
  corresponding Latent Gold output