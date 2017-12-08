The function 'create_folders.R' created folders for every condition in the simulation. 

The code in '11_TOP_BCH2.R' applies the BCH correction method (the improved version of the
method, using the qpsolve() function). A corrected matrix is created for every bootstrap
sample in every iteration, and this for each simulation condition. The files are named
'BCH_amatrix_<iterationnr1-1000><bootstrapnr1-5>.txt'.