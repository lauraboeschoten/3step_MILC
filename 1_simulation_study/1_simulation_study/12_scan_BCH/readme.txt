The code in the file '12_TOP_scan.R' runs through all matrices produced in step 11, 
and investigates whether any inadmissable solutions are produced. An inadmissable 
solution occurs if for a specific combination of scores on the covariates and latent 
variable, there is a zero probability of being in any of the latent classes.
With such an outcome, it is not possible to estimate posterior membership probailities
for this combination of scores and therefore the results are inadmissable.