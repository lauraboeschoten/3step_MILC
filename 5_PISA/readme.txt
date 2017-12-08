Every folder contains an R script that can be run to obtain the results corresponding
to the latent gold / three-step model specified in the folder. The ML and BCH folder can
only be run after the 'excluded' model is run, since some of the output of the 'excluded'
model is needed. 

Note that when latent class analysis is applied in this dataset, the problem of label 
switching can arise. If you run it for yourself, the order of the latent classes can be different
compared to the models that were used here. The R could to create imputations based on the 
posteriors should then be adjusted correspondingly. This has to be checked by hand. 

All folders contain the relevant latent gold output that is needed to obtain the results as were
discussed in the paper.

For this example, we use two SPSS files that are linked by person ID. 
The first one, "NL_kids.sav" contains a subset of the original "CY6_MS_CMB_STU_QQQ.sav"
file that can be downloaded from 
http://vs-web-fs-1.oecd.org/pisa/PUF_SPSS_COMBINED_CMB_STU_QQQ.zip with only the
Dutch students.

The second one, "PISA_NL_book43.sav" contains a subset of the original 
"CY6_MS_CMB_STU_COG.sav.sav" file that can be downloaded from
http://vs-web-fs-1.oecd.org/pisa/PUF_SPSS_COMBINED_CMB_STU_COG.zip
and only contains the questions from booklet 43 and responses from the Dutch students. 



