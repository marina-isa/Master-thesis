# Genetic predisposition for obesity and relationship with breast cancer - Codes

`Whole MR and Clust_non.st.ftd.R` is the script used for data preparation in 2SMR, and clustering with MR Clust. 

`bianca_files` contains scripts using sensitive data from the UK Biobank, namely: 
-  `bmi_prs.sh`, the bash script to run PLINK
-  `bmi_prs.R`, the R script used in the bash script to calculate polygenic scores (PGS)
-  `bmiprs_lm_plot.R`, the R script for rank-based inverse normal transformation of UK Biobank individuals' BMI, running linear regression analysis on transformed BMI with PGS, and forest plot of beta estimates
-  `bcbmi_prs_plot.R`, the R script for logistic regression of breast cancer in UK Biobank women and PGS, and forest plot of beta estimates
