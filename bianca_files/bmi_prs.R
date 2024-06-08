#script for running prs
#load libraries

library(tidyverse)

### make lists of txt files (of tibbles) (read in, assign to list)

c1=vroom::vroom("/proj/sens2017538/nobackup/marina/data/Clusters_PRS/clust1nsf_forPRS.txt")
c2=vroom::vroom("/proj/sens2017538/nobackup/marina/data/Clusters_PRS/clust2nsf_forPRS.txt")
c3=vroom::vroom("/proj/sens2017538/nobackup/marina/data/Clusters_PRS/clust3nsf_forPRS.txt")
c4=vroom::vroom("/proj/sens2017538/nobackup/marina/data/Clusters_PRS/clust4nsf_forPRS.txt")
cNull=vroom::vroom("/proj/sens2017538/nobackup/marina/data/Clusters_PRS/clustNullnsf_forPRS.txt")

cluster_list = list(clust1nsf_forPRS = c1,
                    clust2nsf_forPRS = c2, 
                    clust3nsf_forPRS = c3,
                    clust4nsf_forPRS = c4, 
                    clustNullnsf_forPRS = cNull)

clusters = names(cluster_list)

#create function for cmd

plink_rsc = "plink --bfile /proj/sens2017538/nobackup/pascal/mukg_merged " # can use same files as Pascal
score = "--score /proj/sens2017538/nobackup/marina/data/Clusters_PRS/" #folder where I have SNP files for cluster
 
cols = "1 2 3 header sum "
out = "--out /proj/sens2017538/nobackup/marina/results/bmi_prs_bcclust_NEW/"

#function
prs_func = function(dat){
  cmd = paste0(plink_rsc,
               score,paste0(dat, ".txt "),cols,
               out,paste0(dat))
  
  system(cmd)
}

#use function to loop over all PRSs

lapply(clusters, function(x) prs_func(x))


## for when running in server
q()
