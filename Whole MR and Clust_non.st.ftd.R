## script same as two sample MR but no steiger filtering, create new dataset to export 

library(TwoSampleMR)
library(ggplot2)
library(tidyverse)
library(readxl)

# set seed 
set.seed(14771998)

setwd("/Users/marma900/Documents/RStudio")

## BMI data
# read in BMI GWAS (and change the p-value that is too small)

BMI_GWAS <- read_excel("Huang2022_BMI_meta extra new.xlsx")

head(BMI_GWAS)

# keep only SNPs p≤5e-8, remove duplicates ##performed 23Apr: remove non-biallelic in BOTH NEA and EA (17apr addition) 
# then format to exposure dataset for TwoSampleMR

BMI_ftd <- filter(BMI_GWAS, P <= 5e-8)

BMI_ftd <- BMI_ftd[nchar(BMI_ftd$NEA) == 1, ]

BMI_ftd <- BMI_ftd[nchar(BMI_ftd$EA) == 1, ]

BMI_ftd <- BMI_ftd[!duplicated(BMI_ftd$rsID), ]

bmi_exp <- format_data(BMI_ftd, type = "exposure", 
                       snp_col = "rsID", 
                       beta_col = "BETA", 
                       se_col = "SE", 
                       
                       samplesize_col = "sample_size",
                       
                       effect_allele_col = "EA", 
                       other_allele_col = "NEA",
                       pval_col = "P", 
                       chr_col = "chr", 
                       pos_col = "pos")
# taken out bc missing in dataset:                   eaf_col = "Freq_Tested_Allele_in_HRS", 






## BC data
# load BC mixed GWAS data, keep meta elements only

BC_GWAS_mix <- vroom::vroom("icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt")
names(BC_GWAS_mix)

BC_meta <- select(BC_GWAS_mix, var_name, Freq.Gwas, Effect.Meta, Baseline.Meta, Beta.meta, sdE.meta, p.meta)

BC_meta <- separate(BC_meta, var_name, into = c("chr", "pos", "NonEffAll", "EffAll"))


## both
# create new column in BMI set: chr_pos, which will contain chrom and pos, to then align to chr_pos
#   column in BC set and only keep SNPs of interest -> make dataset smaller. 
#   also keep rsIDs so can align them and have them in BC set

BMI_exp_mini <- unite(bmi_exp, chr.pos.BMI, chr.exposure, pos.exposure, sep = ".", remove = FALSE)
BMI_exp_mini <-  dplyr::select(BMI_exp_mini, chr.pos.BMI, SNP)

BC_meta <-  unite(BC_meta, chr.pos.BC, chr, pos, sep = ".", remove = FALSE)

BC_meta <- BC_meta %>% inner_join(BMI_exp_mini, by = c("chr.pos.BC" = "chr.pos.BMI"))


# add columns to outcome (BC): N.cases, N.controls, N.total

BC_meta <- mutate(BC_meta, 
                  N.cases = 133384,
                  N.controls = 113789,
                  N.total = 247173)


## prepare BC as outcome data 

bc_outc <- format_data(BC_meta, type = "outcome", 
                       snp_col = "SNP", 
                       beta_col = "Beta.meta", 
                       se_col = "sdE.meta", 
                       eaf_col = "Freq_Tested_Allele_in_HRS",
                       effect_allele_col = "Effect.Meta", 
                       other_allele_col = "Baseline.Meta",
                       
                       ncase_col = "N.cases",
                       ncontrol_col = "N.controls",
                       samplesize_col = "N.total",
                       
                       
                       pval_col = "p.meta", 
                       chr_col = "chr", 
                       pos_col = "pos")

# harmonisation, and ### remove any palindromic SNPs ### (step added 27mar24)

dat_non.st.ftd <- harmonise_data(
    exposure_dat = bmi_exp, 
    outcome_dat = bc_outc)

dat_non.st.ftd <- dat_non.st.ftd[dat_non.st.ftd$mr_keep == TRUE, ]

# save dat
save(dat_non.st.ftd, file = "dat2sMR_non.st.ftd.Rda")


### MR Clust ###

library(mrclust)

## load datasets if starting from here
load("dat2sMR_non.st.ftd.Rda")

# take values from TwoSampleMR "dat" dataframe, name after vectors that will be needed for MRClust

bx <- dat_non.st.ftd$beta.exposure
by <- dat_non.st.ftd$beta.outcome
bxse <- dat_non.st.ftd$se.exposure
byse <- dat_non.st.ftd$se.outcome
snp_names <- dat_non.st.ftd$SNP

ratio_est <- by / bx
ratio_est_se <- byse / abs(bx)

# run MRClust with values needed for it

res_mrc_non.st.ftd <- mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, bx = bx, by = by, bxse = bxse, byse = byse, obs_names = snp_names)

# get clust results data and save
clustdata_non.St.ftd <- res_mrc_non.st.ftd$plots$two_stage$data
save(clustdata_non.St.ftd, file = "clustdata_non.St.ftd.Rda")

# plot all results

plot_BMI_BC_best = res_mrc_non.st.ftd$plots$two_stage +
  ggplot2::xlim(0, max(abs(bx) + 2*bxse)) +
  ggplot2::xlab("Genetic association with BMI") +
  ggplot2::ylab("Genetic association with BC") +
  ggplot2::ggtitle("All SNPs and their association with BMI and BC, non-St ftd")

plot_BMI_BC_best


# only plot clusters with min 4 variants and allocation probability ≥80% 

res80 = mrclust::pr_clust(dta = res_mrc_non.st.ftd$results$best, prob = 0.8, min_obs =  4)

keep80 = which(snp_names %in% res80$observation)
bx80   = bx[keep80]
bxse80 = bxse[keep80]
by80   = by[keep80]
byse80 = byse[keep80]
snp_names80 = snp_names[keep80]

plot.BMI_BC.pr80 = two_stage_plot(res = res80, bx = bx80, by = by80, bxse = bxse80,
                                  byse = byse80, obs_names = snp_names80) + 
  ggplot2::xlim(0, max(abs(bx80) + 2*bxse80)) + 
  ggplot2::xlab("Genetic association with BMI") + 
  ggplot2::ylab("Genetic association with BC") + 
  ggplot2::ggtitle("SNPs assoc. with BMI and BC, ≥80% cluster assoc. prob. and min. 4/cluster, non-St ftd");

plot.BMI_BC.pr80


## get causal estimates for each cluster with TSMR ##

library(TwoSampleMR)

# load data if starting here
load("clustdata_non.St.ftd.Rda")

# or if coming from previious steps: get list of SNPs with clusters
clustdata_non.St.ftd <- res_mrc_non.st.ftd$plots$two_stage$data

# filter data for individual clusters, create dataset for every cluster
clust1_nsf<- clustdata_non.St.ftd[clustdata_non.St.ftd$cluster == 1, ]
clust2_nsf<- clustdata_non.St.ftd[clustdata_non.St.ftd$cluster == 2, ]
clust3_nsf<- clustdata_non.St.ftd[clustdata_non.St.ftd$cluster == 3, ]
clust4_nsf<- clustdata_non.St.ftd[clustdata_non.St.ftd$cluster == 4, ]
clustNull_nsf<- clustdata_non.St.ftd[clustdata_non.St.ftd$cluster == 0, ]


# join clust and data prepared in TSMR for individual clusters, by SNP rsIDs

clust1_nsf_forTSMR <- dat_non.st.ftd %>% inner_join(clust1_nsf, by = c("SNP" = "observation"))
clust2_nsf_forTSMR <- dat_non.st.ftd %>% inner_join(clust2_nsf, by = c("SNP" = "observation"))
clust3_nsf_forTSMR <- dat_non.st.ftd %>% inner_join(clust3_nsf, by = c("SNP" = "observation"))
clust4_nsf_forTSMR <- dat_non.st.ftd %>% inner_join(clust4_nsf, by = c("SNP" = "observation"))
clustNull_nsf_forTSMR <- dat_non.st.ftd %>% inner_join(clustNull_nsf, by = c("SNP" = "observation"))

# run TSMR for each cluster

clust1_TSMR_res_nsf <- mr(clust1_nsf_forTSMR)
clust2_TSMR_res_nsf <- mr(clust2_nsf_forTSMR)
clust3_TSMR_res_nsf <- mr(clust3_nsf_forTSMR)
clust4_TSMR_res_nsf <- mr(clust4_nsf_forTSMR)
clustNull_TSMR_res_nsf <- mr(clustNull_nsf_forTSMR)


### Own MRC plotting attempt based on Pascal's 
## using IVW from TSMR

install.packages("MetBrewer")
library(MetBrewer)

clustdata_non.St.ftd <- res_mrc_non.st.ftd$plots$two_stage$data

# create set that has Inverse Variance Weighted beta from TSMR, to use for slope
slopes <- data.frame(cluster_class = c(1, 2, 3, 4), 
                     slope = c(-0.9983151, 2.310281, 1.2493482, -2.685874))

# plot 

clustdata_non.St.ftd%>%
  filter(cluster_class != "Junk", probability>0.8) %>%
  ggplot(aes(x=bx, y = by, color = cluster_class)) +
  
  geom_point(aes(size = probability), alpha = 0.4) +
  geom_errorbarh(aes(xmin = bx - 1.96 * bxse,
                     xmax = bx + 1.96 * bxse,
                     color = cluster_class),
                 linetype = "solid",alpha = 0.2) +
  
  geom_errorbar(aes(ymin = by - 1.96 * byse,
                    ymax = by + 1.96 * byse,
                    color = cluster_class),
                linetype = "solid", alpha = 0.2) +
  
  labs(title = "Clustered Causal Estimates of BC ~ BMI, P > 0.8", 
       color = "Cluster", size = "Cluster \ninclusion \nprobability") +
  xlab("Genetic effect on BMI") +
  ylab("Genetic effect on BC") +
  scale_color_met_d("Archambault") +
  scale_fill_met_d("Archambault") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold")) +
  geom_abline(data = slopes, aes(slope = slope, intercept = 0, color = as.factor(cluster_class)), linetype = "dotted")



### additional manipulations ###


## export lists of SNP clusters for PRS calculations as txt files

# get list of SNPs with clusters
clustdata_non.St.ftd <- res_mrc_non.st.ftd$plots$two_stage$data

# filter data for individual clusters, create dataset for every cluster
clust1_nsf<- clustdata_non.St.ftd[clustdata_non.St.ftd$cluster == 1, ]
clust2_nsf<- clustdata_non.St.ftd[clustdata_non.St.ftd$cluster == 2, ]
clust3_nsf<- clustdata_non.St.ftd[clustdata_non.St.ftd$cluster == 3, ]
clust4_nsf<- clustdata_non.St.ftd[clustdata_non.St.ftd$cluster == 4, ]
clustNull_nsf<- clustdata_non.St.ftd[clustdata_non.St.ftd$cluster == 0, ]

# get MVP only BMI GWAS

BMI_MVP_GWAS <- read_excel("Huang2022_BMI_MVP only.xlsx")

# join clusts and data prepared in TSMR for individual clusters, by SNP rsIDs
 
clust1_nsf_forPRS <- BMI_MVP_GWAS %>% inner_join(clust1_nsf, by = c("rsID" = "observation"))
clust2_nsf_forPRS <- BMI_MVP_GWAS %>% inner_join(clust2_nsf, by = c("rsID" = "observation"))
clust3_nsf_forPRS <- BMI_MVP_GWAS %>% inner_join(clust3_nsf, by = c("rsID" = "observation"))
clust4_nsf_forPRS <- BMI_MVP_GWAS %>% inner_join(clust4_nsf, by = c("rsID" = "observation"))
clustNull_nsf_forPRS <- BMI_MVP_GWAS %>% inner_join(clustNull_nsf, by = c("rsID" = "observation"))

# select only columns needed for BMI PRS checking
clust1_nsf_forPRS <- select(clust1_nsf_forPRS, rsID, EA, BETA, cluster)
clust2_nsf_forPRS <- select(clust2_nsf_forPRS, rsID, EA, BETA, cluster)
clust3_nsf_forPRS <- select(clust3_nsf_forPRS, rsID, EA, BETA, cluster)
clust4_nsf_forPRS <- select(clust4_nsf_forPRS, rsID, EA, BETA, cluster)
clustNull_nsf_forPRS <- select(clustNull_nsf_forPRS, rsID, EA, BETA, cluster)


# export as txt files
write.table(clust1_nsf_forPRS, "Clusters_PRS/clust1nsf_forPRS.txt", quote = FALSE, row.names = FALSE)
write.table(clust2_nsf_forPRS, "Clusters_PRS/clust2nsf_forPRS.txt", quote = FALSE, row.names = FALSE)
write.table(clust3_nsf_forPRS, "Clusters_PRS/clust3nsf_forPRS.txt", quote = FALSE, row.names = FALSE)
write.table(clust4_nsf_forPRS, "Clusters_PRS/clust4nsf_forPRS.txt", quote = FALSE, row.names = FALSE)
write.table(clustNull_nsf_forPRS, "Clusters_PRS/clustNullnsf_forPRS.txt", quote = FALSE, row.names = FALSE)





