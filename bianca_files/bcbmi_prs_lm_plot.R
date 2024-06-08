### script for BMI PRS analysing and plotting
### PART I: lm  
### PART II: plot lm BMI ~ PRS for all clusters
### PART III: forest plot of lm BMI ~ PRS estimates, all clusters


library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)


# load phenotype data and keep females only
load("/proj/sens2017538/nobackup/marina/data/bc_phenos.Rda")
bc_phenos <- filter(bc_phenos, sex == 0)

# load PRS for each cluster
clust1PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmi_prs_bcclust/clust1nsf_forPRS.profile", header = TRUE)
clust2PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmi_prs_bcclust/clust2nsf_forPRS.profile", header = TRUE)
clust3PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmi_prs_bcclust/clust3nsf_forPRS.profile", header = TRUE)
clust4PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmi_prs_bcclust/clust4nsf_forPRS.profile", header = TRUE)
clustNullPRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmi_prs_bcclust/clustNullnsf_forPRS.profile", header = TRUE)

# join sets
clust1prs_BCpheno <- inner_join(bc_phenos, clust1PRS, by = c("IID" = "IID"))
clust2prs_BCpheno <- inner_join(bc_phenos, clust2PRS, by = c("IID" = "IID"))
clust3prs_BCpheno <- inner_join(bc_phenos, clust3PRS, by = c("IID" = "IID"))
clust4prs_BCpheno <- inner_join(bc_phenos, clust4PRS, by = c("IID" = "IID"))
clustNullprs_BCpheno <- inner_join(bc_phenos, clustNullPRS, by = c("IID" = "IID"))

### PART I: run lm ###

# run lm and save summaries with CI
clust1_lm_bc <- glm(BC ~ SCORESUM, data = clust1prs_BCpheno, family = binomial)
clust1_lm_bc <- tidy(clust1_lm_bc, conf.int = TRUE)

clust2_lm_bc <- glm(BC ~ SCORESUM, data = clust2prs_BCpheno, family = binomial)
clust2_lm_bc <- tidy(clust2_lm_bc, conf.int = TRUE)

clust3_lm_bc <- glm(BC ~ SCORESUM, data = clust3prs_BCpheno, family = binomial)
clust3_lm_bc <- tidy(clust3_lm_bc, conf.int = TRUE)

clust4_lm_bc <- glm(BC ~ SCORESUM, data = clust4prs_BCpheno, family = binomial)
clust4_lm_bc <- tidy(clust4_lm_bc, conf.int = TRUE)

clustNull_lm_bc <- glm(BC ~ SCORESUM, data = clustNullprs_BCpheno, family = binomial)
clustNull_lm_bc <- tidy(clustNull_lm_bc, conf.int = TRUE)




### PART II: plot lm all clusters ###

# Add identifier column to dataframes containing prs and pheno info
clust1prs_BCpheno$clust <- "clust1"
clust2prs_BCpheno$clust <- "clust2"
clust3prs_BCpheno$clust <- "clust3"
clust4prs_BCpheno$clust <- "clust4"
clustNullprs_BCpheno$clust <- "clustNull"

# Combine dataframes
allBCphenos <- rbind(clust1prs_BCpheno, clust2prs_BCpheno, clust3prs_BCpheno, clust4prs_BCpheno, clustNullprs_BCpheno)

# Plot
ggplot(allBCphenos, aes(x = SCORESUM, y = BC, color = clust)) +
  geom_smooth(method = "lm", se = FALSE, fullrange=TRUE) +
  labs(title = "BC ~ PRS, all clusters") +
  xlab("PRS") +
  ylab("BC risk")  +
  theme_classic(base_size = 14) +
  scale_color_met_d("Archambault") 



### PART III: forest plot (from here onwards Null cluster = cluster5?)

# join lm tibbles
allclusts_lm_bc <- bind_rows(clust1_lm_bc, clust2_lm_bc, clust3_lm_bc, clust4_lm_bc, clustNull_lm_bc,
                          .id = "cluster")
allclusts_lm_bc <- filter(allclusts_lm_bc, term=='SCORESUM')

## forest plot alone
bc_fplot <- 
  allclusts_lm_bc %>% 
  ggplot(aes(y = fct_rev(cluster), color = cluster))+
  theme_classic(base_size = 14) +
  
  geom_point(aes(x = estimate), shape = 15, size = 4) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) + 
  
  geom_vline(xintercept = 0, linetype = "dotted") +
  
  labs(x = NULL, y = NULL) +
  
#  labs (title = "BC-PRS association", x = "Beta", y = "Cluster") +
  theme(plot.title=element_text(face="bold")) +
  
  theme(legend.position = "none") + 
  
  scale_color_met_d("Archambault") 


bc_fplot


## attempt using khstats approach, patchwork of ggplots

#add column with SNP nrs
allclusts_lm_bc$SNP_count <- c(391, 38, 146, 18, 1517)

# add row to adapt spaces 
# allclusts_lm_bc <- rbind(allclusts_lm_bc, list(' ', '', '', '', '', '', '', '', 'SNP count'))



library(tidyverse)
library(MetBrewer)

p_main <- 
  allclusts_lm_bc %>% 
  ggplot(aes(y = fct_rev(cluster), color = cluster))+
  theme_classic() +
  
  geom_point(aes(x = estimate), shape = 15, size = 3) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) + 
  
  geom_vline(xintercept = 0, linetype = "dashed") +

    labs (title = "BC-PRS association", x = "Beta", y = "Cluster") +
  theme(plot.title=element_text(face="bold")) +
  
  theme(legend.position = "none") + 
  
  scale_color_met_d("Archambault") 


p_main

# make the right side of the plot: SNP counts

# make extra dataset to not mess up main plot
allclusts_lm_bc_extra <- rbind(allclusts_lm, list(0, '', '', '', '', '', '', '', 'SNP count'))


p_right <-
  allclusts_lm_bc_extra  |>
  ggplot() +
  geom_text(
    aes(x = 0, y = fct_rev(cluster), label = SNP_count),
    hjust = 0,
    fontface = ifelse(allclusts_lm_extra$SNP_count == "SNP count", "bold", "plain")
  ) +
  theme_void()

p_right

# put together
library(patchwork)

layout <- c(area(t = 4, l = 0, b = 30, r = 9),
            area(t = 0, l = 9, b = 30, r = 10))


p_main + p_right + plot_layout(design = layout)













########### might not need #############
# for labels of cluster
# change names
allclusts_lm_bc$cluster <- factor(allclusts_lm_bc$cluster,
                               levels = 1:5,
                               labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Null Cluster"))

# reverse order of clusters so will be 1 to Null once flip coordinates
allclusts_lm_bc$cluster <- factor(allclusts_lm_bc$cluster, levels=rev(allclusts_lm_bc$cluster))
########### might not need ends #############



