### script for BMI PRS analysing and plotting
### PART I: lm  
### PART II: plot lm BMI ~ PRS for all clusters
### PART III: forest plot of lm BMI ~ PRS estimates, all clusters

library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)

install.packages("MetBrewer")
library(MetBrewer)
library(tidyverse)


# load phenotype data and keep females only
load("/proj/sens2017538/nobackup/marina/data/bmi_phenos.Rda")
bmi_phenos <- filter(bmi_phenos, sex == 0)

# create column with rank-based inverse normal transformation
bmi_phenos <- rename(bmi_phenos, BMI.raw = BMI)

rntransf <- function(x) {
  tmp.1 <- (order(order(x[!is.na(x)]))-0.5)/length(x[!is.na(x)])
  tmp.2 <- rep(NA,length(x))
  tmp.2[!is.na(x)] <- tmp.1
  return(qnorm(tmp.2,0,1))
}

BMI <- rntransf(bmi_phenos$BMI.raw)

bmi_phenos$BMI <- BMI


# load PRS for each cluster
clust1PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmi_prs_bcclust/clust1nsf_forPRS.profile", header = TRUE)
clust2PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmi_prs_bcclust/clust2nsf_forPRS.profile", header = TRUE)
clust3PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmi_prs_bcclust/clust3nsf_forPRS.profile", header = TRUE)
clust4PRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmi_prs_bcclust/clust4nsf_forPRS.profile", header = TRUE)
clustNullPRS <- read.table("/proj/sens2017538/nobackup/marina/results/bmi_prs_bcclust/clustNullnsf_forPRS.profile", header = TRUE)

# join sets
clust1prs_pheno <- inner_join(bmi_phenos, clust1PRS, by = c("IID" = "IID"))
clust2prs_pheno <- inner_join(bmi_phenos, clust2PRS, by = c("IID" = "IID"))
clust3prs_pheno <- inner_join(bmi_phenos, clust3PRS, by = c("IID" = "IID"))
clust4prs_pheno <- inner_join(bmi_phenos, clust4PRS, by = c("IID" = "IID"))
clustNullprs_pheno <- inner_join(bmi_phenos, clustNullPRS, by = c("IID" = "IID"))


### PART I: run lm ###

# run lm and save sumamries with CI
clust1_lm <- lm(BMI ~ SCORESUM, data = clust1prs_pheno)
clust1_lm <- tidy(clust1_lm, conf.int = TRUE)

clust2_lm <- lm(BMI ~ SCORESUM, data = clust2prs_pheno)
clust2_lm <- tidy(clust2_lm, conf.int = TRUE)

clust3_lm <- lm(BMI ~ SCORESUM, data = clust3prs_pheno)
clust3_lm <- tidy(clust3_lm, conf.int = TRUE)

clust4_lm <- lm(BMI ~ SCORESUM, data = clust4prs_pheno)
clust4_lm <- tidy(clust4_lm, conf.int = TRUE)

clustNull_lm <- lm(BMI ~ SCORESUM, data = clustNullprs_pheno)
clustNull_lm <- tidy(clustNull_lm, conf.int = TRUE)


### PART II: plot lm all clusters ###

# Add identifier column to dataframes containing prs and pheno info
clust1prs_pheno$clust <- "clust1"
clust2prs_pheno$clust <- "clust2"
clust3prs_pheno$clust <- "clust3"
clust4prs_pheno$clust <- "clust4"
clustNullprs_pheno$clust <- "clustNull"

# Combine dataframes
allphenos <- rbind(clust1prs_pheno, clust2prs_pheno, clust3prs_pheno, clust4prs_pheno, clustNullprs_pheno)

# Plot
ggplot(allphenos, aes(x = SCORESUM, y = BMI, color = clust)) +
  geom_smooth(method = "lm", se = FALSE, fullrange=TRUE) +
  labs(title = "BMI ~ PRS, all clusters") +
  xlab("PRS") +
  ylab(bquote(BMI[INT])) +
#  expand_limits(x = c(0, 3), y = c(20, 35)) + 
  theme_classic(base_size = 14) +
  scale_fill_discrete(name = "Cluster", labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Null cluster")) +
  scale_color_met_d("Archambault") 





### PART III: forest plot (from here onwards Null cluster = cluster5?)

# join lm tibbles
allclusts_lm <- bind_rows(clust1_lm, clust2_lm, clust3_lm, clust4_lm, clustNull_lm,
                         .id = "cluster")
allclusts_lm <- filter(allclusts_lm, term=='SCORESUM')


## forest plot alone
bmi_fplot <- 
  allclusts_lm %>% 
  ggplot(aes(y = fct_rev(cluster), color = cluster))+
  theme_classic(base_size = 16) +
  
  geom_point(aes(x = estimate), shape = 15, size = 4) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) + 
  
  geom_vline(xintercept = 1, linetype = "dotted") +
  
  #geom_vline(xintercept = 5.136203, linetype = "dashed") +
  #annotate("text", x = 5.1, y = 3, label = "BMI SD", angle = 90) +
  
  
  labs(x = NULL, y = NULL) +
  
  #  labs (title = "BMI-PRS association", x = "Beta", y = "Cluster") +
  theme(plot.title=element_text(face="bold")) +
  
  theme(legend.position = "none") + 
  
  scale_color_met_d("Archambault") 


bmi_fplot




## attempt using khstats approach, patchwork of ggplots (not so recommended)


# for labels of cluster
# change names
allclusts_lm$cluster <- factor(allclusts_lm$cluster,
                               levels = 1:5,
                               labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Null Cluster"))

# reverse order of clusters so will be 1 to Null once flip coordinates
allclusts_lm$cluster <- factor(allclusts_lm$cluster, levels=rev(allclusts_lm$cluster))

# for SNP counts
#add column with SNP nrs
allclusts_lm$SNP_count <- c(391, 38, 146, 18, 1517)




p_main <- 
  allclusts_lm %>% 
  ggplot(aes(y = fct_rev(cluster), color = cluster))+
  theme_classic() +
  
  geom_point(aes(x = estimate), shape = 15, size = 3) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) + 
  
  geom_vline(xintercept = 5.136203, linetype = "dashed") +
  annotate("text", x = 5.1, y = 3, label = "BMI SD", angle = 90) +
  labs (title = "BMI-PRS association", x = "Beta", y = "Cluster") +
  theme(plot.title=element_text(face="bold")) +
  
  theme(legend.position = "none") + 
  scale_color_met_d("Archambault") 

  
p_main

# make the right side of the plot: SNP counts

#add column with SNP nrs
allclusts_lm$SNP_count <- c(391, 38, 146, 18, 1517)
#make extra dataset
allclusts_lm_extra <- rbind(allclusts_lm, list(0, '', '', '', '', '', '', '', 'SNP count'))


p_right <-
  allclusts_lm_extra  |>
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





## until here 04Apr24 18:44























# incomprehensible mess for another day (not Thu04Apr24)

# run linear model
clust1_lm <- lm(BMI ~ SCORESUM, data = clust1prs_pheno)
summary(clust1_lm, conf.int = TRUE, conf.level = 0.95)
confint(clust1_lm)

clust2_lm <- lm(BMI ~ SCORESUM, data = clust2prs_pheno)
summary(clust2_lm, conf.int = TRUE, conf.level = 0.95)
confint(clust2_lm)

clust3_lm <- lm(BMI ~ SCORESUM, data = clust3prs_pheno)
summary(clust3_lm, conf.int = TRUE, conf.level = 0.95)
confint(clust3_lm)

clust3_lm <- lm(BMI ~ SCORESUM, data = clust3prs_pheno)
summary(clust3_lm, conf.int = TRUE, conf.level = 0.95)
confint(clust3_lm)

clust4_lm <- lm(BMI ~ SCORESUM, data = clust4prs_pheno)
summary(clust4_lm, conf.int = TRUE, conf.level = 0.95)
confint(clust4_lm)

clustNull_lm <- lm(BMI ~ SCORESUM, data = clustNullprs_pheno)
summary(clustNull_lm, conf.int = TRUE, conf.level = 0.95)
confint(clustNull_lm)


## multiple regression 
# Add an identifier column
clust1prs_pheno$group <- "PRS1"
clust2prs_pheno$group <- "PRS2"
clust3prs_pheno$group <- "PRS3"
clust4prs_pheno$group <- "PRS4"
clustNullprs_pheno$group <- "PRSNull"

# Combine the dataframes
df_mult_lm <- rbind(clust1prs_pheno, clust2prs_pheno, clust3prs_pheno, clust4prs_pheno, clustNullprs_pheno)

# make group column contain PRS per group value?
df_mult_lm <- df_mult_lm %>%
  spread(key = group, value = SCORESUM)

# run lm on all
all_clusts_lm <- lm(BMI ~ PRS1 + PRS2 + PRS3 + PRS4 + PRSNull, data = df_mult_lm )



# plot lm 
### also to be continued ###
## first part until next disclaimer tested on second run and works ##

# Add an identifier column
clust1prs_pheno$clust <- "clust1"
clust2prs_pheno$clust <- "clust2"
clust3prs_pheno$clust <- "clust3"
clust4prs_pheno$clust <- "clust4"
clustNullprs_pheno$clust <- "clustNull"

# Combine the dataframes
df <- rbind(clust1prs_pheno, clust2prs_pheno, clust3prs_pheno, clust4prs_pheno, clustNullprs_pheno)

# Plot
ggplot(df, aes(x = SCORESUM, y = BMI, color = clust)) +
  geom_smooth(method = "lm", se = FALSE, fullrange=TRUE) +
  labs(title = "PRS ~ BMI, all clusters") +
  xlab("PRS") +
  ylab("BMI")


### part that was tested ends here ###

ggplot(clust1prs_pheno, aes(x=SCORESUM, y = BMI)) +
  geom_smooth(method=lm, se=FALSE) +
  #geom_point() +
  labs(title = "BMI ~ PRS, cluster 1") +
  xlab("PRS") +
  ylab("BMI")

ggplot(clust2prs_pheno, aes(x = SCORESUM, y = BMI)) +
  geom_smooth(method = lm, se = FALSE, fullrange = FALSE) +
  #geom_point() +
  labs(title = "BMI ~ PRS, cluster 2") +
  xlab("PRS") +
  ylab("BMI")



ggplot(aes(x = BMI, y = SCORESUM)) + 
  geom_smooth(data = clust1prs_pheno, method = lm, se = FALSE, color = "#81A0E4") +
  geom_smooth(data = clust2prs_pheno, method = lm, se = FALSE, color = "#7B4872") + 
  geom_smooth(data = clust3prs_pheno, method = lm, se = FALSE, color = "#ED958A") + 
  geom_smooth(data = clust4prs_pheno, method = lm, se = FALSE, color = "#E98324") + 
  geom_smooth(data = clustNullprs_pheno, method = lm, se = FALSE, color = "#FFD103") +
  labs(title = "PRS ~ BMI, all clusters") +
  xlab("BMI") +
  ylab("PRS")


c("#81A0E4", "#7B4872", "#ED958A", "#E98324", "#FFD103")



## adapt to make plot with all BMI-PRS lm? ##

clust1prs_pheno%>%
  ggplot(aes(x=clust1prs_pheno$BMI, y = clust1prs_pheno$SCORESUM, color = cluster_class)) +
  
  geom_point(aes(size = probability), alpha = 0.4) +
  geom_errorbarh(aes(xmin = bx - 1.96 * bxse,
                     xmax = bx + 1.96 * bxse,
                     color = cluster_class),
                 linetype = "solid",alpha = 0.2) +
  
  geom_errorbar(aes(ymin = by - 1.96 * byse,
                    ymax = by + 1.96 * byse,
                    color = cluster_class),
                linetype = "solid", alpha = 0.2) +
  
  labs(title = "Clustered Causal Estimates of BMI ~ BC, P > 0.8") +
  xlab("Genetic effect on BMI") +
  ylab("Genetic effect on BC") +
  scale_color_met_d("Archambault") +
  scale_fill_met_d("Archambault") +
  geom_abline(data = slopes, aes(slope = slope, intercept = 0, color = as.factor(cluster_class)), linetype = "dotted") 


# Plot attempt modeled on Pascal's 
# maybe not necessary: 
SNP_count <- c(391, 38, 146, 18, 1517)


# plot, adapted from Pascal
p_main <- 
  allclusts_lm %>% 
  ggplot(aes(x = cluster, y = estimate)) +
  geom_pointrange(aes(ymin = conf.low,
                      ymax = conf.high),
                  position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 5.136203, lty=2) +
  labs(title = "BMI-PRS association",
       x = "Cluster",
       y = "Estimate") +
  coord_flip() +
  theme_void()

p_main + p_right



