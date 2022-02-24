rm(list=ls())
library(phyloseq)
library(vegan)
library(gamlss)

#moved the raw objects to the 16S and ITS data folders of the repo. 
phyloseqD <- readRDS("./data/PhyloseqObjects/16S/merged_simper_baruta.Rdata")

#330 taxa in original
phyloseqD

#Convert into a normal data frame
otu_tab_bac <- data.frame(t(otu_table(phyloseqD)))
metadata_bac <- data.frame(sample_data(phyloseqD))
metadata_bac$Time_herb <- paste(metadata_bac$Time, metadata_bac$Herbicide, sep = "_")

#QC
table(rownames(otu_tab_bac) == rownames(metadata_bac))

#adding sampling info back in and renaming
dat <- cbind(metadata_bac$Time_herb, otu_tab_bac)
names(dat)[names(dat) == 'metadata_bac$Time_herb'] <- 'Time_herb'

#Order by treatment, clean up a bit, and doublecheck that the data look good
dat <- dat[order(dat$Time_herb),]
treatments_full <- dat[,1]
dat <- dat[,-1]
dat[1:5,1:10]

# Standardize (or rarify)
#dat_rare <- data.frame(vegan::drarefy(dat, 500)) #This results in proportion data, which may not model properly with the lm

# Hellinger standardize 
dat_rare <- sqrt(dat / rowSums(dat))

dat_rare_with_treatment <- data.frame(treatments_full, dat_rare)

#Make a few new fields in dat_rare_with_treatment just so we are sure to not mess some indexing up. 
dat_rare_with_treatment$time <- as.factor(gsub("(T[123]).*","\\1", dat_rare_with_treatment$treatments_full))
dat_rare_with_treatment$herbicide <- gsub("T[123]_(.*)","\\1", dat_rare_with_treatment$treatments_full)
dat_rare_with_treatment$herbicide <- relevel(as.factor(dat_rare_with_treatment$herbicide), ref = "Non-Treated")

regressionresults <- data.frame(matrix(ncol = 9, nrow = 1))

#Use a zero-inflated beta regression to deal with our proportion data that lie within the interval [0,1)
# See this nice stackexchange post for more: https://stats.stackexchange.com/questions/309047/zero-inflated-beta-regression-using-gamlss-for-vegetation-cover-data

k <- 1
for(i in names(otu_tab_bac)){
  reg <- gamlss(dat_rare_with_treatment[,names(dat_rare_with_treatment) == i] ~ 
              dat_rare_with_treatment$time + dat_rare_with_treatment$herbicide, family = BEZI)
  regressionresults[k, 1] <- i
  regressionresults[k, 2:9] <- c(coef(reg, what = "mu"),
                                  Rsq(reg)) #generalised R-squared of Nagelkerke (1991) for a GAMLSS mode
  k <- k + 1
}

#NOTE: round up and T1 are the reference conditions. 
names(regressionresults) <- c("taxon", "intercept", "timeT2", "timeT3", "Aatrex", "Clarity","Hand","Powermax", "Rsquared")
max(regressionresults$Rsquared)
summary(regressionresults$Rsquared)
