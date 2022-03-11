rm(list=ls())
library(phyloseq)
library(vegan)
library(gamlss)

#moved the raw objects to the 16S and ITS data folders of the repo. 
phyloseqD <- readRDS("./data/PhyloseqObjects/ITS/merged_simper_baruta_ITS.Rdata")

#66 taxa in original ITS object
phyloseqD

#Convert into a normal data frame
otu_tab_ITS <- data.frame(t(otu_table(phyloseqD)))
metadata_ITS <- data.frame(sample_data(phyloseqD))
tax_ITS <- data.frame(tax_table(phyloseqD))
metadata_ITS$Time_herb <- paste(metadata_ITS$Time, metadata_ITS$Herbicide, sep = "_")

#QC
table(rownames(otu_tab_ITS) == rownames(metadata_ITS))

#adding sampling info back in and renaming
dat <- cbind(metadata_ITS$Time_herb, otu_tab_ITS)
names(dat)[names(dat) == 'metadata_ITS$Time_herb'] <- 'Time_herb'

#Order by treatment, clean up a bit, and doublecheck that the data look good
dat <- dat[order(dat$Time_herb),]
treatments_full <- dat[,1]
dat <- dat[,-1]
dat[1:5,1:10]

# Standardize (or rarify)
#dat_rare <- data.frame(vegan::drarefy(dat, 500)) #This results in proportion data, which may not model properly with the lm

# Hellinger standardize 
dat_rare <- sqrt(dat / rowSums(dat))

#Or calculate the probability a taxon is in a sample of 500...a multinomial 
#dat_rare <- data.frame(vegan::drarefy(dat, 500))


dat_rare_with_treatment <- data.frame(treatments_full, dat_rare)

#Make a few new fields in dat_rare_with_treatment just so we are sure to not mess some indexing up. 
dat_rare_with_treatment$time <- as.factor(gsub("(T[123]).*","\\1", dat_rare_with_treatment$treatments_full))
dat_rare_with_treatment$herbicide <- gsub("T[123]_(.*)","\\1", dat_rare_with_treatment$treatments_full)
dat_rare_with_treatment$herbicide <- relevel(as.factor(dat_rare_with_treatment$herbicide), ref = "Non-Treated")

regressionresults <- data.frame(matrix(ncol = 9, nrow = 1))

#Use a zero-inflated beta regression to deal with our proportion data that lie within the interval [0,1)
# See this nice stackexchange post for more: https://stats.stackexchange.com/questions/309047/zero-inflated-beta-regression-using-gamlss-for-vegetation-cover-data

k <- 1
for(i in names(otu_tab_ITS)){
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

plot(table(tax_ITS$Order), las = 2)


taxa_and_reg_results<-cbind(regressionresults, tax_ITS)
write.csv(taxa_and_reg_results, "/Users/gordoncuster/Desktop/Git_Projects/Herbicide_Microbes_PT1/data/PhyloseqObjects/ITS/RegressionresultsITS.csv")
