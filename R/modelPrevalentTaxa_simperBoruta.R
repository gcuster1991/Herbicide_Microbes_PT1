rm(list=ls())
library(phyloseq)
library(vegan)
library(gamlss)

#moved the raw objects to the 16S and ITS data folders of the repo. 
phyloseqD <- readRDS("./data/PhyloseqObjects/ITS/merged_simper_baruta_ITS.Rdata")

#330 taxa in original 16S
#66 taxa in original 16S
phyloseqD

#Convert into a normal data frame
otu_tab_bac <- data.frame(t(otu_table(phyloseqD)))
metadata_bac <- data.frame(sample_data(phyloseqD))
metadata_bac$Time_herb <- paste(metadata_bac$Time, metadata_bac$Herbicide, sep = "_")


#QC
table(rownames(otu_tab_bac) == rownames(metadata_bac))

#adding sampling info back in and renaming
dat <- cbind(metadata_bac$Time_herb, metadata_bac$Total_Weed_Veg, otu_tab_bac)
names(dat)[names(dat) == 'metadata_bac$Time_herb'] <- 'Time_herb'
names(dat)[names(dat) == 'metadata_bac$Total_Weed_Veg'] <- 'Total_Weed_Veg'
#Order by treatment, clean up a bit, and doublecheck that the data look good
dat <- dat[order(dat$Time_herb),]
treatments_full <- dat[,c(1:2)]
dat <- dat[,-c(1,2)]

dat[1:5,1:10]

# Standardize (or rarify)
#dat_rare <- data.frame(vegan::drarefy(dat, 500)) #This results in proportion data, which may not model properly with the lm

# Hellinger standardize 
dat_rare <- sqrt(dat / rowSums(dat))

#Or calculate the probability a taxon is in a sample of 500...a multinomial 
#dat_rare <- data.frame(vegan::drarefy(dat, 500))

#remerge 
dat_rare_with_treatment <- data.frame(treatments_full, dat_rare)

#Make a few new fields in dat_rare_with_treatment just so we are sure to not mess some indexing up. 
dat_rare_with_treatment$time <- as.factor(gsub("(T[123]).*","\\1", dat_rare_with_treatment$Time_herb))
dat_rare_with_treatment$herbicide <- gsub("T[123]_(.*)","\\1", dat_rare_with_treatment$Time_herb)
dat_rare_with_treatment$herbicide <- relevel(as.factor(dat_rare_with_treatment$herbicide), ref = "Non-Treated")

dat_rare_with_treatment <-na.omit(dat_rare_with_treatment)

regressionresults <- data.frame(matrix(ncol = 10, nrow = 1))

#Use a zero-inflated beta regression to deal with our proportion data that lie within the interval [0,1)
# See this nice stackexchange post for more: https://stats.stackexchange.com/questions/309047/zero-inflated-beta-regression-using-gamlss-for-vegetation-cover-data

k <- 1
for(i in names(otu_tab_bac)){
  reg <- gamlss(dat_rare_with_treatment[,names(dat_rare_with_treatment) == i] ~ 
              dat_rare_with_treatment$time + dat_rare_with_treatment$herbicide + dat_rare_with_treatment$Total_Weed_Veg, family = BEZI)
  regressionresults[k, 1] <- i
  regressionresults[k, 2:10] <- c(coef(reg, what = "mu"),
                                  Rsq(reg)) #generalised R-squared of Nagelkerke (1991) for a GAMLSS mode
  k <- k + 1
}

#NOTE: non-treated and T1 are the reference conditions. 
names(regressionresults) <- c("taxon", "intercept", "timeT2", "timeT3", "Atrazine-Mesotrione", "Dicamba","Handweeded","Glyphosate", "Total Weedy Veg", "Rsquared")
max(regressionresults$Rsquared)
summary(regressionresults$Rsquared)


tax_table_from_regression_results <- data.frame(tax_table(phyloseqD))
tax_table_from_regression_results$Taxonomy<- paste(tax_table_from_regression_results[,1],tax_table_from_regression_results[,2],tax_table_from_regression_results[,3],
                                                   tax_table_from_regression_results[,4], tax_table_from_regression_results[,5],tax_table_from_regression_results[,6],
                                                   tax_table_from_regression_results[,7], sep = ";")
#check random rows to make sure otus line up
rownames(tax_table_from_regression_results)[22]
regressionresults$taxon[22]
rownames(tax_table_from_regression_results)[45]
regressionresults$taxon[45]

regressionresults<-cbind(regressionresults, tax_table_from_regression_results$Taxonomy)

regressionresults <- regressionresults[order(regressionresults$Rsquared, decreasing = T),]
#write.csv(regressionresults, "/Users/gordoncuster/Desktop/Git_Projects/Herbicide_Microbes_PT1/Writing/Supplementary_Materials/Final Supp Materials April 22/16S_regression.csv")
write.csv(regressionresults, "/Users/gordoncuster/Desktop/Git_Projects/Herbicide_Microbes_PT1/Writing/Supplementary_Materials/Final Supp Materials April 22/ITS_regression.csv")


#extract top 20 taxa from R2
extract_names<-regressionresults$taxon[1:20]
extract_names<- str_replace(extract_names, pattern = "[.]",replacement = "=")

phyloseqD_filt<-subset_taxa(phyloseqD, taxa_names(phyloseqD) %in% extract_names)
phyloseqD_filt_trt_time <- merge_samples(phyloseqD_filt, group = "herb_time")

micrUBIfuns::plot_taxa_heatmap(phyloseqD_filt_trt_time, rm_na = F,scale_by = "taxa", cluster_rows=F, cluster_columns = T , row_dend_reorder =TRUE)


