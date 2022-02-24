rm(list=ls())
library(phyloseq)

#moved the raw objects to the 16S and ITS data folders of the repo. 
load("./Herb16S_PS_Orig.RData")

#21300 taxa in original
HerbPt1_PS_Orig
#remove all non-bacterial reads
#21028 taxa in bacteria only ps object. 272 non-bacterial taxa 
ps_bac<-subset_taxa(HerbPt1_PS_Orig, Kingdom == "k:Bacteria")
#create seperate dataframe of archea
#101 archeal taxa
ps_archea<-subset_taxa(HerbPt1_PS_Orig, Kingdom=="k:Archaea")


otu_tab_bac<-data.frame(t(otu_table(ps_bac)))
metadata_bac<-data.frame(sample_data(ps_bac))
metadata_bac$Time_herb<-paste(metadata_bac$Time, metadata_bac$Herbicide, sep = "_")

rownames(otu_tab_bac) == rownames(metadata_bac)

dat<-cbind(metadata_bac$Time_herb, otu_tab_bac)
names(dat)[names(dat) == 'metadata_bac$Time_herb'] <- 'Time_herb'
dat[1:10]

#Order by treatment
dat <- dat[order(dat$Time_herb),]
treatments_full <- dat[,1]
dat <- dat[,-1]
dat[1:5,1:10]

dat <- dat[,rev(order(colSums(dat[,1:length(dat)])))]
microbenames <- names(dat)

#############
# Important #
#############
# Rarify or not here. I am calculating probability of occurrence given n observations

dat_rare <- data.frame(vegan::drarefy(dat, 500))
names(dat_rare)

dat_rare_with_treatment <- data.frame(treatments_full, dat_rare)

T1 <- dat_rare_with_treatment[grep("T1", dat_rare_with_treatment$treatments_full),]
T1 <- T1[,order(colSums(T1[,2:length(T1)]), decreasing = T)]

T2 <- dat_rare_with_treatment[grep("T2", dat_rare_with_treatment$treatments_full),]
T2 <- T2[,order(colSums(T2[,2:length(T2)]), decreasing = T)]

T3 <- dat_rare_with_treatment[grep("T3", dat_rare_with_treatment$treatments_full),]
T3 <- T3[,order(colSums(T3[,2:length(T3)]), decreasing = T)]

#Figure out the median counts per treatment for the top n
medianFinder <- function(x, n){
  median_top10 <- data.frame(matrix(ncol = n, nrow = 5))
  
  for(i in 2:(n+1)){
    median_top10[,i-1] <- aggregate(x[,i] ~ 
                                      as.factor(x$treatments_full), 
                                    FUN=median)[,2]
  }
  median_top10 <- data.frame(unique(x[,1]), median_top10)
  names(median_top10) <- names(x)[1:(n+1)]
  return(median_top10)
}

medianFinder(x = T1, 25)
medianFinder(T2, 25)
medianFinder(T3, 25)

#Next need to build models for each taxon and output the results

#Get the unique abundant taxa
tomodel <- unique(c(names(medianFinder(x = T1, 25)),
              names(medianFinder(T2, 25)),
                names(medianFinder(T3, 25))))
tomodel <- tomodel[-1] #remove treatments_full

#Make a few new fields in dat_rare_with_treatment just so we are sure to not mess some indexing up. 
dat_rare_with_treatment$time <- gsub("(T[123]).*","\\1", dat_rare_with_treatment$treatments_full)
dat_rare_with_treatment$herbicide <- gsub("T[123]_(.*)","\\1", dat_rare_with_treatment$treatments_full)
dat_rare_with_treatment$herbicide <- relevel(as.factor(dat_rare_with_treatment$herbicide), ref = "Non-Treated")

regressionresults <- data.frame(matrix(ncol = 10, nrow = 1))

k <- 1
for(i in tomodel){
  #This is not quite right if the response is between 0 and 1
  reg <- lm(dat_rare_with_treatment[,names(dat_rare_with_treatment) == i]~ 
     dat_rare_with_treatment$time + dat_rare_with_treatment$herbicide)
     regressionresults[k, 1] <- i
     regressionresults[k, 2:10] <- c(summary(reg)$coefficients[,1],
                              summary(reg)$r.squared, 
                              summary(reg)$adj.r.squared)
     k <- k + 1
}

names(regressionresults) <- c("taxon", "intercept", "timeT2", "timeT3", "Aatrex", "Clarity","Hand","Round Up", "R squared", "Adj. R squared")
#Note that R2 is less thatn 10% for ALL taxa. 
regressionresults
