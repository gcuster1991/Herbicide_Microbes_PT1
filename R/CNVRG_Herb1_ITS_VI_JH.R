

library(phyloseq)
library(CNVRG)
load("./HerpPT1_ITS_PS_Orig.RData")

HerbPt1_PS_Orig
#remove all non-fungal reads
ps_fungal<-subset_taxa(HerbPt1_PS_Orig, Kingdom == "d:Fungi")


otu_tab_fungal<-data.frame(t(otu_table(ps_fungal)))
metadata_fungal<-data.frame(sample_data(ps_fungal))
metadata_fungal$Time_herb<-paste(metadata_fungal$Time, metadata_fungal$Herbicide, sep = "_")

#QC
rownames(otu_tab_fungal) == rownames(otu_tab_fungal)

dat<-cbind(metadata_fungal$Time_herb, otu_tab_fungal)
names(dat)[names(dat) == 'metadata_fungal$Time_herb'] <- 'Time_herb'
dat[1:10]

#Order by treatment
dat <- dat[order(dat$Time_herb),]

#check for zero counts and rows with no data
table(rowSums(dat[,3:length(dat)]) > 0)
table(colSums(dat[,3:length(dat)]) > 0)

#Some OTUs are not present in these data. Remove them. 
dat2 <- dat[, c(TRUE, TRUE, colSums(dat[,3:length(dat)]) > 0)]

#check
table(colSums(dat2[,3:length(dat2)]) > 0 )
dat2[1:3,1:5]

#For ease, extracting treatment vector
treatments <- dat2$Time_herb
#remove treatment vector from dataframe
dat2 <- dat2[,-2]

#The data are extremely skewed, as expected.
summary(colSums(dat2[,2:length(dat2)] ))
#hist(colSums(dat2[,2:length(dat2)] ))


#Lets try a threshold of 50. This cuts down the number of things we have to model by three quarters ish
table(colSums(dat2[,2:length(dat2)]) > 50)
#table(colSums(dat2[,2:length(dat2)]) > 100)

#Make a new dataframe and add a one to it for modeling.
dat3 <- dat2[, c(1, 1+ which(colSums(dat2[,2:length(dat2)]) > 100))]
dat3[,2:length(dat3)] <-  1 + dat3[,2:length(dat3)]
dim(dat3)

modelOut <- cnvrg_VI(
  countData = dat3,
  starts = indexer(treatments)$starts,
  ends = indexer(treatments)$ends,
  #algorithm = "NUTS",
  #chains = 2,
  #burn = 500,
  #samples = 1500,
  #thinning_rate = 2,
  output_samples = 250,
  #  cores = 16,
  params_to_save = c("pi","p")
)
save.image("//project/microbiome/data/seq/combined/gcuster/otu/HerbPt1/ITS/CNVRG_Herbpt1_ITS_VI_model_100reads.RData")

ests <- extract_point_estimate(model_out = modelOut, countData = dat3,params = c("p","pi"))
forExport <- data.frame(treatments, dat3[,1],ests$pointEstimates_p)
names(forExport)[2] <- "sample"

save.image("//project/microbiome/data/seq/combined/gcuster/otu/HerbPt1/ITS/CNVRG_Herbpt1_ITS_VI_model_w_param_est.RData")
