
library(phyloseq)
library(CNVRG)
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

#Subset to those with 10 or more reads across samples
dat3 <- dat2[, c(1, 1+ which(colSums(dat2[,2:length(dat2)]) > 10))]
#table(colSums(dat3[, 2:length(dat3)]) > 3)

#Add a one to avoid taking the log of zero during modeling
dat3[,2:length(dat3)] <-  1 + dat3[,2:length(dat3)]

dim(dat3) #The y dimension here is way to big for this to run in a reasonable amount of time using HMC. VI might work, but
#the easiest thing to do is to subset the data so that taxa need to have more reads to be considered during modeling. 

#Determine how many taxa have low read counts. First, I will subtract one to undue the addition directly above. 
#then make some plots to help decide on a threshold, then subset the data. 
dat3[,2:length(dat3)] <-  dat3[,2:length(dat3)] -1

#The data are extremely skewed, as expected.
summary(colSums(dat3[,2:length(dat3)] ))
hist(colSums(dat3[,2:length(dat3)] ))
#Lets try a threshold of 50. This cuts down the number of things we have to model by half.
table(colSums(dat3[,2:length(dat3)]) > 50)
table(colSums(dat3[,2:length(dat3)]) > 100)
table(colSums(dat3[,2:length(dat3)]) > 50)

#remake dat3 and try modeling.
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


save.image("//project/microbiome/data/seq/combined/gcuster/otu/HerbPt1/16S/CNVRG_Herbpt1_16S_VI_model_100reads.RData")

ests <- extract_point_estimate(model_out = modelOut, countData = dat3,params = c("p","pi"))

forExport <- data.frame(treatments, dat3[,1],ests$pointEstimates_p)
names(forExport)[2] <- "sample"

save.image("//project/microbiome/data/seq/combined/gcuster/otu/HerbPt1/16S/CNVRG_Herbpt1_16S_VI_model_w_param_est.RData")
