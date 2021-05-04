rm(list=ls())
set.seed(666)
#install.packages("CNVRG")
library(rstan)
library(CNVRG)

#This takes STDIN so we call this script via the command line for different files
#We pass in a file to STDIN. Rscript CNVRG_modedling_script.R yourfile
#inargs[1] is then the path to your file. 

inargs <- commandArgs(trailingOnly = TRUE)

dat <- read.csv(inargs[1], 
                fill = T, header = T, stringsAsFactors = F)

#Order by treatment
dat <- dat[order(dat$treatment),]

#check for zero counts and rows with no data, then add one to all data
table(rowSums(dat[,3:length(dat)]) > 0)
table(colSums(dat[,3:length(dat)]) > 0)
#Some OTUs are not present in these data. Remove them. 
dat2 <- dat[, c(TRUE, TRUE, colSums(dat[,3:length(dat)]) > 0)]

#check
table(colSums(dat2[,3:length(dat2)]) > 0 )
dat2[1:3,1:5]

#For ease, extracting treatment vector
treatments <- dat2$treatment
#remove treatment vector from dataframe
dat2 <- dat2[,-2]

dat2[,2:length(dat2)] <-  1 + dat2[,2:length(dat2)]

#Commented out options are useful for cnvrg_HMC
modelOut <- cnvrg_VI(
  countData = dat2,
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
ests <- extract_point_estimate(modelOut = modelOut, countData = dat2,
                               treatments = length(unique(treatments)))
forExport <- data.frame(treatments, dat2[,1],ests$pointEstimates_p)
names(forExport)[2] <- "sample"

write.csv(forExport, file = paste(inargs[1], "_p_estimates.csv", sep = ""))
save.image(file = paste(inargs[1], ".Rdata", sep = ""))


