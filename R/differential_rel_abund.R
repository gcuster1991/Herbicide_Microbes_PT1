#Calculate the number of taxa that differ in relative abundance
#among treatment groups. Cycle through each time step and locus. 
rm(list=ls())

library(CNVRG)

inputs <- list.files(".", pattern = "*HMC\\.Rdata")

summary_by_group <- data.frame(matrix(nrow = 20, ncol = 2))

for(i in 1:length(inputs)){
  load(inputs[i])
  transformed <- isd_transform(model_out = modelOut, countData = dat3, isd_index = which(names(dat3)=="ISD"))
  diffs <- diff_abund(model_out = transformed, countData = dat3)
  save(diffs, file = paste(inputs[i], "_diff_estimates.Rdata", sep = ""))
}
