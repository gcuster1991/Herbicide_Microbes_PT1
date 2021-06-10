#Calculate the number of taxa that differ in relative abundance
#among treatment groups. Cycle through each time step and locus. 
rm(list=ls())

library(CNVRG)

inputs <- list.files(".", pattern = "*HMC\\.Rdata")

summary_by_group <- data.frame(matrix(nrow = 20, ncol = 2))

for(i in 1:length(inputs)){
  load(inputs[i])
  diffs <- diff_abund(model_out = modelOut, countData = dat3)
  write.csv(diffs, file = paste(gsub("Rdata","",inputs[i]), "_diff_estimates.csv", sep = ""), row.names = F)
  for(j in 1:length(diffs[,1])){
    summary_by_group[j,2] <- length(which(diffs[j,2:length(diffs)] > 0.95 | diffs[j,2:length(diffs)] < 0.05))
  }
  summary_by_group[,1] <- diffs[,1]
  names(summary_by_group) <- c("comparison", "num_taxa")
  write.csv(summary_by_group, file = paste(gsub("Rdata","",inputs[i]), "_diff_Summary.csv", sep = ""), row.names = F)
  rm(modelOut, dat3)
}
