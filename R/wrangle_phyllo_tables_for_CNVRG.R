otus <- read.csv("./otutabITS_cleaned", stringsAsFactors = F)
head(otus)

coligos <- read.table("./coligoITs_17sep20", header = T)
coligos2 <- read.table("./coligoits_harrison_17sep20", header = T)

#Figure out the corn ones
coligos2 <- coligos2[, -(grep("e[pn]", names(coligos2)))]
coligos2 <- coligos2[, -(grep("_\\d+_\\d+_", names(coligos2)))]
#Doesn't appear to be any for its in this one
rm(coligos2)

names(coligos) %in% names(otus)

otus[1 + length(otus[,1]), ] <- 0
for(i in 1:length(coligos)){
  if(names(coligos)[i] %in% names(otus)){
    otus[length(otus[,1]), which(names(otus) == names(coligos)[i])] <- 
      coligos[2, i]
  }
}

#QC
otus[length(otus[,1]), ]

#Ok, split into time ?
test <- read.csv("forModeling_ITS_otuTables/ITS_time1_otu.csv")
test[1:3,1:3]
test$X
test$treatment

#Need to spit into times and add treatment