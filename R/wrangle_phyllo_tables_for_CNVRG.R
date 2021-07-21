rm(list=ls())
otus <- read.csv("./data/otutabITS_cleaned_phyllosphere", stringsAsFactors = F)
head(otus)

#The purpose of bringing in the coligo dfs is to get the ISD counts
coligos <- read.table("./data/coligo_phyllosphere//coligoITs_17sep20", header = T)
coligos2 <- read.table("./data/coligo_phyllosphere//coligoits_harrison_17sep20", header = T)

#Figure out the corn samples
coligos2 <- coligos2[, -(grep("e[pn]", names(coligos2)))]
coligos2 <- coligos2[, -(grep("_\\d+_\\d+_", names(coligos2)))]
#Doesn't appear to be any for its in this one
rm(coligos2)

names(coligos) %in% names(otus)

#Exdtract ISD counts by sample
otus[1 + length(otus[,1]), ] <- 0
for(i in 1:length(coligos)){
  if(names(coligos)[i] %in% names(otus)){
    otus[length(otus[,1]), which(names(otus) == names(coligos)[i])] <- 
      coligos[2, i]
  }
}

#QC, this shows us the ISD counts for each sample
otus[length(otus[,1]), ]

#Ok, split into time and add treatment?
test <- read.csv("./data/Herb_Pt_1_IC_anion_cation_import_to_R.csv", stringsAsFactors = F)
test[1:3,1:3]

test$Sample_ID
names(otus) <- gsub('ITS_\\w+_\\w+_',"",names(otus))
test$Sample_ID <- gsub("[_ ]", "", test$Sample_ID)
test$Sample_ID <- gsub("S", "P", test$Sample_ID)

test$Sample_ID %in% names(otus)

otus_t <- data.frame(names(otus)[2:length(otus)], t(otus[,2:length(otus)]))
names(otus_t) <- c("samples", otus$OTUID)
head(otus_t)
names(otus_t)[length(otus_t)] <- "ISD"
otus_t$samples <- gsub("p", "P", otus_t$samples)
otus_t$samples <- gsub("103B", "103", otus_t$samples)
#otus_t$samples <- gsub("(103P[123])", "\\1T2", otus_t$samples)

#What is this first sample that is called "ITS"
#Bring in the file names, extract sample names and figure out which are missing
#filenames made like this: 
#/project/microbiome/data/seq/psomagen_17sep20_novaseq2/tfmergedreads/ITS/corn/corn.* > filenames 
nms <- read.table("filenames", stringsAsFactors = F)
nms$V2 <- gsub(".*corn\\.(.*)\\.ITS.*Harrison.ITS.*","\\1" , nms$V1)
head(nms$V2)

otus_t$samples %in%  nms$V2 # all good except for the wacky "ITS"
#Therefore, the missing thing must be from this library: cd /project/microbiome/data/seq/psomagen_9oct20_novaseq3/tfmergedreads/ITS/sm18-herb/
#The only thing in there is 101p1, is this what is in the samples but not in the names?

setdiff(otus_t$samples, test$Sample_ID)

otus_t$samples[grep("101",otus_t$samples)] #101p1 is missing and has a weird name in the 9oct library, so I think this is it.

otus_t$samples[otus_t$samples == "ITS"] <- "101P1"

#What is in the metadata but not hte seq data
setdiff(test$Sample_ID, otus_t$samples)
grep("203",otus_t$samples, value =T)
#So where are these samples?
#are there odd dupes in the metadata?
grep("203",test$Sample_ID, value =T)

#should be 2 for PCR duplicates
#101p1 got combined bc original had same names for files I think, during otu table making. 
table(otus_t$samples) == 2

#Combine PCR duplicates, then add treatment info, then split by time


combo_dataframe <- data.frame(matrix(nrow = length(otus_t$samples)/2, ncol = length(otus_t)))
k <- 1
for(i in unique(otus_t$samples)){
  combo_dataframe[k, ] <- c(i, colSums(otus_t[otus_t$samples == i, 2:length(otus_t)]))
  k <- k + 1
}

dim(combo_dataframe)
#remove NA rows
combo_dataframe <- combo_dataframe[!is.na(combo_dataframe$X1),]
names(combo_dataframe) <- names(otus_t)

combo_dataframe$samples %in% test$Sample_ID

mergeddat <- merge(test, combo_dataframe, by.x = "Sample_ID", by.y = "samples")
dim(mergeddat)
#Ok, now we break up by treatment time

t2 <- mergeddat[grep("T2",mergeddat$Sample_ID),]
t3 <- mergeddat[grep("T3",mergeddat$Sample_ID),]
t1 <- mergeddat[-c(grep("T3",mergeddat$Sample_ID), grep("T2",mergeddat$Sample_ID)),]

dim(t1)
dim(t2)
dim(t3)

write.csv(t1, file = "./forModeling_ITS_otuTables/phyllo_ITS_time1_otu.csv", row.names = F)
write.csv(t2, file = "./forModeling_ITS_otuTables/phyllo_ITS_time2_otu.csv", row.names = F)
write.csv(t3, file = "./forModeling_ITS_otuTables/phyllo_ITS_time3_otu.csv", row.names = F)

########################################################
# Do over for 16s

rm(list=ls())
otus <- read.csv("./data/otutab16s_cleaned_phyllosphere", stringsAsFactors = F)
head(otus)

#The purpose of bringing in the coligo dfs is to get the ISD counts
coligos <- read.table("./data/coligo_phyllosphere//coligo16s_17sep20", header = T)
coligos2 <- read.table("./data/coligo_phyllosphere//coligo16s_harrison_17sep20", header = T)

#Figure out the corn samples
coligos2 <- coligos2[, -(grep("e[pn]", names(coligos2)))]
coligos2 <- coligos2[, -(grep("_\\d+_\\d+_", names(coligos2)))]
names(coligos2) #101p2, 301p3

names(coligos) %in% names(otus)

#Exdtract ISD counts by sample
otus[1 + length(otus[,1]), ] <- 0
for(i in 1:length(coligos)){
  if(names(coligos)[i] %in% names(otus)){
    otus[length(otus[,1]), which(names(otus) == names(coligos)[i])] <- 
      coligos[1, i]
  }else if(names(coligos2)[i] %in% names(otus)){
    otus[length(otus[,1]), which(names(otus) == names(coligos2)[i])] <- 
      coligos2[20, i]
  }
}

grep("105", names(coligos), value = T)
grep("105", names(otus), value = T)

#QC, this shows us the ISD counts for each sample
otus[length(otus[,1]), ]

#Ok, split into time and add treatment?
test <- read.csv("./data/Herb_Pt_1_IC_anion_cation_import_to_R.csv", stringsAsFactors = F)
test[1:3,1:3]

test$Sample_ID
names(otus) <- gsub('rna16S_\\w+_\\w+_',"",names(otus))
test$Sample_ID <- gsub("[_ ]", "", test$Sample_ID)
test$Sample_ID <- gsub("S", "P", test$Sample_ID)

test$Sample_ID %in% names(otus)
#a few dont match, maybe T1s?
setdiff(test$Sample_ID, names(otus))
setdiff(names(otus), test$Sample_ID)
names(otus) <- gsub("p", "P", names(otus))
names(otus) <- gsub("103B", "103", names(otus))

otus_t <- data.frame(names(otus)[2:length(otus)], t(otus[,2:length(otus)]))
names(otus_t) <- c("samples", otus$OTUID)
head(otus_t)
names(otus_t)[length(otus_t)] <- "ISD"

#What is this first sample that is called "16S"
#Bring in the file names, extract sample names and figure out which are missing
#filenames made like this: 
#/project/microbiome/data/seq/psomagen_17sep20_novaseq2/tfmergedreads/16S/corn/corn.* > filenames16 
# nms <- read.table("filenames16", stringsAsFactors = F)
# nms$V2 <- gsub(".*corn\\.(.*)\\.16S.*Harrison.16S.*","\\1" , nms$V1)
# head(nms$V2)
# 
# otus_t$samples %in%  nms$V2
# 
# setdiff(otus_t$samples, test$Sample_ID)
# setdiff(test$Sample_ID, otus_t$samples)
# 
# otus_t$samples[grep("101",otus_t$samples)] #101p1 is missing and has a weird name in the 9oct library, so I think this is it.
otus_t$samples <- as.character(otus_t$samples)
otus_t$samples[otus_t$samples == "X16S"] <- "101P1"

#should be 2 for PCR duplicates
#101p1 got combined bc original had same names for files I think, during otu table making. 
table(otus_t$samples) == 2

#Combine PCR duplicates, then add treatment info, then split by time

combo_dataframe <- data.frame(matrix(nrow = length(otus_t$samples)/2, ncol = length(otus_t)))
k <- 1
for(i in unique(otus_t$samples)){
  combo_dataframe[k, ] <- c(i, colSums(otus_t[otus_t$samples == i, 2:length(otus_t)]))
  k <- k + 1
}

dim(combo_dataframe)
#remove NA rows
combo_dataframe <- combo_dataframe[!is.na(combo_dataframe$X1),]
names(combo_dataframe) <- names(otus_t)

combo_dataframe$samples %in% test$Sample_ID

mergeddat <- merge(test, combo_dataframe, by.x = "Sample_ID", by.y = "samples")
dim(mergeddat)
#Ok, now we break up by treatment time

t2 <- mergeddat[grep("T2",mergeddat$Sample_ID),]
t3 <- mergeddat[grep("T3",mergeddat$Sample_ID),]
t1 <- mergeddat[-c(grep("T3",mergeddat$Sample_ID), grep("T2",mergeddat$Sample_ID)),]

dim(t1)
dim(t2)
dim(t3)

write.csv(t1, file = "./forModeling_16S_otuTables/phyllo_16S_time1_otu.csv", row.names = F)
write.csv(t2, file = "./forModeling_16S_otuTables/phyllo_16S_time2_otu.csv", row.names = F)
write.csv(t3, file = "./forModeling_16S_otuTables/phyllo_16S_time3_otu.csv", row.names = F)

