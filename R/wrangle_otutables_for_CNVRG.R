dat <- read.csv("./data/otutable16s.esv_withISD_and_mtDNA_cpDNA_etc.csv")
head(dat)

time1 <- read.csv("./forModeling_16s_otuTables/16S_time1_otu.csv", header = T)
time2 <- read.csv("./forModeling_16s_otuTables/16S_time2_otu.csv", header = T)
time3 <- read.csv("./forModeling_16s_otuTables/16S_time3_otu.csv", header = T)

t1 <- as.character(time1$X)
t2 <- as.character(time2$X)
t3 <- as.character(time3$X)

names(dat) <- gsub(".*_([A-Za-z0-9]*$)","\\1",names(dat))
write.csv(dat[,names(dat) %in% t1],
          file = "./forModeling_16s_otuTables/16S_time1_otu.csv", 
          row.names = F)
write.csv(dat[,names(dat) %in% t2],
          file = "./forModeling_16s_otuTables/16S_time2_otu.csv", 
          row.names = F)
write.csv(dat[,names(dat) %in% t3],
          file = "./forModeling_16s_otuTables/16S_time3_otu.csv", 
          row.names = F)
#Do for ITS
rm(list=ls())
dat <- read.csv("./data/otutableISD.esv.ITS.ISD_mt_cpDNA_etc.csv")
head(dat)

time1 <- read.csv("./forModeling_ITS_otuTables/ITS_time1_otu.csv", header = T)
time2 <- read.csv("./forModeling_ITS_otuTables/ITS_time2_otu.csv", header = T)
time3 <- read.csv("./forModeling_ITS_otuTables/ITS_time3_otu.csv", header = T)

t1 <- as.character(time1$X)
t2 <- as.character(time2$X)
t3 <- as.character(time3$X)

names(dat) <- gsub(".*_([A-Za-z0-9]*$)","\\1",names(dat))
write.csv(dat[,names(dat) %in% t1],
          file = "./forModeling_ITS_otuTables//ITS_time1_otu.csv", 
          row.names = F)
write.csv(dat[,names(dat) %in% t2],
          file = "./forModeling_ITS_otuTables/ITS_time2_otu.csv", 
          row.names = F)
write.csv(dat[,names(dat) %in% t3],
          file = "./forModeling_ITS_otuTables/ITS_time3_otu.csv", 
          row.names = F)
