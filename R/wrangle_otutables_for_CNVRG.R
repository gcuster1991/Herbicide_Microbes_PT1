rm(list=ls())
dat <- read.csv("./data/otutable16s.esv_withISD_and_mtDNA_cpDNA_etc.csv")
head(dat)

time1 <- read.csv("~/Desktop/Herbicide_Microbes_PT1//forModeling_16s_otuTables/16S_time1_otu.csv", header = T)
time2 <- read.csv("~/Desktop/Herbicide_Microbes_PT1//forModeling_16s_otuTables/16S_time2_otu.csv", header = T)
time3 <- read.csv("~/Desktop/Herbicide_Microbes_PT1//forModeling_16s_otuTables/16S_time3_otu.csv", header = T)

t1 <- as.character(time1$X)
t2 <- as.character(time2$X)
t3 <- as.character(time3$X)

#sanity check to see if the time points had duplicated sample names
intersect(t1, t2)
intersect(t1, t3)
intersect(t3, t2)

names(dat) <- gsub(".*_([A-Za-z0-9]*$)","\\1",names(dat))
duplicated(names(dat))

#Ok., dat has samples that were sequenced multiple times and need to be combined
dat2 <- data.frame(matrix(nrow = length(dat[,1]), ncol = 1))
k <- 1
for(i in unique(names(dat)[2:length(names(dat))])){
  if(length(which(names(dat) == i)) > 1){
    dat2[,k] <- rowSums(dat[,names(dat) == i])
  }else{
    dat2[,k] <- dat[,names(dat) == i]
  }
  names(dat2)[k] <- i
  k <- k + 1
}

t1dat <- data.frame(dat2[,names(dat2) %in% t1])
setdiff(names(t1dat), t1)

t2dat <- data.frame(dat2[,names(dat2) %in% t2])
setdiff(names(t2dat), t2)

t3dat <- data.frame(dat2[,names(dat2) %in% t3])
setdiff(names(t3dat), t3)

#transpose and then add treatement and otu names
t1dat <- data.frame(t(t1dat))
t2dat <- data.frame(t(t2dat))
t3dat <- data.frame(t(t3dat))
names(t1dat) <- dat$OTUID
names(t2dat) <- dat$OTUID
names(t3dat) <- dat$OTUID

t1dat$treatment <- time1$treatment
t2dat$treatment <- time2$treatment
t3dat$treatment <- time3$treatment

###

write.csv(t1dat,
          file = "./forModeling_16s_otuTables/16S_time1_otu.csv", 
          row.names = T)
write.csv(t2dat,
          file = "./forModeling_16s_otuTables/16S_time2_otu.csv", 
          row.names = T)
write.csv(t3dat,
          file = "./forModeling_16s_otuTables/16S_time3_otu.csv", 
          row.names = T)

#Do for ITS
rm(list=ls())
dat <- read.csv("./data/otutableISD.esv.ITS.ISD_mt_cpDNA_etc.csv")
head(dat)

time1 <- read.csv("~/Desktop/Herbicide_Microbes_PT1///forModeling_ITS_otuTables/ITS_time1_otu.csv", header = T)
time2 <- read.csv("~/Desktop/Herbicide_Microbes_PT1///forModeling_ITS_otuTables/ITS_time2_otu.csv", header = T)
time3 <- read.csv("~/Desktop/Herbicide_Microbes_PT1///forModeling_ITS_otuTables/ITS_time3_otu.csv", header = T)

t1 <- as.character(time1$X)
t2 <- as.character(time2$X)
t3 <- as.character(time3$X)

#sanity check to see if the time points had duplicated sample names
intersect(t1, t2)
intersect(t1, t3)
intersect(t3, t2)

names(dat) <- gsub(".*_([A-Za-z0-9]*$)","\\1",names(dat))
duplicated(names(dat))

#Ok., I think dat has samples that were sequenced multiple times and need to be combined
dat2 <- data.frame(matrix(nrow = length(dat[,1]), ncol = 1))
k <- 1
for(i in unique(names(dat)[2:length(names(dat))])){
  if(length(which(names(dat) == i)) > 1){
    dat2[,k] <- rowSums(dat[,names(dat) == i])
  }else{
    dat2[,k] <- dat[,names(dat) == i]
  }
  names(dat2)[k] <- i
  k <- k + 1
}

t1dat <- data.frame(dat2[,names(dat2) %in% t1])
setdiff(names(t1dat), t1)

t2dat <- data.frame(dat2[,names(dat2) %in% t2])
setdiff(names(t2dat), t2)

t3dat <- data.frame(dat2[,names(dat2) %in% t3])
setdiff(names(t3dat), t3)

#transpose and then add treatement and otu names
t1dat <- data.frame(t(t1dat))
t2dat <- data.frame(t(t2dat))
t3dat <- data.frame(t(t3dat))
names(t1dat) <- dat$OTUID
names(t2dat) <- dat$OTUID
names(t3dat) <- dat$OTUID

t1dat$treatment <- time1$treatment
t2dat$treatment <- time2$treatment
t3dat$treatment <- time3$treatment



write.csv(t1dat,
          file = "./forModeling_ITS_otuTables/ITS_time1_otu.csv", 
          row.names = T)
write.csv(t2dat,
          file = "./forModeling_ITS_otuTables/ITS_time2_otu.csv", 
          row.names = T)
write.csv(t3dat,
          file = "./forModeling_ITS_otuTables/ITS_time3_otu.csv", 
          row.names = T)

