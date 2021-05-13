dat1 <- read.csv("post_modeling_ML_estimates/16S_time1_otu.csv_p_estimates.csv",
                 stringsAsFactors = F, header = T)
dat1[1:3,1:5]

dat2 <- read.csv("post_modeling_ML_estimates/16S_time1_otu.csv_p_estimatesHMC.csv",
                 stringsAsFactors = F, header = T)
dat2[1:3,1:5]

k <- 1
est <- NA
for(i in names(dat1)[4:length(dat1)]){
  if(i %in% names(dat2)){
   c <- cor.test(dat1[,names(dat1)==i], dat2[, names(dat2)==i])
   est[k] <- c$estimate
   k  <- k + 1
  }
}
hist(est)
min(est) #wow
length(est)
