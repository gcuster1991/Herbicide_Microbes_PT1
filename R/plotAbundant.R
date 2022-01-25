
library(phyloseq)
library(CNVRG)
library(tidyverse)

#moved the raw objects to the 16S and ITS data folders of the repo. 
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
treatments_full <- dat[,1]
dat <- dat[,-1]
dat[1:5,1:10]

dat <- dat[,rev(order(colSums(dat[,1:length(dat)])))]
microbenames <- names(dat)

#############
# Important #
#############
# Rarify or not here. I am calculating probability of occurrence. This is the max likelihood multinomial parameters

dat <- vegan::drarefy(dat, 500)

#Figure out the median counts per treatment for the top n
n <- 25 #pick how many taxa to plot

median_top10 <- data.frame(matrix(ncol = n, nrow = 15))

for(i in 1:n){
  median_top10[,i] <- aggregate(dat[,i]~treatments_full, FUN=median)[,2]
}

median_top10 <- data.frame(aggregate(dat[,i]~treatments_full, FUN=median)[,1], 
  median_top10)

names(median_top10) <- c("treatment", microbenames[1:n])
#Assign colors for plotting
median_top10$col <- "red"

median_top10$col[grep("Aatrex",median_top10$treatment)] <- "blue"
median_top10$col[grep("Clarity",median_top10$treatment)] <- "green"
median_top10$col[grep("Hand",median_top10$treatment)] <- "orange"
median_top10$col[grep("Non-Treated",median_top10$treatment)] <- "black"
median_top10$col[grep("Round",median_top10$treatment)] <- "purple"


median_top10$herbicide <- gsub("T[1-3]_(.*)", "\\1", median_top10$treatment)
median_top10$time <- gsub("(T[1-3])_.*", "\\1", median_top10$treatment)

pdf(width = 6, height = 8, file = "abundantTaxa.pdf")
par(mfrow = c(3,2))
for(j in unique(median_top10$herbicide)){
  plot(median_top10[median_top10$herbicide == j,2],
  type = 'b',
  ylim = c(0,1), 
  ylab = j, 
  las = 2,
  xlab = "Time Point",
  frame.plot = F,
  xpd = NA,
  col = median_top10$col[grep(j,median_top10$treatment)] ,
  xaxt = "n")
axis(side = 1, at = c(1,2,3), labels = c(1,2,3))

  for(i in 3:11){
    lines(x = c(1,2,3),
        y = median_top10[median_top10$herbicide == j,i],
        type = 'b', 
        col = median_top10$col[grep(j,median_top10$treatment)] ,
        xpd = NA)
  }
}
dev.off()

#Which is the most variable? This may not work if n changes from 25.

apply(median_top10[,2:(n+1)], MARGIN =2, FUN = sd)

#Looks like the 4th is the most variable. check. Note that I am using 4+1 here to index.
names(median_top10)[5]

par(mfrow = c(2,3))
for(j in unique(median_top10$herbicide)){
  plot(median_top10[median_top10$herbicide == j,5],
       type = 'b',
       ylim = c(0,1), 
       ylab = j, 
       las = 2,
       xlab = "Time Point",
       frame.plot = F,
       xpd = NA,
       col = median_top10$col[grep(j,median_top10$treatment)] ,
       xaxt = "n")
  axis(side = 1, at = c(1,2,3), labels = c(1,2,3))

}
#Yup centroid 8 is very variable.


#further examination into taxonomy of those top taxa
centroid_list<-names(median_top10)[2:(ncol(median_top10)-3)]
centroid_list<-str_replace(centroid_list, pattern = "[[.]]", replacement = "=")

table(rownames(tax_table(HerbPt1_PS_Orig)) %in% centroid_list)
top25_abund<-subset_taxa(physeq = HerbPt1_PS_Orig,  rownames(tax_table(HerbPt1_PS_Orig)) %in% centroid_list)
top25_abund_tax<-data.frame(tax_table(top25_abund))
write.csv(top25_abund_tax, "../../Top25_16S_tax.csv")

#look at diversity within these plots. 
sort(sample_sums(top25_abund))
sample_data(top25_abund)$herb_time <- paste(sample_data(top25_abund)$Herbicide,sample_data(top25_abund)$Time, sep = "_" )
top25_abund_rare<-rarefy_even_depth(top25_abund, sample.size = 500)

plot_richness(top25_abund_rare, x = "herb_time", measures = c("Shannon", "Simpson", "Observed"), color = "Herbicide") + geom_boxplot()
ggsave("../../../Figures/16S_top25_abund.pdf")

table(data.frame(tax_table(top25_abund_rare))[,2])
table(data.frame(tax_table(top25_abund_rare))[,3])
table(data.frame(tax_table(top25_abund_rare))[,4])
table(data.frame(tax_table(top25_abund_rare))[,5])
table(data.frame(tax_table(top25_abund_rare))[,6])

