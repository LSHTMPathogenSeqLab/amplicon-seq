library(dplyr)
library(scales)
library(data.table)
library("stringr")
library(tidyverse)
library(ggplot2)
library(janitor)


# load in depth_info.txt output from amplicon sequencing pipeline
depth <- read.table("/home/ashley/Documents/amplicon_sequencing/Ngodhe_MDA/Ngodhe_8/depth_info.txt", comment="", header=TRUE)
depth<-as.data.frame(depth)

## convert depth txt file into correct csv format

# filter by desired P. falciparum resistance-associated snp positions (from DHFR, DHPS, MDR1, CRT, AP2MU, AND K13)
# can use any snp positions you want for any species (can also use gene, chrom, or sample ID)

depth<-filter(depth, pos == 748239 | pos == 748262 | pos == 748410 | pos == 748577 | pos == 549681 | pos == 549685 | pos == 549993 | pos == 550117 | pos == 550212 | 
                      pos == 958071 | pos == 958145 | pos == 958440 | pos == 958484 | pos == 960702 | pos == 960990 | pos == 961013 | pos == 961104 | pos == 961625 | 
                      pos == 403613 |pos == 403617 | pos == 403620 | pos == 403621 | pos == 403625 | 
                      pos == 718391 | pos == 718433 | 
                      pos == 1725661 | pos == 1725626 | pos == 1725571 | pos == 1725520 | pos == 1725382 | 
                      pos == 1725370 | pos == 1725340 | pos == 1725316 | pos == 1725259 | pos ==1725266)
is.na(depth) <- depth=="."
depth[is.na(depth)] <- 0

depth <- depth %>% select(-c(chrom, pos, gene))

depth<-as.data.frame(t(depth))
depth<- depth %>% row_to_names(row_number = 1)
rownames(depth) <- 1:nrow(depth)

write.csv(depth, "/home/ashley/Documents/amplicon_sequencing/Ngodhe_MDA/Ngodhe_8/depth_info.csv", row.names=FALSE, quote=FALSE)

#load in new csv
depth<-read.csv("/home/ashley/Documents/amplicon_sequencing/Ngodhe_MDA/Ngodhe_8/depth_info.csv",header=T,sep=',')

# plot depth of resistance snps without k13 snps (tend to have a much higher read depth and can skew the plot)
# create data label list - we will plot the box plot without data labels then add them in manually so we ca adjust size and angle of the text

lablist<-as.vector(c("N51I", "C59R", "S108N", "I164L", "S436H", "G437A", "K540E", "A581G", "A613S",
                     "F61Y", "N86Y", "Y184F", "T199S", "F938Y", "S1034C", "N1042D", "F1072Y", "D1246Y",
                     "C72S","V73-", "M74I", "N75E", "K76T", "R146K", "S160N"))

pdf("/home/ashley/Documents/amplicon_sequencing/Ngodhe_MDA/plots/read_depth_no_k13_NG8.pdf")
boxplot(depth %>% select(c(Asn51, Cys59, Ser108, Ile164, Ser436, Gly437, Lys540, Ala581, Ala613,
                         Phe61, Asn86, Tyr184, Thr199, Phe938, Ser1034, Asn1042, Phe1072, Asp1246, Met74,
                         Cys72, Val73, Asn75, Lys76, Arg146, Ser160)), pars  =  list(xaxt = "n"), ylab = "Read Depth", 
        col = c('green', 'green', 'green', 'green', "yellow", "yellow", "yellow", "yellow", "yellow", 
                "orange", "orange", "orange", "orange", "orange", "orange", "orange", "orange", "orange", 
                "blue", "blue", "blue", "blue", "blue", "purple", "purple"),
        pch=20)

# change 1:25 depending on the number of snps you would like to plot ( I am plotting 25 snps here)
axis(1, at=c(1:25), labels = FALSE)

# can adjust size of label font with cex, offset will move the labels up and down, srt will change the angle of the text - currently set to diagonal
text(c(1:25), par("usr")[3] - 1, labels = lablist, srt = 55, pos = 1, offset= 1.6, xpd = TRUE, cex = 0.9)

legend("topleft",inset=c(0,0), col=c('green', "yellow", "orange", "blue", "purple"),pch=20,leg=c("PfDHFR","PfDHPS", "PfMDR1", "PfCRT", "PfAP2MU"), cex=0.9)
dev.off()



# k13 only
lablist<-as.vector(c("F446I", "N458Y", "M476I", "Y493H", "R539T", "I543T", "P553L", "R561H", "A578S", "C580Y"))

pdf("/home/ashley/Documents/amplicon_sequencing/Ngodhe_MDA/plots/read_depth_k13_NG8.pdf")
boxplot(depth %>% select(c(Phe446, Asn458, Met476, Tyr493, Arg539, Ile543, Pro553, Arg561, Ala578, Cys580)), pars  =  list(xaxt = "n"), ylab = "Read Depth", 
        col = c("red", "red", "red", "red", "red", "red", "red", "red", "red", "red"),
        pch=20)

axis(1, at=c(1:10), labels = FALSE)
text(c(1:10), par("usr")[3] - 1, labels = lablist, srt = 55, pos = 1, offset= 1.6, xpd = TRUE, cex = 0.8)
legend("topleft",inset=c(0,0), col= "red",pch=20,leg="PfK13", cex=0.9)
dev.off()
