library(dplyr)
library(scales)
library(data.table)
library("stringr")
library(tidyverse)
library(ggplot2)
library(janitor)


###### merge read depth data ######
setwd('/home/ashley/Documents/amplicon_sequencing/human_genotyping/CP1_round_1/pool_1/')
pool_1<- read.table("depth_info_pool_1.txt", comment="", header=TRUE)
pool_1<-as.data.frame(pool_1)

setwd('/home/ashley/Documents/amplicon_sequencing/human_genotyping/CP1_round_1/pool_2/')
pool_2<- read.table("depth_info_pool_2.txt", comment="", header=TRUE)
pool_2<-as.data.frame(pool_2)

setwd('/home/ashley/Documents/amplicon_sequencing/human_genotyping/CP1_round_1/pool_3')
pool_3<- read.table("depth_info_pool_3.txt", comment="", header=TRUE)
pool_3<-as.data.frame(pool_3)

setwd('/home/ashley/Documents/amplicon_sequencing/human_genotyping/CP1_round_1/pool_4')
pool_4<- read.table("depth_info_pool_4.txt", comment="", header=TRUE)
pool_4<-as.data.frame(pool_4)

setwd('/home/ashley/Documents/amplicon_sequencing/human_genotyping/CP1_round_1/pool_5')
pool_5<- read.table("depth_info_pool_5.txt", comment="", header=TRUE)
pool_5<-as.data.frame(pool_5)

setwd('/home/ashley/Documents/amplicon_sequencing/human_genotyping/CP1_round_2/pool_6')
pool_6<- read.table("depth_info_pool_6.txt", comment="", header=TRUE)
pool_6<-as.data.frame(pool_6)

setwd('/home/ashley/Documents/amplicon_sequencing/human_genotyping/CP1_round_2/pool_7')
pool_7<- read.table("depth_info_pool_7.txt", comment="", header=TRUE)
pool_7<-as.data.frame(pool_7)

setwd('/home/ashley/Documents/amplicon_sequencing/human_genotyping/CP1_round_2/pool_8')
pool_8<- read.table("depth_info_pool_8.txt", comment="", header=TRUE)
pool_8<-as.data.frame(pool_8)

setwd('/home/ashley/Documents/amplicon_sequencing/human_genotyping/CP1_round_2/pool_9')
pool_9<- read.table("depth_info_pool_9.txt", comment="", header=TRUE)
pool_9<-as.data.frame(pool_9)

setwd('/home/ashley/Documents/amplicon_sequencing/human_genotyping/CP1_round_2/pool_10')
pool_10<- read.table("depth_info_pool_10.txt", comment="", header=TRUE)
pool_10<-as.data.frame(pool_10)

depth <-merge(pool_1,pool_2,by=c("chrom","pos"))
depth <-merge(depth,pool_3,by=c("chrom","pos"))
depth <-merge(depth,pool_4,by=c("chrom","pos"))
depth <-merge(depth,pool_5,by=c("chrom","pos"))
depth <-merge(depth,pool_6,by=c("chrom","pos"))
depth <-merge(depth,pool_7,by=c("chrom","pos"))
depth <-merge(depth,pool_8,by=c("chrom","pos"))
depth <-merge(depth,pool_9,by=c("chrom","pos"))
depth <-merge(depth,pool_10,by=c("chrom","pos"))

#filter based on position

depth_1<-filter(depth, pos == 159204893 |pos == 143781321 | pos == 154536002 | pos == 154535277 | pos == 154533025 | pos == 154534440 |
                  pos == 154534125 | pos == 154533122 | 
                  pos == 5227002 | pos == 5227003 | pos == 5226943 | pos == 5226963 | pos == 5226925 )

is.na(depth_1) <- depth_1=="."
depth_1[is.na(depth_1)] <- 0
