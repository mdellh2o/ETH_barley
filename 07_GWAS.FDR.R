#===============================================================================
# authors: Matteo Dell'Acqua and Leonardo Caproni
# date: 08/2022
#-------------------------------------------------------------------------------
# Reference: DOI
#-------------------------------------------------------------------------------
# Description: Get significant SNPs using False Discovery Rate
#===============================================================================

#preliminary
rm(list=ls())
main <- "~/barley"
setwd(main)

#load libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(ape)
library("qvalue")
library(plyr)

# list files
files=list.files(path= "output/GWAS_MAF005/results/PC.10", pattern="*FarmCPU.csv", full.names = TRUE) 
# files
pattern="*FarmCPU.csv"
pat=".FarmCPU.csv"

# q values
thrfdr<-0.05 #define FDR thr

outlist <- list()
minsig<-c()

for(i in 1:length(files)){ #i=1
  tmp <- read.csv(files[i], header=T)
  
  colnames(tmp)[8]<-"pvalue"
  tmp$qvalue<-qvalue(tmp[,8])$qvalue
  
  hit<-tmp[which(tmp$qvalue<thrfdr),]
  
  
  if(nrow(hit)>0){
    #store mininum sig
    minsig[i]<-max(hit$pvalue)
    names(minsig)[i]<-files[i]
    
    out<-data.frame(trait=files[i], hit)
    outlist[[i]] <- out
  }#if hit
}

minsig<-as.data.frame(tibble::enframe(minsig))
colnames(minsig)[1] <- "trait"

#minsig <- tibble::rownames_to_column(minsig, "trait")

minsig[,1]<-sub("output/GWAS_MAF005/results/PC.10/", "", minsig[,1])
minsig[,1]<-sub(".FarmCPU.csv", "", minsig[,1])

write.table(minsig, file = "output/GWAS_MAF005/PC.10.min.sig.FDR.005.txt", quote = F, row.names = F, sep = "\t")

#Significat association after FDR
sigFDR<-do.call(rbind, outlist)
dim(sigFDR)
table(sigFDR[,1])

#fix values
sigFDR[,1]<-sub("output/GWAS_MAF005/results/PC.10/", "", sigFDR[,1])
sigFDR[,1]<-sub(".FarmCPU.csv", "", sigFDR[,1])


#split environmental and phenotypic associations
bio.res <- subset(sigFDR,grepl("bio",sigFDR$trait))
phen.res <- subset(sigFDR,grepl("D",sigFDR$trait))

length(unique(bio.res$SNP))
length(unique(phen.res$SNP))

############################################################

write.table(sigFDR, file="output/supplementary.table.4.REV1.MAF005.PC.10.fdr005.txt", quote = F, row.names = F, sep = "\t")