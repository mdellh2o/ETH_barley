#===============================================================================
# authors: Matteo Dell'Acqua and Leonardo Caproni
# date: 08/2022
#-------------------------------------------------------------------------------
# Reference: DOI
#-------------------------------------------------------------------------------
# Description: 	Genome-Wide Association Studies of non-collinear biolclimatic
#				variables and phenotypes.
#===============================================================================

# load libraries 
library(data.table)
library(tidyverse)
library(rMVP)

#GWAS
rm(list=ls())
options(stringsAsFactors = F)
wd<-"~/barley"
setwd(wd)

#get phenos, genos
phe <-"output/phenotypes.txt"
vcf <- "input/barley.snps.noWild.FINAL.AF005.vcf" 
vcf.pruned <- "input/barley.snps.noWild.FINAL.AF005.436.pruned_150_5_05.vcf" 

# make an output GWAS directory
GWAS <-"output/GWAS_MAF005"
dir.create(GWAS)

#2.1) Prepare MVP datasets for full SNPs data
MVP.Data(fileVCF=vcf, 
         filePhe=phe, 
         sep.phe="\t", 
         fileKin=F, 
         filePC=F,
         out="output/GWAS_MAF005/mvp.out")
genotype <- attach.big.matrix("output/GWAS_MAF005/mvp.out.geno.desc")
phenotype <- read.table("output/GWAS_MAF005/mvp.out.phe",head=TRUE)
map <- read.table("output/GWAS_MAF005/mvp.out.geno.map" , head = TRUE)
#K <- MVP.K.VanRaden(genotype)
#CV <- MVP.PCA(genotype, pcs.keep = 30)

map[1:5,1:5]

#2.1.1) Get K and PC for full SNPs dataset to add as covariates
#K <- MVP.K.VanRaden(genotype)
#CV <- MVP.PCA(genotype, pcs.keep=10, K=K)
#reg.covs <- list(K, CV)

#2.2) Prepare MVP datasets for pruned SNPs data. The "teff_pruned.hmp.txt"
MVP.Data(fileVCF=vcf.pruned, 
         filePhe=phe, sep.phe="\t", 
         fileKin=FALSE, 
         filePC=FALSE, 
         out="output/GWAS_MAF005/mvp.pruned")

genotype_pruned <- attach.big.matrix("output/GWAS_MAF005/mvp.pruned.geno.desc")
K.p <- MVP.K.VanRaden(genotype_pruned)
CV.p <- MVP.PCA(genotype_pruned, pcs.keep = 30)
#MVP.PCAplot(CV.p, plot3D = F)

#3) Run GWAS.
#3.1) With pruned SNPs dataset covariates
results <-"output/GWAS_MAF005/results"
dir.create(results)
setwd(results)

 #Define threshold
 n.markers <- 30515 # 47492 #n of markers for GWAS
 alpha <- 0.05 #Bonferroni alpha
 hap.blocks <- 2064 #n of non-associated markers r2=0.3 in the set
 thr <- (alpha/hap.blocks)*n.markers #0.512

#set a vector for errors
errs<-c()

#Run for loop MVP
for (j in c(2,3,4)){ 
  
  dirname<-paste0("PC.", j)
  dir.create(paste0("./",dirname))
  setwd(dirname)
  
  for(i in c(2:ncol(phenotype)-3)){    ##(2:5,50:52,68:70,116:118,161:165,197:201
    tryCatch({ #this function is used to handle errors (in any)
      tmpMVP <- 
        MVP(
          phe=phenotype[,c(1,i)], 
          geno=genotype,
          map=map,
          K=K.p,
          #CV.GLM=CV.p[,c(1:j)],
          CV.MLM=CV.p[,c(1:j)],
          CV.FarmCPU=CV.p[,c(1:j)],
          #nPC.FarmCPU=2, # "If pcs have been added in covariate files, PLEASE DO NOT assign value to nPC.GLM, nPC.MLM, nPC.FarmCPU"
          #nPC.GLM=2,
          #nPC.MLM=2,
          #perc=1,
          #priority="speed",
          ncpus=3,
          vc.method="EMMA",
          maxLoop=10,
          method.bin="FaST-LMM",
          permutation.threshold=F,
          file.output=T,
          #permutation.rep=500,
          threshold=0.05,
          method=c("FarmCPU","MLM"),
          file="jpg",
          dpi=300,
          col = c("grey30","grey60") ### M plot colors
        )
      save(tmpMVP, file=paste0(colnames(phenotype)[i],".rMVP.output.Rdata"))
      gc()
      
    }, error=function(e){errs<-c(errs, i)}) #trycatch ends #for i loop ends
  }#for i
  #go back one level
  setwd("..")
} # for jloop ends

setwd(wd)
