#===============================================================================
# author: Leonardo Caproni
# date: 08/2022
#-------------------------------------------------------------------------------
# Reference: DOI
#-------------------------------------------------------------------------------
# Description: 	Bioclimatic and phenological diversity analysis.
#===============================================================================

#PRELIMINARY
wd<-"~/barley/"
setwd(wd)
rm(list=ls())
options(stringsAsFactors = F)

# load libraries
library(vcfR)
library(data.table)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(adegenet)
library(car)
library(pegas)
library(ape)
library(seqinr)
library(ade4)
library(factoextra)
library(ggsci)
library(genetics)
library(EMMREML)
library(compiler)
library(hierfstat)
library(poppr)
library(hierfstat)

##########################################################################################
#DIVERSITY ANALYSIS PCA & DAPC
############################################################################################
# load vcf
myvcf <- read.vcfR ("input/barley.snps.noWild.FINAL.AF005.436.pruned_150_5_05.vcf", verbose = FALSE)
#myvcf <- read.vcfR ("input/barley.snps.noWild.FINAL.AF005.vcf", verbose = FALSE)
genind <- vcfR2genind(myvcf)

genind_scale <- scaleGen(genind, NA.method="mean") # Impute missing data using mean. This genotype dataset is used for PCA

# 3 PCA
## 3.1 Diversity analysis using PCA
geno.pca<- dudi.pca(genind_scale,cent=TRUE,scale=TRUE,scannf=FALSE, nf=3)
geno.pca_scores <- data.frame(geno.pca$li)# order scores
geno.pca_scores <- data.frame(scale(geno.pca_scores)) #scale scores
fviz_eig(geno.pca) #visualize scree plot
get_eigenvalue(geno.pca) #Visualize loadings statistics (variance explained etc...)

## DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS
set.seed(10101)

grp <- find.clusters(genind 
                     # n.pca = 435,
                     # max.n.clust=50
                     # n.clust = 11
                     ) #kept 435 PCs
dapc1 <- dapc(genind, grp$grp, n.pca = 3, 
              n.da = 5,
              ) #kept 3 discriminant functions

bic <- as.data.frame(grp$Kstat)

write.table(bic,file = "output/BIC.vs.nCLUST.11.txt", quote = F, sep = "\t", row.names = T)

# dapc1

scatter(dapc1, cell = 0, pch = 18:23, cstar = 2, mstree = FALSE, lwd = 2, lty = 2, posi.da="topright")

##Plot PCA by DAPC assig
assig <- as.data.frame(grp$grp)
geno <- rownames(assig)
rownames(assig) <- NULL
assig <- cbind(geno,assig)
colnames(assig)<-c("geno", "Cluster")

geno1 <- rownames(geno.pca_scores)
rownames(geno.pca_scores) <- NULL
geno.pca_scores<- cbind(geno1,geno.pca_scores)

geno.pca_scores[1:4,1:3]

geno.pca_scores_assig <- merge(x= assig, y= geno.pca_scores, by.x = "geno", by.y= "geno1" )
geno.pca_scores_assig[1:4,1:5]

## save outputs
save.image(file="diversity.analysis.step.1.11clust.Rdata")

# load AEZs of Ethiopia
aez<- readOGR(dsn = "input/AEZ_32/Agro_ecology.shp")

#revert to dataframe
passbio <- as.data.frame(passbio)

#get coordinates
coord <- passbio[, which(colnames(passbio) %in% c("ID","region","zone", "LAT","LON"))]

#get coord of genotyped accessions & omit NAs
coord <- coord[coord$ID %in% mIDs$IDs,]
coord <- na.omit(coord)

#now include DAPC assicgnation
coord <- merge(coord, assig, by.x="ID", by.y="geno")
coord <- merge(coord, passbio[,c(1,11,12)], by="ID")

#get intersection bw points and AEZs
pts<-st_as_sf(coord, coords = c("LON", "LAT"))
pts$LON<-coord$LON
pts$LAT<-coord$LAT

aez2 <- st_as_sf(aez)
st_crs(pts)<-st_crs(aez2)

inter<-st_intersection(aez2, pts)

#get relevant zones
dftmp<-data.frame(inter)
hitzones<-unique(dftmp[,"AEZ31"])

# add a factor for plotting
aez2$AEZ<-as.character(aez2$AEZ31)
aez2$AEZ[which(!aez2$AEZ %in% hitzones)]<-"ZNR"
aez2$AEZ<-as.factor(aez2$AEZ)

# keep agroecologies with at least 2 hits
table(dftmp[,"AEZ31"]) 
aez2$AEZ<-as.factor(aez2$AEZ)
table(aez2$AEZ)

#include a rougher AEZ definition
aez2$AEZ_type<-as.factor(sub("[0-9]$", "", aez2$AEZ))
aez2$AEZ_type

############################################################################################
#export metadata information
############################################################################################
metaout<-dftmp[,c("ID", "AEZ_type", "AEZ31", "region", "zone", "LON", "LAT")]

#check uniqueness before merging
unique(dftmp$ID[duplicated(dftmp$ID)])

out<-merge(metaout, geno.pca_scores_assig, by.x="ID", by.y="geno", all.y=T)
#out<- out[,c(1:6,8,7,9:ncol(out))]
unique(out$ID[duplicated(out$ID)])
#dlete dulplicated B_220 B_450
out <- out[-c(125,326),]

#now include all info
out <- merge(passbio [,c(1,11,12)], out,  by = "ID", all.x=T)

save(out, file="metadata.diversity.analysis.NEW.11clust.Rdata")
