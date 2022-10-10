#===============================================================================
# author: Leonardo Caproni
# date: 08/2022
#-------------------------------------------------------------------------------
# Reference: DOI
#-------------------------------------------------------------------------------
# Description: 	partial Redudancy Analysis (pRDA). Variance partitionig and 
#				outlier detection. Adapted from Capblancq and Forester 2021
#				DOI: https://doi.org/10.1111/2041-210X.13722
#===============================================================================

#PRELIMINARY
wd<-"~/barley/"
setwd(wd)
rm(list=ls())
options(stringsAsFactors = F)

# load libraries
library(pegas)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(cowplot)
library(corrplot)
library(rgeos)
library(stringr)
library(CMplot)

#load environmental data
load(file = "bioclim.barley.383.Rdata") # load env
passport <- read.delim(file = "output/out.RES.11clust.txt")
hulls <- read.delim(file="output/barley.grain.hulls.txt")
pheno <- read.delim(file = "output/phenotypes.OK.txt", header = T)
#gen <- read.table("input/snp.barley.forGF", header = T, row.names = 1)
gen <- read.table("input/snp.barley.pruned.forGF", header = T, row.names = 1)


# start with geno data, check missingness
perc.miss <- sum(is.na(gen))/(nrow(gen)*ncol(gen)) * 100
perc.miss

# simple imputation
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
gen.imp <- as.data.frame(gen.imp)

#if you want to add row type
passport <- merge(passport, hulls, by="ID")
env<- merge(env, passport [, c(1:3,14,4:7,10:13)], by="ID")
#env <-na.omit(env)

colnames(env) [32:34] <- c("PC1","PC2","PC3")

# need 239 geno and env data, same order
ord<-env[,1]
gen.imp<-subset(gen.imp, rownames(gen.imp) %in% ord)
gen.imp<-gen.imp[order(match(rownames(gen.imp), ord)), , drop = FALSE]

# confirm that genotypes and environmental data are in the same order
identical(rownames(gen.imp), env[,1]) 

# extract geo features
coord <- env[,c("LON","LAT")]
coord <-na.omit(coord)
pcnm <- pcnm(dist(coord))  #this generates the PCNMs, you could stop here if you want all of them
keep <- ceiling(length(which(pcnm$value > 0))/2)

# subset just non-colinear environmental predictors + altitude
non.collinear <-c ("bio4", "bio2", "bio15", "bio3","bio12", "bio19", "bio18", "bio14", "bio9") # as from VIF

# prepare data for RDA
pred <- env[, c(non.collinear)]
pred.alt <-  env[, c(non.collinear,"altitude")]
row.type <- as.factor(env[,c("row.type")])
grain.hull <- as.factor(env[,c("grain.hulls")])
regions <- as.factor(env[,"region"])
aezs <- as.factor(env[,"AEZ31"])
dapc.clust <- as.factor(env[,"Cluster"])
geo <- scores(pcnm)[,1:keep]
geo.red <- scores(pcnm)[,1:10]
pred.pcagen <- env[, c(non.collinear, "PC1","PC2","PC3")]

row.hull <-cbind(row.type, grain.hull)
geo.regions <-  cbind(regions, aezs,geo.red)

# extract phenology data
phenotypes <- subset(pheno, pheno$ID %in% ord)
phenotypes <- phenotypes[, 11:13]

# also extract genetic PCs
gen.pcs <- env[,32:34]

# save and load data
rm(gen)
save.image(file="output/RDA.str.metadata.1.Rdata")

###############################################
# Redundancy Analysis 
###############################################

#loda metadata
load(file="output/RDA.str.metadata.1.Rdata")
# should try with pruned

names(gen.imp) <- gsub("X", "", names(gen.imp))
names(gen.imp) <- gsub("H", "H_", names(gen.imp))

## Pure neutral population structure model  
RDA_env <- rda(gen.imp ~ bio4 + bio2 + bio15 + bio3 + bio12 + bio19 + bio18 + bio14 + bio9 + Condition(PC1 + PC2 + PC3),  env)
RDA_env

RsquareAdj(RDA_env)
screeplot(RDA_env, main="Eigenvalues of constrained axes")

## load Function rdadapt
#### Function to conduct a RDA based genome scan
rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

# use RDA adapt
rdadapt_env<-rdadapt(RDA_env, 2)

## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_env$p.values)
-log10(thres_env)

outliers <- data.frame(Loci = colnames(gen.imp)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], chr = unlist(lapply(strsplit(colnames(gen.imp)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))
outliers <- outliers[order(outliers$chr, outliers$p.value),]

## List of outlier names
outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$chr)])
Outliers <- rep("Neutral", length(colnames(gen.imp)))
Outliers[colnames(gen.imp)%in%outliers$Loci] <- "Outliers"
TAB_manhatan <- data.frame(pos = 1:length(colnames(gen.imp)), 
                           pvalues = rdadapt_env$p.values, 
                           Outliers = Outliers)

# make a Manhattan plot
name.loc<- colnames(gen.imp)
m.plot <- cbind(name.loc, TAB_manhatan)


m.plot[c('CHR', 'POS')] <- str_split_fixed(m.plot$name.loc, '_', 2)
names(m.plot)
m.plot <- m.plot[,c(1,5,6,3)]
  
colnames(m.plot)[1]<-"SNP"
colnames(m.plot)[4]<-"RDA_env"  

CMplot(m.plot, 
       plot.type='m', band=1,
       cex = 0.6,
       col=c("grey30","grey60"),
       threshold= thres_env, 
       threshold.col=c("darkgreen","darkred"),
       signal.line=1, 
       signal.col="red", 
       amplify=T, signal.cex = 0.8, signal.pch = 19, 
       ylim=NULL, LOG10 = TRUE,
       width=12,height=4,
       file = "pdf",
       memo="full")

outliers[c('CHR', 'POS')] <- str_split_fixed(outliers$Loci, '_', 2)
colnames(outliers) [1]<- "SNP"
outliers <- outliers[,c(1,4,5,2)]

write.table(outliers, file="output/RDA.outliers.pruned.txt", sep = '\t', quote = F, row.names = F)



## Null model
RDA0 <- rda(gen.imp ~ 1,  env[,c(2:22,32:34)]) 

#RDA full model
RDAfull <- rda(gen.imp ~ PC1 + PC2 + PC3 + LON + LAT + bio1 + 
                 bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + 
                 bio8 + bio9 + bio10 + bio11 + bio12 + bio13 +
                 bio14 + bio15 + bio16 + bio17 + bio19,  env)

mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)

# Variance Partitioning: Partial RDA

## Full model
pRDAfull <- rda(gen.imp ~ PC1 + PC2 + PC3 + LON + LAT + bio4 + bio2 + bio15 + bio3 + bio12 + bio19 + bio18 + bio14 + bio9,  env)
pRDAfull
RsquareAdj(pRDAfull)
aov.full <- anova(pRDAfull)

## Pure climate model
pRDAclim <- rda(gen.imp ~ bio4 + bio2 + bio15 + bio3 + bio12 + bio19 + bio18 + bio14 + bio9 +  Condition(LON + LAT + PC1 + PC2 + PC3),  env)
pRDAclim
RsquareAdj(pRDAclim)
aov.clim <- anova(pRDAclim)

##Pure geography model
pRDAgeog <- rda(gen.imp ~ LON + LAT + Condition(bio4 + bio2 + bio15 + bio3 + bio12 + bio19 + bio18 + bio14 + bio9  +  PC1 + PC2 + PC3),  env)
pRDAgeog
RsquareAdj(pRDAgeog)
aov.geog <- anova(pRDAgeog)

## Pure neutral population structure model  
pRDAstruct <- rda(gen.imp ~ PC1 + PC2 + PC3 + Condition(LON + LAT + bio4 + bio2 + bio15 + bio3 + bio12 + bio19 + bio18 + bio14 + bio9),  env)
RsquareAdj(pRDAstruct)
aov.struct <- anova(pRDAstruct)

## Pure geno
pRDAgeno <- rda(gen.imp ~ PC1 + PC2 + PC3 ,  env)
pRDAgeno
RsquareAdj(pRDAgeno)
aov.geno <- anova(pRDAgeno)

## Pure neutral population structure model  
RDA_env <- rda(gen.imp ~ bio4 + bio2 + bio15 + bio3 + bio12 + bio19 + bio18 + bio14 + bio9 + Condition(PC1 + PC2 + PC3),  env)
RDA_env
RsquareAdj(RDA_env)
aov.env <- anova(RDA_env)


save.image(file="output/RDA.results.Rdata")