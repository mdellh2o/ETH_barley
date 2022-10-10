#===============================================================================
# authors: Matteo Dell'Acqua and Leonardo Caproni
# date: 08/2022
#-------------------------------------------------------------------------------
# Reference: DOI
#-------------------------------------------------------------------------------
# Description: 	Gradient Forest analysis and estimation of Genomic Vulnerability.
#				Adapted from Rhoné et al., 2020.
#				DOI: https://doi.org/10.1038/s41467-020-19066-4
#===============================================================================

rm(list=ls())
options(stringsAsFactors = F)

#PRELIMINARY
maindir <- "~/barley/"
setwd(maindir)
#load libraries
library(ggplot2)
library(ggcorrplot)
library(ggfortify)
library(patchwork)
library(tidyverse)
library(ggsci)
library(gghighlight)
library(hrbrthemes)
library(hrbrthemes)

#define DAPC palette 

pal.cl <- c("#d2982e", "#7065bb", "#a9b23f", "#b54f90", "#7db754",
            "#ba4a4e", "#43c9b0", "#bb5c2d", "#59bc77", "#b79046",
            "#537326")

# load data
load(file="metadata.diversity.analysis.NEW.11clust.Rdata")
load(file="output/passport.bioclim.data.Rdata")

#only the ones with coordinates
dftmp <- out[,c("ID", "LAT", "LON","Cluster")]
dftmp<-na.omit(dftmp)
dftmp1 <- as.data.frame(subset(passbio, passbio$ID %in% dftmp$ID))

biovar<-dftmp1[,c(1, grep("^bio", colnames(dftmp1)))]
biovar<-na.omit(biovar) #drop NAs if any

# keep only non-colinear variables
tokeep <-c ( "bio2", "bio3", "bio4", "bio9", "bio12","bio14", "bio15",  "bio18", "bio19")
biovar1 <- biovar[,c("ID",tokeep)]

#Scale bioclimatic variables
for(i in c(2,5)){
  biovar1[,i] <- (biovar1[,i]/10)
}

for(j in c(3,4)){
  biovar1[,j] <- (biovar1[,j]/100)
}


# Pricipal component analysis of Bioclimatic variables of Ethiopian landraces
pcbio<-prcomp(biovar1[,-1], scale. =T, center=T)
str(pcbio)

pcbioplot <- autoplot(pcbio, 
                      scale. = TRUE, 
                      data = dftmp, 
                      colour = 'Cluster', 
                      size = 1.2,
                      loadings = TRUE,
                      loadings.colour = 'black', 
                      loadings.label = TRUE, 
                      loadings.label.size = 4,
                      alpha = 0.9,
                      loadings.label.repel=T
                      # frame = TRUE,
                      # frame.type = 'norm',
                      # frame.colour = 'Cluster'
                      ) +
                #ggsci::scale_fill_d3() +
                # ggsci::scale_color_d3() +
                scale_color_manual(values=pal.cl) +
                geom_vline(linetype = "dashed", xintercept = 0, color="gray10") +
                geom_hline(linetype = "dashed", yintercept = 0, color="gray10") +
                ggplot2::theme_light() +
                theme(legend.position="none") +
                labs(title = "Bioclimatic Diversity" )

#correlations
pcout<-pcbio$x[,1:3]
pcout<-data.frame(ID=biovar[,1], pcout)
colnames(pcout)<-paste0(colnames(pcout), "_bio")

#bring it back to passbio
biopc<-merge(biovar1, pcout, by.x="ID", by.y="ID_bio", all.x=T)
rawbio<-biopc[,grep("^bio", colnames(biopc))]
rawpc<-biopc[,grep("_bio$", colnames(biopc))]
colnames(rawpc) <- sub("_bio", "", colnames(rawpc))

#plot COR
corrbioPC<-ggcorrplot(cor(rawpc,rawbio, use="pairwise.complete"), 
                      #method = "circle",
                      outline.col = "gray",
                      ggtheme = ggplot2::theme_light) 
                      #ggsci::scale_color_nejm()
  
load(file="BLUPs.quantitative.traits.barley.predicted.Rdata")
############################################################

raw.phenology<- blups.met[,-grep("ArsiNegele|Holeta|2017|2016", colnames(blups.met))]
raw.phenology <- raw.phenology [,c(1,4:6)]

write.table(raw.phenology, file="publication/supplementary tables/supplementary.table.X.txt", sep = '\t', row.names = F, quote = F)
pheno.clust <- merge(raw.phenology, out, by.x="ID", all.x =T)

# PCA phenology, highlight DAPC clusters
pcpheno<-prcomp(raw.phenology[,-1], scale. =T, center=T)

pcphenoplot <- autoplot(pcpheno, 
                      scale. = TRUE, 
                      data = pheno.clust, 
                      colour = 'Cluster',
                      shape = 'type',
                      size = 'type',
                      loadings = TRUE,
                      loadings.colour = 'black', 
                      loadings.label = TRUE, 
                      loadings.label.size = 4,
                      alpha = 0.9,
                      loadings.label.vjust = 0.5,
                      loadings.label.hjust = -0.5) +
              scale_color_manual(values=pal.cl, na.translate=FALSE) +
              geom_vline(linetype = "dashed", xintercept = 0, color="gray10") +
              geom_hline(linetype = "dashed", yintercept = 0, color="gray10") +
              scale_shape_manual(values = c(17,16,15), na.translate = F) +
              scale_size_manual(values=c(1.5,1.3,1.3), na.translate = F) +
              ggplot2::theme_light() +
              guides(colour=guide_legend(title="DAPC \ncluster", order = 1),
                     shape=guide_legend(title="Status", order = 2),
                     size="none") +
              labs(title = "Phenological Diversity")
              
pcout.pheno<-pcpheno$x[,1:3]
pcout.pheno<-data.frame(ID=raw.phenology[,1], pcout.pheno)
colnames(pcout.pheno)<-paste0(colnames(pcout.pheno), "_pheno")


###add first 3 PCs for ph (only combined measures, separated by metric traits and farmer traitas)
comb<-merge(out, biopc, by.x="ID", by.y="ID", all.x = TRUE)
comb<-merge(comb, pcout.pheno, by.x="ID", by.y="ID_pheno", all.x = TRUE)

# merge all passport info, pheno and bio in a big dataset
pheno<-merge(comb, raw.phenology, by="ID", all.x=T)
names(pheno)

# make sure they are factors
f<-c("AEZ_type","AEZ31", "region", "zone", "type","row.type","Cluster")

pheno[f]<-lapply(pheno[f], as.factor)

#get out names by clusters
byclust<-split(pheno, pheno[,"Cluster"])
idclust<-lapply(byclust, function(x) x[,1])  


#########################################################################
#Non-parametric Kruskas-Wallis Test
#########################################################################

names(pheno)
formulae2 <- lapply(colnames(pheno)[11:ncol(pheno)], function(x) as.formula(paste0(x, " ~ Cluster")))

res.kw <- lapply(formulae2, function(x) kruskal.test(x, data = pheno))
names(res.kw) <- format(formulae2)

pvals.kw <- lapply(formulae2, function(x) kruskal.test(x, data = pheno)[["p.value"]])
names(pvals.kw) <- format(formulae2)


###Bonferroni correction
alpha=0.05

thr<-alpha/length(pvals.kw)

hits.kw<-which(pvals.kw<thr)
length(pvals.kw[hits.kw])
sig_kw <- t(as.data.frame(pvals.kw[hits.kw]))

############################################################################################
#export data 
############################################################################################
names(pheno)
pheno <- pheno[,c(1:13,23:28,14:22,29:31)]

gwas.ph <- pheno [,c(1,20:31)]
row.t.df <- pheno[,c(1,3)]
row.t.df$row.type <- as.character(row.t.df$row.type)


row.t.df <- row.t.df %>%
  dplyr::mutate(Two.rowed = case_when(startsWith(row.type, "Two") ~ "1",
                                     !startsWith(row.type, "Two") ~ "0"))
row.t.df <- row.t.df %>%
  dplyr::mutate(Six.rowed = case_when(startsWith(row.type, "Six") ~ "1",
                                      !startsWith(row.type, "Six") ~ "0"))
row.t.df <- row.t.df %>%
  dplyr::mutate(Irregular = case_when(startsWith(row.type, "Irr") ~ "1",
                                      !startsWith(row.type, "Irr") ~ "0"))

row.t.df <- row.t.df [,c(1,3:5)]

gwas.ph <- merge(gwas.ph, row.t.df, by="ID")

# write table for GWAS
write.table(gwas.ph, file = "output/phenotypes.txt", quote = FALSE, sep="\t", row.names = FALSE)