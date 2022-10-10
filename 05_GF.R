#===============================================================================
# author: Leonardo Caproni
# date: 08/2022
#-------------------------------------------------------------------------------
# Reference: DOI
#-------------------------------------------------------------------------------
# Description: 	Gradient Forest analysis and estimation of Genomic Vulnerability.
#				Adapted from Rhoné et al., 2020.
#				DOI: https://doi.org/10.1038/s41467-020-19066-4
#===============================================================================

#PRELIMINARY
options(stringsAsFactors = F)
rm(list=ls())

maindir <- "~/003_SSSUP/barley"
inputdir <- "~/003_SSSUP/barley/input"
setwd(maindir)

# load libraries
library(gradientForest)
library(rgdal)
library(raster)

#################################################################################################
load("GF_MAF005/metadata.beforeGF.new.data.Rdata")
load("GF_pruned/GF.function.Rdata")
load("modelEnv.ok.Rdata")
#################################################################################################
#load(file="metadata.afterGF.pruned.set.Rdata")

# create a folder to store results
GF.pruned <-"output/GF_pruned_REV1"
dir.create(GF.pruned)

# Look at the GF function
summary(gf.mem)

#arrange variables by importance
by.imp.mem <- names(importance(gf.mem))


##############################################################
# DEFINE CROPPING AREA
##############################################################
## project vulnerability only on the cropping area by AEZ
aez<- readOGR(dsn = "input/AEZ_32/Agro_ecology.shp")
summary(aez)

#only AEZ > 10 hits
table(out$AEZ31)>10

ETH_AEZs <- subset(aez, AEZ31 == "H3"|AEZ31 == "H4"|
                     AEZ31 == "M3"|AEZ31 == "M4"|
                     AEZ31 == "SH3"|AEZ31 == "SH4"|
                     AEZ31 == "SM3"|AEZ31 == "SM4")

# now crop and mask all model envs
clipped.modelEnv<-crop(modelEnv, extent(ETH_AEZs))
clipped.modelEnv<-mask(clipped.modelEnv, ETH_AEZs)
modelEnv<-clipped.modelEnv

clipped.modelEnv<-crop(modelFutureEnv_4550, extent(ETH_AEZs))
clipped.modelEnv<-mask(clipped.modelEnv, ETH_AEZs)
modelFutureEnv_4550<-clipped.modelEnv

clipped.modelEnv<-crop(modelFutureEnv_4570, extent(ETH_AEZs))
clipped.modelEnv<-mask(clipped.modelEnv, ETH_AEZs)
modelFutureEnv_4570<-clipped.modelEnv

clipped.modelEnv<-crop(modelFutureEnv_8550, extent(ETH_AEZs))
clipped.modelEnv<-mask(clipped.modelEnv, ETH_AEZs)
modelFutureEnv_8550<-clipped.modelEnv

clipped.modelEnv<-crop(modelFutureEnv_8570, extent(ETH_AEZs))
clipped.modelEnv<-mask(clipped.modelEnv, ETH_AEZs)
modelFutureEnv_8570<-clipped.modelEnv

##############################################################
# MAKE PREDICTIONS
##############################################################
# Predict allelic turnover actual environment
clim.land <- extract(modelEnv, 1:ncell(modelEnv), df = TRUE) #clim.layer.crop
clim.land <- na.omit(clim.land)
pred.mem <- predict(gf.mem, clim.land[,-1])

# For projections with CMIP5, rcp45 and year 2050
clim.land.4550 <- raster::extract(modelFutureEnv_4550, 1:ncell(modelFutureEnv_4550), df = TRUE) #clim.layer.crop
clim.land.4550 <- na.omit(clim.land.4550)
proj.mem.4550 <- predict(gf.mem, clim.land.4550[,-1]) 

# For projections with CMIP5, rcp45 and year 2070
clim.land.4570 <- extract(modelFutureEnv_4570, 1:ncell(modelFutureEnv_4570), df = TRUE) 
clim.land.4570 <- na.omit(clim.land.4570)
proj.mem.4570 <- predict(gf.mem, clim.land.4570[,-1]) 

# For projections with CMIP5, rcp85 and year 2050
clim.land.8550 <- extract(modelFutureEnv_8550, 1:ncell(modelFutureEnv_8550), df = TRUE) #clim.layer.crop
clim.land.8550 <- na.omit(clim.land.8550)
proj.mem.8550 <- predict(gf.mem, clim.land.8550[,-1]) 

# For projections with CMIP5, rcp85 and year 2070
clim.land.8570 <- extract(modelFutureEnv_8570, 1:ncell(modelFutureEnv_8570), df = TRUE) 
clim.land.8570 <- na.omit(clim.land.8570)
proj.mem.8570 <- predict(gf.mem, clim.land.8570[,-1]) 

#############################################################################################################
#Plot allelic turnover
#############################################################################################################
PCs <- prcomp(pred.mem, center=T, scale=F) #For pred.mem
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 +a2 -a1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
mask<-modelEnv$bio19 #Precipitation of Coldest Quarter
mask[]<-as.numeric(mask[]>0)
rastR <- rastG <- rastB <- mask
rastR[clim.land$ID] <- r
rastG[clim.land$ID] <- g
rastB[clim.land$ID] <- b
rgb.rast <- stack(rastR, rastG, rastB)

# #For projection rcp45 year 2050
PCs.proj.4550 <- prcomp(proj.mem.4550, center=T, scale.=F) #For proj.mem
a1.proj.4550 <- PCs.proj.4550$x[, 1]
a2.proj.4550 <- PCs.proj.4550$x[, 2]
a3.proj.4550 <- PCs.proj.4550$x[, 3]
r.proj.4550 <- a1.proj.4550 + a2.proj.4550
g.proj.4550 <- -a2.proj.4550
b.proj.4550 <- a3.proj.4550 +a2.proj.4550 -a1.proj.4550
r.proj.4550 <- (r.proj.4550 - min(r.proj.4550))/(max(r.proj.4550) - min(r.proj.4550)) * 255
g.proj.4550 <- (g.proj.4550 - min(g.proj.4550))/(max(g.proj.4550) - min(g.proj.4550)) * 255
b.proj.4550 <- (b.proj.4550 - min(b.proj.4550))/(max(b.proj.4550) - min(b.proj.4550)) * 255
mask.proj.4550<-modelFutureEnv_4550$bio19 #Precipitation of Coldest Quarter
mask.proj.4550[]<-as.numeric(mask.proj.4550[]>0)
rastR.proj.4550 <- rastG.proj.4550 <- rastB.proj.4550 <- mask.proj.4550
rastR.proj.4550[clim.land.4550$ID] <- r.proj.4550
rastG.proj.4550[clim.land.4550$ID] <- g.proj.4550
rastB.proj.4550[clim.land.4550$ID] <- b.proj.4550
rgb.rast.proj.4550 <- stack(rastR.proj.4550, rastG.proj.4550, rastB.proj.4550)

# #For projection rcp45 year 2070
PCs.proj.4570 <- prcomp(proj.mem.4570, center=T, scale.=F) #For proj.mem
r.proj.4570 <- PCs.proj.4570$x[, 1]
g.proj.4570 <- PCs.proj.4570$x[, 2]
b.proj.4570 <- PCs.proj.4570$x[, 3]
r.proj.4570 <- (r.proj.4570 - min(r.proj.4570))/(max(r.proj.4570) - min(r.proj.4570)) * 255
g.proj.4570 <- (g.proj.4570 - min(g.proj.4570))/(max(g.proj.4570) - min(g.proj.4570)) * 255
b.proj.4570 <- (b.proj.4570 - min(b.proj.4570))/(max(b.proj.4570) - min(b.proj.4570)) * 255
mask.proj.4570<-modelFutureEnv_4570$bio19 #Precipitation of Coldest Quarter
mask.proj.4570[]<-as.numeric(mask.proj.4570[]>0)
rastR.proj.4570 <- rastG.proj.4570 <- rastB.proj.4570 <- mask.proj.4570
rastR.proj.4570[clim.land.4570$ID] <- r.proj.4570
rastG.proj.4570[clim.land.4570$ID] <- g.proj.4570
rastB.proj.4570[clim.land.4570$ID] <- b.proj.4570
rgb.rast.proj.4570 <- stack(rastR.proj.4570, rastG.proj.4570, rastB.proj.4570)

# #For projection rcp85 year 2050
PCs.proj.8550 <- prcomp(proj.mem.8550, center=T, scale.=F) #For proj.mem
a1.proj.8550 <- PCs.proj.8550$x[, 1]
a2.proj.8550 <- PCs.proj.8550$x[, 2]
a3.proj.8550 <- PCs.proj.8550$x[, 3]
r.proj.8550 <- a1.proj.8550 + a2.proj.8550
g.proj.8550 <- -a2.proj.8550
b.proj.8550 <- a3.proj.8550 +a2.proj.8550 -a1.proj.8550
r.proj.8550 <- (r.proj.8550 - min(r.proj.8550))/(max(r.proj.8550) - min(r.proj.8550)) * 255
g.proj.8550 <- (g.proj.8550 - min(g.proj.8550))/(max(g.proj.8550) - min(g.proj.8550)) * 255
b.proj.8550 <- (b.proj.8550 - min(b.proj.8550))/(max(b.proj.8550) - min(b.proj.8550)) * 255
mask.proj.8550<-modelFutureEnv_8550$bio19 #Precipitation of Coldest Quarter
mask.proj.8550[]<-as.numeric(mask.proj.8550[]>0)
rastR.proj.8550 <- rastG.proj.8550 <- rastB.proj.8550 <- mask.proj.8550
rastR.proj.8550[clim.land.8550$ID] <- r.proj.8550
rastG.proj.8550[clim.land.8550$ID] <- g.proj.8550
rastB.proj.8550[clim.land.8550$ID] <- b.proj.8550
rgb.rast.proj.8550 <- stack(rastR.proj.8550, rastG.proj.8550, rastB.proj.8550)

# #For projection rcp85 year 2070
PCs.proj.8570 <- prcomp(proj.mem.8570, center=T, scale.=F) #For proj.mem
r.proj.8570 <- PCs.proj.8570$x[, 1]
g.proj.8570 <- PCs.proj.8570$x[, 2]
b.proj.8570 <- PCs.proj.8570$x[, 3]
r.proj.8570 <- (r.proj.8570 - min(r.proj.8570))/(max(r.proj.8570) - min(r.proj.8570)) * 255
g.proj.8570 <- (g.proj.8570 - min(g.proj.8570))/(max(g.proj.8570) - min(g.proj.8570)) * 255
b.proj.8570 <- (b.proj.8570 - min(b.proj.8570))/(max(b.proj.8570) - min(b.proj.8570)) * 255
mask.proj.8570<-modelFutureEnv_8570$bio19 #Precipitation of Coldest Quarter
mask.proj.8570[]<-as.numeric(mask.proj.8570[]>0)
rastR.proj.8570 <- rastG.proj.8570 <- rastB.proj.8570 <- mask.proj.8570
rastR.proj.8570[clim.land.8570$ID] <- r.proj.8570
rastG.proj.8570[clim.land.8570$ID] <- g.proj.8570
rastB.proj.8570[clim.land.8570$ID] <- b.proj.8570
rgb.rast.proj.8570 <- stack(rastR.proj.8570, rastG.proj.8570, rastB.proj.8570)

#############################################################################################################
# Estimate genomic offset, so-called "genomic vulnerability"
#############################################################################################################
temp.4550 <- vector("numeric", length = nrow(proj.mem.4550))
for (i in 1:ncol(proj.mem.4550)) {
  temp.4550 <- temp.4550 + (proj.mem.4550[,i]-pred.mem[,i])^2
}

temp.4570 <- vector("numeric", length = nrow(proj.mem.4570))
for (i in 1:ncol(proj.mem.4570)) {
  temp.4570 <- temp.4570 + (proj.mem.4570[,i]-pred.mem[,i])^2
}

temp.8550 <- vector("numeric", length = nrow(proj.mem.8550))
for (i in 1:ncol(proj.mem.8550)) {
  temp.8550 <- temp.8550 + (proj.mem.8550[,i]-pred.mem[,i])^2
}

temp.8570 <- vector("numeric", length = nrow(proj.mem.8570))
for (i in 1:ncol(proj.mem.8570)) {
  temp.8570 <- temp.8570 + (proj.mem.8570[,i]-pred.mem[,i])^2
}

##############################################################
GenVuln.4550 <- data.frame(sqrt(temp.4550))
GenVuln.4570 <- data.frame(sqrt(temp.4570))
GenVuln.8550 <- data.frame(sqrt(temp.8550))
GenVuln.8570 <- data.frame(sqrt(temp.8570))

GenVuln <- cbind(clim.land[,c(1)], GenVuln.4550, GenVuln.4570, GenVuln.8550, GenVuln.8570 )

colnames(GenVuln)[1] <- "cell_ID"
colnames(GenVuln)[2] <- "RCP 4.5 2050"
colnames(GenVuln)[3] <- "RCP 4.5 2070"
colnames(GenVuln)[4] <- "RCP 8.5 2050"
colnames(GenVuln)[5] <- "RCP 8.5 2070"
summary(GenVuln)

# assign coordinates to each of the pixels
clim.land2 <- extract(modelEnv, 1:ncell(modelEnv), df = TRUE)
clim.land2 <- cbind(coordinates(modelEnv), clim.land2)
clim.land2 <- na.omit(clim.land2)
genVuln <- cbind(clim.land2[,c(1,2)], GenVuln)

# Make GenVuln rasters
coordinates(genVuln) <- ~ x + y
gridded(genVuln) <- TRUE

genVuln.rast.4550 <- raster(genVuln, "RCP 4.5 2050")
names(genVuln.rast.4550)<-"Offset"

genVuln.rast.4570 <- raster(genVuln, "RCP 4.5 2070")
names(genVuln.rast.4570)<-"Offset"

genVuln.rast.8550 <- raster(genVuln, "RCP 8.5 2050")
names(genVuln.rast.8550)<-"Offset"

genVuln.rast.8570 <- raster(genVuln, "RCP 8.5 2070")
names(genVuln.rast.8570)<-"Offset"

#save output
save(genVuln, GenVuln, pred.mem, proj.mem.4550,
     proj.mem.4570, proj.mem.8550, proj.mem.8570,
     file="output/genVuln.Rdata")