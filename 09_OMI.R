#===============================================================================
# authors: 	Leonardo Caproni and Matteo Dell'Acqua
# date: 08/2022
#-------------------------------------------------------------------------------
# Reference: DOI
#-------------------------------------------------------------------------------
# Description: Optimal migration Inference; adapted from Rhoné et al. 2020
#				DOI: https://doi.org/10.1038/s41467-020-19066-4
#===============================================================================

#PRELIMINARY
rm(list=ls())
maindir <- "~/barley"
inputdir <- "~/barley/input"
setwd(maindir)
options(stringsAsFactors = F)

# load libraires
library(fields)
library(geosphere)
library(dbscan)
library(maps)
library(shape)
library(RColorBrewer)
library(raster)
library(rgdal)
library(tidyverse)

#load data of GF and GV
load(file="output/genVuln.Rdata")

# set colors parameters for figures
colramp = colorRampPalette(c("turquoise", "blue"))(30)

#get ethiopian extension
eth <- getData("GADM", country="ETH", level=0)
eth1 <- getData("GADM", country="ETH", level=1)

model.extent<-extent(eth)

# revert as data frame with coordinates
df.gv <- as.data.frame(genVuln, xy=T)

#format as the script requires
#PixelID <- rownames(df.gv)
genVuln<- df.gv[,c(6,7,1,2:5)]
colnames(genVuln)[1] <- "Longitude.DD"
colnames(genVuln)[2] <- "Latitude.DD"
colnames(genVuln)[3] <- "PixelID"

# Prepare all need objects for the loop and save them in different folders
coord <- cbind(genVuln$Longitude.DD, genVuln$Latitude.DD)
listModL <- names(genVuln)[4:ncol(genVuln)]

# Prepare data for the analysis
# 1 make foleders
dir.create("output/optimal_migration")
for (j in 1:length(listModL)){
  folder<-dir.create(paste0("output/optimal_migration/",listModL[j]))
}

# 2 make Rdata to run the loop
Predict.gf <- cbind(genVuln[,c(1,2,3)], pred.mem)

Proj.gf <- cbind(genVuln[,c(1,2,3)], proj.mem.4550)
save(Proj.gf, file=paste0("output/optimal_migration/",listModL[1],"/GFres.RData"))

Proj.gf <- cbind(genVuln[,c(1,2,3)], proj.mem.4570)
save(Proj.gf, file=paste0("output/optimal_migration/",listModL[2],"/GFres.RData"))

Proj.gf <- cbind(genVuln[,c(1,2,3)], proj.mem.8550)
save(Proj.gf, file=paste0("output/optimal_migration/",listModL[3],"/GFres.RData"))

Proj.gf <- cbind(genVuln[,c(1,2,3)], proj.mem.8570)
save(Proj.gf, file=paste0("output/optimal_migration/",listModL[4],"/GFres.RData"))

#clean memory
rm(Proj.gf, proj.mem.4550, proj.mem.4570,proj.mem.8550, proj.mem.8570, df.gv)

# Set list of non-collinear variables as from VIF
importantVAR<-c("bio2", "bio3", "bio4", "bio9", "bio12", "bio14", "bio15", "bio18",  "bio19")

#######################################################
# Infer and plot migration per cliamte model
#######################################################
## this may take a while

RES_Migr <- data.frame()

for (modL in listModL){ # modL <- "RCP.4.5.2050"
  
  datModL <- genVuln[, c("Longitude.DD", "Latitude.DD", "PixelID", modL)]
  
  # identify the 5% highest vulnerable pixels
  get5pc<-quantile(datModL[,4], .95)
  dat_maxVuln <- datModL[which(datModL[,4] > get5pc),]

  # Calc the distance between the 10% most vulnerable pixels for clustering
  dist_pts_maxVuln <- distm(dat_maxVuln[,c(1,2)], fun=distGeo)/1000
  
  # dbscan clustering
  db <- dbscan(dist_pts_maxVuln, eps = 500, minPts = 16) #Density-Based Spatial Clustering of Applications with Noise
  dat_maxVuln$groups <- db$cluster
  
  # Extract the higher value of GV by cluster
  Extract <- do.call(rbind, lapply(split(dat_maxVuln,dat_maxVuln$groups), 
                                   function(x) {return(x[which.max(x[,modL]),])}))
  Extract <- Extract[c(Extract$groups != 0),] #remove isolated pixels (unclustered pixels)
  selPixel <- Extract$PixelID
  
  # plot Edist to current conditions
  colnames(Extract)[4] <- "GenVulnMax"
  Extract$model <- modL
  
  load(paste0("output/optimal_migration/", modL, "/GFres.RData"))
  
  for (Pix in selPixel){ #Pix=26213 
    # Extract the projected genomic composition (PGC) at the selected pixel in the future climate 
    ProjmaxVuln.gf <- Proj.gf[Proj.gf$PixelID==Pix,]
    
    temp <- vector("numeric", length = nrow(Predict.gf))
    
    # Calculate the euclidian distance between the PGC at the selected pixel in the future climate  
    # and the PGC over all pixels in the current climate
    for (i in importantVAR) {
      temp <- temp + (ProjmaxVuln.gf[,i]-Predict.gf[,i])^2
    }
    
    Dist <- cbind(Predict.gf[,c(1,2,3)], sqrt(temp))
    colnames(Dist)[4] <- "EDist_to_maxVuln"
    
    # Extract the 5% smallest ED
    EDmin_3pc <- tail(sort(Dist$EDist_to_maxVuln, decreasing=TRUE), round(0.05*nrow(Dist),
                                                                          digits=0))
    dat_EDmin_3pc <- Dist[Dist$EDist_to_maxVuln %in% EDmin_5pc,]
    
    # Extract the 1% smallest ED
    EDmin_1pc <- tail(sort(Dist$EDist_to_maxVuln, decreasing=TRUE), round(0.01*nrow(Dist),
                                                                          digits=0))
    dat_EDmin_1pc <- Dist[Dist$EDist_to_maxVuln %in% EDmin_1pc,]
    
    
    # Compute the distance between the vulnerable pixel and the 3% closest ED
    Dist2Pix <- NULL
    for (i in 1:nrow(dat_EDmin_5pc)){
      Dist2Pix[i] <- distGeo(c(Dist[Dist$PixelID==Pix,]$Longitude.DD,
                               Dist[Dist$PixelID==Pix,]$Latitude.DD), 
                             c(dat_EDmin_3pc$Longitude.DD[i], dat_EDmin_3pc$Latitude.DD[i]))/1000
    }
    dat_EDmin_3pc <- cbind(dat_EDmin_3pc, Dist2Pix)
    
    # Compute the distance between the vulnerable pixel and the 1% closest ED
    Dist2Pix <- NULL
    for (i in 1:nrow(dat_EDmin_1pc)){
      Dist2Pix[i] <- distGeo(c(Dist[Dist$PixelID==Pix,]$Longitude.DD,
                               Dist[Dist$PixelID==Pix,]$Latitude.DD),
                             c(dat_EDmin_1pc$Longitude.DD[i], dat_EDmin_1pc$Latitude.DD[i]))/1000
    }
    dat_EDmin_1pc <- cbind(dat_EDmin_1pc, Dist2Pix)
    
    
    
    # Save the selected vulnerable pixels per group and ED min pixel
    ResUniPixel <- as.data.frame(c(Extract[Extract$PixelID==Pix,],
                                   Dist[which.min(Dist$EDist_to_maxVuln),],
                                   distGeo(c(Extract[Extract$PixelID==Pix,]$Longitude.DD, 
                                             Extract[Extract$PixelID==Pix,]$Latitude.DD),
                                           c(Dist[which.min(Dist$EDist_to_maxVuln),]$Longitude.DD,
                                             Dist[which.min(Dist$EDist_to_maxVuln),]$Latitude.DD))/1000,
                                   dat_EDmin_5pc[which.min(dat_EDmin_5pc$Dist2Pix),],
                                   dat_EDmin_1pc[which.min(dat_EDmin_1pc$Dist2Pix),]))
    
    colnames(ResUniPixel) <- c("lon.VulnPix", "lat.VulnPix", "VulnPixID", 
                               "GenVulnMax", "groups","model", 
                               "lon.MinVul", "lat.MinVul", "MinVulnPixID", "EDist_to_maxVuln", 
                               "geoDist2Vuln",
                               "lon.Clos5pc_MinVul", "lat.Clos3pc_MinVul", "Clos3pc_MinVulnPixID", 
                               "EDist_to_maxVuln_Closest_3pc", "geoDist2Vuln_Closest_3pc",
                               "lon.Clos1pc_MinVul", "lat.Clos1pc_MinVul", "Clos1pc_MinVulnPixID", 
                               "EDist_to_maxVuln_Closest_1pc", "geoDist2Vuln_Closest_1pc")
    
    RES_Migr <- rbind(RES_Migr, ResUniPixel)
  }
}

write.table(RES_Migr, "output/optimal_migration/clustering&Adaptation_results_allRCPs_2050_2070_500.txt", quote=FALSE, row.names = FALSE,
            sep="\t")

############### restart from here 
### Plot migration for all climates models

#set colors
colramp = colorRampPalette(c("turquoise", "blue"))(30)

#choose Km
Km=500
new.dir <- paste0("~/barley/output/optimal_migration/opt.migr.",Km)
dir.create(new.dir)
RES_Migr <- read.table(file = paste("output/optimal_migration/clustering&Adaptation_results_allRCPs_2050_2070_",Km,"Km.txt", sep=""), header=TRUE, sep="\t")

RES_Migr <- RES_Migr[, c("lon.VulnPix", "lat.VulnPix", "VulnPixID", "GenVulnMax", "groups","model", 
                         "lon.MinVul", "lat.MinVul", "MinVulnPixID", "EDist_to_maxVuln",
                         "geoDist2Vuln")]

#unique(RES_Migr$model)#check

for (modL in listModL) { # modL <- "RCP.4.5.2070"
  pdf(file = paste0(new.dir,"/", modL,".pdf"), width = 12, height = 3.5)
  
  par(mfrow=c(1,3), mgp=c(1.5, 0.5, 0), mar=c(3, 3, 1, 1), oma=c(0.5, 0.5, 2, 0.5))
  temp <- RES_Migr[RES_Migr$model == modL,]
  
  qtn<-quantile(temp$EDist_to_maxVuln)
  qtn<-round(qtn,5)
  
  temp$Class <- colramp[30]
  temp$Class[temp$EDist_to_maxVuln < qtn[4]] <- colramp[15]
  temp$Class[temp$EDist_to_maxVuln < qtn[2]] <- colramp[1]
  #### Migration plot
  library(shape) #(library to define arrows style)
  plot(temp[,c(1,2)], pch=".", xlim=c(model.extent@xmin, model.extent@xmax), 
       ylim=c(model.extent@ymin, model.extent@ymax), xlab="", ylab="", asp=1)
  plot(eth1, add=T)
  for (i in c(1:nrow(temp))){
    lines(c(temp$lon.MinVul[i], temp$lon.VulnPix[i]),
          c(temp$lat.MinVul[i], temp$lat.VulnPix[i]),
          col = temp$Class[i], lwd=1.5, lty=1)
    points(temp[,c(1,2)], pch=20, col="grey26")
  }
  
  # add legend
  legend("topright", legend=c(paste("ED <", qtn[2]), paste(qtn[2], "< ED <", qtn[4]), paste("ED >", qtn[4])),
         col=c(colramp[1], colramp[15], colramp[30]),
         lty=1, lwd=1.5, cex=0.8, y.intersp=0.8)
  
  ### Histogram migration load
  Class_fac <- factor(as.factor(temp$Class), levels = unique(temp$Class))
  plot(Class_fac, col=c(colramp[1], colramp[15], colramp[30]), border="white", 
       names.arg=c(paste("ED <", qtn[2]), paste(qtn[2], "< ED <", qtn[4]), paste("ED >", qtn[4])), # < 0.0023", "ED > 0.0023"), ##### Check
       cex.names=0.7, ylab="Counts")
  
  #### Histogram of geographic distances
  hist(temp$geoDist2Vuln, breaks=10, col="grey", border=FALSE, main="", xlab="Geographical distance of migration (km)")
  mtext(modL, outer = TRUE,  cex=1.2)

  dev.off()
  
}

