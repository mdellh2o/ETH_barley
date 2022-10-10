#===============================================================================
# authors: Leonardo Caproni and Matteo Dell'Acqua
# date: 08/2022
#-------------------------------------------------------------------------------
# Reference: DOI
#-------------------------------------------------------------------------------
# Description: 	Genetic PCA and Discriminant anlysis of Principal Components 
#				(DAPC). It also derives AEZ for each sampling point.
#===============================================================================

# Preliminary
rm(list=ls())
inputdir<-"~/003_SSSUP/barley/input"
setwd(inputdir)

library(raster)
library(maptools)
library(rgdal)
library(data.table)

# Read the .csv file containing all info about teff accessions used in the experiment
file<-"passport.txt"
passport <- fread(file)
#str(passport)
names(passport)
summary(passport)

# Import row.type information
row.type  <-fread("row.type.txt")
names(row.type)
row.type <- row.type[,c(3,6)]
colnames(row.type)[2] <- "row.type"

row.type$row.type <- as.factor(row.type$row.type)

passport <- merge(passport, row.type, by.x= "sample name", by.y="DNA_Code", all.x=T)

summary(passport)

# check uniqueness of entries
unique(passport$`sample name`[duplicated(passport$`sample name`)])

#extract mapped features and put them in a SPDF
pasmap<-passport[which(!is.na(passport[,"lat_dec"])),]
pass<-SpatialPointsDataFrame(data=pasmap, coords =pasmap[,c("lon_dec", "lat_dec")],
                             proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))  

plot(pass, main="Map of Teff")

#extract extension and attributes
minmax<-extent(pass)
att<-pass@data 

#get elevation
alt<-getData('alt', country='ETH')

#get additional data
et1<-getData('GADM' , country="ETH", level=1)

#get bioclim at max resolution using tef extension to get the right tile
#bcdata<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/TEFF.distribution/analysis/0.load.shp.add.layer"
#bioc <- getData("worldclim",var="bio",lon = minmax[1], lat=minmax[3],res=0.5)

# load historical data ERA-Interim
datafiles <- Sys.glob("historical/*.tif") #Or whatever identifies your files
bioc <- stack(datafiles)

#make a sample plot
plot(bioc[[1]])
plot(et1, add=T)
plot(pass, add=T)

#extract values
altitude<-raster::extract(alt,pass)
values <- raster::extract(bioc,pass)
biovar <- cbind.data.frame(ID=att[,1],altitude, values)
head(biovar)


#bring together with passport
names(biovar)
#edit names biovar
colnames(biovar) <- sub("CHELSA_","", colnames(biovar))
colnames(biovar) <- sub("10_","", colnames(biovar))
summary(biovar)

names(passport)
passport <- passport[,-20] #drop old altitude

passbio<-merge(passport, biovar, by.x = "sample name", by.y = "ID", all.x=T)


#check number of accessions mapped
head(passbio)
length(which(!is.na(passbio$lat_dec)))
length(which(!is.na(passbio$lon_dec)))

#rearrange the dataset
names(passbio)
passbio<-passbio[,c(1:3,9,12:15,18,19,22:25,36:43,26:35)]
unique(passbio$Country)
passbio <- filter(passbio, (Country %in% c("Ethiopia", NA)))

# fix colnames
colnames(passbio) [1] <- "ID"
colnames(passbio) [9] <- "LAT"
colnames(passbio) [10] <- "LON"

unique(passbio$ID[duplicated(passbio$ID)])

summary(passbio)

#save object for 
save(passbio, file="../output/passport.bioclim.data.OK.Rdata")
write.csv(passbio, file="../output/passport.bioclim.data.OK.csv", row.names=F, quote=F)

