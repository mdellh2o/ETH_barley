#this script performs penotypic analysis
#derives BLUPs and h2

options(stringsAsFactors = F)
library(corrplot)
library(ggplot2)
library(vioplot)
library(asreml)
library(reshape2)
library(ggcorrplot)
library(patchwork)
library(bestNormalize)
library(tidyr)
library(tidyverse)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/BARLEY/output"
setwd(wd)

## state functions to derive BLUPs, BLUEs, and variances
getBLUPs = function(m){
  res = data.frame(coef(m)$rand)
  res = res[grepl("ID",row.names(res)) & !grepl(":",row.names(res)) ,  , drop=FALSE]
  ##res = res[grepl("GY:ID_",row.names(res)), , drop=FALSE]
  return(res)
}

getBLUEs = function(m){
  res = data.frame(coef(m)$fix)
  res = res[grepl("ID",row.names(res)) & !grepl(":",row.names(res)) ,  , drop=FALSE]
  ##res = res[grepl("GY:ID_",row.names(res)), , drop=FALSE]
  return(res)
}

getVariance = function(asr, comp){
  var = summary(asr)$varcomp
  idx = which(rownames(var)==comp)
  v = var$component[idx]
  print(paste("variance component", v))
  return(v)
}


#######################
# get original datasets from basazen
########################

met<-read.delim("phenology_traits.txt")
head(met)

colnames(met)[1:10]<-c("LOCATION", "YEAR", "PLOTN", "REP", "BLK", "ROW", "COLUMN", "ACCESSION", "TRT", "ID")

metinfo<-met[,c("LOCATION","YEAR", "ID", "PLOTN", "ACCESSION", "TRT", "REP", "BLK", "ROW", "COLUMN")]
mettraits<-met[,c(11:ncol(met))]
met<-cbind(metinfo, mettraits)

names(met)
names(met)[1:10]<-c("LOCATION", "YEAR", "ID", "PLOT", "ACCESSION", "TRT", "REP", "BLK", "ROW", "COL")

# fix values
met$YEAR[met$YEAR == "1"]<-2016
met$YEAR[met$YEAR == "2"]<-2017

met$LOCATION[met$LOCATION == "ARSINEG"]<-"ArsiNegele"
met$LOCATION[met$LOCATION == "HOLLETA"]<-"Holeta"

#sort it
met<-met[order(met[, "LOCATION"],met[, "YEAR"], met[, "ID"], met[,"REP"]),]

#fix types of data
f = colnames(met)[1:10]
mtraits =  colnames(met)[!colnames(met) %in% f]
head(met[f])

met[,f] = lapply(met[,f], as.factor)
sapply(met[,f], levels)

head(met[,mtraits])
met[,mtraits] = lapply(met[,mtraits], as.numeric)
summary(met[,mtraits])

#save the raw data in tabular format
write.csv(met, file="barley.metric.traits.raw.txt", row.names=F, quote=F)

#save the raw data in tabular format
# write.csv(farm, file="barley.PVS.traits.raw.txt", row.names=F, quote=F)

##############################
#drop problematic phenotypes
met<-met[,-which(colnames(met) %in% c("NB_AUDP", "NB", "PM_AUDP", "PM", "SCALD_AUDP", "SCALD",
                                      "BCT", "HCT", "MCT"))]
names(met)

######################

###################
# check correlations
###################

## check data classes and assign factors
head(met)
f = colnames(met)[1:10]
traits =  colnames(met)[!colnames(met) %in% f]

sapply(met, class)

#remove missing IDs
dim(met)
if(length(which(is.na(met$ID)))>0){
  met<-met[!is.na(met$ID),] 
}
dim(met)

#plot traits
for (t in traits){ ## t = traits[1]
  print(t)
  tmp<-met[,c(f, t)]
  png(paste0("./traits/", t, ".vioplot.png"))
    vioplot(tmp[,ncol(tmp)] ~  tmp[,"REP"] + tmp[,"LOCATION"]+ tmp[,"YEAR"], col='lightgrey', 
          ylab=t, xlab="Location:Replica:year")
  dev.off()
}

#get updated list of traits
f = colnames(met)[1:10]
traits =  colnames(met)[!colnames(met) %in% f]

#ff = colnames(farm)[1:11]
#ftraits =  colnames(farm)[!colnames(farm) %in% ff]

#create subsets by locations
an<-subset(met, LOCATION=="ArsiNegele")
anmean<-aggregate(GY ~  ID, data=an, FUN = mean, na.rm = TRUE)

ho<-subset(met, LOCATION=="Holeta")
homean<-aggregate(GY ~  ID, data=ho, FUN = mean, na.rm = TRUE)

cor(anmean[,2], homean[,2])

an[traits] = lapply(an[traits], as.numeric)
ho[traits] = lapply(ho[traits], as.numeric)

###check correlations within location
anr1<-subset(an, REP==1)
anr2<-subset(an, REP==2)
cran<-cor(anr1[,traits], anr2[,traits], use="complete")

hor1<-subset(ho, REP==1)
hor2<-subset(ho, REP==2)
crho<-cor(hor1[,traits], hor2[,traits], use="complete")

anplot<-ggcorrplot(cran) + labs(title="ArsiNegele")
hoplot<-ggcorrplot(crho) + labs(title="Holetta")

crpl<-anplot | hoplot
crpl
ggsave(crpl, filename="correlation.replicas.png", width = 12, height = 6)


##############
## make BLUPs for metric data
#get updated list of traits
f = colnames(met)[1:10]
traits =  colnames(met)[!colnames(met) %in% f]

## get list of traits and number of trials
nLoc = 2
nYear = 2
nRep = 2

# set up dataframe for holding blups and h2
blups.met = data.frame(ID = unique(met$ID), row.names = unique(met$ID))
head(blups.met)
h2.met = data.frame()

# loop through traits with model for combined year, within year, within location, within location & year
for(t in traits){  ## t = traits[1]
  
  ## model for across year analysis  
  # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(LOCATION) + id(YEAR) + id(ID):id(LOCATION) + id(LOCATION):id(YEAR) + id(REP):id(LOCATION):id(YEAR), data=met, maxit=60, start.values=TRUE)
  # iv = asr$vparameters.table
  # iv$Constraint[]='U'
  # 
  asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
             random= ~ id(ID) 
             + id(LOCATION) 
             + id(YEAR) 
             + id(ID):id(LOCATION) 
             + id(ID):id(YEAR) 
             + id(LOCATION):id(YEAR) 
             + id(ID):id(LOCATION):id(YEAR) 
             + id(REP):id(LOCATION):id(YEAR)
             , data=met, maxit=100) ##, G.param = iv)
  
  print(summary(asr)$varcomp)
  
  
  b = predict(asr, "ID")
  b = b$pvals
  b.tmp = data.frame(b$ID, b$predicted.value)
  names(b.tmp) = c("ID", t)
  head(b.tmp)
  
  blups.met = merge(blups.met, b.tmp, by="ID")
  head(blups.met)
  
 # b = getBLUPs(asr)
#  b.tmp = data.frame(row.names(b), b$effect); names(b.tmp) = c("ID", t)
#  head(b.tmp)
#  b.tmp$ID<-sub("^ID_", "", b.tmp$ID)
  
 # blups.met = merge(blups.met, b.tmp, by="ID")
#  head(blups.met)
  
  Vg = getVariance(asr, 'ID')
  Vgl = getVariance(asr, 'ID:LOCATION')
  Vgy = getVariance(asr, 'ID:YEAR')
  Vgyl = getVariance(asr, 'ID:LOCATION:YEAR')
  Ve = getVariance(asr, 'units!R')
  
  h2 = Vg/(Vg + Vgl/nLoc + Vgy/nYear + Vgyl/(nLoc*nYear) + Ve/(nRep*nLoc*nYear)) ## h2 calculation for multiple locations in multiple years
  print(paste("H2 for", t, ":", h2))
  
  h2.met = rbind(h2.met, data.frame(trait=t, year="ALL", location="ALL", h2=h2))
  
  asr=NA
  
  ## models within year, across locations ##
  for(y in levels(met$YEAR)){  ## y=2016
    
    writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', y))
    
    # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(ID):id(LOCATION) + id(LOCATION) + id(REP):id(LOCATION), data=met, subset= YEAR==y, maxit=20, start.values=TRUE)
    # iv = asr$vparameters.table
    # iv$Constraint[]='U'
    
    asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
               random= ~ id(ID) 
               + id(LOCATION) 
               + id(ID):id(LOCATION) 
               + id(REP):id(LOCATION), 
               data=met, subset= YEAR==y, maxit=60) ##, G.param = iv)
    
    print(summary(asr)$varcomp)
    
    name = paste(t, y, sep=".")
    b = predict(asr, "ID")
    b = b$pvals
    b.tmp = data.frame(b$ID, b$predicted.value)
    names(b.tmp) = c("ID", name)
    head(b.tmp)
    
    blups.met = merge(blups.met, b.tmp, by="ID")
    head(blups.met)
    
    Vg = getVariance(asr, 'ID')
    Vgl = getVariance(asr, 'ID:LOCATION')
    Ve = getVariance(asr, 'units!R')
    
    h2 = Vg/(Vg + Vgl/nLoc + Ve/(nRep*nLoc))  ## calculation of H2, with one year at two or more locations
    print(paste("H2 for", t, "in year", y, ":", h2))
    
    h2.met = rbind(h2.met, data.frame(trait=t, year=y, location="ALL", h2=h2))
    
    asr=NA
  }   
  
  
  ## models for one location across years ##
  for(l in levels(met$LOCATION)){  ## l="ArsiNegele"
    
    writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', l))
    
    # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(ID):id(LOCATION) + id(LOCATION) + id(REP):id(LOCATION), data=met, subset= YEAR==y, maxit=20, start.values=TRUE)
    # iv = asr$vparameters.table
    # iv$Constraint[]='U'
    
    asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
               random= ~ id(ID) 
               + id(YEAR) 
               + id(ID):id(YEAR) 
               + id(REP):id(YEAR), 
               data=met, subset= LOCATION==l, maxit=60) ##, G.param = iv)
    
    print(summary(asr)$varcomp)
    
    name = paste(t, l, sep=".")
    b = predict(asr, "ID")
    b = b$pvals
    b.tmp = data.frame(b$ID, b$predicted.value)
    names(b.tmp) = c("ID", name)
    head(b.tmp)
    
    blups.met = merge(blups.met, b.tmp, by="ID")
    head(blups.met)
    
    Vg = getVariance(asr, 'ID')
    Vgy = getVariance(asr, 'ID:YEAR')
    Ve = getVariance(asr, 'units!R')
    
    h2 = Vg/(Vg + Vgy/nYear + Ve/nRep) ## H2 calculation for one location in two or more years 
    print(paste("H2 for", t, "in", l, ":", h2))
    
    h2.met = rbind(h2.met, data.frame(trait=t, year="ALL", location=l, h2=h2))
    
    asr=NA
  }   
  
  
  ## models for single location-year ##
  for(y in levels(met$YEAR)){  ## t='DB'; y=2016; 
    for(l in levels(met$LOCATION)){ #l='ArsiNegele'
      
      writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', y, 'in', l))
      
      asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
                 random= ~ id(ID) 
                 + id(REP), 
                 data=met, subset=c(YEAR==y & LOCATION==l), maxit=20)
      
      print(summary(asr)$varcomp)
      
      name = paste(t, y, l, sep=".")
      b = predict(asr, "ID")
      b = b$pvals
      b.tmp = data.frame(b$ID, b$predicted.value)
      names(b.tmp) = c("ID", name)
      head(b.tmp)
      
      blups.met = merge(blups.met, b.tmp, by="ID")
      head(blups.met)
      
      
      Vg = getVariance(asr, 'ID')
      Ve = getVariance(asr, 'units!R')
      
      h2 = Vg/(Vg + Ve/(nRep))  ## H2 calculation for one location in one year
      print(paste("H2 for", t, "in", l, "in year", y, ":", h2))
      
      h2.met = rbind(h2.met, data.frame(trait=t, year=y, location=l, h2=h2))
      
      asr=NA
      
    }
  }
}

h2.met[,4]<-round(h2.met[,4],3)
h2.met

dim(blups.met)
rownames(blups.met)<-blups.met$ID

#save model output
save(met, h2.met, blups.met, file="BLUPs.quantitative.traits.barley.predicted.Rdata")

#make plots
for (t in traits){ ## t = traits[1]
  tmp<-blups.met[,c(1, grep(t, colnames(blups.met)))]
  png(paste0("./traits/", t, ".BLUP.vioplot.png"))
  vioplot(tmp[,2:ncol(tmp)], col='lightgrey', 
          ylab=t, xlab="Location:Replica")
  dev.off()
}


#correlate among locations
anblup<-blups.met[,grep("ArsiNegele", colnames(blups.met))]
hoblup<-blups.met[,grep("Holeta", colnames(blups.met))]
cranho<-cor(anblup, hoblup, use="pairwise.complete")

anhoplot<-ggcorrplot(cranho, title = "ArsiNegele VS Holeta")
anhoplot
ggsave(anhoplot, filename="BLUPs.correlation.by.location.png", height=12, width=12)

#correlate among years
blup16<-blups.met[,grep("2016", colnames(blups.met))]
blup17<-blups.met[,grep("2017", colnames(blups.met))]
cr1617<-cor(blup16, blup17, use="pairwise.complete")

plot1617<-ggcorrplot(cr1617, title = "2016 Vs 2017")
plot1617
ggsave(plot1617, filename="BLUPs.correlation.by.year.png", height=12, width=12)


#check cors across phenotypes
combined<-blups.met[,-grep("ArsiNegele|Holeta|20", colnames(blups.met))]

p.mat <- cor_pmat(combined[,2:ncol(combined)], use="complete")
corplot<-ggcorrplot(cor(combined[,2:ncol(combined)]), p.mat=p.mat, type = "lower")
corplot
ggsave(corplot, filename="BLUPs.correlation.png", height=12, width=12)


###################
# farmer traits
###################

#farmer traits now
#Arsineghelle
an.farm<-read.delim("../input/Individual scoring of the farmers at Arsi Negelle.txt")
an.farm[,1]<-"ArsiNegele"
colnames(an.farm)
names(an.farm)<-toupper(names(an.farm))
head(an.farm)

###check correlations within location
#sort
an.farm<-an.farm[order(an.farm$REP, an.farm$DNA_CODE),]
anr1<-subset(an.farm, REP==1)
anr2<-subset(an.farm, REP==2)
cran<-cor(anr1[,10:ncol(anr1)], anr2[,10:ncol(anr2)], use="complete")
anplot<-ggcorrplot(cran) + labs(title="ArsiNegele")
anplot

#remove all computed averages
an.farm<-an.farm[,-grep("_AVR", names(an.farm))]
colnames(an.farm)[1:9]<-c("LOCATION", "PLOTN", "REP", "BLK", "ROW", "COLUMN", "ACCESSION", "TRT", "ID")

#make the dataframe in a long format
tmp<-an.farm
tmp.long<- melt(tmp, id.vars = colnames(tmp)[1:9])

tmp.long$F_TYPE<-tmp.long$variable
tmp.long$F_ID<-tmp.long$variable

tmp.long$variable<-sub("_.*_F.*$", "", tmp.long$variable)
tmp.long$F_TYPE<-sub("_F.*$", "", tmp.long$F_TYPE)
tmp.long$F_TYPE<-sub("^.*_", "", tmp.long$F_TYPE)
tmp.long$F_TYPE<-sub("F", "", tmp.long$F_TYPE)
tmp.long$F_ID<-sub("^.*_.*_", "", tmp.long$F_ID)

#fix order
tmp.long<-tmp.long[,c(1:9,12:13,10:11)]
head(tmp.long)

#make it into wide
tmp.wide<- tmp.long %>% pivot_wider(names_from = variable, values_from = value)
an.farm.wide<-tmp.wide

#Holetta
ho.farm<-read.delim("../input/Individual scoring of the farmers at Holetta.txt")
ho.farm[,1]<-"Holeta"
colnames(ho.farm)
names(ho.farm)<-toupper(names(ho.farm))
head(ho.farm)

###check correlations within location
#sort
ho.farm<-ho.farm[order(ho.farm$REP, ho.farm$DNACODE),]

hor1<-subset(ho.farm, REP==1)
hor2<-subset(ho.farm, REP==2)
crho<-cor(hor1[,10:ncol(hor1)], hor2[,10:ncol(hor2)], use="complete")
hoplot<-ggcorrplot(crho) + labs(title="Holeta")
hoplot

#remove all computed averages
ho.farm<-ho.farm[,-grep("_AVR", names(ho.farm))]
colnames(ho.farm)[1:9]<-c("LOCATION", "PLOTN", "REP", "BLK", "ROW", "COLUMN", "ACCESSION", "TRT", "ID")

#make the dataframe in a long format
tmp<-ho.farm
tmp.long<- melt(tmp, id.vars = colnames(tmp)[1:9])

tmp.long$F_TYPE<-tmp.long$variable
tmp.long$F_ID<-tmp.long$variable

tmp.long$variable<-sub("_.*_F.*$", "", tmp.long$variable)
tmp.long$F_TYPE<-sub("_F.*$", "", tmp.long$F_TYPE)
tmp.long$F_TYPE<-sub("^.*_", "", tmp.long$F_TYPE)
tmp.long$F_TYPE<-sub("F", "", tmp.long$F_TYPE)
tmp.long$F_ID<-sub("^.*_.*_", "", tmp.long$F_ID)

#fix order
tmp.long<-tmp.long[,c(1:9,12:13,10:11)]
head(tmp.long)

#make it into wide
tmp.wide<- tmp.long %>% pivot_wider(names_from = variable, values_from = value)
ho.farm.wide<-tmp.wide

#combine all data!!!
farm<-rbind(an.farm.wide, ho.farm.wide)
head(farm)

#fix the code for farmer appreciation values
names(farm)[which(names(farm) == "MAT")]<-"MA"
names(farm)[which(names(farm) == "GY")]<-"YA"
names(farm)[which(names(farm) == "BIO")]<-"BA"
names(farm)[which(names(farm) == "OVR")]<-"OA"

#fix types of data
f = colnames(farm)[1:11]
ftraits =  colnames(farm)[!colnames(farm) %in% f]
head(farm[f])

farm[,f] = lapply(farm[,f], as.factor)
sapply(farm[,f], levels)

head(farm[,ftraits])
farm[,ftraits] = lapply(farm[,ftraits], as.numeric)
summary(farm[,ftraits])

#save the raw data in tabular format
write.csv(farm, file="barley.farmer.traits.raw.txt", row.names=F, quote=F)

#collapse by farmers, getting an average per farmer group and check correlations
head(farm)
mean.farm<-farm %>% group_by(LOCATION, REP, ACCESSION, F_TYPE) %>% 
      summarise(MA= mean(MA))

an<-subset(mean.farm, LOCATION=="ArsiNegele")
anr1<-subset(an, REP==1)
anr2<-subset(an, REP==2)
stopifnot(all(anr1$ACCESSION == anr2$ACCESSION))
cor(anr1$MA, anr2$MA)

ho<-subset(mean.farm, LOCATION=="Holeta")
hor1<-subset(ho, REP==1)
hor2<-subset(ho, REP==2)
stopifnot(all(hor1$ACCESSION == hor2$ACCESSION))
cor(hor1$MA, hor2$MA)

ggcorrplot(cor(data.frame(hor1$MA, hor2$MA, anr1$MA, anr2$MA)))

#subset to IDs also present in metric data
#farm<-farm[farm[,"ID"] %in% rownames(blups.met),]
dim(farm)

## test of full model ## 
t = 'MA'
asr=asreml(fixed = as.formula(paste(t,'~ 1')),
           random= ~ id(ID) 
           + id(LOCATION) 
           # + id(ID):id(LOCATION)
           #+ id(BLK):id(REP):id(LOCATION)
           #+ id(F_TYPE)
           + id(F_TYPE) 
           + id(F_TYPE):id(LOCATION) , #one of the groups is not present in both localities
           data=farm, maxit=60)

print(summary(asr)$varcomp)

#now go on trait by trait
nLoc=2
nRep=2
#nBlk=25
#nGrp=3
#nTyp=3
nTyp=2

## set up dataframe for holding blups and h2
blups.farm = data.frame(ID = unique(farm$ID), row.names = unique(farm$ID))
head(blups.farm)
h2.farm = data.frame()

## loop through all traits
for(t in ftraits){  ## t = ftraits[3]

  # 
  asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
             random= ~ id(ID) 
             + id(ID)
             + id(F_TYPE)
             + id(ID):id(LOCATION)
             + id(ID):id(F_TYPE)
             + id(F_TYPE):id(LOCATION),
             data=farm, maxit=100) 
  
  print(summary(asr)$varcomp)

  b = predict(asr, "ID")
  b = b$pvals
  b.tmp = data.frame(b$ID, b$predicted.value)
  names(b.tmp) = c("ID", t)
  head(b.tmp)
  
  blups.farm = merge(blups.farm, b.tmp, by="ID")
  head(blups.farm)
  
  Vg = getVariance(asr, 'ID')
  Vgl = getVariance(asr, 'ID:LOCATION')
  Vt = getVariance(asr, 'F_TYPE')
  Vtl = getVariance(asr, 'F_TYPE:LOCATION')
  Ve = getVariance(asr, 'units!R')
  
  ## 
  h2 = Vg/(Vg + Vgl/nLoc + Vt/nTyp + Vtl/(nTyp*nLoc) + Ve/(nRep*nLoc*nTyp))  
  print(paste("H2 for", t, ":", h2))
  
  h2.farm = rbind(h2.farm, data.frame(trait=t, type='ALL', location="ALL", h2=h2))
  head(h2.farm)
  
  asr=NA
  
  ## BLUPs and H2 calculation for each location ##
  for(l in levels(farm$LOCATION)){  ## l='Holeta'
    
    writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', l))
    
    # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(REP) , data=farm2, subset= LOCATION==l, maxit=60) ##  start.values=TRUE)
    # iv = asr$vparameters.table
    # iv$Constraint[]='U'
    
    asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
               random= ~ id(ID) 
               + id(REP)
               + id(F_TYPE), 
               data=farm, subset= LOCATION==l, maxit=60) ##, G.param = iv)

    
    print(summary(asr)$varcomp)
    
    Vg = getVariance(asr, 'ID')
    Vr = getVariance(asr, 'REP')
    Vt = getVariance(asr, 'F_TYPE')
    Ve = getVariance(asr, 'units!R')
    
    h2 = Vg/(Vg + Vr/nRep + Vt/nTyp + Ve/(nRep*nTyp))  
    
    print(paste("H2 for", t, 'in', l, ":", h2))
    
    h2.farm = rbind(h2.farm, data.frame(trait=t, type='ALL', location=l, h2=h2))
    
    ## blups 
    name = paste(t, l, sep=".")
    b = predict(asr, "ID")
    b = b$pvals
    b.tmp = data.frame(b$ID, b$predicted.value)
    names(b.tmp) = c("ID", name)
    head(b.tmp)
    
    blups.farm = merge(blups.farm, b.tmp, by="ID")
    head(blups.farm)
    
    asr=NA
    
  }   
  
  ## BLUPs and H2 calculation for each type ##
  for(g in levels(farm$F_TYPE)){  ## g='M'
    
    writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', g))
    
    # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(REP) , data=farm2, subset= LOCATION==l, maxit=60) ##  start.values=TRUE)
    # iv = asr$vparameters.table
    # iv$Constraint[]='U'
    
    asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
               random= ~ id(ID) 
               + id(LOCATION)
               + id(ID):id(LOCATION),
               #+ id(REP):id(LOCATION),
               data=farm, subset= F_TYPE==g, maxit=60) ##G.param=iv)
    
    print(summary(asr)$varcomp)
    
    Vg = getVariance(asr, 'ID')
    Vgl = getVariance(asr, 'ID:LOCATION')
    Ve = getVariance(asr, 'units!R')
    
   # h2 = Vg/(Vg + Ve/(nLoc*nRep)) 
    h2 = Vg/(Vg + Vgl/nLoc + Ve/(nLoc*nRep)) 
   
    print(paste("H2 for", t, 'for', g, ":", h2))
    
    h2.farm = rbind(h2.farm, data.frame(trait=t, type=g, location='ALL', h2=h2))
    
    ## blups 
    name = paste(t, g, sep=".")
    b = predict(asr, "ID")
    b = b$pvals
    b.tmp = data.frame(b$ID, b$predicted.value)
    names(b.tmp) = c("ID", name)
    
    head(b.tmp)
    
    blups.farm = merge(blups.farm, b.tmp, by="ID")
    head(blups.farm)
    
    asr=NA
    
  }     
  
  
  ## BLUPs and H2 calculation for each gender at each location ##
  for(l in levels(farm$LOCATION)){  ## l='Holeta'; g='M'
    for(g in levels(farm$F_TYPE)){
      
      writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', g, 'in', l))
      
      # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(REP) , data=farm2, subset= LOCATION==l, maxit=60) ##  start.values=TRUE)
      # iv = asr$vparameters.table
      # iv$Constraint[]='U'
      
      asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
                 random= ~ id(ID) 
                 +id(REP),
                 data=farm, subset= F_TYPE==g & LOCATION==l, maxit=60) ##, G.param = iv)
      
      print(summary(asr)$varcomp)
      
      Vg = getVariance(asr, 'ID')
      Vr = getVariance(asr, 'REP')
      Ve = getVariance(asr, 'units!R')
      
      h2 = Vg/(Vg + Vr/(nRep) + Ve/(nRep))
      
      print(paste("H2 for", t, 'for', g, 'in', l,  ":", h2))
      
      h2.farm = rbind(h2.farm, data.frame(trait=t, type=g, location=l, h2=h2))
      
      ## blups 
      name = paste(t, l, g, sep=".")
      b = predict(asr, "ID")
      b = b$pvals
      b.tmp = data.frame(b$ID, b$predicted.value)
      names(b.tmp) = c("ID", name)
      
      head(b.tmp)
      
      blups.farm = merge(blups.farm, b.tmp, by="ID")
      head(blups.farm)
   
      asr=NA
    }
  }     
}  # end of BLUP calculations

head(blups.farm)
dim(blups.farm)
rownames(blups.farm)<-blups.farm[,1]

h2.farm[,4]<-round(h2.farm[,4],3)
h2.farm

#make plots
for (t in ftraits){ ## t = traits[1]
  tmp<-blups.farm[,c(1, grep(t, colnames(blups.farm)))]
  png(paste0("./traits/", t, ".BLUP.vioplot.png"))
  vioplot(tmp[,2:ncol(tmp)], col='lightgrey', ylab=t, las=2)
  dev.off()
}


#correlate among locations
anblup<-blups.farm[,grep("ArsiNegele", colnames(blups.farm))]
hoblup<-blups.farm[,grep("Holeta", colnames(blups.farm))]
cranho<-cor(anblup, hoblup, use="complete")

hoanplot<-ggcorrplot(cranho, title = "ArsiNegele VS Holeta")
hoanplot
ggsave(hoanplot, filename="farmer.BLUPs.correlation.by.location.png", height=12, width=12)

#correlate among  genders
mblup<-blups.farm[,grep("W$", colnames(blups.farm))]
wblup<-blups.farm[,grep("M$", colnames(blups.farm))]
crmw<-cor(mblup, wblup, use="complete")

mwplot<-ggcorrplot(crmw, title = "Men VS Women")
mwplot
ggsave(mwplot, filename="farmer.BLUPs.correlation.by.gender.png", height=12, width=12)


#mhoe a correlation plot on combined values
combined<-blups.farm[,-grep("Ars|Hol", colnames(blups.farm))]

p.mat <- cor_pmat(combined[,2:ncol(combined)], use="complete")
corplot<-ggcorrplot(cor(combined[,2:ncol(combined)]), p.mat=p.mat, type = "lower")
corplot

ggsave(corplot, filename="farmer.BLUPs.correlation.png", height=12, width=12)

#save output
save(farm, h2.farm, blups.farm, file="../output/pvs.BLUPs.barley.Rdata")

#####################################
#correlate metric traits with farmer traits

for (i in c("Hol", "Arsi")){
  combmet<-blups.met[,2:ncol(blups.met)]
  combmet<-combmet[,grep(i, colnames(combmet))]
  
  combfarm<-blups.farm[,2:ncol(blups.farm)]
  combfarm<-combfarm[,grep(i, colnames(combfarm))]
  #combfarm<-combfarm[,-grep("[MW]$", colnames(combfarm))]
  
  cr<-cor(combmet, combfarm, use="complete")
  
  combcorr<-ggcorrplot(cr)
  combcorr
  
  ggsave(combcorr, filename=paste0("farmer.metric.BLUPs.",i,".correlation.png"), height=12, width=12)
}#for i


#combine outputs in a trait file useful for GWAS
pheno<-merge(blups.met, blups.farm, by="ID")
rownames(pheno)<-pheno[,1]
write.table(pheno, file="../output/all.BLUPs.barley.txt", sep="\t", quote=F, row.names=F)

#write out h2
write.table(h2.farm, file="../output/h2.farm.txt", sep="\t", quote=F, row.names=F)
write.table(h2.met, file="../output/h2.met.txt", sep="\t", quote=F, row.names=F)

############################################
#
#           RESTART FROM HERE
#
##########################################

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/BARLEY/output"
setwd(wd)

library(ggplot2)
library(tidyr)
library(ggExtra)
library(patchwork)

load("pvs.BLUPs.barley.Rdata")
load("BLUPs.quantitative.traits.barley.predicted.Rdata")

head(blups.farm)
head(blups.met)

#make scatterplot of farmers scores
p1<-ggplot(blups.farm, aes(x=MA.ArsiNegele.M, y=MA.ArsiNegele.W)) + geom_point()+
  theme_bw() + xlim(1, 5)+ ylim(1, 5)  + geom_abline(intercept = 0, slope = 1) + geom_point(color='darkred')+ coord_fixed()
p1

p2<-ggplot(blups.farm, aes(x=MA.Holeta.M, y=MA.Holeta.W)) + geom_point()+
  theme_bw()+ xlim(1, 5)+ ylim(1, 5)  + geom_abline(intercept = 0, slope = 1)+ geom_point(color='navyblue')+ coord_fixed()
p2 

combplot<-p1|p2
combplot

ggsave(combplot, file="gender.correlation.plot.pdf", height=4, width=9)

#define ranking by farmers in both locations
perc<-.75
Anperc<-quantile(blups.farm[,"OA.ArsiNegele"], perc)
Hoperc<-quantile(blups.farm[,"OA.Holeta"], perc)

selAn<-blups.farm[which(blups.farm[,"OA.ArsiNegele"]>Anperc),"ID"]
selHo<-blups.farm[which(blups.farm[,"OA.Holeta"]>Hoperc),"ID"]

length(selAn)
length(selHo)

which(selAn %in% selHo)

#make a pca on phentoypes and plot farmer selection
tokeep<-names(blups.met)
tokeep<-tokeep[-grep("ID|An|Ho|Ord", tokeep)]
tokeep
traitpc<-prcomp(blups.met[,tokeep])

#add rank by farmers
evaltraits<-c("OA", "OA.ArsiNegele", "OA.Holeta", "OA.M", "OA.W")
toplot<-cbind(data.frame(traitpc$x[,1:5]), blups.farm[, evaltraits])
#in each evalutation column, just put the top percentile
perc<-.9
for (i in evaltraits){ # i = "OA"
  tmp<-quantile(toplot[,i], perc)
  top<-rownames(toplot)[which(toplot[,i] >= tmp)]
  toplot[which(!rownames(toplot) %in% top),i]<-NA
}
toplot

farmersel<-ggplot(toplot, aes(x=PC1, y=PC2)) + geom_point(alpha = 0.7, col="gray") + 
    theme_bw() +
    geom_point(data = subset(toplot, !is.na(OA.ArsiNegele)), colour = "red", alpha=0.2, cex=2) +
    geom_point(data = subset(toplot, !is.na(OA.Holeta)), colour = "blue", alpha=0.2, cex=2)
farmersel
ggsave(farmersel, file="PCA.farmers.selection.pdf", height=4, width=4)

#move from wide to tall format
gathercols <- c("OA.Holeta", "OA.ArsiNegele", "OA.Holeta.M", "OA.Holeta.W", "OA.ArsiNegele.M", "OA.ArsiNegele.W")

#from wide to tall
farmoa<-gather(blups.farm, key, score, gathercols, factor_key=TRUE)
farmoa<-farmoa[,c("ID", "key", "score")]
head(farmoa)  
  
#plot
pl1 <- ggplot(farmoa, aes(var1, var2)) + geom_point(alpha=0.3)
ggMarginal(pl1, type = "histogram")
  
  
  
