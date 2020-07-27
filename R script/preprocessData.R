rm(list=ls())
gc()

### Load Required Packages  ####
library(compositions)
library(data.table)
library(reshape2)
library(doParallel)
library(igraph)
library(caret)
library(tidyverse)
library(cvAUC)
library(glmnet)
library(xgboost)
library(vegan)
library(knitr)


# Read csv
micro = read_csv("Data/mouseMicrobiome_noAduj.csv")
#get otu names
cnames = data.frame(col = colnames(micro))
cnames = cnames %>% 
  filter(str_detect(string = col,pattern = "k__"))
cnames = separate(cnames,1,sep = ";",remove = F,
                  into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
#Rename Feture for ease of use
cnames$ID = paste("V",1:nrow(cnames),sep = "")
featurePosition = which(str_detect(colnames(micro),pattern = "k__")==T)

#OTU data
otu = data.frame(micro[,featurePosition])
colnames(otu) = cnames$ID
#get metdata
metadata = data.frame(micro[,-featurePosition])
#add cage_mouse_strain identifier
metadata$cage_mouse = paste(metadata$Cage,metadata$Mouse,metadata$Strain,sep = "_")
metadata$mouse_cage = paste(metadata$Mouse,metadata$Cage,sep = "_")
rownames(otu) = micro$index

#examine otu table
rs = rowSums(otu)
hist(rs,breaks = 25)
summary(rs)

#add total reads to metadata
metadata$totalReads = rs

#remove samples with too few reads
minReads = 10000 #min total reads per sample
sum(rs<minReads)
rownames(otu)[rs<minReads]
removeSamples = metadata[rs<minReads,] #samples where total reads are below min. number of reads
#outlier sampleID - cage_mouse
outlierSamples = removeSamples$cage_mouse

#Filter Metadata - remove outlier samples
keepMetaData = metadata %>% 
  filter(!cage_mouse %in% outlierSamples )
keepOTU = otu[rownames(otu)%in%keepMetaData$index,]
keepOTU = data.frame(index = rownames(keepOTU),keepOTU)

#join Metadata and OTU
rawData = left_join(keepMetaData,keepOTU)
#remove due to missing paired sample
pairedMissing = names(table(rawData$cage_mouse)[table(rawData$cage_mouse)<2])
#filter on paired missing 
rawData = rawData %>% 
  filter(!cage_mouse%in%pairedMissing)

#data overview
table(rawData$Sensitization);(table(metadata$Sensitization))
table(rawData$Strain);(table(metadata$Strain))
table(rawData$Reactors);(table(metadata$Reactors))
clo(table(rawData$IgE));clo(table(metadata$IgE))


#data parms
mouse_cages = as.character(unique(rawData$mouse_cage))
strains = as.character(unique(rawData$Strain))
treatments = as.character(unique(rawData$Sensitization))
ids = as.character(unique(rawData$cage_mouse))

#get features
features = rawData %>% 
  dplyr::select(starts_with(match = "V"))




###########################
## CC027H Anlaysis ###
###########################
cc027_dat = rawData %>% 
  filter(Strain == strains[2])
features = cc027_dat %>% 
  dplyr::select(starts_with(match = "V"))
#pre
cc027_pre = cc027_dat %>% 
  filter(Time=="pre")
features.pre = cc027_pre %>% 
  dplyr::select(starts_with(match = "V"))
cc027_features.pre = clo(features.pre + 5e-07 ) #add pseudo count
rownames(cc027_features.pre) = cc027_pre$mouse_cage
#post
cc027_post = cc027_dat %>% 
  filter(Time=="post") %>% 
  filter(mouse_cage %in% cc027_pre$mouse_cage)
features.post = cc027_post %>% 
  dplyr::select(starts_with(match = "V"))
cc027_features.post = clo(features.post + 5e-07)
rownames(cc027_features.post) = cc027_post$mouse_cage
#QC
sum(!cc027_post$mouse_cage == cc027_pre$mouse_cage) #should be 0
#compute diff in pre post composition
post_pre.pertubation = clo(cc027_features.post/cc027_features.pre)
post_pre.pertubation = data.frame(mouse_cage = rownames(post_pre.pertubation),post_pre.pertubation)
post_pre.pertubation = left_join(cc027_post[,1:23],post_pre.pertubation)
features.pert = post_pre.pertubation %>% 
  dplyr::select(starts_with(match = "V"))

##  perturb Data ###
#By Sens
id_sen = paste(post_pre.pertubation$Strain,post_pre.pertubation$Sensitization,sep = "_")
df = data.frame(Sensitzation = id_sen,features.pert)
cc027_perturb = df
cc027_perturb_meta = cc027_pre %>% 
  dplyr::select(-starts_with(match = "V"))
##########################################################





###########################
## C3H Anlaysis ###
###########################
C3H_dat = rawData %>% 
  filter(Strain == strains[1])
features = C3H_dat %>% 
  dplyr::select(starts_with(match = "V"))
#pre
C3H_pre = C3H_dat %>% 
  filter(Time=="pre")
features.pre = C3H_pre %>% 
  dplyr::select(starts_with(match = "V"))
C3H_features.pre = clo(features.pre + 5e-07 ) #add pseudo count
rownames(C3H_features.pre) = C3H_pre$mouse_cage
#post
C3H_post = C3H_dat %>% 
  filter(Time=="post") %>% 
  filter(mouse_cage %in% C3H_pre$mouse_cage)
features.post = C3H_post %>% 
  dplyr::select(starts_with(match = "V"))
C3H_features.post = clo(features.post + 5e-07)
rownames(C3H_features.post) = C3H_post$mouse_cage
#QC
sum(!C3H_post$mouse_cage == C3H_pre$mouse_cage) #should be 0
#compute diff in pre post composition
post_pre.pertubation = clo(C3H_features.post/C3H_features.pre)
post_pre.pertubation = data.frame(mouse_cage = rownames(post_pre.pertubation),post_pre.pertubation)
post_pre.pertubation = left_join(C3H_post[,1:23],post_pre.pertubation)
features.pert = post_pre.pertubation %>% 
  dplyr::select(starts_with(match = "V"))

##  perturb Data ###
#By Sens
id_sen = paste(post_pre.pertubation$Strain,post_pre.pertubation$Sensitization,sep = "_")
df = data.frame(Sensitzation = id_sen,features.pert)
C3H_perturb = df
C3H_perturb_meta = C3H_pre %>% 
  dplyr::select(-starts_with(match = "V"))
##########################################################




##########################################################
### Pre Composition CC027 vs C3h
##########################################################
C3H.pre = data.frame(Strain = C3H_pre$Strain,C3H_features.pre)
CC027.pre = data.frame(Strain = cc027_pre$Strain,cc027_features.pre)
### Data
df = rbind(C3H.pre,CC027.pre)
df.metadata = rbind(C3H_pre,cc027_pre)
#write output data 
write_csv(x = df,path = "Output/preCompositionByStrain.csv")
write_csv(x = df.metadata,path = "Output/preCompositionByStrainMetadata.csv")
##########################################################



###################################
### Pertubation CC027 and C3h #####
###################################
df = rbind(cc027_perturb,C3H_perturb)
df.metadata = rbind(cc027_perturb_meta,C3H_perturb_meta)
#write output data 
write_csv(x = df,path = "Output/pertubationByStrain.csv")
write_csv(x = df.metadata,path = "Output/pertubationByStrainMetadata.csv")
###################################
