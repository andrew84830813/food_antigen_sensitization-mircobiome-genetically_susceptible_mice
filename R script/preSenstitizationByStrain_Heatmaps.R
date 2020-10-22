# Author:   Andrew Hinton
# Email:    andrew84@email.unc.edu
# github:   https://github.com/andrew84830813/food_antigen_sensitization-mircobiome-genetically_susceptible_mice.git




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
library(ggsci)

### start parallel cluster
ncores = detectCores()
clus <- makeCluster(ncores-1) 
registerDoParallel(clus) 



### Read Data
dat = read_csv("Output/preCompositionByStrain.csv")
df.metadata = read_csv("Output/preCompositionByStrainMetadata.csv")


#Load Helper Functions
functions_Path = "Functions/"
setwd(functions_Path)
file_names <- dir() 
for(i in 1:length(file_names)){
  message(i)
  source(file = file_names[i] )
}
setwd("..")


#weighted Analysis
weighted = T
weightALRs = F

# ### compute clr matrix
clr_mat = data.frame(easyCODA::CLR(data = dat[,-1],weight = weighted)$LR)
d = parallelDist::parDist(as.matrix(clr_mat))

### compute plr matrix
# plr_mat = data.frame(easyCODA::LR(data = dat[,-1],weight = weighted)$LR)
# d = parallelDist::parDist(as.matrix(plr_mat))

### by sample pairwise distance matrix
dat = data.frame(dat)
labels = dat[,1]


##########################################################-
### Calculate Procrustes(Greenacre) and PseduoF for all denoms in Additive Log-ratio (ALR) transform ####
##########################################################-
cn = colnames(dat[,-1])
alr.procrustes = foreach(i = 1:ncol(dat[,-1]),.combine = rbind,.packages = c("compositions","vegan"))%dopar%{
  # compute alr transformation using the ith feature as denom
  if(weightALRs ==T){
    #apply weights
    ph.alr = easyCODA::ALR(data = dat[,-1],denom = i,weight = weighted)
    wts.w = ph.alr$LR.wt
    ph.alr  = ph.alr$LR
    colnm = colnames(ph.alr)
    ph.alr = sapply(1:length(wts.w), function(x) ph.alr[,x]*wts.w[x])
    colnames(ph.alr) = colnm
  }else{
    ph.alr = easyCODA::ALR(data = dat[,-1],denom = i,weight = weighted)$LR
  }
 
  
  ph.d = parallelDist::parDist(as.matrix(ph.alr))
  # Compute procrustes correlation between alr and full logratio configurations distance matices 
  proc = vegan::protest(X = ph.d,Y = d,permutations = 1)
  # compuet pseudo f with adonis function
  a.df = data.frame(Strain = labels)
  pmv = adonis2(ph.d~Strain,data = a.df,permutations = 1)
  #Store data
  data.frame(Denom = cn[i],i = i,
             Correlation = proc$t0,
             pval = proc$signif,
             pseudoF = pmv$F[1],
             permanova_pval = pmv$`Pr(>F)`[1])
}
##########################################################-



##########################################################-
### Select top AlR denom, that is ALR denom. that has the highest correlation to the full logratio configuration
##########################################################-
topALRs = top_n(alr.procrustes,n = 1,wt = Correlation)
i = topALRs$i[1]
### compute weights (avg proportions) of each taxa 
weights = data.frame(Denom = names(colMeans(dat[,-1])), weight = colMeans(dat[,-1]))

### Apply ALR transformation using the top denom from above
#procrusted table
alr.procrustes1 = left_join(alr.procrustes,weights)
colnames(alr.procrustes1)[1] = "ID"
alr.procrustes1 = left_join(cnames,alr.procrustes1)
reference = alr.procrustes1[i,]
##########################################################-


##########################################################-
### Procrustes Correlation Visualization
##########################################################-
#compute alr transformation using the ith feature as denom
if(weightALRs ==T){
  #apply weights
  ph.alr = easyCODA::ALR(data = dat[,-1],denom = i,weight = weighted)
  wts.w = ph.alr$LR.wt
  ph.alr  = ph.alr$LR
  colnm = colnames(ph.alr)
  ph.alr = sapply(1:length(wts.w), function(x) ph.alr[,x]*wts.w[x])
  colnames(ph.alr) = colnm
}else{
  ph.alr = easyCODA::ALR(data = dat[,-1],denom = i,weight = weighted)$LR
}
plot(alr.procrustes1$weight,alr.procrustes1$Correlation)
#compute pairwise sample distance using ith alr tranformation
ph.d = dist(ph.alr)
#Compute procrustes correlation between alr and full logratio configurations distance matices 
proc = vegan::protest(X = ph.d,Y = d,permutations = 100)
dvec = as.vector(as.matrix(d)[lower.tri(as.matrix(d),diag = F)])
ph.vec = as.vector(as.matrix(ph.d)[lower.tri(as.matrix(ph.d),diag = F)])
#Plot Output
tiff(filename = "Figures/Supplemental/preSentization_topDenomProcrutes.tiff",width = 4,height = 3,units = "in",res = 300)
data.frame(d = dvec,alr.dist = ph.vec) %>%
  ggplot(aes(d,alr.dist)) +
  geom_point(alpha = .5)+
  theme_classic()+
  xlab("Logratio distances")+
  ylab(paste("ALR:V",reference$i, " Distances",sep = ""))+
  labs(caption = paste("Procrustes Correlation = ",round(proc$t0,digits = 4)))+
  ggtitle("Pre-sensitization",subtitle =  paste(reference$ID,": ",reference$Order,reference$Family,reference$Genus,reference$Species,sep = ""))+
  theme(plot.subtitle = element_text(size = 8),plot.title = element_text(face = "bold"))
dev.off()
#Write Supplemental PseudoF and Procrustes Table
write_csv(x = alr.procrustes1,path = "Output/supplementalTable_preSentitization_procrustesAnalysis.csv")

## Visualize Range og Proc. Corr.
tiff(filename = "Figures/Supplemental/preSentization_ProcrutesOPT.tiff",width = 4,height = 3,units = "in",res = 300)
alr.procrustes1$Label = if_else(alr.procrustes1$ID==reference$ID,reference$ID,NULL)
ggplot(alr.procrustes1,aes(reorder(ID,Correlation),Correlation,group = 1,label = Label))+
  geom_point(alpha = 1,size = 1.5,shape = 24)+
  geom_label(nudge_y = -.1,nudge_x = -12.1)+
  geom_segment(x = 230 , y = 0 , xend = 230 , yend = reference$Correlation  ,col = "red",lty = "dashed")+
  geom_line()+
  geom_hline(yintercept = reference$Correlation,lty = "dashed",col= "red")+
  ggpubr::theme_classic2()+
  theme(axis.text.x = element_blank())+
  ylab("Procrustes Correlation")+
  scale_y_continuous(limits = c(0.25,1))+  
  xlab(expression(paste(Taxa[i],sep="")))
dev.off()
###############################################################3


#######################################################
#### <ML_Anlalysis> Output already stores in output (Runtime can be lengthy depending on cpu)
#######################################################
mlData = data.frame(React = dat[,1],ph.alr) 
nrepeats =  10
nfolds = 10
#cvPerformance = read_csv("Output/presen_ALR_cvTrainPerformance.csv")
#nullPerformance = read_csv("Output/presen_ALR_nullTrainPerformance.csv")

#=============================
## Empirical Performance
#=============================
permuteLabels = F
cvPerformance = data.frame()
preds.df = data.frame()
for(s in 1:nrepeats){
  foldData = kfoldDataPartition(df = mlData,
                                kfold = nfolds,
                                permuteLabel = permuteLabels,
                                seed = s)
  #=====================================================-
  ## Train Model
  #=====================================================-
  set.seed(s)
  modelPerf = foreach(f = 1:nfolds,.combine = rbind)%do%{
    suppressMessages(suppressWarnings({
    ##  select fold data for train and test splilt
    trainData = (foldData[[f]]$xtrain_combinedFolds[,-1])
    ytrain = factor(foldData[[f]]$xtrain_combinedFolds[,1])
    classes = as.character(unique(ytrain))
    ## retrieve test data from kfold partition
    testData = (foldData[[f]]$xtest_kthFold[,-1])
    ytest = factor(foldData[[f]]$xtest_kthFold[,1])
    #boruta Feature Selection
    b = Boruta::Boruta(x = trainData,
                       y = as.factor(ytrain),
                       doTrace = 0,maxRun = 100)
    dec = data.frame(Ratio = names(b$finalDecision),Decision = b$finalDecision)
    keep = dec %>% 
      filter(Decision!="Rejected")
    trainData = subset(trainData,select = keep$Ratio)
    testData = subset(testData,select = keep$Ratio)
    ## Train Model
    glm.mdl1=train(x = trainData,
                   y =ytrain,
                   #method = "glmnet",
                   method = "rf",importance = T,proximity = F,ntree = 5000,tuneGrid = expand.grid(mtry = 2),
                   metric = "ROC",
                   trControl = train_control
    )
    
   
    ## performance
    preds = predict.train(glm.mdl1,testData,type = "prob")[,1]
    rocobj = pROC::roc(ytest,preds)
    ## AUC
    dir = rocobj$direction
    thrhold = pROC::coords(rocobj,x="best",input = "threshold",best.method = "youden")
    if(dir == "<"){
      preds1 = as.factor(if_else(preds<thrhold[1][,1],rocobj$levels[1],rocobj$levels[2]))
      AC = confusionMatrix(preds1,factor(ytest,levels = c(rocobj$levels[1],rocobj$levels[2])))
    }else{
      preds1 = as.factor(if_else(preds>thrhold[1][,1],rocobj$levels[1],rocobj$levels[2]))
      AC = confusionMatrix(preds1,factor(ytest,levels = c(rocobj$levels[1],rocobj$levels[2])))
    }
    ## Procrustes Analysis
    rows = as.numeric(foldData[[f]]$xtrain_IDs$ID)
    procrCorrTrain = pwDistProcCorrelation(baseMatrix = dat[rows,-1],stretchMatrix = trainData)
    rows = as.numeric(foldData[[f]]$xtest_IDs$ID)
    procrCorrTest = pwDistProcCorrelation(baseMatrix = dat[rows,-1],stretchMatrix = testData)
    auroc = pROC::auc(ytest,preds)
    }))
    message("seed-",s,", fold-",f)
    preds.df = rbind(preds.df,data.frame(Fold = f,Seed = s,Prob = preds,Label = ytest))
    
    ## Performance
    data.frame(fold = f,seed = s,
               AUC =auroc,
                t(AC$overall),
               Train_procrustedCorr = procrCorrTrain,
               Test_procrustedCorr = procrCorrTest,
               t(AC$byClass),
               numFeatures = ncol(trainData),
               getTrainPerf(glm.mdl1),glm.mdl1$bestTune,
               positiveClass = AC$positive)
    
  }
  #write sth repeat kfold test performance
  cvPerformance  = rbind(cvPerformance,modelPerf)
}


## performance
cvPerformance %>% 
  summarise_all(.funs = mean)
cvPerformance$Labels = 'True'
#write_csv(x = cvPerformance,path = "Output/presen_ALR_cvTrainPerformance.csv")
#write_csv(x = preds.df,path = "Output/presen_ALR_TrainPredictions.csv")

#=============================
## Null Performance (Permuted Labels) 
#=============================
permuteLabels = T
nrepeats = 100
nullPerformance = data.frame()
null.preds.df = data.frame()
for(s in 1:nrepeats){
  foldData = kfoldDataPartition(df = mlData,
                                kfold = nfolds,
                                permuteLabel = permuteLabels,
                                seed = s)
  #=====================================================-
  ## Train Model
  #=====================================================-
  set.seed(s)
  modelPerf = foreach(f = 1:nfolds,.combine = rbind)%do%{
    preds = c()
   suppressMessages(suppressWarnings({
      ##  select fold data for train and test splilt
      trainData = (foldData[[f]]$xtrain_combinedFolds[,-1])
      ytrain = factor(foldData[[f]]$xtrain_combinedFolds[,1])
      classes = as.character(unique(ytrain))
      ## retrieve test data from kfold partition
      testData = (foldData[[f]]$xtest_kthFold[,-1])
      ytest = factor(foldData[[f]]$xtest_kthFold[,1])
      #boruta Feature Selection
      b = Boruta::Boruta(x = trainData,
                         y = as.factor(ytrain),
                         doTrace = 0,maxRun = 100)
      dec = data.frame(Ratio = names(b$finalDecision),Decision = b$finalDecision)
      keep = dec %>% 
        filter(Decision!="Rejected")
      #train with all features if boruta returns no features
      if(nrow(keep)>0){
        trainData = subset(trainData,select = keep$Ratio)
        testData = subset(testData,select = keep$Ratio) 
      }
      ## Train Model
      glm.mdl1=train(x = trainData,
                     y =ytrain,
                     #method = "glmnet",
                     method = "rf",importance = T,proximity = F,ntree = 5000,tuneGrid = expand.grid(mtry = 2),
                     metric = "AUC",
                     trControl = train_control
      )
      ## performance
      preds = predict.train(glm.mdl1,testData,type = "prob")[,1]
      rocobj = pROC::roc(ytest,preds)
      ## AUC
      dir = rocobj$direction
      thrhold = pROC::coords(rocobj,x="best",input = "threshold",best.method = "youden")
      if(dir == "<"){
        preds1 = as.factor(if_else(preds<thrhold[1][,1],rocobj$levels[1],rocobj$levels[2]))
        AC = confusionMatrix(preds1,factor(ytest,levels = c(rocobj$levels[1],rocobj$levels[2])))
      }else{
        preds1 = as.factor(if_else(preds>thrhold[1][,1],rocobj$levels[1],rocobj$levels[2]))
        AC = confusionMatrix(preds1,factor(ytest,levels = c(rocobj$levels[1],rocobj$levels[2])))
      }
      ## Procrustes Analysis
      # Train
      rows = as.numeric(foldData[[f]]$xtrain_IDs$ID)
      errorStatus = catchError(expr =pwDistProcCorrelation(baseMatrix = dat[rows,-1],stretchMatrix = trainData) )
      if(errorStatus==T){
        procrCorrTrain = NA
      }else{
        procrCorrTrain = pwDistProcCorrelation(baseMatrix = dat[rows,-1],stretchMatrix = trainData) 
      }
      # Test
      rows = as.numeric(foldData[[f]]$xtest_IDs$ID)
      errorStatus = catchError(expr = pwDistProcCorrelation(baseMatrix = dat[rows,-1],stretchMatrix = testData) )
      if(errorStatus==T){
        procrCorrTest = NA
      }else{
        procrCorrTest = pwDistProcCorrelation(baseMatrix = dat[rows,-1],stretchMatrix = testData)  
      }
      ## class perf
      auroc = pROC::auc(ytest,preds)
    }))
    message("seed-",s,", fold-",f)
    null.preds.df = rbind(null.preds.df,data.frame(Fold = f,Seed = s,Prob = preds,Label = ytest))
    
    ## Model Performance
    data.frame(fold = f,seed = s,
               AUC = auroc,
               t(AC$overall),
               Train_procrustedCorr = procrCorrTrain,
               Test_procrustedCorr = procrCorrTest,
               t(AC$byClass),
               numFeatures = ncol(trainData),
               getTrainPerf(glm.mdl1),glm.mdl1$bestTune,
               positiveClass = AC$positive)
    
  }
  #write sth repeat kfold test performance
  nullPerformance  = rbind(nullPerformance,modelPerf)
}
## performance
np_ = nullPerformance %>% 
  group_by(seed) %>% 
  summarise_all(.funs = mean)
nullPerformance$Labels = 'Permuted'
#write_csv(x = nullPerformance,path = "Output/presen_ALR_nullTrainPerformance.csv")
#write_csv(x = null.preds.df,path = "Output/presen_ALR_nullTrainPredictions.csv")

#Visualize
permutationPerf = nullPerformance %>% 
 group_by(seed,Labels) %>% 
  summarise_all(.funs = mean)

#permutation pval
pval = (sum(mean(cvPerformance$AUC)<permutationPerf$AUC)+1) / (nrow(permutationPerf)+1)
tiff(filename = "Figures/Main/prePertubationMicrobiome_permutationTesting.tiff",width = 4,height = 4,units = "in",res = 300)
ggplot(permutationPerf,aes(AUC))+
  geom_histogram(aes(fill = Labels),col = "white")+
  #scale_fill_manual(values = "grey40")+
  geom_vline(aes(xintercept = mean(cvPerformance$AUC),col = "Classification AUROC"),
             lty = "dashed",size = 1,
             show.legend = T)+
  geom_vline(aes(xintercept = mean(nullPerformance$AUC),col = "Chance AUROC"),
             lty = "dashed",size = 1,
             show.legend = T)+
  scale_color_manual(values = c("Black","Green"))+
  scale_fill_manual(values = "grey40" )+
  theme_classic()+
  scale_x_continuous(limits = c(0.5,1))+
  ggtitle("Predicting strain using pre-senstization gut micobiome")+
  xlab("10-fold Cross Validated AUROC")+
  theme(axis.text = element_text(hjust = 1,size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.key.size = unit(x = .1,units = "in"),legend.box = "vertical",
        legend.spacing.x = unit(x = .05,units = "in"),legend.margin = margin(),
        legend.spacing.y = unit(x = 0.00001,units = "in"),
        legend.key.width = unit(x = .1,units = "in"),
        legend.title = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = 8,face = "bold"),
        plot.caption = element_text(size = 8,face = "italic"))+
  labs(caption =  paste("Permutation Test(n=100), pval =",round(pval,5)))
dev.off()

## ROC Curve
#True lables
nrepeats = 10
ddf = data.frame()
for(s in 1:nrepeats){
  foldData = kfoldDataPartition(df = mlData,
                                kfold = nfolds,
                                permuteLabel = permuteLabels,
                                seed = s)
  
  for(f in 1:nfolds){
    ph = data.frame(Seed = s,ytest = foldData[[f]]$xtest_kthFold[,1],foldData[[f]]$xtest_IDs)  
    ddf = rbind(ddf,ph)
  }
  
}
test = cbind(ddf,preds.df[,-2])
test = test %>% 
  group_by(ytest,ID) %>% 
  summarise(meanProb = mean(Prob))
roco = roc(test$ytest,test$meanProb,ci=TRUE, plot=FALSE)
ciobj <- ci.se(roco)
curve.label = paste0("True Labels(Mean AUC = ",round(roco$auc,3),")")
dat.ci.true <- data.frame(FPR = 1-as.numeric(rownames(ciobj)),TPR = ciobj[, 1],upper = ciobj[, 3],Trial = curve.label)
y = roco$sensitivities
x = 1 - roco$specificities
true.label = data.frame(TPR = y,FPR = x) 
true.label = true.label %>% 
  arrange(FPR,TPR)
true.label$Trial = curve.label


#Null Labels
nrepeats = 100
ddf.null = data.frame()
permuteLabels = T
for(s in 1:nrepeats){
  foldData = kfoldDataPartition(df = mlData,
                                kfold = nfolds,
                                permuteLabel = permuteLabels,
                                seed = s)
  
  for(f in 1:nfolds){
    ph = data.frame(Seed = s,ytest = foldData[[f]]$xtest_kthFold[,1],foldData[[f]]$xtest_IDs)  
    ddf.null = rbind(ddf.null,ph)
  }
  
}
test.null = cbind(ddf.null,null.preds.df[,-2])
test.null = test.null %>% 
  group_by(ytest,ID) %>% 
  summarise(meanProb = mean(Prob))
roco.null = roc(test.null$ytest,test.null$meanProb,ci=TRUE, plot=FALSE)
ciobj <- ci.se(roco.null)
curve.label = paste0("Permuted Labels(Mean AUC = ",round(roco.null$auc,3),")")
dat.ci <- data.frame(FPR = 1-as.numeric(rownames(ciobj)),TPR = ciobj[, 1],upper = ciobj[, 3],Trial = curve.label )
y = roco.null$sensitivities
x = 1 - roco.null$specificities
null.curve = data.frame(TPR = y,FPR = x) 
null.curve = null.curve %>% 
  arrange(FPR,TPR)
null.curve$Trial = curve.label

#Final ROC Curves
roc.curves = rbind(true.label,null.curve)
tiff(filename = "Figures/Main/preSen_rocCurveByStrain.tiff",width = 3.5,height = 2.75,units = "in",res = 300)
ggplot(roc.curves,aes(FPR,TPR,col = Trial,fill = Trial))+
  geom_line(alpha = 1,size = 1)+
  theme_classic()+
  #geom_ribbon(data = dat.ci, aes(x = FPR, ymin = TPR, ymax = upper), alpha= 0.2,colour = NA)+
  geom_abline(slope=1, intercept = 0, alpha=1, color = "black",size = .75) +
  theme(legend.position = c(.65,.15))+
  xlab("False Positive Rate")+
  ylab("True Positive Rate")+
  theme(axis.text = element_text(hjust = 1,size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.margin = margin(),legend.key.size = unit(.1,units = "in"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        plot.caption = element_text(size = 8,face = "italic"))
dev.off()
######################################################################################################






#######################################################
### Identify Min Signature
#######################################################
ph.alr = data.frame(ph.alr)
set.seed(08272008)
b = Boruta::Boruta(x = ph.alr,
                   y = as.factor(dat[,1]),
                   doTrace = 0,maxRun = 100)
dec = data.frame(Ratio = names(b$finalDecision),Decision = b$finalDecision)
keep = dec %>% 
  filter(Decision!="Rejected")
alr.df = subset(ph.alr,select = keep$Ratio)
#######################################################







###############################################################
### Statistical Analysis ####
###############################################################
alr.d =  parallelDist::parDist(as.matrix(alr.df))
#permanova
a.df = data.frame(Type = labels)
pmv = adonis2(alr.d~Type,data = a.df,permutations = 5000)
(pmv)
#dispersion test
mod = betadisper(as.dist(ph.d),group = labels)
anova(mod)
plot(mod)
mod.HSD <- TukeyHSD(mod)
dis = data.frame(Treatment = labels,Dist = mod$distances)
#tiff(filename = "Figures/Main/pre-postPertubationMicrobiome_betadisp.tiff",width = 3,height = 2.25,units = "in",res = 300)
ggplot(dis,aes(Treatment,Dist,fill = Treatment))+
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot()+
  geom_jitter(width = .1,alpha = .6)+
  theme_classic()+
  #scale_fill_brewer(palette = "Paired")+
  #labs(caption = paste('PERMDISP2 F=',round(mod$`F value`[1],5),", p=",round(mod$`Pr(>F)`[1],5),sep = ""))+
  theme(axis.text = element_text(hjust = 1,size = 8,colour = "black"),legend.position = "none",
        axis.title.x = element_blank(),axis.title.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45),plot.caption = element_text(size = 8))+
  ylab("Distance to Centroid")
###############################################################




###############################################################
### Principle Coordinate Analysis with using ALR with pseudo F maximizing denom ####
###############################################################
mds = cmdscale(ph.d,eig = T)
mds.df = data.frame(Strain = dat[,1],mds$points)
#Compute PCoA variation
pcoa_var = mds$eig[mds$eig>0]
pcoa_var = (pcoa_var / sum(pcoa_var))*100
#Write Figure Output
tiff(filename = "Figures/Main/presensitizationMicrobiome_withALR-PCoA.tiff",width = 3,height = 2.75,units = "in",res = 300)
ggplot(mds.df,aes(X1,X2,fill = Strain))+
  geom_point(size = 3,pch = 21,alpha = .7)+
  theme_classic()+
  labs(caption = paste('PERMANOVA F=',round(pmv$F,5),", p=",round(pmv$`Pr(>F)`,5),sep = ""))+
  theme(axis.text = element_text(hjust = 1,size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.position = "none",
        legend.margin = margin(),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        plot.caption = element_text(size = 8,face = "italic"))+
  xlab(paste("PCo1 (",round(pcoa_var[1],digits = 2),"% Variation)",sep = ""))+
  ylab(paste("PCo2 (",round(pcoa_var[2],digits = 2),"% Variation)",sep = ""))+
  scale_fill_manual(values = c("Red","Blue"))
dev.off()
##############################################################



### Presenstization Heatmap PreProcess #####
keepRatios = data.frame(Group = labels,alr.df)
ncolors = n_distinct(keepRatios[,1])
colors = c("Red","Blue")
color.df = data.frame(Color = colors,
                      Group = as.character(unique(keepRatios[,1])))
color.df = left_join(keepRatios,color.df)
mat = as.matrix((keepRatios[,-1]))
rownames(mat) = paste(keepRatios[,1],1:nrow(keepRatios),sep = "YlGn")
RColorBrewer::brewer.pal.info
col <- colorRampPalette(RColorBrewer::brewer.pal(10, "PRGn"))(1000)
col <- colorRampPalette(RColorBrewer::brewer.pal(10, "PRGn")[10:1])(1000)
mat = t(mat)
rownames(mat) = str_replace(string = rownames(mat),pattern = "\\.",replacement = "/")
sc = color.df$Color
# heatmap microbiome composition
tiff(filename = "Figures/Main/presensitizationMicrobiome_withALR-Heatmap.tiff",width = 9,height = 6,units = "in",res = 300)
gplots::heatmap.2( mat,
                   col = col,#viridis::viridis(n = 1000,direction = -1,option = "D"),#gplots::redblue(n = 1000) ,#viridis::viridis(n = 1000,option = "D"),#,
                   Rowv = TRUE,
                   margins = c(2, 5.5), labCol = FALSE,labRow = rownames(mat),
                   hclustfun = function(x) hclust(x, method = "ward.D2"),
                   distfun = function(x) parallelDist::parallelDist(x,method = "euclidean"),
                   key=TRUE, symkey=TRUE, 
                   scale = "none",
                   sepwidth=c(0.01,0.01),
                   sepcolor="black",colsep=1:ncol(mat),rowsep=1:nrow(mat),
                   density.info="density", trace="none",
                   ColSideColors = sc,xlab = 'Sample',ylab = "Log-ratio",
                   #main = "Presensitization Gut Mircobiome",
                   #colCol = "white",
                   key.title = "Log-ratio",key.xlab = paste("log( v_i / ",reference$ID,")",sep = "" ),
                   #RowSideColors = cc,
                   cexRow = 0.75,cexCol = .75
                   
)
dev.off()


#=====================================-
### Kruskall Test ####
#=====================================-
keepRatios = data.frame(Group = dat[,1],alr.df)

gather.alr2 = keepRatios %>% 
  gather(key = "Ratio",value = "value",2:ncol(keepRatios)) %>% 
  separate(2,into = c("ID","denom"),sep = "\\.",remove = F)
unqTaxa = unique(gather.alr2$Ratio)

krus.pval = data.frame()
for(i in 1:length(unqTaxa)){
  ph = gather.alr2 %>% 
    filter(Ratio == unqTaxa[i])
  kt = wilcox.test(value ~ Group, data = ph)
  ph.results = data.frame(Ratio = unqTaxa[i],pval = kt$p.value)
  #bind data
  krus.pval = rbind(krus.pval,ph.results)
}
krus.pval$p.adj = p.adjust(krus.pval$pval,method = "BH")
krus.pval$signf =   if_else(krus.pval$pval<=0.05,1,0)
sum(krus.pval$signf,na.rm = T)

#top N signficant ratios by pvalues
krus.pval_topN = krus.pval %>% 
  top_n(12,-pval) %>% 
  arrange(desc(-pval))
krus.pval_topN = separate(krus.pval_topN,1,into = c("ID","denom"),sep = "\\.",remove = F)
krus.pval_topN  = left_join(krus.pval_topN,cnames) %>% 
  arrange(desc(-pval))
krus.pval_topN$ofg = paste(krus.pval_topN$Order,krus.pval_topN$Family,krus.pval_topN$Genus,krus.pval_topN$Species,sep = "|")

#visualize difference between strain boxplot
keepRatios = data.frame(Strain = dat[,1],ph.alr)
keyRatios = keepRatios %>% 
  dplyr::select(Strain,one_of(krus.pval_topN$Ratio))
keyRatios = gather(keyRatios,"Ratio","value",2:ncol(keyRatios))
keyRatios = separate(keyRatios,2,into = c("ID","denom"),sep = "\\.",remove = F)
keyRatios$ID = factor(keyRatios$ID,levels = krus.pval_topN$ID)
keyRatios = left_join(keyRatios,cnames)
keyRatios$ofg = paste(keyRatios$Order,keyRatios$Family,keyRatios$Genus,keyRatios$Species,sep = "|")
keyRatios$ofg = factor(keyRatios$ofg,levels = krus.pval_topN$ofg)
keyRatios$ofg = str_replace_all(keyRatios$ofg,pattern = "__",replacement = "_")

ggplot(keyRatios,aes(Strain,value,col = Strain)) +
  geom_boxplot(aes(fill = Strain),outlier.shape = NA,col = "black",alpha = .5)+
  geom_jitter(aes(fill = Strain),pch = 21,width = .1,alpha = 1,size = 3,col = "black")+
  facet_wrap(.~ofg,nrow = 3)+
  scale_fill_manual(values = c("red","blue"))+
  scale_color_manual(values = c("red","blue"))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(caption = paste("Ref. Denom. = V",reference$i," | ",reference$col,sep = ""))+
  ylab("Weighted Log-ratio")
############################################################




############################################################
### Mean LR ####
############################################################
library(ggpubr)
library(rstatix)
keepRatios = data.frame(ID =df.metadata$index,Group = dat[,1],alr.df)
gather.alr1 = keepRatios %>% 
  gather(key = "Ratio",value = "value",3:ncol(keepRatios)) %>% 
  #filter(Ratio%in% krus.pval_topN$Ratio) %>% 
  mutate(Ratio = str_replace(string = Ratio,pattern = "\\.",replacement = " / "))
#Enriched Analysis
enrich = keepRatios %>% 
  gather(key = "Ratio",value = "value",3:ncol(keepRatios)) %>% 
  #filter(Ratio%in% krus.pval_topN$Ratio) %>% 
  group_by(Group,Ratio) %>% 
  summarise(meanLR = mean(value)) %>% 
  mutate(Ratio = str_replace(string = Ratio,pattern = "\\.",replacement = " / ")) %>% 
  arrange(desc(meanLR)) 
enrich1 = enrich %>% 
  tidyr::spread("Group","meanLR") %>% 
  mutate(Enriched_C3H = if_else(C3H>CC027,1,0))
c3h_enriched = enrich1$Ratio[enrich1$Enriched_C3H==1]
cc027_enriched = enrich1$Ratio[enrich1$Enriched_C3H!=1]
#=========================
#C3H Enriched
#=========================
gather.alr = gather.alr1 %>% 
  filter(Ratio %in%c3h_enriched)
#stat test
stat.test <- gather.alr %>%
  group_by(Ratio) %>%
  rstatix::t_test(value ~ Group) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj")
# Create a bar plot with error bars (mean +/- sd)
bp <- ggpubr::ggbarplot(
  gather.alr, x = "Ratio", y = "value", add = c("mean_ci"),
   size = .25,add.params = list(size=.25, position = position_dodge(0.8)),
  fill= "Group", palette = c("Red", "Blue"),ylab = "log ( Ratio )",
  position = position_dodge(0.8),
  ggtheme = theme(legend.position = "top",
                  axis.line = element_line(),
                  axis.title = element_text(face = "bold"),
                  panel.background = element_blank())
)
# Add p-values onto the bar plots
stat.test <- stat.test %>%
  rstatix::add_xy_position(fun = "mean_ci", x = "Ratio",) 
p = bp + ggpubr::stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0.01,y.position = 11.5,size = 2
)
tiff(filename = "Figures/Supplemental/preSentization_logratioSignatureC3HEnriched.tiff",width = 6,height = 3,units = "in",res = 300)
ggpubr::ggpar(p,x.text.angle = 45,font.main = 8,font.xtickslab  = 8)
dev.off()
#=========================
#CC027 Enriched
#=========================
gather.alr = gather.alr1 %>% 
  filter(Ratio %in%cc027_enriched)
#stat test
stat.test <- gather.alr %>%
  group_by(Ratio) %>%
  rstatix::t_test(value ~ Group) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj")
# Create a bar plot with error bars (mean +/- sd)
bp <- ggpubr::ggbarplot(
  gather.alr, x = "Ratio", y = "value", add = c("mean_ci"),
  size = .25,add.params = list(size=.25, position = position_dodge(0.8)),
  fill= "Group", palette = c("Red", "Blue"),ylab = "log ( Ratio )",
  position = position_dodge(0.8),
  ggtheme = theme(legend.position = "top",
                  axis.line = element_line(),
                  axis.title = element_text(face = "bold"),
                  panel.background = element_blank())
)
# Add p-values onto the bar plots
stat.test <- stat.test %>%
  rstatix::add_xy_position(fun = "mean_ci", x = "Ratio",) 
p = bp + ggpubr::stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0.01,y.position = 11.5,size = 2
)
tiff(filename = "Figures/Supplemental/preSentization_logratioSignatureCC027Enriched.tiff",width = 6,height = 3,units = "in",res = 300)
ggpubr::ggpar(p,x.text.angle = 45,font.main = 8,font.xtickslab  = 8,font.caption = 6,
      caption = paste("Ref. Denom. = V",reference$i," | ",reference$col,sep = ""))
dev.off()





############################################################
###  Train Final Model to get feature importance
############################################################
set.seed(08272008)
finalModel = train(x = alr.df,
                   y =dat[,1],
                   method = "rf",importance = T,proximity = F,ntree = 5000,tuneGrid = expand.grid(mtry = 2),
                   metric = "AUC",
                   trControl = train_control
)
vi = varImp(finalModel)
vi.df = data.frame(Ratio = rownames(vi$importance),Importance = vi$importance[,2] )
vi.df$Ratio = str_replace(string = vi.df$Ratio,pattern = "\\.",replacement = " / ")
vi.df = left_join(vi.df,enrich1)
vi.df= vi.df %>% 
  mutate(Enriched = factor(if_else(Enriched_C3H==1,"C3H_enriched","CC027_enriched")))
tiff(filename = "Figures/Supplemental/preSentization_logratioSignatureFeatureImportance.tiff",width = 5,height = 6,units = "in",res = 300)
  ggpubr::ggdotchart(vi.df, x = "Ratio", y = "Importance",ylab = "Random Forest Feature Importance",
             color = "Enriched",                                # Color by groups
             palette = c("Red", "Blue"), # Custom color palette
             sorting = "descending",                       # Sort value in descending order
             add = "segments", add.params = list(color = "black",size = .75),                            # Add segments from y = 0 to dots
             rotate = TRUE,                                # Rotate vertically
             group = "Enriched",  
             x.text.col = F,                              # Order by groups
             dot.size = 3,                             # Large dot size
             #label = round(dfm$mpg),                        # Add mpg values as dot labels
             ggtheme = theme(legend.position = "top",
                             axis.line = element_line(),axis.text = element_text(size = 8),
                             axis.title = element_text(face = "bold"),
                             panel.background = element_blank())+
               ggpubr::theme_cleveland()
  )
            
dev.off()
