# Author:   Andrew Hinton
# Email:    andrew84@email.unc.edu
# github:   https://github.com/andrew84830813/food_antigen_sensitization-mircobiome-genetically_susceptible_mice.git



### Read Data
dat = read_csv("Output/pertubationByStrain.csv")
df.metadata = read_csv("Output/pertubationByStrainMetadata.csv")

#Load Helper Functions
functions_Path = "Functions/"
setwd(functions_Path)
file_names <- dir() 
for(i in 1:length(file_names)){
  message(i)
  source(file = file_names[i] )
}
setwd("..")


### start parallel cluster
stopCluster(clus)
ncores = detectCores()
clus <- makeCluster(ncores-1) 
registerDoParallel(clus) 


#weighted Analysis
weighted = T # weight CLR matrix
weightALRs = F

### compute clr matrix
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
### Select top AlR denom, that is ALR denom. that has the highest correlation to the full logratio configuration ########
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
### Procrustes Correlation Visualization #####
##########################################################-
#compute alr transformation using the ith feature as denom
coda.alr = easyCODA::ALR(data = dat[,-1],denom = i,weight = weighted)
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
tiff(filename = "Figures/Supplemental/pre-postPertubation_topDenom_PseudoF.tiff",width = 4,height = 3,units = "in",res = 300)
data.frame(d = dvec,alr.dist = ph.vec) %>%
  ggplot(aes(d,alr.dist)) +
  geom_point(alpha = .5)+
  theme_classic()+
  xlab("Logratio distances")+
  ylab(paste("ALR:V",reference$i, " Distances",sep = ""))+
  labs(caption = paste("Procrustes Correlation = ",round(proc$t0,digits = 4)))+
  ggtitle("Pertubation",subtitle = paste(reference$ID,": ",reference$Order,reference$Family,reference$Genus,reference$Species,sep = ""))+
  theme(plot.subtitle = element_text(size = 8),plot.title = element_text(face = "bold"))
dev.off()
#Write Supplemental PseudoF and Procrustes Table
write_csv(x = alr.procrustes1,path = "Output/supplementalTable_pre-postPertubation_procrustes-pseudoF_Analysis.csv")

## Visualize Range og Proc. Corr.
tiff(filename = "Figures/Supplemental/perturbation_ProcrutesOPT.tiff",width = 4,height = 3,units = "in",res = 300)
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
#### <ML_Anlalysis>  ######
#######################################################
targetVariable = as.factor(dat[,1])
mlData = data.frame(React = targetVariable,ph.alr) 
nrepeats =  10
nfolds = 5
#Metrics
usemultiCLassACC = T
usebinACC = F
usemultiClass.AUC = T
useAUPRC = F
usebinaryAUROC = F
#AUPC
minorityClass = names(which.min(table(mlData[,1])))
majorityClass = as.character(unique(mlData[,1][mlData[,1]!=minorityClass]))
#=============================
## Empirical Performance
#=============================
permuteLabels = F
cvPerformance = data.frame()
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
    #intialize metric
    multiCLassACC = NA
    binACC = NA
    multiClass.AUC = NA
    AUPRC = NA
    binaryAUROC = NA
    
    suppressMessages(suppressWarnings({
      ##  select fold data for train and test splilt
      trainData = (foldData[[f]]$xtrain_combinedFolds[,-1])
      ytrain = factor(foldData[[f]]$xtrain_combinedFolds[,1],levels = levels(targetVariable))
      classes = as.character(unique(ytrain))
      ## retrieve test data from kfold partition
      testData = (foldData[[f]]$xtest_kthFold[,-1])
      ytest = factor(foldData[[f]]$xtest_kthFold[,1],levels = levels(targetVariable))

      #RFE Feature Selection
      fs = rfeSelection.ByMetric(train_ratio = trainData,test_ratio = testData,
                                 metric = "multiClass.ROC",
                                 ytrain = ytrain,
                                 sets = 10,
                                 impMeasure = "impurity_corrected",
                                 nreps = 10,
                                 minPercentFeatReturn = .1)
        trainData = fs$reducedTrainRatios
        testData = fs$reducedTestRatio
      
     
      ## Train Model
      glm.mdl1=train(x = trainData,
                     y =ytrain,
                     method = "rf",importance = T,proximity = F,ntree = 5000,tuneGrid = expand.grid(mtry = c(2,sqrt(ncol(trainData)))),
                     metric = "AUC",
                     trControl = train_control
      )
      
      #multiCLass Accuracy
      if(usemultiCLassACC==T){
        preds = factor(predict.train(glm.mdl1,testData),levels = levels(targetVariable))
        AC = confusionMatrix(preds,ytest)
        multiCLassACC =colMeans(AC$byClass,na.rm = T)
      }
      
      #Binary Accuracy
      if(usebinACC==T){
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
        binACC = cbind(data.frame(t(AC$overall)),data.frame(t(AC$byClass)))
      }
     
      ## Mulitclass ROC
      if(usemultiClass.AUC==T){
        preds = predict.train(glm.mdl1,testData,type = "prob")
        rocobj = pROC::multiclass.roc(ytest,preds)
        multiClass.AUC = pROC::auc(rocobj)
      }
      
      ## AUPRC
      if(useAUPRC==T){
        #area under precision recall curve
        negClass.col = which(colnames(preds)==majorityClass)
        posClass = preds[,negClass.col]
        negClass = preds[,negClass.col]
        posClass = posClass[ytest==minorityClass]
        negClass = negClass[ytest==majorityClass]
        AUPRC = PRROC::pr.curve(scores.class0 = negClass,scores.class1 = posClass,curve = TRUE)$auc.integral  
      }
      
      #Binary AUROC
      if(usebinaryAUROC==T){
        preds = predict.train(glm.mdl1,testData,type = "prob")[,1]
        binaryAUROC = pROC::auc(ytest,preds)
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
    }))
    
    message("seed-",s,", fold-",f)
    ## Performance
    data.frame(fold = f,seed = s,
               t(multiCLassACC),
               binACC,
               multiClass.AUC,
               AUPRC,
               binaryAUROC,
               Train_procrustedCorr = procrCorrTrain,
               Test_procrustedCorr = procrCorrTest,
               numFeatures = ncol(trainData),
               trainAUC = getTrainPerf(glm.mdl1)[2],glm.mdl1$bestTune)
    
  }
  #write sth repeat kfold test performance
  cvPerformance  = rbind(cvPerformance,modelPerf)
}
## performance
cvPerformance %>% 
  summarise_all(.funs = mean)
cvPerformance$Labels = 'True'
write_csv(x = cvPerformance,path = "Output/perturbation_ALR_cvTrainPerformance.csv")
#=============================
## Null Performance
#=============================
nrepeats =  100
nfolds = 5
permuteLabels = T
nullPerformance = data.frame()
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
      
      #RFE Feature Selection
      # fs = rfeSelection.ByMetric(train_ratio = trainData,test_ratio = testData,metric = "multiClass.ROC",ytrain = ytrain,
      #                       sets = 10,impMeasure = "impurity_corrected",nreps = 5,minPercentFeatReturn = .1)
      # trainData = fs$reducedTrainRatios
      # testData = fs$reducedTestRatio
      
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
                     #method = "svmRadial",
                     method = "rf",importance = T,proximity = F,ntree = 5000,tuneGrid = expand.grid(mtry = c(2,sqrt(ncol(trainData)))),
                     metric = "Mean_Balanced_Accuracy",
                     trControl = train_control
      )
      
      #Balanced ACC
      preds = predict.train(glm.mdl1,testData)
      AC = confusionMatrix(preds,ytest)
      AC = data.frame(Group =  rownames(AC$byClass), AC$byClass)
      AC$fold = f
      AC$seed = s
      meanBalancedACC = mean(AC$Balanced.Accuracy)
      ## Mulit ROC
      preds = predict.train(glm.mdl1,testData,type = "prob")
      rocobj = pROC::multiclass.roc(ytest,preds)
      multiClass.AUC = pROC::auc(rocobj)
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
    }))
    
    message("seed-",s,", fold-",f)
    ## Performance
    data.frame(fold = f,seed = s,
               multiCLassAUC =multiClass.AUC,
               meanBalancedACC =meanBalancedACC,
               Train_procrustedCorr = procrCorrTrain,
               Test_procrustedCorr = procrCorrTest,
               numFeatures = ncol(trainData),
               trainAUC = getTrainPerf(glm.mdl1)[2],glm.mdl1$bestTune)
    
  }
  #write sth repeat kfold test performance
  nullPerformance  = rbind(nullPerformance,modelPerf)
}
## performance
nullPerformance %>% 
  summarise_all(.funs = mean)
nullPerformance$Labels = 'Permuted'
write_csv(x = nullPerformance,path = "Output/perturbation_ALR_nullTrainPerformance.csv")
#Visualize
permutationPerf = nullPerformance %>% 
  group_by(seed,Labels) %>% 
  summarise_all(.funs = mean)

#permutation pval
pval = (sum(mean(cvPerformance$multiCLassAUC)<permutationPerf$multiCLassAUC)+1) / (nrow(permutationPerf)+1)
tiff(filename = "Figures/Supplemental/PertubationMicrobiome_permutationTesting.tiff",width = 4,height = 4,units = "in",res = 300)
ggplot(permutationPerf,aes(multiCLassAUC))+
  geom_histogram(aes(fill = Labels),col = "white")+
  scale_fill_brewer(palette = "Dark2")+
  geom_vline(aes(xintercept = mean(cvPerformance$multiCLassAUC),col = "Classification  multi-class AUROC"),
             lty = "dashed",size = 1,
             show.legend = T)+
  geom_vline(aes(xintercept = mean(nullPerformance$multiCLassAUC),col = "Chance  multi-class AUROC"),
             lty = "dashed",size = 1,
             show.legend = T)+
  annotate(geom="label", 
           x=.8, y=20, 
           label=paste("AUROC = ",round(mean(cvPerformance$multiCLassAUC),3)),
           color="black",fill = "white")+
  scale_color_manual(values = c("Black","Red"))+
  scale_fill_manual(values = "purple4" )+
  theme_classic()+
  xlab("5-fold Cross Validated  multi-class AUROC")+
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
        plot.caption = element_text(size = 8,face = "italic"))+
  labs(caption =  paste("Permutation Test(repeats=100), pval =",round(pval,5)))
dev.off()
#####################################################




#######################################################
### Identify Min Signature ############
#######################################################
targetVariable = as.factor(dat[,1])
ph.alr = data.frame(ph.alr)
set.seed(08272008)
## RFE
fs = rfeSelection.ByMetric(train_ratio = ph.alr,test_ratio = ph.alr,
                           metric = "multiClass.ROC",
                           ytrain = targetVariable,
                           sets = 10,
                           impMeasure = "impurity_corrected",
                           nreps = 10,
                           minPercentFeatReturn = .1)
alr.df = fs$reducedTrainRatios
## final model
glm.mdl1=train(x = ph.alr,
               y =targetVariable,
               method = "rf",importance = T,proximity = F,ntree = 5000,tuneGrid = expand.grid(mtry = c(2,sqrt(ncol(ph.alr)))),
               metric = "AUC",
               trControl = train_control
)

rfModel = glm.mdl1$finalModel
randomForest::importance(rfModel)
#var imp
vi  =varImp(glm.mdl1)
vi = vi$importance
vi = data.frame(ID = rownames(vi),vi)
vi = separate(vi,col = 1,into = c("ID","denom"),sep = "\\.")
vi.gather = gather(vi,"Group","Imp",3:ncol(vi))
vi.means = data.frame(Ratio = rownames(vi),meanImp = rowMeans(vi[,-2:-1]))
#######################################################




##########################################################-
### Statistical Analysis of best alr transformed data ####
##########################################################-
#permanova
ph.d = dist(alr.df)
labels = dat[,1]
a.df = data.frame(Type = targetVariable)
pmv = adonis2(ph.d~Type,data = a.df,permutations = 5000)
(pmv)
#dispersion test
mod = betadisper(as.dist(ph.d),group = labels)
anova(mod)
mod.HSD <- TukeyHSD(mod)
plot(mod,label = F)
##########################################################-



##########################################################-
### Principle Coordinate Analysis with using ALR with pseudo F maximizing denom ####
##########################################################-
mds = cmdscale(ph.d,eig = T)
mds.df = data.frame(Strain_Treatment = labels,mds$points)
#Compute PCoA variation
pcoa_var = mds$eig[mds$eig>0]
pcoa_var = (pcoa_var / sum(pcoa_var))*100
#compute centroids
centroid = mds.df %>% 
  gather(key = "Dim",value = "Val",2:ncol(mds.df)) %>% 
  group_by(Strain_Treatment,Dim) %>% 
  summarise(mn = mean(Val)) %>% 
  spread(key = "Dim",value = "mn") %>% 
  mutate(Strain = if_else(str_detect(Strain_Treatment,"CC027"),"CC027","C3H"))
disp = left_join(mds.df,centroid,by = 'Strain_Treatment')
disp = separate(disp,col = 1,into = c("Strain","Sensistization"),sep = "_",remove = F)
mds.df = separate(mds.df,col = 1,into = c("Strain","Sensistization"),sep = "_",remove = F)
unqLabels = unique(labels)[c(1,6,2,7,3,8,4,9,5,10)]
## Relevel Factor
mds.df$Strain_Treatment = factor(mds.df$Strain_Treatment,levels = unqLabels)
disp$Strain_Treatment = factor(disp$Strain_Treatment,levels = unqLabels)
centroid$Strain_Treatment = factor(centroid$Strain_Treatment,levels = unqLabels)
centroid = separate(centroid,col = 1,into = c("Strain","Sensistization"),sep = "_",remove = F)
mds.df$Sensistization = factor(mds.df$Sensistization,level = c("Peanut","Walnut","Milk","Egg","PBS"))
centroid$Sensistization = factor(centroid$Sensistization,level = c("Peanut","Walnut","Milk","Egg","PBS"))
disp$Sensistization = factor(disp$Sensistization,level = c("Peanut","Walnut","Milk","Egg","PBS"))

#PCoA Plot
tiff(filename = "Figures/Main/pre-postPertubationMicrobiome_withALR-PCoA.tiff",width = 5,height = 3,units = "in",res = 300)
ggplot(data = mds.df,aes(X1,X2,col = Sensistization,fill = Sensistization))+
  #geom_hline(yintercept = 0,alpha = .5)+
  #geom_vline(xintercept = 0,alpha = .5)+
  xlab(paste("PCo1 (",round(pcoa_var[1],digits = 2),"% of total variation)",sep = ""))+
  ylab(paste("PCo2 (",round(pcoa_var[2],digits = 2),"% of total variation)",sep = ""))+
  geom_point(aes(shape = Strain),alpha = 1,size = 1)+
  geom_segment(aes(x = X1.y , y = X2.y , xend = X1.x , yend = X2.x ,linetype =  Strain, colour = Sensistization),
               data = disp,alpha = .6,size = .75)+
  geom_point(data = centroid,aes(X1,X2,shape = Strain),size = 3)+
  scale_color_jco()+#(palette = "RdYlBu",direction = 1)+
  scale_fill_jco()+#(palette = "RdYlBu")+
  theme_bw()+
  labs(caption = paste('PERMANOVA F=',round(pmv$F,5),", p=",round(pmv$`Pr(>F)`,5),sep = ""))+
  theme(axis.text = element_text(hjust = 1,size = 7),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.text.x = element_text(size = 7),
        #legend.position = c(.9,.2),
        legend.key.size = unit(.1,units = "in"),
        plot.caption = element_text(size = 6,face = "italic"))+
  theme(panel.grid = element_blank())
dev.off()
##########################################################-



##########################################################-
###  Beta Disp ####
##########################################################-
interactions = factor(labels,levels = unqLabels)
mod = betadisper(as.dist(ph.d),group = interactions)
anv = anova(mod)
permutest(mod)
mod.HSD <- TukeyHSD(mod)
dis = data.frame(Treatment = interactions,Dist = mod$distances)
dis = separate(dis,col = 1,into = c("Strain","Sensistization"),sep = "_",remove = F)
tiff(filename = "Figures/Main/pre-postPertubationMicrobiome_betadisp.tiff",width = 6,height = 5,units = "in",res = 300)
ggplot(dis,aes(Treatment,Dist,fill = Sensistization))+
  #stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha =.5,outlier.shape = NA)+
  theme_classic()+
  geom_jitter(width = .15,pch = 21,size = 1,color = "black")+
  scale_fill_jco()+
  labs(caption = paste('Betadispersion F=',round(anv$`F value`[1],5),", p=",round(anv$`Pr(>F)`[1],5),sep = ""))+
  theme(axis.text.x = element_text(hjust = 1,size = 7,angle = 45),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = "bold"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        legend.key.size = unit(.1,units = "in"),
        plot.caption = element_text(size = 6,face = "italic"))+
  ylab("Distance to Centroid")
dev.off()
##########################################################-



##########################################################-
### Post hoc  #####
##########################################################-
interaction.df = data.frame(Interaction = dat[,1],alr.df)
interactions = as.character(dat[,1])
classes = as.character(unique(interactions))
c = combinat::combn2(classes)
pw_diffMean = foreach(i = 1:nrow(c),.combine = rbind,.packages = "vegan")%do%{
  classes = c[i,]
  ph = interaction.df %>% 
    filter(Interaction %in% classes)
  phd = parallelDist::parDist(as.matrix(ph[,-1]))
  a.df = data.frame(Type = ph$Interaction)
  pmv = adonis2(phd~.,data = a.df,permutations = 100000)
  message(i)
  data.frame(T1 = classes[1],T2 = classes[2],
             pseudoF = pmv$F[1],pval = pmv$`Pr(>F)`[1])
}
pw_diffMean$False_Discovery_Rate = p.adjust(pw_diffMean$pval,method = "BH")
pw_diffMean = pw_diffMean %>% 
  arrange(desc(-pval))
treatment_byStrain = pw_diffMean %>% 
  mutate(keep = if_else(str_split(T1,"_",simplify = T)[,2]==str_split(T2,"_",simplify = T)[,2],1,0)) %>% 
  filter(keep==1) %>% 
  dplyr::select(-keep) %>% 
  arrange(desc(-pval))
write_csv(pw_diffMean,path = "Output/supp1-perturbationMulComp.csv")

#spread pairwise
spread.pw = pw_diffMean[,c(1,2,5)] %>% 
  spread("T2",value = "False_Discovery_Rate")
classes.df = data.frame(T1 = classes)
spread.pw = right_join(spread.pw,classes.df) %>% 
    dplyr::select(T1,one_of(classes))
#Visulaize differences
el = pw_diffMean %>% 
  filter(False_Discovery_Rate<=0.05)
g = graph_from_edgelist(el = as.matrix(el[,1:2]),directed = F)
E(g)$weight = 1/el$pseudoF^2
plot(g,layout = layout.fruchterman.reingold,vertex.size = strength(graph = g,weights = 1/(E(g)$weight)^.5 ))
##########################################################-




##########################################################-
### mean LR analysis ########
##########################################################-
#average results by strain-sensitization
interaction.df = data.frame(Interaction = dat[,1],alr.df)
gather.alr = interaction.df %>% 
  gather(key = "Ratio",value = "value",2:ncol(alr.df)) %>% 
  group_by(Interaction,Ratio) %>% 
  summarise(meanLR = mean(value)) %>% 
  spread(key = "Ratio",value = "meanLR") 
gather.alr = data.frame(gather.alr)
##########################################################-




##########################################################-
# perform cluster analysis ##############
##########################################################-
#React Information
sampleNames = gather.alr[,1]
react = data.frame( Interaction = dat[,1],
                    Strain = df.metadata$Strain,
                    React = df.metadata$Reactors,
                    Sens = df.metadata$Sensitization,
                    IgE = df.metadata$IgE
)

#Graph
#Disim Based on Pseudo F
spread.pw = pw_diffMean %>% 
  dplyr::select(T1,T2,pseudoF,False_Discovery_Rate) %>% 
  mutate(weight = if_else(False_Discovery_Rate>0.05,pseudoF,0))
el = spread.pw
g = graph_from_edgelist(as.matrix(el[,1:2]),directed = F)
E(g)$weight = el$weight
g = igraph::simplify(g,remove.multiple = T,remove.loops = T)
adMar = get.adjacency(g,attr = "weight",sparse = F)
#remove egg distconnected component
adMar = adMar
g = graph_from_adjacency_matrix(adMar,mode = "undirected",weighted = T)


#Inverse Weighting
E(g)$weight = (1/E(g)$weight^2)

#Extract KNN Adj Matrix
adjacency_matrix <- igraph::as_adjacency_matrix(graph = g,sparse = F,
                                                attr = "weight"
)

## Maximize React Entropy
# Resolution Parm
sampleNames = rownames(adjacency_matrix)
set.seed(08272008)
x = seq(.5,3,length.out = 100)
#React Cluster Entropy
entropy_ = c()
numComm = c()
j = 1
for(i in x){
  partition <- leiden::leiden(adjacency_matrix,resolution_parameter = i)
  dust = data.frame(Interaction = sampleNames,Memb = partition)
  dust = separate(data = dust,col = 1,into = c("strain","treatment"))
  
  #React Entropy
  c = data.frame(Interaction = rownames(adjacency_matrix),Cluster = partition)
  ncolors = max(c$Cluster)
  c$Cluster = as.factor(c$Cluster)
  clus.metadata = left_join(react, c)
  p = clus.metadata %>% 
    group_by(Cluster,React) %>% 
    summarise(n = n()) %>% 
    spread(key = "React",value = "n",fill = 0) %>% 
    mutate(reactPercent = reactor / (reactor + non))
  numComm[j] = n_distinct(partition)
  entropy_[j] = vegan::diversity(p$reactor)
  
  j = j+1
}

minEntrop = min(entropy_)
results = data.frame(ResParm = x, Entropy = entropy_) %>% 
  filter(Entropy==minEntrop) %>% 
  arrange(desc(ResParm))
resParm = results[1,1]

#Visualize
## Num Communities
df = data.frame(resolutionParameter = x,numberCommunities = numComm)
size = if_else(x==resParm,3,1)
color = if_else(x==resParm,"Red","Black")
tiff(filename = "Figures/Supplemental/pre-postPertubationMicrobiome_withALR-NumCommunities.tiff",width =4 ,height = 3,units = "in",res = 300)
ggplot(df,aes(resolutionParameter,numberCommunities))+
  geom_line()+
  theme_classic()+
  #annotate(x = 2.7,y = .35,geom = "text",label = paste0("(",round(resParm,2),",",round(minEntrop,2),")"))+
  geom_point(col = color,size = size)+
  ylab("Number of Communities")+
  xlab("Resolution Parameter")+
  theme(axis.text.x = element_text(hjust = 1,size = 7,angle = 45),
        #axis.title.x = element_blank(),
        axis.title = element_text(size = 8,face = "bold"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        legend.key.size = unit(.1,units = "in"),
        plot.caption = element_text(size = 6,face = "italic"))
dev.off()
## Entropy
df = data.frame(resolutionParameter = x,reactEntropp = entropy_)
size = if_else(x==resParm,3,1)
color = if_else(x==resParm,"Red","Black")
tiff(filename = "Figures/Supplemental/pre-postPertubationMicrobiome_withALR-ResParmMixEntropyMin.tiff",width =4 ,height = 3,units = "in",res = 300)
ggplot(df,aes(resolutionParameter,reactEntropp))+
  geom_line()+
  theme_classic()+
  annotate(x = 2.7,y = .35,geom = "text",label = paste0("(",round(resParm,2),",",round(minEntrop,2),")"))+
  geom_point(col = color,size = size)+
  ylab("Entropy_Reactor | Clustering(resolutionParm)")+
  xlab("Resolution Parameter")+
  theme(axis.text.x = element_text(hjust = 1,size = 7,angle = 45),
        #axis.title.x = element_blank(),
        axis.title = element_text(size = 8,face = "bold"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        legend.key.size = unit(.1,units = "in"),
        plot.caption = element_text(size = 6,face = "italic"))
dev.off()

#Cluster given tune resolution parameter
cluster <- leiden::leiden(adjacency_matrix,resolution_parameter = resParm) 
c = data.frame(Interaction = rownames(adjacency_matrix),Cluster = cluster)
ncolors = n_distinct(c$Cluster)
c$Cluster = as.factor(c$Cluster)
c = c %>% 
  arrange(desc(Interaction)) %>% 
  arrange(desc(Cluster))
colors = ggsci::pal_npg()(ncolors)
color.df = data.frame(Color = colors,
                      Cluster = as.character(unique(c$Cluster)))
color.df = left_join(c,color.df)
color.df = color.df %>% 
  filter(Interaction %in%  gather.alr[,1])
vcolor = data.frame(Interaction = V(g)$name)
vcolor = left_join(vcolor,color.df)
#vis graph
tiff(filename = "Figures/Supplemental/pre-postPertubationMicrobiome_withALR-clusterLeiden.tiff",width = 6,height = 6,units = "in",res = 300)
plot(g,vertex.color = vcolor$Color,
     edge.width=2.5*E(g)$weight,
    layout = layout_with_fr , edge.curved = .25,
    vertex.label.cex=0.7, vertex.label.dist=2.5,
     #vertex.shape = if_else(df.metadata$Strain=="C3H","square","circle"),
     vertex.size = 15,margin = c(0,0,0,0),rescale = T
     )
dev.off()

#Visualize
clus.metadata = left_join(react, c)
clus.metadata$Cluster
p = clus.metadata %>% 
  group_by(Cluster,React) %>% 
  summarise(n = n()) %>% 
  spread(key = "React",value = "n",fill = 0) %>% 
  mutate(reactPercent = reactor / (reactor + non))
testStat = vegan::diversity(p$reactor)

tiff(filename = "Figures/Supplemental/pre-postPertubationMicrobiome_withALR-clusterReactor.tiff",width = 4,height = 4,units = "in",res = 300)
ggplot(p,aes(Cluster,reactor))+
  geom_bar(stat = "identity",col = "black",fill = colors[ncolors:1])+
  coord_flip()+
  theme_classic()+
  ylab("Number of Anaphylactic Mice")+
  labs(caption = paste("Entropy = ",round(testStat,5)))+
  theme(text = element_text(size = 8))+
  theme(axis.text.x = element_text(hjust = 1,size = 7,angle = 45),
        #axis.title.x = element_blank(),
        axis.title = element_text(size = 8,face = "bold"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        legend.key.size = unit(.1,units = "in"),
        plot.caption = element_text(size = 6,face = "italic"))
dev.off()

entropySim = function(){
  c$Cluster = sample(c$Cluster)
  clus.metadata = left_join(react, c)
  clus.metadata$Cluster
  p = clus.metadata %>% 
    group_by(Cluster,React) %>% 
    summarise(n = n()) %>% 
    spread(key = "React",value = "n",fill = 0) %>% 
    mutate(reactPercent = reactor / (reactor + non))
  vegan::diversity(p$reactor)
}
#Simulate null entropy distribution
set.seed(08272008)
ntrials = 5000
reactEntropy.null = replicate(ntrials,entropySim(),simplify = T)
hist(reactEntropy.null)
pval.entropy = (sum(testStat>reactEntropy.null)+1)/ (ntrials+1)
permEntropy = data.frame(Pemuted_Entropy = reactEntropy.null)
tiff(filename = "Figures/Supplemental/pre-postPertubationMicrobiome_withALR-reactClusterentropy.tiff",width =4,height = 4,units = "in",res = 300)
ggplot(permEntropy,aes(Pemuted_Entropy))+
  geom_histogram(bins = 10,col = "white")+
  xlab("Permuted Entropy (n = 5000)")+
  theme_classic()+
  geom_vline(xintercept = testStat,col = "red",lty = "dashed")+
  annotate(geom = "label",x = .35,y=250,size = 2,label = paste("Test Entropy\n",round(testStat,3),sep = ""))+
  labs(caption = paste("p = ",round(pval.entropy,digits = 5),sep = ""),title = "Between Cluster 'React' Entropy ")+
  theme(text = element_text(size = 8))
dev.off()
##########################################################-





##########################################################-
###  Heatmap PreProcess #####
##########################################################-
sc = color.df$Color
# log ratio matrix
mat = left_join(data.frame(Interaction=c$Interaction),gather.alr)
mat = as.matrix(mat[,-1])
rownames(mat) = c$Interaction
colnames(mat) = str_replace(string = colnames(mat),pattern = ".V79",replacement = "")
col <- colorRampPalette(RColorBrewer::brewer.pal(10, "PuOr"))(1000)
mat = t(mat)
# heatmap microbiome composition
tiff(filename = "Figures/Supplemental/pre-postPertubationMicrobiome_withALR-interactionHeatmap.tiff",width = 8,height = 5.7,units = "in",res = 300)
gplots::heatmap.2( mat,
                   col = col,#gplots::redblue(n = 1000), #viridis::viridis(n = 1000,direction = 1,option = "C"),#gplots::redblue(n = 1000) ,#viridis::viridis(n = 1000,option = "D"),#,
                   Rowv = TRUE,Colv = F,
                   margins = c(6,4), 
                   sepwidth=c(0.01,0.01),
                   sepcolor="black",colsep=1:ncol(mat),rowsep=1:nrow(mat),
                   hclustfun = function(x) hclust(x, method = "complete"),
                   distfun = function(x) parallelDist::parallelDist(x,method = "euclidean"),
                   key=TRUE, symkey=TRUE, 
                   scale = "none",
                   density.info="density", trace="none",
                   ColSideColors = sc,
                   #main = "Presensitization Gut Mircobiome",
                   key.title = "Log-ratio",key.xlab = paste("log( v_i / v_",reference$i,")",sep = "" ),
                   #RowSideColors = cc,
                   cexRow = .75,cexCol = .75,
                   #adjCol = c(1,1)
                   
)
dev.off()
######################




########################################-
### PCA Analysis ####
########################################-
#Contribution PCA Biplot (adapted from Greenacre easyCODA package)
obj = easyCODA::PCA(data =gather.alr[,-1],weight = F, )
dim = c(1, 2)
axes.inv = c(1, 1)
obj.rpc <- obj$rowcoord[, dim] %*% diag(obj$sv[dim] * axes.inv)
obj.csc <- obj$colcoord[, dim] %*% diag(axes.inv)
obj.cpc <- obj.csc %*% diag(obj$sv[dim])
obj.ccc <- obj.csc * sqrt(obj$colmass)
obj.crd <- obj.ccc
rescale = 1
cex = c(0.8, 0.8)

#Gathered PCA
gather.alr = left_join(c,gather.alr)[,-2]
pc = prcomp(gather.alr[,-1],center = T)
pc_var = (pc$sdev/sum(pc$sdev))*100
pc_loading = data.frame(Ratio = rownames(pc$rotation),0.95 * rescale * obj.crd) #contribution biplot weightings
loading.df = data.frame(Ratio = rownames(pc$rotation),
                        x.start = 0, 
                        x.end = pc_loading[,2],
                        y.start = 0, 
                        y.end = pc_loading[,3])
loading.df$Ratio = str_replace_all(string = loading.df$Ratio,pattern = "\\.",replacement = " / ") 
ismap.df = data.frame(Group = c$Interaction,pc$x,Cluster = c$Cluster)
#convex hull
hull_cyl <- ismap.df %>%
  group_by(Cluster) %>%
  dplyr::slice(chull(PC1, PC2))
scalingFactor = 7

tiff(filename = "Figures/Main/pre-postPertubationMicrobiome_withALR-interactionPCA.tiff",width = 4.5,height = 3,units = "in",res = 300)
ggplot(loading.df)+
  geom_segment(data = loading.df,aes(x = x.start,y = y.start,xend = scalingFactor*x.end, yend = scalingFactor*y.end),
               alpha = .35,col = "darkblue",size = .5,
               arrow = arrow(length = unit(0.03, "npc")))+
  geom_text(aes(label=Ratio,x = scalingFactor*x.end,y = scalingFactor*y.end), 
            size=2.25,
            check_overlap = T)+
  geom_text(data = ismap.df,
            aes(PC1,PC2,col = Cluster,fill = Cluster,label = Group),
            nudge_y = -.25,size = 2.25)+
  geom_point(data =ismap.df,aes(PC1,PC2,col = Cluster,fill = Cluster),
             size = 2.5,pch = 21,col = "black")+
  geom_polygon(data = hull_cyl, aes(PC1,PC2,col = Cluster,fill = Cluster),alpha = 0.15)+
  theme_bw()+
  scale_x_continuous(limits = c(-8,10))+
  scale_fill_manual(values = colors[ncolors:1])+
  scale_color_manual(values = colors[ncolors:1])+
  geom_vline(xintercept = 0,size = .25)+
  geom_hline(yintercept = 0,size = .25)+
  xlab(paste("PC1 (",round(pc_var[1],digits = 2),"% of total variation)",sep = ""))+
  ylab(paste("PC2 (",round(pc_var[2],digits = 2),"% of total variation)",sep = ""))+
  theme(axis.text = element_text(hjust = 1,size = 8),panel.grid = element_blank(),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        legend.text = element_text(size = 8),
        #legend.direction = "horizontal",
        legend.title = element_text(size = 8),
        #legend.box.background = element_rect(colour = "black"),legend.key = element_rect(colour = "black"),
        legend.position = "right",legend.key.size = unit(.15,units = "in"),
        #legend.background = element_rect(colour = "black"),
        axis.text.x = element_text(size = 8),plot.caption = element_text(size = 8))
dev.off()





#####################################-
### Kruskall Test ####
#####################################-
interaction.df = data.frame(Interaction = dat[,1],alr.df)
mat = left_join(c,interaction.df)
gather.alr2 = mat %>%
  gather(key = "Ratio",value = "value",3:ncol(mat)) %>%
  separate(3,into = c("ID","denom"),sep = "\\.") %>%
  group_by(Cluster,ID) %>%
  summarise(meanLR = mean(value),n = n(),sd = sd(value),CI = (1-.96*sd)/sqrt(n),lb = meanLR - CI, ub = meanLR + CI)
gather.alr2 = right_join(cnames,gather.alr2)
gather.alr2$Family_Genus = paste(gather.alr2$Order,gather.alr2$Family,gather.alr2$Genus,gather.alr2$Species,sep = " ")
unqTaxa = unique(gather.alr2$Family_Genus) # from presen script


#Krus test
interaction.df = data.frame(Interaction = dat[,1],alr.df)
mat = left_join(c,interaction.df)
gather.alr1 = mat %>%
  gather(key = "Ratio",value = "value",3:ncol(mat)) %>%
  separate(3,into = c("ID","denom"),sep = "\\.") 
gather.alr1 = right_join(cnames,gather.alr1)
gather.alr1$Family_Genus = paste(gather.alr1$Order,gather.alr1$Family,gather.alr1$Genus,gather.alr1$Species,sep = " ")
krus.df = data.frame()
for( i in 1:length(unqTaxa)){
  ph = gather.alr1 %>% 
    filter(Family_Genus==unqTaxa[i])
  #Anova
  md = lm(value~Cluster,data = ph)
  anv =anova(md)
  #Krus
  kt = kruskal.test(x = ph$value,g = ph$Cluster)
  ph = data.frame(Taxa = unqTaxa[i],pval.kruskal = kt$p.value,anova.pval = anv$`Pr(>F)`[1])
  krus.df = rbind(krus.df,ph)
}
krus.df$FDR = p.adjust(p = krus.df$pval,method = "BH")
krus.df$FDR.anova = p.adjust(p = krus.df$anova.pval,method = "BH")
krus.df_ = krus.df %>% 
  filter(Taxa!="__ __ __ __")
### Supplement 3
write_csv(x = krus.df_,path = "Output/supp3-perturbationKruskallTest.csv")



krus.df.final = krus.df %>% 
  filter(FDR<0.05)

## GG blot
df = gather.alr1 %>% 
  filter(Family_Genus %in% krus.df.final$Taxa)
stat.test <- df %>%
  group_by(Family_Genus) %>%
  rstatix::wilcox_test(value ~ Cluster) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj")
bp <- ggpubr::ggbarplot(
  df, x = "Family_Genus", y = "value", fill = "Cluster",
  palette = "npg", add = "mean_ci",
  position = position_dodge(0.8),ggtheme = theme(axis.text.x = element_text(size = 8,angle = 90,hjust = 1))
)
stat.test <- stat.test %>%
  rstatix::add_xy_position(x = "Family_Genus", fun = "mean_ci", dodge = 0.8)
bp + 
  ggpubr::stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0.00,
    bracket.nudge.y = -2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))



#order by cluster
df = gather.alr2 %>% 
  filter(Family_Genus %in% krus.df.final$Taxa)
orderCluster = 4
clusN.alr = df %>% 
  filter(Cluster==orderCluster) %>% 
  filter(Family_Genus!="__ __ __ __") %>% 
  arrange(desc(meanLR))
#reorder data
df$Family_Genus = factor(df$Family_Genus,levels = clusN.alr$Family_Genus)
refFamGen = paste(reference$Order,reference$Family,reference$Genus,reference$Species,sep = " ")
df$Cluster = paste("Cluster",df$Cluster)
df$Family_Genus = str_replace_all(string = df$Family_Genus,pattern = "o__",replacement = "")
df$Family_Genus = str_replace_all(string = df$Family_Genus,pattern = "f__",replacement = "")
df$Family_Genus = str_replace_all(string = df$Family_Genus,pattern = "g__",replacement = "")
df$Family_Genus = str_replace_all(string = df$Family_Genus,pattern = "s__",replacement = "")
df$Family_Genus = str_replace_all(string = df$Family_Genus,pattern = "__",replacement = "")
df =df %>% 
  filter(Family_Genus!="")
df$Family_Genus = paste(df$Family_Genus," (",df$ID,")",sep = "")

#visualize pertubation signature
p = ggplot(df,aes(Family_Genus,meanLR,fill = meanLR))+
  geom_col(col = "black",width = .6,size = .5)+
  geom_errorbar(aes(ymin = lb , ymax = ub ), width=0.2,size = .5)+
  facet_grid(.~Cluster)+
  coord_flip()+
  scale_fill_distiller(palette = "RdBu",direction = 1)+
  scale_color_brewer(palette = "Set3")+
  scale_y_continuous(breaks = c(-2.5,0,2.5))+
  #ggridges::theme_ridges()+
  theme_minimal()+
  geom_hline(yintercept = 0,col = "black",size = .5)+
  labs(caption = paste("Top ALR Denominator = ",refFamGen," (",reference$ID,")",sep = ""))+
  theme(legend.position = "none",text = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        strip.background.x =  element_rect(fill = colors),
        #strip.text.x = element_text(colour = colors),
        axis.line.x = element_line(colour = "black"),
        strip.text = element_text(size = 8,face = "bold",colour = "white"),
        axis.text.x = element_text(size = 8),
        plot.caption = element_text(size = 6,face = "italic"),
        axis.text.y = element_text(size = 8),
        axis.ticks.x = element_line(),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.y = element_blank())+
  ylab("Mean logratio")


g <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
fills <- colors[ncolors:1]
k <- 1
for (i in strip_t) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

tiff(filename = "Figures/Main/pre-postPertubationMicrobiome_withALR-ratioSignature.tiff",width = 10,height = 3.75,units = "in",res = 300)
grid::grid.draw(g)
dev.off()







#Variable Importance
vi.gather = right_join(vi.gather,cnames)
vi.gather$Family_Genus = paste(vi.gather$Order,vi.gather$Family,vi.gather$Genus,vi.gather$Species,sep = " ")
vi.gather = vi.gather %>% 
  filter(Family_Genus %in% clusN.alr$Family_Genus)
vi.gather$Family_Genus = str_replace_all(string = vi.gather$Family_Genus,pattern = "o__",replacement = "")
vi.gather$Family_Genus = str_replace_all(string = vi.gather$Family_Genus,pattern = "f__",replacement = "")
vi.gather$Family_Genus = str_replace_all(string = vi.gather$Family_Genus,pattern = "g__",replacement = "")
vi.gather$Family_Genus = str_replace_all(string = vi.gather$Family_Genus,pattern = "s__",replacement = "")
vi.gather$Family_Genus = str_replace_all(string = vi.gather$Family_Genus,pattern = "__",replacement = "")

#visualize pertubation signature
tiff(filename = "Figures/Supplemental/pre-postPertubationMicrobiome_withALR-FeatureImportance.tiff",width = 6,height = 6.5 ,units = "in",res = 300)
ggplot(vi.gather,aes(Group,Imp,fill = Group))+
  geom_bar(stat = "identity",col = "black")+
  facet_wrap(.~Family_Genus,labeller = label_wrap_gen(width=10))+
  scale_fill_jco()+
  ylab("Random Forest Feature importance")+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none",axis.title.y = element_blank(),
        strip.text = element_text(size = 6),axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))
dev.off()
