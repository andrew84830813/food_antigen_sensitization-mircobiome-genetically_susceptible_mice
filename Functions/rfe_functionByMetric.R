#######################################################################
########     <<<Random Forest Recursive Selection Function>>>    ########
#######################################################################
###   train_ratio = logratio matrix (samples x features/taxa) to be reduced
###   test_ratio = logratio matrix (samples x features/taxa) to be calulated for testing outside of fold
###   impMeasure = random forest metric to rank features
###   sets = number of rounds/interals of feature selection
###   ntrees = number of trees to grow at each round 
###   minPercentFeatReturn = 0.3
###   nreps = number of iterations of training for each round/interval of features selection
#######################################################################



rfeSelection.ByMetric = function(train_ratio,test_ratio,ytrain,
                                 #defaultParms =defaultParms,useXGB = F,
                                 metric = "AUROC",
                                 impMeasure = "impurity_corrected",
                                 sets = 10,ntrees = 2000,mtry_ = 2,
                                 minPercentFeatReturn = .3,nreps = 5){ 
          
  
 
  #test matrix to select reduced features from
  train_ratio = data.frame(train_ratio)
  baseRatios = train_ratio
  xtst.pc = data.frame(test_ratio)
  
  #AUPC
  minorityClass = names(which.min(table(ytrain)))
  majorityClass = unique(ytrain[ytrain!=minorityClass])
  
  
  #prevent from training with 0 features 
  subsets = round(seq(round(ncol(train_ratio)*minPercentFeatReturn),ncol(train_ratio),length.out = sets))
  if(min(subsets)<=1){
    pos = which(subsets<=1)
    subsets = subsets[-pos]
  }
  
  sets = length(subsets)
  subsets = subsets[sets:1]
  subsetAUC = data.frame()
  importanceList = list()
  
  #random forest recirsive feature selection 
  importance.df =c()
  glm.trainAUC = c()
  glm.auprc = c()
  multiClass.roc = c()
  rfeCV = foreach(r=1:nreps,.combine = rbind ,
                  .packages=c("ranger","pROC","tidyverse","caret","compositions")) %do% {
                    
                      n=mtry_
                      ph.df = data.frame(Status = factor(ytrain),train_ratio)
               
                      rf.ranger = ranger::ranger(formula = Status~.,data = ph.df,
                                                 num.trees = ntrees,mtry = n,
                                                 max.depth = 0,
                                                 min.node.size = 1,
                                                 write.forest = T,probability = T,
                                                 #splitrule = "extratrees",
                                                 sample.fraction = 1,
                                                 importance = impMeasure,
                                                 scale.permutation.importance = T)
                      message(r)
                      
                      #AUC
                      glm.trainAUC[r] = pROC::auc(ytrain,rf.ranger$predictions[,1])
                     
                      #AUPRC
                      posClass = rf.ranger$predictions[,1]
                      negClass = rf.ranger$predictions[,1]
                      posClass = posClass[ytrain==minorityClass]
                      negClass = negClass[ytrain==majorityClass]
                      AUPRC = PRROC::pr.curve(scores.class0 = negClass,scores.class1 = posClass,curve = TRUE)$auc.integral
                      glm.auprc[r] = AUPRC
                      
                      if(metric=="multiClass.ROC"){
                        #multiclass ROC
                        rocobj = pROC::multiclass.roc(ytrain,rf.ranger$predictions)
                        multiClass.roc[r] = pROC::auc(rocobj)  
                      }
                      
                      #output
                      data.frame(Rep = r,Ratio = names(ranger::importance(rf.ranger)),
                                 Imp = ranger::importance(rf.ranger))
                    }
                    
                    
  vim = rfeCV %>% 
    group_by(Ratio) %>% 
    dplyr::summarise(meanImp = mean(Imp))
  
  #select subset ny metric
  if(metric=="AUROC"){
    ph = data.frame(Feats = ncol(train_ratio),AUC = mean(glm.trainAUC),sdAUC = sd(glm.trainAUC))
    subsetAUC = rbind(subsetAUC,ph)
    importanceList[[1]] = vim
  }else if(metric=="AUPRC"){
    ph = data.frame(Feats = ncol(train_ratio),AUC = mean(glm.auprc),sdAUC = sd(glm.auprc))
    subsetAUC = rbind(subsetAUC,ph)
    importanceList[[1]] = vim
  }else if(metric=="multiClass.ROC"){
    ph = data.frame(Feats = ncol(train_ratio),AUC = mean(multiClass.roc),sdAUC = sd(multiClass.roc))
    subsetAUC = rbind(subsetAUC,ph)
    importanceList[[1]] = vim
  }
  
  
  for(i in 2:sets){
    #subset
    sn = top_n(importanceList[[i-1]],n = subsets[i],wt = meanImp)
    train_ratio = data.frame(subset(train_ratio, select=as.character(sn$Ratio)))
    importance.df =c()
    glm.trainAUC = c()
    glm.auprc = c()
    multiClass.roc = c()
    rfeCV = foreach(r=1:nreps,.combine = rbind ,
                    .packages=c("ranger","pROC","tidyverse","caret","compositions")) %do% {
                      
                      n=mtry_
                      ph.df = data.frame(Status = factor(ytrain),train_ratio)
                      
                      rf.ranger = ranger::ranger(formula = Status~.,data = ph.df,
                                                 num.trees = ntrees,mtry = n,
                                                 max.depth = 0,
                                                 min.node.size = 1,
                                                 write.forest = T,probability = T,
                                                 #splitrule = "extratrees",
                                                 sample.fraction = 1,
                                                 importance = impMeasure,
                                                 scale.permutation.importance = T)
                      message(r)
                      
                      #AUC
                      glm.trainAUC[r] = pROC::auc(ytrain,rf.ranger$predictions[,1])
                      
                      #AUPRC
                      posClass = rf.ranger$predictions[,1]
                      negClass = rf.ranger$predictions[,1]
                      posClass = posClass[ytrain==minorityClass]
                      negClass = negClass[ytrain==majorityClass]
                      AUPRC = PRROC::pr.curve(scores.class0 = negClass,scores.class1 = posClass,curve = TRUE)$auc.integral
                      glm.auprc[r] = AUPRC
                      
                      if(metric=="multiClass.ROC"){
                        #multiclass ROC
                        rocobj = pROC::multiclass.roc(ytrain,rf.ranger$predictions)
                        multiClass.roc[r] = pROC::auc(rocobj)  
                      }
                      
                      #output
                      data.frame(Rep = r,Ratio = names(ranger::importance(rf.ranger)),
                                 Imp = ranger::importance(rf.ranger))
                      }
                      
                      
                   
    vim = rfeCV %>% 
      group_by(Ratio) %>% 
      dplyr::summarise(meanImp = mean(Imp))
    
    #select subset ny metric
    if(metric=="AUROC"){
      ph = data.frame(Feats = ncol(train_ratio),AUC = mean(glm.trainAUC),sdAUC = sd(glm.trainAUC))
      subsetAUC = rbind(subsetAUC,ph)
      importanceList[[i]] = vim
    }else if(metric=="AUPRC"){
      ph = data.frame(Feats = ncol(train_ratio),AUC = mean(glm.auprc),sdAUC = sd(glm.auprc))
      subsetAUC = rbind(subsetAUC,ph)
      importanceList[[i]] = vim
    }else if(metric=="multiClass.ROC"){
      ph = data.frame(Feats = ncol(train_ratio),AUC = mean(multiClass.roc),sdAUC = sd(multiClass.roc))
      subsetAUC = rbind(subsetAUC,ph)
      importanceList[[i]] = vim
    }
    
    message(paste("Subset-",subsets[i]," ",i,"/",sets," Calculated............",sep = ""))
  }
  
  #select best features
  i = which.max(subsetAUC$AUC)
  varImp_ = importanceList[[i]]
  train_ratio = data.frame(subset(baseRatios, select=as.character(varImp_$Ratio)))
  xtst.pc1 = data.frame(subset(xtst.pc, select=as.character(varImp_$Ratio)))
  subsetAUC$CI = (1.96*(subsetAUC$sdAUC))/sqrt(nreps)
  subsetAUC$lower = subsetAUC$AUC - subsetAUC$CI
  subsetAUC$upper = subsetAUC$AUC + subsetAUC$CI
  
  
  
  return(list(reducedTrainRatios = train_ratio,
              reducedTestRatio = xtst.pc1,
              trainPerformance = subsetAUC)
  )
}



