train_control <- trainControl(method="repeatedcv",
                              repeats = 5,
                              number = 5,
                              
                              seeds = NULL,
                              classProbs = TRUE,
                              savePredictions = F,
                              allowParallel = TRUE,
                              summaryFunction = caret::multiClassSummary
)




#Functions
kfoldDataPartition = function(df,kfold,seed = 08272008,permuteLabel = F){ 
  set.seed(seed)
  
  if(permuteLabel==T){
    df[,1] = sample(df[,1])
  }
  f = createFolds(df[,1],k = kfold,list = FALSE)
  df$fold = f
  df_ = df
  df_ = data.frame(ID = rownames(df_),df_)
  
  #Partition Data
  foldData = list()
  ph_list = list()
  
  for (j in 1:kfold){
    xtrain = df%>%
      filter(fold != j)%>%
      dplyr::select(-fold)
    
    xtrain_ID = df_%>%
      filter(fold != j)%>%
      dplyr::select(ID,1,fold)
    
    xtest = df%>%
      filter(fold == j)%>%
      dplyr::select(-fold)
    
    xtest_ID = df_%>%
      filter(fold == j)%>%
      dplyr::select(ID,1,fold)
    
    ph_list[[1]] = xtrain
    names(ph_list)[1] = "xtrain_combinedFolds"
    ph_list[[2]] = xtest
    names(ph_list)[2] = "xtest_kthFold"
    ph_list[[3]] = xtrain_ID
    names(ph_list)[3] = "xtrain_IDs"
    ph_list[[4]] = xtest_ID
    names(ph_list)[4] = "xtest_IDs"
    
    foldData[[j]] = ph_list
    names(foldData)[j] = paste("fold_",j,sep = "")
  }
  
  return(foldData)
}



#Functions
catchError <- function(expr){
  tryCatch(expr,
           error = function(e){
             T
           },
           warning = function(w){
             F
           },
           finally = {
             F
           })
  
}



fastImputeZeroes = function(df.cdata2, impFactor = 1e-11){
  cn = colnames(df.cdata2)
  rn = rownames(df.cdata2)
  df.cdata2 = as.matrix(clo(df.cdata2))
  
  if(is.null(impFactor)){
    impFactor =min(df.cdata2[df.cdata2>0])/10  
  }
  
  df_imputed = foreach(i  = 1:nrow(df.cdata2),.combine = rbind)%dopar%{
    sampleData = df.cdata2[i,]
    nz = sum(sampleData==0)
    sampleData[sampleData==0] = impFactor
    sampleData[sampleData!=0] =  sampleData[sampleData!=0]*(1-impFactor*nz)
  }
  
  colnames(df_imputed) = cn
  rownames(df_imputed) = rn
  return(df_imputed)
}





pwDistProcCorrelation = function(baseMatrix,stretchMatrix,weighted = T){
  ### compute clr matrix
  clr.base = data.frame(easyCODA::CLR(data = baseMatrix,weight = weighted)$LR)
  ### distance matrices 
  base.d= parallelDist::parDist(as.matrix(clr.base))
  stretch.d= parallelDist::parDist(as.matrix(stretchMatrix))
  
  #proc. Correlation
  proc = vegan::protest(X = stretch.d,Y = base.d,permutations = 1)
  
  return(proc$t0)
}

