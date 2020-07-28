# Author:   Andrew Hinton
# Email:    andrew84@email.unc.edu
# github:   https://github.com/andrew84830813/food_antigen_sensitization-mircobiome-genetically_susceptible_mice.git



### Read Data
dat = read_csv("Output/preCompositionByStrain.csv")

### start parallel cluster
ncores = detectCores()
clus <- makeCluster(ncores-1) 
registerDoParallel(clus) 


### compute clr matrix
clr_mat = data.frame(clr(dat[,-1]))

### by sample pairwise distance matrix
dat = data.frame(dat)
labels = dat[,1]
d = parallelDist::parDist(as.matrix(clr_mat))

### Calculate Procrustes(Greenacre) and PseduoF for all denoms in Additive Log-ratio (ALR) transform ####
cn = colnames(dat[,-1])
alr.procrustes = foreach(i = 1:ncol(dat[,-1]),.combine = rbind,.packages = c("compositions","vegan"))%dopar%{
  # compute alr transformation using the ith feature as denom
  ph.alr = data.frame(alr(dat[,-1],ivar = i))
  # Compute procrustes correlation between alr and full logratio configurations distance matices 
  proc = vegan::protest(X = ph.alr,Y = d,permutations = 100)
  
  # compuet pseudo f with adonis function
  ph.d = parallelDist::parDist(as.matrix(ph.alr))
  a.df = data.frame(Strain = labels)
  pmv = adonis2(ph.d~Strain,data = a.df,permutations = 5000)
  data.frame(Denom = cn[i],i = i,Correlation = proc$t0,
             pval = proc$signif,
             pseudoF = pmv$F[1],
             permanova_pval = pmv$`Pr(>F)`[1])
}



### Select top AlR denom, that is ALR denom. that has the highest correlation to the full logratio configuration
topALRs = top_n(alr.procrustes,n = 1,wt = pseudoF)
i = topALRs$i[1]
### compute weights (avg proportions) of each taxa 
weights = data.frame(Denom = names(colMeans(dat[,-1])), weight = colMeans(dat[,-1]))

### Apply ALR transformation using the top denom from above
### Transformed data
keepRatios =  data.frame(Strain = dat[,1],alr(dat[,-1],ivar = i))
#procrusted table
alr.procrustes1 = left_join(alr.procrustes,weights)
colnames(alr.procrustes1)[1] = "ID"
alr.procrustes1 = left_join(cnames,alr.procrustes1)


### Procrustes Correlation Visualization
#compute alr transformation using the ith feature as denom
ph.alr = data.frame(alr(dat[,-1],ivar = i))
#compute pairwise sample distance using ith alr tranformation
ph.d = dist(ph.alr)
#Compute procrustes correlation between alr and full logratio configurations distance matices 
proc = vegan::protest(X = ph.d,Y = d,permutations = 100)
dvec = as.vector(as.matrix(d)[lower.tri(as.matrix(d),diag = F)])
ph.vec = as.vector(as.matrix(ph.d)[lower.tri(as.matrix(ph.d),diag = F)])
#Plot Output
tiff(filename = "Figures/Supplemental/preSentization_topDenomProcrutes.tiff",width = 5,height = 5,units = "in",res = 300)
data.frame(d = dvec,alr.dist = ph.vec) %>%
  ggplot(aes(d,alr.dist)) +
  geom_point(alpha = .5)+
  theme_classic()+
  xlab("Logratio distances")+
  ylab(paste("ALR:V",i, " Distances",sep = ""))+
  labs(caption = paste("Procrustes Correlation = ",round(proc$t0,digits = 4)))+
  ggtitle("Pre-sensitization")
dev.off()
#Write Supplemental PseudoF and Procrustes Table
write_csv(x = alr.procrustes1,path = "Output/supplementalTable_preSentitization_procrustesAnalysis.csv")
###############################################################3




### Principle Coordinate Analysis with using ALR with pseudo F maximizing denom ####
mds = cmdscale(ph.d,eig = T)
mds.df = data.frame(Strain = dat[,1],mds$points)
#Compute PCoA variation
pcoa_var = mds$eig[mds$eig>0]
pcoa_var = (pcoa_var / sum(pcoa_var))*100
#Write Figure Output
tiff(filename = "Figures/Supplemental/presensitizationMicrobiome_withALR-PCoA.tiff",width = 9,height = 5,units = "in",res = 300)
ggplot(mds.df,aes(X1,X2,fill = Strain))+
  geom_point(size = 5,pch = 21,alpha = .7)+
  theme_bw()+
  xlab(paste("PCo1 (",round(pcoa_var[1],digits = 2),"% of total variation)"))+
  ylab(paste("PCo2 (",round(pcoa_var[2],digits = 2),"% of total variation)"))+
  scale_fill_manual(values = c("Red","Blue"))
dev.off()



### Presenstization Heatmap PreProcess #####
ncolors = n_distinct(keepRatios$Strain)
colors = c("Red","Blue")
color.df = data.frame(Color = colors,
                      Strain = as.character(unique(keepRatios$Strain)))
color.df = left_join(keepRatios,color.df)
mat = as.matrix((keepRatios[,-1]))
rownames(mat) = keepRatios[,1]
col <- colorRampPalette(RColorBrewer::brewer.pal(10, "PuOr"))(1000)
mat = t(mat)
sc = color.df$Color
# heatmap microbiome composition
tiff(filename = "Figures/Main/presensitizationMicrobiome_withALR-Heatmap.tiff",width = 9,height = 5,units = "in",res = 300)
gplots::heatmap.2( mat,
                   col = col,#viridis::viridis(n = 1000,direction = -1,option = "D"),#gplots::redblue(n = 1000) ,#viridis::viridis(n = 1000,option = "D"),#,
                   Rowv = TRUE,
                   #margins = c(2, 2),
                   hclustfun = function(x) hclust(x, method = "average"),
                   distfun = function(x) parallelDist::parallelDist(x,method = "euclidean"),
                   key=TRUE, symkey=TRUE, 
                   scale = "none",
                   density.info="density", trace="none",
                   ColSideColors = sc,main = "Presensitization Gut Mircobiome",
                   colCol = "white",key.title = "Log-ratio",key.xlab = paste("log( v_i / v_",i,")",sep = "" ),
                   #RowSideColors = cc,
                   cexRow = 0.5,
                   
)
dev.off()









