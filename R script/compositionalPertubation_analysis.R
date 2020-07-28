# Author:   Andrew Hinton
# Email:    andrew84@email.unc.edu
# github:   https://github.com/andrew84830813/food_antigen_sensitization-mircobiome-genetically_susceptible_mice.git



### Read Data
dat = read_csv("Output/pertubationByStrain.csv")
df.metadata = read_csv("Output/pertubationByStrainMetadata.csv")


### start parallel cluster
ncores = detectCores()
clus <- makeCluster(ncores-1) 
registerDoParallel(clus) 


### compute clr matrix
clr_mat = data.frame(clr(dat[,-1]))
clr_mat = data.frame(easyCODA::CLR(data = dat[,-1])$LR)


### by sample pairwise distance matrix
dat = data.frame(dat)
labels = dat[,1]
d = parallelDist::parDist(as.matrix(clr_mat))


##########################################################-
### Calculate Procrustes(Greenacre) and PseduoF for all denoms in Additive Log-ratio (ALR) transform ####
##########################################################-
cn = colnames(dat[,-1])
alr.procrustes = foreach(i = 1:ncol(dat[,-1]),.combine = rbind,.packages = c("compositions","vegan"))%dopar%{
  # compute alr transformation using the ith feature as denom
  #ph.alr = data.frame(alr(dat[,-1],ivar = i))
  ph.alr = easyCODA::ALR(data = dat[,-1],denom = i)$LR 
  ph.d = parallelDist::parDist(as.matrix(ph.alr))
  # Compute procrustes correlation between alr and full logratio configurations distance matices 
  proc = vegan::protest(X = ph.d,Y = d,permutations = 100)
  # compuet pseudo f with adonis function
  a.df = data.frame(Strain = labels)
  pmv = adonis2(ph.d~Strain,data = a.df,permutations = 5000)
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
### Transformed data
keepRatios =  data.frame(Strain = dat[,1],alr(dat[,-1],ivar = i))
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
#ph.alr = data.frame(alr(dat[,-1],ivar = i))
ph.alr = easyCODA::ALR(data = dat[,-1],denom = i)$LR 
#compute pairwise sample distance using ith alr tranformation
ph.d = dist(ph.alr)
#Compute procrustes correlation between alr and full logratio configurations distance matices 
proc = vegan::protest(X = ph.d,Y = d,permutations = 100)
dvec = as.vector(as.matrix(d)[lower.tri(as.matrix(d),diag = F)])
ph.vec = as.vector(as.matrix(ph.d)[lower.tri(as.matrix(ph.d),diag = F)])

#Plot Output
tiff(filename = "Figures/Supplemental/pre-postPertubation_topDenom_PseudoF.tiff",width = 5,height = 5,units = "in",res = 300)
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
write_csv(x = alr.procrustes1,path = "Output/supplementalTable_pre-postPertubation_procrustes-pseudoF_Analysis.csv")
###############################################################3



##########################################################-
### Statistical Analysis of best alr transformed data ####
##########################################################-
#permanova
a.df = data.frame(Type = labels)
pmv = adonis2(ph.d~Type,data = a.df,permutations = 5000)
(pmv)
#dispersion test
mod = betadisper(as.dist(ph.d),group = labels)
anova(mod)
mod.HSD <- TukeyHSD(mod)
##########################################################-



##########################################################-
### Principle Coordinate Analysis with using ALR with pseudo F maximizing denom ####
##########################################################-
mds = cmdscale(ph.d,eig = T)
mds.df = data.frame(Strain_Treatment = dat[,1],mds$points)
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
shape = if_else(str_detect(labels,"CC027"),"CC027","C3H")
mds.df$Strain = shape
disp$Strain = shape
unqLabels = unique(labels)[c(1,6,2,7,3,8,4,9,5,10)]
#Relevel Factor
mds.df$Strain_Treatment = factor(mds.df$Strain_Treatment,levels = unqLabels)
disp$Strain_Treatment = factor(disp$Strain_Treatment,levels = unqLabels)
centroid$Strain_Treatment = factor(centroid$Strain_Treatment,levels = unqLabels)
#PCoA Plot
tiff(filename = "Figures/Main/pre-postPertubationMicrobiome_withALR-PCoA.tiff",width = 9,height = 5,units = "in",res = 300)
ggplot(data = mds.df,aes(X1,X2,col = Strain_Treatment,fill = Strain_Treatment))+
  geom_hline(yintercept = 0,alpha = .5)+
  geom_vline(xintercept = 0,alpha = .5)+
  xlab(paste("PCo1 (",round(pcoa_var[1],digits = 2),"% of total variation)",sep = ""))+
  ylab(paste("PCo2 (",round(pcoa_var[2],digits = 2),"% of total variation)",sep = ""))+
  geom_point(aes(shape = Strain),alpha = 1,size = 1)+
  geom_segment(aes(x = X1.y , y = X2.y , xend = X1.x , yend = X2.x , colour = Strain_Treatment),
               data = disp,alpha = .6,size = .5)+
  geom_point(data = centroid,aes(X1,X2,shape = Strain),size = 6)+
  scale_color_brewer(palette = "Paired",direction = 1)+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  labs(caption = paste('PERMANOVA F=',round(pmv$F,5),", p=",round(pmv$`Pr(>F)`,5),sep = ""))+
  theme(axis.text = element_text(hjust = 1,size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),legend.text = element_text(size = 8),legend.title = element_text(size = 8),
        axis.text.x = element_text(size = 8),plot.caption = element_text(size = 8))+
  theme(panel.grid = element_blank())
dev.off()
##########################################################-



##########################################################-
#Beta Disp
##########################################################-
interactions = factor(labels,levels = unqLabels)
mod = betadisper(as.dist(ph.d),group = interactions)
anv = anova(mod)
permutest(mod)
mod.HSD <- TukeyHSD(mod)
dis = data.frame(Treatment = interactions,Dist = mod$distances)
tiff(filename = "Figures/Main/pre-postPertubationMicrobiome_betadisp.tiff",width = 3,height = 2.25,units = "in",res = 300)
ggplot(dis,aes(Treatment,Dist,fill = Treatment))+
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot()+
  theme_classic()+
  scale_fill_brewer(palette = "Paired")+
  labs(caption = paste('PERMDISP2 F=',round(anv$`F value`[1],5),", p=",round(anv$`Pr(>F)`[1],5),sep = ""))+
  theme(axis.text = element_text(hjust = 1,size = 8,colour = "black"),legend.position = "none",
        axis.title.x = element_blank(),axis.title.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45),plot.caption = element_text(size = 8))+
  ylab("Distance to Centroid")
dev.off()
##########################################################-



##########################################################-
### Post hoc  ###
##########################################################-
alr.df = data.frame(Interaction = interactions,ph.alr)
interactions = as.character(interactions)
classes = as.character(unique(interactions))
c = combinat::combn2(classes)
pw_diffMean = foreach(i = 1:nrow(c),.combine = rbind,.packages = "vegan")%do%{
  classes = c[i,]
  ph = alr.df %>% 
    filter(Interaction %in% classes)
  phd = parallelDist::parDist(as.matrix(ph[,-1]))
  a.df = data.frame(Type = ph$Interaction)
  pmv = adonis2(phd~.,data = a.df,permutations = 100000)
  message(i)
  data.frame(T1 = classes[1],T2 = classes[2],
             pseudoF = pmv$F[1],pval = pmv$`Pr(>F)`[1])
}
pw_diffMean$pval.adj = p.adjust(pw_diffMean$pval,method = "BH")
pw_diffMean = pw_diffMean %>% 
  arrange(desc(-pval.adj))
##########################################################-



##########################################################-
### identify min signature
##########################################################-
#select key features
ph.alr = data.frame(ph.alr)
nrepeats = 10
keep.df  =data.frame()
for(s in 1:nrepeats){
  set.seed(s)
  b = Boruta::Boruta(x = ph.alr,
                     y = as.factor(interactions),doTrace = 0,maxRun = 300)
  dec = data.frame(Ratio = names(b$finalDecision),Decision = b$finalDecision)
  keep = dec %>% 
    filter(Decision!="Rejected")
  keep.df = rbind(keep.df,data.frame(seed = s,keep))
  message(s)
}
#aggregate boruta feature selection result across nrepeats
agg_keep.df = keep.df %>% 
  group_by(Ratio) %>% 
  summarise(n = n()/nrepeats)
#threshold for number of times feature selected as important
threshold = 0.8
#select important ratio given selection threshold
agg_keep.df = agg_keep.df %>% 
  filter(n>=threshold)
#select final Ratios
trainRatios = subset(ph.alr,select = as.character(agg_keep.df$Ratio))
alr.df = data.frame(Interaction = interactions,trainRatios)
#average results by interactoin
gather.alr = alr.df %>% 
  gather(key = "Ratio",value = "value",2:ncol(alr.df)) %>% 
  group_by(Interaction,Ratio) %>% 
  summarise(meanLR = mean(value)) %>% 
  spread(key = "Ratio",value = "meanLR") 
gather.alr = data.frame(gather.alr)
##########################################################-


##########################################################-
# perform cluster analysis
##########################################################-
d =  parallelDist::parallelDist(as.matrix(gather.alr[,-1]),method = "euclidean")
cls = foreach(k = 2:9,.combine = rbind)%do%{
  #Hie.Clustering
  # c = factoextra::hcut(x = d,isdiss = T,k = k,hc_method = "ward.D2")
  # data.frame(K = k,Sil = c$silinfo$avg.width)
  
  #PAM Clustering
  sil = mean(cluster::silhouette(cluster::pam(x = d,diss = T,k = k))[,3])
  data.frame(K = k,Sil = sil)
  
}
plot(cls)
optClus = cls$K[which.max(cls$Sil)]
ggplot(cls,aes(K,Sil,group = 1))+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks = cls$K)+
  theme_classic()+
  geom_vline(xintercept = optClus,lty = "dashed")+
  ylab("Avg silhouette width")+
  xlab('Number of clusters k')

mat = as.matrix(gather.alr[,-1])
rownames(mat) = gather.alr[,1]
pam.res3 <- cluster::pam(x = mat,k = optClus,stand = F,)
factoextra::fviz_cluster(object =pam.res3 , ellipse.type = "euclid", 
                         star.plot = TRUE, 
                         repel = TRUE, 
                         ggtheme = theme_bw() )
sil <- cluster::silhouette(cluster::pam(x = d,diss = T,k = optClus))
factoextra::fviz_silhouette(sil)
##########################################################-






###  Heatmap PreProcess #####
#Clustering
c = data.frame(Interaction = names(pam.res3$clustering),Cluster = pam.res3$clustering)
ncolors = max(c$Cluster)
c$Cluster = as.factor(c$Cluster)

colors = RColorBrewer::brewer.pal(n = ncolors,name = "Dark2")[1:ncolors]
color.df = data.frame(Color = colors,
                      Cluster = as.character(unique(c$Cluster)))
color.df = left_join(c,color.df)
color.df = color.df %>% 
  filter(Interaction %in% gather.alr[,1])
sc = color.df$Color
# log ratio matrix
mat = as.matrix(gather.alr[,-1])
rownames(mat) = gather.alr[,1]
col <- colorRampPalette(RColorBrewer::brewer.pal(10, "PuOr"))(1000)
mat = t(mat)
# heatmap microbiome composition
#tiff(filename = "Figures/Main/pre-postPertubationMicrobiome_withALR-interactionHeatmap.tiff",width = 9,height = 5,units = "in",res = 300)
gplots::heatmap.2( mat,
                   col = col,#viridis::viridis(n = 1000,direction = -1,option = "D"),#gplots::redblue(n = 1000) ,#viridis::viridis(n = 1000,option = "D"),#,
                   Rowv = TRUE,
                   #margins = c(2, 2),
                   hclustfun = function(x) hclust(x, method = "average"),
                   distfun = function(x) parallelDist::parallelDist(x,method = "euclidean"),
                   key=TRUE, symkey=TRUE, 
                   scale = "none",
                   density.info="density", trace="none",
                   ColSideColors = sc,
                   #main = "Presensitization Gut Mircobiome",
                   key.title = "Log-ratio",key.xlab = paste("log( v_i / v_",i,")",sep = "" ),
                   #RowSideColors = cc,
                   cexRow = 0.5,cexCol = .75,adjCol = c(.85,1)
                   
)
dev.off()


#Visualize Distance Matrix
d =  parallelDist::parallelDist(as.matrix(gather.alr[,-1]),method = "euclidean")
cc = data.frame(Interaction = names(pam.res3$clustering),Cluster = pam.res3$clustering) %>% 
  arrange(desc(Cluster))
cc$Cluster = as.factor(cc$Cluster)
cc = left_join(cc,color.df)
dmat = as.matrix(d)
dmat = 1/dmat^2
diag(dmat) = 0
diag(dmat) = max(dmat)
colnames(dmat) = gather.alr[,1]
rownames(dmat) = gather.alr[,1]
dmat = dmat[cc$Interaction,cc$Interaction]
dmat = data.frame(T1 = cc$Interaction,dmat) 
dmat =  gather(dmat,key = "T2",value = "Dist",2:ncol(dmat))
dmat$T1 =factor(dmat$T1,levels = cc$Interaction)
dmat$T2 =factor(dmat$T2,levels = cc$Interaction)
th = .35
th.val = quantile(dmat$Dist,probs = th)
ggplot(dmat,aes(T1,T2,fill = Dist))+
  geom_tile()+
  scale_fill_distiller(palette = "BuPu",direction = 1)+
  theme_bw()+
  theme(axis.text = element_text(colour = cc$Color,face = "bold"))
dmat = dmat %>% 
  mutate(Dist = if_else(Dist>th.val,0,1))
dmat$Dist = as.factor(dmat$Dist)
ggplot(dmat,aes(T1,T2,fill = Dist))+
  geom_tile(col = "white")+
  scale_fill_manual(values = c("white","grey40"))+
  theme_bw()+
  theme(axis.text = element_text(colour = cc$Color,face = "bold"))




#Compute average LR values by cluster
clus_alr = right_join(c,alr.df)

gather.alr = clus_alr %>% 
  gather(key = "Ratio",value = "value",3:ncol(clus_alr)) %>% 
  group_by(Cluster,Ratio) %>% 
  summarise(meanLR = mean(value),n = n(),sd = sd(value),CI = (1-.96*sd)/sqrt(n),lb = meanLR - CI, ub = meanLR + CI) %>% 
  separate(2,into = c("ID","denom"),sep = "\\.")
gather.alr1 = right_join(cnames,gather.alr)
gather.alr1$Family_Genus = paste(gather.alr1$Family,gather.alr1$Genus,gather.alr1$Species,sep = "|") 
gather.alr1 = gather.alr1 %>% 
  filter(Family_Genus!="f__|g__|s__")
#visualize pertubation signature
ggplot(gather.alr1,aes(Family_Genus,meanLR,fill = meanLR))+
  geom_col(col = "black")+
  geom_errorbar(aes(ymin = lb , ymax = ub ), width=0.2)+
  facet_wrap(.~Cluster,nrow = 1)+
  coord_flip()+
  scale_fill_distiller(palette = "RdBu",direction = 1)+
  scale_color_brewer(palette = "Set3")+
  ggridges::theme_ridges()+
  theme(legend.position = "none",
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 7),
        axis.title.x = element_text(hjust = 0.5,size = 8),
        axis.title.y = element_blank())+
  geom_hline(yintercept = 0,col = "black",size = 1)+
  ylab("Mean logratio")



