#Read Data
fecal_IgA = read_csv("Data/Total Fecal IgA_pre sens.csv")
metadata = read_csv("Output/preCompositionByStrainMetadata.csv")
features = read_csv("Output/preCompositionByStrain.csv") #relative adundances
features = data.frame(features)
labels = features[,1]
features = features[,-1]
counts = read_csv("Output/preSenCountsByStrain.csv")[,-1]





#Alpha Diveristy
shannon <- diversity(features, "shannon")
simp = diversity(features, "simpson")
invsimp <- diversity(features, "inv")

#Chao Estimates
chao = estimateR(counts)
diversity.df = data.frame(Group = labels,H = shannon,Simpson = simp,invSimpson = invsimp,t(chao))


######################################
#Shannon Entropy
#######################################
diversityIndex = "H"
div.gathered = diversity.df %>% 
  gather(key = "Diversity",value = "Index",2:ncol(diversity.df)) %>% 
  filter(Diversity == diversityIndex)
div.gathered[,1] = factor(div.gathered[,1],levels =  c("CC027","C3H"))
my_comparisons = list(c("C3H","CC027"))
#kruskall walis test
kruskal.test(Index ~ Group, data = div.gathered)
tiff(filename = "Figures/Main/preSenMicrobiome_Entropy.tiff",width = 3,height = 2.8,units = "in",res = 300)
p = ggpubr::ggboxplot(div.gathered, x = "Group", y = "Index",
                      outlier.shape = NA,
                      color = "Group",
                      palette =c("Blue", "Red"),xlab = "Strain",ylab = "Shannon Entropy",
                      add = "jitter",add.params = list(size = 3,shape = 21,linetype = "black",alpha = .5,width = .1),
                      ggtheme = theme(legend.position = "none",
                                  axis.line = element_line(),
                                  axis.text = element_text(face = "bold",size = 10),
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_text(size = 10,face = "bold"),
                                  panel.background = element_blank()
                  ))+
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label = "p.signif",label.y = 3.3,
                             tip.length = 0)# Add pairwise comparisons p-value
p = ggpubr::ggpar(p,ylim = c(1.75,3.4))
p
dev.off()



######################################
# Richness Species
######################################
diversityIndex = "S.obs"
div.gathered = diversity.df %>% 
  gather(key = "Diversity",value = "Index",2:ncol(diversity.df)) %>% 
  filter(Diversity == diversityIndex)
div.gathered$Group = factor(div.gathered$Group,levels =  c("CC027","C3H"))
#kruskall walis test
kruskal.test(Index ~ Group, data = div.gathered)
my_comparisons = list(c("C3H","CC027"))
#plots
ggplot(div.gathered,aes(Group,Index,col = Group))+
  geom_boxplot(outlier.shape = NA)+
  scale_color_manual(values = c("Red","Blue"))+
  geom_jitter(aes(fill = Group),width = .1,size = 3,alpha = 1,shape = 21,col = "black")+
  scale_fill_manual(values = c("Red","Blue"))+
  theme_classic()+
  theme(legend.position = "none")+
  ylab("Taxa Observed")
tiff(filename = "Figures/Main/preSenMicrobiome_Richness.tiff",width = 3,height = 2.8,units = "in",res = 300)
p = ggpubr::ggboxplot(div.gathered, x = "Group", y = "Index",
                      outlier.shape = NA,
                      color = "Group",
                      palette =c("Blue", "Red"),xlab = "Strain",ylab = "# of Unique Taxa Assignments",
                      add = "jitter",add.params = list(size = 3,shape = 21,linetype = "black",alpha = .5,width = .1),
                      ggtheme = theme(legend.position = "none",
                                      axis.line = element_line(),
                                      axis.text = element_text(face = "bold",size = 10),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_text(size = 10,face = "bold"),
                                      panel.background = element_blank()
                      ))+
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label = "p.signif",
                             tip.length = 0)# Add pairwise comparisons p-value
p
dev.off()
#############################################



#############################################
## Fecal IgA
#############################################
fecal_IgA = data.frame(fecal_IgA)
tiff(filename = "Figures/Main/preSenMicrobiome_FecalIgA.tiff",width = 3,height = 2.8,units = "in",res = 300)
p = ggpubr::ggboxplot(fecal_IgA, x = "Strain", y = "Fecal_IgA",
                      outlier.shape = NA,
                      color = "Strain",
                      palette =c("Blue", "Red"),xlab = "Strain",
                      ylab = "Total Fecal IgA (\u03BCg/mL)",
                      add = "jitter",add.params = list(size = 3,shape = 21,linetype = "black",alpha = .5,width = .1),
                      ggtheme = theme(legend.position = "none",
                                      axis.line = element_line(),
                                      axis.text = element_text(face = "bold",size = 10),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_text(size = 10,face = "bold"),
                                      panel.background = element_blank()
                      ))+
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label = "p.signif",label.y = 550,
                             tip.length = 0)# Add pairwise comparisons p-value
p = ggpubr::ggpar(p,ylim = c(0,575))
p
dev.off()
#########################################################################################



#################################
u = str_split(string = unqTaxa,n = 2,pattern = "\\.",simplify = T)
u = u[,1]
u = c(u,reference$ID)

#Relative Abunance
ra = data.frame(index = df.metadata$index,dat[,-1])

#Order
Order = data.frame(Taxa = cnames$ID,Order = cnames$Order)
Order = Order %>% 
  filter(Taxa %in% u) %>% 
  mutate(Order = if_else(Order%in%c("__","o__"),"Other",Order)) %>% 
  mutate(Order = str_replace(Order,"o__",""))
#Summarise RA across non unique Genera
ra =melt(data = ra, 
         id.vars = "index",
         variable.name = "Taxa", 
         measure.vars = colnames(ra)[-1])

Order = left_join(Order,ra)
Order = Order %>% 
  group_by(index,Order) %>% 
  summarise(RA = sum(value)) %>% 
  spread(key = "Order",value = "RA")
#gather
cc027 = Order[dat[,1]=="CC027",]
C3H = Order[dat[,1]=="C3H",]

#Compute Mean
prob = 0.65
cc027.mean = data.frame(Group = "CC027",
                        Order = names((mean.acomp(acomp(cc027[,-1])))),
                        RA = as.vector(mean.acomp(acomp(cc027[,-1]))))
th = quantile(cc027.mean$RA,probs = prob)
cc027.mean = cc027.mean %>% 
  mutate(Label = if_else(RA>=th,Order,""))

C3H.mean = data.frame(Group = "C3H",
                      Order = names((mean.acomp(acomp(C3H[,-1])))),
                      RA = as.vector(mean.acomp(acomp(C3H[,-1]))))
th = quantile(C3H.mean$RA,probs = prob)
C3H.mean = C3H.mean %>% 
  mutate(Label = if_else(RA>=th,Order,""))

#combine data
means = rbind(cc027.mean,C3H.mean)
# means = means %>% 
#   filter(Order!="Other")
#order by relative abundance
ra.order = means %>% 
  group_by(Order) %>% 
  summarise(meanra = mean(RA)) %>% 
  arrange(desc(-meanra))
means$Order = factor(means$Order,levels = ra.order$Order)
means$Group = factor(means$Group,levels = c("CC027","C3H"))


col <- colorRampPalette(RColorBrewer::brewer.pal(10, "Paired"))(n_distinct(means$Order))

tiff(filename = "Figures/Main/preSenMicrobiome_RA.tiff",width = 3,height = 2.8,units = "in",res = 300)
ggplot(means,aes(Group,RA,fill = Order,label = Label))+
  geom_col(width = .9,
           col = "black"
           )+
  #geom_bar(stat="identity") +
  #geom_text(position = position_stack(vjust = .5),check_overlap = F)+
  scale_fill_manual(values = col[length(col):1])+
  theme_minimal()+
  ylab("Mean Relative Abundance")+
  xlab("Strain")+
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(colour = "black",size = 1),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.height = unit(.1,units = "in"),
        axis.text.x =  element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 8),
        plot.caption = element_text(size = 8,face = "italic"))

dev.off()





#################################-
#################################-
#################################
u = str_split(string = unqTaxa,n = 2,pattern = "\\.",simplify = T)
u = u[,1]
u = c(u,reference$ID)

#Relative Abunance
ra = data.frame(index = df.metadata$index,dat[,-1])
#Genus
Genus = data.frame(Taxa = cnames$ID,Genus = cnames$Genus)
Genus = Genus %>% 
  filter(Taxa %in% u) %>% 
  mutate(Genus = if_else(Genus%in%c("__","o__"),"Other",Genus)) %>% 
  mutate(Genus = str_replace(Genus,"o__",""))
#Summarise RA across non unique Genera
ra =melt(data = ra, 
         id.vars = "index",
         variable.name = "Taxa", 
         measure.vars = colnames(ra)[-1])

Genus = left_join(Genus,ra)
Genus = Genus %>% 
  group_by(index,Genus) %>% 
  summarise(RA = sum(value)) %>% 
  spread(key = "Genus",value = "RA")
#gather
cc027 = Genus[dat[,1]=="CC027",]
C3H = Genus[dat[,1]=="C3H",]

#Compute Mean
prob = 0.65
cc027.mean = data.frame(Group = "CC027",
                        Genus = names((mean.acomp(acomp(cc027[,-1])))),
                        RA = as.vector(mean.acomp(acomp(cc027[,-1]))))
th = quantile(cc027.mean$RA,probs = prob)
cc027.mean = cc027.mean %>% 
  mutate(Label = if_else(RA>=th,Genus,""))

C3H.mean = data.frame(Group = "C3H",
                      Genus = names((mean.acomp(acomp(C3H[,-1])))),
                      RA = as.vector(mean.acomp(acomp(C3H[,-1]))))
th = quantile(C3H.mean$RA,probs = prob)
C3H.mean = C3H.mean %>% 
  mutate(Label = if_else(RA>=th,Genus,""))

#combine data
means = rbind(cc027.mean,C3H.mean)
# means = means %>% 
#   filter(Genus!="Other")
#Genus by relative abundance
ra.Genus = means %>% 
  group_by(Genus) %>% 
  summarise(meanra = mean(RA)) %>% 
  arrange(desc(-meanra))
means$Genus = factor(means$Genus,levels = ra.Genus$Genus)
means$Group = factor(means$Group,levels = c("CC027","C3H"))


col <- colorRampPalette(RColorBrewer::brewer.pal(10, "Paired"))(n_distinct(means$Genus))

tiff(filename = "Figures/Main/preSenMicrobiome_RA.tiff",width = 3,height = 2.8,units = "in",res = 300)
ggplot(means,aes(Group,RA,fill = Genus,label = Label))+
  geom_col(width = .9,
           col = "black"
  )+
  #geom_bar(stat="identity") +
  #geom_text(position = position_stack(vjust = .5),check_overlap = F)+
  #scale_fill_manual(values = col[length(col):1])+
  theme_minimal()+
  ylab("Mean Relative Abundance")+
  xlab("Strain")+
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(colour = "black",size = 1),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.height = unit(.1,units = "in"),
        axis.text.x =  element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 8),
        plot.caption = element_text(size = 8,face = "italic"))

dev.off()










#################################-
#################################-
#################################
u = str_split(string = unqTaxa,n = 2,pattern = "\\.",simplify = T)
u = u[,1]
u = c(u,reference$ID)

#Relative Abunance
ra = data.frame(index = df.metadata$index,dat[,-1])
#Order
Order = data.frame(Taxa = cnames$ID,Order = cnames$Order)
Order = Order %>% 
  filter(Taxa %in% u) %>% 
  mutate(Order = if_else(Order%in%c("__","o__"),"Other",Order)) %>% 
  mutate(Order = str_replace(Order,"o__",""))
#Summarise RA across non unique Genera
ra =melt(data = ra, 
         id.vars = "index",
         variable.name = "Taxa", 
         measure.vars = colnames(ra)[-1])

Order = left_join(Order,ra)
Order = Order %>% 
  group_by(index,Order) %>% 
  summarise(RA = sum(value)) %>% 
  spread(key = "Order",value = "RA")
#gather
cc027 = Order[dat[,1]=="CC027",]
C3H = Order[dat[,1]=="C3H",]

#Compute Mean
prob = 0.65
cc027.mean = data.frame(Group = "CC027",
                        Order = names((mean.acomp(acomp(cc027[,-1])))),
                        RA = as.vector(mean.acomp(acomp(cc027[,-1]))))
th = quantile(cc027.mean$RA,probs = prob)
cc027.mean = cc027.mean %>% 
  mutate(Label = if_else(RA>=th,Order,""))

C3H.mean = data.frame(Group = "C3H",
                      Order = names((mean.acomp(acomp(C3H[,-1])))),
                      RA = as.vector(mean.acomp(acomp(C3H[,-1]))))
th = quantile(C3H.mean$RA,probs = prob)
C3H.mean = C3H.mean %>% 
  mutate(Label = if_else(RA>=th,Order,""))

#combine data
means = rbind(cc027.mean,C3H.mean)
# means = means %>% 
#   filter(Order!="Other")
#Order by relative abundance
ra.Order = means %>% 
  group_by(Order) %>% 
  summarise(meanra = mean(RA)) %>% 
  arrange(desc(-meanra))
means$Order = factor(means$Order,levels = ra.Order$Order)
means$Group = factor(means$Group,levels = c("CC027","C3H"))


col <- colorRampPalette(RColorBrewer::brewer.pal(10, "Paired"))(n_distinct(means$Order))

tiff(filename = "Figures/Main/preSenMicrobiome_RA.tiff",width = 3,height = 2.8,units = "in",res = 300)
ggplot(means,aes(Group,RA,fill = Order,label = Label))+
  geom_col(width = .9,
           col = "black"
  )+
  #geom_bar(stat="identity") +
  #geom_text(position = position_stack(vjust = .5),check_overlap = F)+
  #scale_fill_manual(values = col[length(col):1])+
  theme_minimal()+
  ylab("Mean Relative Abundance")+
  xlab("Strain")+
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(colour = "black",size = 1),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.height = unit(.1,units = "in"),
        axis.text.x =  element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 8),
        plot.caption = element_text(size = 8,face = "italic"))

dev.off()

