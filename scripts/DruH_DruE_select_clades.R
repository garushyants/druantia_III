library(ggtree)
library(ggplot2)
library(ape)
library(phytools)
library(stringr)
library(dplyr)
library(scales)
library(treeio)
library(this.path)

mainpath<-paste0(dirname(this.path()),"/../")
setwd(mainpath)

#Read trees
DruH3tree <- read.iqtree("./data/phylogenetic_trees/DruH3_mmseqs98.trimmed.modelselection.IQTree.treefile")
DruE3tree <- read.iqtree("./data/phylogenetic_trees/DruE3_mmseqs98.trimmed.modelselection.IQTree.treefile")
#reroot in the middle
DruH3TreeMidRoot<-midpoint.root(DruH3tree@phylo)
DruE3TreeMidRoot<-midpoint.root(DruE3tree@phylo)


##Test clades assignment from TreeCluster and pick the best set

DruE3CladesMed<-read.csv("./data/treecluster/DruE3_treecluster_med_clade_3",
                         sep="\t", header=T)
############
#Drawing basic DruE3 tree
DruE3BasicTreePlot<-ggtree(DruE3TreeMidRoot,layout = 'circular', open.angle=40, 
                           size=0.2,
                           color="#636363")%<+% DruE3BootstrapValuesA80 +
  geom_nodepoint(aes(size = UFboot),
                 color = '#4292c6',
                 alpha=.3)+
  scale_size_continuous(range = c(0.01,1))+
  geom_treescale(y=1, x=3, fontsize=3, linesize=0.7, offset=1)
DruE3BasicTreePlot

####Adding tippoint with clusters to understand what have to be merged
ClNum<-length(unique(DruE3CladesMed$ClusterNumber))
colors <- distinctColorPalette(ClNum)
DruE3BasicTreePlot %<+% DruE3CladesMed +
  geom_tippoint(aes(color=as.factor(ClusterNumber)), size=2) +
  scale_color_manual(values=colors)
# #And clade 6 here include two additinal ones that have to be put in the separate clusters because they are separate on DruH tree
# TempTreeWithInternalNodes<-ggtree(DruE3tree@phylo, 
#                                   size=0.2,
#                                   color="#636363") %<+% DruE3tree@data +
#   geom_text2(aes(subset = !isTip, label = node), hjust = -0.3, size=2)
# TempTreeWithInternalNodes
# ggsave("DruE3_temp_internal_nodes.png", plot = TempTreeWithInternalNodes,
#        path =FigDir, height=35, width =10, units="cm",dpi=300)
###From that I know
##Cluster 6 in this initial assignment has to be split in 3 with mrca:968,1057,1081
cluster2nodes<-c(968,1057,1081)
names(cluster2nodes)<-c(2,3,4)
Cluster2RestructuredDf<- lapply(cluster2nodes, function(n) {
  subtree <- extract.clade(DruE3TreeMidRoot, node = n)
  tips <- subtree$tip.label   # get tip numbers
  data.frame(
    ClAdj = rep(names(cluster2nodes)[cluster2nodes == n],length(tips)),
    SequenceName =tips,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
#I want to merge -1,1-5 into one cluster and then it seems like a reasonable assignment
DruE3CladesMedF<-merge(DruE3CladesMed,Cluster2RestructuredDf,by="SequenceName", all=T)

DruE3CladesMedF$ClNumAdj<-ifelse(DruE3CladesMedF$ClusterNumber < 6, 1,
                                 ifelse(DruE3CladesMedF$ClusterNumber == 6,
                                        as.integer(DruE3CladesMedF$ClAdj),
                                        DruE3CladesMedF$ClusterNumber-2))
##I get 15 clusters out of it that all seems reasonable
#Saving final clusters to a file
EClustToSave<-DruE3CladesMedF[,c(1,4)]
names(EClustToSave)<-c("Sequence","Cluster")
EClustToSave<-EClustToSave[order(EClustToSave$Cluster),]
write.table(EClustToSave,
            file = "./data/phylogenetic_trees/DruE3_med_clade_3_final_clades.tsv",
            sep="\t",
            quote = F,
            row.names = F)
##

myclustercolors<-c("#a6cee3","#1f78b4","#b2df8a",
                   "#33a02c","#fb9a99","#e31a1c",
                   "#fdbf6f","#ff7f00","#cab2d6",
                   "#6a3d9a","#ffff99","#8c510a",
                   "#00441b","#4d4d4d","#dfc27d")
names(myclustercolors)<-unique(DruE3CladesMedF$ClNumAdj)

##get common ancestors for all clusters so I can do the blocks on the tree

#getting the common ancestor
DruE3ClusterAncestryNodeOfInterest<-DruE3CladesMedF %>% group_by(ClNumAdj) %>%
  summarise(ClMRCA=getMRCA(DruE3TreeMidRoot, SequenceName))
##
#Checking that it works as expected
DruE3TreeWithCladesSimple<-DruE3BasicTreePlot +
  geom_hilight(data = DruE3ClusterAncestryNodeOfInterest,
               mapping = aes(node=ClMRCA,
                             fill = as.factor(ClNumAdj)),
               alpha=0.2,
               extend=.05)+
  scale_fill_manual(values=myclustercolors,name="Cluster")
DruE3TreeWithCladesSimple

#####Now do the the same for the DruH3 clusters
DruH3CladesMed<-read.csv("./data/treecluster//DruH3_treecluster_med_clade_2.4",
                         sep="\t", header=T)
###There are problems with this clustering because I miss the stem of of DruE 14 and the whole DruE 8,
###because some of the sequences are assigned as singletons
###Do the same procedure as above for DruE
HTempTreeWithInternalNodes<-ggtree(DruH3tree@phylo, 
                                   size=0.2,
                                   color="#636363") %<+% DruE3tree@data +
  geom_text2(aes(subset = !isTip, label = node), hjust = -0.3, size=2)
HTempTreeWithInternalNodes
###
#do some restructuring
Hcluster2nodes<-c(601,1065,583,746)
names(Hcluster2nodes)<-c(54,48,43,45)
HCluster2RestructuredDf<- lapply(Hcluster2nodes, function(n) {
  subtree <- extract.clade(DruH3TreeMidRoot, node = n)
  tips <- subtree$tip.label   # get tip numbers
  data.frame(
    ClAdj = rep(names(Hcluster2nodes)[Hcluster2nodes == n],length(tips)),
    SequenceName =tips,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
DruH3CladesMedF<-merge(DruH3CladesMed,HCluster2RestructuredDf,by="SequenceName", all=T)
###################
##Aligning it with Clusters for DruE3
DruH3CladesMedF$ClNumAdj<-ifelse(is.na(DruH3CladesMedF$ClAdj),
                                 ifelse(DruH3CladesMedF$ClusterNumber %in% c(27,28),27,#5
                                              ifelse(DruH3CladesMedF$ClusterNumber %in% c(22,11,17,16), 11,#15
                                               ifelse(DruH3CladesMedF$ClusterNumber %in% c(6,9,10), 6,#9
                                                      DruH3CladesMedF$ClusterNumber))), DruH3CladesMedF$ClAdj)
##Renumbering clusters so they would match the corresponding ones in DruE3
ClDruH3DF<-data.frame(DruH3ClustOld = c(13,18,26,23,27,12,8,48,6,4,7,5,43,54,45),
                      DruE3Clust = seq(1:15))

DruH3CladesMedFi<-merge(DruH3CladesMedF,ClDruH3DF, by.x="ClNumAdj", by.y = "DruH3ClustOld", all.x = T)
#After that I get 15 clusters as expected
##Saving final list to a file
HClustToSave<-DruH3CladesMedFi[,c(2,5)]
names(HClustToSave)<-c("Sequence","Cluster")
HClustToSave<-HClustToSave[order(HClustToSave$Cluster),]
write.table(HClustToSave,
            file = "./data/phylogenetic_trees/DruH3_med_clade_2.4_final_clades.tsv",
            sep="\t",
            quote = F,
            row.names = F)
