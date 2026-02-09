library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ape)
library(phytools)
library(randomcoloR)
library(stringr)
library(tidyr)
library(dplyr)
library(scales)
library(treeio)
library(readr)
library(ggnewscale)
library(this.path)

###################
mainpath<-paste0(dirname(this.path()),"/../")
setwd(mainpath)

FigDir<-"figures"

if (!dir.exists(FigDir)){
  dir.create(FigDir)
} else {
  print("Directory already exists!")
}
##################
#Read trees
DruH3tree <- read.iqtree("./data/phylogenetic_trees/DruH3_mmseqs98.trimmed.modelselection.IQTree.treefile")
DruE3tree <- read.iqtree("./data/phylogenetic_trees/DruE3_mmseqs98.trimmed.modelselection.IQTree.treefile")
#reroot in the middle
DruH3TreeMidRoot<-midpoint.root(DruH3tree@phylo)
DruE3TreeMidRoot<-midpoint.root(DruE3tree@phylo)

#DruH3TreeMidRoot$tip.label

##############################
##geting BootstrapValues
DruE3BootstrapValues<-get.data(DruE3tree)
#subset UFboot > 80
DruE3BootstrapValuesA80<-subset(DruE3BootstrapValues,
                                DruE3BootstrapValues$UFboot >= 80)

DruH3BootstrapValues<-get.data(DruH3tree)
#subset UFboot > 80
DruH3BootstrapValuesA80<-subset(DruH3BootstrapValues,
                                DruH3BootstrapValues$UFboot >= 80)
##############################
#Read in the data on the clusters
DruE3CladesAssignments<-read.csv("./data/phylogenetic_trees/DruE3_med_clade_3_final_clades.tsv", sep="\t")
DruH3CladesAssignments<-read.csv("./data/phylogenetic_trees/DruH3_med_clade_2.4_final_clades.tsv", sep="\t")

GetNodesOfInterest<-function(Df, tree)
{
  AncestryNodesOfInterest<-Df%>% 
    filter(!is.na(Cluster)) %>%
    group_by(Cluster) %>%
    summarise(ClMRCA=getMRCA(tree, Sequence))
  return(AncestryNodesOfInterest)
}

DruE3ClusterAncestryNodeOfInterest<-GetNodesOfInterest(DruE3CladesAssignments,
                                                       DruE3TreeMidRoot)
DruH3ClusterAncestryNodeOfInterest<-GetNodesOfInterest(DruH3CladesAssignments,
                                                       DruH3TreeMidRoot)

##set up clades colors
myclustercolors<-c("#a6cee3","#1f78b4","#b2df8a",
                   "#33a02c","#fb9a99","#e31a1c",
                   "#fdbf6f","#ff7f00","#cab2d6",
                   "#6a3d9a","#ffff99","#8c510a",
                   "#00441b","#dfc27d","#4d4d4d")
names(myclustercolors)<-unique(DruE3ClusterAncestryNodeOfInterest$Cluster)

##Draw simple tree with clusters

DruE3BasicTreePlot<-ggtree(DruE3TreeMidRoot,layout = 'circular', open.angle=40, 
                           size=0.2,
                           color="#636363")%<+% DruE3BootstrapValuesA80 +
  geom_nodepoint(aes(size = UFboot),
                 color = '#4292c6',
                 alpha=.3)+
  scale_size_continuous(range = c(0.01,1))+
  geom_treescale(y=1, x=3, fontsize=3, linesize=0.7, offset=1)
DruE3BasicTreePlot
DruE3TreeWithCladesSimple<-DruE3BasicTreePlot +
  geom_hilight(data = DruE3ClusterAncestryNodeOfInterest,
               mapping = aes(node=ClMRCA,
                             fill = as.factor(Cluster)),
               alpha=0.2,
               extend=.05)+
  scale_fill_manual(values=myclustercolors,name="Cluster")
DruE3TreeWithCladesSimple










##Test clades assignment from TreeCluster and pick the best set

DruE3CladesMed<-read.csv("../20250717_DruE_DruH_treecluster/DruE3_treecluster/DruE3_treecluster_med_clade_3",
                         sep="\t", header=T)
#This version contains 16 clusters, but let's look at them more accurately
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
#Actually this med clustering looks reasonable, aside for one clade with not very good support
#And clade 6 here include two additinal ones that have to be put in the separate clusters because they are separate on DruH tree
TempTreeWithInternalNodes<-ggtree(DruE3tree@phylo, 
       size=0.2,
       color="#636363") %<+% DruE3tree@data +
  geom_text2(aes(subset = !isTip, label = node), hjust = -0.3, size=2)
TempTreeWithInternalNodes
ggsave("DruE3_temp_internal_nodes.png", plot = TempTreeWithInternalNodes,
       path =FigDir, height=35, width =10, units="cm",dpi=300)
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
            file = "../20250717_DruE_DruH_treecluster/DruE3_med_clade_3_final_clades.tsv",
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
#do the grouping to color the branches
DruE3TreeWithCladesSimple<-DruE3BasicTreePlot +
  geom_hilight(data = DruE3ClusterAncestryNodeOfInterest,
               mapping = aes(node=ClMRCA,
                             fill = as.factor(ClNumAdj)),
               alpha=0.2,
               extend=.05)+
  scale_fill_manual(values=myclustercolors,name="Cluster")
DruE3TreeWithCladesSimple
ggsave(filename="DruE3_15_clades_simple_tree.png",
       plot=DruE3TreeWithCladesSimple,
       path=FigDir,
       dpi=300,
       width = 35,
       height =35,
       units = "cm")

#####Now do the the same for the DruH3 clusters
DruH3CladesMed<-read.csv("../20250717_DruE_DruH_treecluster/DruH3_treecluster/DruH3_treecluster_med_clade_2.4",
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
##cluster 14 is node 600
##cluster 8 is 1065
##cluster 13 is 583
##I am alos adding 15 so it will be easy to get rid of singletons
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
DruH3CladesMedF$ClNumAdj<-ifelse(is.na(DruH3CladesMedF$ClAdj),
                                 ifelse(DruH3CladesMedF$ClusterNumber %in% c(27,28),27,#5
                                        #ifelse(DruH3CladesMedF$ClusterNumber %in% c(1,2,3), 1,#13
                                        #ifelse(DruH3CladesMedF$ClusterNumber %in% c(14,20,19,15,25,24,21), 14,#14
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
            file = "../20250717_DruE_DruH_treecluster/DruH3_med_clade_2.4_final_clades.tsv",
            sep="\t",
            quote = F,
            row.names = F)
##Get MRCAs for DruH3
DruH3ClusterAncestryNodeOfInterest<-DruH3CladesMedFi %>% group_by(DruE3Clust) %>%
  summarise(ClMRCA=getMRCA(DruH3TreeMidRoot, SequenceName))
DruH3ClusterAncestryNodeOfInterest<-subset(DruH3ClusterAncestryNodeOfInterest,!is.na(DruH3ClusterAncestryNodeOfInterest$DruE3Clust))
DruH3BasicTreePlot<-ggtree(DruH3TreeMidRoot,layout = 'circular', open.angle=40, 
                           size=0.2,
                           color="#636363")%<+% DruH3BootstrapValuesA80 +
  geom_nodepoint(aes(size = UFboot),
                 color = '#4292c6',
                 alpha=.3)+
  scale_size_continuous(range = c(0.01,1))+
  geom_treescale(y=1, x=3, fontsize=3, linesize=0.7, offset=1)

DruH3TreeWithCladesSimple<-DruH3BasicTreePlot +
  geom_hilight(data = DruH3ClusterAncestryNodeOfInterest,
               mapping = aes(node=ClMRCA,fill = as.factor(DruE3Clust)),
               alpha=0.2,
               extend=.05)+
  scale_fill_manual(values=myclustercolors,name="Cluster")
DruH3TreeWithCladesSimple

##save plot
ggsave(filename="DruH3_15_clades_simple_tree.png",
       plot=DruH3TreeWithCladesSimple,
       path=FigDir,
       dpi=300,
       width = 35,
       height =35,
       units = "cm")
##############################
##############################

##############################
#########Get metadata

##Read padloc data for genome counts
DruH3padlocDf<-read.csv("./DruH3_padloc20_refseq.nopseudo.cutoff07.csv",
                   header=F)
DruE3padlocDf<-read.csv("./DruE3_padloc20_refseq.nopseudo.cutoff07.withDruH.csv",
                        header=F)
#adjust genome IDs
DruH3padlocDf$GenomeID<-sapply(strsplit(DruH3padlocDf$V1, ".csv:", fixed=TRUE), 
                                       head, 1)
DruE3padlocDf$GenomeID<-sapply(strsplit(DruE3padlocDf$V1, ".csv:", fixed=TRUE), 
                               head, 1)
##get info on clusters
DruH3Clusters<-read.csv("./DruH3_mmseqs98_cluster.withref.tsv",
                        header=F, sep="\t")
DruE3Clusters<-read.csv("./DruE3_mmseqs98_cluster.withref.tsv",
                        header=F, sep="\t")


######Read assembly summary
AssemblySummary<-read.csv("Dru3_assembly_summary.tsv", 
                          header = F, sep="\t",quote = "", 
                          row.names = NULL, 
                          stringsAsFactors = FALSE)

##get representative genomes
Dru3ReprGenomes<-read.csv("../20250424_DruantiIII_neighborhoods_genomad/Dru3_representative_genomes.txt", header=F)


####Adding cluster info to padloc df
DruH3padlocDfWithCl<-merge(DruH3padlocDf,DruH3Clusters, by.x ="V4", by.y="V2", all.x =T)
DruE3padlocDfWithCl<-merge(DruE3padlocDf,DruE3Clusters, by.x ="V4", by.y="V2", all.x =T)


#####

AssemblySummary$genus<-word(AssemblySummary$V8,1)

CommonGenus<-AssemblySummary %>% group_by(genus) %>%
  count()
CommonGenus<-CommonGenus[order(-CommonGenus$n),]
sum(pull(CommonGenus[c(1:7),2]))/length(AssemblySummary$V1)
TopGenus<-CommonGenus$genus[1:7]
AssemblySummary$TopOnly<-ifelse(AssemblySummary$genus %in% TopGenus, AssemblySummary$genus, "Else")

#picking colors
TreeColors<-c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17","#969696")
names(TreeColors)<-c(TopGenus,"Else")

########################
###Draw circular trees
########################
DrawCircularTree<-function(name,DF, treegraph)#, bootDF)
{
  # DF<-DruH3padlocDfWithCl
  # tree<-DruH3TreeMidRoot
  # name<-"DruH"
  # bootDF<-DruH3BootstrapValuesA80
  DruH3MetadataPADLOCLong<-merge(DF[,c(20,21,1,4,7,8)],
                                 AssemblySummary, by.y="V1", by.x="GenomeID")
  
  DruH3MetadataPADLOCbyWp<-DruH3MetadataPADLOCLong %>%
    group_by(V1.y,TopOnly)%>%
    count()
  names(DruH3MetadataPADLOCbyWp)<-c("ID","Genus","Count")
  
  #adding a small value so it will be visually seen which genus is there
  DruH3MetadataPADLOCbyWp$log10Count<-0.1+log10(DruH3MetadataPADLOCbyWp$Count)
  DruH3MetadataPADLOCbyWp$Genus<-factor(DruH3MetadataPADLOCbyWp$Genus, 
                                        levels=c(TopGenus, 'Else'))
  # 
  # 
  # #Plot
  # #Draw tree
  # DruH3BasicTreePlot<-ggtree(tree,layout = 'circular', open.angle=40, 
  #                            size=0.2,
  #                            color="#636363")%<+% bootDF +
  #   geom_nodepoint(aes(size = UFboot),
  #                  color = '#4292c6',
  #                  alpha=.3)+
  #   scale_size_continuous(range = c(0.01,1))+
  #   geom_treescale(y=1, x=3, fontsize=3, linesize=0.7, offset=1)
  # 
  # #+
  #   #geom_text(aes(label=B), hjust=-.5)
  # DruH3BasicTreePlot
  
  DruH3WithBarPlot<-treegraph +#DruH3BasicTreePlot +
    new_scale_fill() +
    geom_fruit(data = DruH3MetadataPADLOCbyWp, 
               geom = geom_bar, 
               mapping = aes(y=ID,
                             x=log10Count,
                             fill=Genus),
               pwidth=0.5, 
               orientation="y", 
               stat="identity",
               offset = 0.35,
               axis.params=list(axis="x",
                                text.size=3),
               grid.params=list())+
    scale_fill_manual(values=TreeColors)+
    ggtitle(name)
  DruH3WithBarPlot
  
  #save
  ggsave(filename=paste0(name,"_mmseqs98.trimmed.IQTree.tree.png"),
         plot=DruH3WithBarPlot,
         path=FigDir,
         dpi=300,
         width = 35,
         height =35,
         units = "cm")
  ggsave(filename=paste0(name,"_mmseqs98.trimmed.IQTree.tree.svg"),
         plot=DruH3WithBarPlot,
         path=FigDir,
         dpi=300,
         width = 35,
         height =35,
         units = "cm")
  return(DruH3WithBarPlot)
}

#DruH3
DrawCircularTree("DruH",DruH3padlocDfWithCl, DruH3TreeWithCladesSimple)#DruH3TreeMidRoot, DruH3BootstrapValuesA80)
DruEBaseTree<-DrawCircularTree("DruE",DruE3padlocDfWithCl, DruE3TreeWithCladesSimple)#DruE3TreeMidRoot, DruE3BootstrapValuesA80)

#DruENoMidRoot<-DrawCircularTree("DruE_no_rerooting",DruE3padlocDfWithCl, DruE3tree@phylo, DruE3BootstrapValuesA80)

##########################################################################
##############Draw connected trees

#I only want to draw connections in 555 representative genomes
DruH3padlocDfWithClRepr<-subset(DruH3padlocDfWithCl,
                                DruH3padlocDfWithCl$GenomeID %in% Dru3ReprGenomes$V1)
DruE3padlocDfWithClRepr<-subset(DruE3padlocDfWithCl,
                                DruE3padlocDfWithCl$GenomeID %in% Dru3ReprGenomes$V1)

DruE3padlocDfWithClRepr$SystemID<-sapply(strsplit(DruE3padlocDfWithClRepr$V1.x, ".csv:", fixed=TRUE), 
                               tail, 1)
DruH3padlocDfWithClRepr$SystemID<-sapply(strsplit(DruH3padlocDfWithClRepr$V1.x, ".csv:", fixed=TRUE), 
                                         tail, 1)

GetConnections<-merge(DruE3padlocDfWithClRepr[,c(20,21,22)],
      DruH3padlocDfWithClRepr[,c(20,21,22)],
      by=c("GenomeID","SystemID"))
colnames(GetConnections)[3:4]<-c("druE","druH")


#renaming DruH leaves according to connections above
druHtipsDF<-data.frame(druH=DruH3TreeMidRoot$tip.label)
#I only able to keep one occurence so some of the connections will be lost
druHtipsmerged<-merge(druHtipsDF, aggregate(druE ~ druH, data=GetConnections, head, 1), by="druH", all.x = T) 
druHtipsmerged$newlables<-ifelse(is.na(druHtipsmerged$druE),druHtipsmerged$druH,druHtipsmerged$druE)
#It is essential to order tips in the right order
druHtipsmerged<-druHtipsmerged[match(DruH3TreeMidRoot$tip.label, druHtipsmerged$druH),]

DruH3Mod<-DruH3TreeMidRoot
DruH3Mod$tip.label<-druHtipsmerged$newlables
########################################################
###Let's plot trees now


p1 <- ggtree(DruE3TreeMidRoot,
             size=0.5,
             color="#636363")

p2<-ggtree(DruH3Mod,
           size=0.5,
           color="#636363")

d1 <- p1$data
d2 <- p2$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1 + geom_tree(data=d2)

dd <- bind_rows(d1, d2) %>% 
  filter(!is.na(label))
ddtips<-subset(dd, dd$isTip =='TRUE')

DruEvsDruH<-pp + geom_line(aes(x, y, group=label), data=ddtips, color='lightblue', linewidth=.2, alpha =.8)
DruEvsDruH

ggsave(filename="DruE3_vs_DruH3_mmseqs98.trimmed.IQTree.tree.png",
       plot=DruEvsDruH,
       path=FigDir,
       dpi=300,
       width = 60,
       height =25,
       units = "cm")  

########################################
###Plotting Genomad data on the DruE phylogenetic tree

##Plasmids
pathtogenomadfolder<-"../20250424_DruantiIII_neighborhoods_genomad/genomad_output/"
filelistgenomadplasmid = list.files(pattern="\\plasmid_summary.tsv$",
                                    recursive = T,
                                    path = pathtogenomadfolder)
setwd(paste0(mainpath,"/",pathtogenomadfolder))
GenomadResultsPlasmids<-readr::read_tsv(filelistgenomadplasmid, id="file_name")
setwd(mainpath)
GenomadResultsPlasmids<-separate(data =  GenomadResultsPlasmids,
                                 col = file_name,
                                 into=c("GenomeID",NA,NA),
                                 sep="/")
colnames(GenomadResultsPlasmids)[2]<-"ID"

DruWithGenomadPlasmid<-merge(DruE3padlocDfWithClRepr,
                             GenomadResultsPlasmids,
                             by.x= c("GenomeID","V2"),
                             by.y= c("GenomeID","ID"),
                             all.x =T)
nrow(subset(DruWithGenomadPlasmid, !is.na(DruWithGenomadPlasmid$length)))
DruGenomadPlSh<-DruWithGenomadPlasmid[,c(1:3,5,8,23,27)]
colnames(DruGenomadPlSh)[7]<-"GenomadPlasmidScore"

#########################
#viruses
filelistgenomadvirus = list.files(pattern="\\virus_summary.tsv$",
                                  recursive = T,
                                  path = pathtogenomadfolder)
setwd(paste0(mainpath,"/",pathtogenomadfolder))
GenomadResultsViruses<-readr::read_tsv(filelistgenomadvirus, id="file_name")
setwd(mainpath)

GenomadResultsViruses<-separate(data = GenomadResultsViruses,
                                col = file_name,
                                into=c("GenomeID",NA,NA),
                                sep="/")
GenomadResultsViruses<-separate(data = GenomadResultsViruses,
                                col = seq_name,
                                into=c("ID","Provirus"),
                                sep="\\|", remove = F)
GenomadResultsViruses<-separate(data = GenomadResultsViruses,
                                col = coordinates,
                                into=c("Start","End"),
                                sep="\\-")

#it is essential to convert to numbers here, because otherwise it is interpreted as string
GenomadResultsViruses$Start<-as.integer(GenomadResultsViruses$Start)
GenomadResultsViruses$End<-as.integer(GenomadResultsViruses$End)

GenomadResultsVirusesCons<-subset(GenomadResultsViruses,
                                  GenomadResultsViruses$virus_score > .8)


PreGenomadVirusesOnTmnContigs<-merge(DruE3padlocDfWithClRepr,
                                     GenomadResultsVirusesCons, by.x =c("GenomeID","V2"),
                                     by=c("GenomeID","ID"))
PreGenomadVirusesOnTmnContigs$InProphage<-ifelse((PreGenomadVirusesOnTmnContigs$V12 >= PreGenomadVirusesOnTmnContigs$Start) & 
                                                   (PreGenomadVirusesOnTmnContigs$V13 <= PreGenomadVirusesOnTmnContigs$End) &
                                                   (PreGenomadVirusesOnTmnContigs$V13 <= PreGenomadVirusesOnTmnContigs$End) & 
                                                   (PreGenomadVirusesOnTmnContigs$V13 >= PreGenomadVirusesOnTmnContigs$Start),
                                                 1,0)

GenomadVirusesOnDruContigs<-subset(PreGenomadVirusesOnTmnContigs,
                                   InProphage == 1)


##merging plasmid and Viral data
DruVirPlas<-merge(DruGenomadPlSh,
                  GenomadVirusesOnDruContigs[,c(1:3,35,31,26:28)],
                  all.x =T,
                  by = c("GenomeID","V2","V4"))
#####
DruVirPlas$InProphage<-ifelse(DruVirPlas$virus_score>0, 1, 0)
DruVirPlas$InPlasmid<-ifelse(DruVirPlas$length>0, 4, 0)
MGEDataForPlot<-DruVirPlas[,c("V4","InPlasmid","InProphage")]
colnames(MGEDataForPlot)[2:3]<-c("Plasmid","Prophage")
MGEDataForPlot[is.na(MGEDataForPlot)]<-0
MGEDataForPlotLong<-gather(MGEDataForPlot,
                           key = "Location", value = "Prediction", Plasmid:Prophage)



DruEWithPlasmid<-DruEBaseTree +
new_scale_fill()+
  geom_fruit(data=MGEDataForPlotLong,
             geom = geom_tile,
             mapping = aes(y=V4, x= Location,
                           fill = as.character(Prediction)),
             pwidth=.1,
             offset = 0.05,
             axis.params=list(axis="x",
                              text.size=2,
                              text.angle=60,
                              vjust=0,
                              hjust=1))+
  scale_fill_manual(values=c("#ffffff","#bebada","#8dd3c7"), guide="none")
DruEWithPlasmid

ggsave(filename="DruE3_mmseqs98.tree.withMGE.svg",
       plot=DruEWithPlasmid,
       path=FigDir,
       dpi=300,
       width = 35,
       height =30,
       units = "cm")

ggsave(filename="DruE3_mmseqs98.tree.withMGE.png",
       plot=DruEWithPlasmid,
       path=FigDir,
       dpi=300,
       width = 35,
       height =30,
       units = "cm")

##################
##Compare my clades with DruH and DruE profiles
DruE3HMMProf<-unique(DruE3padlocDf[,c("V4","V3","V6","V5")])
DruH3HMMProf<-unique(DruH3padlocDf[,c("V4","V3","V6","V5")])

DruE3TreeWithCladesAndPADLOCHMM<-DruE3TreeWithCladesSimple %<+% DruE3HMMProf +
  geom_tippoint(aes(color = V6)) +
  scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb"),
                     name = "PADLOC HMM")
DruE3TreeWithCladesAndPADLOCHMM
DruH3TreeWithCladesAndPADLOCHMM<-DruH3TreeWithCladesSimple %<+% DruH3HMMProf +
  geom_tippoint(aes(color = V6)) +
  scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
                              "#fb9a99","#e31a1c","#fdbf6f","#ff7f00",
                              "#cab2d6","#6a3d9a","#ffff99","#b15928",
                              "#bfbfbf"),
                     name = "PADLOC HMM")
DruH3TreeWithCladesAndPADLOCHMM

ggsave(filename="DruE3_mmseqs98.tree.withPADLOCHMM.png",
       plot=DruE3TreeWithCladesAndPADLOCHMM,
       path=FigDir,
       dpi=300,
       width = 35,
       height =30,
       units = "cm")

ggsave(filename="DruH3_mmseqs98.tree.withPADLOCHMM.png",
       plot=DruH3TreeWithCladesAndPADLOCHMM,
       path=FigDir,
       dpi=300,
       width = 35,
       height =30,
       units = "cm")





