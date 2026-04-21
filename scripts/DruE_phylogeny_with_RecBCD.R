library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ape)
library(phytools)
library(treeio)
library(this.path)
library(dplyr)
library(ggnewscale)

mainpath<-paste0(dirname(this.path()),"/../")
setwd(mainpath)

FigDir<-"figures/DruE_RecBCD_analysis"

if (!dir.exists(FigDir)){
  dir.create(FigDir,recursive = TRUE)
} else {
  print("Directory already exists!")
}
###
DruEtree <- read.iqtree("./data/phylogenetic_trees/DruE3_mmseqs98.trimmed.modelselection.IQTree.treefile")
DruETreeMidRoot<-midpoint.root(DruEtree@phylo)
##
AdditionalDruData<-read.csv("./data/Supplementary_table_1_dataset_info.tsv", sep ="\t", header =T)
AdditionalDruDataRepr<-subset(AdditionalDruData, AdditionalDruData$representative_genome == "Y" &
                                AdditionalDruData$Protein == "DruE3")
##
GetNodesOfInterest<-function(Df, tree)
{
  AncestryNodesOfInterest<-Df%>% 
    filter(!is.na(Cluster)) %>%
    group_by(Cluster) %>%
    summarise(ClMRCA=getMRCA(tree, TreeRepresentative))
  return(AncestryNodesOfInterest)
}
DruE3ClusterAncestryNodeOfInterest<-GetNodesOfInterest(AdditionalDruDataRepr,
                                                       DruETreeMidRoot)
##set up clades colors
myclustercolors<-c("#a6cee3","#1f78b4","#b2df8a",
                   "#33a02c","#fb9a99","#e31a1c",
                   "#8c510a","#ff7f00","#cab2d6",
                   "#ffff99","#6a3d9a","#fdbf6f",
                   "#00441b","#dfc27d","#4d4d4d")
names(myclustercolors)<-unique(DruE3ClusterAncestryNodeOfInterest$Cluster)

###Basic tree
BasicTreePlot<-ggtree(DruETreeMidRoot,layout = 'fan', open.angle=10, 
                      size=0.2,
                      color="#636363")+
  geom_treescale(y=1, x=3.2, fontsize=3, linesize=0.7, offset=1)
TreeWithCladesSimple<-BasicTreePlot +
  geom_hilight(data = DruE3ClusterAncestryNodeOfInterest,
               mapping = aes(node=ClMRCA,
                             fill = as.factor(Cluster)),
               alpha=0.2,
               extend=.05)+
  scale_fill_manual(values=myclustercolors,name="Cluster")

#########Get HMMer Rec search outputs
RecBraw<-read.csv("./data/RecBCD_hmmer/RecB_representative_genomes_hmmer_2060421.tab", header = F, sep="\t")
RecCraw<-read.csv("./data/RecBCD_hmmer/RecC_representative_genomes_hmmer_2060421.tab", header = F, sep="\t")
RecDraw<-read.csv("./data/RecBCD_hmmer/RecD_representative_genomes_hmmer_2060421.tab", header = F, sep="\t")
###RecC hits are the most reliable ones with good hits
#get best record for genome
RecCbest <- RecCraw %>%
  group_by(V1) %>%
  slice_min(V8, n = 1, with_ties = FALSE) %>%
  ungroup()
RecDbest <- RecDraw %>%
  group_by(V1) %>%
  slice_min(V8, n = 1, with_ties = FALSE) %>%
  ungroup()
RecBbest <- RecBraw %>%
  group_by(V1) %>%
  slice_min(V8, n = 1, with_ties = FALSE) %>%
  ungroup()

RecForPlot<-merge(RecCbest[,c(1,4,6,7)], AdditionalDruDataRepr[,c(1,3)], all.y = T, by.x ="V1", by.y = "RefseqID")
RecDForPlot<-merge(RecDbest[,c(1,4,6,7)], AdditionalDruDataRepr[,c(1,3)], all.y = T, by.x ="V1", by.y = "RefseqID")
RecBForPlot<-merge(RecBbest[,c(1,4,6,7)], AdditionalDruDataRepr[,c(1,3)], all.y = T, by.x ="V1", by.y = "RefseqID")
##add it to the tree
DruETreeWithRecBCD<-TreeWithCladesSimple +
  new_scale_fill()+
  geom_fruit(data=RecBForPlot,
             geom = geom_tile,
             mapping = aes(y=TreeRepresentative, x= V4,
                           fill = V7),
             color="#bdbdbd",
             pwidth=.1,
             offset = 0.1)+
  scale_fill_gradient(high="#08519c", low ="#f7fbff", name ="RecB hmmsearch bitscore", na.value = "#ffffff")+
  new_scale_fill()+
  geom_fruit(data=RecForPlot,
             geom = geom_tile,
             mapping = aes(y=TreeRepresentative, x= V4,
                           fill = V4),
             color="#bdbdbd",
             pwidth=.1,
             offset = 0.022)+
  scale_fill_manual(values=c("#6a3d9a"), name ="", na.value = NA,na.translate = FALSE)+
  new_scale_fill()+
  geom_fruit(data=RecDForPlot,
             geom = geom_tile,
             mapping = aes(y=TreeRepresentative, x= V4,
                           fill = V7),
             color="#bdbdbd",
             pwidth=.1,
             offset = 0.022)+
  scale_fill_gradient(high="#7a0177", low ="#fff7f3", name ="RecD hmmsearch bitscore", na.value = "#ffffff")

DruETreeWithRecBCD
###add info on where Zorya is also found
readPADLOC<-function(file)
{
  df<-read.csv(file,header = F)
  df$RefseqID<-sub("\\.csv.*", "", df$V1)
  dfessential<-df[,c(4,20,2,3,7,6,12:14)]
  names(dfessential)<-c("ProteinID","RefseqID","Contig","System","Protein","PADLOC_HMM","Start","End","Strand")
  return(dfessential)
}
ZorEPADLOC<-readPADLOC("./data/ZorE2_padloc20_refseq.nopseudo.activenucl.csv")

ZorEinReprGenomes<-merge(ZorEPADLOC, AdditionalDruDataRepr[,c(1,3)], by ="RefseqID")
#Plot
DruETreeWithRecCDZorE<-DruETreeWithRecBCD + new_scale_fill()+
  geom_fruit(data=ZorEinReprGenomes,
             geom = geom_tile,
             mapping = aes(y=TreeRepresentative, x= Protein,
                           fill = Protein),
             color="#bdbdbd",
             pwidth=.1,
             offset = 0.025)+
  scale_fill_manual(values=c("#006837"), name ="", na.value = "#ffffff")
DruETreeWithRecCDZorE

ggsave("DruE_tree_with_RecBCD_ZorE.pdf",
       plot=DruETreeWithRecCDZorE,
       path=FigDir,
       width=40,
       height=30,
       dpi=300,
       units="cm")