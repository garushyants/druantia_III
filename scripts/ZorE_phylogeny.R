library(ggtree)
library(ggplot2)
library(ape)
library(phytools)
library(treeio)
library(this.path)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(ggtreeExtra)
library(Biostrings)

###################
mainpath<-paste0(dirname(this.path()),"/../")
setwd(mainpath)

FigDir<-"figures/ZorE_DruE_analysis"

if (!dir.exists(FigDir)){
  dir.create(FigDir,recursive = TRUE)
} else {
  print("Directory already exists!")
}
##################
#Read trees
ZorEtree <- read.iqtree("./data/phylogenetic_trees/ZorE2_mmseqs98.activenucl.trimmed.IQTree.treefile")
DruEtree <- read.iqtree("./data/phylogenetic_trees/DruE3_mmseqs98.trimmed.modelselection.IQTree.treefile")
#reroot in the middle
ZorETreeMidRoot<-midpoint.root(ZorEtree@phylo)
DruETreeMidRoot<-midpoint.root(DruEtree@phylo)

##################
###get taxonomy and genome counts
##############################
##geting BootstrapValues
DruEBootstrapValues<-get.data(DruEtree)
#subset UFboot > 80
DruEBootstrapValuesA80<-subset(DruEBootstrapValues,
                                DruEBootstrapValues$UFboot >= 80)

ZorEBootstrapValues<-get.data(ZorEtree)
#subset UFboot > 80
ZorEBootstrapValuesA80<-subset(ZorEBootstrapValues,
                               ZorEBootstrapValues$UFboot >= 80)
##############################
###MMseqs clusters
readMMseqs<-function(file)
{
  MMseqsdf<-read.csv(file, sep ="\t", header = F)
  names(MMseqsdf)<-c("TreeRepresentative","ProteinID")
  return(MMseqsdf)
}
#
DruE3MMseqs<-readMMseqs("./data/mmseqs/DruE3_mmseqs98_cluster.withref.tsv")
ZorEMMseqs<-readMMseqs("./data/mmseqs/ZorE2_mmseqs98_cluster.activenucl.tsv")
##PADLOC data
readPADLOC<-function(file)
{
  df<-read.csv(file,header = F)
  df$RefseqID<-sub("\\.csv.*", "", df$V1)
  dfessential<-df[,c(4,20,2,3,7,6,12:14)]
  names(dfessential)<-c("ProteinID","RefseqID","Contig","System","Protein","PADLOC_HMM","Start","End","Strand")
  return(dfessential)
}
#
DruE3PADLOC<-readPADLOC("./data/DruE3_padloc20_refseq.nopseudo.cutoff07.withDruH.csv")
ZorEPADLOC<-readPADLOC("./data/ZorE2_padloc20_refseq.nopseudo.activenucl.csv")

####################
###Assembly summary
readAssemblySummary<-function(file)
{
  ASdf<-read.csv(file, sep="\t", header =F)
  ASessential<-ASdf[,c(1,7,8,9,12)]
  names(ASessential)<-c("RefseqID","Taxid",
                                     "Name",
                                     "Strain","AssemblyStatus")
  return(ASessential)
}
DruEAS<-readAssemblySummary("./data/DruE3_assembly_summary_refseq_20250415.txt")
ZorEAS<-readAssemblySummary("./data/ZorE_assembly_summary_refseq_20250415.txt")
###
###Taxonomic information
readtaxonomy<-function(taxpath)
{
  TaxselectedFiles<-list.files(pattern="\\.tsv$",
                               path = taxpath)
  setwd(paste0(mainpath,"/",taxpath))
  Taxdata<-readr::read_tsv(TaxselectedFiles, id="file_name")
  setwd(mainpath)
  TaxdataEssential<-unique(Taxdata[,c(2,11,12,14,16,18,20,22,24,26)])
  names(TaxdataEssential)[1]<-"Taxid"
  return(TaxdataEssential)
}
DruEtaxinfo<-readtaxonomy("./data/taxonomy_info_datasets/")
ZorEtaxinfo<-readtaxonomy("./data/ZorE_taxonomy_info_datasets/")

#############################
MergeMetadata<-function(PADLOC,AssemblySum,Tax,MMseqs)
{
  DfwithGenomeInfo<-merge(PADLOC,
                          AssemblySum,
                          by="RefseqID", all.x = T)
  DfwithTaxonomy<-merge(DfwithGenomeInfo,
                        Tax,
                        by="Taxid",
                        all.x = T)
  DfTaxMMseqs<-merge(DfwithTaxonomy,
                     MMseqs,
                     by="ProteinID",
                     all.x = T)
  return(DfTaxMMseqs)
}
ZorEMetadataall<-MergeMetadata(ZorEPADLOC,ZorEAS,ZorEtaxinfo,ZorEMMseqs)
DruEMetadataall<-MergeMetadata(DruE3PADLOC,DruEAS,DruEtaxinfo,DruE3MMseqs)
#############################
##Plot ZorE phylogenetic tree
ZorEBasicTree<-ggtree(ZorETreeMidRoot,
                      layout = 'fan', open.angle=10, 
                      size=0.2,
                      color="#636363")%<+% ZorEBootstrapValuesA80 +
  geom_nodepoint(aes(size = UFboot),
                 color = '#4292c6',
                 alpha=.3)+
  scale_size_continuous(range = c(0.01,1))+
  geom_treescale(y=1, x=3.2, fontsize=3, linesize=0.7, offset=1)
ZorEBasicTree
############################
##Adding more info on taxonomy
ZorECommonGenusInSet<-ZorEMetadataall %>% group_by(`Genus name`) %>%
  count()
ZorECommonGenusInSet<-ZorECommonGenusInSet[order(-ZorECommonGenusInSet$n),]

sum(pull(ZorECommonGenusInSet[c(1:11),2]))/length(ZorEMetadataall$ProteinID)
#I am taking 11 most common genus (>30 genomes), because they together account for ~88% of all genomes
ZorEGenusToKeep<-c(pull(ZorECommonGenusInSet[c(1:11),1]))
ZorEMetadataall$TopGenus<-ifelse(ZorEMetadataall$`Genus name` %in% ZorEGenusToKeep,
                                 ZorEMetadataall$`Genus name`, "Other")

ZorEHitsCountsPerLeaf<-ZorEMetadataall %>%
  group_by(TreeRepresentative)%>%
  count()
#transform to log10
ZorEHitsCountsPerLeaf$logCount<-log10(ZorEHitsCountsPerLeaf$n)+.1

ZorEGenusCountsPerLeaf<-ZorEMetadataall %>%
  group_by(TreeRepresentative, TopGenus) %>% 
  count(TopGenus) %>%
  group_by(TreeRepresentative) %>%
  mutate(percent = n / sum(n) * 100)
ZorEGenusCountsPerLeaf$TopGenus<-factor(ZorEGenusCountsPerLeaf$TopGenus,
                                       levels=c(ZorEGenusToKeep,"Other"))

ZorECommonClassInSet<-ZorEMetadataall %>% group_by(`Class name`) %>%
  count()%>%
  arrange(desc(n))
ZorEClassDf<-ZorEMetadataall[,c("TreeRepresentative","Class name")]
ZorEClassDf$CommonClass<-ifelse(ZorEClassDf$`Class name` %in% ZorECommonClassInSet$`Class name`[1:7],
                            ZorEClassDf$`Class name`,
                            "Other")
#I pick all Classes with > 10 genomes
ZorEClassCountsPerLeaf<-ZorEClassDf %>%
  group_by(TreeRepresentative, CommonClass) %>% 
  count(CommonClass) %>%
  group_by(TreeRepresentative) %>%
  mutate(percent = n / sum(n) * 100)
unique(ZorEClassCountsPerLeaf$CommonClass)
ZorEClassCountsPerLeaf$CommonClass<-factor(ZorEClassCountsPerLeaf$CommonClass,
                                       levels=c(ZorECommonClassInSet$`Class name`[1:7],
                                                "Other"))
##Also adding info for ZoryaII found in the same genome with Druantia III
ZorEWithDruall<-merge(ZorEMetadataall,
                      DruEMetadataall[,c(3,4,6,8:10,23)],
                      by = "RefseqID",
                      all.x =T)
ZorEWithDruall$distance<-ifelse(!is.na(ZorEWithDruall$Start.y),
                                ifelse(ZorEWithDruall$Contig.x == ZorEWithDruall$Contig.y,
                                       abs((ZorEWithDruall$Start.y+ZorEWithDruall$End.y-ZorEWithDruall$Start.y)/2 -
                                             (ZorEWithDruall$Start.x+ZorEWithDruall$End.x-ZorEWithDruall$Start.x)/2),
                                       Inf),
                                NA)
ZorEWithDruall$Dru3<-ifelse(!is.na(ZorEWithDruall$Start.y),
                            ifelse(ZorEWithDruall$distance < 20000,
                              ifelse(ZorEWithDruall$Strand.x == ZorEWithDruall$Strand.y,
                                     ifelse(ZorEWithDruall$distance < 6435,
                                            "adjucent",
                                            "neighborhood (20k)"),
                                     ifelse(ZorEWithDruall$distance < 3551,
                                            "adjucent","neighborhood (20k)")),
                              "same genome"),
                              NA)
ZorEWithDru<- ZorEWithDruall[,c(23,6,32)] %>% 
  group_by(TreeRepresentative.x,Protein.x) %>%
  summarise(
    Dru3 = if (all(is.na(Dru3))) NA else dplyr::first(sort(Dru3[!is.na(Dru3)])),
    .groups = "drop"
  )
###############################
genuscolors<-c("#e31a1c","#ff7f00","#b2df8a","#33a02c",
               "#fb9a99","#fdbf6f","#a6cee3","#1f78b4",
               "#ffff99","#cab2d6","#6a3d9a","#d9d9d9")
names(genuscolors)<-c(ZorEGenusToKeep,"Other")

classcolors<-c("#fbb4ae","#b3cde3","#ccebc5",
               "#decbe4","#fed9a6","#ffffcc","#fddaec","#d9d9d9")
names(classcolors)<-c(ZorECommonClassInSet$`Class name`[1:7],"Other")
############################
ZorETreeToSave<-ZorEBasicTree+
  new_scale_fill()+
  geom_fruit(data=ZorEWithDru,
             geom = geom_tile,
             mapping =aes(y = TreeRepresentative.x,
                          x=Protein.x,
                          fill = Dru3),
             offset =0.01, pwidth =.1)+
  scale_fill_manual(values = c("#542788","#80cdc1","#dfc27d"),
                    na.value = "white",
                    na.translate = FALSE,
                    name = "Druantia III location")+
  new_scale_fill()+
  geom_fruit(data=ZorEClassCountsPerLeaf,
             geom = geom_col,
             mapping = aes(y=TreeRepresentative,
                           fill = CommonClass,
                           x = percent),
             offset =0.02, pwidth =.02)+
  scale_fill_manual(values=classcolors, name = "Class")+
  guides(fill = guide_legend(nrow = 4))+
  new_scale_fill()+
  geom_fruit(data=ZorEGenusCountsPerLeaf,
             geom = geom_col,
             mapping = aes(y=TreeRepresentative,
                           fill = TopGenus,
                           x = percent),
             offset =0.01, pwidth =.06)+
  scale_fill_manual(values=genuscolors, name = "Genus")+
  guides(fill = guide_legend(nrow = 4))+
  geom_fruit(data = ZorEHitsCountsPerLeaf,
             geom = geom_col,
             mapping = aes(y=TreeRepresentative,
                           x=logCount), fill= "#878787",
             pwidth =.4,
             axis.params=list(axis="x",
                              text.size=3,
                              line.size=.3),
             grid.params = list(size=.3,
                                alpha=.3))
ZorETreeToSave           
ggsave("ZorE_tree_with_Tax_and_DruE_location.pdf",
       plot=ZorETreeToSave,
       path=FigDir,
       width=40,
       height=30,
       dpi=300,
       units="cm")
############################
###Read DruE clusters and plotbasic tree
DruE3CladesAssignments<-read.csv("./data/phylogenetic_trees/DruE3_med_clade_3_final_clades.tsv", sep="\t")

GetNodesOfInterest<-function(Df, tree)
{
  AncestryNodesOfInterest<-Df%>% 
    filter(!is.na(Cluster)) %>%
    group_by(Cluster) %>%
    summarise(ClMRCA=getMRCA(tree, Sequence))
  return(AncestryNodesOfInterest)
}
DruE3ClusterAncestryNodeOfInterest<-GetNodesOfInterest(DruE3CladesAssignments,
                                                       DruETreeMidRoot)
myclustercolors<-c("#a6cee3","#1f78b4","#b2df8a",
                   "#33a02c","#fb9a99","#e31a1c",
                   "#8c510a","#ff7f00","#cab2d6",
                   "#ffff99","#6a3d9a","#fdbf6f",
                   "#00441b","#dfc27d","#4d4d4d")
names(myclustercolors)<-unique(DruE3ClusterAncestryNodeOfInterest$Cluster)

DruEBasicTreePlot<-ggtree(DruETreeMidRoot,
                      size=0.5,
                      color="#525252") + #%<+% DruEBootstrapValuesA80 +
  # geom_nodepoint(aes(size = UFboot),
  #                color = '#4292c6',
  #                alpha=.6)+
  scale_size_continuous(range = c(0.01,1))+
  geom_treescale(y=1, x=3.2, fontsize=3, linesize=0.7, offset=1)
DruETreeWithCladesSimple<-DruEBasicTreePlot +
  geom_hilight(data = DruE3ClusterAncestryNodeOfInterest,
               mapping = aes(node=ClMRCA,
                             fill = as.factor(Cluster)),
               alpha=0.2,
               extend=.05)+
  scale_fill_manual(values=myclustercolors,name="Cluster")
DruETreeWithCladesSimple

####Connection data
Connectionsall<-ZorEWithDruall[,c(1,4,23,30)]
Connections<-unique(subset(Connectionsall,
                    !is.na(Connectionsall$TreeRepresentative.y)))
names(Connections)[c(3,4)]<-c("ZorE","DruE")
Connections$group<-seq(1,length(Connections$RefseqID))
ConnectionsWithClades<-merge(Connections,DruE3CladesAssignments,
                             by.x = "DruE",
                             by.y ="Sequence", all.x =T)

Connectionslong<-ConnectionsWithClades %>%
  pivot_longer(cols = c("ZorE","DruE"),
               names_to = "gene",
               values_to= "label")

#####
ZorETree<-ggtree(ZorETreeMidRoot) +
  scale_size_continuous(range = c(0.01,1))

Edata <- DruETreeWithCladesSimple$data
Zdata <- ZorETree$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
Zdata$x <- max(Zdata$x) - Zdata$x + max(Edata$x) + 1
Zdata$y <- Zdata$y*(max(Edata$y)/max(Zdata$y))
#Combine trees together
TreesSideBySide <- DruETreeWithCladesSimple + geom_tree(data=Zdata,
                                              size=0.5,
                                              color="#525252") 
##get coordinates for all tips on both trees
dd <- bind_rows(Edata, Zdata) %>% 
  filter(!is.na(label))
ddtips<-subset(dd, dd$isTip =='TRUE')

ConnectionslongWithCoord<-merge(Connectionslong,
                                ddtips, by = "label",all.x =T)

DruEZorETree<-TreesSideBySide+ geom_line(aes(x, y, group=group, color=as.factor(Cluster)), ConnectionslongWithCoord,
                           linewidth =.4)+
  scale_color_manual(values=myclustercolors,name="Cluster")+
  guides(color = "none")
DruEZorETree

ggsave("ZorE_DruE_tree_with_connections.pdf",
       plot=DruEZorETree,
       path=FigDir,
       width=30,
       height=21,
       dpi=300,
       units="cm")

########################
########################
#Looking at ZorE most C-terminal domain
# read alignment (FASTA format assumed)
ZorEAln <- readAAStringSet("./data/ZorE2_mmseqs98.MUSCLE5.nuclease_motif.faa")
# extract positions 915–1040
ZorECterminus <- subseq(ZorEAln, start = 915, end = 1040)
# count gaps ("-") per sequence
Ctermgap_counts <- vcountPattern("-", ZorECterminus)
# build dataframe
Ctermgap_df <- data.frame(
  sequence = names(ZorECterminus),
  gap_count = as.numeric(Ctermgap_counts),
  row.names = NULL
)
Ctermgap_df$Cterm<-ifelse(Ctermgap_df$gap_count/(1040-915+1) > 0.7, 0,1)

Ctermgap_df <- Ctermgap_df %>%
  mutate(
    Protein_ID = sub(" .*", "", sequence)
  )
Ctermgap_df$line<-rep("ZorE C-terminus", length(Ctermgap_df$sequence))
###Let's now plot it on the tree

ZorETreeToSaveWithCTerm<-ZorEBasicTree+
  new_scale_fill()+
  geom_fruit(data=Ctermgap_df,
             geom = geom_tile,
             mapping =aes(y = Protein_ID,
                          x=line,
                          fill = as.factor(Cterm)),
             offset =0.01, pwidth =.07)+
  scale_fill_manual(values = c("white","#4393c3"),
                    na.translate = FALSE,
                    name = "ZorE C-terminus")+
  new_scale_fill()+
  geom_fruit(data=ZorEWithDru,
             geom = geom_tile,
             mapping =aes(y = TreeRepresentative.x,
                          x=Protein.x,
                          fill = Dru3),
             offset =0.015, pwidth =.1)+
  scale_fill_manual(values = c("#542788","#80cdc1","#dfc27d"),
                    na.value = "white",
                    na.translate = FALSE,
                    name = "Druantia III location")+
  new_scale_fill()+
  geom_fruit(data=ZorEClassCountsPerLeaf,
             geom = geom_col,
             mapping = aes(y=TreeRepresentative,
                           fill = CommonClass,
                           x = percent),
             offset =0.02, pwidth =.02)+
  scale_fill_manual(values=classcolors, name = "Class")+
  guides(fill = guide_legend(nrow = 4))+
  new_scale_fill()+
  geom_fruit(data=ZorEGenusCountsPerLeaf,
             geom = geom_col,
             mapping = aes(y=TreeRepresentative,
                           fill = TopGenus,
                           x = percent),
             offset =0.01, pwidth =.06)+
  scale_fill_manual(values=genuscolors, name = "Genus")+
  guides(fill = guide_legend(nrow = 4))+
  geom_fruit(data = ZorEHitsCountsPerLeaf,
             geom = geom_col,
             mapping = aes(y=TreeRepresentative,
                           x=logCount), fill= "#878787",
             pwidth =.4,
             axis.params=list(axis="x",
                              text.size=3,
                              line.size=.3),
             grid.params = list(size=.3,
                                alpha=.3))
ZorETreeToSaveWithCTerm     
ggsave("ZorE_tree_with_Tax_and_DruE_location_ZorECterm.pdf",
       plot=ZorETreeToSaveWithCTerm,
       path=FigDir,
       width=40,
       height=30,
       dpi=300,
       units="cm")
