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

FigDir<-"figures/phylogenetic_trees"

if (!dir.exists(FigDir)){
  dir.create(FigDir,recursive = TRUE)
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
###MMseqs clusters
DruE3MMseqs<-read.csv("./data/mmseqs/DruE3_mmseqs98_cluster.withref.tsv", sep ="\t", header = F)
DruH3MMseqs<-read.csv("./data/mmseqs/DruH3_mmseqs98_cluster.withref.tsv", sep ="\t", header = F)
MMseqs<-rbind(DruE3MMseqs,DruH3MMseqs)
names(MMseqs)<-c("TreeRepresentative","ProteinID")
###Read PADLOC data
DruE3PADLOC<-read.csv("./data/DruE3_padloc20_refseq.nopseudo.cutoff07.withDruH.csv",header = F)
DruH3PADLOC<-read.csv("./data/DruH3_padloc20_refseq.nopseudo.cutoff07.csv",header = F)
#
TransformPADLOCdata<-function(df){
  df$RefseqID<-sub("\\.csv.*", "", df$V1)
  dfessential<-df[,c(4,20,2,7,6,12:14)]
  names(dfessential)<-c("ProteinID","RefseqID","Contig","Protein","PADLOC_HMM","Start","End","Strand")
  return(dfessential)
}
#get simple data
DruE3PADLOCEssential<-TransformPADLOCdata(DruE3PADLOC)
DruH3PADLOCEssential<-TransformPADLOCdata(DruH3PADLOC)
PADLOCEssential<-rbind(DruE3PADLOCEssential,DruH3PADLOCEssential)
####################
###Assembly summary
AssemblySummary<-read.csv("./data/DruE3_assembly_summary_refseq_20250415.txt", sep="\t", header =F)
AssemblySummaryEssential<-AssemblySummary[,c(1,7,8,9,12)]
names(AssemblySummaryEssential)<-c("RefseqID","Taxid",
                                   "Name",
                                   "Strain","AssemblyStatus")
###Taxonomic information
TaxPath<-"./data/taxonomy_info_datasets/"
TaxselectedFiles<-list.files(pattern="\\.tsv$",
                              path = TaxPath)
setwd(paste0(mainpath,"/",TaxPath))
Taxdata<-readr::read_tsv(TaxselectedFiles, id="file_name")

setwd(mainpath)
TaxdataEssential<-unique(Taxdata[,c(2,11,12,14,16,18,20,22,24,26)])
names(TaxdataEssential)[1]<-"Taxid"
###Merge PADLOCdata with genome info 
PADLOCWithGenomeInfo<-merge(PADLOCEssential,AssemblySummaryEssential,
      by="RefseqID",
      all.x =T)
###Add full taxonomy
PADLOCwithTaxonomy<-merge(PADLOCWithGenomeInfo,
                          TaxdataEssential,
                          by="Taxid",
                          all.x = T)
###Add MMseqs info
DruantiaTaxMMseqs<-merge(PADLOCwithTaxonomy,
                         MMseqs,
                         by="ProteinID",
                         all.x = T)
###Representative genomes
RepresentativeGenomes<-read.csv("./data/Dru3_representative_genomes.txt", header=F)
names(RepresentativeGenomes)<-"RefseqID"
RepresentativeGenomes$representative_genome<-rep("Y",length(RepresentativeGenomes$RefseqID))
###
DruantiaTaxMMseqsRep<-merge(DruantiaTaxMMseqs,RepresentativeGenomes,
                            by ="RefseqID",
                            all.x =T)
##############################
#Read in the data on the clusters
DruE3CladesAssignments<-read.csv("./data/phylogenetic_trees/DruE3_med_clade_3_final_clades.tsv", sep="\t")
DruH3CladesAssignments<-read.csv("./data/phylogenetic_trees/DruH3_med_clade_2.4_final_clades.tsv", sep="\t")

CladesAssignments<-rbind(DruE3CladesAssignments,
                         DruH3CladesAssignments)

DruantiaAllInfo<-merge(DruantiaTaxMMseqsRep,
                       CladesAssignments,
                       by.x = "TreeRepresentative",
                       by.y = "Sequence",
                       all.x = T)
##################################################################
##Read Genomad output 
##Plasmids
pathtogenomadfolder<-"./data/genomad_output/"
filelistgenomadplasmid = list.files(pattern="\\plasmid_summary.tsv$",
                                    recursive = T,
                                    path = pathtogenomadfolder)
setwd(paste0(mainpath,"/",pathtogenomadfolder))
GenomadResultsPlasmids<-readr::read_tsv(filelistgenomadplasmid, id="file_name")
setwd(mainpath)
GenomadResultsPlasmids<-separate(data =  GenomadResultsPlasmids,
                                 col = file_name,
                                 into=c("RefseqID",NA,NA),
                                 sep="/")
colnames(GenomadResultsPlasmids)[2]<-"Contig"
GenomadPlasmidsSimple<-GenomadResultsPlasmids[c(1,2,7)]
#get Plasmids in data
GenomadPlasmidsSimpleIndata<-subset(GenomadPlasmidsSimple,
                                    GenomadPlasmidsSimple$RefseqID %in% DruantiaAllInfo$RefseqID &
                                      GenomadPlasmidsSimple$Contig %in% DruantiaAllInfo$Contig)
# 
# DruWithGenomadPlasmid<-merge(DruE3padlocDfWithClRepr,
#                              GenomadResultsPlasmids,
#                              by.x= c("GenomeID","V2"),
#                              by.y= c("GenomeID","ID"),
#                              all.x =T)
# nrow(subset(DruWithGenomadPlasmid, !is.na(DruWithGenomadPlasmid$length)))
# DruGenomadPlSh<-DruWithGenomadPlasmid[,c(1:3,5,8,23,27)]
# colnames(DruGenomadPlSh)[7]<-"GenomadPlasmidScore"

######
#extracting phage info
filelistgenomadvirus = list.files(pattern="\\virus_summary.tsv$",
                                  recursive = T,
                                  path = pathtogenomadfolder)
setwd(paste0(mainpath,"/",pathtogenomadfolder))
GenomadResultsViruses<-readr::read_tsv(filelistgenomadvirus, id="file_name")
setwd(mainpath)

GenomadResultsViruses<-separate(data = GenomadResultsViruses,
                                col = file_name,
                                into=c("RefseqID",NA,NA),
                                sep="/")
GenomadResultsViruses<-separate(data = GenomadResultsViruses,
                                col = seq_name,
                                into=c("Contig","Provirus"),
                                sep="\\|", remove = F)
GenomadResultsViruses<-separate(data = GenomadResultsViruses,
                                col = coordinates,
                                into=c("Start","End"),
                                sep="\\-")

GenomadResultsViruses$Start<-as.integer(GenomadResultsViruses$Start)
GenomadResultsViruses$End<-as.integer(GenomadResultsViruses$End)
#add coordinates where prophage is a whole contig
GenomadResultsViruses$Start<-ifelse(is.na(GenomadResultsViruses$Start),
                                    1,
                                    GenomadResultsViruses$Start)
GenomadResultsViruses$End<-ifelse(is.na(GenomadResultsViruses$End),
                                  as.integer(GenomadResultsViruses$length),
                                  GenomadResultsViruses$End)

GenomadResultsVirusesCons<-subset(GenomadResultsViruses,
                                  GenomadResultsViruses$virus_score > .8)
##now get intersection with data

PreGenomadVirusesInDru<-merge(GenomadResultsVirusesCons,
                                     DruantiaAllInfo[,c(1,2,3,5,8:10)], by =c("RefseqID","Contig"))

PreGenomadVirusesInDru$InProphage<-ifelse((PreGenomadVirusesInDru$Start.y >= PreGenomadVirusesInDru$Start.x) & 
                                                   (PreGenomadVirusesInDru$End.y <= PreGenomadVirusesInDru$End.x) &
                                                   (PreGenomadVirusesInDru$Start.y <= PreGenomadVirusesInDru$End.x) & 
                                                   (PreGenomadVirusesInDru$End.y >= PreGenomadVirusesInDru$Start.x),
                                                 1,0)

GenomadVirusesOnDruContigs<-subset(PreGenomadVirusesInDru,
                                   InProphage == 1)
GenomadVirusesSimple<-GenomadVirusesOnDruContigs[,c(1,2,16,11)]

####Add information on GEI from TreasureIsland#################
pathtotifolder<-"./data/treasureisland_results/"
filelistti = list.files(pattern="\\combined.csv$",
                                  path = pathtotifolder)
setwd(paste0(mainpath,"/",pathtotifolder))
#filter out empty files
filelisttifi<-filelistti[file.info(filelistti)$size > 0]
TiResults<-readr::read_csv(filelisttifi, id="file_name", col_names = F)
setwd(mainpath)
TiResults$RefseqID<-str_replace(TiResults$file_name,"_combined.csv","")
names(TiResults)[3:5]<-c("Contig","GEI_start","GEI_end")
PreTiResultsInDru<-merge(TiResults,
                              DruantiaAllInfo[,c(1,2,3,5,8:10)], by =c("RefseqID","Contig"))
#....and get intersection
PreTiResultsInDru$GEI<-ifelse((PreTiResultsInDru$Start >= PreTiResultsInDru$GEI_start) & 
                                  (PreTiResultsInDru$End <= PreTiResultsInDru$GEI_end) &
                                  (PreTiResultsInDru$Start <= PreTiResultsInDru$GEI_end) & 
                                  (PreTiResultsInDru$End >= PreTiResultsInDru$GEI_start),
                                "GEI",NA)
GEIOnDruContigs<-subset(PreTiResultsInDru,
                                   GEI == "GEI")
GEIOnDruContigsSimple<-GEIOnDruContigs[,c(1,2,8,13)]
###################################################
###Add MGE data to the main dataset
DruantiaAllWithPlasmids<-merge(DruantiaAllInfo,
                               GenomadPlasmidsSimpleIndata,
                               by = c("RefseqID","Contig"),
                               all.x =T)
DruantiaAllWithPlasmidsVirus<-merge(DruantiaAllWithPlasmids,
                          GenomadVirusesSimple,
                          by = c("RefseqID","Contig","TreeRepresentative"),
                          all.x = T)
DruantiaAllWithMGE<-merge(DruantiaAllWithPlasmidsVirus,
                          GEIOnDruContigsSimple,
                                    by = c("RefseqID","Contig","TreeRepresentative"),
                                    all.x = T)
# ####Save all metadata about dataset to file
# write.table(DruantiaAllWithMGE, file = "./data/Supplementary_table_1_dataset_info.tsv",
#             row.names = F, sep="\t", quote =F)
##################################################################
####Exploring gene order
DruHDruEOrder<-DruantiaAllWithMGE %>%
  group_by(RefseqID, Contig) %>%
  summarize(
    start_E = Start[Protein == "DruE3"][1],
    end_E   = End[Protein == "DruE3"][1],
    start_H = Start[Protein == "DruH3"][1],
    end_H   = End[Protein == "DruH3"][1],
    strand  = Strand[1],
    order = case_when(
      strand == "+" & start_E < start_H ~ "druEdruH",
      strand == "-" & start_H < start_E ~ "druEdruH",
      TRUE ~ "druHdruE"
    ),
    distance = case_when(
      strand == "+" ~ start_E - end_H +1,  
      strand == "-" ~ end_E - start_H +1,
      TRUE ~ NA_real_
    ),
    .groups = "drop"
  )
unique(DruHDruEOrder$order)
DruHDruEOrder<-subset(DruHDruEOrder, DruHDruEOrder$order !="druEdruH") #there is one case of incorrect positioning
#In all cases the order is preserved aside from one glitch in grouping, but
#checking of initial files suggest that the gene order is the same
#also check the intergenic distance
ggplot(DruHDruEOrder) +
  geom_histogram(aes(x=distance),
                 bins =300)+
  geom_vline(xintercept = mean(DruHDruEOrder$distance))+
  xlim(-500,500)
##################################################################

####Getting nodes for clusters for later vizualization

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
                   "#8c510a","#ff7f00","#cab2d6",
                   "#ffff99","#6a3d9a","#fdbf6f",
                   "#00441b","#dfc27d","#4d4d4d")
names(myclustercolors)<-unique(DruE3ClusterAncestryNodeOfInterest$Cluster)


DrawCircularTreeWithData<-function(tree,clusternodes, bootstraps)
{
  tree<-DruE3TreeMidRoot
  clusternodes<-DruE3ClusterAncestryNodeOfInterest
  bootstraps<- DruE3BootstrapValuesA80
  BasicTreePlot<-ggtree(tree,layout = 'fan', open.angle=10, 
                    size=0.2,
                    color="#636363")%<+% bootstraps +
    geom_nodepoint(aes(size = UFboot),
                   color = '#4292c6',
                   alpha=.3)+
    scale_size_continuous(range = c(0.01,1))+
    geom_treescale(y=1, x=3.2, fontsize=3, linesize=0.7, offset=1)
  TreeWithCladesSimple<-BasicTreePlot +
    geom_hilight(data = clusternodes,
                 mapping = aes(node=ClMRCA,
                               fill = as.factor(Cluster)),
                 alpha=0.2,
                 extend=.05)+
    scale_fill_manual(values=myclustercolors,name="Cluster")
  ##metadata for the tree
  Metadata<-subset(DruantiaAllWithMGE,
                   DruantiaAllWithMGE$TreeRepresentative %in% tree$tip.label)
  #############MGE data
  MGEdata<-subset(Metadata[,c("TreeRepresentative",
                              "plasmid_score",
                              "virus_score",
                              "GEI")],
                  !(is.na(Metadata$plasmid_score) &
                      is.na(Metadata$virus_score) & 
                      is.na(Metadata$GEI))
                  )
  MGEdata$Plasmid<-ifelse(!is.na(MGEdata$plasmid_score),"plasmid",NA)
  MGEdata$Prophage<-ifelse(!is.na(MGEdata$virus_score),"prophage",NA)
  MGEdatalong<-MGEdata[,c(1,4:6)] %>% pivot_longer(cols = -TreeRepresentative,
                                        names_to = "variable",
                                        values_to = "Location",
                                        values_drop_na = TRUE)
  ########################################
  ##get taxonomy for outer circles
  CommonGenusInSet<-Metadata %>% group_by(`Genus name`) %>%
    count()
  CommonGenusInSet<-CommonGenusInSet[order(-CommonGenusInSet$n),]
  
  sum(pull(CommonGenusInSet[c(1:12),2]))/length(Metadata$ProteinID)
  #I am taking 12 most common genus, because they together account for ~92% of all genomes, and tmn there found at least in >50 genomes per genus
  GenusToKeep<-c(pull(CommonGenusInSet[c(1:12),1]))
  Metadata$TopGenus<-ifelse(Metadata$`Genus name` %in% GenusToKeep,
                              Metadata$`Genus name`, "Other")
  
  HitsCountsPerLeaf<-Metadata %>%
    group_by(TreeRepresentative)%>%
    count()
  #transform to log10
  HitsCountsPerLeaf$logCount<-log10(HitsCountsPerLeaf$n)+.1
  
  TaxonomyCountsPerLeaf<-Metadata %>%
    group_by(TreeRepresentative, TopGenus) %>% 
    count(TopGenus) %>%
    group_by(TreeRepresentative) %>%
    mutate(percent = n / sum(n) * 100)
  TaxonomyCountsPerLeaf$TopGenus<-factor(TaxonomyCountsPerLeaf$TopGenus,
                                         levels=c(GenusToKeep,"Other"))
  
  CommonClassInSet<-Metadata %>% group_by(`Class name`) %>%
    count()%>%
    arrange(desc(n))
  ClassDf<-Metadata[,c("TreeRepresentative","Class name")]
  ClassDf$CommonClass<-ifelse(ClassDf$`Class name` %in% CommonClassInSet$`Class name`[1:5],
                              ClassDf$`Class name`,
                              "Other")
  #I pick all Classes with > 15 genomes
  ClassCountsPerLeaf<-ClassDf %>%
    group_by(TreeRepresentative, CommonClass) %>% 
    count(CommonClass) %>%
    group_by(TreeRepresentative) %>%
    mutate(percent = n / sum(n) * 100)
  unique(ClassCountsPerLeaf$CommonClass)
  ClassCountsPerLeaf$CommonClass<-factor(ClassCountsPerLeaf$CommonClass,
                                         levels=c(CommonClassInSet$`Class name`[1:5],
                                                  "Other"))
  ########################################
  #create colors
  genuscolors<-c("#e31a1c","#33a02c","#1f78b4","#6a3d9a",
                 "#ff7f00","#b15928","#ffff99","#fb9a99",
                 "#fdbf6f","#b2df8a","#a6cee3","#cab2d6", "#d9d9d9")
  names(genuscolors)<-c(GenusToKeep,"Other")
  
  classcolors<-c("#fbb4ae","#b3cde3","#ccebc5",
                 "#decbe4","#fed9a6","#d9d9d9")
  names(classcolors)<-c(CommonClassInSet$`Class name`[1:5],"Other")
  #########################################
  ####Vizualize all the info
  TreeWithMGE<-TreeWithCladesSimple +
    new_scale_fill()+
    geom_fruit(data=MGEdatalong,
               geom = geom_tile,
               mapping = aes(y=TreeRepresentative, x= Location,
                             fill = Location),
               color="#bdbdbd",
               pwidth=.04,
               offset = 0.1,
               axis.params=list(axis="x",
                                text.size=2,
                                text.angle=60,
                                vjust=0,
                                hjust=1))+
    scale_fill_manual(values=c("#6a3d9a","#238443","#0570b0"), guide="none", na.value = "#ffffff")
  TreeWithMGEAndTax<-TreeWithMGE +
    new_scale_fill()+
    geom_fruit(data=ClassCountsPerLeaf,
               geom = geom_col,
               mapping = aes(y=TreeRepresentative,
                             fill = CommonClass,
                             x = percent),
               offset =0.01, pwidth =.02)+
    scale_fill_manual(values=classcolors, name = "Class")+
    guides(fill = guide_legend(nrow = 4))+
    new_scale_fill()+
    geom_fruit(data=TaxonomyCountsPerLeaf,
               geom = geom_col,
               mapping = aes(y=TreeRepresentative,
                             fill = TopGenus,
                             x = percent),
               offset =0.01, pwidth =.1)+
    scale_fill_manual(values=genuscolors, name = "Genus")+
    guides(fill = guide_legend(nrow = 4))
  
  MainTreeToSave<-TreeWithMGEAndTax +
    geom_fruit(data = HitsCountsPerLeaf,
               geom = geom_col,
               mapping = aes(y=TreeRepresentative,
                             x=logCount), fill= "#878787",
               pwidth =.4,
               axis.params=list(axis="x",
                                text.size=3,
                                line.size=.3),
               grid.params = list(size=.3,
                                  alpha=.3))+
    theme(legend.position = "bottom")
  return(MainTreeToSave)
}
################################
#####Draw circular tree plots
DruE3CircularPlot<-DrawCircularTreeWithData(DruE3TreeMidRoot,
                                            DruE3ClusterAncestryNodeOfInterest,
                                            DruE3BootstrapValuesA80)
DruH3CircularPlot<-DrawCircularTreeWithData(DruH3TreeMidRoot,
                                            DruH3ClusterAncestryNodeOfInterest,
                                            DruH3BootstrapValuesA80)
####save circular plots
ggsave("DruE_tree_with_MGE_and_Tax.pdf",
       plot=DruE3CircularPlot,
       path=FigDir,
       width=30,
       height=30,
       dpi=300,
       units="cm")

ggsave("DruH_tree_with_MGE_and_Tax.pdf",
       plot=DruH3CircularPlot,
       path=FigDir,
       width=30,
       height=30,
       dpi=300,
       units="cm")


#########################################
#########################################
###Draw interactions between DruE and DruH

InfoSubset<-DruantiaAllInfo[,c("RefseqID","Contig","Protein",
                   "TreeRepresentative")]
DruEInfoSubset<-subset(InfoSubset, InfoSubset$Protein == "DruE3")
DruHInfoSubset<-subset(InfoSubset, InfoSubset$Protein == "DruH3")

DruParallel<-merge(DruEInfoSubset,
                   DruHInfoSubset,
                   by = c("RefseqID","Contig"))
names(DruParallel)[c(4,6)]<-c("druE","druH")
#get all unique pairs
DruAlluniquepairs<-unique(DruParallel[c(4,6)])
DruAlluniquepairs$group<-seq(1,length(DruAlluniquepairs$druE))
DruAllpairslong<-DruAlluniquepairs %>%
  pivot_longer(cols = starts_with("dru"),
               names_to = "gene",
               values_to= "label")

#Draw basic trees
#DruE
DruETree<-ggtree(DruE3TreeMidRoot,
                 size=0.2,
                 color="#636363")

DruETreeWithCl<-DruETree +
  geom_hilight(data = DruE3ClusterAncestryNodeOfInterest,
               mapping = aes(node=ClMRCA,
                             fill = as.factor(Cluster)),
               alpha=0.2)+
  scale_fill_manual(values=myclustercolors,name="Cluster")
DruETreeWithCl

#DruH
DruHTree<-ggtree(DruH3TreeMidRoot,
                 size=0.2,
                 color="#636363") +
  geom_nodepoint(aes(size = UFboot),
                 color = '#4292c6',
                 alpha=.3)+
  scale_size_continuous(range = c(0.01,1))

#get plot data
Edata <- DruETreeWithCl$data
Hdata <- DruHTree$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
Hdata$x <- max(Hdata$x) - Hdata$x + max(Edata$x) + 1
#Combine trees together
TreesSideBySide <- DruETreeWithCl + geom_tree(data=Hdata,
                                 size=0.2,
                                 color="#636363") 
###get coordinates of all tips
dd <- bind_rows(Edata, Hdata) %>% 
  filter(!is.na(label))
ddtips<-subset(dd, dd$isTip =='TRUE')
#merge with all possible pairs
DruAllpairslongWithCoord<-merge(DruAllpairslong,
      ddtips, by = "label",all.x =T)

#Get nodes that connect more than with one node
DruAllPairs<-DruAllpairslong %>% pivot_wider(id_cols = group,
                                             names_from = gene,
                                             values_from = label)
DruEcounts<-as.data.frame(table(DruAllPairs$druE)) %>% arrange(desc(Freq))
DruHcounts<-as.data.frame(table(DruAllPairs$druH)) %>% arrange(desc(Freq))
EdgeCountsAll<-rbind(DruEcounts,DruHcounts)
names(EdgeCountsAll)<-c("label","n")

DruAllpairslongWithCoordCounts<-merge(DruAllpairslongWithCoord,
                                      EdgeCountsAll,
                                      by="label")
DruAllpairslongWithCoordCounts<-DruAllpairslongWithCoordCounts %>% group_by(group) %>%
  mutate(vertex_degree = max(n))


VertexDegreeDf<-DruAllpairslongWithCoordCounts[,c(1:3,12)] %>%
  pivot_wider(id_cols = group,values_from = c(label,n),names_from = gene) %>% arrange(desc(n_druE),desc(n_druH))
# #Saving the table with vertex degrees fo the furture use
# write.table(VertexDegreeDf, file = "./data/Dru_vertex_degree_summary.tsv",
#             sep="\t", row.names = F, quote = F)

###draw final tree with connections
DruEvsDruH<-TreesSideBySide + geom_line(aes(x, y, group=group, color=vertex_degree), DruAllpairslongWithCoordCounts,
                                        alpha =.8, linewidth =.1) +
  scale_color_gradient2(low = 'lightblue', mid = 'lightblue', high = '#d73027', midpoint =1,
                        breaks = seq(1,13,3))
DruEvsDruH

#save
ggsave(filename="DruE3_vs_DruH3_tree_with_connections.pdf",
       plot=DruEvsDruH,
       path=FigDir,
       dpi=300,
       width = 45,
       height =20,
       units = "cm") 

DruPairsHighDegree<-subset(DruAllpairslongWithCoordCounts, DruAllpairslongWithCoordCounts$vertex_degree >9)
DruEvsDruHHighDegree<-TreesSideBySide + geom_line(aes(x, y, group=group, color=vertex_degree), DruPairsHighDegree,
                                        alpha =.8, linewidth =.1) +
  scale_color_gradient2(low = 'lightblue', mid = 'lightblue', high = '#d73027', midpoint =1,
                        breaks = seq(1,13,3))
DruEvsDruHHighDegree

ggsave(filename="DruE3_vs_DruH3_tree_with_high_degree_connections.pdf",
       plot=DruEvsDruHHighDegree,
       path=FigDir,
       dpi=300,
       width = 45,
       height =20,
       units = "cm") 