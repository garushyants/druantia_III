library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(gggenes)
library(ggtree)
library(ggpubr)
library(aplot)
library(ape)
library(phytools)
library(scales)
library(ggnewscale)
library(tidyr)
library(this.path)
library(RColorBrewer)
library(ggh4x)
###################
mainpath<-paste0(dirname(this.path()),"/../")
setwd(mainpath)

Window = 20000 #the size of the neighborhood to extract

FigDir<-"figures/DruantiaIII_neighborhoods"

if (!dir.exists(FigDir)){
  dir.create(FigDir,recursive = TRUE)
} else {
  print("Directory already exists!")
}
###############
##In this version I am combining DefenseFinder and PADLOC outputs to get the more detailed picture
###############

#Starting by reading GFFs from representative genomes
GFFPath<-"./data/Druantia_neigborhoods/representative_genomes"
GFFFiles<-list.files(pattern="\\_genomic.gff.gz$",
                                       path = GFFPath)
setwd(paste(mainpath,GFFPath,sep="/"))
#Reading archived GFFs takes some time
GFFdata<-readr::read_tsv(GFFFiles, id="file_name", skip = 9, col_names = F)
#I select here both CDS and some non-coding genes
GFFdataCDS<-subset(GFFdata, !(GFFdata$X3 %in% c("gene","pseudogene","exon","region")) & !is.na(GFFdata$X2))
GFFdataCDS$GenomeID<-str_replace(GFFdataCDS$file_name,"_genomic.gff.gz","")

setwd(mainpath)

###Extracting fields
#######Select the required annotation fields
extractField<-function(pattern, column){
  r<-regexpr(pattern,column)
  out <- rep(NA,length(column))
  out[r!=-1] <- regmatches(column, r)
  out<-str_replace_all(out,";","")
  return(out)
}
#####
GFFdataCDS$X9<-paste0(GFFdataCDS$X9,";")#I do that in case the product is last record
GFFdataCDS$ID<-extractField('ID=[^-]+-([^;]+)\\;',GFFdataCDS$X9)
GFFdataCDS$ID<-str_split(GFFdataCDS$ID,"-",simplify = T, n = 5)[,2]
GFFdataCDS$gene<-extractField('\\;gene=([^;]+)\\;',GFFdataCDS$X9)
GFFdataCDS$gene<-str_replace(GFFdataCDS$gene,"gene=","")
GFFdataCDS$note<-extractField('\\;Note=([^;]+)\\;',GFFdataCDS$X9)
GFFdataCDS$note<-str_replace(GFFdataCDS$note,"Note=","")
GFFdataCDS$product<-extractField('\\;product=([^;]+)\\;',GFFdataCDS$X9)
GFFdataCDS$product<-str_replace(GFFdataCDS$product,"product=","")
names(GFFdataCDS)[2]<-"seqid"

#############################
###read DruE tree
DruE3tree <- read.tree("./data/phylogenetic_trees/DruE3_mmseqs98.trimmed.modelselection.IQTree.treefile")

#reroot
DruE3TreeMidRoot<-midpoint.root(DruE3tree)

###load clades
 DruE3clades<-read.csv("./data/phylogenetic_trees/DruE3_med_clade_3_final_clades.tsv", header=T, sep="\t")
 #get mrca
 DruE3NodeOfInterest<-DruE3clades %>% group_by(Cluster) %>%
   summarise(ClMRCA=getMRCA(DruE3TreeMidRoot, Sequence))

######################
#Reading In PADLOC data to find the borders
######################
##Read PADLOC files
PADLOCpath<-"./data/Druantia_neigborhoods/PADLOC_representative"
PADLOCselectedFiles<-list.files(pattern="\\.csv$",
                                recursive = T,
                                path = PADLOCpath)

setwd(paste0(mainpath,"/",PADLOCpath))
PADLOCdata<-readr::read_csv(PADLOCselectedFiles, id="file_name")
PADLOCdata$GenomeID<-str_replace(PADLOCdata$file_name,".csv","")
setwd(mainpath)

##subset Druantia III
PADLOCDruIII<-subset(PADLOCdata,
                     PADLOCdata$system == "druantia_type_III")
tmp<-PADLOCDruIII %>% group_by(seqid, system.number, GenomeID) %>%
  mutate(system.start = min(c(start,end)),
         system.end = max(c(start,end)))
DruEBorders<-subset(tmp, tmp$protein.name == "DruE3")
#get borders
DruEBorders$regionleft<-ifelse((DruEBorders$system.start-Window)>0, 
                                          DruEBorders$system.start-Window, 0)
DruEBorders$regionright<-DruEBorders$system.end+Window
#I have to remove the ones are DruE but not the ones on the tree
DruE3ReprSystemsOnTree<-subset(DruEBorders, DruEBorders$target.name %in% DruE3tree$tip.label)
#And then remove hits from one genome to get rid of duplicates
#WP_021699417.1 removing GCF_900455475.1
#WP_021702843.1 removing GCF_900455475.1

UniqDru3Borders<-subset(DruE3ReprSystemsOnTree,
                              DruE3ReprSystemsOnTree$GenomeID != 'GCF_900455475.1')
##########################
###############subset PADLOC within borders
PADLOCMainDF<-merge(PADLOCdata, UniqDru3Borders[,c(21,24,25,3)], 
              by = c("GenomeID",
                     "seqid"))

PADLOCMergeWithRegionsDF<-PADLOCMainDF %>% mutate(ROI = case_when(
  start >= regionleft & start <= regionright &
    end >= regionleft & end <= regionright ~ "Y"
))
##also remove pseudogene
PADLOCMergeWithRegionsDF<-subset(PADLOCMergeWithRegionsDF,
                           !startsWith(PADLOCMergeWithRegionsDF$target.name,"pseudo"))
###and partial domains (???)
PADLOCMergeWithRegionsDF$protein.name<-ifelse(PADLOCMergeWithRegionsDF$hmm.coverage < .5,
                                              paste(PADLOCMergeWithRegionsDF$protein.name,"partial",sep="_"),
                                              PADLOCMergeWithRegionsDF$protein.name)

PADLOCSelectedRegionsDF<-subset(PADLOCMergeWithRegionsDF,
                          PADLOCMergeWithRegionsDF$ROI =="Y")

###########################
##Subsetting GFF to the list of genes of interest

GFFMainDF<-merge(GFFdataCDS, UniqDru3Borders[c(3,21,24,25)],
              by = c("GenomeID",
                     "seqid"))
#this is essntial because otherwise selected regions are incorrect
GFFMainDF$X4<-as.numeric(GFFMainDF$X4)
GFFMainDF$X5<-as.numeric(GFFMainDF$X5)
GFFMergeWithRegionsDF<-GFFMainDF %>% mutate(ROI = case_when(
  X4 >= regionleft & X4 <= regionright &
    X5 >= regionleft & X5 <= regionright ~ "Y"
))
######################
#This is a table of all GFF records that are used for the downstream PFAM annotation
GFFSelectedRegionsDF<-subset(GFFMergeWithRegionsDF,
                          GFFMergeWithRegionsDF$ROI =="Y")

# ###Saving this to run the PFAM annotation later on
# GFFIDsToSave<-GFFSelectedRegionsDF[,c("GenomeID","ID")]
# write.table(GFFIDsToSave, file = "PFAM_GeneIDs_20000.tsv",
#             sep="\t",quote = F, row.names = F, col.names = F)

##########################################
##Reading PFAM hmmscan output
PFAMpath<-"data/Druantia_neigborhoods/PFAM_search_20251121/"
PFAMselectedFiles<-list.files(pattern="\\.csv$",
                                recursive = T,
                                path = PFAMpath)

setwd(paste0(mainpath,"/",PFAMpath))
PFAMdata<-readr::read_csv(PFAMselectedFiles, id="file_name")
PFAMdata$GenomeID<-str_replace(PFAMdata$file_name,"_pfamscan.csv","")
setwd(mainpath)

PFAMdata$PFAMID<-str_split_i(PFAMdata$hmm_acc,"\\.",1)
##
PFAMdatashort<-PFAMdata %>% group_by(GenomeID,seq_id) %>%
  summarise(PFAM_IDs = paste(sort(unique(PFAMID)), collapse = ";"),
            PHMM_names = paste(unique(hmm_name), collapse = ";"),
            .groups = "drop")

##Merge PFAM with GFF
GFFwPFAMofInterest<-merge(GFFSelectedRegionsDF[,c(1,2,5,6,7,9,12,13,14,15)],
                          PFAMdatashort,
                          by.x = c("GenomeID","ID"),
                          by.y = c("GenomeID","seq_id"),
                          all.x=T)

############################################
##Read DefenseFinder genes files
DefenseFinderpath<-"data/Druantia_neigborhoods/DefenseFinder_output/"
#this is way faster than doing a recursive search for files
DefenseFinderselectedFolders<-list.files(path = DefenseFinderpath)
DefenseFinderselectedFiles<-paste0(DefenseFinderselectedFolders, "/",
                                   DefenseFinderselectedFolders,
                                   "_protein.faa_defense_finder_genes.tsv")
setwd(paste0(mainpath,"/",DefenseFinderpath))

DefenseFinderdata<-readr::read_tsv(DefenseFinderselectedFiles, id="file_name")
DefenseFinderdata$GenomeID<-str_split_i(DefenseFinderdata$file_name,"/",1)

setwd(mainpath)

###
##Merge DefenseFinder results with GFF selected regions
DFwithCoord<-merge(GFFwPFAMofInterest,
                   DefenseFinderdata, by.x = c("GenomeID","ID"),
                   by.y = c("GenomeID","hit_id"), all.x =T)

DFwithCoordSh<-DFwithCoord[,c(1:12,15,35:37)]
##################################
##Add geNomad Annotations
#I do this annotation on the same selected proteoomes as PFAM search
GenomadAnnotation<-read.csv("./data/Druantia_neigborhoods/GenomadAnnotateResults_20260325.tsv",
                            sep="\t", header=F)
#Filter results
GenomadAnnotationFiltered<-subset(GenomadAnnotation,
                                  GenomadAnnotation$V3 < 0.001 & #E-value filter as in default genomad
                                    GenomadAnnotation$V6 > 0.75) #V6 is target coverage
##
GenomadAnnotationMeta<-read.csv("./data/Druantia_neigborhoods/genomad_metadata_v1.9.tsv", 
                                sep="\t")
GenomadAnnotationFull<-merge(GenomadAnnotationFiltered, 
                             GenomadAnnotationMeta,
                             by.x="V2",
                             by.y="MARKER",
                             all.x =T)
###Merge with other
collapse_safe <- function(x) {
  x <- sort(unique(na.omit(x)))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

GenomadAnnoCollapsed<-GenomadAnnotationFull %>% group_by(V1) %>%
  summarise(
    annotation_accessions        = collapse_safe(ANNOTATION_ACCESSIONS),
    annotation_descriptions      = collapse_safe(ANNOTATION_DESCRIPTION),
    genomad_hmm_ids               = collapse_safe(V2),
    genomad_specificity_class     = collapse_safe(SPECIFICITY_CLASS),
    .groups = "drop"
  )
########Merge with GFF and PFAM above
GFFwPGofInterest<-merge(DFwithCoordSh,
                        GenomadAnnoCollapsed,
                        by.x = c("ID"),
                        by.y = c("V1"),
                        all.x=T)
##################################
#####Create the final annotation dataframe
##Merge With PADLOC
AllFunctionalAnno<-merge(GFFwPGofInterest, 
                         PADLOCSelectedRegionsDF[,c(1:6,9,14:17)],
                         by.x=c("GenomeID","seqid","ID"),
                         by.y=c("GenomeID","seqid","target.name"), all=T)

AllFunctionalAnno$Gstart<-ifelse(is.na(AllFunctionalAnno$start),
                                 AllFunctionalAnno$X4, AllFunctionalAnno$start)
AllFunctionalAnno$Gend<-ifelse(is.na(AllFunctionalAnno$end),
                                 AllFunctionalAnno$X5, AllFunctionalAnno$end)
AllFunctionalAnno$Gstrand<-ifelse(is.na(AllFunctionalAnno$strand),
                               AllFunctionalAnno$X7, AllFunctionalAnno$strand)

#Do final round of merging annotations
AllFunctionalAnno$Agene<-ifelse(is.na(AllFunctionalAnno$protein.name),
                                ifelse(is.na(AllFunctionalAnno$gene_name),
                                       ifelse(is.na(AllFunctionalAnno$name),
                                              " ",
                                              AllFunctionalAnno$name),
                                       AllFunctionalAnno$gene_name),
                                AllFunctionalAnno$protein.name)
AllFunctionalAnno$Agene<-str_replace(AllFunctionalAnno$Agene,"Druantia__DruE_3","DruE3_DF")
AllFunctionalAnno$Agene<-str_replace(AllFunctionalAnno$Agene,"Druantia__DruH","DruH3")
AllFunctionalAnno$ggstrand<-ifelse(AllFunctionalAnno$Gstrand == '+',
                                   T,F)

AllFunctionalAnno<-AllFunctionalAnno[order(AllFunctionalAnno$seqid, AllFunctionalAnno$Gstart),]
###There are overlaping genes that I have to remove
AllFunctionalAnno$duplic<-ifelse((AllFunctionalAnno$Gend == lag(AllFunctionalAnno$Gend, default=-100) |
                                     AllFunctionalAnno$Gend == lead(AllFunctionalAnno$Gend, default=-100) | 
                                    AllFunctionalAnno$Gstart == lag(AllFunctionalAnno$Gstart, default=-100) |
                                       AllFunctionalAnno$Gstart == lead(AllFunctionalAnno$Gstart, default=-100)) &
                                   is.na(AllFunctionalAnno$protein.name),
                                  F,T)
AllFunctionalAnnoNoDupl<-subset(AllFunctionalAnno, AllFunctionalAnno$duplic)
#############################################################################
###############Merging annotations to show on the final plot
AllFunctionalAnnoNoDupl$HMM_Annot_IDs<-ifelse(!is.na(AllFunctionalAnnoNoDupl$annotation_accessions),
                                  AllFunctionalAnnoNoDupl$annotation_accessions,
                                  ifelse(!is.na(AllFunctionalAnnoNoDupl$PFAM_IDs),
                                         AllFunctionalAnnoNoDupl$PFAM_IDs,
                                         NA))
AllFunctionalAnnoNoDupl$product_description<-ifelse(!is.na(AllFunctionalAnnoNoDupl$protein.name),
                                                    AllFunctionalAnnoNoDupl$protein.name,
                                                    ifelse(!is.na(AllFunctionalAnnoNoDupl$gene_name),
                                                           AllFunctionalAnnoNoDupl$gene_name,
                                                           ifelse(!is.na(AllFunctionalAnnoNoDupl$annotation_descriptions),
                                                                  AllFunctionalAnnoNoDupl$annotation_descriptions,
                                                                  ifelse(is.na(AllFunctionalAnnoNoDupl$product),
                                                                         ifelse(is.na(AllFunctionalAnnoNoDupl$note),
                                                                                AllFunctionalAnnoNoDupl$X3,
                                                                                AllFunctionalAnnoNoDupl$note),
                                                                         AllFunctionalAnnoNoDupl$product))))
####created combined annotation
AllFunctionalAnnoNoDupl$PlotLabel<-ifelse(!is.na(AllFunctionalAnnoNoDupl$protein.name),
                                          AllFunctionalAnnoNoDupl$protein.name,
                                          ifelse(!is.na(AllFunctionalAnnoNoDupl$gene_name),
                                                 AllFunctionalAnnoNoDupl$gene_name,
                                                 ifelse(!is.na(AllFunctionalAnnoNoDupl$gene),
                                                        AllFunctionalAnnoNoDupl$gene,
                                                        ifelse(!is.na(AllFunctionalAnnoNoDupl$PFAM_IDs),
                                                               AllFunctionalAnnoNoDupl$PFAM_IDs,
                                                               ifelse(!is.na(AllFunctionalAnnoNoDupl$annotation_descriptions),
                                                                      AllFunctionalAnnoNoDupl$annotation_accessions,
                                                                      ifelse(AllFunctionalAnnoNoDupl$X3 != "CDS",
                                                                             ifelse(AllFunctionalAnnoNoDupl$X3 =="sequence_feature",
                                                                                    AllFunctionalAnnoNoDupl$note,
                                                                                    AllFunctionalAnnoNoDupl$X3),
                                                                             AllFunctionalAnnoNoDupl$product))))))

#Add viral things
AllFunctionalAnnoNoDupl$Integrase<-grepl("integrase|recombinase", AllFunctionalAnnoNoDupl$product_description, ignore.case=T)
#Add fill
AllFunctionalAnnoNoDupl$fill<-ifelse(AllFunctionalAnnoNoDupl$PlotLabel %in% c("DruH3","DruE3"),
                                     "druantia_type_III",
                                     ifelse(AllFunctionalAnnoNoDupl$PlotLabel == "ZorE2",
                                            "ZorE",
                                            ifelse(AllFunctionalAnnoNoDupl$PlotLabel == "WYL",
                                                   "WYL",
                                            ifelse(AllFunctionalAnnoNoDupl$X3 != "CDS",
                                                   ifelse(AllFunctionalAnnoNoDupl$X3 =="sequence_feature",
                                                          NA,
                                                          AllFunctionalAnnoNoDupl$X3),
                                                   ifelse(!is.na(AllFunctionalAnnoNoDupl$system.number),
                                                          ifelse(startsWith(AllFunctionalAnnoNoDupl$system,"PDC"),
                                                                 NA,AllFunctionalAnnoNoDupl$system),
                                                          ifelse(!is.na(AllFunctionalAnnoNoDupl$activity) &
                                                                   AllFunctionalAnnoNoDupl$activity == "Antidefense",
                                                                 "Antidefense",
                                                                 ifelse(!is.na(AllFunctionalAnnoNoDupl$subtype),AllFunctionalAnnoNoDupl$subtype,
                                                                 ifelse(AllFunctionalAnnoNoDupl$Integrase,
                                                                        "integrase",
                                                                        ifelse(!is.na(AllFunctionalAnnoNoDupl$genomad_specificity_class),
                                                                               ifelse(grepl("VV",AllFunctionalAnnoNoDupl$genomad_specificity_class),
                                                                                      "phage",
                                                                                      ifelse(grepl("PP",AllFunctionalAnnoNoDupl$genomad_specificity_class),
                                                                                             "plasmid",
                                                                                             NA)
                                                                               ),
                                                                               NA)
                                                                 ))))))))
AllFunctionalAnnoNoDupl$X3<-ifelse(is.na(AllFunctionalAnnoNoDupl$X3),
                                  "CDS",AllFunctionalAnnoNoDupl$X3)

##############################################################################
###Now add original protein ID
AllFuncWithMol<-merge(AllFunctionalAnnoNoDupl,UniqDru3Borders[,c(3,5,21)],
                      by=c("GenomeID","seqid"))
#Add clade ID
AllFuncWithMolClade<-merge(AllFuncWithMol,DruE3clades, by.x = "target.name", by.y ="Sequence")

#ARMADA adjustment
AllFuncWithMolClade <- AllFuncWithMolClade %>%
  group_by(GenomeID, seqid) %>%
  mutate(
    arm_match = case_when(
      (grepl("ZorD_partial", PlotLabel) |
        grepl("PF00176;PF00271", PlotLabel)) ~ "ArmA",
      grepl("REase_MTase_IIG", PlotLabel) ~ "ArmB",
      grepl("DruE_partial", PlotLabel) ~ "ArmC",
      grepl("PF00580;PF13361", PlotLabel) ~ "ArmD",
      TRUE ~ NA_character_
    ),
    
    # count DISTINCT ARM components
    arm_count = n_distinct(arm_match, na.rm = TRUE)
  ) %>%
  mutate(
    PlotLabel = ifelse(
      Cluster == 11 & !is.na(arm_match) & arm_count >= 2,
      arm_match,
      PlotLabel
    ),
    
    fill = ifelse(
      Cluster == 11 & !is.na(arm_match) & arm_count >= 2,
      "ARMADA",
      fill
    )
  ) %>%
  ungroup() %>%
  select(-arm_match, -arm_count)

# ##Saving annotations as a supplementary table
# SuppDataDFToSave<-AllFuncWithMolClade[,c(1:5,9,30:32,36:38,40,41)]
# colnames(SuppDataDFToSave)[c(1,3,4,5,7:9)]<-c("Tree_representative","ContigID","ProteinID","gene_type","start","end","strand")
# write.table(SuppDataDFToSave, "./data/Supplementary_table_2_DruantiaIII_genomic_contexts_20k.tsv",
#             sep="\t",
#             row.names = F,
#             quote = F)

############################################################################
###Let's draw figures per cluster for each neighborhood
#cluster colors
myclustercolors<-c("#a6cee3","#1f78b4","#b2df8a",
                   "#33a02c","#fb9a99","#e31a1c",
                   "#8c510a","#ff7f00","#cab2d6",
                   "#ffff99","#6a3d9a","#fdbf6f",
                   "#00441b","#dfc27d","#4d4d4d")
names(myclustercolors)<-unique(DruE3clades$Cluster)

##
TipColorDF<-data.frame(label=DruE3TreeMidRoot$tip.label)
TipColorDF$color<-ifelse(TipColorDF$label == "WP_020219138.1", "ECOR19", "")

TipColorScheme<-c("#cb181d", "#636363")
names(TipColorScheme)<-c("ECOR19","")

DrawClusterContexts<-function(cluster)
{
  #cluster<-11
  #subset tree
  Node<-subset(DruE3NodeOfInterest, DruE3NodeOfInterest$Cluster ==cluster)$ClMRCA
  subtree<-extract.clade(DruE3TreeMidRoot, node = Node)
  DruE3BasicWBt<-ggtree(subtree,
                        size=.5,
                        color="#636363") +
    #geom_tiplab()+
    geom_hilight(node=getMRCA(subtree,subtree$tip.label),fill = myclustercolors[cluster],
                 alpha=0.2,
                 extend=.35) +
    scale_fill_manual(values=myclustercolors,name="Cluster") +
    new_scale_fill()
  
  DruE3BasicTreePlot<- DruE3BasicWBt %<+%  TipColorDF +
    geom_tippoint(aes(colour=color),size=1)+
    scale_colour_manual(values=TipColorScheme, guide="none")+
    theme_tree() 
  
  leaf_order<-DruE3BasicTreePlot$data %>%
    filter(isTip) %>% arrange (y)
  
  ###
  GenesToPlotCluster<-subset(AllFuncWithMolClade,AllFuncWithMolClade$target.name %in% subtree$tip.label)
  
  ##
  essential_colors<-c("#377eb8","#4daf4a","#c51b8a","#c51b8a","#ec7014","#e41a1c","#c6dbef")
  names(essential_colors)<-c("druantia_type_III","ZorE","tmRNA","tRNA","phage","integrase","plasmid")
  allfillgroups<-unique(GenesToPlotCluster$fill)
  other_groups <- setdiff(allfillgroups, names(essential_colors))
  base_palette<-brewer.pal(9,"Set3")
  extended_palette <- colorRampPalette(base_palette)(length(other_groups))
  names(extended_palette)<-other_groups
  combinedcolors<-c(essential_colors,extended_palette)
  ##
  #trying to create compatible coordinates between samples
  #the idea is to align everything on DruE
  DruE3CoordFP<-subset(GenesToPlotCluster,
                       GenesToPlotCluster$protein.name == "DruE3")
  #DruE3CoordFP$middle<-DruE3CoordFP$Gstart + (DruE3CoordFP$Gend-DruE3CoordFP$Gstart)/2
  #I want to flip the ones that have other orientation
  DruE3CoordFP<-DruE3CoordFP%>%group_by(target.name) %>%
    mutate(middle = Gstart +(Gend-Gstart)/2)
  GenesToPlotNC<-merge(GenesToPlotCluster, DruE3CoordFP[,c(1,34,42)], by= "target.name")
  GenesToPlotNC$target.name<-factor(GenesToPlotNC$target.name, levels = leaf_order$label)
  GenesToPlotNC$pstart<-ifelse(GenesToPlotNC$ggstrand.y,
                               GenesToPlotNC$Gstart - GenesToPlotNC$middle,
                               GenesToPlotNC$middle - GenesToPlotNC$Gend)
  GenesToPlotNC$pend<-ifelse(GenesToPlotNC$ggstrand.y,
                             GenesToPlotNC$Gend - GenesToPlotNC$middle,
                             GenesToPlotNC$middle - GenesToPlotNC$Gstart)
  GenesToPlotNC$pstrand<-ifelse(GenesToPlotNC$ggstrand.y,
                                GenesToPlotNC$ggstrand.x,
                                !GenesToPlotNC$ggstrand.x)
  colnames(GenesToPlotNC)[1]<-"molecule"
  
  ###finally ploting
  GenesPlot<-ggplot(data = GenesToPlotNC,
                    aes(y =  molecule,
                        xmin = pstart,
                        xmax = pend))+ 
    geom_hline(aes(yintercept = molecule),
               linewidth =.5, color ="#bdbdbd")+
    geom_gene_arrow(aes(fill = fill,
                        forward=pstrand),
                    arrowhead_height = unit(4, "mm"),
                    arrow_body_height = unit(4, "mm"),
                    arrowhead_width = unit(1, "mm"))+
    geom_gene_label(aes(label=PlotLabel))+
    scale_fill_manual(values = combinedcolors, na.value="white", name ="")+
    theme_tree()+
    theme(legend.position = "right")
    
p<-GenesPlot %>% insert_left(DruE3BasicTreePlot, width =.1)

  
  ggsave(paste0("DruE_Cl",cluster,"_neighborhoods_and_HGT_20000.pdf"),
         plot=p,
         path = "./figures/DruantiaIII_neighborhoods/",
         width=50,
         height = ifelse(length(leaf_order$label)/2 <=10,
                         12,
                         length(leaf_order$label)/2),
         limitsize = F,
         units="cm",
         dpi=300)
}
####Save plots for all clusters
for(i in unique(AllFuncWithMolClade$Cluster))
{
  DrawClusterContexts(i)
}


###########################################
####Do summary and domain analysis
################################
#####Let's start with summary by cluster
#I use fill variable to describe essential classes

cluster_counts <- AllFuncWithMolClade %>%
  group_by(Cluster) %>%
  summarise(n_target_unique = n_distinct(target.name), .groups = "drop")

DrawCategoriesByCluster<-unique(AllFuncWithMolClade[,c("target.name","Cluster","GenomeID","seqid","fill")])%>% 
  filter(!is.na(fill)) %>%
  group_by(Cluster, fill) %>%
  summarise(
    fill_n = n(),
    .groups = "drop"
  ) %>%
  mutate(fill = gsub("_", " ", fill))%>%
  left_join(cluster_counts, by = "Cluster")
DrawCategoriesByCluster$percentage<-DrawCategoriesByCluster$fill_n*100/DrawCategoriesByCluster$n_target_unique
#for plot I select only the ones that occur in 20% of genomes per cluster
PercCutOff<-20
allfillgroups<-unique(subset(DrawCategoriesByCluster,
                             DrawCategoriesByCluster$percentage> PercCutOff)$fill)
DrawCategoriesByClusterForPlot<-subset(DrawCategoriesByCluster,
                                       DrawCategoriesByCluster$fill %in% allfillgroups)
DrawCategoriesByClusterForPlot$header<-paste0("Cluster ",DrawCategoriesByClusterForPlot$Cluster," (n=",
                                              DrawCategoriesByClusterForPlot$n_target_unique,
                                              ")")
DrawCategoriesByClusterForPlot$header<-factor(DrawCategoriesByClusterForPlot$header,
                                                 levels = unique(
                                                   DrawCategoriesByClusterForPlot[order(DrawCategoriesByClusterForPlot$Cluster),]$header))
###repeat the coloring scheme to match the ones on the long plots
essential_colors<-c("#377eb8","#4daf4a","#c51b8a","#c51b8a","#ec7014","#e41a1c","#c6dbef")
names(essential_colors)<-c("druantia type III","ZorE","tmRNA","tRNA","phage","integrase","plasmid")
other_groups <- setdiff(allfillgroups, names(essential_colors))
base_palette<-brewer.pal(9,"Set3")
extended_palette <- colorRampPalette(base_palette)(length(other_groups))
names(extended_palette)<-other_groups
Fillcombinedcolors<-c(essential_colors,extended_palette)
#plot summary for defense systems
SystemsSummary20k<-ggplot(data=DrawCategoriesByClusterForPlot, aes(x=fill, y =percentage, fill = fill))+
  geom_col()+
  facet_wrap(~header, ncol =5)+
  scale_fill_manual(values = Fillcombinedcolors, name ="")+
  ylab("% of genomes")+
  xlab("")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,hjust =1, vjust =.5, size=12),
        axis.text.y = element_text(size=12),
        strip.text = element_text(size = 12),
        legend.position = "none")
SystemsSummary20k  
#save summary figure
ggsave("DruantiaIII_neighborhoods_Summary_20k_column.pdf",
       plot=SystemsSummary20k,
       path = "./figures/DruantiaIII_neighborhoods/",
       width=35,
       height = 25,
       units="cm",
       dpi=300)

####################################################
#Draw representative neighborhoods for each cluster

RepresentativeGenomesList<-c("GCF_900044045.1","GCF_003730115.1","GCF_000785575.1",
                             "GCF_001005725.1",
                         "GCF_001557765.1","GCF_001586155.1","GCF_002174125.1","GCF_002875665.1",
                         "GCF_900168205.1",
                         "GCF_003047145.2",
                         "GCF_008806955.1",
                         "GCF_000380185.1",
                         "GCF_020694005.1", "GCF_014652135.1","GCF_020447145.1")

RepresentativeGenomes<-subset(AllFuncWithMolClade,
                              AllFuncWithMolClade$GenomeID %in% RepresentativeGenomesList)
RepresentativeGenomes$fill<-gsub("_"," ",RepresentativeGenomes$fill)
###now let's recalculate gene starts so everything aligns well
allfillgroupsrepr<-unique(RepresentativeGenomes$fill[!is.na(RepresentativeGenomes$fill)])
other_groupsrepr <- setdiff(allfillgroupsrepr, names(essential_colors))
base_palette<-brewer.pal(9,"Set3")
extended_paletterepr <- colorRampPalette(base_palette)(length(other_groupsrepr))
names(extended_paletterepr)<-other_groupsrepr
combinedcolorsrepr<-c(essential_colors,extended_paletterepr)

RepresentativeGenomes$header<-paste0("Cluster ", RepresentativeGenomes$Cluster)
RepresentativeGenomes$header<-factor(RepresentativeGenomes$header,
                                     levels = unique(
                                       RepresentativeGenomes[order(RepresentativeGenomes$Cluster),]$header))

###finally ploting
ReprGenesPlot<-ggplot(data = RepresentativeGenomes,
                  aes(y =  GenomeID,
                      xmin = Gstart,
                      xmax = Gend))+ 
  geom_hline(aes(yintercept = GenomeID),
             linewidth =.5, color="#bdbdbd")+
  geom_gene_arrow(aes(fill = fill,
                      forward=ggstrand),
                  arrowhead_height = unit(6, "mm"),
                  arrow_body_height = unit(6, "mm"),
                  arrowhead_width = unit(2, "mm"))+
  geom_gene_label(aes(label=PlotLabel))+
  facet_wrap2(~ header, scales = "free", ncol=1,
              strip = strip_themed(
                background_x = elem_list_rect(fill = scales::alpha(myclustercolors, 0.4))
              ))+
  ylab("")+
  scale_fill_manual(values = combinedcolorsrepr, 
                    na.value="white", name ="",
                    guide = guide_legend(ncol = 1),
                    breaks = names(combinedcolorsrepr))+
  theme_minimal()+
  theme(legend.position = "right",
        axis.text.y = element_text(size = 12, color ="black"))
ReprGenesPlot

#save
ggsave("DruantiaIII_neighborhoods_20k_representative_genomes.pdf",
       plot=ReprGenesPlot,
       path = "./figures/DruantiaIII_neighborhoods/",
       width=45,
       height = 32,
       units="cm",
       dpi=300)


###################
#Update domain story that I was showing before






#In reality I have to do that by clade and by not merging the domains
##Merge PFAM with GFF
PFAMinSelectedRegions<-merge(AllFuncWithMol[,c(1:2,4:8,11)],
                             PFAMdata[,c(2:18)],
                             by.x =  c("GenomeID","target.name"), by.y = c("GenomeID","seq_id"), all.x = T)
PFAMinSelectedRegionsWClades<-merge(PFAMinSelectedRegions,
                                    DruE3clades, by.x = "target.name", by.y = "Sequence", all.x = T)

DruEClusterSizes<-DruE3clades |> group_by(Cluster) |> summarise(ClSize = n())
#Removing domains found in DruE3 and DruH3, and only keeping neighboring ones
PFAMinSelectedRegionsWCladesNoDru<-subset(PFAMinSelectedRegionsWClades,
                                          PFAMinSelectedRegionsWClades$system != "druantia_type_III")
DomainRankingsByClade<-PFAMinSelectedRegionsWCladesNoDru |> group_by(Cluster,hmm_name,PFAMID, clan) |>
  summarise (n = n()) |>
  arrange(desc(n)) 

DomainRankingsByCladeNoNA<-DomainRankingsByClade  %>% drop_na()

DomainRankingsByCladeWSize<-merge(DomainRankingsByCladeNoNA, DruEClusterSizes, by = "Cluster", all.x =T)

DruEMostCommonDomainsByClade<-subset(DomainRankingsByCladeWSize, DomainRankingsByCladeWSize$n > DomainRankingsByCladeWSize$ClSize*.4)
###I want to add empty clusters
DruEMostCommonDomainsByCladeAllClusters<-merge(DruEMostCommonDomainsByClade, DruEClusterSizes,
                                               by=c("Cluster","ClSize"), all.y = T)
DruEMostCommonDomainsByCladeAllClusters$perc<-DruEMostCommonDomainsByCladeAllClusters$n*100/DruEMostCommonDomainsByCladeAllClusters$ClSize
DruEMostCommonDomainsByCladeAllClusters$hmm_name[is.na(DruEMostCommonDomainsByCladeAllClusters$hmm_name)] <- ""
DruEMostCommonDomainsByCladeAllClusters$ClusterHeader<-paste0("Clade ",DruEMostCommonDomainsByCladeAllClusters$Cluster, " n=", 
                                                              DruEMostCommonDomainsByCladeAllClusters$ClSize)
DruEMostCommonDomainsByCladeAllClusters$ClusterHeader<-factor(DruEMostCommonDomainsByCladeAllClusters$ClusterHeader,
                                                              levels= unique(DruEMostCommonDomainsByCladeAllClusters$ClusterHeader))

#Plot  

DruantiaIIIDomainsInNeib<-ggplot(data=DruEMostCommonDomainsByCladeAllClusters,
                                 aes(x = hmm_name, y = perc, fill = clan)) +
  geom_col() +
  facet_wrap(~ ClusterHeader, ncol =5)+#, scales = "free_y")+
  theme_classic()+
  ylab("%")+
  theme(axis.text.x = element_text(angle =90, vjust =.5, hjust =1))+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
                             "#fb9a99","#e31a1c","#fdbf6f","#ff7f00",
                             "#cab2d6","#6a3d9a"), na.translate=FALSE)

DruantiaIIIDomainsInNeib
ggsave("DruE_domains_in_neighborhood_40p_20000.pdf",
       path=FiguresAndDataFolder,
       plot=DruantiaIIIDomainsInNeib,
       width=24, height=15,
       limitsize = F,
       units="cm",
       dpi=300)
####Also I want to figure out how many of those has some kinds of tRNAs nearby...
####For the ones that I have it, I will pull all the genomes to look at the gene contexts and try to reconstruct whole the cargo that is there
NonCDSgenesIntegrase<-subset(GenesToPlotNC, GenesToPlotNC$X3 !="CDS" |
                               GenesToPlotNC$phage == "Integrase")

NonCDSgenesWClade<-merge(NonCDSgenesIntegrase, DruE3clades, by.x = "molecule", by.y = "Sequence")
NonCDSgenesWCladeCounts<-merge(NonCDSgenesWClade, DruEClusterSizes, by = "Cluster")
NonCDSgenesWCladeCounts$ClusterHeader<-paste0("Clade ",NonCDSgenesWCladeCounts$Cluster, " n=", 
                                              NonCDSgenesWCladeCounts$ClSize)
NonCDSgenesWCladeCounts$ClusterHeader<-factor(NonCDSgenesWCladeCounts$ClusterHeader,
                                              levels= unique(NonCDSgenesWCladeCounts$ClusterHeader))
NonCDSgenesWCladeCounts$GeneLabels<-ifelse(NonCDSgenesWCladeCounts$X3 == "riboswitch",
                                           "riboswitch",
                                           ifelse(NonCDSgenesWCladeCounts$X3 =="CDS",
                                                  "Integrase",
                                                  NonCDSgenesWCladeCounts$Agene)
)
NonCDSgenesWCladeCounts$GeneLabelsSimp<-ifelse(NonCDSgenesWCladeCounts$X3 == "sequence_feature",
                                               NonCDSgenesWCladeCounts$Agene,
                                               ifelse(NonCDSgenesWCladeCounts$X3 =="CDS",
                                                      "Integrase",
                                                      NonCDSgenesWCladeCounts$X3)
)
#plot
DruIIIMobilitySimp<-ggplot(NonCDSgenesWCladeCounts, aes(x = GeneLabelsSimp, fill = color))+
  geom_histogram(stat = "count") +
  facet_wrap(~ClusterHeader, scales = "free_y")+
  theme_classic()+
  ylab("count")+
  xlab("")+
  theme(axis.text.x = element_text(angle =90, vjust =.5, hjust =1))
DruIIIMobilitySimp
ggsave("DruIII_integrase_RNAs_in_neighborhood_simple_20000.pdf",
       path=FiguresAndDataFolder,
       plot=DruIIIMobilitySimp,
       width=24, height=15,
       limitsize = F,
       units="cm",
       dpi=300)
###So indeed there are more promising classes aside from Clade 11
###Clade1 (has tRNAs and tmRNAs),Clade2, Clade6 (t and tm), 
###Clade12, Clade13 (also has rRNAs), Clade14 (has tRNAs in less than half cases)
###Clade1 & 12 are less spread out and had mostly representatives from Vibrio
###Clade2 is also multiple species, so it might be interesting
###Clade14 can be interesting because it is soil bacteria, and they might have something very different

DruIIIMobility<-ggplot(NonCDSgenesWCladeCounts, aes(x = GeneLabels, fill = color))+
  geom_histogram(stat = "count") +
  facet_wrap(~ClusterHeader, scales = "free_y")+
  theme_classic()+
  ylab("count")+
  xlab("")+
  theme(axis.text.x = element_text(angle =90, vjust =.5, hjust =1),
        legend.position = "bottom")
DruIIIMobility
ggsave("DruIII_integrase_RNAs_in_neighborhood_20000.pdf",
       path=FiguresAndDataFolder,
       plot=DruIIIMobility,
       width=35, height=19,
       limitsize = F,
       units="cm",
       dpi=300)














































# AllFuncWithMol<-merge(AllFunctionalAnnoNoDupl[,c(1,2,3,9,17,18,23:27)], UniqDru3Borders[,c(3,5,21)],
#       by=c("GenomeID","seqid"))
# GenesToPlot<-AllFuncWithMol[,c(11,1:10)]
# names(GenesToPlot)[1]<-"molecule"
# GenesToPlot$fill<-ifelse(is.na(GenesToPlot$system), GenesToPlot$phage,
#                          ifelse(startsWith(GenesToPlot$system,"PD"),
#                                 GenesToPlot$phage,
#                                 GenesToPlot$system))
# 
# GenesColors<-hue_pal()(78)#rainbow(69)
# GenesColors<-append(GenesColors, "#525252", after=2)
# allSysNames<-unique(GenesToPlot$fill[!is.na(GenesToPlot$fill)])
# names(GenesColors)<-allSysNames[order(allSysNames)]
# 
# 
# 
# ##################################################################
# #####Getting all other types of INFO for plotting
# ##################################################################
# 
# ############Taxonomic distribution
# DruE3padlocDf<-read.csv("../20250423_Druantia3_whole_systems/DruE3_padloc20_refseq.nopseudo.cutoff07.withDruH.csv",
#                         header=F)
# #adjust genome IDs
# DruE3padlocDf$GenomeID<-sapply(strsplit(DruE3padlocDf$V1, ".csv:", fixed=TRUE), 
#                                head, 1)
# ##get info on clusters
# DruE3Clusters<-read.csv("../20250423_Druantia3_whole_systems/DruE3_mmseqs98_cluster.withref.tsv",
#                         header=F, sep="\t")
# 
# 
# ######Read assembly summary
# AssemblySummary<-read.csv("../20250423_Druantia3_whole_systems/Dru3_assembly_summary.tsv", 
#                           header = F, sep="\t",quote = "", 
#                           row.names = NULL, 
#                           stringsAsFactors = FALSE)
# 
# ##get representative genomes
# Dru3ReprGenomes<-read.csv("../20250424_DruantiIII_neighborhoods_genomad/Dru3_representative_genomes.txt", header=F)
# 
# 
# ####Adding cluster info to padloc df
# DruE3padlocDfWithCl<-merge(DruE3padlocDf,DruE3Clusters, by.x ="V4", by.y="V2", all.x =T)
# 
# AssemblySummary$genus<-word(AssemblySummary$V8,1)
# 
# CommonGenus<-AssemblySummary %>% group_by(genus) %>%
#   count()
# CommonGenus<-CommonGenus[order(-CommonGenus$n),]
# sum(pull(CommonGenus[c(1:7),2]))/length(AssemblySummary$V1)
# TopGenus<-CommonGenus$genus[1:7]
# AssemblySummary$TopOnly<-ifelse(AssemblySummary$genus %in% TopGenus, AssemblySummary$genus, "Else")
# 
# #picking colors
# TreeColors<-c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17","#969696")
# names(TreeColors)<-c(TopGenus,"Else")
# 
# ##########
# DruE3MetadataPADLOCLong<-merge(DruE3padlocDfWithCl[,c(20,21,1,4,7,8)],
#                                AssemblySummary, by.y="V1", by.x="GenomeID")
# 
# DruE3MetadataPADLOCbyWp<-DruE3MetadataPADLOCLong %>%
#   group_by(V1.y,TopOnly)%>%
#   count()
# names(DruE3MetadataPADLOCbyWp)<-c("molecule","Genus","Count")
# 
# #adding a small value so it will be visually seen which genus is there
# DruE3MetadataPADLOCbyWp$log10Count<-0.1+log10(DruE3MetadataPADLOCbyWp$Count)
# DruE3MetadataPADLOCbyWp$Genus<-factor(DruE3MetadataPADLOCbyWp$Genus, 
#                                       levels=c(TopGenus, 'Else'))
# 
# # #################################################################
# # ########################################
# # ###Plotting Genomad data on the DruE phylogenetic tree
# # DruE3padlocDfWithClRepr<-subset(DruE3padlocDfWithCl,
# #                                 DruE3padlocDfWithCl$GenomeID %in% Dru3ReprGenomes$V1)
# # 
# # DruE3padlocDfWithClRepr$SystemID<-sapply(strsplit(DruE3padlocDfWithClRepr$V1.x, ".csv:", fixed=TRUE), 
# #                                          tail, 1)
# # ##Plasmids
# # pathtogenomadfolder<-"../20250424_DruantiIII_neighborhoods_genomad/genomad_output/"
# # filelistgenomadplasmid = list.files(pattern="\\plasmid_summary.tsv$",
# #                                     recursive = T,
# #                                     path = pathtogenomadfolder)
# # setwd(paste0(mainpath,"/",pathtogenomadfolder))
# # GenomadResultsPlasmids<-readr::read_tsv(filelistgenomadplasmid, id="file_name")
# # setwd(mainpath)
# # GenomadResultsPlasmids<-separate(data =  GenomadResultsPlasmids,
# #                                  col = file_name,
# #                                  into=c("GenomeID",NA,NA),
# #                                  sep="/")
# # colnames(GenomadResultsPlasmids)[2]<-"ID"
# # 
# # DruWithGenomadPlasmid<-merge(DruE3padlocDfWithClRepr,
# #                              GenomadResultsPlasmids,
# #                              by.x= c("GenomeID","V2"),
# #                              by.y= c("GenomeID","ID"),
# #                              all.x =T)
# # nrow(subset(DruWithGenomadPlasmid, !is.na(DruWithGenomadPlasmid$length)))
# # DruGenomadPlSh<-DruWithGenomadPlasmid[,c(1:3,5,8,23,27)]
# # colnames(DruGenomadPlSh)[7]<-"GenomadPlasmidScore"
# # 
# # #########################
# # #viruses
# # filelistgenomadvirus = list.files(pattern="\\virus_summary.tsv$",
# #                                   recursive = T,
# #                                   path = pathtogenomadfolder)
# # setwd(paste0(mainpath,"/",pathtogenomadfolder))
# # GenomadResultsViruses<-readr::read_tsv(filelistgenomadvirus, id="file_name")
# # setwd(mainpath)
# # 
# # GenomadResultsViruses<-separate(data = GenomadResultsViruses,
# #                                 col = file_name,
# #                                 into=c("GenomeID",NA,NA),
# #                                 sep="/")
# # GenomadResultsViruses<-separate(data = GenomadResultsViruses,
# #                                 col = seq_name,
# #                                 into=c("ID","Provirus"),
# #                                 sep="\\|", remove = F)
# # GenomadResultsViruses<-separate(data = GenomadResultsViruses,
# #                                 col = coordinates,
# #                                 into=c("Start","End"),
# #                                 sep="\\-")
# # 
# # #it is essential to convert to numbers here, because otherwise it is interpreted as string
# # GenomadResultsViruses$Start<-as.integer(GenomadResultsViruses$Start)
# # GenomadResultsViruses$End<-as.integer(GenomadResultsViruses$End)
# # 
# # GenomadResultsVirusesCons<-subset(GenomadResultsViruses,
# #                                   GenomadResultsViruses$virus_score > .8)
# # 
# # 
# # PreGenomadVirusesOnTmnContigs<-merge(DruE3padlocDfWithClRepr,
# #                                      GenomadResultsVirusesCons, by.x =c("GenomeID","V2"),
# #                                      by=c("GenomeID","ID"))
# # PreGenomadVirusesOnTmnContigs$InProphage<-ifelse((PreGenomadVirusesOnTmnContigs$V12 >= PreGenomadVirusesOnTmnContigs$Start) & 
# #                                                    (PreGenomadVirusesOnTmnContigs$V13 <= PreGenomadVirusesOnTmnContigs$End) &
# #                                                    (PreGenomadVirusesOnTmnContigs$V13 <= PreGenomadVirusesOnTmnContigs$End) & 
# #                                                    (PreGenomadVirusesOnTmnContigs$V13 >= PreGenomadVirusesOnTmnContigs$Start),
# #                                                  1,0)
# # 
# # GenomadVirusesOnDruContigs<-subset(PreGenomadVirusesOnTmnContigs,
# #                                    InProphage == 1)
# # 
# # 
# # ##merging plasmid and Viral data
# # DruVirPlas<-merge(DruGenomadPlSh,
# #                   GenomadVirusesOnDruContigs[,c(1:3,35,31,26:28)],
# #                   all.x =T,
# #                   by = c("GenomeID","V2","V4"))
# # #####
# # DruVirPlas$InProphage<-ifelse(DruVirPlas$virus_score>0, 1, 0)
# # DruVirPlas$InPlasmid<-ifelse(DruVirPlas$length>0, 4, 0)
# # MGEDataForPlot<-DruVirPlas[,c("V4","InPlasmid","InProphage")]
# # colnames(MGEDataForPlot)[2:3]<-c(1,2)
# # MGEDataForPlot[is.na(MGEDataForPlot)]<-0
# # MGEDataForPlotLong<-gather(MGEDataForPlot,
# #                            key = "Location", value = "Prediction", 2:3)
# # names(MGEDataForPlotLong)[1]<-'molecule'
# # MGEDataForPlotLong$Location<-as.integer(MGEDataForPlotLong$Location)
# 
# TipColorDF<-data.frame(label=DruE3TreeMidRoot$tip.label)
# TipColorDF$color<-ifelse(TipColorDF$label == "WP_020219138.1", "ECOR19", "")
# 
# TipColorScheme<-c("#cb181d", "#636363")
# names(TipColorScheme)<-c("ECOR19","")
# 
# #cluster colors
# myclustercolors<-c("#a6cee3","#1f78b4","#b2df8a",
#                    "#33a02c","#fb9a99","#e31a1c",
#                    "#8c510a","#ff7f00","#cab2d6",
#                    "#ffff99","#6a3d9a","#fdbf6f",
#                    "#00441b","#dfc27d","#4d4d4d")
# names(myclustercolors)<-unique(DruE3clades$Cluster)
# #################################################################
# #Missing from this plot is Pfam annotations by contig
# 
# #################################################################
# #I plot each cluster independently
# DrawClusterContexts<-function(cluster)
# {
#   #cluster<-11
#   #subset tree
#   Node<-subset(DruE3NodeOfInterest, DruE3NodeOfInterest$Cluster ==cluster)$ClMRCA
#   subtree<-extract.clade(DruE3TreeMidRoot, node = Node)
#   DruE3BasicWBt<-ggtree(subtree,
#                         size=1,
#                         color="#636363") +
#     geom_nodelab(size=3) +
#     geom_tiplab()+
#     scale_x_continuous(limits = c(0,3.6))+
#     geom_hilight(node=getMRCA(subtree,subtree$tip.label),fill = myclustercolors[cluster],
#                  alpha=0.2,
#                  extend=.05) +
#     scale_fill_manual(values=myclustercolors,name="Cluster") +
#     new_scale_fill()
#   
#   DruE3BasicTreePlot<- DruE3BasicWBt %<+%  TipColorDF +
#     geom_tippoint(aes(colour=color),size=3)+
#     scale_colour_manual(values=TipColorScheme, guide="none")+
#     theme_tree2() 
#   leaf_order<-DruE3BasicTreePlot$data %>%
#     filter(isTip) %>% arrange (y)
#   
#   MGEDataCluster<-subset(MGEDataForPlotLong,
#                          MGEDataForPlotLong$molecule %in% subtree$tip.label)
#   MGEDataCluster$molecule<-factor(MGEDataCluster$molecule,
#                                       levels=leaf_order$label)
#   HGTPLot<-ggplot(data = MGEDataCluster, 
#                   aes(y=molecule,
#                       x=Location,
#                       fill = as.character(Prediction),
#                       color= as.character(Prediction)))+
#     geom_tile()+
#     scale_color_manual(values=c("#ffffff","#bebada","#8dd3c7"), guide="none")+
#     scale_fill_manual(values=c("#ffffff","#bebada","#8dd3c7"), guide="none")+
#     theme_tree2()+
#     theme(legend.position = 'none')
#   ##
#   GenesToPlotCluster<-subset(GenesToPlot,GenesToPlot$molecule %in% subtree$tip.label)
#   GenesToPlotCluster$molecule<-factor(GenesToPlotCluster$molecule, levels = leaf_order$label)
#   #trying to create compatible coordinates between samples
#   #the idea is to align everything on DruE
#   DruE3CoordFP<-subset(GenesToPlotCluster,
#                        GenesToPlotCluster$protein.name == "DruE3")
#   #DruE3CoordFP$middle<-DruE3CoordFP$Gstart + (DruE3CoordFP$Gend-DruE3CoordFP$Gstart)/2
#   #I want to flip the ones that have other orientation
#   DruE3CoordFP<-DruE3CoordFP%>%group_by(molecule) %>%
#     mutate(middle = Gstart +(Gend-Gstart)/2)
#   GenesToPlotNC<-merge(GenesToPlotCluster, DruE3CoordFP[,c(1,11,13)], by= "molecule")
#   
#   GenesToPlotNC$pstart<-ifelse(GenesToPlotNC$ggstrand.y,
#                                GenesToPlotNC$Gstart - GenesToPlotNC$middle,
#                                GenesToPlotNC$middle - GenesToPlotNC$Gend)
#   GenesToPlotNC$pend<-ifelse(GenesToPlotNC$ggstrand.y,
#                              GenesToPlotNC$Gend - GenesToPlotNC$middle,
#                              GenesToPlotNC$middle - GenesToPlotNC$Gstart)
#   GenesToPlotNC$pstrand<-ifelse(GenesToPlotNC$ggstrand.y,
#                                 GenesToPlotNC$ggstrand.x,
#                                 !GenesToPlotNC$ggstrand.x)
#   
#   ###finally ploting
#   GenesPlot<-ggplot(data = GenesToPlotNC,
#                     aes(y =  molecule,
#                         xmin = pstart,
#                         xmax = pend))+ 
#     geom_hline(aes(yintercept = molecule),
#                linewidth =.5, color ="#bdbdbd")+
#     geom_gene_arrow(aes(fill = fill,
#                         forward=pstrand))+
#     geom_gene_label(aes(label=Agene))+
#     scale_fill_manual(values = GenesColors, na.value="white")+
#     theme_tree2()+
#     theme(legend.position = "right")
#   
#   p<- ggarrange(DruE3BasicTreePlot,
#                 HGTPLot,
#                 GenesPlot, nrow =1,
#                 widths = c(0.5,.05,1),
#                 legend = "none",
#                 align='hv')
#   p
#   
#   ggsave(paste0("DruE_Cl",cluster,"_neighborhoods_and_HGT_20000.pdf"),
#          plot=p,
#          path = "../20260211_DruE_neighborhood_by_cluster",
#          width=50,
#          height = ifelse(length(leaf_order$label)/2 <=10,
#                          12,
#                          length(leaf_order$label)/2),
#          limitsize = F,
#          units="cm",
#          dpi=300)
# }
# 
# DrawClusterContexts(11)
# DrawClusterContexts(13)
# DrawClusterContexts(2)
# DrawClusterContexts(8)
# #################################################################
# DruE3BasicWBt<-ggtree(DruE3TreeMidRoot,
#                            size=1,
#                            color="#636363") +
#   geom_nodelab(size=3) +
#   geom_tiplab()+
#   scale_x_continuous(limits = c(0,3.6))+
#   geom_hilight(data = DruE3NodeOfInterest,
#                mapping = aes(node=ClMRCA,fill = as.factor(Cluster)),
#                alpha=0.2,
#                extend=.05) +
#   scale_fill_manual(values=myclustercolors,name="Cluster") +
#   new_scale_fill()
# 
# DruE3BasicTreePlot<- DruE3BasicWBt %<+%  TipColorDF +
#   geom_tippoint(aes(colour=color),size=3)+
#   scale_colour_manual(values=TipColorScheme, guide="none")+
#   theme_tree2() 
# 
# ##get leaf order
# leaf_order<-DruE3BasicTreePlot$data %>%
#   filter(isTip) %>% arrange (y)
# ##
# ###
# #doing species plot
# DruE3MetadataPADLOCbyWp$molecule<-factor(DruE3MetadataPADLOCbyWp$molecule,
#                                          levels=leaf_order$label)
# 
# SpeciesPlot<-ggplot(data=DruE3MetadataPADLOCbyWp,
#                     aes(x= log10Count,
#                     y = molecule,
#                     fill = Genus))+
#   geom_bar(stat='identity')+
#   scale_fill_manual(values=TreeColors)+
#   theme_tree2()+
#   theme(legend.position = 'none')
# 
# ###
# #doing HGT plot
# MGEDataForPlotLong$molecule<-factor(MGEDataForPlotLong$molecule,
#                                          levels=leaf_order$label)
# HGTPLot<-ggplot(data = MGEDataForPlotLong, 
#                 aes(y=molecule,
#                     x=Location,
#                     fill = as.character(Prediction),
#                     color= as.character(Prediction)))+
#   geom_tile()+
#   scale_color_manual(values=c("#ffffff","#bebada","#8dd3c7"), guide="none")+
#   scale_fill_manual(values=c("#ffffff","#bebada","#8dd3c7"), guide="none")+
#   theme_tree2()+
#   theme(legend.position = 'none')
# 
# ###
# #creating genes plot
# GenesToPlot$molecule<-factor(GenesToPlot$molecule, levels = leaf_order$label)
# #trying to create compatible coordinates between samples
# #the idea is to align everything on DruE
# DruE3CoordFP<-subset(GenesToPlot,
#                      GenesToPlot$protein.name == "DruE3")
# #DruE3CoordFP$middle<-DruE3CoordFP$Gstart + (DruE3CoordFP$Gend-DruE3CoordFP$Gstart)/2
# #I want to flip the ones that have other orientation
# DruE3CoordFP<-DruE3CoordFP%>%group_by(molecule) %>%
#   mutate(middle = Gstart +(Gend-Gstart)/2)
# GenesToPlotNC<-merge(GenesToPlot, DruE3CoordFP[,c(1,11,13)], by= "molecule")
# 
# GenesToPlotNC$pstart<-ifelse(GenesToPlotNC$ggstrand.y,
#                              GenesToPlotNC$Gstart - GenesToPlotNC$middle,
#                              GenesToPlotNC$middle - GenesToPlotNC$Gend)
# GenesToPlotNC$pend<-ifelse(GenesToPlotNC$ggstrand.y,
#                           GenesToPlotNC$Gend - GenesToPlotNC$middle,
#                           GenesToPlotNC$middle - GenesToPlotNC$Gstart)
# GenesToPlotNC$pstrand<-ifelse(GenesToPlotNC$ggstrand.y,
#                               GenesToPlotNC$ggstrand.x,
#                               !GenesToPlotNC$ggstrand.x)
# 
# ###finally ploting
# GenesPlot<-ggplot(data = GenesToPlotNC,
#        aes(y =  molecule,
#            xmin = pstart,
#            xmax = pend))+ 
#   geom_gene_arrow(aes(fill = fill,
#                       forward=pstrand))+
#   geom_gene_label(aes(label=Agene))+
#   scale_fill_manual(values = GenesColors, na.value="white")+
#   theme_tree2()+
#   theme(legend.position = "none")
# 
# 
# ####################
# ##Arranging All plots together
# p<- ggarrange(DruE3BasicTreePlot,
#               SpeciesPlot,
#               HGTPLot,
#           GenesPlot, nrow =1,
#           widths = c(0.5,.3,.05,1),
#           legend = "none",
#           align='hv')
# 
# ######################
# ggsave("DruE_long_with_neighborhoods_and_HGT_20000.pdf",
#        plot=p,
#        path = "../20260211_DruE_neighborhood_by_cluster",
#        limitsize = FALSE,
#        width =190,
#        height=170,
#        units="cm",
#        dpi=1000)
#   
# ###############################
# ################################
# #The simple
# DomainRankings<-AllFunctionalAnnoNoDupl |> group_by(name) |>
#   summarise (n = n()) |>
#   arrange(desc(n))
# #methylases are common
# 
# #In reality I have to do that by clade and by not merging the domains
# ##Merge PFAM with GFF
# PFAMinSelectedRegions<-merge(AllFuncWithMol[,c(1:2,4:8,11)],
#                              PFAMdata[,c(2:18)],
#                              by.x =  c("GenomeID","target.name"), by.y = c("GenomeID","seq_id"), all.x = T)
# PFAMinSelectedRegionsWClades<-merge(PFAMinSelectedRegions,
#                                     DruE3clades, by.x = "target.name", by.y = "Sequence", all.x = T)
# 
# DruEClusterSizes<-DruE3clades |> group_by(Cluster) |> summarise(ClSize = n())
# #Removing domains found in DruE3 and DruH3, and only keeping neighboring ones
# PFAMinSelectedRegionsWCladesNoDru<-subset(PFAMinSelectedRegionsWClades,
#                                           PFAMinSelectedRegionsWClades$system != "druantia_type_III")
# DomainRankingsByClade<-PFAMinSelectedRegionsWCladesNoDru |> group_by(Cluster,hmm_name,PFAMID, clan) |>
#   summarise (n = n()) |>
#   arrange(desc(n)) 
# 
# DomainRankingsByCladeNoNA<-DomainRankingsByClade  %>% drop_na()
# 
# DomainRankingsByCladeWSize<-merge(DomainRankingsByCladeNoNA, DruEClusterSizes, by = "Cluster", all.x =T)
# 
# DruEMostCommonDomainsByClade<-subset(DomainRankingsByCladeWSize, DomainRankingsByCladeWSize$n > DomainRankingsByCladeWSize$ClSize*.4)
# ###I want to add empty clusters
# DruEMostCommonDomainsByCladeAllClusters<-merge(DruEMostCommonDomainsByClade, DruEClusterSizes,
#                                                by=c("Cluster","ClSize"), all.y = T)
# DruEMostCommonDomainsByCladeAllClusters$perc<-DruEMostCommonDomainsByCladeAllClusters$n*100/DruEMostCommonDomainsByCladeAllClusters$ClSize
# DruEMostCommonDomainsByCladeAllClusters$hmm_name[is.na(DruEMostCommonDomainsByCladeAllClusters$hmm_name)] <- ""
# DruEMostCommonDomainsByCladeAllClusters$ClusterHeader<-paste0("Clade ",DruEMostCommonDomainsByCladeAllClusters$Cluster, " n=", 
#                                                               DruEMostCommonDomainsByCladeAllClusters$ClSize)
# DruEMostCommonDomainsByCladeAllClusters$ClusterHeader<-factor(DruEMostCommonDomainsByCladeAllClusters$ClusterHeader,
#                                                               levels= unique(DruEMostCommonDomainsByCladeAllClusters$ClusterHeader))
# 
# #Plot  
# 
# DruantiaIIIDomainsInNeib<-ggplot(data=DruEMostCommonDomainsByCladeAllClusters,
#                                  aes(x = hmm_name, y = perc, fill = clan)) +
#   geom_col() +
#   facet_wrap(~ ClusterHeader, ncol =5)+#, scales = "free_y")+
#   theme_classic()+
#   ylab("%")+
#   theme(axis.text.x = element_text(angle =90, vjust =.5, hjust =1))+
#   scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
#                              "#fb9a99","#e31a1c","#fdbf6f","#ff7f00",
#                              "#cab2d6","#6a3d9a"), na.translate=FALSE)
# 
# DruantiaIIIDomainsInNeib
# ggsave("DruE_domains_in_neighborhood_40p_20000.pdf",
#        path=FiguresAndDataFolder,
#        plot=DruantiaIIIDomainsInNeib,
#        width=24, height=15,
#        limitsize = F,
#        units="cm",
#        dpi=300)
# ####Also I want to figure out how many of those has some kinds of tRNAs nearby...
# ####For the ones that I have it, I will pull all the genomes to look at the gene contexts and try to reconstruct whole the cargo that is there
# NonCDSgenesIntegrase<-subset(GenesToPlotNC, GenesToPlotNC$X3 !="CDS" |
#                                GenesToPlotNC$phage == "Integrase")
# 
# NonCDSgenesWClade<-merge(NonCDSgenesIntegrase, DruE3clades, by.x = "molecule", by.y = "Sequence")
# NonCDSgenesWCladeCounts<-merge(NonCDSgenesWClade, DruEClusterSizes, by = "Cluster")
# NonCDSgenesWCladeCounts$ClusterHeader<-paste0("Clade ",NonCDSgenesWCladeCounts$Cluster, " n=", 
#                                               NonCDSgenesWCladeCounts$ClSize)
# NonCDSgenesWCladeCounts$ClusterHeader<-factor(NonCDSgenesWCladeCounts$ClusterHeader,
#                                               levels= unique(NonCDSgenesWCladeCounts$ClusterHeader))
# NonCDSgenesWCladeCounts$GeneLabels<-ifelse(NonCDSgenesWCladeCounts$X3 == "riboswitch",
#                                            "riboswitch",
#                                            ifelse(NonCDSgenesWCladeCounts$X3 =="CDS",
#                                                   "Integrase",
#                                                   NonCDSgenesWCladeCounts$Agene)
# )
# NonCDSgenesWCladeCounts$GeneLabelsSimp<-ifelse(NonCDSgenesWCladeCounts$X3 == "sequence_feature",
#                                                NonCDSgenesWCladeCounts$Agene,
#                                                ifelse(NonCDSgenesWCladeCounts$X3 =="CDS",
#                                                       "Integrase",
#                                                       NonCDSgenesWCladeCounts$X3)
# )
# #plot
# DruIIIMobilitySimp<-ggplot(NonCDSgenesWCladeCounts, aes(x = GeneLabelsSimp, fill = color))+
#   geom_histogram(stat = "count") +
#   facet_wrap(~ClusterHeader, scales = "free_y")+
#   theme_classic()+
#   ylab("count")+
#   xlab("")+
#   theme(axis.text.x = element_text(angle =90, vjust =.5, hjust =1))
# DruIIIMobilitySimp
# ggsave("DruIII_integrase_RNAs_in_neighborhood_simple_20000.pdf",
#        path=FiguresAndDataFolder,
#        plot=DruIIIMobilitySimp,
#        width=24, height=15,
#        limitsize = F,
#        units="cm",
#        dpi=300)
# ###So indeed there are more promising classes aside from Clade 11
# ###Clade1 (has tRNAs and tmRNAs),Clade2, Clade6 (t and tm), 
# ###Clade12, Clade13 (also has rRNAs), Clade14 (has tRNAs in less than half cases)
# ###Clade1 & 12 are less spread out and had mostly representatives from Vibrio
# ###Clade2 is also multiple species, so it might be interesting
# ###Clade14 can be interesting because it is soil bacteria, and they might have something very different
# 
# DruIIIMobility<-ggplot(NonCDSgenesWCladeCounts, aes(x = GeneLabels, fill = color))+
#   geom_histogram(stat = "count") +
#   facet_wrap(~ClusterHeader, scales = "free_y")+
#   theme_classic()+
#   ylab("count")+
#   xlab("")+
#   theme(axis.text.x = element_text(angle =90, vjust =.5, hjust =1),
#         legend.position = "bottom")
# DruIIIMobility
# ggsave("DruIII_integrase_RNAs_in_neighborhood_20000.pdf",
#        path=FiguresAndDataFolder,
#        plot=DruIIIMobility,
#        width=35, height=19,
#        limitsize = F,
#        units="cm",
#        dpi=300)
# 
# 
# 
# 
# 
# 
# 
# 
#   
#   
# 
# 
# 
# 
# 
# 
# 
