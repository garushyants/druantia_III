library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(gggenes)
library(ggtree)
library(ggpubr)
library(ape)
library(phytools)
library(scales)
#library(ggnewscale)
library(tidyr)

mainpath<-"/Volumes/garushyants/druantia_zorya/20250506_DruIII_defense_neighborhoods"
setwd(mainpath)
###############
##
##In this version I am combining DefenseFinder and PADLOC outputs to get the more detailed picture
##I don't think that in reality it improves things by much, but I prefer to have it all together
##And also add info from general annotation in GFFs
##2025/09/25 And more importantly I add clades that I selected from DruE/DruH trees
###############

##The important thing that I need to get genome coordinates
##And genes of interest for downstream PFAM analysis
##And for that I need to read info from GFF files
GFFPath<-"./representative_genomes"
GFFFiles<-list.files(pattern="\\_genomic.gff.gz$",
                                       path = GFFPath)
setwd(paste(mainpath,GFFPath,sep="/"))
#Reading archived GFFs takes some time
GFFdata<-readr::read_tsv(GFFFiles, id="file_name", skip = 9, col_names = F)
GFFdataCDS<-subset(GFFdata, GFFdata$X3 == "CDS")
GFFdataCDS$GenomeID<-str_replace(GFFdataCDS$file_name,"_genomic.gff.gz","")
setwd(mainpath)

#####Right way to extract fields from GFF
extractField<-function(pattern, column){
  #pattern<-'\\;product=([^;]+)\\;'
  #column<-TmnNeigGFF$V10
  r<-regexpr(pattern,column)
  out <- rep(NA,length(column))
  out[r!=-1] <- regmatches(column, r)
  out<-str_replace_all(out,";","")
  return(out)
}
#####
GFFdataCDS$ID<-extractField('ID=cds-([^;]+)\\;',GFFdataCDS$X9)
GFFdataCDS$ID<-str_replace(GFFdataCDS$ID,"ID=cds-","")
GFFdataCDS$product<-extractField('\\;product=([^;]+)\\;',GFFdataCDS$X9)
GFFdataCDS$product<-str_replace(GFFdataCDS$product,"product=","")
names(GFFdataCDS)[2]<-"seqid"

#############################
###read DruE tree
DruE3tree <- read.tree("../20250423_Druantia3_whole_systems/DruE3_mmseqs98.trimmed.modelselection.IQTree.treefile")

#reroot
DruE3TreeMidRoot<-midpoint.root(DruE3tree)

###load clades
 DruE3clades<-read.csv("../20250717_DruE_DruH_treecluster/DruE3_med_clade_3_final_clades.tsv", header=T, sep="\t")
 #get mrca
 DruE3NodeOfInterest<-DruE3clades %>% group_by(Cluster) %>%
   summarise(ClMRCA=getMRCA(DruE3TreeMidRoot, Sequence))

######################
#Reading In PADLOC data to find the borders
######################
##Read PADLOC files
PADLOCpath<-"./PADLOC_representative"
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
Window = 20000
DruEBorders$regionleft<-ifelse((DruEBorders$system.start-Window)>0, 
                                          DruEBorders$system.start-Window, 0)
DruEBorders$regionright<-DruEBorders$system.end+Window
#I have to remove the ones are DruE but not the ones on the tree
DruE3ReprSystemsOnTree<-subset(DruEBorders, DruEBorders$target.name %in% DruE3tree$tip.label)
#And then remove hits from one genome to get rid of duplicates
#I have duplicates for the following
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
#This is a table of 20519 records that I can use for the PFAM run later on, but also allows to deal with DefenseFinder
GFFSelectedRegionsDF<-subset(GFFMergeWithRegionsDF,
                          GFFMergeWithRegionsDF$ROI =="Y")


###Saving this to run the PFAM annotation later on
GFFIDsToSave<-GFFSelectedRegionsDF[,c("GenomeID","ID")]
write.table(GFFIDsToSave, file = "PFAM_GeneIDs_15000.tsv",
            sep="\t",quote = F, row.names = F, col.names = F)

##########################################
##Reading PFAM hmmscan output
PFAMpath<-"./PFAM_search_20250515/"
PFAMselectedFiles<-list.files(pattern="\\.csv$",
                                recursive = T,
                                path = PFAMpath)

setwd(paste0(mainpath,"/",PFAMpath))
PFAMdata<-readr::read_csv(PFAMselectedFiles, id="file_name")
PFAMdata$GenomeID<-str_replace(PFAMdata$file_name,"_pfamscan.csv","")
setwd(mainpath)

PFAMdata$PFAMID<-str_split_i(PFAMdata$hmm_acc,"\\.",1)
###PFAM phage associated domains
PFAMPhage<-read.csv("PFAM_list_of_Phage_domains_fromRoman.tsv",
                    sep="\t", header=F)
PFAMPhage$IsPhage<-rep("Ph", length(PFAMPhage$V1))
##
PFAMdataWPhage<-merge(PFAMdata,PFAMPhage[,c(1,4)],
                      by.x="PFAMID", by.y = "V1", all.x =T)

PFAMdataMDom<-PFAMdataWPhage %>% group_by(GenomeID,seq_id) %>%
  fill(IsPhage) %>%
  summarise(name=str_c(unique(hmm_name), collapse="/"),
            phage = unique(IsPhage))
names(PFAMdataMDom)[2]<-"ID"


######

##Merge PFAM with GFF
GFFwithPFAM<-merge(GFFSelectedRegionsDF[,c(1,2,6,7,9,12,13)],
                   PFAMdataMDom,
                   by =  c("GenomeID","ID"), all.x = T)

############################################
##Read DefenseFinder genes files
DefenseFinderpath<-"./DefenseFinder_output_20250514"
DefenseFinderselectedFiles<-list.files(pattern="\\_defense_finder_genes.tsv$",
                                recursive = T,
                                path = DefenseFinderpath)
setwd(paste0(mainpath,"/",DefenseFinderpath))

DefenseFinderdata<-readr::read_tsv(DefenseFinderselectedFiles, id="file_name")
DefenseFinderdata$GenomeID<-str_split_i(DefenseFinderdata$file_name,"/",1)

setwd(mainpath)

###
##Merge DefenseFinder results with GFF selected regions
DFwithCoord<-merge(GFFwithPFAM,
                   DefenseFinderdata, by.x = c("GenomeID","ID"),
                   by.y = c("GenomeID","hit_id"), all.x =T)

DFwithCoordSh<-DFwithCoord[,c(1:8,9,12,15,32:34)]
##################################
#####Create the final annotation file
##Merge With PADLOC
AllFunctionalAnno<-merge(DFwithCoordSh, 
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
######


AllFuncWithMol<-merge(AllFunctionalAnnoNoDupl[,c(1,2,9,17,18,23:27)], UniqDru3Borders[,c(3,5,21)],
      by=c("GenomeID","seqid"))
GenesToPlot<-AllFuncWithMol[,c(11,1:10)]
names(GenesToPlot)[1]<-"molecule"
GenesToPlot$fill<-ifelse(is.na(GenesToPlot$system), GenesToPlot$phage,
                         ifelse(startsWith(GenesToPlot$system,"PD"),
                                GenesToPlot$phage,
                                GenesToPlot$system))

GenesColors<-hue_pal()(78)#rainbow(69)
GenesColors<-append(GenesColors, "#525252", after=2)
allSysNames<-unique(GenesToPlot$fill[!is.na(GenesToPlot$fill)])
names(GenesColors)<-allSysNames[order(allSysNames)]



##################################################################
#####Getting all other types of INFO for plotting
##################################################################

############Taxonomic distribution
DruE3padlocDf<-read.csv("../20250423_Druantia3_whole_systems/DruE3_padloc20_refseq.nopseudo.cutoff07.withDruH.csv",
                        header=F)
#adjust genome IDs
DruE3padlocDf$GenomeID<-sapply(strsplit(DruE3padlocDf$V1, ".csv:", fixed=TRUE), 
                               head, 1)
##get info on clusters
DruE3Clusters<-read.csv("../20250423_Druantia3_whole_systems/DruE3_mmseqs98_cluster.withref.tsv",
                        header=F, sep="\t")


######Read assembly summary
AssemblySummary<-read.csv("../20250423_Druantia3_whole_systems/Dru3_assembly_summary.tsv", 
                          header = F, sep="\t",quote = "", 
                          row.names = NULL, 
                          stringsAsFactors = FALSE)

##get representative genomes
Dru3ReprGenomes<-read.csv("../20250424_DruantiIII_neighborhoods_genomad/Dru3_representative_genomes.txt", header=F)


####Adding cluster info to padloc df
DruE3padlocDfWithCl<-merge(DruE3padlocDf,DruE3Clusters, by.x ="V4", by.y="V2", all.x =T)

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

##########
DruE3MetadataPADLOCLong<-merge(DruE3padlocDfWithCl[,c(20,21,1,4,7,8)],
                               AssemblySummary, by.y="V1", by.x="GenomeID")

DruE3MetadataPADLOCbyWp<-DruE3MetadataPADLOCLong %>%
  group_by(V1.y,TopOnly)%>%
  count()
names(DruE3MetadataPADLOCbyWp)<-c("molecule","Genus","Count")

#adding a small value so it will be visually seen which genus is there
DruE3MetadataPADLOCbyWp$log10Count<-0.1+log10(DruE3MetadataPADLOCbyWp$Count)
DruE3MetadataPADLOCbyWp$Genus<-factor(DruE3MetadataPADLOCbyWp$Genus, 
                                      levels=c(TopGenus, 'Else'))

#################################################################
########################################
###Plotting Genomad data on the DruE phylogenetic tree
DruE3padlocDfWithClRepr<-subset(DruE3padlocDfWithCl,
                                DruE3padlocDfWithCl$GenomeID %in% Dru3ReprGenomes$V1)

DruE3padlocDfWithClRepr$SystemID<-sapply(strsplit(DruE3padlocDfWithClRepr$V1.x, ".csv:", fixed=TRUE), 
                                         tail, 1)
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
colnames(MGEDataForPlot)[2:3]<-c(1,2)
MGEDataForPlot[is.na(MGEDataForPlot)]<-0
MGEDataForPlotLong<-gather(MGEDataForPlot,
                           key = "Location", value = "Prediction", 2:3)
names(MGEDataForPlotLong)[1]<-'molecule'
MGEDataForPlotLong$Location<-as.integer(MGEDataForPlotLong$Location)

TipColorDF<-data.frame(label=DruE3TreeMidRoot$tip.label)
TipColorDF$color<-ifelse(TipColorDF$label == "WP_020219138.1", "ECOR19", "")

TipColorScheme<-c("#cb181d", "#636363")
names(TipColorScheme)<-c("ECOR19","")

#cluster colors
myclustercolors<-c("#1f78b4","#a6cee3","#8c510a","#00441b",
                   "#b2df8a","#cab2d6","#4d4d4d","#dfc27d",
                   "#fb9a99","#fdbf6f","#6a3d9a","#33a02c",
                   "#e31a1c","#ff7f00","#ffff99")
names(myclustercolors)<-unique(DruE3clades$Cluster)

#################################################################
#plot long tree
DruE3BasicWBt<-ggtree(DruE3TreeMidRoot,
                           size=1,
                           color="#636363") +
  geom_nodelab(size=3) +
  geom_tiplab()+
  scale_x_continuous(limits = c(0,3.6))+
  geom_hilight(data = DruE3NodeOfInterest,
               mapping = aes(node=ClMRCA,fill = as.factor(Cluster)),
               alpha=0.2,
               extend=.05) +
  scale_fill_manual(values=myclustercolors,name="Cluster") +
  new_scale_fill()

DruE3BasicTreePlot<- DruE3BasicWBt %<+%  TipColorDF +
  geom_tippoint(aes(colour=color),size=3)+
  scale_colour_manual(values=TipColorScheme, guide="none")+
  theme_tree2() 

##get leaf order
leaf_order<-DruE3BasicTreePlot$data %>%
  filter(isTip) %>% arrange (y)
##
###
#doing species plot
DruE3MetadataPADLOCbyWp$molecule<-factor(DruE3MetadataPADLOCbyWp$molecule,
                                         levels=leaf_order$label)

SpeciesPlot<-ggplot(data=DruE3MetadataPADLOCbyWp,
                    aes(x= log10Count,
                    y = molecule,
                    fill = Genus))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=TreeColors)+
  theme_tree2()+
  theme(legend.position = 'none')

###
#doing HGT plot
MGEDataForPlotLong$molecule<-factor(MGEDataForPlotLong$molecule,
                                         levels=leaf_order$label)
HGTPLot<-ggplot(data = MGEDataForPlotLong, 
                aes(y=molecule,
                    x=Location,
                    fill = as.character(Prediction),
                    color= as.character(Prediction)))+
  geom_tile()+
  scale_color_manual(values=c("#ffffff","#bebada","#8dd3c7"), guide="none")+
  scale_fill_manual(values=c("#ffffff","#bebada","#8dd3c7"), guide="none")+
  theme_tree2()+
  theme(legend.position = 'none')

###
#creating genes plot
GenesToPlot$molecule<-factor(GenesToPlot$molecule, levels = leaf_order$label)
#trying to create compatible coordinates between samples
#the idea is to align everything on DruE
DruE3CoordFP<-subset(GenesToPlot,
                     GenesToPlot$protein.name == "DruE3")
#DruE3CoordFP$middle<-DruE3CoordFP$Gstart + (DruE3CoordFP$Gend-DruE3CoordFP$Gstart)/2
#I want to flip the ones that have other orientation
DruE3CoordFP<-DruE3CoordFP%>%group_by(molecule) %>%
  mutate(middle = Gstart +(Gend-Gstart)/2)
GenesToPlotNC<-merge(GenesToPlot, DruE3CoordFP[,c(1,11,13)], by= "molecule")

GenesToPlotNC$pstart<-ifelse(GenesToPlotNC$ggstrand.y,
                             GenesToPlotNC$Gstart - GenesToPlotNC$middle,
                             GenesToPlotNC$middle - GenesToPlotNC$Gend)
GenesToPlotNC$pend<-ifelse(GenesToPlotNC$ggstrand.y,
                          GenesToPlotNC$Gend - GenesToPlotNC$middle,
                          GenesToPlotNC$middle - GenesToPlotNC$Gstart)
GenesToPlotNC$pstrand<-ifelse(GenesToPlotNC$ggstrand.y,
                              GenesToPlotNC$ggstrand.x,
                              !GenesToPlotNC$ggstrand.x)

###finally ploting
GenesPlot<-ggplot(data = GenesToPlotNC,
       aes(y =  molecule,
           xmin = pstart,
           xmax = pend))+ 
  geom_gene_arrow(aes(fill = fill,
                      forward=pstrand))+
  geom_gene_label(aes(label=Agene))+
  scale_fill_manual(values = GenesColors, na.value="white")+
  theme_tree2()+
  theme(legend.position = "none")


####################
##Arranging All plots together
p<- ggarrange(DruE3BasicTreePlot,
              SpeciesPlot,
              HGTPLot,
          GenesPlot, nrow =1,
          widths = c(0.5,.3,.05,1),
          legend = "none",
          align='hv')

######################
ggsave("DruE_long_with_neighborhoods_and_HGT_20000.pdf",
       plot=p,
       limitsize = FALSE,
       width =190,
       height=170,
       units="cm",
       dpi=1000)
  









  
  







