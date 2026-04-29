library(tidyr)
library(dplyr)
library(Biostrings)
library(ggplot2)
library(purrr)
library(ggpubr)
library(xml2)
library(this.path)

mainpath<-paste0(dirname(this.path()),"/../../")
setwd(mainpath)

###read chromosome
# load fasta
fasta <- readDNAStringSet("./data/endopore/pUC19_Cm.fasta")
#get pUC19
seq <- fasta[[1]]  # DNAString
seq_name <- names(fasta)[1]
###############

####analyze DruE###############################################
DruEData<-read.csv("./data/endopore/ZB9BNH_4_drue_breakpoints.tab", header = F, sep="\t")

colstosep<-c("V2","V3","V4","V5")

for (col in colstosep) {
  DruEData <- separate(DruEData, !!col, into = paste0(col, c("_chrom", "_coords")), sep = ":")
}

#Filtering the cases that have the correct orientation
DruEdataCmOnEnds<-subset(DruEData,
                         DruEData$V2_chrom == "Cm" & DruEData$V5_chrom == "Cm" &
                           DruEData$V3_chrom == "pUC19" & DruEData$V4_chrom == "pUC19")

###getting most common pUC19 cut coordinates

pUCcolstosep<-c("V3_coords","V4_coords")

for (col in pUCcolstosep) {
  DruEdataCmOnEnds <- separate(DruEdataCmOnEnds, !!col, into = paste0(col, c("_start", "_end")), sep = "-")
}

###get most common ends
DruEMostCommonEnds<-DruEdataCmOnEnds %>%
  pivot_longer(
    cols = c(V3_coords_start, V4_coords_end),
    names_to = "source",
    values_to = "coord"
  ) %>%
  count(source, coord, sort = TRUE) %>%
  filter(n > 100)

DruEPositionsToExtract<-as.numeric(unique(DruEMostCommonEnds$coord))

###extracting sequences
window = 40

DruEneighborhoods <- sapply(DruEPositionsToExtract, function(pos) {
  start <- max(1, pos - window)
  end   <- min(length(seq), pos + window)
  as.character(subseq(seq, start, end))
})
#builddf
DruEseqdf<- data.frame(
  coord = DruEPositionsToExtract,
  neighborhood = DruEneighborhoods
)
#create fasta objects
DruEseqstosave <- DNAStringSet(DruEseqdf$neighborhood)
names(DruEseqstosave) <- paste0("coord_", DruEseqdf$coord)
#save neighborhoods
writeXStringSet(DruEseqstosave, filepath = "./data/endopore/ZB9BNH_4_drue_neighborhoods.fasta", format = "fasta")
#loading those to MEME Suite, but I am not able to find any obvious motifs around the region
#no strong specific motif is found
################
#Let's vizualize most common cuts
DruEEndsForPlot <- DruEdataCmOnEnds %>%
  pivot_longer(
    cols = c(V3_coords_start, V4_coords_end),
    values_to = "coord"
  ) %>%
  count(coord, sort = TRUE)

DruEtop_coords <- DruEMostCommonEnds %>%
  arrange(desc(n)) %>%
  slice_head(n = 15) %>%
  mutate(coord_num = as.numeric(coord))

DruECutPlot1<-ggplot(data = DruEEndsForPlot) +
  geom_col(aes(x = as.numeric(coord), y = n), fill = "#7570b3",
           width = 3) +
  annotate("rect", xmin = 1, xmax = length(seq), ymin = 50000, ymax = 60000,
           fill = "#1b9e77", alpha = 0.9)+
  annotate("text", x = (1 + length(seq)) / 2,  
           y = 42000,                          
           label = "pUC19",
           color = "black", size = 5, fontface = "bold")+
  # geom_rect(xmin=1,xmax=length(seq),ymin=30000,ymax=40000,
  #           fill = "#1b9e77") +
  geom_point(data = DruEtop_coords,
               aes(x = coord_num, 
                   y = 68000),
               fill = "#e6ab02",
               shape=25, size =6) +
  scale_fill_manual(values= c("#7570b3","#e7298a","#66a61e","#e6ab02"),
                    guide ="none")+
  scale_y_log10(
    expand=c(0,0),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_continuous(expand = c(0, 0),
    breaks = seq(0, length(seq), by = 100)
  ) +
  labs(x = "Coordinate", y = "log10(read count)") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black", size =14)
  ) +
  ggtitle("DruE")
DruECutPlot1
ggsave(file = "./figures/endopore/DruE_most_common_breakpoints_20260416_15s30s.pdf",
       plot = DruECutPlot1,
       width =40, height =18, units ="cm",dpi=300)

#####################
####Draw distances between cuts
DruEDistForPlot<-subset(DruEdataCmOnEnds, as.numeric(DruEdataCmOnEnds$V3_coords_end) > 2600 &
                          as.numeric(DruEdataCmOnEnds$V4_coords_start) < 100)
DruEDistForPlot$distance<-as.numeric(DruEDistForPlot$V3_coords_start) - as.numeric(DruEDistForPlot$V4_coords_end)
DruEDistForPlot$fill<-ifelse(as.numeric(DruEDistForPlot$distance) < 1, 
                             "potential primary",
                             "secondary")

DruEDistAll<-ggplot(data = DruEDistForPlot) +
  geom_histogram(aes(x=distance,
                     fill = fill), binwidth = 10)+
  theme_classic()+
  ylab("log10(read count)")+
  scale_x_continuous(limits = c(-100,2000))+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  scale_fill_manual(values = c("#7570b3",
                               "#1b9e77"), name = "")+
  theme(
    axis.text = element_text(color = "black", size =14),
    legend.text = element_text(color = "black", size =14),
    axis.title = element_text(color = "black", size =14)
  )
DruEDistAll
ggsave(file = "./figures/endopore/DruE_cleavage_distances_all.pdf",
       plot = DruEDistAll,
       width =30, height =20, units ="cm",dpi=300)
##zoom in 
##but it seems that mostly data is on secondary cuts
DruEDistZoom<-ggplot(data = DruEDistForPlot) +
  geom_histogram(aes(x=distance,
                     fill = fill
                     #y = after_stat(count*100 / sum(count))
  ), binwidth = 1)+
  theme_classic()+
  scale_x_continuous(limits = c(-50,100),
                     breaks = c(seq(-50,100, by =10)))+
  scale_y_continuous(limits = c(0,500),name = "read count")+
  scale_fill_manual(values = c("#7570b3",
                               "#1b9e77"), guide = "none")+
  theme(
    axis.text = element_text(color = "black", size =14),
    axis.title = element_text(color = "black", size =14)
  )
DruEDistZoom
# ggsave(file = "DruE_cleavage_distances_zoom.pdf",
#        plot = DruEDistZoom,
#        width =30, height =20, units ="cm",dpi=300)
#but this plot is not very useful because it counts the same cleavage location multiple time
#####################
####Draw most common cleavage scenarios
DruECleavageDistCases<-DruEDistForPlot %>% group_by(V3_coords_start,V4_coords_end) %>% 
  summarize (count = n()) %>%
  arrange(-count)
DruECleavageDistCasesTop20<-DruECleavageDistCases[c(1:20),]
DruECleavageDistCasesTop20<-DruECleavageDistCasesTop20 %>% arrange(as.numeric(V3_coords_start))
DruECleavageDistCasesTop20$rank <- seq(20,1)

DruEtmptop<-DruECleavageDistCasesTop20[,c(1,3,4)]
colnames(DruEtmptop)[1]<-"start"
DruEtmptop$end<-rep("2686",length(DruEtmptop$rank))
DruEtmptop$ycoord<-DruEtmptop$rank+0.5
DruEtmpbottom<-DruECleavageDistCasesTop20[,c(2:4)]
colnames(DruEtmpbottom)[1]<-"end"
DruEtmpbottom$start<-rep("1",length(DruEtmpbottom$rank))
DruEtmpbottom$ycoord<-DruEtmpbottom$rank
DruECleavageTopForPlot<-rbind(DruEtmptop,DruEtmpbottom)

DruEMostCommonCases<-ggplot(DruECleavageTopForPlot) + geom_rect(aes(xmin = as.numeric(start),
                                               xmax = as.numeric(end), 
                                               ymin = ycoord - sqrt(count/(5*max(count))), 
                                               ymax =ycoord+sqrt(count/(5*max(count))),
                                               fill= as.factor(rank)))+
  scale_x_continuous(expand =c(0,0))+
  # scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
  #                            "#fb9a99","#e31a1c","#fdbf6f","#ff7f00",
  #                            "#cab2d6","#6a3d9a"),
  #                   guide = "none")+
  guides(fill="none")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
#DruEMostCommonCases

DruECombinedPlot<-ggarrange(DruECutPlot1,
          DruEMostCommonCases, 
          ncol = 1,
          heights = c(1,0.3),
          align = "v")
DruECombinedPlot

ggsave(file = "./figures/endopore/DruE_most_common_breakpoints_20260416_15s30s_comb.pdf",
       plot = DruECombinedPlot,
       width =40, height =22, units ="cm",dpi=300)


#####Save filtered tables
SaveFilteredData<-function(DF,name){
  DataToSave <- DF %>%
    pmap_chr(~ {
      vals <- c(...)
      
      id <- vals[1]
      
      pairs <- c(
        # 2:3  → Cm 435-867
        paste0(vals[2], ":", vals[3]),
        
        # 4:5-6 → pUC19 2280 2686
        paste0(vals[4], ":", vals[5], "-", vals[6]),
        
        # 7:8-9
        paste0(vals[7], ":", vals[8], "-", vals[9]),
        
        # 10:11
        paste0(vals[10], ":", vals[11])
      )
      
      paste(c(id, pairs), collapse = "\t")
    })
  writeLines(DataToSave, paste0("./data/endopore/",name,"_breakpoints_filtered.tab"))
}

SaveFilteredData(DruEdataCmOnEnds,"DruE")