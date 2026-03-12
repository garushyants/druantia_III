library(dplyr)
library(readr)
library(ggplot2)
library(this.path)

###################
mainpath<-paste0(dirname(this.path()),"/../")
setwd(mainpath)

FigDir<-"figures/habitat_info"

if (!dir.exists(FigDir)){
  dir.create(FigDir,recursive = TRUE)
} else {
  print("Directory already exists!")
}

KimGenomeMetadata<-read.csv("data/habitat_data/spire_v1_genome_metadata.tsv", sep="\t")
KimMicrontology<-read.csv("data/habitat_data/spire_v1_microntology_mod.tsv", sep="\t", header =F)
DruantiaInfo<-read.csv("data/Supplementary_table_1_dataset_info.tsv", sep="\t")

###collapsing by the longest annotation
KimMicrontologyCollapsed<-KimMicrontology %>%
  group_by(V1) %>%
  mutate(colon_count = str_count(V3, ":")) %>%
  slice_max(colon_count, n = 1, with_ties = FALSE) %>%
  select(-colon_count)

KimGenomeMetadataWSpecies<- KimGenomeMetadata %>% 
  left_join(
    KimMicrontologyCollapsed %>%
      distinct(V1, .keep_all = TRUE),
    by = c("derived_from_sample" = "V1")
  )

#Select only the ones that have info up to species
HabitatDataShort<-subset(KimGenomeMetadataWSpecies[,c(27,28,32)],
                         KimGenomeMetadataWSpecies$species !="")

DruantiaDruEInfo<-subset(DruantiaInfo, 
                         DruantiaInfo$Protein == "DruE3")

###trying to merge with habitats
##first at species level
DruantiaDruEInfoWHabitats<-DruantiaDruEInfo %>% 
  left_join(
    HabitatDataShort[,c(2,3)] %>%
              distinct(species, .keep_all = TRUE),
            by = c("Species.name" = "species")
  )
#then on genus level
DruantiaDruEInfoWHabitatsGenus<-DruantiaDruEInfoWHabitats %>% 
  left_join(
    HabitatDataShort[,c(1,3)] %>%
      distinct(genus, .keep_all = TRUE),
    by = c("Genus.name" = "genus")
  )
#merging those with priority on species level
DruantiaDruEInfoWHabitatsGenus$Habitat<-ifelse(!is.na(DruantiaDruEInfoWHabitatsGenus$V3.x),
                                               DruantiaDruEInfoWHabitatsGenus$V3.x,
                                               DruantiaDruEInfoWHabitatsGenus$V3.y)

###Now let's group by tree representative 
DruEOnlyRepresentatives_habitat <- DruantiaDruEInfoWHabitatsGenus %>%
  group_by(TreeRepresentative) %>%
  summarise(
    Cluster = first(Cluster),
    Habitat_combined = first(unique(na.omit(Habitat))),
    Habitat_species_only=first(sort(na.omit(V3.x)))
  )

#Let's recalculate in percentage
DruEOnlyRepresentatives_habitat_percent <- DruEOnlyRepresentatives_habitat %>%
  group_by(Cluster)%>%
  count(Habitat_combined) %>%
  mutate(perc = n / sum(n) * 100)
####
clustercolors<-c("#a6cee3","#1f78b4","#b2df8a",
                 "#33a02c","#fb9a99","#e31a1c",
                 "#8c510a","#ff7f00","#cab2d6",
                 "#ffff99","#6a3d9a","#fdbf6f",
                 "#00441b","#dfc27d","#4d4d4d")

################
#Plot data for tree representatives only
#remove NAs
DruERepr_habitat_nona<-subset(DruEOnlyRepresentatives_habitat_percent,
       !is.na(DruEOnlyRepresentatives_habitat_percent$Habitat_combined))
###
ReprHabitatPlot<-ggplot(data = DruERepr_habitat_nona) +
  geom_col(aes(x = Habitat_combined,
               y = perc,
               fill = as.factor(Cluster)),
           color="#bdbdbd",
           linewidth=.2)+
  # facet_wrap(~Cluster, ncol =15,
  #            scales = "free")+
  facet_grid(. ~ Cluster, 
             scales = "free_x", 
             space = "free_x")+
  scale_fill_manual(values = clustercolors, name = "Cluster" ) +
  scale_y_continuous(breaks = seq(0,100,by=10))+
  xlab("habitat")+
  ylab("%")+
  theme_classic()+
  theme(axis.text.x = element_text(angle =90,hjust=1,vjust=.5, size=11),
        axis.text.y =element_text(size=12),
        legend.text = element_text(size=12),
        legend.title =element_text(size=12),
        strip.background = element_blank(),
        strip.text = element_blank() )
ReprHabitatPlot

ggsave("Druantia_III_habitats_spire_representatives.pdf",
       path = FigDir,
       plot = ReprHabitatPlot,
       width = 20, height=8,
       dpi=300)

# ####################
# #Plot data for all genomes
# #I am not using this one because it is obviously very unbalanced
# DruEAllHabitats_percent<-DruantiaDruEInfoWHabitatsGenus %>%
#   group_by(Cluster)%>%
#   count(Habitat) %>%
#   mutate(perc = n / sum(n) * 100)
# #remove NAs
# DruEAllHabitats_percent_nona<-subset(DruEAllHabitats_percent,
#                               !is.na(DruEAllHabitats_percent$Habitat))
# ##
# ggplot(data = DruEAllHabitats_percent_nona) +
#   geom_col(aes(x = Habitat,
#                y = perc,
#                fill = as.factor(Cluster)),
#            color="#bdbdbd",
#            linewidth=.2,)+
#   # facet_wrap(~Cluster, ncol =15,
#   #            scales = "free")+
#   facet_grid(. ~ Cluster,
#              scales = "free_x",
#              space = "free_x")+
#   scale_fill_manual(values = clustercolors, name = "Cluster" ) +
#   xlab("habitat")+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle =90,hjust=1,vjust=.5))
