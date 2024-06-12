#The code for creating the 16S cyanobacterial tree and grep for 18S sequence matches are also in here

####Packages and other universal data for analyses####
library(biomformat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggtree)
library(ggnewscale)
library(forcats)
library(stringr)
library(mgsub)

#colour
col <- species_colours <- c((brewer.pal(8, "Set2")[c(1:8)]), (brewer.pal(12, "Paired")[c(1:12)]), "black", "white", (brewer.pal(12, "Set3")[c(1:12)]))

####16S dataframe####
####16S Feature-table####
OTU_reads <- read_biom("../tax_data/16S_feature-table.biom")
OTU_table_16S <- as.data.frame(as.matrix(biom_data(OTU_reads)))
OTU_table_16S$Feature.ID <- rownames(OTU_table_16S)
OTU_table_16S <- OTU_table_16S %>% select("S240", "S241", "S421", "S422", "S430", "S439", "S440", "S441", "Feature.ID", "E4-neg")
#colnames(OTU_table_16S)
OTU_table_16S <- pivot_longer(OTU_table_16S, cols = 1:8, names_to = "Fish.ID", values_to = "Reads")
OTU_table_16S <- filter(OTU_table_16S, Reads > 0)
OTU_table_16S$Reads_nrm <- OTU_table_16S$Reads - OTU_table_16S$`E4-neg`

####16S Taxonomy####
#Cyanoseq for Cyanobacteria
CyanoSeq <- read.delim("../tax_data/16S_all_taxonomy_CyanoSeq.tsv")
CyanoSeq <- separate(CyanoSeq, col = "Taxon", into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = F) #will say discarded data, but it is just NAs/blanks - some strings have ';' at the end and so extra column is created
CyanoSeq <- CyanoSeq %>% mutate(kingdom = rep("k_spk_Bacteria"))
CyanoSeq <- CyanoSeq[,c(1:3, 11, 4:10)]

#PR2 for Chloroplast sequences
PR2 <- read.delim("../tax_data/16S_all_taxonomy_PR2.tsv")

#SILVA for everything else
SILVA <- read.delim("../tax_data/16S_all_taxonomy_SILVA.tsv")
SILVA <- separate(SILVA, col = "Taxon", into = c("domain", "superkingdom", "kingdom", "subkingdom", "major_clade", "superphylum", "phylum", "subphylum", "infraphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", "superfamily", "family", "subfamily", "genus", "species"), sep = ";", remove = F)

##############################Cyanoseq check##############################
#There were differences in the number of ASVs for S riv assigned to CyanoSeq (higher) than SILVA. The reason in SILVA they are assigned as Chloroplast, but in CyanoSeq they do not go below phylum level.
#CyanoSeq_only <- anti_join(left_join(OTU_table_16S, CyanoSeq) %>% filter(grepl("Cyanobacteriota", Taxon)==T) %>% filter(grepl("Chloroplast", Taxon)==F)  %>% select(Feature.ID) %>% unique(), left_join(OTU_table_16S, SILVA) %>% filter(grepl("Cyanobacteria", Taxon)==T) %>% filter(grepl("Chloroplast", Taxon)==F) %>% select(Feature.ID) %>% unique())
#SILVA %>% filter(Feature.ID %in% CyanoSeq_only$Feature.ID) %>% View
#CyanoSeq %>% filter(Feature.ID %in% CyanoSeq_only$Feature.ID) %>% View
###########################################################################

#to extract sequences from Pr2
Chloroplast_Feat <- SILVA %>% filter(grepl("Chloroplast", Taxon)) %>% select(Feature.ID) %>% unique()

PR2_Chloro <- PR2 %>% filter(Feature.ID %in% Chloroplast_Feat$Feature.ID)
PR2_Chloro <- separate(PR2_Chloro, col = "Taxon", into = c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = F)

#Do this once df_16S has been created
#Chloroplast spp from SILVA for S riv
SILVA %>% filter(Feature.ID %in% OTU_table_16S$Feature.ID) %>% filter(Feature.ID %in% df_16S$Feature.ID) %>% filter(grepl("Cyanobac",Taxon)==T & order == "Chloroplast") %>% select(species) %>% unique()

x <- SILVA %>% filter(Feature.ID %in% OTU_table_16S$Feature.ID) %>% filter(Feature.ID %in% df_16S$Feature.ID) %>% filter(grepl("Cyanobac",Taxon)==T & order == "Chloroplast") %>% select(Feature.ID, species, Confidence) %>% rename(SILVA_species = species, con_SILVA=Confidence)

SILVA %>% filter(Feature.ID %in% OTU_table_16S$Feature.ID) %>% filter(Feature.ID %in% df_16S$Feature.ID) %>% filter(grepl("Cyanobac",Taxon)==T & order == "Chloroplast") %>% select(Feature.ID) %>% left_join(., PR2_Chloro) %>% select(genus) %>% unique()

y <- SILVA %>% filter(Feature.ID %in% OTU_table_16S$Feature.ID) %>% filter(Feature.ID %in% df_16S$Feature.ID) %>% filter(grepl("Cyanobac",Taxon)==T & order == "Chloroplast") %>% select(Feature.ID) %>% left_join(., PR2_Chloro) %>% select(Feature.ID, family, genus, species, Confidence) %>% rename(PR2_species = species, PR2_con = Confidence)

full_join(x,y) %>% View

#to extract sequences from Cyanoseq#########################
Cyanobac_feat <- SILVA %>% filter(grepl("Cyanobacteria", Taxon)==T) %>% filter(grepl("Chloroplast", Taxon)==F) %>% .[,c("Feature.ID")] %>% unique()

CyanoSeq_Cyano <- CyanoSeq %>% filter(Feature.ID %in% Cyanobac_feat)

SILVA_rem <- SILVA %>% filter(! Feature.ID %in% CyanoSeq_Cyano$Feature.ID) %>% filter(! Feature.ID %in% PR2_Chloro$Feature.ID)

#Do this once df_16S has been created
#Cyanobacteria spp from SILVA for S riv
SILVA %>% filter(Feature.ID %in% OTU_table_16S$Feature.ID) %>% filter(Feature.ID %in% df_16S$Feature.ID) %>% filter(grepl("Cyanobac",Taxon)==T & order != "Chloroplast") %>% select(species) %>% unique()

x <- SILVA %>% filter(Feature.ID %in% OTU_table_16S$Feature.ID) %>% filter(Feature.ID %in% df_16S$Feature.ID) %>% filter(grepl("Cyanobac",Taxon)==T & order != "Chloroplast") %>% select(Feature.ID, species, Confidence) %>% rename(SILVA_species = species, con_SILVA=Confidence)

SILVA %>% filter(Feature.ID %in% OTU_table_16S$Feature.ID) %>% filter(Feature.ID %in% df_16S$Feature.ID) %>% filter(grepl("Cyanobac",Taxon)==T & order != "Chloroplast") %>% select(Feature.ID) %>% left_join(., CyanoSeq_Cyano) %>% select(genus) %>% unique()

y <- SILVA %>% filter(Feature.ID %in% OTU_table_16S$Feature.ID) %>% filter(Feature.ID %in% df_16S$Feature.ID) %>% filter(grepl("Cyanobac",Taxon)==T & order != "Chloroplast") %>% select(Feature.ID) %>% left_join(., CyanoSeq_Cyano) %>% select(Feature.ID, family, genus, species, Confidence) %>% rename(Cyano_species = species, Cyano_con = Confidence)

full_join(x,y) %>% View

#Make df_16S
df_16S_tax <- rbind(PR2_Chloro, CyanoSeq_Cyano, SILVA_rem[,c("Feature.ID", "Taxon", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "Confidence")])

df_16S <- left_join(OTU_table_16S, df_16S_tax)

#change Fish IDs to more reader friendly
fish_name_change <- data.frame(Fish.ID = c("S240", "S241", "S421", "S422", "S430", "S439", "S441", "S440"),
                               Fish.ID2 = c("SP1", "SP2", "MR1", "MR2", "ND1", "CR1", "CR3", "CR2")) 

df_16S <- df_16S %>% mutate(Fish.ID = mgsub(Fish.ID, fish_name_change$Fish.ID, fish_name_change$Fish.ID2))
# A tibble: 1 x 2
# ASVs  Reads
# <int>  <dbl>
#   1  2089 116674

df_16S <- filter(df_16S, Reads_nrm > 5)
df_16S <- df_16S %>% filter(Confidence >= 0.97)
df_16S %>% summarise(ASVs = length(unique(Feature.ID)), Reads = sum(Reads_nrm))
# A tibble: 1 x 2
#ASVs  Reads
#<int>  <dbl>
#  1  1686 115198

#Supplementary Table S5 - df_16S
#write.table(df_16S, file = "./16S_ASV_tax.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

df_16S %>% summarise(ASVs = length(unique(Feature.ID)), Reads = sum(Reads_nrm))

####16S tree####
Cyano_tree_seq <- df_16S %>% filter(grepl("Cyanobacteriota", Taxon)==T) %>% filter(grepl("Chloroplast", Taxon)==F) %>% .[,c("Feature.ID")] %>% unique()

#Get sequences from 16S dna-sequences.fasta output from Qiime2
#write.table(x = Cyano_tree_seq, file = "../Cyano_tree/Cyan_filt.txt", quote = F, row.names = F, col.names = F)
#copy and paste into nano bash file - doesn't seem to work if directly uploaded
#grep -A 1 -Fwf Cyano_filt.txt /path_to/dna-sequences.fasta > Cyano_seq.fasta
#sed -i '/--/d' Cyano_seq.fasta
#module load QIIME2/2021.2
#qiime tools import --input-path ./Cyano_seq.fasta --output-path ./Cyano_seq.qza --type 'FeatureData[Sequence]'
#qiime phylogeny align-to-tree-mafft-fasttree --i-sequences Cyano_seq.qza --o-alignment Cyano_aligned --o-masked-alignment Cyano_masked_aligned --o-tree Cyano_unrooted.tree --o-rooted-tree Cyano_rooted.tree --verbose
#unzip Cyano_unrooted.tree.qza -d rm.tree
#cp rm.tree/*/data/*.nwk ./Cyano_Sriv_tree.nwk
#rm -r rm.tree

tree <- read.tree("../Cyano_tree/Cyano_Sriv_tree.nwk")

tip.label.order <- tree[["tip.label"]] %>% as.data.frame() %>% rename(., tip.label=1)

tree_meta <- 
  df_16S[,c("Feature.ID", "Fish.ID", "Reads_nrm")] %>% 
  group_by(Fish.ID) %>% 
  mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>%
  filter(Feature.ID %in% Cyano_tree_seq$Feature.ID) %>%
  left_join(., CyanoSeq %>% select(Feature.ID, order, genus) %>% unique()) %>% 
  mutate(order = if_else(is.na(order)==T, "Unassigned Cyanophyceae", order)) %>%
    mutate(genus = if_else(is.na(genus)==T, "", genus)) %>%
  ungroup() %>% 
  group_by(Feature.ID) %>% 
  mutate(Fish.no = as.character(length(Fish.ID))) %>% 
  ungroup() %>%
  select(-Fish.ID, -Reads_nrm) %>%
  aggregate(RRA~., ., mean) %>%
  rename(., tip.label=Feature.ID) %>% 
  unique() %>% 
  left_join(tip.label.order, .) 

rownames(tree_meta) <- tree_meta$tip.label

tree_meta$node <-  rep(1:nrow(tree_meta))

circ <- ggtree(tree, layout = "circular", branch.length = "branch.length", size = 0.01)

#to genus names as tip label
circ$data$tip_genus <- c(tree_meta$genus, rep(NA, times=91))

#add unicellular or filamentous information
uni_fil <- readxl::read_xlsx("./Sriv_uni_fila.xlsx", sheet = 1,col_names = T)
uni_fil <- rename(uni_fil, genus=Genus)

test <- c(tree_meta$genus, rep(NA, times=91))
test <- as.data.frame(test) %>% rename(., genus = test)

uni_fil2 <- left_join(test, uni_fil[,c("genus", "Type")])

circ$data$uniORfil <- uni_fil2$Type

circ$data$uniORfilCOLOUR <- if_else(grepl(pattern = "Filamentous", x = circ$data$uniORfil, ignore.case = T), "Filamentous", if_else(grepl(pattern = "Unicellular", x = circ$data$uniORfil, ignore.case = T), "Unicellular", NA))

#had a problem with align = true and error said it was due to colour
#Googled: and suggested to update ggtree
#version used in Oct 2023 3.1.4.992
#library(devtools)
#install_github("YuLab-SMU/ggtree") 3.9.1
#library(ggtree)

circ2 <- circ+
  geom_tiplab(aes(label=tip_genus, color = uniORfilCOLOUR), align=TRUE, linesize=.5, size = 8*0.35, offset = 0.1)

p1 <- 
  gheatmap(circ2, tree_meta %>% select(Fish.no), offset=0, width=.1, colnames = F, color = "black") +
  scale_fill_brewer(palette = "Greys")+
  labs(fill="Samples (n)")+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8, colour = "black"),
        legend.key.size = unit(0.3, "cm"),
        legend.position = "left",
        legend.box = 'vertical',
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt"))

p2 <- p1+
    new_scale_fill()

p3 <- gheatmap(p2, tree_meta %>% select(RRA), offset=0.035, width=.1, colnames = F, color = "black")+
  scale_fill_gradientn(colours = c("dark blue", "light blue", "white", "orange", "red"))+
  labs(fill="RRA (%)")+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8, colour = "black"),
        legend.key.size = unit(0.3, "cm"),
        legend.position = "left",
        legend.box = 'vertical',
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt"))

p4 <- p3+
  new_scale_fill()+
  geom_hilight(data = tree_meta, mapping=aes(node=node, fill=order), type = "rect", align = "right", alpha = 0.8)+
  # scale_fill_manual(values = c("#D9D9D9", #light grey #Aegeococcocales
  #                              "#E6AB02", #light brown #Chroococcales
  #                              "#666666", #grey #Coleofasciculales
  #                              "#FFED6F", #light yellow #Desertifilales
  #                              "#66A61E",#dark green #Gastranaerophilales
  #                              "#FCCDE5", #blue #Leptolyngbyales
  #                              "#CCEBC5",#light green #Nodosilineales
  #                              "#BEBADA", #light orange #Nostocales
  #                              "#80B1D3", #light purple #Oculatellales
  #                              "#B3DE69", #purple #Oscillatoriales
  #                              "#A6761D", #brown #Prochlorothrichales
  #                              "#D95F02", #orange-red #Spirulinales
  #                              "#FB8072",#red #Synechococcales
  #                              "white", #yellow-brown #Unassigned cyano
  #                              "#1B9E77",#Vampirovibrionales
  #                              "#8DD3C7" #grey #Verrucomicrobia
  # ))+
scale_fill_manual(values = c("#A6CEE3", #light grey #Aegeococcocales
                             "#1F78B4", #light brown #Chroococcales
                             "#B2DF8A", #grey #Coleofasciculales
                             "#33A02C", #light yellow #Desertifilales
                             "#FB9A99",#dark green #Gastranaerophilales
                             "#E31A1C", #blue #Leptolyngbyales
                             "#FDBF6F",#light green #Nodosilineales
                             "#FF7F00", #light orange #Nostocales
                             "#CAB2D6", #light purple #Oculatellales
                             "#6A3D9A", #purple #Oscillatoriales
                             "#FFFF99", #brown #Prochlorothrichales
                             "#FFD92F", #orange-red #Spirulinales
                             "#DFC27D",#red #Synechococcales
                             "white", #yellow-brown #Unassigned cyano
                             "#543005",#Vampirovibrionales
                             "#8DD3C7" #grey #Verrucomicrobia
))+
  geom_point(data = circ$data %>% mutate(bootstrap = if_else(grepl("0.9.*", label), "Bootstrap > 0.9", NA)) %>% filter(is.na(bootstrap)==F), aes(color=bootstrap), size=1)+
  #give colours for bootstrap, and tiplabels
  #scale_color_manual(values = c(brewer.pal(9, "YlOrRd")[c(9)], brewer.pal(11, "PRGn")[c(2)], brewer.pal(11, "PiYG")[c(10)]))+
  scale_color_manual(values = c(brewer.pal(9, "YlOrRd")[c(9)], "black", "grey"))+
  labs(fill="Cyanobacteria order", color = "")+
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, "cm"),
        legend.position = "right",
        legend.box = 'vertical',
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt")
  )

p4

#84 mm/174 mm wide x 234 mm
#ggsave(filename = "CyanoSeq_tree_Sriv.pdf", plot = p4, path = "./", device = "pdf", width = 234, height = 174, units = "mm", dpi = 300)

######

#16S bar figures
df_16S <- df_16S %>% mutate(phylum = if_else(class=="Bacillariophyta", "Bacillariophyta", phylum))

df_16S %>% select(phylum, Feature.ID, Fish.ID) %>% unique() %>% group_by(phylum) %>% summarise(ASVs = length(unique(Feature.ID)), FOO = length(unique(Fish.ID))) %>% write.csv(x = ., file = "./Phyla_16S.txt", quote = F, row.names = F)


col_16S <- c("#66C2A5", #Actinobacteriota
             "#FDBF6F", #Bacillariophyta
             "#8DA0CB", #Bacteroidota
             "#FC8D62", #Campylobacterota
             "#A6D854", #Chlorophyta
             "#FFD92F", #Cyanobacteria
             "#E5C494", #Firmicutes
             "#33A02C", #Fusobacteriota
             "#A6CEE3", #Myxococcota
             "#1F78B4", #NB1-j
             "#964B00", #Ochrophyta
             "#B3B3B3", #Other
             "#FB9A99", #Planctomycetota
             "#E78AC3", #Proteobacteria
             "#E31A1C", #Rhodophyta
             "#FF7F00" #Verrucomicrobiota
             )

Phyla_16S <- 
  df_16S %>% 
  #filter(grepl("Cyano", phylum)==T) %>%
    mutate(phylum = if_else(grepl("Rhodophyceae", Taxon)==T, "Rhodophyta", if_else(grepl("Bacillariophyta", Taxon)==T, "Bacillariophyceae", if_else(grepl("Phaeophyceae", Taxon)==T, "Phaeophyceae", if_else(grepl("Ochrophyta", Taxon)==T, "Other Ochrophyta", if_else(grepl("Chlorophyta", Taxon)==T, "Chlorophyta", phylum)))))) %>%
  filter(is.na(phylum)==F) %>% group_by(Fish.ID) %>% mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>% aggregate(RRA~phylum+Fish.ID, ., sum) %>% mutate(cat = if_else(RRA < 3, "Other", phylum)) %>%
  mutate(Fish.ID = gsub("^S", "", Fish.ID)) %>%
  ggplot()+
  geom_bar(stat = "identity", aes(x = Fish.ID, y = RRA, fill = cat), colour = "black")+
  #facet_grid(~., space = "free", scales = "free")+
  scale_fill_manual(values = col_16S)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7, colour = "black"),
    axis.text.x = element_text(size = 7, colour = "black"), #, angle = 90, vjust = 0.5, hjust=1
    axis.title = element_text(size = 8, colour = "black"),
    legend.title = element_text(size = 7, colour = "black"),
    legend.text = element_text(size = 7, colour = "black"),
    legend.key.size = unit(0.3, "cm"),
    strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
    strip.background = element_rect(fill = "white", color = "black")
  )+
  labs(fill = "Phylum")+
  ylab("Relative read abundance (%)")+
  xlab("Fish ID")

#ggsave(filename = "Phyla_16S_otherOchr.jpeg", plot = Phyla_16S, path = "./", device = "jpeg", width = 174, height = 60, units = "mm", dpi = 300)

#number of fish with algae
df_16S %>% filter(is.na(phylum)==F) %>% group_by(Fish.ID) %>% mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>% aggregate(RRA~phylum+Fish.ID, ., sum) %>% mutate(cat = if_else(RRA < 1, "Other", phylum)) %>% filter(grepl("phyta", phylum)) %>% group_by(phylum) %>% summarise(n = length(Fish.ID))

#bubble plot for algae
BP_16S <- 
  df_16S %>% 
  group_by(Fish.ID) %>% mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>% ungroup() %>% 
  mutate(Fish.ID = gsub("^S", "", Fish.ID)) %>%
  filter(grepl("Rhodophyta|Bacillariophyta|Ochrophyta|Chlorophyta|Dinoflagellata", phylum)==T) %>%
    mutate(phylum = if_else(grepl("Rhodophyceae", Taxon)==T, "Rhodophyta", if_else(grepl("Bacillariophyta", Taxon)==T, "Bacillariophyceae", if_else(grepl("Phaeophyceae", Taxon)==T, "Phaeophyceae", if_else(grepl("Ochrophyta", Taxon)==T, "Other Ochrophyta", if_else(grepl("Chlorophyta", Taxon)==T, "Chlorophyta", phylum)))))) %>%
  #mutate(phylum = gsub("[a-z].*$", "", phylum)) %>%
  group_by(Fish.ID,genus) %>% mutate(decide = if_else(grepl("uncultured", species)==F & is.na(species)==F, "species", "genus")) %>% ungroup() %>%
  group_by(Fish.ID,decide, genus, species) %>% mutate(ASVs = length(decide)) %>% ungroup() %>%
  group_by(Fish.ID,species) %>% mutate(ASVs_species = length(species)) %>% ungroup() %>%
  mutate(group = if_else(grepl("uncultured", species)==F & is.na(species)==F, species, genus)) %>%
  aggregate(RRA~Fish.ID+group+ASVs+phylum, ., sum) %>%
  #mutate(group = gsub("^[a-z].*_..._", "~", group)) %>%
  filter(grepl("_X", group)==F) %>%
  mutate(group = gsub("_", " ", group)) %>%
  ggplot()+
  geom_point(shape = 21, aes(x=Fish.ID, y=fct_rev(group), fill = RRA, size = ASVs))+
  facet_grid(phylum~., space = "free", scales = "free")+
  scale_fill_gradientn(colours = c("dark blue", "light blue", "white", "orange", "red"))+
  theme_bw()+
  theme(
    panel.spacing = unit(0, "lines"),
    panel.grid.major.y = element_line(colour = "gray"),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 7, colour = "black"),
    axis.text.x = element_text(size = 7, colour = "black"), #, angle = 90, vjust = 0.5, hjust=1
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 8, colour = "black"),
    legend.title = element_text(size = 7, colour = "black"),
    legend.text = element_text(size = 7, colour = "black"),
    legend.key.size = unit(0.3, "cm"),
    strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
    strip.text.x = element_text(angle = 0, size = 7, colour = "black"),
    strip.background = element_rect(fill = "white", color = "black")
  )+
  ylab("Genus or species")+
  xlab("Fish ID")+
  labs(fill = "RRA (%)", size = "ASVs (n)")

#ggsave(filename = "16S_bubble_plot_otherOch.jpeg", plot = BP_16S, path = "./", device = "jpeg", width = 100, height = 60, units = "mm", dpi = 300)

df_16S %>% filter(phylum == "Cyanobacteriota") %>% select(genus, Fish.ID) %>% unique() %>% group_by(genus) %>% summarise(n = length(Fish.ID)) %>% filter(n >= 4)

df_16S %>% filter(phylum == "Cyanobacteriota") %>% select(genus, Fish.ID) %>% unique() %>% filter(is.na(genus)==F) %>% group_by(Fish.ID) %>% summarise(n = length(genus))

#presence of algae for different phyla
df_16S %>% 
  group_by(Fish.ID) %>% mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>% ungroup() %>% 
  mutate(Fish.ID = gsub("^S", "", Fish.ID)) %>%
  filter(grepl("Rhodophyta|Bacillariophyta|Ochrophyta|Chlorophyta|Dinoflagellata", phylum)==T) %>%
  mutate(phylum = if_else(grepl("Rhodophyceae", Taxon)==T, "Rhodophyta", if_else(grepl("Bacillariophyta", Taxon)==T, "Bacillariophyceae", if_else(grepl("Phaeophyceae", Taxon)==T, "Phaeophyceae", if_else(grepl("Ochrophyta", Taxon)==T, "Other Ochrophyta", if_else(grepl("Chlorophyta", Taxon)==T, "Chlorophyta", phylum)))))) %>%
  filter(phylum == "Chlorophyta") %>% #change for algae interested in here
  aggregate(RRA~Fish.ID, ., sum)

#Cyanobacteria genera
Cyano_df_16S <- df_16S %>% group_by(Fish.ID) %>% mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>% ungroup() %>% group_by(Fish.ID, family) %>% mutate(ASVs = length(Feature.ID)) %>% ungroup() %>% filter(grepl("Cyanobacteriota", Taxon)==T) %>% mutate(genus = if_else(genus == "" | is.na(genus)==T, "Unassigned", genus)) %>% filter(family != "") %>% aggregate(RRA~Fish.ID+genus+ASVs+family, ., sum)

  
Cyano_df_16S <-
Cyano_df_16S %>%
ggplot()+
  geom_point(shape = 21, aes(x=Fish.ID, y=fct_rev(genus), fill = RRA, size = ASVs))+
  scale_fill_gradientn(colours = c("dark blue", "light blue", "white", "orange", "red"))+
  facet_grid(family~Fish.ID, space = "free", scales = "free")+
  theme_bw()+
  theme(
    panel.grid.minor.y = element_line(colour = "grey"),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text.y = element_text(size = 7, colour = "black", face = "italic"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 8, colour = "black"),
    legend.title = element_text(size = 7, colour = "black"),
    legend.text = element_text(size = 7, colour = "black"),
    legend.key.size = unit(0.3, "cm"),
    strip.text.y = element_text(angle = 0, size = 7, colour = "black", face = "bold"),
    strip.text.x = element_text(angle = 0, size = 7, colour = "black", face = "bold"),
    strip.background = element_rect(fill = "white", color = "black")
  )+
  ylab("Genus")+
  labs(fill = "RRA (%)", size = "ASVs (n)")

#ggsave(filename = "Cyano_df_16S.jpeg", plot = Cyano_df_16S, path = "./", device = "jpeg", width = 174, height = 200, units = "mm", dpi = 300)

####18S dataframe####
####18S Feature-table####
OTU_FD <- read_biom("../tax_data/feature-table_FD.biom")
OTU_FD <- as.data.frame(as.matrix(biom_data(OTU_FD)))
OTU_FD$Feature.ID <- rownames(OTU_FD)
OTU_FD$Primer <- rep("FD-Bra")

OTU_Zhan <- read_biom("../tax_data/feature-table_Zhan.biom")
OTU_Zhan <- as.data.frame(as.matrix(biom_data(OTU_Zhan)))
OTU_Zhan$Feature.ID <- rownames(OTU_Zhan)
OTU_Zhan$Primer <- rep("Zhan")

OTU_Zim <- read_biom("../tax_data/feature-table_Zim.biom")
OTU_Zim <- as.data.frame(as.matrix(biom_data(OTU_Zim)))
OTU_Zim$Feature.ID <- rownames(OTU_Zim)
OTU_Zim$Primer <- rep("Zimm")

OTU_Stoeck <- read_biom("../tax_data/feature-table_Stoeck.biom")
OTU_Stoeck <- as.data.frame(as.matrix(biom_data(OTU_Stoeck)))
OTU_Stoeck$Feature.ID <- rownames(OTU_Stoeck)
OTU_Stoeck$Primer <- rep("Stoeck")

OTU_table_18S <- rbind(OTU_FD, OTU_Stoeck, OTU_Zhan, OTU_Zim)

OTU_table_18S <- pivot_longer(OTU_table_18S, cols = 1:8, names_to = "Fish.ID", values_to = "Reads")
OTU_table_18S <- filter(OTU_table_18S, Reads > 0)
OTU_table_18S %>% group_by(Primer) %>% summarise(ASVs = length(unique(Feature.ID)), Reads = sum(Reads))
OTU_table_18S$Reads_nrm <- OTU_table_18S$Reads - OTU_table_18S$Neg
OTU_table_18S %>% group_by(Primer) %>% summarise(ASVs = length(unique(Feature.ID)), Reads = sum(Reads_nrm))

####18S Taxonomy####
#SILVA 138.1
tax_FD <- read.delim("../tax_data/FD_taxonomy.tsv")
tax_FD$Primer <- rep("FD-Bra")
tax_Zhan <- read.delim("../tax_data/Zhan_taxonomy.tsv")
tax_Zhan$Primer <- rep("Zhan")
tax_Zim <- read.delim("../tax_data/Zim_taxonomy.tsv")
tax_Zim$Primer <- rep("Zimm")
tax_Stoeck <- read.delim("../tax_data/Stoek_taxonomy.tsv")
tax_Stoeck$Primer <- rep("Stoeck")

SILVA_tax_18S <- rbind(tax_FD, tax_Stoeck, tax_Zhan, tax_Zim)

####Feature-table joined to taxonomy####
df_18S <- left_join(OTU_table_18S, SILVA_tax_18S)
df_18S <- separate(df_18S, col = "Taxon", into = c("domain", "superkingdom", "kingdom", "subkingdom", "major_clade", "superphylum", "phylum", "subphylum", "infraphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", "superfamily", "family", "subfamily", "genus", "species"), sep = ";", remove = F)

#change Fish IDs to more reader friendly
fish_name_change <- data.frame(Fish.ID = c("240", "241", "421", "422", "430", "439", "441", "440"),
                               Fish.ID2 = c("SP1", "SP2", "MR1", "MR2", "ND1", "CR1", "CR3", "CR2")) 

df_18S <- df_18S %>% mutate(Fish.ID = mgsub(Fish.ID, fish_name_change$Fish.ID, fish_name_change$Fish.ID2))

df_18S <- df_18S %>% mutate(Primer=gsub("Zimm", "Zim", Primer))

#Supplementary Table S6 - df_18S
#write.table(df_18S, file = "./18S_ASV_tax.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

####18S stats####
colnames(df_18S)

df_18S %>% group_by(Primer) %>% summarise(ASVs = length(unique(Feature.ID)), Reads = sum(Reads_nrm))
# A tibble: 4 x 3
# Primer  ASVs  Reads
# <chr>  <int>  <dbl>
#   1 FD-Bra   426 113341
# 2 Stoeck   277  48564
# 3 Zhan     293  30547
# 4 Zimm     575  50706

#how many ASVs were not fish
df_18S %>% filter(grepl("Teleostei", Taxon)==F) %>% group_by(Primer) %>% summarise(ASVs = length(unique(Feature.ID)), Reads = sum(Reads_nrm))
# A tibble: 4 x 3
# Primer  ASVs Reads
# <chr>  <int> <dbl>
#   1 FD-Bra   418 13353 #~89
# 2 Stoeck   269  7330 #~85
# 3 Zhan     287  5753 #~80
# 4 Zimm     561 44160 #

#remove fish
df_18S <- df_18S %>% filter(grepl("Teleostei", Taxon)==F)

df_18S <- df_18S %>% filter(Reads_nrm > 5)
df_18S %>% group_by(Primer) %>% summarise(ASVs = length(unique(Feature.ID)), Reads = sum(Reads_nrm))
# A tibble: 4 x 3
# Primer  ASVs Reads
# <chr>  <int> <dbl>
#   1 FD-Bra   283 12863
# 2 Stoeck   166  6979
# 3 Zhan     169  5358
# 4 Zimm     438 43754

df_18S %>% filter(Confidence < 0.97)
#all ASVs were assigned

#18S phyla bar with RRA

#so cholorophyta is green, Rhodo red and Ochro brown
colours <- c("#66C2A5", #Annelida
             "#FC8D62", #Apicomplexa
             "#8DA0CB", #Athropoda
             "#FDBF6F", #Bacillariophyta
             "#A6D854", #Chlorophyta
             "#FFD92F", #Ciliophora
             "#E5C494", #Cnidaria
             "#33A02C", #Dinoflag
             "#A6CEE3", #Mollusca
             "#1F78B4", #Nematoda
             "#B3B3B3", #Other
             "#964B00", #Phaeo
             "#FB9A99", #Platyhelm
             "#E78AC3", #Porifera
             "#E31A1C", #Rhodophyta
             "#FF7F00", #Tunicata
             "#CAB2D6","#6A3D9A","#FFFF99","#B15928","black","white","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")

Phyla_18S <- 
df_18S %>% 
  mutate(phylum = if_else(grepl("Rhodophyceae", Taxon)==T, "Rhodophyta", if_else(grepl("Diatomea", Taxon)==T, "Bacillariophyceae", if_else(grepl("Phaeophyceae", Taxon)==T, "Phaeophyceae", if_else(grepl("Ochrophyta", Taxon)==T, "Other Ochrophyta", if_else(grepl("Chlorophyta", Taxon)==T, "Chlorophyta", phylum)))))) %>%
  filter(is.na(phylum)==F) %>% group_by(Fish.ID, Primer) %>% mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>% aggregate(RRA~phylum+Fish.ID+Primer, ., sum) %>% mutate(cat = if_else(RRA<5, "Other", phylum)) %>%
  #filter(phylum == "Rhodophyta" | phylum == "Bacillariophyta" | phylum == "Ochrophyta" | phylum == "Chlorophyta") %>%
  ggplot()+
  geom_bar(stat = "identity", aes(x = Fish.ID, y = RRA, fill = cat), colour = "black")+
  facet_grid(Primer~., space = "free", scales = "free")+
  scale_fill_manual(values = colours)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7, colour = "black"),
    axis.text.x = element_text(size = 7, colour = "black"), #, angle = 90, vjust = 0.5, hjust=1
    axis.title = element_text(size = 8, colour = "black"),
    legend.title = element_text(size = 7, colour = "black"),
    legend.text = element_text(size = 7, colour = "black"),
    legend.key.size = unit(0.3, "cm"),
    strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
    strip.background = element_rect(fill = "white", color = "black")
  )+
  labs(fill = "Phylum")+
  ylab("Relative read abundance (%)")+
  xlab("Fish ID")

#ggsave(filename = "Phyla_18S_otherOch.jpeg", plot = Phyla_18S, path = "./", device = "jpeg", width = 174, height = 100, units = "mm", dpi = 300)

#mean and median testing
mean_med_18S <- 
df_18S %>% 
  mutate(phylum = if_else(grepl("Rhodophyceae", Taxon)==T, "Rhodophyta", if_else(grepl("Diatomea", Taxon)==T, "Bacillariophyceae", if_else(grepl("Phaeophyceae", Taxon)==T, "Phaeophyceae", if_else(grepl("Ochrophyta", Taxon)==T, "Other Ochrophyta", if_else(grepl("Chlorophyta", Taxon)==T, "Chlorophyta", phylum)))))) %>%
  filter(Primer == "Stoeck") %>%
  filter(is.na(phylum)==F) %>% group_by(Fish.ID, Primer) %>% mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>%
aggregate(RRA~Fish.ID+group, ., sum) 
aggregate(RRA~phylum+Fish.ID+Primer, ., sum) %>% #mutate(phylum = if_else(RRA<5, "Other", phylum)) %>%
  pivot_wider(names_from = Fish.ID, values_from = RRA, values_fill = 0) %>%
  pivot_longer(cols = 3:10, names_to = "Fish.ID", values_to = "RRA") %>%
  group_by(phylum) %>%
  rename(value = RRA) %>%
  select(-Primer)

mean_med_18S <- rbind(mean_med_18S, mean_med_18S %>% summarise(value = mean(value)) %>% mutate(Fish.ID = rep("Mean")) %>% select(phylum, Fish.ID, value), mean_med_18S %>% summarise(value = median(value)) %>% mutate(Fish.ID = rep("Median")) %>% select(phylum, Fish.ID, value))

mean_med_18S_plot <- mean_med_18S %>%
  #get rid of RRA 0s as otherwise they are plotted in the bar graph
  filter(value>0) %>%
  ggplot()+
  geom_bar(stat = "identity", aes(x = factor(Fish.ID, levels = c("CR1","CR2","CR3","MR1","MR2","ND1","SP1","SP2","Median","Mean")), y = value), colour = "black")+
  facet_grid(phylum~., space = "fixed")+ #, scales = "free"
  scale_fill_manual(values = colours)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7, colour = "black"),
    axis.text.x = element_text(size = 7, colour = "black"), #, angle = 90, vjust = 0.5, hjust=1
    axis.title = element_text(size = 8, colour = "black"),
    legend.title = element_text(size = 7, colour = "black"),
    legend.text = element_text(size = 7, colour = "black"),
    legend.key.size = unit(0.3, "cm"),
    strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
    strip.background = element_rect(fill = "white", color = "black")
  )+
  ylab("Relative read abundance (%)")+
  xlab("")

#ggsave(filename = "mean_med_18S.jpeg", plot = mean_med_18S_plot, path = "./", device = "jpeg", width = 174, height = 300, units = "mm", dpi = 300)

#Supplementary Table S7
#overlap of phyla among primers
df_18S %>%
  filter(Primer != "Zim") %>%
  mutate(phylum = if_else(grepl("Rhodophyceae", Taxon)==T, "Rhodophyta", if_else(grepl("Diatomea", Taxon)==T, "Bacillariophyceae", if_else(grepl("Phaeophyceae", Taxon)==T, "Phaeophyceae", if_else(grepl("Ochrophyta", Taxon)==T, "Other Ochrophyta", if_else(grepl("Chlorophyta", Taxon)==T, "Chlorophyta", phylum)))))) %>%
filter(is.na(phylum)==F) %>% select(Fish.ID, Primer, phylum) %>%  unique() %>% group_by(Fish.ID, phylum) %>% 
  mutate(no = length(Primer)) %>% ungroup() %>% select(Fish.ID, phylum, no) %>% unique() %>% pivot_wider(names_from = Fish.ID, values_from = no, values_fill = 0) %>%
  write.table(.,file = "./Sub2-TableS7.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#Supplementary Table S8
#Dinoflag
df_18S %>%
  filter(Primer != "Zim") %>%
  mutate(phylum = if_else(grepl("Rhodophyceae", Taxon)==T, "Rhodophyta", if_else(grepl("Diatomea", Taxon)==T, "Bacillariophyceae", if_else(grepl("Phaeophyceae", Taxon)==T, "Phaeophyceae", if_else(grepl("Ochrophyta", Taxon)==T, "Other Ochrophyta", if_else(grepl("Chlorophyta", Taxon)==T, "Chlorophyta", phylum)))))) %>%
  filter(is.na(phylum)==F) %>% group_by(Fish.ID, Primer) %>% mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>% aggregate(RRA~phylum+Fish.ID+Primer, ., sum) %>% mutate(cat = if_else(RRA<5, "Other", phylum)) %>%
  filter(phylum == "Dinoflagellata") %>%
  select(Fish.ID, Primer, RRA) %>%
  pivot_wider(names_from = Primer, values_from = RRA, values_fill = 0) %>%
  write.table(.,file = "./Sub2-TableS8.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#number of unique phyla per ind/primer
df_18S %>% 
  mutate(phylum = if_else(grepl("Rhodophyceae", Taxon)==T, "Rhodophyta", if_else(grepl("Diatomea", Taxon)==T, "Bacillariophyta", if_else(grepl("Ochrophyta", Taxon)==T, "Ochrophyta", if_else(grepl("Chlorophyta", Taxon)==T, "Chlorophyta", phylum))))) %>%
  filter(is.na(phylum)==F) %>% select(Fish.ID, Primer, phylum) %>%  unique() %>% group_by(Fish.ID, Primer) %>% mutate(no = length(phylum)) %>% ungroup() %>% select(Fish.ID, Primer, no) %>% unique() %>% pivot_wider(names_from = Primer, values_from = no, values_fill = 0) %>% mutate(avg = rowSums(.[,2:4])/3)

#number of individuals per phyla
df_18S %>% 
  mutate(phylum = if_else(grepl("Rhodophyceae", Taxon)==T, "Rhodophyta", if_else(grepl("Diatomea", Taxon)==T, "Bacillariophyta", if_else(grepl("Ochrophyta", Taxon)==T, "Ochrophyta", if_else(grepl("Chlorophyta", Taxon)==T, "Chlorophyta", phylum))))) %>%
  filter(is.na(phylum)==F) %>% select(Fish.ID, Primer, phylum) %>%  unique() %>% group_by(Primer, phylum) %>% mutate(no = length(Fish.ID)/8*100) %>% ungroup() %>% select(phylum, Primer, no) %>% unique() %>% pivot_wider(names_from = Primer, values_from = no, values_fill = 0) %>% View

#Species and genus on one plot
Bubble_GS <- 
df_18S %>% 
  mutate(phylum = if_else(grepl("Rhodophyceae", Taxon)==T, "Rhodophyta", if_else(grepl("Diatomea", Taxon)==T, "Bacillariophyceae", if_else(grepl("Phaeophyceae", Taxon)==T, "Phaeophyceae", if_else(grepl("Ochrophyta", Taxon)==T, "Other Ochrophyta", if_else(grepl("Chlorophyta", Taxon)==T, "Chlorophyta", phylum)))))) %>%
  group_by(Fish.ID, Primer) %>% mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>% ungroup() %>% 
  group_by(Fish.ID, Primer,genus) %>% mutate(decide = if_else(grepl("uncultured", species)==F & is.na(species)==F, "species", "genus")) %>% ungroup() %>%
  group_by(Fish.ID, Primer,decide, genus, species) %>% mutate(ASVs = length(decide)) %>% ungroup() %>%
  group_by(Fish.ID, Primer,species) %>% mutate(ASVs_species = length(species)) %>% ungroup() %>%
  mutate(group = if_else(grepl("uncultured", species)==F & is.na(species)==F, species, genus)) %>% filter(grepl("Rhodophyta|Bacillariophyceae|Other Ochrophyta|Phaeophyceae|Chlorophyta|Dinoflagellata", phylum)==T) %>% 
  aggregate(RRA~Fish.ID+Primer+group+ASVs+phylum, ., sum) %>%
  #mutate(group = gsub("^[a-z].*_..._", "~", group)) %>%
  filter(grepl("^[a-z].*_..._", group)==F) %>%
  mutate(group = gsub("_", " ", group)) %>%
  ggplot()+
  geom_point(shape = 21, aes(x=Primer, y=fct_rev(group), fill = RRA, size = ASVs))+
  facet_grid(phylum~Fish.ID, space = "free", scales = "free")+
  scale_fill_gradientn(colours = c("dark blue", "light blue", "white", "orange", "red"))+
  theme_bw()+
  theme(
    panel.spacing = unit(0, "lines"),
    panel.grid.major.y = element_line(colour = "gray"),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 7, colour = "black"),
    axis.text.x = element_text(size = 7, colour = "black", angle = 90, vjust = 0.5, hjust=1),
    axis.title = element_text(size = 8, colour = "black"),
    legend.title = element_text(size = 7, colour = "black"),
    legend.text = element_text(size = 7, colour = "black"),
    legend.key.size = unit(0.3, "cm"),
    legend.position = "bottom",
    legend.margin = margin(0,0,0,0),
    legend.box.spacing = unit(0, "pt"),
    strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
    strip.text.x = element_text(angle = 0, size = 7, colour = "black"),
    strip.background = element_rect(fill = "white", color = "black")
  )+
  ylab("Genus or species")+
  labs(fill = "RRA (%)", size = "ASVs (n)")

#ggsave(filename = "Bubble_GS_OtherOch.jpeg", plot = Bubble_GS, path = "./", device = "jpeg", width = 174, height = 234, units = "mm", dpi = 300)

#Other taxa
Bubble_GS_nonalgae <- 
df_18S %>% 
  mutate(phylum = if_else(grepl("Rhodophyceae", Taxon)==T, "Rhodophyta", if_else(grepl("Diatomea", Taxon)==T, "Bacillariophyta", if_else(grepl("Ochrophyta", Taxon)==T, "Ochrophyta", if_else(grepl("Chlorophyta", Taxon)==T, "Chlorophyta", phylum))))) %>%
  group_by(Fish.ID, Primer) %>% mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>% ungroup() %>% 
  group_by(Fish.ID, Primer,genus) %>% mutate(decide = if_else(grepl("uncultured", species)==F & is.na(species)==F, "species", "genus")) %>% ungroup() %>%
  group_by(Fish.ID, Primer,decide, genus, species) %>% mutate(ASVs = length(decide)) %>% ungroup() %>%
  group_by(Fish.ID, Primer,species) %>% mutate(ASVs_species = length(species)) %>% ungroup() %>%
  mutate(group = if_else(grepl("uncultured", species)==F & is.na(species)==F, species, genus)) %>%
  filter(grepl("Rhodophyta|Bacillariophyta|Ochrophyta|Chlorophyta|Dinoflagellata", phylum)==F) %>% 
  aggregate(RRA~Fish.ID+Primer+group+ASVs+phylum, ., sum) %>%
  #mutate(group = gsub("^[a-z].*_..._", "~", group)) %>%
  filter(grepl("^[a-z].*_..._", group)==F) %>%
  mutate(group = gsub("_", " ", group)) %>%
  ggplot()+
  geom_point(shape = 21, aes(x=Primer, y=fct_rev(group), fill = RRA, size = ASVs))+
  facet_grid(phylum~Fish.ID, space = "free", scales = "free")+
  scale_fill_gradientn(colours = c("dark blue", "light blue", "white", "orange", "red"))+
  theme_bw()+
  theme(
    panel.spacing = unit(0, "lines"),
    panel.grid.major.y = element_line(colour = "gray"),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 7, colour = "black"),
    axis.text.x = element_text(size = 7, colour = "black", angle = 90, vjust = 0.5, hjust=1),
    axis.title = element_text(size = 8, colour = "black"),
    legend.title = element_text(size = 7, colour = "black"),
    legend.text = element_text(size = 7, colour = "black"),
    legend.key.size = unit(0.3, "cm"),
    strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
    strip.text.x = element_text(angle = 90, size = 7, colour = "black"),
    strip.background = element_rect(fill = "white", color = "black")
  )+
  ylab("Genus or species")+
  labs(fill = "RRA (%)", size = "ASVs (n)")

ggsave(filename = "Bubble_GS_nonalgae.jpeg", plot = Bubble_GS_nonalgae, path = "./", device = "jpeg", width = 180, height = 234, units = "mm", dpi = 300)

####GREP####
#Using grep to investigate sequence matches among 18S primers
#using dna-sequences.fasta_${Primer} for primer.fasta files

#FD v STOECK
#ensure that there is a complete overlap for FD sequences, they need to be trimmed to remove the last three characters
#sed '/^[A-Z]/s/...$//' FD.fasta > FD_rm3char.fasta   
#while IFS='' read -r line; do fastaID=$(echo $line | sed 's/^>.*/YES/' ); if [[ "$fastaID" == "YES" ]]; then ID=$(echo $line); else grep -B 1 $line Stoeck.fasta | grep '^>' | sed 's/$/, /' > match.txt; MATCH=$(cat match.txt | tr -d '\n'); printf "%s\t%s\n" "$ID" "$MATCH" >> FDvSTOECK.txt; fi; done < FD_rm3char.fasta

#FD v ZHAN
#FD lies within ZHAN and so there is no need to trim to assess the overlap
#while IFS='' read -r line; do fastaID=$(echo $line | sed 's/^>.*/YES/' ); if [[ "$fastaID" == "YES" ]]; then ID=$(echo $line); else grep -B 1 $line Zhan.fasta | grep '^>' | sed 's/$/, /' > match.txt; MATCH=$(cat match.txt | tr -d '\n'); printf "%s\t%s\n" "$ID" "$MATCH" >> FDvZHAN.txt; fi; done < FD.fasta

#FD v ZIM
#FD lies within ZIM and so there is no need to trim to assess the overlap
#while IFS='' read -r line; do fastaID=$(echo $line | sed 's/^>.*/YES/' ); if [[ "$fastaID" == "YES" ]]; then ID=$(echo $line); else grep -B 1 $line Zim.fasta | grep '^>' | sed 's/$/, /' > match.txt; MATCH=$(cat match.txt | tr -d '\n'); printf "%s\t%s\n" "$ID" "$MATCH" >> FDvZIM.txt; fi; done < FD.fasta

#Using this output
#Stoeck
FDvStoeck <- read.delim("../grep_exact_match/FDvSTOECK.txt", header = F, "\t")
FDvStoeck_match <- filter(FDvStoeck, V2 != "")
FDvStoeck_nomatch <- filter(FDvStoeck, V2 == "")
FDvStoeck_nomatch <- FDvStoeck_nomatch %>% rename(FD = 1, Match.ID = 2)

FDvStoeck_df <- c()
for (i in 1:nrow(FDvStoeck_match)) {
  x <- FDvStoeck_match[i,]
  match_n = str_count(string = x[,2], pattern = ",")
  
  FD.id = x[,1]
  
  match <- strsplit(x = x[,2], split = ", ")
  match <- data.frame(match) %>% rename(Match.ID = 1)
  
  FDvStoeck_df <- rbind(FDvStoeck_df, cbind(rep(FD.id, match_n), match))
}

FDvStoeck_df <- FDvStoeck_df %>% rename(FD = 1)

FDvStoeck <- rbind(FDvStoeck_nomatch, FDvStoeck_df)

FDvStoeck$Primer <- rep("Stoeck")

#Zhan
FDvZhan <- read.delim("../grep_exact_match/FDvZhan.txt", header = F, "\t")
FDvZhan_match <- filter(FDvZhan, V2 != "")
FDvZhan_nomatch <- filter(FDvZhan, V2 == "")
FDvZhan_nomatch <- FDvZhan_nomatch %>% rename(FD = 1, Match.ID = 2)

FDvZhan_df <- c()
for (i in 1:nrow(FDvZhan_match)) {
  x <- FDvZhan_match[i,]
  match_n = str_count(string = x[,2], pattern = ",")
  
  FD.id = x[,1]
  
  match <- strsplit(x = x[,2], split = ", ")
  match <- data.frame(match) %>% rename(Match.ID = 1)
  
  FDvZhan_df <- rbind(FDvZhan_df, cbind(rep(FD.id, match_n), match))
}

FDvZhan_df <- FDvZhan_df %>% rename(FD = 1)

FDvZhan <- rbind(FDvZhan_nomatch, FDvZhan_df)

FDvZhan$Primer <- rep("Zhan")

#Zim
FDvZim <- read.delim("../grep_exact_match/FDvZim.txt", header = F, "\t")

FDvZim_match <- filter(FDvZim, V2 != "")
FDvZim_nomatch <- filter(FDvZim, V2 == "")
FDvZim_nomatch <- FDvZim_nomatch %>% rename(FD = 1, Match.ID = 2)

FDvZim_df <- c()
for (i in 1:nrow(FDvZim_match)) {
  x <- FDvZim_match[i,]
  match_n = str_count(string = x[,2], pattern = ",")
  
  FD.id = x[,1]
  
  match <- strsplit(x = x[,2], split = ", ")
  match <- data.frame(match) %>% rename(Match.ID = 1)
  
  FDvZim_df <- rbind(FDvZim_df, cbind(rep(FD.id, match_n), match))
}

FDvZim_df <- FDvZim_df %>% rename(FD = 1)

FDvZim <- rbind(FDvZim_nomatch, FDvZim_df)

FDvZim$Primer <- rep("Zim")

#Join the dataframes
Grep <- rbind(FDvStoeck, FDvZhan, FDvZim)

#all feature IDs in 18S dataframe
Feat.ID <- df_18S %>% select(Feature.ID) %>% mutate(Feature.ID = paste(">", Feature.ID, sep = "")) %>% unique()

#retaining only those ASVs which passed filtering
Grep <- Grep %>% filter(FD %in% Feat.ID$Feature.ID)
Grep <- Grep %>% filter(Match.ID %in% Feat.ID$Feature.ID)

Grep <- Grep %>% mutate(Feature.ID = gsub("^>", "", FD))
Grep <- inner_join(Grep, unique(df_18S[,c(2, 7)])) %>% rename(Taxon_FD = Taxon)

Grep <- Grep %>% mutate(Feature.ID = gsub("^>", "", Match.ID))
Grep <- inner_join(Grep, unique(df_18S[,c(2, 7)])) %>% rename(Taxon_Match = Taxon)

Grep <- Grep %>% mutate(Match = if_else(Taxon_FD == Taxon_Match, "T", "F"))

Grep %>% group_by(Primer, FD) %>% mutate(Total = length(Match.ID)) %>% group_by(Primer, FD, Match) %>% mutate(Same = length(Match), Look = if_else(Total == Same, "Ignore", "Look")) %>% aggregate(cbind(Total, Same)~Match+Primer, ., sum)

Grep %>% group_by(Primer, FD) %>% mutate(Total = length(Match.ID)) %>% group_by(Primer, FD, Match) %>% mutate(Same = length(Match), Look = if_else(Total == Same, "Ignore", "Look")) %>% filter(Match == "F") %>% View
