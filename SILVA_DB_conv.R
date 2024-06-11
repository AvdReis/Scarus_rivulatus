library(dplyr)
library(tidyr)

####NOTES####
#taxo.tsv is the tax levels from SILVA
#https://www.arb-silva.de/no_cache/download/archive/release_138.1/Exports/taxonomy/
#tax_slv_ssu_138.1.txt.gz

#Mothur forum
#Some of the code has been copied and modified from
#https://raw.githubusercontent.com/rec3141/diversity-scripts/master/convert_silva_taxonomy.r
############

map.in <- read.table("/path_to/taxo.tsv",header=F,sep="\t",stringsAsFactors=F)
map.in <- map.in[,c(1,3)]
colnames(map.in) <- c("taxlabel","taxlevel")

#fix Escherichia nonsense
map.in$taxlevel[which(map.in$taxlabel=="Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia;")] <- "genus"

map.in$taxlabel <- gsub(";$", "", map.in$taxlabel)
map.in$taxlabel <- gsub(".*;", "", map.in$taxlabel)
map.in$taxlabel <- gsub(" ", "_", map.in$taxlabel)

map.in <- unique(map.in)
map.in %>% group_by(taxlabel) %>% summarise(count = length(taxlabel)) %>% filter(count > 1)

# 1 Incertae Sedis         10
map.in <- filter(map.in, taxlabel != "Incertae_Sedis")
# 5 uncultured              8
map.in <- filter(map.in, taxlabel != "uncultured")

# 2 Labyrinthulomycetes     2
filter(map.in, taxlabel == "Labyrinthulomycetes")
map.in$taxlevel <- if_else(map.in$taxlabel == "Labyrinthulomycetes", "class", map.in$taxlevel)
map.in <- unique(map.in)

# 3 SAR                     2
filter(map.in, taxlabel == "SAR")
map.in$taxlevel <- if_else(map.in$taxlabel == "SAR", "major_clade", map.in$taxlevel)
map.in <- unique(map.in)

# 4 Stramenopiles           2
filter(map.in, taxlabel == "Stramenopiles")
map.in$taxlevel <- if_else(map.in$taxlabel == "Stramenopiles", "kingdom", map.in$taxlevel)
map.in <- unique(map.in)

#change noted in original work
map.in$taxlevel <- if_else(map.in$taxlabel == "RsaHf231", "phylum", map.in$taxlevel)

# bring in the old taxonomic levels from SILVA and remap them using the new levels
#cd /mnt/c/Users/Aimz/Dropbox/PF_NGS/Jan2021/CyanoSeq_SILVA_tax/
#sed 's/ /_/g' SILVA_tax.txt > SILVA_tax_mod.txt
#sed -i 's/;/\t/g' SILVA_tax_mod.txt

tax.in <- read.delim(file = "C:/Users/Aimz/Dropbox/PF_NGS/Jan2021/CyanoSeq_SILVA_tax/SILVA_tax_mod.txt",header=F, sep = "\t")


# library(readxl)
# test <- read_xlsx(path = "C:/Users/Aimz/Dropbox/PF_NGS/Jan2021/CyanoSeq_SILVA_tax/SILVA_tax_mod.xlsx", sheet = 1, col_names = F)
  
tax.in <- 
tax.in %>% pivot_longer(., cols = 2:21, names_to = "order", values_to = "taxlabel") %>% filter(taxlabel != "" & taxlabel != "uncultured" & taxlabel != "Incertae_Sedis") %>% left_join(., map.in)
tax.in <- filter(tax.in, taxlabel != "unidentified")

tax.in <- tax.in %>% mutate(taxlevel = if_else(is.na(taxlevel)==T & grepl("_", taxlabel)==T, "species", taxlevel))
tax.in <- filter(tax.in, is.na(taxlevel)==F)


#selection = c("phylum", "class", "order", "family", "genus", "species")
#tax.in.df <- tax.in %>% filter(., taxlevel %in% selection)
tax.in.df <- tax.in %>% select(-order) %>% unique()
tax.in.df <- aggregate(taxlabel~., tax.in.df, toString)
tax.in.df$taxlabel <- gsub(", ", "#", tax.in.df$taxlabel)

tax.in.df <- tax.in.df %>% pivot_wider(names_from = "taxlevel", values_from = "taxlabel")

colnames(tax.in.df)
                      
tax.in.df <- tax.in.df[,c("V1", "domain", "superkingdom", "kingdom", "subkingdom", "major_clade", "superphylum", "phylum", "subphylum", "infraphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", "superfamily", "family", "subfamily", "genus", "species")]

tax.out <- tax.in.df %>% rename("V1" = V1, "d" = domain, "spk" = superkingdom, "k" = kingdom, "sbk" = subkingdom, "mc" = major_clade, "spp" = superphylum, "p" = phylum, "sbp" = subphylum, "ifp" = infraphylum, "spc" = superclass, "c" = class, "sbc" = subclass, "ifc" = infraclass, "spo" = superorder, "o" = order, "sbo" = suborder, "spf" = superfamily, "f" = family, "sbf" = subfamily, "g" = genus, "s" = species)

#checking no NAs in domain
#unique(tax.out$domain) %>% is.na()

levels_fill <-data.frame(LEVEL_UP = c("d", "spk", "k", "sbk", "mc", "spp", "p", "sbp", "ifp", "spc", "c", "sbc", "ifc", "spo", "o", "sbo", "spf", "f", "sbf", "g"),
                         LEVEL = c("spk", "k", "sbk", "mc", "spp", "p", "sbp", "ifp", "spc", "c", "sbc", "ifc", "spo", "o", "sbo", "spf", "f", "sbf", "g", "s")
)

 test <- tax.out[1:6,]
# test2 <- c()

tax.final <- c()

for(i in 1:nrow(tax.out)) {
  row <- tax.out[i,]
  print(i)
  for (x in levels_fill$LEVEL) {
    if (is.na(row[,c(x)])==T) {
      onedown <- levels_fill %>% filter(LEVEL %in% x) %>% .[1,1]
      row[,c(x)] <- paste(x, "__", row[,c(onedown)])
    }
  }
  tax.final <- rbind(tax.final, row)
}

tax.final <- tax.out
for (i in levels_fill$LEVEL) {
  onedown <- levels_fill %>% filter(LEVEL %in% i) %>% .[1,1]
  tax.final[[i]] <- if_else(is.na(tax.final[[i]])==T, paste(i, tax.final[[onedown]], sep="_"), tax.final[[i]])
}


#For rdp classifier use 6 linnean
SILVA_lin.tax <- tax.final[,c("V1", "p", "c", "o", "f", "g", "s")]
SILVA_lin.tax <- SILVA_lin.tax %>% mutate(tax_all = "")
for (i in 2:ncol(SILVA_lin.tax)) {
  if (i < 7) {
    SILVA_lin.tax$tax_all <- paste(SILVA_lin.tax$tax_all, SILVA_lin.tax[[i]], ";", sep = "")
  }
  
}

#for rpd classifier using all
tax_final_SILVA.tax <- tax.final
tax_final_SILVA.tax <- tax_final_SILVA.tax %>% mutate(tax_all = "")

for (i in 2:ncol(tax_final_SILVA.tax)) {
  if (i < 23) {
    tax_final_SILVA.tax$tax_all <- paste(tax_final_SILVA.tax$tax_all, tax_final_SILVA.tax[[i]], ";", sep = "")
  }
}

write.table(tax_final_SILVA.tax[,c(1,23)], file = "./all_levels_given_SILVA.tax", quote = F, col.names = F, row.names = F)
