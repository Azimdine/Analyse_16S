###Package charging###

library(magrittr)
library(biomformat)
library(lme4)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(grid)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
#devtools::install_github("gauravsk/ranacapa")
library(ranacapa)
library(RColorBrewer)
library(randomcoloR)
library(cowplot)
library(rstatix)
library(RADanalysis)
library(permute)
#remotes::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
library(microbiome)
library(tibble)
library(knitr)
library(nlme)
library(data.table)
library(readstata13)
library(tidyverse)
#library(randomForest)
#library(hrbrthemes)
library(gcookbook)
library(ggConvexHull)
library(ggnewscale)
library(plyr)
library(ggrepel)
library(corrplot)
library(ape)
library(gridExtra)
library(stringr)
library(UpSetR)
library(viridis)
library(ComplexHeatmap)

####First part analysis: charging data, cleaning data and quality controle####

### set the work directory and set working environment###
rm(list=ls(all=TRUE)) # go get a clean workspace

setwd("/Users/ahabib/Documents/B16716_Patrick/partie_1/16S")
#getwd()

taxFile="taxonomy_table.csv"
OTUFile="ASV_table.csv"

OTU_data = read.csv(OTUFile, sep =";") # load the table, with colnames being sample names
# Formating OTU_data to be accepted as input in merge_phyloseq
head(rownames(OTU_data)) # rownames are number, they have to be changed as OTU names
rownames(OTU_data) = OTU_data[,1] # set rownames as OTU names
head(rownames(OTU_data)) # rownames are correct but the 1st column is still present and need to be removed
OTU_data = OTU_data[,-c(1)] # remove OTU names from the table itself

tax_data = read.csv(taxFile, sep =";") # load the table, with colnames being rank
# Formating tax_data to be accepted as input in merge_phyloseq
rownames(tax_data) = tax_data[,1] # set rownames as OTU names
tax_data = tax_data[,-c(1)] # remove OTU names from the table itself
tax_data = as.matrix(tax_data) # transforms as matrix (or else tax_table will do it with a warning)

OTU = otu_table(OTU_data, taxa_are_rows = TRUE)
TAX = tax_table(tax_data)
ps = phyloseq(OTU, TAX)
#sample_names(ps) <- gsub("^X", "", sample_names(ps))
metadata = read.csv("Metadata.csv", sep =";")

# Formating metadata to be accepted as input in merge_phyloseq
metadata[,1] = as.factor(metadata[,1]) # Samples ID's as factor
rownames(metadata) <- metadata[,1] # Assign rownames to be Sample ID's
metadata = metadata[,-c(1)] # Remove the Sample ID's column that is now the name of the rows
metadata = sample_data(metadata) # Transform dataframe into phyloseq class object.

#create phyloseq object#
sample_names(ps)
sample_names(metadata)
phyloseq = merge_phyloseq(ps, metadata)
#The print_ps from microbiomeutilities can give information from the data#
phyloseq
print_ps(phyloseq)

#### Quality control analysis ####


#This will be used to determine further filters.

sort(sample_sums(phyloseq))
colnames(tax_table(phyloseq))

seq_sums <- sort(sample_sums(phyloseq))

# Convertir en DataFrame
seq_sums_df <- data.frame(Sample = names(seq_sums), TotalReads = seq_sums)

# Afficher les premières lignes du tableau
head(seq_sums_df)

# Exporter le tableau sous format CSV
write.csv(seq_sums_df, "sequence_sums_per_sample.csv", row.names = FALSE)

###plot reads per samples###
# Extraire les reads par échantillon
df_reads <- data.frame(
  SampleID = names(sample_sums(phyloseq)),
  Reads = sample_sums(phyloseq)
)

# Ordonner par reads croissants
df_reads$SampleID <- factor(df_reads$SampleID, levels = df_reads$SampleID[order(df_reads$Reads)])

# Plot avec noms visibles
pdf("samplesequeningdepth.pdf", #name of file to print. 
    width = 12, # define plot width and height. completely up to user.
    height = 10)
sampledepth<-ggplot(df_reads, aes(x = SampleID, y = Reads)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) + # rotation et taille
  labs(title = "Nombre de reads par échantillon avant filtration",
       x = "Échantillons",
       y = "Nombre de reads")
sampledepth
dev.off()
ggsave(filename = "sampledepth.png", plot = sampledepth, device = "png", width = 12, height = 10)


#view(tax_data)
# Assumes Kingdom is the first rank and goes as far as it can
#for us, Domain is the first rank
#for us, Domain is the first rank
colnames(tax_table(phyloseq)) <- c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", 
                                   "Species")[1:length(rank_names(phyloseq))] # Renaming taxonomy column names

#Remove unwanted characters, as Eukaryota, Chloroplast and Mitochondria#
Gooddata_rem <- phyloseq %>% subset_taxa(Kingdom != "Archaea"& Family!= "Mitochondria"& Class!="Chloroplast"& Genus!="NA")
Gooddata_rem
print_ps(Gooddata_rem)

###distribution of read/sample
###plot reads per samples###
# Extraire les reads par échantillon
df_reads <- data.frame(
  SampleID = names(sample_sums(Gooddata_rem)),
  Reads = sample_sums(Gooddata_rem)
)

# Ordonner par reads croissants
df_reads$SampleID <- factor(df_reads$SampleID, levels = df_reads$SampleID[order(df_reads$Reads)])

# Plot avec noms visibles
pdf("samplesequeningdepthfiltration.pdf", #name of file to print. 
    width = 12, # define plot width and height. completely up to user.
    height = 10)
sampledepthfil<-ggplot(df_reads, aes(x = SampleID, y = Reads)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) + # rotation et taille
  labs(title = "Nombre de reads par échantillon après filtration",
       x = "Échantillons",
       y = "Nombre de reads")
sampledepthfil
dev.off()
ggsave(filename = "sampledepthfiltration.png", plot = sampledepth, device = "png", width = 12, height = 10)

#joined figures

pdf("Depth_av_ap.pdf", 
    width = 10, 
    height = 12)
plot_grid(sampledepth,sampledepthfil,
          labels=c("A", "B"), ncol = 1, nrow = 2)
dev.off()

###normalization####
##A normalization with rarefaction so each sample will have the same sequence depth###

###methode 1: even depth with the min reads

Gooddata <- rarefy_even_depth(Gooddata_rem, sample.size = min(sample_sums(Gooddata_rem)), rngseed = 123)
Gooddata
print_ps(Gooddata)

##methode 2: CSS pour garder tous les échantillon
# Méthode CSS complète
#library(metagenomeSeq)

# Conversion
#mgseq <- phyloseq_to_metagenomeSeq(Gooddata_rem)

# Calcul du paramètre de normalisation
#p <- cumNormStatFast(mgseq)
#mgseq <- cumNorm(mgseq, p = p)

# Extraction des données normalisées
#norm_counts <- MRcounts(mgseq, norm = TRUE, log = FALSE)

# Creation du nouvel objet phyloseq
#Gooddata_CSS <- Gooddata_rem
#otu_table(Gooddata_CSS) <- otu_table(norm_counts, taxa_are_rows = TRUE)

# Vérification
#sort(sample_sums(Gooddata_CSS))  # Maintenant normalisé !
#Gooddata <-Gooddata_CSS
#print_ps(Gooddata)
#summarize_phyloseq(Gooddata)

###methode 3: cutoff: <= 5000 reads

#samplerar<- rarefy_even_depth(phyloseq, sample.size = 5000, verbose = FALSE, replace = FALSE,rngseed = 123)

#Gooddata <-samplerar

####Second part of the analysis: relative abundance, alpha and beta diversity####

### Gross analysis to see who is there in the dataset###
##look how much taxa compose the full dataset##
dfBacteria = subset_taxa(Gooddata, Kingdom=="Bacteria")
get_taxa_unique(dfBacteria, "Phylum") # to see how many Subdivision, we have 21 Subdivision
get_taxa_unique(dfBacteria, "Class") # to see how many genera, 42 different Class
get_taxa_unique(dfBacteria, "Order") # to see how many genera, 73 different Order
get_taxa_unique(dfBacteria, "Family") # to see how many families, 96 different families
get_taxa_unique(dfBacteria, "Genus") # to see how many genera, 151 different genera
get_taxa_unique(dfBacteria, "Species") # to see how many species, 156 different species/OTU are present

####Relative abundance####
###Prevalence abundance###

## Summarize taxa information ##
#Get taxa summary
#This can be used for entire dataset.
#phylum level
p0 <- Gooddata
tx.sum <- taxa_summary(p0, "Phylum")
tx.sum
write.csv(tx.sum, "prevalence_Phylum.csv")

##Arranging metadata##
sample_data(Gooddata)$Well<- factor((sample_data(Gooddata)$Well), levels=c("3B_8","3B_11","3W_8","3W_11"))
sample_data(Gooddata)$Samples<- factor((sample_data(Gooddata)$Samples), levels=c("B8","B11","W8","W11"))
###Plot creation abundance###

###Phylum###

ps2.rel <- microbiome::transform(Gooddata, "compositional")

# collapse at Phylum level
ps2.Phylum.rel <- tax_glom(physeq = ps2.rel, taxrank = "Phylum", NArm = FALSE)

# add Phylum names (convert to character), replace NA with "Unclassified"
Phylum <- as.character(tax_table(ps2.Phylum.rel)[, "Phylum"])
Phylum[is.na(Phylum)] <- "Unclassified Phylum"
taxa_names(ps2.Phylum.rel) <- make.unique(Phylum)

# create dataframe
ps2.Phylum.rel.df <- data.table(psmelt(ps2.Phylum.rel))

# calculate mean relative abundance by Phylum
ps2.Phylum.rel.df[, mean := mean(Abundance, na.rm = TRUE), by = "OTU"]
write.csv(ps2.Phylum.rel.df, "Phylum_abundance.csv")

# merge less abundant taxa (<0.1%)
#ps2.Phylum.rel.df[(mean <= 0.001), OTU := "Less Abundant Phylum"]

# Calcul de l'abondance totale pour garder les 20 family les plus abondants
top_phylum <- ps2.Phylum.rel.df[, .(total_abundance = sum(Abundance)), by = Phylum]
top_phylum <- top_phylum[order(-total_abundance)][1:10]$Phylum  

# Regroupe tous les autres genres sous "Less Abundant family"
ps2.Phylum.rel.df[!Phylum %in% top_phylum, Phylum := "Less Abundant Phylum"]

# Creating df with summarized lesser abundant taxa abundance
ps2.Phylum.rel.df <- ps2.Phylum.rel.df[, sum(Abundance), by = list(Phylum, Sample, Well,Samples)]
colnames(ps2.Phylum.rel.df)[5] <- "Abundance"

# define color palette
colcodes.Phylum <- distinctColorPalette(length(unique(ps2.Phylum.rel.df$OTU)) + 13)


#write.csv(ps2.Phylum.rel.df, "Phylum_abundance_less.csv")

###Plot creation abundance###

pdf("Phylum-barplot_all.pdf", 
    width = 10, 
    height = 12)
Phylum.bar <- ggplot(data=ps2.Phylum.rel.df, aes(x=Sample, y=Abundance*100, fill=Phylum)) + 
  geom_bar(position="stack", stat="identity", color = "black", linewidth = 0.05, width = 1)  + 
  xlab("Samples") + ylab("Relative abundance of Phylum")   + scale_fill_manual(values = colcodes.Phylum) + 
  labs(fill = "Phylum")+ coord_cartesian(expand = FALSE ) + guides(fill=guide_legend(ncol=1)) + 
  scale_y_continuous(breaks = c(seq(from = 0, to = 100, by = 5)))+theme(legend.text=element_text(size=13),legend.title = element_text(size=13, face="bold"))+
  theme(
    axis.text.x = element_text(size = 10, angle = 0),
    axis.text.y = element_text(size = 10, angle = 45, vjust = 0.5, hjust = 1),
    axis.title = element_text(size=13, face="bold"),
    plot.title = element_text(size=15, face="bold", margin = ggplot2::margin(t = 0, r = 0, b = 10)),
    strip.text.x = element_text(size=10, face="bold"),
    legend.text = element_text(face = "italic")
  )
Phylum.bar
dev.off()
ggsave(filename = "Phylum-barplot_all.png", plot = Phylum.bar, device = "png", width = 12, height = 10)

###Well##

taxa_p <- tax_glom(Gooddata, "Phylum")
taxa_p_melt = merge_samples(taxa_p,"Well")
taxa_p_melt = transform_sample_counts(taxa_p_melt, function(x) 100 * x/sum(x))

taxa_p_s = taxa_p_melt %>%
  psmelt() %>%
  arrange(Phylum)

# CORRECTION : Vérifions la structure des données
print(str(taxa_p_s))

Phylum_totals <- aggregate(Abundance ~ Phylum, data = taxa_p_s, 
                          FUN = sum, na.rm = TRUE)

# Renommer la colonne et trier
names(Phylum_totals)[2] <- "total_abundance"
Phylum_totals <- Phylum_totals[order(-Phylum_totals$total_abundance), ]

print(head(Phylum_totals, 10))

# Sélectionner des Phylum les plus abondants
top_30_Phylum <- head(Phylum_totals$Phylum,10)

print(top_30_Phylum)  # Vérifions la liste

# Créer une nouvelle variable
taxa_p_s_top30 <- taxa_p_s %>%
  mutate(Phylum_taxa = ifelse(Phylum %in% top_30_Phylum, as.character(Phylum), 
                                   "Less Abundant Phylum"))

write.csv(taxa_p_s_top30, "Phylum_less.csv")

# Vérifier le résultat
table(taxa_p_s_top30$Phylum_taxa, useNA = "always")
print(paste("Nombre de catégories:", length(unique(taxa_p_s_top30$Phylum_taxa))))

# Vérifier le nombre de catégories uniques
length(unique(taxa_p_s_top30$Phylum_taxa)) # Devrait être 31 (30 + "Autres")

# Adapter votre palette de couleurs pour 31 catégories
# (Assurez-vous que colcodes.genus a au moins 31 couleurs)

# Modifier le plot
pdf("Phylum-barplot_Well.pdf", width = 10, height = 9)
Phylum_barplot <- ggplot(taxa_p_s_top30, aes(x = Sample, y = Abundance, fill = Phylum_taxa)) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(0.2, "lines")) +
  theme(strip.text = element_text(colour = 'black'), 
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=12),
        legend.text = element_text(face = "italic"), 
        axis.text.y = element_text(angle=90, hjust = 1, size=12), 
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=12, face="bold",
                                    margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(angle = 90, hjust = 0.5, size=12, face="bold", 
                                    margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)), 
        plot.title = element_text(angle = 0, hjust = 0.5, size=12, face="bold", 
                                  margin = ggplot2::margin(t = 0, r = 0, b = 20, l = 0))) +
  geom_bar(stat = "identity", width = 0.8,linewidth = 0.05) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_fill_manual(values = colcodes.Phylum) + # Assurez-vous d'avoir 31 couleurs
  theme(legend.text=element_text(size=12), 
        legend.title = element_text(size=12, face="bold")) + 
  ylab("Relative Abundance") +
  xlab("Well") +
  ggtitle("Relative Abundance of Phylum according to Well")

Phylum_barplot
dev.off()
ggsave(filename = "Phylum-barplot_Wel.png", plot = Phylum_barplot, device = "png", width = 10, height = 9)

###Samples##

taxa_p <- tax_glom(Gooddata, "Phylum")
taxa_p_melt = merge_samples(taxa_p,"Samples")
taxa_p_melt = transform_sample_counts(taxa_p_melt, function(x) 100 * x/sum(x))

taxa_p_s = taxa_p_melt %>%
  psmelt() %>%
  arrange(Phylum)

# CORRECTION : Vérifions la structure des données
print(str(taxa_p_s))

Phylum_totals <- aggregate(Abundance ~ Phylum, data = taxa_p_s, 
                           FUN = sum, na.rm = TRUE)

# Renommer la colonne et trier
names(Phylum_totals)[2] <- "total_abundance"
Phylum_totals <- Phylum_totals[order(-Phylum_totals$total_abundance), ]

print(head(Phylum_totals, 10))

# Sélectionner des Phylum les plus abondants
top_30_Phylum <- head(Phylum_totals$Phylum,10)

print(top_30_Phylum)  # Vérifions la liste

# Créer une nouvelle variable
taxa_p_s_top30 <- taxa_p_s %>%
  mutate(Phylum_taxa = ifelse(Phylum %in% top_30_Phylum, as.character(Phylum), 
                              "Less Abundant Phylum"))

write.csv(taxa_p_s_top30, "Phylum_less_Samples.csv")

# Vérifier le résultat
table(taxa_p_s_top30$Phylum_taxa, useNA = "always")
print(paste("Nombre de catégories:", length(unique(taxa_p_s_top30$Phylum_taxa))))

# Vérifier le nombre de catégories uniques
length(unique(taxa_p_s_top30$Phylum_taxa)) # Devrait être 31 (30 + "Autres")

# Adapter votre palette de couleurs pour 31 catégories
# (Assurez-vous que colcodes.genus a au moins 31 couleurs)

# Modifier le plot
pdf("Phylum-barplot_Samples.pdf", width = 10, height = 9)
Phylum_barplot <- ggplot(taxa_p_s_top30, aes(x = Sample, y = Abundance, fill = Phylum_taxa)) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(0.2, "lines")) +
  theme(strip.text = element_text(colour = 'black'), 
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=12),
        legend.text = element_text(face = "italic"), 
        axis.text.y = element_text(angle=90, hjust = 1, size=12), 
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=12, face="bold",
                                    margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(angle = 90, hjust = 0.5, size=12, face="bold", 
                                    margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)), 
        plot.title = element_text(angle = 0, hjust = 0.5, size=12, face="bold", 
                                  margin = ggplot2::margin(t = 0, r = 0, b = 20, l = 0))) +
  geom_bar(stat = "identity", width = 0.8,linewidth = 0.05) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_fill_manual(values = colcodes.Phylum) + # Assurez-vous d'avoir 31 couleurs
  theme(legend.text=element_text(size=12), 
        legend.title = element_text(size=12, face="bold")) + 
  ylab("Relative Abundance") +
  xlab("Samples") +
  ggtitle("Relative Abundance of Phylum according to Samples")

Phylum_barplot
dev.off()
ggsave(filename = "Phylum-barplot_Samples.png", plot = Phylum_barplot, device = "png", width = 10, height = 9)

### Family ###

# create a object with relative abundance
ps2.rel3 <- microbiome::transform(Gooddata, "compositional")

# define the levels to glom 
ps2.Family.rel <- tax_glom(physeq = ps2.rel3, taxrank = "Family", NArm = FALSE)

# CORRECTION : Gérer les noms de taxons dupliqués
# Au lieu d'utiliser directement les noms de famille, créons des noms uniques

# Extraire les noms de famille
family_names <- tax_table(ps2.Family.rel)[, "Family"]

# Remplacer les NA par "Unclassified Family"
family_names[is.na(family_names)] <- "Unclassified Family"

# Vérifier les doublons
duplicated_families <- duplicated(family_names) | duplicated(family_names, fromLast = TRUE)
if(any(duplicated_families)) {
  cat("Familles avec des doublons détectées:\n")
  print(unique(family_names[duplicated_families]))
}

# Créer des noms uniques en ajoutant un identifiant numérique pour les doublons
unique_family_names <- make.unique(family_names, sep = "_")

# Assigner les noms uniques
taxa_names(ps2.Family.rel) <- unique_family_names

# Créer une colonne avec le nom de famille propre pour l'affichage
tax_table(ps2.Family.rel) <- cbind(tax_table(ps2.Family.rel), 
                                   Display_Family = family_names)

# create dataframe
ps2.Family.rel.df <- data.table(psmelt(ps2.Family.rel))

# group df by taxa and calculate mean rel. abundance
ps2.Family.rel.df[, mean := mean(Abundance, na.rm = TRUE), by = "OTU"]

# Sauvegarder le DataFrame en CSV
write.csv(ps2.Family.rel.df, "Family_abundance.csv")

# Calcul de l'abondance totale pour garder les 20 family les plus abondants
top_family <- ps2.Family.rel.df[, .(total_abundance = sum(Abundance)), by = Family]
top_family <- top_family[order(-total_abundance)][1:10]$Family  

# Regroupe tous les autres genres sous "Less Abundant family"
ps2.Family.rel.df[!Family %in% top_family, Family := "Less Abundant Family"]

# Creating df with summarized lesser abundant taxa abundance
ps2.Family.rel.df <- ps2.Family.rel.df[, sum(Abundance), by = list(Family, Sample, Well,Samples)]
colnames(ps2.Family.rel.df)[5] <- "Abundance"
# get color codes - utiliser Display_Family
colcodes.Family <- distinctColorPalette(length(unique(ps2.Family.rel.df$Display_Family)) + 13)

#plotting stacked bar chart
pdf("Family-barplot_all.pdf", 
    width = 10, 
    height = 8)
Family.bar <- ggplot(data=ps2.Family.rel.df, aes(x=Sample, y=Abundance*100, fill=Family)) + 
  geom_bar(position="stack", stat="identity", color = "black", linewidth = 0.05, width = 1)  + 
  xlab("Samples") + ylab("Relative abundance of Family")   + scale_fill_manual(values = colcodes.Family) + 
  labs(fill = "Family")+ coord_cartesian(expand = FALSE ) + guides(fill=guide_legend(ncol=1)) + 
  scale_y_continuous(breaks = c(seq(from = 0, to = 100, by = 5)))+theme(legend.text=element_text(size=13),legend.title = element_text(size=13, face="bold"))+
  theme(
    axis.text.x = element_text(size = 10, angle = 0),
    axis.text.y = element_text(size = 10, angle = 45, vjust = 0.5, hjust = 1),
    axis.title = element_text(size=13, face="bold"),
    plot.title = element_text(size=15, face="bold", margin = ggplot2::margin(t = 0, r = 0, b = 10)),
    strip.text.x = element_text(size=10, face="bold"),
    legend.text = element_text(face = "italic")
  )
Family.bar
dev.off()
ggsave(filename = "Family-bar-plot_all.png", plot = Family.bar, device = "png", width = 10, height = 8)

###Well##

taxa_p <- tax_glom(Gooddata, "Family")
taxa_p_melt = merge_samples(taxa_p,"Well")
taxa_p_melt = transform_sample_counts(taxa_p_melt, function(x) 100 * x/sum(x))

taxa_p_s = taxa_p_melt %>%
  psmelt() %>%
  arrange(Family)

# CORRECTION : Vérifions la structure des données
print(str(taxa_p_s))

family_totals <- aggregate(Abundance ~ Family, data = taxa_p_s, 
                           FUN = sum, na.rm = TRUE)

# Renommer la colonne et trier
names(family_totals)[2] <- "total_abundance"
family_totals <- family_totals[order(-family_totals$total_abundance), ]

print(head(family_totals, 10))

# Sélectionner les 30 familles les plus abondants
top_30_family <- head(family_totals$Family, 10)

print(top_30_family)  # Vérifions la liste

# Créer une nouvelle variable
taxa_p_s_top30 <- taxa_p_s %>%
  mutate(Family_taxa = ifelse(Family %in% top_30_family, as.character(Family), "Other family"))

write.csv(taxa_p_s_top30, "Family_less_Well.csv")

# Vérifier le résultat
table(taxa_p_s_top30$Family_taxa, useNA = "always")
print(paste("Nombre de catégories:", length(unique(taxa_p_s_top30$Family_taxa))))

# Vérifier le nombre de catégories uniques
length(unique(taxa_p_s_top30$Family_taxa)) # Devrait être 31 (30 + "Autres")

# Adapter votre palette de couleurs pour 31 catégories
# (Assurez-vous que colcodes.genus a au moins 31 couleurs)

# Modifier le plot
pdf("Family-barplot_Well.pdf", width = 8, height = 10)
Family_barplot <- ggplot(taxa_p_s_top30, aes(x = Sample, y = Abundance, fill = Family_taxa)) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(0.2, "lines")) +
  theme(strip.text = element_text(colour = 'black'), 
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=12),
        legend.text = element_text(face = "italic"), 
        axis.text.y = element_text(angle=90, hjust = 1, size=12), 
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=12, face="bold",
                                    margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(angle = 90, hjust = 0.5, size=12, face="bold", 
                                    margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)), 
        plot.title = element_text(angle = 0, hjust = 0.5, size=12, face="bold", 
                                  margin = ggplot2::margin(t = 0, r = 0, b = 20, l = 0))) +
  geom_bar(stat = "identity", width = 0.8,linewidth = 0.05) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_fill_manual(values = colcodes.Family) + # Assurez-vous d'avoir 31 couleurs
  theme(legend.text=element_text(size=12), 
        legend.title = element_text(size=12, face="bold")) + 
  ylab("Relative Abundance") +
  xlab("Well") +
  ggtitle("Relative Abundance of Top 10 Family according to Well")

Family_barplot
dev.off()
ggsave(filename = "Family_barplot_Well.png", plot = Family_barplot, device = "png", width = 8, height = 9)

###Samples##

taxa_p <- tax_glom(Gooddata, "Family")
taxa_p_melt = merge_samples(taxa_p,"Samples")
taxa_p_melt = transform_sample_counts(taxa_p_melt, function(x) 100 * x/sum(x))

taxa_p_s = taxa_p_melt %>%
  psmelt() %>%
  arrange(Family)

# CORRECTION : Vérifions la structure des données
print(str(taxa_p_s))

family_totals <- aggregate(Abundance ~ Family, data = taxa_p_s, 
                           FUN = sum, na.rm = TRUE)

# Renommer la colonne et trier
names(family_totals)[2] <- "total_abundance"
family_totals <- family_totals[order(-family_totals$total_abundance), ]

print(head(family_totals, 10))

# Sélectionner les 30 familles les plus abondants
top_30_family <- head(family_totals$Family, 10)

print(top_30_family)  # Vérifions la liste

# Créer une nouvelle variable
taxa_p_s_top30 <- taxa_p_s %>%
  mutate(Family_taxa = ifelse(Family %in% top_30_family, as.character(Family), "Other family"))

write.csv(taxa_p_s_top30, "Family_less_Samples.csv")

# Vérifier le résultat
table(taxa_p_s_top30$Family_taxa, useNA = "always")
print(paste("Nombre de catégories:", length(unique(taxa_p_s_top30$Family_taxa))))

# Vérifier le nombre de catégories uniques
length(unique(taxa_p_s_top30$Family_taxa)) # Devrait être 31 (30 + "Autres")

# Adapter votre palette de couleurs pour 31 catégories
# (Assurez-vous que colcodes.genus a au moins 31 couleurs)

# Modifier le plot
pdf("Family-barplot_Samples.pdf", width = 8, height = 10)
Family_barplot <- ggplot(taxa_p_s_top30, aes(x = Sample, y = Abundance, fill = Family_taxa)) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(0.2, "lines")) +
  theme(strip.text = element_text(colour = 'black'), 
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=12),
        legend.text = element_text(face = "italic"), 
        axis.text.y = element_text(angle=90, hjust = 1, size=12), 
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=12, face="bold",
                                    margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(angle = 90, hjust = 0.5, size=12, face="bold", 
                                    margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)), 
        plot.title = element_text(angle = 0, hjust = 0.5, size=12, face="bold", 
                                  margin = ggplot2::margin(t = 0, r = 0, b = 20, l = 0))) +
  geom_bar(stat = "identity", width = 0.8,linewidth = 0.05) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_fill_manual(values = colcodes.Family) + # Assurez-vous d'avoir 31 couleurs
  theme(legend.text=element_text(size=12), 
        legend.title = element_text(size=12, face="bold")) + 
  ylab("Relative Abundance") +
  xlab("Samples") +
  ggtitle("Relative Abundance of Top 10 Family according to Samples")

Family_barplot
dev.off()
ggsave(filename = "Family_barplot_Samples.png", plot = Family_barplot, device = "png", width = 8, height = 9)

#####Genus#####

ps2.rel2 <- microbiome::transform(Gooddata, "compositional")
ps2.Genus.rel <- tax_glom(physeq = ps2.rel2, taxrank = "Genus", NArm = FALSE)

#melt to dataframe and convert factor to character
ps2.Genus.rel.df <- data.table(psmelt(ps2.Genus.rel))
write.csv(ps2.Genus.rel.df, "genus_all.csv")

ps2.Genus.rel.df$Phylum <- as.character(ps2.Genus.rel.df$Phylum)
ps2.Genus.rel.df$Class <- as.character(ps2.Genus.rel.df$Class)
ps2.Genus.rel.df$Order <- as.character(ps2.Genus.rel.df$Order)
ps2.Genus.rel.df$Family <- as.character(ps2.Genus.rel.df$Family)
ps2.Genus.rel.df$Genus <- as.character(ps2.Genus.rel.df$Genus)

# Remplacer les NA dans la colonne Genus par "Unassigned"
ps2.Genus.rel.df[is.na(Genus), Genus := "Unassigned"]

# Calcul de la moyenne et du maximum d'abondance
ps2.Genus.rel.df[, mean := ((sum(Abundance, na.rm = TRUE)) / 4), by = "Genus"]
ps2.Genus.rel.df[, maximum := max(Abundance[Abundance > 0], na.rm = TRUE), by = "Genus"]

# Calcul de l'abondance totale pour garder les 20 genres les plus abondants
top_genres <- ps2.Genus.rel.df[, .(total_abundance = sum(Abundance)), by = Genus]
top_genres <- top_genres[order(-total_abundance)][1:10]$Genus  # Garde les 30 genres les plus abondants

# Regroupe tous les autres genres sous "Less Abundant Genera"
ps2.Genus.rel.df[!Genus %in% top_genres, Genus := "Other Genera"]
##verifier si les genres sont uniques##
unique(ps2.Genus.rel.df$Genus)
# Filtrer pour enlever les lignes où Genus est NA et les valeurs vides
ps2.Genus.rel.df <- ps2.Genus.rel.df[!is.na(Genus)]
ps2.Genus.rel.df$Genus[ps2.Genus.rel.df$Genus == ""] <- "Unclassified"
#ps2.Genus.rel.df <- ps2.Genus.rel.df[ps2.Genus.rel.df$Genus != "", ]
# Agréger les données
ps2.Genus.rel.df <- ps2.Genus.rel.df[, sum(Abundance), by = list(Genus, Sample, Well,Samples)]
colnames(ps2.Genus.rel.df)[5] <- "Abundance"

# Génération des couleurs uniques pour chaque genre
colcodes.genus <- distinctColorPalette(length(unique(ps2.Genus.rel.df$Genus)) + 7)

# Visualisation en graphique en barres empilées
pdf("Genus-barplot_all.pdf", width = 10, height = 12)
Genus.bar <- ggplot(data=ps2.Genus.rel.df, aes(x=Sample, y=Abundance*100, fill=Genus)) + 
  geom_bar(position="stack", stat="identity", color = "black", size=0.05, width = 1)  + 
  xlab("Samples") + ylab("Relative abundance of Genera")   + 
  scale_fill_manual(values = colcodes.genus) + labs(fill = "Genus") + 
  coord_cartesian(expand = FALSE, ylim = c(0, 100) ) + guides(fill=guide_legend(ncol=1)) + 
  scale_y_continuous(breaks = c(seq(from = 0, to = 100, by = 5)))+theme(legend.text=element_text(size=13),legend.title = element_text(size=13, face="bold"))+
  theme(
    axis.text.x = element_text(size = 10, angle = 0),
    axis.text.y = element_text(size = 10, angle = 45, vjust = 0.5, hjust = 1),
    axis.title = element_text(size=13, face="bold"),
    plot.title = element_text(size=15, face="bold", margin = ggplot2::margin(t = 0, r = 0, b = 10)),
    strip.text.x = element_text(size=10, face="bold"),
    legend.text = element_text(face = "italic")
  )
Genus.bar
dev.off()
ggsave(filename = "Genus-barplot_all.png", plot = Genus.bar, device = "png", width = 10, height = 12)

##Well##

taxa_p <- tax_glom(Gooddata, "Genus")
taxa_p_melt = merge_samples(taxa_p,"Well")
taxa_p_melt = transform_sample_counts(taxa_p_melt, function(x) 100 * x/sum(x))

taxa_p_s = taxa_p_melt %>%
  psmelt() %>%
  arrange(Genus)

# CORRECTION : Vérifions la structure des données
print(str(taxa_p_s))

genus_totals <- aggregate(Abundance ~ Genus, data = taxa_p_s, 
                          FUN = sum, na.rm = TRUE)

# Renommer la colonne et trier
names(genus_totals)[2] <- "total_abundance"
genus_totals <- genus_totals[order(-genus_totals$total_abundance), ]

print(head(genus_totals, 10))

# Sélectionner les 30 genres les plus abondants
top_30_genus <- head(genus_totals$Genus, 10)

print(top_30_genus)  # Vérifions la liste

# Créer une nouvelle variable
taxa_p_s_top30 <- taxa_p_s %>%
  mutate(Genus_taxa = ifelse(Genus %in% top_30_genus, as.character(Genus), "Other genera"))

write.csv(taxa_p_s_top30, "Genus_less_Well.csv")

# Vérifier le résultat
table(taxa_p_s_top30$Genus_taxa, useNA = "always")
print(paste("Nombre de catégories:", length(unique(taxa_p_s_top30$Genus_taxa))))

# Vérifier le nombre de catégories uniques
length(unique(taxa_p_s_top30$Genus_taxa)) # Devrait être 31 (30 + "Autres")

# Adapter votre palette de couleurs pour 31 catégories
# (Assurez-vous que colcodes.genus a au moins 31 couleurs)

# Modifier le plot
pdf("Genus-barplot_Well.pdf", width = 10, height = 8)
Genus_barplot <- ggplot(taxa_p_s_top30, aes(x = Sample, y = Abundance, fill = Genus_taxa)) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(0.2, "lines")) +
  theme(strip.text = element_text(colour = 'black'), 
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=10),
        legend.text = element_text(face = "italic"), 
        axis.text.y = element_text(angle=90, hjust = 1, size=12), 
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=12, face="bold",
                                    margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(angle = 90, hjust = 0.5, size=12, face="bold", 
                                    margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)), 
        plot.title = element_text(angle = 0, hjust = 0.5, size=12, face="bold", 
                                  margin = ggplot2::margin(t = 0, r = 0, b = 20, l = 0))) +
  geom_bar(stat = "identity", width = 0.8,linewidth = 0.05) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_fill_manual(values = colcodes.genus) + # Assurez-vous d'avoir 31 couleurs
  theme(legend.text=element_text(size=12), 
        legend.title = element_text(size=12, face="bold")) + 
  ylab("Relative Abundance") +
  xlab("Well") +
  ggtitle("Relative Abundance of Top 10 Genera according to Well")

Genus_barplot
dev.off()
ggsave(filename = "Genus-barplot_Well.png", plot = Genus_barplot, device = "png", width = 10, height = 8)

##Samples##

taxa_p <- tax_glom(Gooddata, "Genus")
taxa_p_melt = merge_samples(taxa_p,"Samples")
taxa_p_melt = transform_sample_counts(taxa_p_melt, function(x) 100 * x/sum(x))

taxa_p_s = taxa_p_melt %>%
  psmelt() %>%
  arrange(Genus)

# CORRECTION : Vérifions la structure des données
print(str(taxa_p_s))

genus_totals <- aggregate(Abundance ~ Genus, data = taxa_p_s, 
                          FUN = sum, na.rm = TRUE)

# Renommer la colonne et trier
names(genus_totals)[2] <- "total_abundance"
genus_totals <- genus_totals[order(-genus_totals$total_abundance), ]

print(head(genus_totals, 10))

# Sélectionner les 30 genres les plus abondants
top_30_genus <- head(genus_totals$Genus, 10)

print(top_30_genus)  # Vérifions la liste

# Créer une nouvelle variable
taxa_p_s_top30 <- taxa_p_s %>%
  mutate(Genus_taxa = ifelse(Genus %in% top_30_genus, as.character(Genus), "Other genera"))

write.csv(taxa_p_s_top30, "Genus_less_Samples.csv")

# Vérifier le résultat
table(taxa_p_s_top30$Genus_taxa, useNA = "always")
print(paste("Nombre de catégories:", length(unique(taxa_p_s_top30$Genus_taxa))))

# Vérifier le nombre de catégories uniques
length(unique(taxa_p_s_top30$Genus_taxa)) # Devrait être 31 (30 + "Autres")

# Adapter votre palette de couleurs pour 31 catégories
# (Assurez-vous que colcodes.genus a au moins 31 couleurs)

# Modifier le plot
pdf("Genus-barplot_Samples.pdf", width = 10, height = 8)
Genus_barplot <- ggplot(taxa_p_s_top30, aes(x = Sample, y = Abundance, fill = Genus_taxa)) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(0.2, "lines")) +
  theme(strip.text = element_text(colour = 'black'), 
        axis.text.x = element_text(angle = 0, hjust = 0.5, size=10),
        legend.text = element_text(face = "italic"), 
        axis.text.y = element_text(angle=90, hjust = 1, size=12), 
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=1, size=12, face="bold",
                                    margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(angle = 90, hjust = 0.5, size=12, face="bold", 
                                    margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)), 
        plot.title = element_text(angle = 0, hjust = 0.5, size=12, face="bold", 
                                  margin = ggplot2::margin(t = 0, r = 0, b = 20, l = 0))) +
  geom_bar(stat = "identity", width = 0.8,linewidth = 0.05) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_fill_manual(values = colcodes.genus) + # Assurez-vous d'avoir 31 couleurs
  theme(legend.text=element_text(size=12), 
        legend.title = element_text(size=12, face="bold")) + 
  ylab("Relative Abundance") +
  xlab("Samples") +
  ggtitle("Relative Abundance of Top 10 Genera according to Samples")

Genus_barplot
dev.off()
ggsave(filename = "Genus-barplot_Samples.png", plot = Genus_barplot, device = "png", width = 10, height = 8)

### ALPHA DIVERSITY: CHAO1 AND InVSIMPSON AND BETA DIVERSITY PCOA PLOTS####

## get a general overview of the richness according infestations state and all parasites variables ##
# make a vector with the alpha diversity estimators we want to calculate

# Vérifier le nombre total de reads et de taxons
cat("Nombre total de reads:", sum(sample_sums(Gooddata)), "\n")
cat("Nombre total de taxons:", ntaxa(Gooddata), "\n")
cat("Nombre minimum de reads par échantillon:", min(sample_sums(Gooddata)), "\n")

# Vérifier la présence de singletons
singletons <- taxa_sums(Gooddata) == 1
cat("Nombre de singletons:", sum(singletons), "\n")

###plotting Alpha diversity####


##Well####

alpha_Gooddata = c("Chao1", "Shannon","InvSimpson")

# adding plotting information to object "D"
D <- plot_richness(phyloseq, "Well", measures=alpha_Gooddata, color="Well")

# plot data from "D" as a boxplot with ggplot2 plus wilcox.test statistical analysis
#comp <- make_pairs(sample_data(phyloseq)$Well)

pdf("Alpha-div_Well.pdf", 
    width = 8, 
    height = 10)
Deg<-D + geom_boxplot(data=D$data , aes(x=Well, color=NULL))+ 
  #stat_compare_means(
    #comparisons = comp,
    #label = "p.format",
    #tip.length = 0.05,
    #method = "wilcox.test",
    #p.adjust.methods="holm")+ 
  theme(axis.text=element_text(size=14,angle = 0, hjust = 0.5, vjust=0.5),
        axis.title = element_text(size=14, face="bold"), plot.title =element_text(size=14, face="bold", 
                                                                                  margin =ggplot2:: margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=14, face="bold" )) +
  xlab("Well")+theme(legend.text=element_text(size=14),
                         legend.title = element_text(size=14, face="bold"))+
  theme(legend.position="none")+theme(axis.text.x = element_text(angle=0, hjust = 0.5, vjust=0.5))
Deg
dev.off()
ggsave(filename = "Alpha-div_Well.png", plot = Deg, device = "png",width = 8, height = 10)

##Samples####

alpha_Gooddata = c("Chao1", "Shannon","InvSimpson")

# adding plotting information to object "D"
D <- plot_richness(phyloseq, "Samples", measures=alpha_Gooddata, color="Samples")

# plot data from "D" as a boxplot with ggplot2 plus wilcox.test statistical analysis
#comp <- make_pairs(sample_data(phyloseq)$Samples)

pdf("Alpha-div_Samples.pdf", 
    width = 8, 
    height = 10)
Deg<-D + geom_boxplot(data=D$data , aes(x=Samples, color=NULL))+ 
  #stat_compare_means(
  #comparisons = comp,
  #label = "p.format",
  #tip.length = 0.05,
  #method = "wilcox.test",
  #p.adjust.methods="holm")+ 
  theme(axis.text=element_text(size=14,angle = 0, hjust = 0.5, vjust=0.5),
        axis.title = element_text(size=14, face="bold"), plot.title =element_text(size=14, face="bold", 
                                                                                  margin =ggplot2:: margin(t = 0, r = 0, b = 20, l = 0)), strip.text.x = element_text(size=14, face="bold" )) +
  xlab("Samples")+theme(legend.text=element_text(size=14),
                     legend.title = element_text(size=14, face="bold"))+
  theme(legend.position="none")+theme(axis.text.x = element_text(angle=0, hjust = 0.5, vjust=0.5))
Deg
dev.off()
ggsave(filename = "Alpha-div_Samples.png", plot = Deg, device = "png",width = 8, height = 10)

#####Beta diversity#####

###Exploration of feature prevalence in the dataset###
#which we will define here as the number of samples in which a taxon appears at least once.
# Compute prevalence of each feature, store as data.frame#

##All dataset##
prevdf = apply(X = otu_table(Gooddata),
               MARGIN = ifelse(taxa_are_rows(Gooddata), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(Gooddata),
                    tax_table(Gooddata))


###Prevalence Filtering to make threshold of taxa###
##Taxa filtering##
##les fonctions ci-dessous ont éliminé bcp de taxa; c'est justement l'idée de cette approche. 
#The following uses five percent of all samples as the prevalence threshold.

# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(Gooddata)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps = prune_taxa(keepTaxa, Gooddata)
ps
Gooddata2<-ps
Gooddata2
print_ps(Gooddata2)

####Unconstrained Ordinations- Beta diversity####

###Samples###

# Ordinate using Principal Coordinate analysis
Gooddata_pcoa <- ordinate(
  physeq = Gooddata2, 
  method = "PCoA", 
  distance = "bray"
)

# Méthode manuelle - extraction des données
# Extraire les scores PCoA
pcoa_scores <- as.data.frame(Gooddata_pcoa$vectors[, 1:2])
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")

# Extraire les métadonnées
metadata <- data.frame(sample_data(Gooddata2))

# Combiner les données
plot_data <- cbind(pcoa_scores, Samples = metadata$Samples)

# Calculer les pourcentages de variance expliquée
variance_explained <- round(100 * Gooddata_pcoa$values$Eigenvalues[1:2] / sum(Gooddata_pcoa$values$Eigenvalues), 1)
col_vector <- c("B8" = "blue","B11" = "red","W8"="orange","W11"="green")
# Créer le graphique manuellement
pdf("PcoA_Samples.pdf", 
    width = 10, 
    height = 8)
Deg <- ggplot(plot_data, aes(x = PCoA1, y = PCoA2, color = as.factor(Samples))) +
  geom_point(size = 4) +
  scale_color_manual(
    values = col_vector,
    name = "Samples"
  ) +
  labs(
    x = paste0("PCoA1 (", variance_explained[1], "%)"),
    y = paste0("PCoA2 (", variance_explained[2], "%)"),
    title = "Bacterial community according to Samples",
    color = "Samples"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14, angle = 0, hjust = 0.5, vjust = 0.5), 
    axis.title = element_text(size = 14, face = "bold"), 
    plot.title = element_text(size = 14, face = "bold", 
                              margin = margin(t = 0, r = 0, b = 20, l = 0)),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold")
  ) 

print(Deg)
dev.off()

ggsave(filename = "PcoA.png", plot = Deg, device = "png", width = 10, height = 8)
# Optionnel: Test PERMANOVA pour significativité

bray_dist <- phyloseq::distance(Gooddata2, method = "bray")
adonis_result <- adonis2(bray_dist ~ Samples, 
                         data = metadata)

cat("Résultat PERMANOVA:\n")
print(adonis_result)

# Ajouter la valeur p au graphique
p_value <- adonis_result$`Pr(>F)`[1]
pdf("PcoA_Samples.pdf", 
    width = 10, 
    height = 8)
Deg <- Deg + 
  annotate("text", x = min(plot_data$PCoA1), y = max(plot_data$PCoA2),
           label = paste("PERMANOVA p =", round(p_value, 4)),
           hjust = 0, vjust = 1, size = 5, fontface = "bold")

print(Deg)
dev.off()

###AVEC PERMANOVA

# 1. Extraire les coordonnées PCoA (axes 1 et 2)
# ----------------------------------------------------------------------------
pcoa_points <- as.data.frame(Gooddata_pcoa$vectors[, 1:2])
colnames(pcoa_points) <- c("Axis.1", "Axis.2")

# Métadonnées
metadata <- data.frame(sample_data(Gooddata2))
metadata$Samples <- as.factor(metadata$Samples)

# Ajouter la variable Samples aux points
pcoa_points$Samples <- metadata$Samples

# ----------------------------------------------------------------------------
# 2. Extraire la matrice d'abondance et la table taxonomique
# ----------------------------------------------------------------------------
abund_matrix <- as.matrix(otu_table(Gooddata2))
taxa_table <- as.data.frame(tax_table(Gooddata2))

# Mettre les échantillons en lignes
if(taxa_are_rows(Gooddata2)) {
  abund_matrix <- t(abund_matrix)
}
# Maintenant : lignes = échantillons, colonnes = ASV

# Vérifier que les noms de colonnes correspondent aux noms de lignes de taxa_table
if(!all(colnames(abund_matrix) %in% rownames(taxa_table))) {
  stop("Les noms des ASV dans la matrice d'abondance ne correspondent pas à la table taxonomique.")
}

# ----------------------------------------------------------------------------
# 3. Projeter les taxons sur l'ordination avec envfit
# ----------------------------------------------------------------------------
set.seed(123)
fit_taxa <- envfit(Gooddata_pcoa$vectors[, 1:2] ~ ., 
                   data = as.data.frame(abund_matrix), 
                   permutations = 999)

# Extraire les coordonnées des flèches
taxa_scores <- as.data.frame(scores(fit_taxa, display = "vectors"))
colnames(taxa_scores) <- c("Axis.1", "Axis.2")
taxa_scores$ASV <- rownames(taxa_scores)
taxa_scores$pval <- fit_taxa$vectors$pvals
taxa_scores$r2 <- fit_taxa$vectors$r

# ----------------------------------------------------------------------------
# 4. Joindre la taxonomie pour obtenir un label lisible (Genre)
# ----------------------------------------------------------------------------
# Récupérer le genre pour chaque ASV
taxa_scores$Genus <- taxa_table[taxa_scores$ASV, "Genus"]

# Nettoyer : remplacer les NA ou chaînes vides par l'ASV (ou par le niveau supérieur)
taxa_scores$Label <- ifelse(is.na(taxa_scores$Genus) | taxa_scores$Genus == "",
                            taxa_scores$ASV,
                            taxa_scores$Genus)

# Éviter les doublons de labels (plusieurs ASV avec le même genre) – on garde celui avec le meilleur R²
taxa_scores <- taxa_scores[order(taxa_scores$r2, decreasing = TRUE), ]
taxa_scores <- taxa_scores[!duplicated(taxa_scores$Label), ]

# Filtrer les taxons significatifs (p < 0.05)
taxa_signif <- subset(taxa_scores, pval < 0.05)

# Si aucun significatif, prendre les 5 avec le meilleur R²
if(nrow(taxa_signif) == 0) {
  message("Aucun taxon significatif (p < 0.05). Affichage des 5 avec le R² le plus élevé.")
  taxa_signif <- head(taxa_scores[order(taxa_scores$r2, decreasing = TRUE), ], 5)
}

# ----------------------------------------------------------------------------
# 5. Variance expliquée par les axes
# ----------------------------------------------------------------------------
eigenvalues <- Gooddata_pcoa$values$Eigenvalues
total_variance <- sum(eigenvalues)
contribution_axis_1 <- eigenvalues[1] / total_variance * 100
contribution_axis_2 <- eigenvalues[2] / total_variance * 100

# ----------------------------------------------------------------------------
# 6. PERMANOVA (avec avertissement si peu de réplicats)
# ----------------------------------------------------------------------------
bray_dist <- phyloseq::distance(Gooddata2, method = "bray")

if(any(table(metadata$Samples) < 2)) {
  cat("ATTENTION : Certains niveaux de 'Samples' n'ont qu'un seul échantillon.\n",
      "Le test PERMANOVA peut donner un R² = 1 et une p-value non fiable.\n",
      "Envisagez de regrouper en condition (B vs W).\n")
}

adonis_result <- adonis2(bray_dist ~ Samples, data = metadata, permutations = 999)
p_value <- adonis_result$`Pr(>F)`[1]
r_squared <- adonis_result$R2[1]

if(is.na(p_value)) {
  p_label <- "p = NA"
  significance <- ""
} else {
  significance <- ifelse(p_value < 0.05, "*", "ns")
  p_label <- sprintf("p = %.3f", p_value)
}

# ----------------------------------------------------------------------------
# 7. Palette de couleurs (automatique ou manuelle)
# ----------------------------------------------------------------------------
samples_levels <- levels(metadata$Samples)
col_vector <- setNames(hue_pal()(length(samples_levels)), samples_levels)
# Si vous préférez la palette d'origine :
# col_vector <- c("B8" = "blue","B11" = "red","W8"="orange","W11"="green")  # jaune -> orange

# ----------------------------------------------------------------------------
# 8. Graphique final avec labels de genre
# ----------------------------------------------------------------------------
A <- ggplot(pcoa_points, aes(x = Axis.1, y = Axis.2)) + 
  geom_point(aes(color = Samples), size = 4, alpha = 0.8) +
  scale_color_manual(values = col_vector, name = "Samples") +
  
  # Flèches des taxons
  geom_segment(data = taxa_signif,
               aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgray", alpha = 0.7, linewidth = 0.6) +
  
  # Étiquettes des taxons : on utilise maintenant le champ "Label"
  geom_text_repel(data = taxa_signif,
                  aes(x = Axis.1, y = Axis.2, label = Label),
                  size = 3.5, fontface = "italic",
                  box.padding = 0.5, point.padding = 0.3,
                  max.overlaps = 20, segment.color = "gray50") +
  
  # Annotation PERMANOVA
  annotate("text",
           x = min(pcoa_points$Axis.1, na.rm = TRUE),
           y = max(pcoa_points$Axis.2, na.rm = TRUE),
           label = sprintf("PERMANOVA\n%s\nR² = %.3f %s", p_label, r_squared, significance),
           hjust = 0, vjust = 1, size = 4.5, fontface = "bold", color = "darkred") +
  
  xlab(paste("Axis 1 (", round(contribution_axis_1, 2), "%)", sep = "")) +
  ylab(paste("Axis 2 (", round(contribution_axis_2, 2), "%)", sep = "")) +
  ggtitle("Eukaryotic community according to Samples") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", alpha = 0.5) +
  
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20)),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  )

print(A)

# Sauvegarder en PDF
pdf("PcoA_Samples_driving_PERMANOVA.pdf", width = 12, height = 14)
print(A)
dev.off()

# Sauvegarder les résultats du PERMANOVA dans un fichier texte
sink("PERMANOVA_results_Degradation.txt")
cat("=== RÉSULTATS PERMANOVA ===\n\n")
cat("Date de l'analyse:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Variable explicative: Location\n\n")
print(adonis_result)
cat("\n")
cat(sprintf("Valeur p: %.4f\n", p_value))
cat(sprintf("R²: %.3f (%.1f%% de variance expliquée)\n", r_squared, r_squared * 100))
cat(sprintf("Significativité: %s\n", significance))
cat("\n\n")
cat("Interprétation:\n")
if(p_value < 0.05) {
  cat("Il existe une différence significative dans la composition des communautés bactériennes\n")
  cat("entre les groupes de patients (1 vs 0).\n")
  cat(sprintf("Le statut tumoral explique %.1f%% de la variance observée.\n", r_squared * 100))
} else {
  cat("Aucune différence significative n'a été détectée dans la composition des communautés\n")
  cat("bactériennes entre les groupes de patients.\n")
}
sink()

cat("Les résultats du PERMANOVA ont été sauvegardés dans 'PERMANOVA_results.txt'\n")

#####Third part of analysis: Heatmap

library(BiocManager)
library(DESeq2)

###Genus###
###Samples###

###top N taxon##

# Étape 1 : Facteur pour la variable "Samples"
sample_data(Gooddata2)$Samples <- as.factor(sample_data(Gooddata2)$Samples)

# Étape 2 : Agrégation au niveau du genre (Genus) et suppression explicite des NA
ps.taxa <- tax_glom(Gooddata2, taxrank = "Genus", NArm = FALSE)

# Supprimer les NA dans la colonne "Genus" de la table taxonomique
ps.taxa <- prune_taxa(!is.na(as.vector(tax_table(ps.taxa)[, "Genus"])), ps.taxa)

# Vérifiez le nombre de taxons avant et après le filtrage
print(ntaxa(ps.taxa))  # Nombre de taxons restants

# Étape 3 : Sous-échantillonnage des données selon les valeurs de "Samples"
ps.taxa.sub <- subset_samples(ps.taxa, Samples %in% c("B8", "B11","W8","W11"))

# Étape 4 : Filtrage des taxons rares (> 90 % de zéros)
ps.taxa.pse.sub <- prune_taxa(
  rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9,
  ps.taxa.sub
)
# Vérifiez le nombre de taxons après le filtrage
print(ntaxa(ps.taxa.pse.sub))

# Transformation en matrice
otu_mat <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Vérifier orientation
if (taxa_are_rows(ps.taxa.pse.sub)) {
  otu_mat <- t(otu_mat)
}

# Ajouter un pseudocount pour éviter log(0)
otu_mat <- otu_mat + 1

# Liste des comparaisons
comparisons <- list(
  c("B8", "B11"),
  c("W8", "W11"), 
  c("B8", "W8"),
  c("B8", "W11"),
  c("B11", "W11"), 
  c("B11", "W8")
)

top_n <- 10

# Boucle
for (comp in comparisons) {
  
  group1 <- comp[1]
  group2 <- comp[2]
  comp_name <- paste0(group1, "_vs_", group2)
  
  cat("\n=== Analyse simplifiée pour", comp_name, "===\n")
  
  # Extraire les échantillons
  sample_ids <- sample_data(ps.taxa.pse.sub)$Samples
  
  s1 <- rownames(sample_data(ps.taxa.pse.sub))[sample_ids == group1]
  s2 <- rownames(sample_data(ps.taxa.pse.sub))[sample_ids == group2]
  
  # Moyennes (ici = valeur unique)
  mean1 <- colMeans(otu_mat[s1, , drop = FALSE])
  mean2 <- colMeans(otu_mat[s2, , drop = FALSE])
  
  # Log2 Fold Change
  log2FC <- log2(mean1 / mean2)
  
  results_df <- data.frame(
    ASV = names(log2FC),
    log2FoldChange = log2FC
  )
  
  # Ajouter taxonomie
  taxa_info <- as.data.frame(tax_table(ps.taxa.pse.sub))
  taxa_info$ASV <- rownames(taxa_info)
  
  results_df <- merge(results_df, taxa_info, by = "ASV")
  
  # Trier par effet absolu
  results_df <- results_df[order(abs(results_df$log2FoldChange), decreasing = TRUE), ]
  
  # Top taxons
  barplot_data <- results_df[1:min(top_n, nrow(results_df)), ]
  
  print(barplot_data[, c("Genus", "log2FoldChange")])
  
  # Sauvegarde
  write.csv(results_df, paste0("RESULTS_noStats_", comp_name, ".csv"), row.names = FALSE)
  
  # BARPLOT
  bar_plot <- ggplot(barplot_data, 
                     aes(x = reorder(Genus, log2FoldChange), 
                         y = log2FoldChange)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    labs(title = paste("Top", top_n, "-", comp_name),
         x = "Genus",
         y = "log2 Fold Change")
  
  pdf(paste0("Barplot_", comp_name, "_noStats.pdf"), width = 10, height = 6)
  print(bar_plot)
  dev.off()
  png(paste0("Barplot_", comp_name, ".png"), width = 12, height = 8, units = "in", res = 300)
  print(bar_plot)
  dev.off()
}
###Regroupement des échantillons 2 par 2 pour pouvoir utilisé DESeq2
##W vs B

library(DESeq2)

# Regroupement B vs W
sample_data(Gooddata2)$Condition <- ifelse(grepl("^B", sample_data(Gooddata2)$Samples), "B", "W")
sample_data(Gooddata2)$Condition <- as.factor(sample_data(Gooddata2)$Condition)

# Agrégation au niveau du genre
ps.taxa <- tax_glom(Gooddata2, taxrank = "Genus", NArm = FALSE)
ps.taxa <- prune_taxa(!is.na(as.vector(tax_table(ps.taxa)[, "Genus"])), ps.taxa)

# Filtrage des taxons rares (présents dans >=2 échantillons)
ps.taxa.filt <- prune_taxa(rowSums(otu_table(ps.taxa) > 0) >= 2, ps.taxa)

# DESeq2
ps_ds <- phyloseq_to_deseq2(ps.taxa.filt, ~ Condition)
ps_ds <- estimateSizeFactors(ps_ds, type = "poscounts")
ds <- DESeq(ps_ds, test = "Wald", fitType = "parametric", minReplicatesForReplace = Inf)
res <- results(ds, contrast = c("Condition", "W", "B"), alpha = 0.05)
res_ordered <- res[order(res$padj, na.last = NA), ]

# Tableau des résultats avec noms de genre
res_df <- as.data.frame(res_ordered)
res_df$Feature <- rownames(res_df)
tax_info <- as.data.frame(tax_table(ps.taxa.filt))
res_df$Genus <- tax_info[res_df$Feature, "Genus"]
res_df$Genus <- ifelse(is.na(res_df$Genus) | res_df$Genus == "", res_df$Feature, res_df$Genus)
write.csv(res_df, "DESeq2_B_vs_W_all_results.csv", row.names = FALSE)

# Top 10
top_n <- 10
top_features <- rownames(res_ordered[1:min(top_n, nrow(res_ordered)), ])

# Volcano plot
volcano_data <- res_df
volcano_data$Significatif <- ifelse(volcano_data$padj < 0.05 & !is.na(volcano_data$padj), "Oui", "Non")
volcano_data$Label <- ifelse(volcano_data$Feature %in% top_features, volcano_data$Genus, "")

p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = Significatif, label = Label)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Non" = "gray70", "Oui" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(size = 3.5, max.overlaps = 20, box.padding = 0.5, fontface = "italic") +
  labs(title = "Top 10 differentially abundant genera : B (n=2) vs W (n=2)", x = "Log2 Fold Change (W vs B)", y = "-Log10(padj)") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave("Volcano_B_vs_W.pdf", p_volcano, width = 8, height = 6)
print(p_volcano)

# Barplot sans légende
bar_data <- res_df[res_df$Feature %in% top_features, ]
bar_data <- bar_data[order(bar_data$log2FoldChange, decreasing = TRUE), ]

p_bar <- ggplot(bar_data, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange,
                              fill = ifelse(log2FoldChange > 0, "Enrichi dans W", "Enrichi dans B"))) +
  geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c("Enrichi dans W" = "steelblue", "Enrichi dans B" = "salmon")) +
  coord_flip() +
  labs(title = paste("Top", top_n, "Top 10 differentially abundant genera (B vs W)"),
       x = "Genre", y = "Log2 Fold Change (W vs B)") +
  theme_minimal() + theme(axis.text.y = element_text(face = "italic", size = 10))
ggsave("Barplot_B_vs_W.pdf", p_bar, width = 8, height = 6)
print(p_bar)

png(paste0("Barplot_B_vs_W.png", comp_name, ".png"), width = 12, height = 8, units = "in", res = 300)
print(p_bar)
dev.off()

# Optionnel : heatmap des comptes normalisés pour ces taxons
library(pheatmap)
# Extraire les comptes normalisés
# Extraire les comptes normalisés
norm_counts <- counts(ds, normalized = TRUE)

# Sélectionner les top taxons
top_counts <- norm_counts[top_features, ]

# Ajouter les noms de genres
tax_info <- as.data.frame(tax_table(ps.taxa.filt))
rownames(top_counts) <- tax_info[top_features, "Genus"]

# Annotation des colonnes
col_annot <- data.frame(Condition = sample_data(ps.taxa.filt)$Condition)
rownames(col_annot) <- colnames(top_counts)

# Heatmap
pheatmap(top_counts,
         scale = "row",
         annotation_col = col_annot,
         main = "Top 10 differentially abundant genera",
         fontsize_row = 10,
         angle_col = 0,
         filename = "Heatmap_B_vs_W.pdf")


####Forth part of analysis: statistical analysis and variable contributor for dissimilarity####

### Statistic analysis ###
##While PcoA gives us a visual of beta-diversity, 
#it does not test for statistical differences. 
#We do this with permutational analysis of variance (PERMANOVA) or analysis of similarity (ANOSIM). 
#These test whether the overall microbial community differs by my variable 
#of interest.
#We use ADONIS (PERMANOVA test) when we have continuous variables
# we perform permutational ANOVA test to see if any of the available information
#in the dataset is indicative of community structure.

#Recalculating distance and running an adonis
set.seed(6)
Gooddata2

#Gooddata2<-subset_samples(Gooddata2, Well!="",NArm=TRUE)
#Gooddata2<-subset_samples(Gooddata2, Samples!="",NArm=TRUE)

project_bray <- phyloseq::distance(Gooddata2, method= "bray",NArm = TRUE)
sample_df <- data.frame(sample_data(Gooddata2))
colnames(sample_df)
#now the adonis test to see if there is a signficant difference according to parasite species variables
res.adonis <- adonis2(project_bray ~ Samples+Well+Information+CYTO, data=sample_df, method="bray",by = "term")

res.adonis 
results<-data.frame(as.data.frame(res.adonis))
write.csv(results, "PERMANOVATEST_allvariable2.csv")

results$variable <- row.names(res.adonis)
results$Contribution <- results$R2 * 100

# ajouter colonne significativité
results$Signif <- ifelse(results$Pr..F. < 0.05, "Significatif", "Non significatif")

pdf("contribution_all_2.pdf", width = 10, height = 8)

Per<-ggplot(results, aes(x = reorder(variable, Contribution), 
                    y = Contribution, 
                    fill = Signif)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(Pr..F., 3)), 
            hjust = -0.1, size = 4) +
  coord_flip() +
  scale_fill_manual(values = c("Significatif" = "steelblue",
                               "Non significatif" = "grey")) +
  theme_bw() +
  xlab("") +
  ylab("% Contribution") +
  ggtitle("Contribution des variables (PERMANOVA)") +
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"))
Per
dev.off()


