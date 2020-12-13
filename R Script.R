#MICB 447 Team 9 Supplementary R Script

# ------- ALPHA DIVERSITY PLOTS FOR DIETARY PROTEIN ---------

setwd("/Users/KristiMacBookAir2015/Desktop/MICB447/alpha")

library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggpubr)

# ----- Setting up -----

# Shannon's diversity function
shannons = function(x){
  present = x[x>0]
  p = present/sum(present)
  -sum(p*log(p))
}

# Pielou's evennenss function
evenness = function(x){
  present = x[x>0]
  p = present/sum(present)
  -sum(p*log(p))/log(sum(x>0))
}

# richness function
richness = function(x){
  return(sum(x>0))
}

# Load data
biom = import_biom("table-with-taxonomy3.biom")
taxa_table = otu_table(biom)
taxonomy = tax_table(biom)
metadata = read.table("dog_metadata.tsv",sep="\t",header=T,row.names = 1)

microbial_samples = colnames(taxa_table)
metadata_samples = rownames(metadata)
which_metadata = c()
for (i in 1:dim(taxa_table)[2]){
  which_metadata = c(which_metadata,which(metadata_samples ==
                                            microbial_samples[i]))
}
metadata = metadata[which_metadata,]

# Calculate alpha diversity
metadata$richness = apply(taxa_table,2,richness)
metadata$shannons = apply(taxa_table,2,shannons)
metadata$evenness = apply(taxa_table,2,evenness)


# ----- PROTEIN SOURCE -----

# Filter to remove "other" sample and blank cells
shannon_filtered_data <- metadata %>% filter(protein_source != "") %>% 
  filter(protein_source != "other")


### Discrete variable, Shannon
shannon_protein_source_plot <- ggplot(shannon_filtered_data,aes(x=protein_source,y=shannons)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitter(width=0.1)) +
  xlab("Protein Source") +
  ylab("Shannon Diversity Index")

shannon_protein_source_plot


# Kruskal-Wallis test
kruskal.test(shannons ~ protein_source, data =
               filter(shannon_filtered_data, protein_source!=""))


# ----- % CRUDE PROTEIN -----

# Filter to remove blank cells
crude_filtered_data <- metadata %>% filter(perc_crude_protein_min_group != "")

### Discrete variable, Shannon
shannon_perc_crude_plot <- ggplot(crude_filtered_data,aes(x=perc_crude_protein_min_group,y=shannons)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitter(width=0.1)) +
  scale_x_discrete(labels = c('14-20%', '21-25%', '26-35%')) +
  xlab("Percent Crude Protein") +
  ylab("Shannon Diversity Index")

shannon_perc_crude_plot


# Kruskal-Wallis test, with blank cells removed
kruskal.test(shannons ~ perc_crude_protein_min_group, data =
               filter(metadata,perc_crude_protein_min_group!=""))



# ---------- BETA DIVERSITY PLOTS FOR DIETARY PROTEIN --------

setwd("/Users/KristiMacBookAir2015/Desktop/MICB447/beta")

library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggpubr)


# Helper functions
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# Importing all data files
biom_file <- import_biom("table-with-taxonomy2.biom")
metadata <- import_qiime_sample_data("dog_metadata.tsv")
tree <- read_tree_greengenes("tree.nwk")


# Combine all objects into phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)
# Overview of phyloseq object
physeq

# Use sample_data() to look at metadata of phyloseq object
# Look only at first 6 lines of table with head()
head(sample_data(physeq))

# Set set of random numbers
set.seed(711)

# taxonomic rank names
rank_names(physeq)
# Rename column names for taxonomic ranks
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class",
                                 "Order", "Family", "Genus",
                                 "Species")
rank_names(physeq)

# -------------------- PROTEIN SOURCE!!!!! ------------------

# BETA DIVERSITY
# Remove "other" and empty rows
# Filter out data based on metadata category with !=
no_other <- subset_samples(physeq, protein_source != "other", protein_source != "") 
no_other

no_other_empty <- subset_samples(no_other, protein_source != "")
no_other_empty

# Exclude samples with less than 15,000 reads
no_other_empty_15000 <- prune_samples(sample_sums(no_other_empty) >= 15000, no_other_empty)
no_other_empty_15000

# Beta diversity PCoA plot
# Diversity requires rarefied taxa tables
# We know from Qiime that we should rarefy to 1250 so right palm
# samples level out
physeq_rar <- rarefy_even_depth(no_other_empty_15000, sample.size = 4000)

# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

# We define the type of analysis with ordinate() and setting the
# method = PCoA, setting distance to type of beta diversity analysis
# you can change analysis to jaccard, bray-curtis and so on


# ----- WEIGHTED UNIFRAC, protein source -----
ord_wunifrac_prot_source <- ordinate(physeq_rar_RA, method = "PCoA",
                distance = "wunifrac")

# Plot data
plot_ordination(physeq_rar_RA,
                ord_wunifrac_prot_source,
                type = "sample",
                color = "protein_source") +
  # Adding text to plot
  annotate(geom = "text",
           label = "PCoA plot of weighted UniFrac by protein source",
           size = 4) + 
  # Manually adjust colours for points
  scale_colour_manual(values = c("blue", "red",
                                 "green", "pink"),
                      labels = c("Chicken", "Hydrolyzed",
                                 "Lamb",
                                 "Fish")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("Protein source")) + 
  theme_bw(base_size = 14)


# ----- BRAY CURTIS, protein source -----

ord_bray_prot_source <- ordinate(physeq_rar_RA, method = "PCoA",
                                     distance = "bray")

# Plot data
plot_ordination(physeq_rar_RA,
                ord_bray_prot_source,
                type = "sample",
                color = "protein_source",
                title = "PCoA of Bray-Curtis distances by protein source") +
  # Adding text to plot
  annotate(geom = "text",
           label = "PCoA plot of weighted UniFrac by protein source",
           size = 4) + 
  # Manually adjust colours for points
  scale_colour_manual(values = c("blue", "red",
                                 "green", "pink"),
                      labels = c("Chicken", "Hydrolyzed",
                                 "Lamb",
                                 "Fish")) +
  guides(colour = guide_legend("Protein source")) + 
  theme_bw(base_size = 14)



# -------------------- % CRUDE PROTEIN!!!!! ------------------

# Remove "other" and empty rows
# Filter out data based on metadata category with !=
no_empty <- subset_samples(physeq, perc_crude_protein_min_group != "") 
no_empty

# Exclude samples with less than 15,000 reads
no_empty_15000 <- prune_samples(sample_sums(no_empty) >= 15000, no_empty)
no_empty_15000

# Beta diversity PCoA plot
# Diversity requires rarefied taxa tables
# We know from Qiime that we should rarefy to 1250 so right palm
# samples level out
physeq_rar <- rarefy_even_depth(no_empty_15000, sample.size = 4000)

# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

# We define the type of analysis with ordinate() and setting the
# method = PCoA, setting distance to type of beta diversity analysis
# you can change analysis to jaccard, bray-curtis and so on


# ----- WEIGHTED UNIFRAC, % crude protein -----
ord_wunifrac_crude <- ordinate(physeq_rar_RA, method = "PCoA",
                                     distance = "wunifrac")

# Plot data
plot_ordination(physeq_rar_RA,
                ord_wunifrac_crude,
                type = "sample",
                color = "perc_crude_protein_min_group") +
  # Adding text to plot
  annotate(geom = "text",
           label = "PCoA plot of weighted UniFrac by % crude protein",
           size = 4) + 
  # Manually adjust colours for points
  scale_colour_manual(values = c("deeppink", "darkseagreen",
                                 "darkmagenta"),
                      labels = c("14-20%", "21-26%",
                                 "26-35%")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("Percent crude protein")) + 
  theme_bw(base_size = 14)


# ----- BRAY CURTIS, protein source -----

ord_bray_crude <- ordinate(physeq_rar_RA, method = "PCoA",
                                 distance = "bray")

# Plot data
plot_ordination(physeq_rar_RA,
                ord_bray_crude,
                type = "sample",
                color = "perc_crude_protein_min_group",
                title = "PCoA of Bray-Curtis distances by % crude protein") +
  # Adding text to plot
  annotate(geom = "text",
           label = "PCoA plot of Bray-Curtis by % crude protein",
           size = 4) + 
  # Manually adjust colours for points
  scale_colour_manual(values = c("deeppink", "darkseagreen",
                                 "darkmagenta"),
                      labels = c("14-20%", "21-26%",
                                 "26-35%")) +
  stat_ellipse(type = "norm", size = 1) + 
  guides(colour = guide_legend("Percent crude protein")) + 
  theme_bw(base_size = 14)


