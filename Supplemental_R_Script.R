#MICB 447 Team 9 Supplemental R Script 
#Joshua Calalang, Honor Cheung, Kristi Lichimo, Bonny So
#16 December 2020

#--------ALPHA AND BETA DIVERSITY FOR IBD Vs HEALTHY SAMPLES---------
library(ggplot2)

### Alpha diversity functions
# Shannon's diversity
shannons = function(x){
  present = x[x>0]
  p = present/sum(present)
  -sum(p*log(p))
}

### Load data
biom = import_biom("table-with-taxonomy.biom")
taxa_table = otu_table(biom)
taxonomy = tax_table(biom)
metadata = read.table("dogs_metadata_at_risk.tsv",sep="\t",header=T,row.names = 1)

#select only samples with metadata
microbial_samples = colnames(taxa_table)
metadata_samples = rownames(metadata)
which_metadata = c()
for (i in 1:dim(taxa_table)[2]){
  which_metadata = c(which_metadata,which(metadata_samples ==
                                            microbial_samples[i]))
}
metadata = metadata[which_metadata,]


### Calculate alpha diversity
metadata$shannons = apply(taxa_table,2,shannons)

#plot data FIGURE 1A
### Discrete variable
graph_disease_stat<-ggplot(metadata,aes(x=disease_stat,y=shannons)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitter(width=0.1))
graph_disease_stat + labs(x = "disease stat")

#Stats FIGURE 1A
wilcox.test(shannons ~ disease_stat, data=metadata)

#Beta diversity
# Loading packages
library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)

# Helper functions
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
# Combine all objects into phyloseq object
physeq <- merge_phyloseq(biom, metadata, tree)
# Overview of phyloseq object
physeq

# taxonomic rank names
rank_names(physeq)
# Rename column names for taxonomic ranks
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class",
                                 "Order", "Family", "Genus",
                                 "Species")
rank_names(physeq)

# Beta diversity PCoA plot
# Diversity requires rarefied taxa tables
physeq_rar <- rarefy_even_depth(physeq, sample.size = 15000)
  

# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

# We define the type of analysis with ordinate() and setting the
# method = PCoA, setting distance to type of beta diversity analysis
# you can change analysis to jaccard, bray-curtis and so on
ord <- ordinate(physeq_rar_RA, method = "PCoA",
                distance = "wunifrac")

# Plot data for IBD healthy FIGURE 1B
plot_ordination(physeq_rar_RA,
                ord,
                type = "sample",
                color = "disease_stat",
                title = "PCoA (Weighted Unifrac)") +
  
  # Manually adjust colours for points
  scale_colour_manual(values = c("#00AFBB", "#E7B800"),
                      labels = c("healthy", "IBD")) 



#---ALPHA DIVERSITY ANALYSIS FOR AT-RISK VS NON-AT-RISK BREED GROUPS--------------

#plot data FIGURE 2
#2 Way ANOVA  
#Visualise data  
graph1<-ggboxplot(metadata, x = "at_risk", y = "shannons", color = "disease_stat",
            add = c("mean_se", "dotplot"),
            palette = c("#00AFBB", "#E7B800")) 
graph1 + labs(x = "at risk")

#Compute 2-Way ANOVA for unbalanced designs FIGURE 2
library(car)
my_anova <- aov(shannons ~ at_risk * disease_stat, data = metadata)
Anova(my_anova, type = "III")
TukeyHSD(my_anova)
                                        
 #--------DIFFERENTIAL ABUNDANCE ANALYSIS FOR AT-RISK VS NON-AT-RISK BREED GROUPS----------
                                        
# Keep only abundant ASVs
# First determine counts of ASV across all samples
total_counts <- taxa_sums(physeq)

# Calculate relative abundance for each ASV
relative_abundance <- calculate_RA(total_counts)

# Determine which ASVs are more abundant than 0.1%
# Change this if you want a different cutoff (0.001 = 0.1%)
abundant <- relative_abundance > 0.001 
abundant_taxa <- prune_taxa(abundant, physeq)
abundant_taxa_filtered<-prune_samples(sample_sums(abundant_taxa)>=5015,abundant_taxa)
#Throws out all false values and keeps the ones that are true
# Check the resulting new phyloseq object with much fewer taxa
abundant_taxa_filtered
#Double check filtered features 
sort(sample_sums(abundant_taxa_filtered),decreasing = T)
# now your phyloseq object is called "Family"
family <- tax_glom(abundant_taxa_filtered, taxrank = "Family", NArm = FALSE)
family

# Only keep nonrisk samples
nonrisk<- subset_samples(family, at_risk == "no")
nonrisk

# Let's limit to mixed breed samples that we filtered above
deseq_feature <- phyloseq_to_deseq2(nonrisk, ~ disease_stat)
# now the rest of the deseq2 code:
geo_means <- apply(counts(deseq_feature), 1, gm_mean)
deseq_feature <- estimateSizeFactors(deseq_feature, geoMeans = geo_means)
deseq_feature <- DESeq(deseq_feature, fitType="local")

feature_diff_abund <- results(deseq_feature)
# Define your cutoff for the adjusted p-value
alpha <- 0.05
# Reformat information as data frame including feature as variable
significant_feature <- as_tibble(feature_diff_abund, rownames = "feature")
# Keep only significant results and sort by adjusted p-value
significant_feature <- filter(significant_feature, padj < alpha) #want only rows that are below alpha
significant_feature <- arrange(significant_feature, padj)#arrange by smallest to largest p-value; sorted by padj

# Get the taxonomic information as a data frame
taxa_df <- as_tibble(as.data.frame(tax_table(physeq)), rownames = "feature")

# Combine the significant features with taxonomic classification
significant_feature <- inner_join(significant_feature, taxa_df)
dim(significant_feature)

# Plot differential abundance
significant_feature %>%
  ggplot(aes(x = log2FoldChange, y = Family)) +
  geom_col() 


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
biom_file <- import_biom("table-with-taxonomy3.biom")
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

# ---------------------- PROTEIN SOURCE -----------------------

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



# ------------------------- % CRUDE PROTEIN --------------------------

# Remove "other" and empty rows
# Filter out data based on metadata category with !=
no_empty <- subset_samples(physeq, perc_crude_protein_min_group != "") 
no_empty

# Exclude samples with less than 15,000 reads
no_empty_15000 <- prune_samples(sample_sums(no_empty) >= 15000, no_empty)
no_empty_15000

# Beta diversity PCoA plot
# Diversity requires rarefied taxa tables

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
                                         
                                         
 
# ------- ALPHA DIVERSITY PLOTS FOR NEUTER STATUS ---------

# Load CRAN packages
library(tidyverse)
library(vegan)

# Load Bioconductor packages
library(phyloseq)
library(DESeq2)
library(ggplot2)

### Alpha diversity functions
# Shannon's diversity
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


# ------- INTACT + HEALTHY ---------                                         
                                         
### Load data
biom = import_biom("/Users/joshuacalalang/Desktop/MICB447_R/Project2/intact_healthy/table-with-taxonomy.biom")
taxa_table = otu_table(biom)
taxonomy = tax_table(biom)
metadata = read.table("/Users/joshuacalalang/Desktop/MICB447_R/Project2/dog_metadata.tsv",sep="\t",header=T,row.names = 1)

microbial_samples = colnames(taxa_table)
metadata_samples = rownames(metadata)
which_metadata = c()
for (i in 1:dim(taxa_table)[2]){
  which_metadata = c(which_metadata,which(metadata_samples ==
                                            microbial_samples[i]))
}
metadata = metadata[which_metadata,]

### Calculate alpha diversity
metadata$richness = apply(taxa_table,2,richness)
metadata$shannons = apply(taxa_table,2,shannons)
metadata$evenness = apply(taxa_table,2,evenness)

### Discrete variable
ggplot(metadata,aes(x=sex,y=shannons)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitter(width=0.1)) +
  xlab("Sex") +
  ylab("Shannon Diversity")

kruskal.test(shannons ~ sex, data = filter(metadata, sex!=""))
                                         
                                         
# ------- INTACT + IBD ---------
                                         
### Load data
biom = import_biom("/Users/joshuacalalang/Desktop/MICB447_R/Project2/intact_IBD/table-with-taxonomy.biom")
taxa_table = otu_table(biom)
taxonomy = tax_table(biom)
metadata = read.table("/Users/joshuacalalang/Desktop/MICB447_R/Project2/dog_metadata.tsv",sep="\t",header=T,row.names = 1)

microbial_samples = colnames(taxa_table)
metadata_samples = rownames(metadata)
which_metadata = c()
for (i in 1:dim(taxa_table)[2]){
  which_metadata = c(which_metadata,which(metadata_samples ==
                                            microbial_samples[i]))
}
metadata = metadata[which_metadata,]

### Calculate alpha diversity
metadata$richness = apply(taxa_table,2,richness)
metadata$shannons = apply(taxa_table,2,shannons)
metadata$evenness = apply(taxa_table,2,evenness)

### Discrete variable
ggplot(metadata,aes(x=sex,y=shannons)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitter(width=0.1)) +
  xlab("Sex") +
  ylab("Shannon Diversity")

kruskal.test(shannons ~ sex, data = filter(metadata, sex!=""))                                         
                                         
                                    
# ------- NEUTERED + HEALTHY ---------   
                                         
### Load data
biom = import_biom("/Users/joshuacalalang/Desktop/MICB447_R/Project2/N_healthy/table-with-taxonomy.biom")
taxa_table = otu_table(biom)
taxonomy = tax_table(biom)
metadata = read.table("/Users/joshuacalalang/Desktop/MICB447_R/Project2/dog_metadata.tsv",sep="\t",header=T,row.names = 1)

microbial_samples = colnames(taxa_table)
metadata_samples = rownames(metadata)
which_metadata = c()
for (i in 1:dim(taxa_table)[2]){
  which_metadata = c(which_metadata,which(metadata_samples ==
                                            microbial_samples[i]))
}
metadata = metadata[which_metadata,]

### Calculate alpha diversity
metadata$richness = apply(taxa_table,2,richness)
metadata$shannons = apply(taxa_table,2,shannons)
metadata$evenness = apply(taxa_table,2,evenness)

### Discrete variable
ggplot(metadata,aes(x=sex,y=shannons)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitter(width=0.1)) +
  xlab("Sex") +
  ylab("Shannon Diversity")

kruskal.test(shannons ~ sex, data = filter(metadata, sex!="")) 
                                         

# ------- NEUTERED + IBD ---------  
                                         
### Load data
biom = import_biom("/Users/joshuacalalang/Desktop/MICB447_R/Project2/N_IBD/table-with-taxonomy.biom")
taxa_table = otu_table(biom)
taxonomy = tax_table(biom)
metadata = read.table("/Users/joshuacalalang/Desktop/MICB447_R/Project2/dog_metadata.tsv",sep="\t",header=T,row.names = 1)

microbial_samples = colnames(taxa_table)
metadata_samples = rownames(metadata)
which_metadata = c()
for (i in 1:dim(taxa_table)[2]){
  which_metadata = c(which_metadata,which(metadata_samples ==
                                            microbial_samples[i]))
}
metadata = metadata[which_metadata,]

### Calculate alpha diversity
metadata$richness = apply(taxa_table,2,richness)
metadata$shannons = apply(taxa_table,2,shannons)
metadata$evenness = apply(taxa_table,2,evenness)

### Discrete variable
ggplot(metadata,aes(x=sex,y=shannons)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitter(width=0.1)) +
  xlab("Sex") +
  ylab("Shannon Diversity")

kruskal.test(shannons ~ sex, data = filter(metadata, sex!=""))
                                         
 
# ------- SUPPLEMENTAL BETA DIVERSITY PLOTS FOR NEUTER STATUS ---------                                         
                                         
# Loading packages
library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)

# Helper functions
gm_mean <- function(x, na.rm = TRUE) {  
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x)) }

# ------- INTACT + HEALTHY ---------                                         
                                         
# Importing all data files
biome_file <- import_biom("/Users/joshuacalalang/Desktop/MICB447_R/Project2/intact_healthy/table-with-taxonomy.biom") # recall: need to use quotation marks to specify file path
metadata <- import_qiime_sample_data("/Users/joshuacalalang/Desktop/MICB447_R/Project2/dog_metadata.tsv")
tree <- read_tree_greengenes("/Users/joshuacalalang/Desktop/MICB447_R/Project2/intact_healthy/tree.nwk")

# Combining all objects into phyloseq object
physeq <- merge_phyloseq(biome_file, metadata, tree)

# Overview of phyloseq object
physeq 

# head function allows us to only look at first 6 lines
head(sample_data(physeq)) # sample_data() allows us to look at the metadata of phyloseq object

# Set set of random numbers (to make downstream analysis reproducible)
set.seed(711)

# taxonomic rank names
rank_names(physeq)

# Shows us first 6 lines [optional]
head(tax_table(physeq))

# Renaming the column names 
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class",
                                 "Order", "Family", "Genus", "Species")

# Beta diversity PCoA plot
# Diversity analysis requires rarefied datataxa tables
# rarefaction:
physeq_rar <- rarefy_even_depth(physeq, sample.size = 3000) 

# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

# We define the type of analysis using " ordinate() " and setting the 
  # method = PCoA, and also setting the distance to type of beta diversity analysis
# you can change analysis to jaccard, bray-curtis, and so on
ord <- ordinate(physeq_rar_RA, method = "PCoA", distance = "wunifrac")

# Plot data
plot_ordination(physeq_rar_RA, 
                ord, 
                type = "samples", 
                color = "sex", 
                title = "PcoA (Weighted Unifrac)") +
  
  # Adding text to plot
  annotate(geom = "text", 
           label = "",
           x = -0.08,
           y = 0.18,
           size = 4) + 

# Adding elipses
stat_ellipse(type = "norm", size = 0.5) + # the size value specifies thickness of the ellipses line
  
guides(colour = guide_legend("Sex")) + # changing legend title
  
  theme_bw(base_size = 14) # changing background colour

                                         
# ------- INTACT + IBD ---------                                           
                                         
# Importing all data files
biome_file <- import_biom("/Users/joshuacalalang/Desktop/MICB447_R/Project2/intact_IBD/table-with-taxonomy.biom") # recall: need to use quotation marks to specify file path
metadata <- import_qiime_sample_data("/Users/joshuacalalang/Desktop/MICB447_R/Project2/dog_metadata.tsv")
tree <- read_tree_greengenes("/Users/joshuacalalang/Desktop/MICB447_R/Project2/intact_IBD/tree.nwk")

# Combining all objects into phyloseq object
physeq <- merge_phyloseq(biome_file, metadata, tree)

# Overview of phyloseq object
physeq 

# head function allows us to only look at first 6 lines
head(sample_data(physeq)) # sample_data() allows us to look at the metadata of phyloseq object

# Set set of random numbers (to make downstream analysis reproducible)
set.seed(711)

# taxonomic rank names
rank_names(physeq)

# Shows us first 6 lines [optional]
head(tax_table(physeq))

# Renaming the column names 
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class",
                                 "Order", "Family", "Genus", "Species")

# Filter data based on metadata category with " == " [useful!!!]
spn <- subset_samples(physeq, spayed_neutered ==  "yes")

# Format data as table with as.data.frame()
as.data.frame(sample_names(spn))

# Show read count for each sample using " sample_sums() "
as.data.frame(sample_sums(spn))

# Exclude samples with less than 15000 reads
spn_15000 <- prune_samples(sample_sums(spn) >= 15000, spn)
spn_15000

# Beta diversity PCoA plot
# Diversity analysis requires rarefied datataxa tables
# rarefaction:
physeq_rar <- rarefy_even_depth(physeq, sample.size = 3000) 

# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

# We define the type of analysis using " ordinate() " and setting the 
  # method = PCoA, and also setting the distance to type of beta diversity analysis
# you can change analysis to jaccard, bray-curtis, and so on
ord <- ordinate(physeq_rar_RA, method = "PCoA", distance = "wunifrac")

# Plot data
plot_ordination(physeq_rar_RA, 
                ord, 
                type = "samples", 
                color = "sex", 
                title = "PcoA (Weighted Unifrac)") +
  
  # Adding text to plot
  annotate(geom = "text", 
           label = "",
           x = -0.08,
           y = 0.18,
           size = 4) + 

# Adding elipses
stat_ellipse(type = "norm", size = 0.5) + # the size value specifies thickness of the ellipses line
  
guides(colour = guide_legend("Sex")) + # changing legend title
  
  theme_bw(base_size = 14) # changing background colour

                                         
# ------- NEUTERED + HEALTHY ---------   
# Importing all data files
biome_file <- import_biom("/Users/joshuacalalang/Desktop/MICB447_R/Project2/N_healthy/table-with-taxonomy.biom") # recall: need to use quotation marks to specify file path
metadata <- import_qiime_sample_data("/Users/joshuacalalang/Desktop/MICB447_R/Project2/dog_metadata.tsv")
tree <- read_tree_greengenes("/Users/joshuacalalang/Desktop/MICB447_R/Project2/N_healthy/tree.nwk")

# Combining all objects into phyloseq object
physeq <- merge_phyloseq(biome_file, metadata, tree)

# Overview of phyloseq object
physeq 

# head function allows us to only look at first 6 lines
head(sample_data(physeq)) # sample_data() allows us to look at the metadata of phyloseq object

# Set set of random numbers (to make downstream analysis reproducible)
set.seed(711)

# taxonomic rank names
rank_names(physeq)

# Shows us first 6 lines [optional]
head(tax_table(physeq))

# Renaming the column names 
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class",
                                 "Order", "Family", "Genus", "Species")

# Filter data based on metadata category with " == " [useful!!!]
spn <- subset_samples(physeq, spayed_neutered ==  "yes")

# Format data as table with as.data.frame()
as.data.frame(sample_names(spn))

# Show read count for each sample using " sample_sums() "
as.data.frame(sample_sums(spn))

# Exclude samples with less than 15000 reads
spn_15000 <- prune_samples(sample_sums(spn) >= 15000, spn)
spn_15000

# Beta diversity PCoA plot
# Diversity analysis requires rarefied datataxa tables
# rarefaction:
physeq_rar <- rarefy_even_depth(physeq, sample.size = 3000) 

# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

# We define the type of analysis using " ordinate() " and setting the 
# method = PCoA, and also setting the distance to type of beta diversity analysis
# you can change analysis to jaccard, bray-curtis, and so on
ord <- ordinate(physeq_rar_RA, method = "PCoA", distance = "wunifrac")

# Plot data
plot_ordination(physeq_rar_RA, 
                ord, 
                type = "samples", 
                color = "sex", 
                title = "PcoA (Weighted Unifrac)") +
  
  # Adding text to plot
  annotate(geom = "text", 
           label = "",
           x = -0.08,
           y = 0.18,
           size = 4) + 
  
  # Adding elipses
  stat_ellipse(type = "norm", size = 0.5) + # the size value specifies thickness of the ellipses line
  
  guides(colour = guide_legend("Sex")) + # changing legend title
  
  theme_bw(base_size = 14) # changing background colour
                                         
                                         
# ------- NEUTERED + IBD ---------                                        
# Importing all data files
biome_file <- import_biom("/Users/joshuacalalang/Desktop/MICB447_R/Project2/N_IBD/table-with-taxonomy.biom") # recall: need to use quotation marks to specify file path
metadata <- import_qiime_sample_data("/Users/joshuacalalang/Desktop/MICB447_R/Project2/dog_metadata.tsv")
tree <- read_tree_greengenes("/Users/joshuacalalang/Desktop/MICB447_R/Project2/N_IBD/tree.nwk")

# Combining all objects into phyloseq object
physeq <- merge_phyloseq(biome_file, metadata, tree)

# Overview of phyloseq object
physeq 

# head function allows us to only look at first 6 lines
head(sample_data(physeq)) # sample_data() allows us to look at the metadata of phyloseq object

# Set set of random numbers (to make downstream analysis reproducible)
set.seed(711)

# taxonomic rank names
rank_names(physeq)

# Shows us first 6 lines [optional]
head(tax_table(physeq))

# Renaming the column names 
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class",
                                 "Order", "Family", "Genus", "Species")

# Filter data based on metadata category with " == " [useful!!!]
spn <- subset_samples(physeq, spayed_neutered ==  "yes")

# Format data as table with as.data.frame()
as.data.frame(sample_names(spn))

# Show read count for each sample using " sample_sums() "
as.data.frame(sample_sums(spn))

# Exclude samples with less than 15000 reads
spn_15000 <- prune_samples(sample_sums(spn) >= 15000, spn)
spn_15000

# Beta diversity PCoA plot
# Diversity analysis requires rarefied datataxa tables
# rarefaction:
physeq_rar <- rarefy_even_depth(physeq, sample.size = 3000) 

# Convert to RA (relative abundance)
physeq_rar_RA <- transform_sample_counts(physeq_rar, function(x) x/sum(x))

# We define the type of analysis using " ordinate() " and setting the 
# method = PCoA, and also setting the distance to type of beta diversity analysis
# you can change analysis to jaccard, bray-curtis, and so on
ord <- ordinate(physeq_rar_RA, method = "PCoA", distance = "wunifrac")

# Plot data
plot_ordination(physeq_rar_RA, 
                ord, 
                type = "samples", 
                color = "sex", 
                title = "PcoA (Weighted Unifrac)") +
  
  # Adding text to plot
  annotate(geom = "text", 
           label = "",
           x = -0.08,
           y = 0.18,
           size = 4) + 
  
  # Adding elipses
  stat_ellipse(type = "norm", size = 0.5) + # the size value specifies thickness of the ellipses line
  
  guides(colour = guide_legend("Sex")) + # changing legend title
  
  theme_bw(base_size = 14) # changing background colour                                         
                                         
                                         
