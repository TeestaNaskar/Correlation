
###run regular correlations
#turn off warnings. To turn back on change value to 0
options(warn = -1)

library(openxlsx)
library(ppcor)
library(igraph)
library(impute)
library(qgraph)
library(DescTools)
library(Hmisc)
library(dplyr)
library(tidyr)
library(corrplot)
library(stringr)
library(Hmisc)
library(pheatmap)

setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/partial correlation immune lipids")
##I downloaded the full lists of the genes associated with term GOBP_PHOSPHOLIPID_METABOLIC_PROCESS
#total lipids and cytokines to be correlated
#download all the genes for the significant modules named OBP_PHOSPHOLIPID_METABOLIC_PROCESS
phospholipid = read.xlsx("GOBP_PHOSPHOLIPID_METABOLIC_PROCESS.v2023.2.Hs.xlsx", sheet = 2)
# Extract the GENE_SYMBOLS column as a single string
gene_symbols_string <- phospholipid$GENE_SYMBOLS[1]
# Split the string into individual gene names, removing any leading or trailing whitespace
genes1 <- strsplit(gene_symbols_string, ",")[[1]]
genes1 <- trimws(genes1)
#genes = phospholipid$GENE_SYMBOLS
#phospholipid = unlist(strsplit(genes, split = ","))
choline = read.xlsx("GOBP_PHOSPHOLIPID_METABOLIC_PROCESS.v2023.2.Hs.xlsx", sheet = 4)
# Extract the GENE_SYMBOLS column as a single string
gene_symbols_string <- choline$GENE_SYMBOLS[1]
# Split the string into individual gene names, removing any leading or trailing whitespace
genes2 <- strsplit(gene_symbols_string, ",")[[1]]
genes2 <- trimws(genes2)

chemokines = read.xlsx("GOBP_PHOSPHOLIPID_METABOLIC_PROCESS.v2023.2.Hs.xlsx", sheet = 6)
# Extract the GENE_SYMBOLS column as a single string
gene_symbols_string <- chemokines$GENE_SYMBOLS[1]
# Split the string into individual gene names, removing any leading or trailing whitespace
genes3 <- strsplit(gene_symbols_string, ",")[[1]]
genes3 <- trimws(genes3)

defense = read.xlsx("GOBP_PHOSPHOLIPID_METABOLIC_PROCESS.v2023.2.Hs.xlsx", sheet = 8)
# Extract the GENE_SYMBOLS column as a single string
gene_symbols_string <- defense$GENE_SYMBOLS[1]
# Split the string into individual gene names, removing any leading or trailing whitespace
genes4 <- strsplit(gene_symbols_string, ",")[[1]]
genes4 <- trimws(genes4)

#combine all four list of genes
list_of_genes <- list(genes1, genes2, genes3, genes4)
# Use Reduce to combine all lists and keep only unique items
combined_list <- Reduce(union, list_of_genes)
print(combined_list)
#choose your group of interest and run the partial corrs
#for male, load first vst counts, subsets male and then select only lipids and immune genes
#load metadata
meta = read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/RNAseq/data/HUMAN.PLACENTA.METADATA/WorkingMetadata.updatedbyTeesta.consolidatedinfo.Anissa.Greg.Placenta.inventory.xlsx", sheet = 3)
meta = meta[1:93,]
#male control
male.control <- meta[meta$Group=='Control' & meta$CSEX =='male',]
rownames(male.control) = male.control$Placenta_Seq_ID
#male cannabis
male.cannabis <- meta[meta$Group=='Cannabis' & meta$CSEX =='male',]
rownames(male.cannabis) = male.cannabis$Placenta_Seq_ID
##female control
female.control <- meta[meta$Group=='Control' & meta$CSEX =='female',]
rownames(female.control) = female.control$Placenta_Seq_ID
#female cannabis
female.cannabis <- meta[meta$Group=='Cannabis' & meta$CSEX =='female',]
rownames(female.cannabis) = female.cannabis$Placenta_Seq_ID
#load vst counts
vst = read.csv("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/RNAseq/data/VSTcountdata.human.96subjwithoutSANDY/VST_counts.abovebasemean5.csv")
#subsetting the male vst
rownames(vst) = vst[,1]
vst = vst[,2:ncol(vst)]
# Subset and transpose the malevst data for control and cannabis conditions
male.control.vst <- vst[,colnames(vst) %in% rownames(male.control)]
male.control.vst.sub <- male.control.vst[rownames(male.control.vst) %in% combined_list,]
male.control.vst.sub <- t(male.control.vst.sub)

male.cannabis.vst <- vst[,colnames(vst) %in% rownames(male.cannabis)]
male.cannabis.vst.sub <- male.cannabis.vst[rownames(male.cannabis.vst) %in% combined_list,]
male.cannabis.vst.sub <- t(male.cannabis.vst.sub)


# Calculate correlation matrices for control and cannabis conditions
control_corrs <- rcorr(as.matrix(male.control.vst.sub), type = "pearson")
cannabis_corrs <- rcorr(as.matrix(male.cannabis.vst.sub), type = "pearson")

# Extract correlation matrices
control_corr_matrix <- control_corrs$r
cannabis_corr_matrix <- cannabis_corrs$r

# Perform hierarchical clustering on the upper triangle of the control correlation matrix
dist_control <- as.dist((1 - control_corr_matrix) / 2)
hclust_control <- hclust(dist_control, method = "complete")
order_control <- hclust_control$order

# Reorder the control and cannabis correlation matrices according to the clustering order
reordered_control_corr_matrix <- control_corr_matrix[order_control, order_control]
reordered_cannabis_corr_matrix <- cannabis_corr_matrix[order_control, order_control]

# Initialize the combined matrix with the reordered control correlation matrix
combined_corr_matrix <- reordered_control_corr_matrix

# Combine the upper triangle of the reordered control matrix with the lower triangle of the reordered cannabis matrix
combined_corr_matrix[lower.tri(combined_corr_matrix)] <- reordered_cannabis_corr_matrix[lower.tri(reordered_cannabis_corr_matrix)]

# Visualize the combined correlation matrix using pheatmap
pheatmap(combined_corr_matrix, cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("dodgerblue", "white", "magenta"))(100),
         show_rownames = FALSE,
         show_colnames = FALSE)

# Subset and transpose the femalevst data for control and cannabis conditions
female.control.vst <- vst[,colnames(vst) %in% rownames(female.control)]
female.control.vst.sub <- female.control.vst[rownames(female.control.vst) %in% combined_list,]
female.control.vst.sub <- t(female.control.vst.sub)

female.cannabis.vst <- vst[,colnames(vst) %in% rownames(female.cannabis)]
female.cannabis.vst.sub <- female.cannabis.vst[rownames(female.cannabis.vst) %in% combined_list,]
female.cannabis.vst.sub <- t(female.cannabis.vst.sub)

# Calculate correlation matrices for control and cannabis conditions
control_corrs <- rcorr(as.matrix(female.control.vst.sub), type = "pearson")
cannabis_corrs <- rcorr(as.matrix(female.cannabis.vst.sub), type = "pearson")

# Extract correlation matrices
control_corr_matrix <- control_corrs$r
cannabis_corr_matrix <- cannabis_corrs$r

# Perform hierarchical clustering on the upper triangle of the control correlation matrix
dist_control <- as.dist((1 - control_corr_matrix) / 2)
hclust_control <- hclust(dist_control, method = "complete")
order_control <- hclust_control$order

# Reorder the control and cannabis correlation matrices according to the clustering order
reordered_control_corr_matrix <- control_corr_matrix[order_control, order_control]
reordered_cannabis_corr_matrix <- cannabis_corr_matrix[order_control, order_control]

# Initialize the combined matrix with the reordered control correlation matrix
combined_corr_matrix <- reordered_control_corr_matrix

# Combine the upper triangle of the reordered control matrix with the lower triangle of the reordered cannabis matrix
combined_corr_matrix[lower.tri(combined_corr_matrix)] <- reordered_cannabis_corr_matrix[lower.tri(reordered_cannabis_corr_matrix)]

# Visualize the combined correlation matrix using pheatmap
pheatmap(combined_corr_matrix, cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("blue", "white", "magenta"))(100),
         show_rownames = FALSE,
         show_colnames = FALSE)
