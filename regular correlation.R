###run regular correlations
setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/Rat_placenta/correlation/")
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

##I downloaded the full lists of the genes associated with term GOBP_PHOSPHOLIPID_METABOLIC_PROCESS
#total lipids and cytokines to be correlated
#download all the genes for the significant modules named OBP_PHOSPHOLIPID_METABOLIC_PROCESS
phospholipid = read.xlsx("MSigDB_GOBP_MICE_DATA.xlsx", sheet = 2)
# Extract the GENE_SYMBOLS column as a single string
gene_symbols_string <- phospholipid$GENE_SYMBOLS[1]
# Split the string into individual gene names, removing any leading or trailing whitespace
genes1 <- strsplit(gene_symbols_string, ",")[[1]]
genes1 <- trimws(genes1)
#genes = phospholipid$GENE_SYMBOLS
#phospholipid = unlist(strsplit(genes, split = ","))
choline = read.xlsx("MSigDB_GOBP_MICE_DATA.xlsx", sheet = 4)
# Extract the GENE_SYMBOLS column as a single string
gene_symbols_string <- choline$GENE_SYMBOLS[1]
# Split the string into individual gene names, removing any leading or trailing whitespace
genes2 <- strsplit(gene_symbols_string, ",")[[1]]
genes2 <- trimws(genes2)

chemokines = read.xlsx("MSigDB_GOBP_MICE_DATA.xlsx", sheet = 6)
# Extract the GENE_SYMBOLS column as a single string
gene_symbols_string <- chemokines$GENE_SYMBOLS[1]
# Split the string into individual gene names, removing any leading or trailing whitespace
genes3 <- strsplit(gene_symbols_string, ",")[[1]]
genes3 <- trimws(genes3)

defense = read.xlsx("MSigDB_GOBP_MICE_DATA.xlsx", sheet = 8)
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

#choose your group of interest and run the partial corrs
#for male, load first vst counts, subsets male and then select only lipids and immune genes
#load metadata

meta = read.table("/Users/teestanaskar/Dropbox/Teesta/Placenta/Rat_placenta/RNAseq/data/metadata.txt")
colnames(meta) = meta[1,]
meta = meta[2:nrow(meta),]
rownames(meta) = meta$ID
#rownames(meta) = gsub("X", "", rownames(meta))
male = meta[meta$Group == 'VEH' & meta$Sex =='Male',]
#load vst counts
vst = read.csv("data/../DEG_STARalignedcounts/bothsex.rat.VST_counts.csv")
#subsetting the male vst
rownames(vst) = vst[,1]
vst = vst[,2:ncol(vst)]
#rownames = rownames(male)
#modified_row_names <- sub("[A-Za-z]+$", "", rownames)
#rownames(male) = modified_row_names
malevst = vst[,colnames(vst) %in% rownames(male)]
#subset male vst for selected lipids and immune genes
#total lipids and cytokines to be correlated
#lipids_immune = c("Trpm8","Trpv1", "Trpv3","Fasn","Daglb","Ifngr1","Phospho1","Abhd4", "Tgfb1l1", "Acsl5", "Il33", "Cnr2", "Ifnar1",  
                 # "Il18", "Cds1", "Ptgs1", "Ccl2", "Il1rapl2", "Il17f", "Ccr1l1", "Il36b", "Pltp", "Cxcr1")
male.control.vst.sub = malevst[rownames(malevst) %in% combined_list,]
male.control.vst.sub = t(male.control.vst.sub)
corrs = rcorr(as.matrix(male.control.vst.sub), type = "pearson")
#corrs = rcorr(as.matrix(lipids_immune[,male.control.vst.sub]),type="pearson")
rownames(corrs$r)=colnames(male.control.vst.sub)
colnames(corrs$r)=colnames(male.control.vst.sub)
rownames(corrs$P)=colnames(male.control.vst.sub)
colnames(corrs$P)=colnames(male.control.vst.sub)
control_male_corrs = corrs
control_male_ps = control_male_corrs$P
control_male_ps[is.nan(control_male_ps)]=1
#next corr for female control
female = meta[meta$Group == 'VEH' & meta$Sex =='Female',]
#rownames(female) = female$Placenta_Seq_ID
#load vst counts
#subsetting the female vst
femalevst = vst[,colnames(vst) %in% rownames(female)]
#subset female vst for selected lipids and immune genes
female.control.vst.sub = femalevst[rownames(femalevst) %in% lipids_immune,]
female.control.vst.sub = t(female.control.vst.sub)
corrs = rcorr(as.matrix(female.control.vst.sub),type="pearson")
rownames(corrs$r)=colnames(female.control.vst.sub)
colnames(corrs$r)=colnames(female.control.vst.sub)
rownames(corrs$P)=colnames(female.control.vst.sub)
colnames(corrs$P)=colnames(female.control.vst.sub)
control_female_corrs = corrs
control_female_ps = control_female_corrs$P
control_female_ps[is.nan(control_female_ps)]=1

#total correlation datasets from controls 
#1) control_male_corrs 2)control_female_corrs

## same for cannabis##########################
male = meta[meta$Group=='THC_CBD' & meta$Sex =='Male',]
#rownames(male) = male$Placenta_Seq_ID
#subsetting vst count
malevst = vst[,colnames(vst) %in% rownames(male)]
#subset male vst for selected lipids and immune genes
#total lipids and cytokines to be correlated
male.cannabis.vst.sub = malevst[rownames(malevst) %in% lipids_immune,]
male.cannabis.vst.sub = t(male.cannabis.vst.sub)
corrs = rcorr(as.matrix(male.cannabis.vst.sub), type = "pearson")
#corrs = rcorr(as.matrix(lipids_immune[,male.cannabis.vst.sub]),type="pearson")
rownames(corrs$r)=colnames(male.cannabis.vst.sub)
colnames(corrs$r)=colnames(male.cannabis.vst.sub)
rownames(corrs$P)=colnames(male.cannabis.vst.sub)
colnames(corrs$P)=colnames(male.cannabis.vst.sub)
cannabis_male_corrs = corrs
cannabis_male_ps = cannabis_male_corrs$P
cannabis_male_ps[is.nan(cannabis_male_ps)]=1
#next corr for female cannabis/THC-CBD
female = meta[meta$Group=='THC_CBD' & meta$Sex =='Female',]
#rownames(female) = female$Placenta_Seq_ID
#load vst counts
#subsetting the female vst
femalevst = vst[,colnames(vst) %in% rownames(female)]
#subset female vst for selected lipids and immune genes
female.cannabis.vst.sub = femalevst[rownames(femalevst) %in% lipids_immune,]
female.cannabis.vst.sub = t(female.cannabis.vst.sub)
corrs = rcorr(as.matrix(female.cannabis.vst.sub),type="pearson")
rownames(corrs$r)=colnames(female.cannabis.vst.sub)
colnames(corrs$r)=colnames(female.cannabis.vst.sub)
rownames(corrs$P)=colnames(female.cannabis.vst.sub)
colnames(corrs$P)=colnames(female.cannabis.vst.sub)
cannabis_female_corrs = corrs
cannabis_female_ps = cannabis_female_corrs$P
cannabis_female_ps[is.nan(cannabis_female_ps)]=1

#total correlation datasets from cannabis 
#1) cannabis_male_corrs 2) cannabis_female_corrs
#four correlations sets for merging plots
#1) control_male_corrs 2)control_female_corrs 3) cannabis_male_corrs 4) cannabis_female_corrs
#############################create plots#########################################
#merged matrix plots
control_merged_male.corr_matrix = control_male_corrs$r
control_corr_UT = upper.tri(control_merged_male.corr_matrix)
control_merged_male.corr_matrix[control_corr_UT] = cannabis_male_corrs$r[control_corr_UT]

control_merged_p_matrix = control_male_ps
control_p_UT = upper.tri(control_merged_p_matrix)
control_merged_p_matrix[control_p_UT] = cannabis_male_ps[control_p_UT]

tiff('control.tiff', 
     width=10, 
     height = 10, 
     units = 'cm', 
     compression ='lzw', res=600)
corrplot(control_merged_male.corr_matrix,
         method='circle',
         is.corr = F,
         p.mat=control_merged_p_matrix,
         insig= "label_sig",
         col= COL2('RdBu', 100),
         order='original',
         diag=F,
         mar=c(1,1,1,1),
         title='Rat VEH THC male',
         tl.col="black",
         tl.cex= 1,
         pch.col='white',
         pch.cex = 2,
         cl.cex=1,
         cl.pos = 'r',
         col.lim= c(-1,1))
dev.off()
#############################create plots#########################################
#merged matrix plots
control_merged_female.corr_matrix = control_female_corrs$r
control_corr_UT = upper.tri(control_merged_female.corr_matrix)
control_merged_female.corr_matrix[control_corr_UT] = cannabis_female_corrs$r[control_corr_UT]

control_merged_female_p_matrix = control_female_ps
control_p_UT = upper.tri(control_merged_female_p_matrix)
control_merged_female_p_matrix[control_p_UT] = cannabis_female_ps[control_p_UT]

tiff('control.tiff', 
     width=10, 
     height = 10, 
     units = 'cm', 
     compression ='lzw', res=600)
corrplot(control_merged_female.corr_matrix,
         method='circle',
         is.corr = F,
         p.mat=control_merged_female_p_matrix,
         insig= "label_sig",
         col= COL2('RdBu', 100),
         order='original',
         diag=F,
         mar=c(1,1,1,1),
         title='Rat Veh THC female',
         tl.col="black",
         tl.cex= 1,
         pch.col='white',
         pch.cex = 2,
         cl.cex=1,
         cl.pos = 'r',
         col.lim= c(-1,1))
dev.off()
##**********************************************************************##
#### Perform hierarchical clustering
# for male after the merging matrix step (marged matrix plots)
row_hclust <- hclust(as.dist(1 - control_merged_male.corr_matrix))
col_hclust <- hclust(as.dist(1 - t(control_merged_male.corr_matrix)))

# Reorder the correlation matrix based on clustering
row_order <- row_hclust$order
col_order <- col_hclust$order

# Reorder the correlation matrix and p-value matrix
control_merged_male.corr_matrix <- control_merged_male.corr_matrix[row_order, col_order]
control_merged_p_matrix <- control_merged_p_matrix[row_order, col_order]

# Plot the reordered correlation matrix
tiff('control.tiff', 
     width = 10, 
     height = 10, 
     units = 'cm', 
     compression = 'lzw', 
     res = 600)

corrplot(control_merged_male.corr_matrix,
         method = 'circle',
         is.corr = FALSE,
         p.mat = control_merged_p_matrix,
         insig = "label_sig",
         col = COL2('RdBu', 100),
         order = 'original',
         diag = FALSE,
         mar = c(1, 1, 1, 1),
         title = 'control cannabis male',
         tl.col = "black",
         tl.cex = 1,
         pch.col = 'white',
         pch.cex = 2,
         cl.cex = 1,
         cl.pos = 'r',
         col.lim = c(-1, 1))

dev.off()
#### hierarchial clustering for female
### after merging matrix step (merged matrix plots)
row_hclust <- hclust(as.dist(1 - control_merged_female.corr_matrix))
col_hclust <- hclust(as.dist(1 - t(control_merged_female.corr_matrix)))

# Reorder the correlation matrix based on clustering
row_order <- row_hclust$order
col_order <- col_hclust$order

# Reorder the correlation matrix and p-value matrix
control_merged_female.corr_matrix <- control_merged_female.corr_matrix[row_order, col_order]
control_merged_female_p_matrix <- control_merged_female_p_matrix[row_order, col_order]

# Plot the reordered correlation matrix
tiff('control.tiff', 
     width = 10, 
     height = 10, 
     units = 'cm', 
     compression = 'lzw', 
     res = 600)

corrplot(control_merged_female.corr_matrix,
         method = 'circle',
         is.corr = FALSE,
         p.mat = control_merged_female_p_matrix,
         insig = "label_sig",
         col = COL2('RdBu', 100),
         order = 'original',
         diag = FALSE,
         mar = c(1, 1, 1, 1),
         title = 'control cannabis female',
         tl.col = "black",
         tl.cex = 1,
         pch.col = 'white',
         pch.cex = 2,
         cl.cex = 1,
         cl.pos = 'r',
         col.lim = c(-1, 1))
