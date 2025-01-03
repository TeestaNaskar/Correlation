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
male = meta[meta$Group=='Control' & meta$CSEX =='male',]
rownames(male) = male$Placenta_Seq_ID
#load vst counts
vst = read.csv("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/RNAseq/data/VSTcountdata.human.96subjwithoutSANDY/VST_counts.abovebasemean5.csv")
#subsetting the male vst
rownames(vst) = vst[,1]
vst = vst[,2:ncol(vst)]
#rownames = rownames(male)
#modified_row_names <- sub("[A-Za-z]+$", "", rownames)
#rownames(male) = modified_row_names
malevst = vst[,colnames(vst) %in% rownames(male)]
#subset male vst for selected lipids and immune genes

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
female = meta[meta$Group=='Control' & meta$CSEX =='female',]
rownames(female) = female$Placenta_Seq_ID
#load vst counts
#subsetting the female vst
femalevst = vst[,colnames(vst) %in% rownames(female)]
#subset female vst for selected lipids and immune genes

female.control.vst.sub = femalevst[rownames(femalevst) %in% combined_list,]
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
male = meta[meta$Group=='Cannabis' & meta$CSEX =='male',]
rownames(male) = male$Placenta_Seq_ID
#subsetting vst count
malevst = vst[,colnames(vst) %in% rownames(male)]
#subset male vst for selected lipids and immune genes
#total lipids and cytokines to be correlated

male.cannabis.vst.sub = malevst[rownames(malevst) %in% combined_list,]
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
#next corr for female control
female = meta[meta$Group=='Cannabis' & meta$CSEX =='female',]
rownames(female) = female$Placenta_Seq_ID
#load vst counts
#subsetting the female vst
femalevst = vst[,colnames(vst) %in% rownames(female)]
#subset female vst for selected lipids and immune genes
female.cannabis.vst.sub = femalevst[rownames(femalevst) %in% combined_list,]
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

tiff('controlvscannabis.male.tiff', 
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
         title='control cannabis male',
         tl.col="black",
         tl.cex= 1,
         pch.col='white',
         pch.cex = 2,
         cl.cex=1,
         cl.pos = 'r',
         col.lim= c(-1,1))
dev.off()
##the above corrplot didn't work and since its going too big I have to use complexheatmap function for that
Heatmap(control_merged_male.corr_matrix,
        name = "Correlation",
        col = colorRamp2(c(-1, 0, 1), c("dodgerblue", "white", "red")),
        show_row_names = FALSE,
        show_column_names = FALSE,
        #cluster_rows = F,
        #cluster_columns = F,
        top_annotation = HeatmapAnnotation(lines = anno_lines(control_merged_p_matrix)))

#or
heatmap_result <- pheatmap(control_merged_male.corr_matrix,
         color = colorRampPalette(c("blue", "white", "magenta"))(100),
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         #cluster_cols = TRUE,
         #annotation_row = control_merged_p_matrix,
         #annotation_col = control_merged_p_matrix,
         main = "Control Cannabis Male Correlation Matrix")
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
         title='control vs cannabis female',
         tl.col="black",
         tl.cex= 1,
         pch.col='white',
         pch.cex = 2,
         cl.cex=1,
         cl.pos = 'r',
         col.lim= c(-1,1))
dev.off()
##**********************************************************************##
Heatmap(control_merged_female.corr_matrix,
        name = "Correlation",
        col = colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "magenta")),
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        top_annotation = HeatmapAnnotation(lines = anno_lines(control_merged_female_p_matrix)))
#or
# Install pheatmap if you haven't already
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)
tiff('control_cannabis_male.tiff', 
     width=10, 
     height = 10, 
     units = 'cm', 
     compression ='lzw', res=600)
# Use pheatmap to plot the correlation matrix
pheatmap(control_merged_female.corr_matrix,
         color = colorRampPalette(c("blue", "white", "magenta"))(100),
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         #annotation_row = control_merged_p_matrix,
         #annotation_col = control_merged_p_matrix,
         main = "Control Cannabis Female Correlation Matrix")

pheatmap(control_merged_male.corr_matrix,
         color = colorRampPalette(c("dodgerblue3", "white", "magenta"))(100),
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         #annotation_row = control_merged_p_matrix,
         #annotation_col = control_merged_p_matrix,
         main = "Control Cannabis Male Correlation Matrix")
#### if i want to modify the colors in prizm then I need to open them in prizm for that exporting the data from R to open in prizm
clustered_rows <- heatmap_result$tree_row
clustered_cols <- heatmap_result$tree_col

# Reorder rows and columns in the original dataframe
df_clustered <- control_merged_male.corr_matrix[clustered_rows, clustered_cols]

# Save the clustered dataframe as a CSV file
write.csv(df_clustered, "clustered_dataframe.csv", row.names = TRUE)

row_hclust <- hclust(as.dist(1 - control_merged_male.corr_matrix))
col_hclust <- hclust(as.dist(1 - t(control_merged_male.corr_matrix)))

# Reorder the correlation matrix based on clustering
row_order <- row_hclust$order
col_order <- col_hclust$order

control_merged_male.corr_matrix <- control_merged_male.corr_matrix[row_order, col_order]
