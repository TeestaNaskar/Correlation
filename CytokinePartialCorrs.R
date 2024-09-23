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

setwd("/Users/teestanaskar/Dropbox/Cannabis Meeting and Projects/pnTHC/Anissa's projects/Prenatal Cannabis/Prenatal THC+CBD/Molecular/Cytokines levels")
setwd("/Users/teestanaskar/Downloads/")
#load and subset data
behavior = read.csv('CondProbByRat.csv')
metadata = subset(behavior[,1:5])
behav_only = behavior[,6:ncol(behavior)] #pull all lipids for easy selection

#choose whether you want to correlate NacSh central lipids or peripheral

veh_behav = behavior[behavior$Group == 'VEH',]
THC_behav = behavior[behavior$Group == 'THC+CBD',]


#choose your group of interest and run the partial corrs
data = veh_behav
corrs = pcor(data[,6:ncol(data)])
rownames(corrs$estimate)=colnames(data)[6:ncol(data)]
colnames(corrs$estimate)=colnames(data)[6:ncol(data)]
rownames(corrs$p.value)=colnames(data)[6:ncol(data)]
colnames(corrs$p.value)=colnames(data)[6:ncol(data)]
veh_corrs = corrs
veh_ps = veh_corrs$p.value
veh_ps[is.nan(veh_ps)]=1

data = THC_behav
corrs = pcor(data[,6:ncol(data)])
rownames(corrs$estimate)=colnames(data)[6:ncol(data)]
colnames(corrs$estimate)=colnames(data)[6:ncol(data)]
rownames(corrs$p.value)=colnames(data)[6:ncol(data)]
colnames(corrs$p.value)=colnames(data)[6:ncol(data)]
THC_corrs = corrs
THC_ps = THC_corrs$p.value
THC_ps[is.nan(THC_ps)]=1



#############################create plots#########################################
#merged matrix plots
merged_corr_matrix = veh_corrs$estimate
corr_UT = upper.tri(merged_corr_matrix)
merged_corr_matrix[corr_UT] = THC_corrs$estimate[corr_UT]

merged_p_matrix = veh_ps
p_UT = upper.tri(merged_p_matrix)
merged_p_matrix[p_UT] = THC_ps[p_UT]

# tiff('neutralrats_NacShLipids_mergedMatrix.tiff', 
#      width=10, 
#      height = 10, 
#      units = 'cm', 
#      compression ='lzw', res=600)
corrplot(merged_corr_matrix,
         method='circle',
         is.corr = F,
         p.mat=merged_p_matrix,
         insig= "label_sig",
         col= COL2('PRGn', 100),
         order='original',
         diag=F,
         mar=c(1,1,1,.5) +0.3,
         title='PartialCors: BottomLeft-Veh, TopRight-THC',
         tl.col="black",
         tl.cex= 1,
         pch.col='white',
         pch.cex = 2,
         cl.cex=.8,
         cl.pos = 'b')
#dev.off()



###run regular correlations
#choose your group of interest and run the partial corrs
data = veh_behav
corrs = rcorr(as.matrix(data[,6:ncol(data)]),type="pearson")
rownames(corrs$r)=colnames(data)[6:ncol(data)]
colnames(corrs$r)=colnames(data)[6:ncol(data)]
rownames(corrs$P)=colnames(data)[6:ncol(data)]
colnames(corrs$P)=colnames(data)[6:ncol(data)]
veh_reg_corrs = corrs
veh_reg_ps = veh_reg_corrs$P
veh_reg_ps[is.nan(veh_reg_ps)]=1

data = THC_behav
corrs = rcorr(as.matrix(data[,6:ncol(data)]),type="pearson")
rownames(corrs$r)=colnames(data)[6:ncol(data)]
colnames(corrs$r)=colnames(data)[6:ncol(data)]
rownames(corrs$P)=colnames(data)[6:ncol(data)]
colnames(corrs$P)=colnames(data)[6:ncol(data)]
THC_reg_corrs = corrs
THC_reg_ps = THC_reg_corrs$P
THC_reg_ps[is.nan(THC_reg_ps)]=1


#############################create plots#########################################
#merged matrix plots
merged_reg_corr_matrix = veh_reg_corrs$r
corr_UT = upper.tri(merged_reg_corr_matrix)
merged_reg_corr_matrix[corr_UT] = THC_reg_corrs$r[corr_UT]

merged_reg_p_matrix = veh_reg_ps
p_UT = upper.tri(merged_p_matrix)
merged_p_matrix[p_UT] = THC_reg_ps[p_UT]

# tiff('neutralrats_NacShLipids_mergedMatrix.tiff', 
#      width=10, 
#      height = 10, 
#      units = 'cm', 
#      compression ='lzw', res=600)
corrplot(merged_reg_corr_matrix,
         method='circle',
         is.corr = F,
         p.mat=merged_reg_p_matrix,
         insig= "label_sig",
         col= COL2('PRGn', 100),
         order='original',
         diag=F,
         mar=c(1,1,1,.5) +0.3,
         title='RegularCors: BottomLeft-Veh, TopRight-THC',
         tl.col="black",
         tl.cex= 1,
         pch.col='white',
         pch.cex = 2,
         cl.cex=.8,
         cl.pos = 'b')
#dev.off()






