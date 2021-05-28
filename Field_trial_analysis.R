## Analysis of bacterial and fungal communities in field trial
## Reena Debray
## April 2021

# Load packages
library(phyloseq)
library(vegan)
library(indicspecies)
library(igraph)
library(Hmisc)
library(philentropy)

### Bacterial alpha and beta diversity

# Subset bacterial phyloseq object to field samples
ag_ps_field<-subset_samples(ag_ps_nocontam_allsamp,Group=="Field")
# Remove taxa that do not occur in any of the field samples
nosamp_taxa<-taxa_names(ag_ps_field)[apply(ag_ps_field@otu_table,1,function(x){mean(x)==0})]
ag_ps_rmnosamp<-prune_taxa(setdiff(taxa_names(ag_ps_field),nosamp_taxa),ag_ps_field)

# Calculate bacterial alpha diversity
ag_ps_rmnosamp@sam_data$observed<-unlist(estimate_richness(ag_ps_rmnosamp,measures="Observed"))
ag_ps_rmnosamp@sam_data$shannon<-unlist(estimate_richness(ag_ps_rmnosamp,measures="Shannon"))
ag_ps_rmnosamp@sam_data$evenness<-ag_ps_rmnosamp@sam_data$shannon/log(ag_ps_rmnosamp@sam_data$observed)

# Calculate bacterial beta diversity within drought treatments
# get_pairwise_dist() function can be found in the file Functions.R
dist_bray<-as.matrix(phyloseq::distance(ag_ps_rmnosamp,method="bray"))
FW_names<-rownames(ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Water.Regiment=="full water",])
D_names<-rownames(ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Water.Regiment=="drought",])
drought_dispersion_analysis<-get_pairwise_dist(dist_bray,list(FW_names,D_names),c("full water","drought"))

# Calculate bacterial beta diversity within genotypes
WT_names<-rownames(ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Genotype=="76R",])
rmc_names<-rownames(ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Genotype=="rmc",])
genotype_dispersion_analysis<-get_pairwise_dist(dist_bray,list(WT_names,rmc_names),c("76R","rmc"))

# Calculate bacterial beta diversity within P levels
highP_names<-rownames(ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Phosphorous.Level==3,])
medP_names<-rownames(ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Phosphorous.Level==2,])
lowP_names<-rownames(ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Phosphorous.Level==1,])
P_dispersion_analysis<-get_pairwise_dist(dist_bray,list(highP_names,medP_names,lowP_names),c("high","medium","low"))

### Fungal alpha and beta diversity

# Subset fungal phyloseq object to field samples
ps_field_noNA <- subset_samples(ps_noNA_ITSx_nocontam, Group=="Field")
# Remove taxa that do not occur in any of the field samples
nosamp_taxa<-taxa_names(ps_field_noNA)[apply(ps_field_noNA@otu_table,1,function(x){mean(x)==0})]
ps_field_noNA<-prune_taxa(setdiff(taxa_names(ps_field_noNA),nosamp_taxa),ps_field_noNA)

# Calculate fungal alpha diversity
ps_field_noNA@sam_data$observed<-unlist(estimate_richness(ps_field_noNA,measures="Observed"))
ps_field_noNA@sam_data$shannon<-unlist(estimate_richness(ps_field_noNA,measures="Shannon"))
ps_field_noNA@sam_data$evenness<-ps_field_noNA@sam_data$shannon/log(ps_field_noNA@sam_data$observed)

# Calculate fungal beta diversity within drought treatments
dist_bray_ITS<-as.matrix(phyloseq::distance(ps_field_noNA,method="bray"))
FW_names<-rownames(ps_field_noNA@sam_data[ps_field_noNA@sam_data$Water.Regiment=="full water",])
D_names<-rownames(ps_field_noNA@sam_data[ps_field_noNA@sam_data$Water.Regiment=="drought",])
drought_dispersion_ITS_analysis<-get_pairwise_dist(dist_bray_ITS,list(FW_names,D_names),c("full water","drought"))

# Calculate fungal beta diversity within genotypes
WT_names<-rownames(ps_field_noNA@sam_data[ps_field_noNA@sam_data$Genotype=="76R",])
rmc_names<-rownames(ps_field_noNA@sam_data[ps_field_noNA@sam_data$Genotype=="rmc",])
genotype_dispersion_ITS_analysis<-get_pairwise_dist(dist_bray_ITS,list(WT_names,rmc_names),c("Wild-type 76R","Reduced\nmycorrhizal"))

# Calculate fungal beta diversity within P levels
highP_names<-rownames(ps_field_noNA@sam_data[ps_field_noNA@sam_data$Phosphorous.Level==3,])
medP_names<-rownames(ps_field_noNA@sam_data[ps_field_noNA@sam_data$Phosphorous.Level==2,])
lowP_names<-rownames(ps_field_noNA@sam_data[ps_field_noNA@sam_data$Phosphorous.Level==1,])
P_dispersion_ITS_analysis<-get_pairwise_dist(dist_bray_ITS,list(highP_names,medP_names,lowP_names),c("High","Medium","Low"))

### Indicator species analysis

# Indicator species analysis for water regimen on wild-type background
abund<-t(subset_samples(ag_ps_rmnosamp,Genotype=="76R")@otu_table)
water<-subset_samples(ag_ps_rmnosamp,Genotype=="76R")@sam_data$Water.Regiment
inv_water = as.matrix(multipatt(abund,water,func="IndVal.g",control = how(nperm=999)))[,1]$sign

# Indicator species analysis for genotype on full-water background
abund<-t(subset_samples(ag_ps_rmnosamp,Water.Regiment=="full water")@otu_table)
genotype<-subset_samples(ag_ps_rmnosamp,Water.Regiment=="full water")@sam_data$Genotype
inv_genotype = as.matrix(multipatt(abund,genotype,func="IndVal.g",control = how(nperm=999)))[,1]$sign

# Identify ASVs that were significant hits for either water regimen or genotype (or both)
all_indic<-unique(c(rownames(inv_water[inv_water$p.value<0.05 & !is.na(inv_water$p.value),]),rownames(inv_genotype[inv_genotype$p.value<0.05 & !is.na(inv_genotype$p.value),])))
water_sig<-inv_water[all_indic,]
genotype_sig<-inv_genotype[all_indic,]
colnames(water_sig)=c("Var1","Var2","index","stat","pval")
water_sig$ASV<-rownames(water_sig)
water_sig$treatment<-rep("water",nrow(water_sig))
colnames(genotype_sig)=c("Var1","Var2","index","stat","pval") 
genotype_sig$ASV<-rownames(genotype_sig)
genotype_sig$treatment<-rep("genotype",nrow(genotype_sig))
indic_asvs<-rbind(water_sig,genotype_sig)

### Fungal-bacterial correlations

# Richness
ITS_16S_df<-data.frame(ps_field_noNA@sam_data$observed,ps_field_noNA@sam_data$shannon,ag_ps_rmnosamp@sam_data$observed,ag_ps_rmnosamp@sam_data$shannon,ps_field_noNA@sam_data$extraction_batch,ps_field_noNA@sam_data$processing_batch,ps_field_noNA@sam_data$Genotype,ps_field_noNA@sam_data$Water.Regiment,ps_field_noNA@sam_data$Phosphorous.Level)
colnames(ITS_16S_df)=c("fungal_richness","fungal_shannon","bac_richness","bac_shannon","extraction_batch","processing_batch","Genotype","Water.Regiment","Phosphorous.Level")
cor.test(ITS_16S_df$fungal_richness,ITS_16S_df$bac_richness)

# Richness (partial R-squared, controlling for study variables)
full<-lm(bac_richness~fungal_richness+Water.Regiment+Genotype+Phosphorous.Level+extraction_batch+processing_batch,ITS_16S_df)
reduced<-lm(bac_richness~Water.Regiment+Genotype+Phosphorous.Level+extraction_batch+processing_batch,ITS_16S_df)
rsq.partial(full,reduced)

# Community composition
dist_bray<-as.matrix(phyloseq::distance(ag_ps_rmnosamp,method="bray"))
dist_bray_ITS<-as.matrix(phyloseq::distance(ps_field_noNA,method="bray"))
mantel(dist_bray,dist_bray_ITS)

# Community composition (partial Mantel, controlling for study variables)
treatments<-ag_ps_rmnosamp@sam_data[,c("Phosphorous.Level","Genotype","Water.Regiment","processing_batch","extraction_batch")]
treatments$Genotype<-as.numeric(as.factor(treatments$Genotype))
treatments$Water.Regiment<-as.numeric(as.factor(treatments$Water.Regiment))
treatments$processing_batch<-as.numeric(as.factor(treatments$processing_batch))
treatments$extraction_batch<-as.numeric(as.factor(treatments$extraction_batch))
treatment_matrix<-philentropy::distance(data.frame(treatments),method="euclidean")
mantel.partial(dist_bray,dist_bray_ITS,treatment_matrix)

### Co-occurrence network

# Convert to presence-absence and combine top 500 bacterial and top 500 fungal ASVs
# bacteria
ag_ps_rmnosamp@sam_data$Treatment<-paste(ag_ps_rmnosamp@sam_data$Phosphorous.Level,paste(ag_ps_rmnosamp@sam_data$Genotype,ag_ps_rmnosamp@sam_data$Water.Regiment,sep="_"),sep="P_")
ps_nw_bacteria<-phyloseq_standardize_otu_abundance(ag_ps_rmnosamp, method = "pa")
topNbac<-names(sort(taxa_sums(ps_nw_bacteria),TRUE)[1:500])
ps_nw_bacteria<-prune_taxa(topNbac,ps_nw_bacteria)
taxa_names(ps_nw_bacteria)=paste("B_",taxa_names(ps_nw_bacteria),sep="")
# fungi
ps_field_noNA@sam_data$Treatment<-paste(ps_field_noNA@sam_data$Phosphorous.Level,paste(ps_field_noNA@sam_data$Genotype,ps_field_noNA@sam_data$Water.Regiment,sep="_"),sep="P_")
ps_nw_fungi<-phyloseq_standardize_otu_abundance(ps_field_noNA, method = "pa")
topNfun<-names(sort(taxa_sums(ps_nw_fungi),TRUE)[1:500])
ps_nw_fungi<-prune_taxa(topNfun,ps_nw_fungi)
taxa_names(ps_nw_fungi)=paste("F_",taxa_names(ps_nw_fungi),sep="")
# combine
ps_nw_otu<-as.matrix(rbind(ps_nw_bacteria@otu_table,ps_nw_fungi@otu_table))

# Use edge_metrics and node_metrics functions (see Functions.R file) to calculate network metrics across water-genotype treatment combinations
nw_sums<-data.frame(matrix(nrow=0,ncol=14))
node_sums<-data.frame(matrix(nrow=0,ncol=4))
nw_list<-list()

i=1
for (treatment in unique(paste(ps_nw_bacteria@sam_data$Water.Regiment,ps_nw_bacteria@sam_data$Genotype))){
  sample_names<-rownames(ps_nw_bacteria@sam_data[ps_nw_bacteria@sam_data$Water_Genotype==treatment,])

  ps_nw_subset<-ps_nw_otu[,sample_names]
  ps_nw_subset<-ps_nw_subset[which(rowSums(ps_nw_subset)!=0),]

  ps_nw_subset_corr<-rcorr(t(ps_nw_subset))
  tmp_r<-ps_nw_subset_corr$r
  tmp_p<-ps_nw_subset_corr$P

  #Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction ("FDR-BH")
  tmp_padjust <- p.adjust(tmp_p, method="BH")
  p.cutoff=0.1
  tmp_r[tmp_r==1]=0
  tmp_r[is.na(tmp_r)]=0
  tmp_r[which(tmp_padjust>p.cutoff)]=0
  
  # delete those rows and columns with sum = 0 (note positive and negative corrs might cancel each other out)
  tmp_r<-tmp_r[which(rowSums(tmp_r)!=0),]
  tmp_r<-tmp_r[,which(colSums(tmp_r)!=0)]
  
  #calculate node metrics
  degrees<-node_metrics(tmp_r)
  treatment_df<-data.frame(rep(treatment,nrow(degrees)),degrees)
  node_sums<-rbind(node_sums,treatment_df)
  
  #calculate network metrics
  nw_sums[i,1]=treatment
  nw_sums[i,2:14]=edge_metrics(tmp_r)
  
  #save graph
  nw_list[[i]]=graph.adjacency(tmp_r,weight=T,mode="undirected")
  
  i=i+1
}

colnames(node_sums)=c("treatment","degree","degree_norm","betweenness")
node_sums$asv<-rownames(node_sums)
node_sums$kingdom<-substr(node_sums$asv,1,1)

colnames(nw_sums)=c("treatment","num_edges","num_nodes","corr_pos_edges","corr_neg_edges","bac_edges","bac_pos_edges","bac_neg_edges","fun_edges","fun_pos_edges","fun_neg_edges","mean_bac_degree","mean_fun_degree","clustering_coef")

#calculate bacteria-fungi connections
nw_sums$BF_edges<-nw_sums$num_edges-nw_sums$bac_edges-nw_sums$fun_edges
nw_sums$BF_pos_edges<-nw_sums$corr_pos_edges-nw_sums$bac_pos_edges-nw_sums$fun_pos_edges
nw_sums$BF_neg_edges<-nw_sums$corr_neg_edges-nw_sums$bac_neg_edges-nw_sums$fun_neg_edges


