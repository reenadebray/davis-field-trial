## Analysis of bacterial and fungal communities in field trial
## Reena Debray
## April 2021

# Load packages
library(phyloseq)
library(vegan)
library(indicspecies)
library(igraph)
library(Hmisc)

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

# Indicator species analysis for water regime on wild-type background
abund<-t(subset_samples(ag_ps_rmnosamp,Genotype=="76R")@otu_table)
water<-subset_samples(ag_ps_rmnosamp,Genotype=="76R")@sam_data$Water.Regiment
inv_water = as.matrix(multipatt(abund,water,func="IndVal.g",control = how(nperm=999)))[,1]$sign

# Indicator species analysis for genotype on full-water background
abund<-t(subset_samples(ag_ps_rmnosamp,Water.Regiment=="full water")@otu_table)
genotype<-subset_samples(ag_ps_rmnosamp,Water.Regiment=="full water")@sam_data$Genotype
inv_genotype = as.matrix(multipatt(abund,genotype,func="IndVal.g",control = how(nperm=999)))[,1]$sign

# Identify ASVs that were significant hits for either water regime or genotype (or both)
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
