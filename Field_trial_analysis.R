## Analysis of bacterial and fungal communities in field trial
## Reena Debray
## April 2021

# Load packages
library(phyloseq)
library(ggplot2)
library(gridExtra)

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

### Figure 1: effects of drought on bacterial and fungal diversity

# Figure 1A: drought & bacterial richness 
ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Water.Regiment=="drought","water2"]="50% deficit"
ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Water.Regiment=="full water","water2"]="Full water"
ag_ps_rmnosamp@sam_data$water2<-factor(ag_ps_rmnosamp@sam_data$water2,levels=c("Full water","50% deficit"))
ggplot(ag_ps_rmnosamp@sam_data,aes(water2,observed))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Observed bacterial richness")

# Figure 1B: drought & fungal richness
ps_field_noNA@sam_data[ps_field_noNA@sam_data$Water.Regiment=="drought","water2"]="50% deficit"
ps_field_noNA@sam_data[ps_field_noNA@sam_data$Water.Regiment=="full water","water2"]="Full water"
ps_field_noNA@sam_data$water2<-factor(ps_field_noNA@sam_data$water2,levels=c("Full water","50% deficit"))
ggplot(ps_field_noNA@sam_data,aes(water2,observed))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Observed fungal richness")

# Figure 1C: drought & bacterial shannon diversity
ggplot(ag_ps_rmnosamp@sam_data,aes(water2,shannon))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Bacterial diversity (Shannon index)")

# Figure 1D: drought & fungal shannon diversity
ggplot(ps_field_noNA@sam_data,aes(water2,shannon))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Bacterial diversity (Shannon index)")

# Figure 1E: bacterial beta diversity within drought treatments
drought_dispersion_analysis[drought_dispersion_analysis$Water.Regiment=="drought","water2"]="50% deficit"
drought_dispersion_analysis[drought_dispersion_analysis$Water.Regiment=="full water","water2"]="Full water"
drought_dispersion_analysis$water2<-factor(drought_dispersion_analysis$water2,levels=c("Full water","50% deficit"))
ggplot(drought_dispersion_analysis,aes(water2,avg_BC))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Bacterial beta diversity within group")

# Figure 1F: fungal beta diversity within drought treatments
drought_dispersion_ITS_analysis[drought_dispersion_ITS_analysis$Water.Regiment=="drought","water2"]="50% deficit"
drought_dispersion_ITS_analysis[drought_dispersion_ITS_analysis$Water.Regiment=="full water","water2"]="Full water"
drought_dispersion_ITS_analysis$water2<-factor(drought_dispersion_ITS_analysis$water2,levels=c("Full water","50% deficit"))
ggplot(drought_dispersion_ITS_analysis,aes(water2,avg_BC))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Fungal beta diversity within group")

### Figure 2: effects of plant genotype on bacterial and fungal diversity

# Figure 2A: genotype & bacterial richness
ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Genotype=="76R","genotype2"]="Wild-type 76R"
ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Genotype=="rmc","genotype2"]="Reduced\nmycorrhizal"
ag_ps_rmnosamp@sam_data$genotype2<-factor(ag_ps_rmnosamp@sam_data$genotype2,levels=c("Wild-type 76R","Reduced\nmycorrhizal"))
ggplot(ag_ps_rmnosamp@sam_data,aes(genotype2,observed))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Observed bacterial richness")

# Figure 2B: genotype & fungal richness
ps_field_noNA@sam_data[ps_field_noNA@sam_data$Genotype=="76R","genotype2"]="Wild-type 76R"
ps_field_noNA@sam_data[ps_field_noNA@sam_data$Genotype=="rmc","genotype2"]="Reduced\nmycorrhizal"
ps_field_noNA@sam_data$genotype2<-factor(ps_field_noNA@sam_data$genotype2,levels=c("Wild-type 76R","Reduced\nmycorrhizal"))
ggplot(ps_field_noNA@sam_data,aes(genotype2,observed))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Observed fungal richness")

# Figure 2C: genotype & bacterial shannon diversity
ggplot(ag_ps_rmnosamp,aes(genotype2,shannon))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Bacterial diversity (Shannon index)")

# Figure 2D: genotype & fungal shannon diversity
ggplot(ps_field_noNA,aes(genotype2,shannon))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Fungal diversity (Shannon index)")

# Figure 2E: bacterial beta diversity within genotypes
genotype_dispersion_analysis[genotype_dispersion_analysis$Genotype=="76R","genotype2"]="Wild-type 76R"
genotype_dispersion_analysis[genotype_dispersion_analysis$Genotype=="rmc","genotype2"]="Reduced\nmycorrhizal"
genotype_dispersion_analysis$genotype2<-factor(genotype_dispersion_analysis$genotype2,levels=c("Wild-type 76R","Reduced\nmycorrhizal"))
ggplot(genotype_dispersion_analysis,aes(genotype2,avg_BC))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Bacterial beta diversity within group")

# Figure 2F: fungal beta diversity within genotypes
genotype_dispersion_ITS_analysis[genotype_dispersion_ITS_analysis$Genotype=="76R","genotype2"]="Wild-type 76R"
genotype_dispersion_ITS_analysis[genotype_dispersion_ITS_analysis$Genotype=="rmc","genotype2"]="Reduced\nmycorrhizal"
genotype_dispersion_ITS_analysis$genotype2<-factor(genotype_dispersion_ITS_analysis$genotype2,levels=c("Wild-type 76R","Reduced\nmycorrhizal"))
ggplot(genotype_dispersion_ITS_analysis,aes(genotype2,avg_BC))+geom_boxplot(width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(width=0.05,size=2,alpha=0.8)+theme_classic(base_size=24)+guides(color=F)+xlab("")+ylab("Fungal beta diversity within group")
