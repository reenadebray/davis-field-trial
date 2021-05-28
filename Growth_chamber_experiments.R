## Analysis of bacterial and fungal communities in growth chamber experiments
## Reena Debray
## April 2021

# Load packages
library(phyloseq)
library(vegan)

## Changes in microbiome composition from the field to the growth chamber
field_G<-subset_samples(ag_ps_nocontam_allsamp,Inoculum=="G sample" | Inoculum=="G inoculum")
field_P<-subset_samples(ag_ps_nocontam_allsamp,Inoculum=="P sample" | Inoculum=="P inoculum")
field_GP<-subset_samples(ag_ps_nocontam_allsamp,Inoculum=="G sample" | Inoculum=="G inoculum" | Inoculum=="P sample" | Inoculum=="P inoculum")
ag_ps_vegan <- veganotu(field_GP)
ag_ps_BC <- vegdist(wisconsin(sqrt(ag_ps_vegan)), method = "bray")

## Distribution of shared and growth-chamber-unique ASVs

# Genotype experiment
genotype<-subset_samples(ag_ps_nocontam_allsamp,Group=="Genotype experiment")
nosamp_taxa<-taxa_names(genotype)[apply(genotype@otu_table,1,function(x){mean(x)==0})]
genotype<-prune_taxa(setdiff(taxa_names(genotype),nosamp_taxa),genotype)
genotype_relabund<-transform_sample_counts(genotype, function(x){x/sum(x)})

GC_unique_ASVs<-c(setdiff(taxa_names(genotype),taxa_names(ag_ps_rmnosamp)))
GC_field_shared<-c(intersect(taxa_names(genotype),taxa_names(ag_ps_rmnosamp)))
G_prop_shared<-data.frame(matrix(nrow=0,ncol=3))
colnames(G_prop_shared)=c("Genotype","Growth.Chamber.Genotype","prop_shared")

#iterate through samples in genotype experiment
for (i in seq(1,nsamples(genotype))){
  sample<-colnames(genotype_relabund@otu_table)[i]
  
  #total relative abundance (from 0 to 1) composed of ASVs shared with field
  G_prop_shared[i,c("sample","Genotype","Growth.Chamber.Genotype")]=c(sample,genotype@sam_data[sample,c("Genotype","Growth.Chamber.Genotype")])
  G_prop_shared[i,"prop_shared"]=sum(genotype_relabund@otu_table[rownames(genotype_relabund@otu_table)%in%GC_field_shared,sample])
}

# Phosphorus experiment
phosphorus<-subset_samples(ag_ps_nocontam_allsamp,Group=="Phosphorous experiment")
nosamp_taxa<-taxa_names(phosphorus)[apply(phosphorus@otu_table,1,function(x){mean(x)==0})]
phosphorus<-prune_taxa(setdiff(taxa_names(phosphorus),nosamp_taxa),phosphorus)
phosphorus_relabund<-transform_sample_counts(phosphorus, function(x){x/sum(x)})

GC_unique_ASVs_P<-setdiff(taxa_names(phosphorus),taxa_names(ag_ps_rmnosamp))
GC_field_shared_P<-intersect(taxa_names(phosphorus),taxa_names(ag_ps_rmnosamp))
P_prop_shared<-data.frame(matrix(nrow=0,ncol=3))
colnames(P_prop_shared)=c("Phosphorous.Level","Growth.Chamber.Phosphorous.Level","prop_shared")

#iterate through samples in genotype experiment
for (i in seq(1,nsamples(phosphorus))){
  sample<-colnames(phosphorus_relabund@otu_table)[i]
  
  #total relative abundance (from 0 to 1) composed of ASVs shared with field
  a<-sum(phosphorus_relabund@otu_table[rownames(phosphorus_relabund@otu_table)%in%GC_field_shared,sample])
  P_prop_shared[i,c("sample","Phosphorous.Level","Growth.Chamber.Phosphorous.Level")]=c(sample,phosphorus@sam_data[sample,c("Phosphorous.Level","Growth.Chamber.Phosphorous.Level")])
  P_prop_shared[i,"prop_shared"]=a
}
P_prop_shared[is.na(P_prop_shared$Phosphorous.Level),"Phosphorous.Level"]="buffer"
#plot shared ASVs
P_prop_shared$Phosphorous.Level<-factor(P_prop_shared$Phosphorous.Level,levels=c("3","1","buffer"))

## Remove growth-chamber-unique ASVs and calculate alpha and beta diversity
# Genotype experiment
genotype_noGCU<-prune_taxa(GC_field_shared,genotype)
# alpha diversity
genotype_noGCU@sam_data$observed<-apply(genotype_noGCU@otu_table,2,function(x){length(x[x>0])})
# beta diversity
dist_bray_GC<-as.matrix(phyloseq::distance(genotype_noGCU,method="bray"))
WT_names<-sample_names(subset_samples(genotype_noGCU, Growth.Chamber.Genotype=="76R"))
rmc_names<-sample_names(subset_samples(genotype_noGCU,Growth.Chamber.Genotype=="rmc"))
GC_dispersion_analysis<-get_pairwise_dist(dist_matrix = dist_bray_GC,grouped_list_of_names = list(WT_names,rmc_names),treatments = c("76R","rmc"))

# Phosphorus experiment
phosphorus_noGCU<-prune_taxa(GC_field_shared_P,phosphorus)
# alpha diversity
phosphorus_noGCU@sam_data$observed<-apply(phosphorus_noGCU@otu_table,2,function(x){length(x[x>0])})
# beta diversity
dist_bray_GCP<-as.matrix(phyloseq::distance(phosphorus_noGCU,method="bray"))
lowP_names<-sample_names(subset_samples(phosphorus_noGCU, Growth.Chamber.Phosphorous.Level=="low"))
highP_names<-sample_names(subset_samples(phosphorus_noGCU,Growth.Chamber.Phosphorous.Level=="high"))
GC_dispersion_analysis_P<-get_pairwise_dist(dist_matrix = dist_bray_GCP,grouped_list_of_names = list(lowP_names,highP_names),treatments = c("low P","high P"))

##Bray-Curtis distances from inoculum to sample
all_plants<-unique(ag_ps_nocontam_allsamp@sam_data$Plant.ID)
g_plants<-all_plants[substr(all_plants,2,2)==8 | substr(all_plants,2,2)==9]
p_plants<-all_plants[substr(all_plants,2,2)!=8 & substr(all_plants,2,2)!=9 & all_plants!=""]
  
# genotype experiment
genotype_exp_beta<-data.frame(matrix(nrow=0,ncol=6))
for (plant in g_plants){
    # for each plant in the field, identify the two plants in the growth chamber that received its inocula
    trio<-subset_samples(ag_ps_nocontam_allsamp@sam_data,Plant.ID==plant)
    inoculum=rownames(trio[trio$Inoculum=="G inoculum",])
    rmc=rownames(trio[trio$Growth.Chamber.Genotype=="rmc" & !is.na(trio$Growth.Chamber.Genotype),])
    WT=rownames(trio[trio$Growth.Chamber.Genotype=="76R" & !is.na(trio$Growth.Chamber.Genotype),])
    
    # calculate dissimilarity between inoculum and rmc sample
    inoculum_rmc<-field_GC_bray[inoculum,rmc]
    plantID<-as.character(trio[inoculum,"Treatment.ID"])
    type<-paste(trio[inoculum,"Genotype"],"rmc",sep="/")
    genotype_exp_beta<-rbind(genotype_exp_beta,c(plantID,type,inoculum_rmc))
    
    # calculate dissimilarity between inoculum and wild-type sample
    inoculum_WT<-field_GC_bray[inoculum,WT]
    plantID<-as.character(trio[inoculum,"Treatment.ID"])
    type<-paste(trio[inoculum,"Genotype"],"76R",sep="/")
    genotype_exp_beta<-rbind(genotype_exp_beta,c(plantID,type,inoculum_WT))
    
    # calculate dissimilarity between rmc sample and wild-type sample
    rmc_WT<-field_GC_bray[rmc,WT]
    samples<-paste(trio[rmc,"Treatment.ID"],trio[WT,"Treatment.ID"])
    genotype_exp_beta<-rbind(genotype_exp_beta,c(samples,"samples",rmc_WT))
}
colnames(genotype_exp_beta)=c("plantID","type","beta_bray")
genotype_exp_beta$beta_bray<-as.numeric(genotype_exp_beta$beta_bray)
genotype_exp_beta[genotype_exp_beta$type=="76R/76R" | genotype_exp_beta$type=="rmc/rmc","class"]="same_genotype"
genotype_exp_beta[genotype_exp_beta$type=="76R/rmc" | genotype_exp_beta$type=="rmc/76R","class"]="different_genotype"
genotype_exp_beta[genotype_exp_beta$type=="samples","class"]="samples"
genotype_exp_beta$inoculum_genotype<-substr(genotype_exp_beta$type,1,3)
genotype_exp_beta$type<-factor(genotype_exp_beta$type,levels=c("76R/76R","rmc/rmc","76R/rmc","rmc/76R"))
genotype_exp_beta$class<-factor(genotype_exp_beta$class,levels=c("same_genotype","different_genotype","samples"))
wilcox.test(beta_bray~class,paired=T,data=genotype_exp_beta[genotype_exp_beta$class!="samples",],alternative="less")

# phosphorous experiment
phosphorous_exp_beta<-data.frame(matrix(nrow=0,ncol=6))
for (plant in p_plants){
    # for each plant in the field, identify the two plants in the growth chamber that received its inocula
    trio<-subset_samples(ag_ps_nocontam_allsamp@sam_data,Plant.ID==plant)
    inoculum=rownames(trio[trio$Inoculum=="P inoculum",])
    low=rownames(trio[trio$Growth.Chamber.Phosphorous.Level=="low" & !is.na(trio$Growth.Chamber.Phosphorous.Level),])
    high=rownames(trio[trio$Growth.Chamber.Phosphorous.Level=="high" & !is.na(trio$Growth.Chamber.Phosphorous.Level),])
    
    #calculate dissimilarity between inoculum and low-P sample
    inoculum_low<-field_GC_bray[inoculum,low]
    plantID<-as.character(trio[inoculum,"Treatment.ID"])
    type<-paste(trio[inoculum,"Phosphorous.Level"],"low",sep="/")
    phosphorous_exp_beta<-rbind(phosphorous_exp_beta,c(plantID,type,inoculum_low))
    
    #calculate dissimilarity between inoculum and high-P sample
    inoculum_high<-field_GC_bray[inoculum,high]
    plantID<-as.character(trio[inoculum,"Treatment.ID"])
    type<-paste(trio[inoculum,"Phosphorous.Level"],"high",sep="/")
    phosphorous_exp_beta<-rbind(phosphorous_exp_beta,c(plantID,type,inoculum_high))
    
    #calculate dissimilarity between high-P sample and low-P sample
    low_high<-field_GC_bray[low,high]
    samples<-paste(trio[low,"Treatment.ID"],trio[high,"Treatment.ID"])
    phosphorous_exp_beta<-rbind(phosphorous_exp_beta,c(samples,"samples",low_high))
    
}
colnames(phosphorous_exp_beta)=c("plantID","type","beta_bray")
phosphorous_exp_beta$beta_bray<-as.numeric(phosphorous_exp_beta$beta_bray)
phosphorous_exp_beta[phosphorous_exp_beta$type=="3/high" | phosphorous_exp_beta$type=="1/low","class"]="same_phosphorous"
phosphorous_exp_beta[phosphorous_exp_beta$type=="3/low" | phosphorous_exp_beta$type=="1/high","class"]="different_phosphorous"
phosphorous_exp_beta[phosphorous_exp_beta$type=="samples","class"]="samples"
phosphorous_exp_beta$type<-factor(phosphorous_exp_beta$type,levels=c("1/low","3/high","1/high","3/low"))
phosphorous_exp_beta$class<-factor(phosphorous_exp_beta$class,levels=c("same_phosphorous","different_phosphorous","samples"))
wilcox.test(beta_bray~class,paired=T,data=phosphorous_exp_beta[phosphorous_exp_beta$class!="samples",],alternative="less")
