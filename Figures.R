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

### Figure 3: Parallel effects of water regime and genotype

# Figure 3A: alpha diversity across factorial water and genotype combinations
ag_ps_rmnosamp@sam_data$Water_Genotype<-paste(ag_ps_rmnosamp@sam_data$Water.Regiment,ag_ps_rmnosamp@sam_data$Genotype)
ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Water_Genotype=="full water rmc","water_genotype_2"]="Reduced\nmycorrhizal \n Full water\n"
ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Water_Genotype=="drought 76R","water_genotype_2"]="Wild-type 76R \n Water deficit\n"
ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Water_Genotype=="drought rmc","water_genotype_2"]="Reduced\nmycorrhizal \n Water deficit\n"
ag_ps_rmnosamp@sam_data[ag_ps_rmnosamp@sam_data$Water_Genotype=="full water 76R","water_genotype_2"]="Wild-type 76R \n Full water\n"
ag_ps_rmnosamp@sam_data$water_genotype_2<-factor(ag_ps_rmnosamp@sam_data$water_genotype_2,levels=c("Wild-type 76R \n Full water\n","Wild-type 76R \n Water deficit\n","Reduced\nmycorrhizal \n Full water\n","Reduced\nmycorrhizal \n Water deficit\n"))
ggplot(ag_ps_rmnosamp@sam_data)+geom_boxplot(aes(water_genotype_2,observed,color=water_genotype_2),width=0.5,size=0.8,outlier.shape=NA)+geom_jitter(aes(water_genotype_2,observed,fill=water_genotype_2),width=0.05,size=3,shape=21)+theme_classic(base_size=21)+scale_fill_manual(values=rev(viridis(4)))+guides(color=F)+scale_color_manual(values=rev(viridis(4)))+guides(fill=F)+xlab("")+ylab("Observed bacterial richness")

# Figure 3B: NMDS ordination across factorial water and genotype combinations
set.seed(123)
plot_ordination(ag_ps_rmnosamp,ordinate(ag_ps_rmnosamp,"NMDS",distance = "bray"),type="samples")+geom_point(size=3,shape=21,aes(fill=water_genotype_2))+theme_classic(base_size=22)+stat_ellipse(aes(color=water_genotype_2),size=1)+scale_color_manual(values=rev(viridis(4)))+scale_fill_manual(values=rev(viridis(4)))+guides(fill=F)

# Figure 3C: indicator species analysis
# Add a "stat" variable with the magnitude of the indicator value, and positive or negative sign depending on the treatment
# In the water treatment, index==1 is drought, index==2 is full water, and index==3 is neither (either they didn't occur in either sample or they were too evenly distributed)
# In the genotype treatment, index==1 is 76R, index==2 is rmc, and index==3 is neither
indic_asvs[indic_asvs$treatment=="water" & indic_asvs$index==1,"stat"]=-1*indic_asvs[indic_asvs$treatment=="water" & indic_asvs$index==1,"stat"]
indic_asvs[indic_asvs$treatment=="water" & indic_asvs$index==3,"stat"]=0
indic_asvs[indic_asvs$treatment=="genotype" & indic_asvs$index==2,"stat"]=-1*indic_asvs[indic_asvs$treatment=="genotype" & indic_asvs$index==2,"stat"]
indic_asvs[indic_asvs$treatment=="genotype" & indic_asvs$index==3,"stat"]=0
indic_asvs$multiplier<-indic_asvs$stat/abs(indic_asvs$stat)
indic_asvs[is.na(indic_asvs$multiplier),"multiplier"]=0
indic_asvs[indic_asvs$pval<0.05 & !is.na(indic_asvs$pval),"label"]="*"
indic_asvs[indic_asvs$pval>=0.05 | is.na(indic_asvs$pval),"label"]=""

# Add family names
for (i in seq(1,nrow(indic_asvs))){asv=indic_asvs[i,"ASV"];indic_asvs[i,"Family"]=as.character(ag_ps_rmnosamp@tax_table[asv,"Family"])}
indic_asvs<-indic_asvs[order(indic_asvs$treatment,indic_asvs$stat),]
indic_asvs$ASV<-factor(indic_asvs$ASV,levels=indic_asvs[indic_asvs$treatment=="water","ASV"])
indic_asvs[is.na(indic_asvs$Family),"Family"]="Unclassified"
indic_asvs[indic_asvs$Family=="Family_XII","Family"]="Bacillales Family XII"
indic_asvs[indic_asvs$Family=="Family_XI","Family"]="Clostridiales Family XI"
indic_asvs$Family_ASV<-paste(indic_asvs$Family," (",indic_asvs$ASV,")",sep="")

# Adjust treatment labels
indic_asvs$treatment<-factor(indic_asvs$treatment,levels=c("water","genotype"))
indic_asvs[indic_asvs$treatment=="water","treatment2"]="Water regime"
indic_asvs[indic_asvs$treatment=="genotype","treatment2"]="Mycorrhizal genotype"
indic_asvs$treatment2<-factor(indic_asvs$treatment2,levels=c("Water regime","Mycorrhizal genotype"))
indic_asvs<-indic_asvs[order(indic_asvs$treatment,indic_asvs$stat),]
indic_asvs$Family_ASV<-factor(indic_asvs$Family_ASV,levels=indic_asvs[indic_asvs$treatment=="water","Family_ASV"])
ggplot(indic_asvs,aes(stat,Family_ASV,fill=(stat<0)))+geom_bar(stat="identity",position="dodge")+theme_classic(base_size=24)+xlab("Indicator value")+ylab("")+guides(fill=F)+scale_fill_manual(values=c("#FCD225FF","#39568CFF"))+facet_wrap(~treatment2)+ theme(panel.spacing = unit(2.5, "lines"),axis.text.y=element_text(size=16),axis.title.x= element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)))+geom_text(aes(stat+0.1*multiplier,Family_ASV,label=label),size=7)

