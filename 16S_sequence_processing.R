## Implementation of DADA2 pipeline on 16S sequences to create phyloseq object for data analysis
## Reena Debray
## April 2021

# Load packages
library(dada2)
library(phyloseq)
library(decontam)

# Initialize "path" variable as directory containing fastq files
ag_fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
ag_fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
ag.sample.names <- sapply(strsplit(basename(ag_fnFs), "_"), `[`, 1)

# Filter reads with Ns or >2 expected errora based on quality scores
# Truncate reads at 240 bp in the forward direction and 140 bp in the reverse direction, or at first position with a quality score <2
ag_filtFs <- file.path(path, "filtered/", paste0(ag.sample.names, "_F_filt.fastq.gz"))
ag_filtRs <- file.path(path, "filtered/", paste0(ag.sample.names, "_R_filt.fastq.gz"))
names(ag_filtFs) <- ag.sample.names
names(ag_filtRs) <- ag.sample.names
ag_out <- filterAndTrim(ag_fnFs, ag_filtFs, ag_fnRs, ag_filtRs, truncLen=c(240,140), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread = TRUE)

# Train core denoising algorithm on error rates
ag_errF <- learnErrors(ag_filtFs, multithread=TRUE)
ag_errR <- learnErrors(ag_filtRs, multithread=TRUE)

# Infer amplicon sequencing variants
ag_dadaFs <- dada(ag_filtFs, err=ag_errF, multithread=TRUE)
ag_dadaRs <- dada(ag_filtRs, err=ag_errR, multithread=TRUE)

# Merge paired reads and remove chimeras
ag_mergers <- mergePairs(ag_dadaFs, ag_filtFs, ag_dadaRs, ag_filtRs, verbose=TRUE)
ag_seqtab <- makeSequenceTable(ag_mergers)
ag_seqtab.nochim <- removeBimeraDenovo(ag_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Assign taxonomy
# We used the Silva taxonomic training dataset located at https://zenodo.org/record/1172783#.YH9ESRNKiL8
# Initialize "training_set" variable as the path to the fasta file
ag_taxa <- assignTaxonomy(ag_seqtab.nochim, training_set, multithread=TRUE)
ag_taxa.print <- ag_taxa 
rownames(ag_taxa.print) <- NULL

# Make ASV table
ag_asv_seqs <- colnames(ag_seqtab.nochim)
 ag_asv_headers <- vector(dim(ag_seqtab.nochim)[2], mode="character")
 for (i in 1:dim(ag_seqtab.nochim)[2]) {
 ag_asv_headers[i] <- paste(">ASV", i, sep="_")
 }
ag_asv_tab <- t(ag_seqtab.nochim)

# Create phyloseq item (using a metadata file called "ag_sample_info")
row.names(ag_asv_tab) <- sub(">", "", ag_asv_headers)
row.names(ag_taxa.print) <- sub(">", "", ag_asv_headers)
ag_OTU = otu_table(ag_asv_tab, taxa_are_rows = TRUE)
ag_TAX = tax_table(ag_taxa)
ag_samples = ag_sample_info
ag_ps <-phyloseq(ag_OTU, ag_TAX, ag_samples)

# Remove chloroplast and mitochondrial sequences & PCR positive control
ag_ps_nochloro <- subset_taxa(ag_ps, (Order!="Chloroplast" | is.na(Order))) #important to specify is.na() here, or phyloseq will remove NA taxa as well
ag_ps_NoCNoM <- subset_taxa(ag_ps_nochloro, (Family!="Mitochondria" | is.na(Family))) #important to specify is.na() here, or phyloseq will remove NA taxa as well
ag_ps_NoCNoM_nocontrols <- subset_samples(ag_ps_NoCNoM, Group=="Field" | Group=="Genotype experiment" | Group=="Phosphorous experiment")
ag_ps_rmpos = subset_samples(ag_ps_NoCNoM, Group != "positive control") #includes sequencing controls, extraction controls, buffer, & fertilizer
ag_samples_df <- as.data.frame(sample_data(ag_ps_rmpos))

# Run decontam function on PCR negative controls
ag_ps_rmpos_negdecontam <- subset_samples(ag_ps_rmpos, (Group!="Extraction control" & Treatment.ID!="(sampling buffer)") | is.na(Treatment.ID))
sample_data(ag_ps_rmpos_negdecontam)$is.neg <- (sample_data(ag_ps_rmpos_negdecontam)$Group == "negative control")
contamdf.prev1 <- isContaminant(ag_ps_rmpos_negdecontam, method="prevalence", neg="is.neg", threshold=0.1)
contam_1<-which(contamdf.prev1$contaminant)

# Run decontam on buffer negative control
ag_ps_rmpos_buffer<- subset_samples(ag_ps_rmpos, (Group!="Extraction control" & !is.na(Treatment.ID)))
sample_data(ag_ps_rmpos_buffer)$is.neg <- (sample_data(ag_ps_rmpos_buffer)$Treatment.ID == "(sampling buffer)")
contamdf.prev2 <- isContaminant(ag_ps_rmpos_buffer, method="prevalence", neg="is.neg", threshold=0.1)
contam_2<-which(contamdf.prev2$contaminant)

# Run decontam on DNA extraction negative controls, blocked by extraction batch
ag_ps_rmpos_extraction<- subset_samples(ag_ps_rmpos, (Treatment.ID != "(sampling buffer)" & !is.na(Treatment.ID)))
sample_data(ag_ps_rmpos_extraction)$is.neg <- (sample_data(ag_ps_rmpos_extraction)$Group == "Extraction control")
sample_data(ag_ps_rmpos_extraction)$batch<-c(rep("B1",24),rep("B2",24),rep("B3",23),rep("B4",24),rep("B5",24),rep("B6",19),rep("B7",24))
contamdf.prev3 <- isContaminant(ag_ps_rmpos_extraction, method="prevalence",batch="batch", neg="is.neg", threshold=0.1)
contam_3<-which(contamdf.prev3$contaminant)

# Merge contaminant lists from different negative controls and filter from data
ag_badTaxa = taxa_names(ag_ps_rmpos)[sort(unique(c(contam_1,contam_2,contam_3)))]
ag_goodTaxa <- setdiff(taxa_names(ag_ps_NoCNoM_nocontrols), ag_badTaxa)
ag_ps_nocontam_allsamp = prune_taxa(ag_goodTaxa, ag_ps_NoCNoM_nocontrols)
