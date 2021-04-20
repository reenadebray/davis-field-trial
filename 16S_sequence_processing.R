## Implementation of DADA2 pipeline on 16S sequences to create phyloseq object for data analysis
## Reena Debray
## April 2021

# Load packages
library(dada2)
library(phyloseq)

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
