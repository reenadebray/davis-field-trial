## Implementation of DADA2 pipeline on ITS sequences to create phyloseq object for data analysis
## Reena Debray
## April 2021

# Load packages
library(dada2)
library(ShortRead)
library(Biostrings)
library(phyloseq)
library(stringr)

# Initialize "path" variable as directory containing fastq files
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Filter reads shorter than 50 bp
# Truncate reads at the first position with a quality score <2
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

# Train core denoising algorithm on error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Dereplicate idnetical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer amplicon sequencing variants
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge paired reads and remove chimeras
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Assign taxonomy
# We used the General Release Fasta at this link: https://unite.ut.ee/repository.php
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
taxa.print <- taxa  #Removing sequence rownames for display only
rownames(taxa.print) <- NULL

# Make ASV table
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_tab <- t(seqtab.nochim)

# Create phyloseq item (using a metadata file called "ag_sample_info")
row.names(asv_tab) <- sub(">", "", asv_headers)
row.names(taxa.print) <- sub(">", "", asv_headers)
OTU = otu_table(asv_tab, taxa_are_rows = TRUE)
TAX = tax_table(taxa.print)
samples = sample_info_ps
ps <-phyloseq(OTU, TAX, samples)

# Identify ASVs that are unclassified at the phylum level
ps_phylumNA <- subset_taxa(ps, is.na(Phylum))
NA_ASV_names <- rownames(tax_table(ps_phylumNA))
NA_ASV_numbers <- as.numeric(str_remove(NA_ASV_names, "ASV_"))
NA_ASV_numbers <- as.vector(NA_ASV_numbers)

# Identify non-fungal eukaryotic sequences using the program ITSx (see instructions in ITSx_sequence_filtering file)
# load ITSx output as "ITSx_output"
# Remove ASVs unclassified at the phylum level that ITSx did not identify as fungi
ITSx_fungal <- subset(ITSx_output, V3 == "F")
ITSx_fungal_numbers <- as.numeric(ITSx_fungal$V1)
nonfungalASVs_ITSx <- NA_ASV_numbers[! NA_ASV_numbers %in% ITSx_fungal_numbers]
nonfungalASVs_ITSx <- paste("ASV_", nonfungalASVs_ITSx, sep="")
fungalTaxa <- setdiff(taxa_names(ps), nonfungalASVs_ITSx)
ps_noNA_ITSx = prune_taxa(fungalTaxa, ps)
