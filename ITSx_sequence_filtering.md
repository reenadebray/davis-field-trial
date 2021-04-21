## Method for identifying fungal from non-fungal sequences in ITS amplicon datasets
## Catherine A. Hernandez
## January 2021

After processing sequencing results and generating a phyloseq object in R, many ASVs were classified as Kingdom “Fungi” but had no Phylum (or lower) assignment. 
I extracted the sequences of these ASVs and BLASTed them – many matched to plants or algae, but some may be unclassified fungi. 
In order to decide which ASVs to retain/remove, I used a program called ITSx (https://microbiology.se/software/itsx/; https://microbiology.se/publ/itsx_users_guide.pdf). 
This program uses hidden Markov models (HMMs) to extract ITS sequences from the amplicons and assign them to putative origins. 
I then compared the origin classifications from ITSx with the top 5 BLAST hits for each ASV, which matched well. 
Moving forward in the analysis, we will retain only the ASVs that ITSx assigned to fungi. 
Assuming the link hasn’t changed, follow the installation and usage instructions in the manual, but below is how I got it to work.

1. Downloaded “Homebrew” package manager for Mac
2. Downloaded HMMER (http://hmmer.org/download.html) to binaries directory on my computer (/usr/local/bin)
3. Downloaded ITSx package to binaries directory.
4. In Mac terminal, check that ITSx download was successful by running: ITSx --help
5. Run ITSx with the following code:./ITSx -i {path to FASTA file of sequences} -o {path to output file}
6. If you want a detailed table of results, add --detailed_results T to the code in line 5
7. “Origin” in the detailed results table lists the putative classification of thesequence origin  (“F” is fungal)
8. The ITSx manual recommends comparing their classifications to BLAST results. For the dataset I worked with (Reena’s field trial experiment) the results matched well.
