# davis-field-trial

Reena Debray, Yvonne Socolar, Griffin Kaulbach, Aidee Guzman, Catherine A. Hernandez, Rose Curley, Alexander Dhond, Timothy Bowles, and Britt Koskella.

## Abstract
Water and nutrient limitation are key stressors that affect plant health and ecosystem function. Yet, the extent to which these stressors act by altering plant-associated microbial communities remains poorly understood, particularly in aboveground plant tissues. Using experimental manipulations in the field and growth chamber, we probe the effects of irrigation, soil fertility, and mycorrhizal associations on bacterial and fungal communities in the tomato phyllosphere (Solanum lycopersicum). Water stress and mycorrhizal disruption reduced bacterial richness and homogenized communities across plants. Fungal communities also shifted under water stress and mycorrhizal disruption, largely driven by reductions in the relative abundance of dominant fungal taxa. We observed substantial parallelism in the individual microbial taxa affected by irrigation and mycorrhizal associations, suggesting that both conditions induce systemic changes in plant physiology that reshape the aboveground microbiome. Our analyses provide novel insights into the indirect impacts of biotic and abiotic conditions on aboveground microbial communities.

## Funding
Support for this work was provided by the Army Research Office (W911NF-17-1-0231 to AG), the National Science Foundation (NSF Graduate Research Fellowships to RD, CH, AG), the University of California, Berkeley (Berkeley Fellowships to YS, CH), the Society for the Study of Evolution (Grant 047408 to RD), the Hellman Fellows Fund (Hellman Fellows Award to BK) and the Marian E. Koshland Integrated Natural Sciences Center (KINSC Summer Scholarship to GK).

## This repository includes:
### Functions.R
Contains all functions written for the analysis.
### 16S_sequence_processing.R
Implementation of DADA2 pipeline and removal of plant chloroplast and mitochondrial sequences.
### ITS_sequence_processing.R
Implementation of DADA2 ITS pipeline and removal of non-fungal eukaryotic sequences.
### ITS_sequence_filtering.md
Describes how to use ITSx to classify ITS sequences as fungal or non-fungal.
### Field_trial_analysis.R
All code used to calculate diversity and compositional metrics in the field trial data..
### Growth_chamber_experiments.R
All code used to calculate diversity and compositional metrics in the growth chamber experiments.
### Figures.R
Code used to generate all figures in the manuscript.
