# Scripts & Data for studying Parallel Adaptation through RNAseq and Soil & Leaf Chemistry data in a Phylogenetic framework

## Differential Gene Expression
This is pairwise differential gene expression analyses in three pairs of closely related _Diospyros_ species growing in a common garden. Each pair consists of one species native to normal soil and one to ultramafic soil.
Although the filtration for genes with low counts (R package edgeR) is done on the entire dataset of six species, the differential expression (R package DESeq2) and GO enrichment (topGO) analysis here is shown for one species pair rev-lab (_D.revolutissima_ vs. _D.labillardierei_) only. Z-scores are also calculated for the enriched GO terms. Of course, this can be easily replicated for the two remaining pairs.
