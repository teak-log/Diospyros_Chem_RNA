# Scripts & Data for studying Parallel Adaptation through RNAseq and Soil & Leaf Chemistry data in a Phylogenetic framework

1. Soil & Leaf Chemistry PCA
1. Differential Gene Expression

## Soil & Leaf Chemistry PCA
Separate PCAs are made for the leaf and soil chemistry datasets to get a broad overview. Each loading is the concentration of an element normalised across the distribution of the variable such that the mean is zero and variance is one.

## Differential Gene Expression
This is pairwise differential gene expression analyses in three pairs of closely related _Diospyros_ species growing in a common garden. Each pair consists of one species native to normal soil and one to ultramafic soil.
Although the filtration for genes with low counts (R package edgeR) is done on the entire dataset of six species, the differential expression (R package DESeq2) and GO enrichment (topGO) analysis here is shown for one species pair rev-lab (_D.revolutissima_ vs. _D.labillardierei_) only. Z-scores are also calculated for the enriched GO terms. This can be easily replicated for the two remaining pairs.
