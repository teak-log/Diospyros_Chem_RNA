# The raw RNAseq bam file was preprocessed using featureCounts to generate the counts file that is used as input for the following R code:
# featureCounts" "-p" "--countReadPairs" "-t" "exon" "-g" "gene_id" "-a" "braker.gtf" "--tmpDir" "/tmp/slurm-8220931" "-T" "4" "-o" "RNAcounts.txt" "RNAinput.bam" >

library(DESeq2)
library(edgeR)
library(ggplot2)
library(topGO)
library(data.table)
library(dplyr)
library(ggtext)


### Load raw counts, species assignments and gene annotations ###
setwd("C:/Users/Admin/Documents/rna/root/")
counts_raw <- read.csv("C:/Users/Admin/Documents/rna/root/rootRNAcounts_ordered", header = TRUE, sep = "\t")
counts <- DGEList(counts = counts_raw[,c(1,7:42)], group =
                    as.factor(rep(c("sed","ult","sed","ult","sed","ult"), rep.int(6,6))))  # native soil-type assignments
species <- as.factor(rep(c("cal","spn","imp","heq","lab","rev"), rep.int(6,6) ))  # species assignments
counts$samples$species <- species
group <- counts$samples$group  # counts dataframe specifying species and soil preferences

annotation <- read.csv("C:/Users/Admin/Documents/rna/amin_redo/vie1167c_OmicsboxAminGene_Fannotation.txt", sep = "\t")
annotation <- annotation[annotation$SeqName %in% counts$genes$Geneid,]
annotation <- annotation[,c("SeqName", "Description","Length","GO.Names","GO.IDs")]  #annotations of genes in counts table


### Filter & normalise ###
dim(counts)  #check initial number of genes
keep.exprs <- filterByExpr(counts, group=species)
counts <- counts[keep.exprs,, keep.lib.sizes=FALSE]  # 10/(median lib size) CPM cutoff met for at least 6 samples (size of a species group)
dim(counts)  #check final number of genes post-filtering


### Set up comparisons and find DEGs ###
dds <- DESeqDataSetFromMatrix(
  countData = counts[["counts"]],
  colData = counts[["samples"]],
  design = ~ species
)

dds <- DESeq(dds)

res_revlab <- results(dds, contrast = c("species", "rev", "lab"))
res_revlab_annot <- merge(
  as.data.frame(res_revlab),
  annotation,
  by.x = "row.names",
  by.y = "row.names",
  all.x = TRUE
)

sig_revlab <- subset(res_revlab_annot, padj<0.05 & abs(log2FoldChange) > 1.5)
write.table(sig_revlab[,c(8,9,11,12,3,6,7,4,5,2,10)], file = "revlab_root_deseq_Ppt05_lfc1pt5.genes", sep = "\t", quote = FALSE, row.names = FALSE)


### Setting up GO and gene universe ###

args <- c("C:/Users/Admin/Documents/rna/amin_redo/topgo_amin.ann","revlab_root_deseq_Ppt05_lfc1pt5.genes")  # Annotation and significant genes list
universeFile <- args[1]
interestingGenesFile <- args[2]

# Read gene2GO mapping and gene universe
geneID2GO <- readMappings(file = universeFile)
geneUniverse <- intersect(names(geneID2GO), counts[["genes"]]$Geneid)
geneID2GO <- geneID2GO[geneUniverse]

genesOfInterest <- read.table(interestingGenesFile, sep = "\t", header = TRUE, quote = "")
genesOfInterest <- as.character(genesOfInterest$SeqName)

geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

setDT(sig_revlab)  # DESeq results with SeqName and log2FoldChange


### Function to run topGO and compute z-scores for one ontology ###
run_topgo_zscore <- function(ontology) {
  
  message("Running topGO for ontology: ", ontology)
  
  myGOdata <- new("topGOdata", description = paste("topGO", ontology),
                  ontology = ontology, allGenes = geneList,
                  annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  resultTopgo <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
  
  numsignif <- sum(attributes(resultTopgo)$score <= 0.05)  # Number of significant GO terms (p ≤ 0.05)
  
  GOresults <- GenTable(myGOdata, topgoFisher = resultTopgo,
                        orderBy = "topgoFisher", topNodes = numsignif)  # Get topGO table
  GOresults$logP <- -log10(as.numeric(GOresults$topgoFisher))
  
  # Calculate z-scores using topGO gene membership
  sigGenes <- names(geneList)[geneList == 1]
  go2genes <- genesInTerm(myGOdata)
  
  zscore_per_GO <- lapply(names(go2genes), function(go) {
    genes_in_go <- intersect(go2genes[[go]], sigGenes)
    logfc <- sig_revlab$log2FoldChange[match(genes_in_go, sig_revlab$SeqName)]
    n <- length(logfc)
    if (n == 0) return(NULL)
    
    data.frame(
      GO.ID = go,
      n = n,
      mean_logFC = mean(logfc, na.rm = TRUE),
      z_score = mean(logfc, na.rm = TRUE) * sqrt(n)
    )
  }) %>% bind_rows()
  
  GOresults <- left_join(GOresults, zscore_per_GO, by = c("GO.ID" = "GO.ID"))  # Join results
  GOresults$Ontology <- ontology  # Add ontology column for clarity
  return(GOresults)
}


### Running GO enrichment for all three ontologies, combining results and plotting bar charts ###
results_BP <- run_topgo_zscore("BP")  # Biological Processes
results_CC <- run_topgo_zscore("CC")  # Cellular Component
results_MF <- run_topgo_zscore("MF")  # Molecular Function

all_results <- bind_rows(results_BP, results_CC, results_MF)  # combine ontologies

write.table(all_results, file = "revlab_root_deseq_padj05_lfc1pt5_zscore.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)  # Save combined results

# Plot bar charts where length of bar = -log(pval), continuous colour scales denote z-score,
# and figures in the bars are of the form no. of DE genes / no. of tested genes in the GO term
GOresults <- read.table("revlab_root_deseq_padj05_lfc1pt5_zscore.txt", sep = "\t", header = TRUE)
GOresults %>% ggplot(aes(x = reorder(Term, z_score), y = logP, fill = z_score)) + geom_col() + geom_text(aes(label = paste(Significant, Annotated, sep = "/")),
                                                                                                         hjust = 1, size = 3) + coord_flip() +
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "grey", name = "z-score") +
  labs(x = "", y = expression(-log[10](p-value)),
       title = "Significantly enriched GO terms from root DEGs recovered by deseq2",
       subtitle = "Comparison: *D. revolutissima* vs *D. labillardierei*") + theme_minimal() + theme(plot.subtitle = element_markdown())
