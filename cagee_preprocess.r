library(ape)
library(phytools)
library(edgeR)
library(ggplot2)
library(ggtree)
library(stringr)
library(dplyr)
library(VennDiagram)
library(gridExtra)

### Generating an ultrametric input tree ###
mytree <- read.tree(file = "astral2_unann.nwk")  # load full phylogeny
plot(mytree)
myprunedtree <- keep.tip(mytree, c("labillardierei", "revolutissima", "spPicNga", "calciphila", "hequetiae", "impolita"))  # selecting species with RNAseq data
plot(myprunedtree)

# plot trees made ultrametric by different methods and choose the one closest to the original tree
par(mfrow=c(4,2), mar=c(0.5,0.5,1,0.5))

plot(myprunedtree, edge.width = 1,
     cex = 0.7,
     main = "A) Non-ultrametric",
     cex.main = 1)

plot(force.ultrametric(myprunedtree, method=c("extend")), edge.width = 1,
     cex = 0.7,
     main = "B) Extension method",
     cex.main = 1)

plot(force.ultrametric(myprunedtree, method=c("nnls")), edge.width = 1,
     cex = 0.7,
     main = "C) NNLS method",
     cex.main = 1)

plot(chronoMPL(myprunedtree), edge.width = 1,
     cex = 0.7,
     main = "D) Ultrametric chronoMPL",
     cex.main = 1)

plot(chronos(myprunedtree), edge.width = 1,
     cex = 0.7,
     main = "E) Ultrametric chronos, λ=1",
     cex.main = 1)
plot(chronos(myprunedtree, lambda = 0.1), edge.width = 1,
     cex = 0.7,
     main = "F) Ultrametric chronos, λ=0.1",
     cex.main = 1)
plot(chronos(myprunedtree, lambda = 0.01), edge.width = 1,
     cex = 0.7,
     main = "G) Ultrametric chronos, λ=0.01",
     cex.main = 1)
plot(chronos(myprunedtree, lambda = 0.001), edge.width = 1,
     cex = 0.7,
     main = "H) Ultrametric chronos, λ=0.001",
     cex.main = 1)

write.tree(chronos(myprunedtree, lambda = 0.1), file = "cagee2_tree.nwk")
                                                    # writing the tree with the lowest PHIIC score with no false convergence to use as input in the CAGEE command


### Clean and Normalise counts ###
counts_raw <- read.csv("C:/Users/Admin/Documents/rna/root/rootRNAcounts_ordered.tsv", header = TRUE, sep = "\t")
counts <- DGEList(counts = counts_raw[,c(1,7:42)], group =
                    as.factor(rep(c("sed","ult","sed","ult","sed", "ult"), rep.int(6,6))))  # native soil-type assignments
species <- as.factor(rep(c("calciphila","spPicNga","impolita","hequetiae","labillardierei","revolutissima"), rep.int(6,6) ))  # species assignments
counts$samples$species <- species
group <- counts$samples$group     # counts dataframe specifying species and soil preferences

annotation <- read.csv("C:/Users/Admin/Documents/rna/amin_redo/vie1167c_OmicsboxAminGene_Fannotation.txt", sep = "\t")
annotation <- annotation[annotation$SeqName %in% counts$genes$Geneid,]  # functional annotations of genes in counts table
annotation <- annotation[,c("SeqName", "Description","Length","GO.Names","GO.IDs")]

dim(counts)  #check initial number of genes
keep.exprs <- filterByExpr(counts, group=species)
counts <- counts[keep.exprs,, keep.lib.sizes=FALSE]  # 10/(median lib size) CPM cutoff met for at least 6 samples (size of a species group)
dim(counts)  #check final number of genes post-filtering

counts <- calcNormFactors(counts, method = "TMM")
counts[["samples"]][["norm.factors"]]  # this set doesn't need a lot of normalisation but just for sake of completion
cpm <- cpm(counts)  # CAGEE does a log conversion itself so we leave its input at the counts per million units


### Per species per gene CPM medians ###
cpm_df <- as.data.frame(cpm)
medians <- sapply(levels(species), function(sp) apply(cpm_df[, species == sp], 1, median))
medians <- as.data.frame(medians)
medians$Gene_ID <- counts[["genes"]]$Geneid
medians <- merge(annotation[,c(1,2)], medians,
                          by.x = "SeqName",   # column in annotation
                          by.y = "Gene_ID",   # column in medians table
                          all.y = TRUE)     # keep all genes that passed filtering
medians <- medians[,c(2,1,3:8)]
names(medians)[1] <- "DESC"
names(medians)[2] <- "Gene_ID"
head(medians)
write.table(medians, file = "cagee2_medians.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
                                                # writing the medians file to use as input (along with previously generated cagee2_tree.nwk) in the CAGEE command:
#cagee --cores 16 --tree cagee2_tree.nwk --infile cagee2_medians.tsv --sigma_tree cagee3_sigmatree.nwk -o C:/Users/Admin/Documents/rna/root/cagee/results2_3sigma
