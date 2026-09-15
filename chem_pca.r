library(ggplot2)
library(dplyr)
library(factoextra)
leafsoil <- read.csv("C:/Users/Admin/Documents/chem/leaf_and_soil_no_amb.tsv",
                     header = TRUE, sep = '\t')


######## LEAF ########
leaf <- dplyr::select(leafsoil, leaf_sample, species, soil_type:S_leaf, Co_leaf:Mo_leaf)  # Pb dropped due to extremely low resolution measurements
rownames(leaf) <- leaf$leaf_sample
leaf <- leaf[!rownames(leaf) %in% c("pan1146","pan1147"),]  # removing two doubtful leaf samples
leaf.pca <- prcomp(leaf[,4:21], scale. = TRUE)  # using all the previously selected elements for the PCA, normalisation: each variable to mean=0 & variance=1
pcadfl <- data.frame(leaf.pca$x)
pcadfl$sample <- rownames(pcadfl)
fviz_pca_biplot(leaf.pca, label="var", col.var = "gray", repel = TRUE) +
  ggtitle("Leaf") + geom_point(aes(shape=leaf$soil_type, col=leaf$species, size=3)) +
  scale_color_manual(values=c("darkorchid","plum","mediumpurple","thistle4","skyblue",
                              "orangered","royalblue","chocolate","brown3","navy","lightseagreen","aquamarine",
                              "khaki3", "seagreen2","palevioletred1", "gold", "darkgoldenrod3", "olivedrab2",
                              "darkolivegreen4", "darkred", "limegreen", "darksalmon", "steelblue4")) +
  geom_text(aes(label = leaf$leaf_sample), size=2) +
  guides(size="none") + guides(color = guide_legend(override.aes = list(size = 4)),
                               shape = guide_legend(override.aes = list(size = 4))) +
  labs(colour = "species", shape = "presumed soil preference")


######## SOIL ########
soil <- dplyr::select(leafsoil, assoc_soil_sample, species, soil_type, C_soil:Co_soil)  # dropped pH columns since there are not independent variables
soil <- soil %>% distinct(assoc_soil_sample, .keep_all = TRUE)  # since soil samples are fewer than leaf samples, occasionally the same soil sample
                                                                # is "associated" to >1 leaf sample, this step removes such duplicates for the PCA
rownames(soil) <- soil$assoc_soil_sample
soil.pca <- prcomp(soil[,4:20], scale. = TRUE)  # using all the previously selected elements for the PCA, normalisation: each variable to mean=0 & variance=1
pcadfs <- data.frame(soil.pca$x)
pcadfs$sample <- rownames(pcadfs)
fviz_pca_biplot(soil.pca, label="var", col.var = "gray", repel = TRUE) +
  ggtitle("Soil") + geom_point(aes(shape=soil$soil_type, col=soil$species, size=3)) +
  scale_color_manual(values=c("darkorchid","plum","mediumpurple","thistle4","skyblue",
                              "orangered","royalblue","chocolate","brown3","navy","lightseagreen","aquamarine",
                              "khaki3", "seagreen2","palevioletred1", "gold", "darkgoldenrod3", "olivedrab2",
                              "darkolivegreen4", "darkred", "limegreen", "darksalmon", "steelblue4")) +
  geom_text(aes(label = soil$assoc_soil_sample), size=2) +
  guides(size="none") + guides(color = guide_legend(override.aes = list(size = 4)),
                               shape = guide_legend(override.aes = list(size = 4))) +
  labs(colour = "species", shape = "presumed soil preference")
