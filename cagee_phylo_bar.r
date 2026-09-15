# The input for this R code is a collated tsv file with both leaf and root analysis results of the CAGEE command run with 3 evolutionary rates:
# cagee --cores 16 --tree cagee2_tree.nwk --infile cagee2_medians.tsv --sigma_tree cagee3_sigmatree.nwk -o C:/Users/Admin/Documents/rna/root/cagee/results2_3sigma

library(dplyr)
library(ggtree)
library(tidyr)
library(ggplot2)


### Set up phylogeny and load CAGEE results ###
tree <- read.tree(text = "(labillardierei:1,(revolutissima:0.8394620094,((impolita:0.4335132768,hequetiae:0.4335132768)<11>:0.2456226085,(calciphila:0.2076888523,spPicNga:0.2076888523)<10>:0.4714470331)<9>:0.1603261241)<8>:0.1605379906)<7>;")
dat <- read.table("3sigma_leaf_root.tsv", header = TRUE, sep = "\t")  # CAGEE results for leaves and roots collated into one file used as input here

min_branch <- min(dat$terminal_branch)
dat$terminal_branch <- dat$terminal_branch / min_branch
dat[c("leaf_in", "leaf_de", "root_in", "root_de")] <- dat[c("leaf_in", "leaf_de", "root_in", "root_de")] / dat$terminal_branch # scaling change by terminal branch length

soil_colors <- c("ultramafic" = "#FF8800", "non-ultramafic" = "#005F60")
dat$soil_color <- soil_colors[dat$soil_type]

bar_colors <- c("leaf_in" = "olivedrab",
                "leaf_de" = "#6B8E2350",
                "root_in" = "#9b7c62",
                "root_de" = "#d8ccc2")

p_tree <- ggtree(tree)  # plot base tree, get tip coordinates

tip_coords <- p_tree$data %>%
  filter(isTip) %>%
  left_join(dat, by = c("label" = "species"))


### Terminal nodes' bar plots setup ###
n_bars <- 4
bar_width <- 0.05
bar_spacing <- 0.06

max_value <- max(tip_coords %>% select(leaf_in, leaf_de, root_in, root_de) %>% unlist())  # max bar height used to adjust bar heights easily later on

bar_df <- tip_coords %>%
  select(label, x, y, leaf_in, leaf_de, root_in, root_de) %>%
  pivot_longer(cols = leaf_in:root_de, names_to = "variable", values_to = "value") %>%
  group_by(label) %>%
  mutate(
    bar_index = row_number() - 1,
    x = x + 0.2 + bar_index * bar_spacing,
    height = value / max_value * 0.5  # bar height
  ) %>%
  ungroup()


### Internal nodes' bar plots setup ###

node_coords <- p_tree$data %>% filter(!isTip) %>% select(node, x, y)  # Get internal node coordinates
node_dat <- read.table("internal_node_values.txt", header = TRUE, sep = "\t")  # Read internal node data

node_dat[c("leaf_in", "leaf_de", "root_in", "root_de")] <- 
  node_dat[c("leaf_in", "leaf_de", "root_in", "root_de")] / min_branch  # Normalize by branch length

node_bar_df <- node_coords %>%
  left_join(node_dat, by = "node") %>%
  pivot_longer(cols = leaf_in:root_de, names_to = "variable", values_to = "value") %>%
  group_by(node) %>%
  mutate(
    bar_index = row_number() - 1,
    x = x + 0.04 + bar_index * bar_spacing,  # smaller offset for nodes
    y = y + 0.0,
    height = value / max_value * 0.5  # bar height
  ) %>%
  ungroup() %>%
  rename(label = node) %>%
  mutate(label = as.character(label))  # convert to character for bind_rows

all_bar_df <- bind_rows(bar_df, node_bar_df)  # Combine tip and node bars

ref_value <- 500
ref_height <- ref_value / max_value * 0.5  # to add a scale for the bar charts


### Plot phylogeny with species and branches coloured by native soil type and evolutionary rate types respectively and bar charts showing number of genes with significant change in expression between two consecutive nodes ###
ggtree(tree, size = 1) %<+% tip_coords +
  geom_tree(aes(color = soil_color), size = 1) +
  geom_tiplab(aes(color = soil_color),
              align = TRUE,
              offset = 0.01,
              fontface = 3) +
  geom_rect(data = all_bar_df,
            aes(xmin = x - bar_width / 2,
                xmax = x + bar_width / 2,
                ymin = y,
                ymax = y + height,
                fill = variable),
            color = NA) +
  scale_color_identity(
    name = expression(atop("Native soil/",
                           paste("Evolutionary rate class (", sigma^2, ")"))),
    breaks = c("#005F60", "#FF8800"),
    labels = c("non-ultramafic", "ultramafic"),
    guide = "legend"
  ) +
  scale_fill_manual(values = bar_colors,
                    labels = c("decreased in leaf","increased in leaf","decreased in root","increased in root"),
                    name = expression(atop("gene expression in tissue", "(#genes scaled by branch length)"))) +
  theme(
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 8)
  ) +
  geom_treescale(x = 0, y = 6.5, width = 0.2, fontsize = 0, linesize = 1) +
  annotate("text", x = 0.1, y = 6.7, label = "unit branch length", size = 2.6) +
  geom_segment(aes(x = max(tip_coords$x) + 0.42,
                   xend = max(tip_coords$x) + 0.42,
                   y = max(tip_coords$y),
                   yend = max(tip_coords$y) + ref_height), linewidth = 1, colour = "darkgray") +
  annotate("text",
           x = max(tip_coords$x) + 0.43,
           y = max(tip_coords$y) + ref_height - 0.1,
           label = sprintf("%s genes/unit \n branch length", ref_value),
           size = 2.6, hjust = 0) +
  xlim(0, max(all_bar_df$x) + 0.1)
