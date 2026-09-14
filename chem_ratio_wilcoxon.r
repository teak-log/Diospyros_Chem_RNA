library(ggplot2)
library(dplyr)
library(ggsignif)
library(grid)
library(gridExtra)
library(gridtext)
library(introdataviz)


### Function to identify and replace outliers with NA ###
replace_outliers_with_na <- function(data, columns) {
  for (col in columns) {
    # Calculate Q1, Q3, and IQR
    Q1 <- quantile(data[[col]], 0.25, na.rm = TRUE)
    Q3 <- quantile(data[[col]], 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    # Define outlier thresholds
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    # Identify outliers
    outliers <- data[[col]] < lower_bound | data[[col]] > upper_bound
    # Display outliers
    if (any(outliers)) {
      cat("Outliers in column", col, ":\n")
      print(data[outliers, col, drop = FALSE])
      cat("\n") # Add a newline for better readability
      } else {
      cat("No outliers in column", col, "\n\n")
    }
    # Replace outliers with NA
    data[[col]][outliers] <- NA
  }
  return(data)
}

leaf_soil_type <- read.csv('C:/Users/Admin/Documents/chem/leaf_and_soil_no_amb_addNA.tsv', header = TRUE, sep = "\t")  # using a dataset filtered of possibly inaccurate measurements as indicated by the measurement facilities
rownames(leaf_soil_type) <- leaf_soil_type$leaf_sample  # each datapoint maps to a unique leaf sample
leaf_soil_type <- leaf_soil_type[!(row.names(leaf_soil_type) %in% c("pan1146","pan1147","eru1109")), ]  # removing a few doubtful samples


### Calculating the samplewise leaf:soil ratio of each element and appending them as new columns at the end of the original dataframe ###
leaf_soil_type_mod <- leaf_soil_type %>% mutate(C_ratio = C_leaf / C_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(N_ratio = N_leaf / N_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(Ca_ratio = Ca_leaf*1000 / Ca_soil)  # equalising measurement units in leaf & soil
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(Mg_ratio = Mg_leaf*1000 / Mg_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(Na_ratio = Na_leaf*1000 / Na_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(K_ratio = K_leaf*1000 / K_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(P_ratio = P_leaf / P_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(Cr_ratio = Cr_leaf / Cr_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(Fe_ratio = Fe_leaf / Fe_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(Mn_ratio = Mn_leaf / Mn_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(Ni_ratio = Ni_leaf / Ni_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(Al_ratio = Al_leaf / Al_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(Cu_ratio = Cu_leaf / Cu_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(Zn_ratio = Zn_leaf / Zn_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(S_ratio = S_leaf / S_soil)
leaf_soil_type_mod <- leaf_soil_type_mod %>% mutate(Co_ratio = Co_leaf / Co_soil)


### Organising ratio dataset & plotting half-violins showing Mann-Whitney U tests ###
anno <- list()  # initiate list of p-values
plot_list <- list()
sig_list <- c(46, 49, 52, 51, 48, 53, 55, 54, 61, 56, 58, 59)  # to help collect the above generated new columns in order of increasing atomic no. of the elements
leaf_soil_type_mod <- replace_outliers_with_na(leaf_soil_type_mod, sig_list)  # applying the outlier function on the ratios

#plotting grid of half violins with inset box plots for sedimentary vs ultramafic samples, showing p-values for Mann-Whitney U tests
for (j in 1:12) {
  i = as.integer(sig_list[j])
  anno[[j]] <- wilcox.test(
    leaf_soil_type_mod[leaf_soil_type_mod$soil_type == "sed", i],
    leaf_soil_type_mod[leaf_soil_type_mod$soil_type == "ult", i],)$p.value  # Mann-Whitney U test
  
  plot_list[[j]] <- ggplot(leaf_soil_type_mod,
                           aes_string(x=1,y=names(leaf_soil_type_mod)[i],
                                      fill="soil_type")) + geom_split_violin() +
    scale_fill_manual(values=c("#005F60", "#FF8800")) +
    geom_boxplot(width=0.2, linewidth=0.1, color="black", alpha=0.1, outlier.shape = NA) +
    guides(fill = "none") + theme_minimal() +
    geom_label(label=paste((sub("_ratio","",colnames(leaf_soil_type_mod)[i]))),
               x=0.6, y=max(as.numeric(na.omit(leaf_soil_type_mod[,i]))),
               color = "black", fill="white") +
    geom_signif(annotation = formatC(anno[[j]], digits = 1),
                y_position=max(as.numeric(na.omit(leaf_soil_type_mod[,i]))) +
                  0.1*max(as.numeric(na.omit(leaf_soil_type_mod[,i]))),
                xmin = 0.9, xmax = 1.1, textsize = 3.2) + labs(x = NULL, y=NULL) +
    coord_cartesian(clip = "off") + theme(axis.text.x = element_blank(),
                                          axis.ticks.x=element_blank())
}
soil_label = richtext_grob(text =
                           '<span style="color:black">Soil type (</span><span style="color:#005F60">sedimentary</span><span style="color:black">/</span><span style="color:#FF8800">ultramafic</span><span style="color:black">)</span>')
grid.arrange(grobs=plot_list, ncol=6, bottom=soil_label,left="leaf:soil concentrations of elements")
