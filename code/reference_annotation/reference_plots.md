Visualizing the Reference Annotations
================
Fallon Ratner
2023-11-23

## Setting up

-   Load libraries
-   Load dataframes
-   Visualize

``` r
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(stringr)
library(tidyr)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
#select own wd
setwd("C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/output/Annotation/Meta_Ref_Sctype/Reference - CellTypist")

#Read in CSV files:
datasets <- c("3 months", "6 months", "34_GA","24_GA",
              "22_GA", "17-18 GW") 
source_values <- c("Velasco_12w" = "3 months",
                   "Velasco_24w" = "6 months",
                  "Polioudakis_17-18GW"= "17-18 GW",
                   "Herring_22GW"= "22_GA",
                   "Herring_24GW"= "24_GA",
                   "Herring_34GW"= "34_GA")

# Initialize an empty dataframe for storing combined data
df_list <- list()
# Loop through each datasets
for (dataset in datasets) {
  # Construct the filename for each datasets
  csv_file <- paste0("proportion_",dataset, ".csv")
  
  # Read the csv with stringsAsFactors set to FALSE
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  # Add source as a column in the df 
  source_name <- names(source_values)[source_values == dataset]
  df$Source <- source_name
  
  # Append tmp_df to df_list
  df_list[[dataset]] <- df
}

# Combine all data.frames in df_list into one
ref <- do.call(rbind, df_list)
```

## Bar plot to see amount of cells labeled as unknown

``` r
#change column labels
ref$predicted_labels[ref$predicted_labels == "Unassigned"] <- "Unknown"
colnames(ref)[colnames(ref) == "predicted_labels"] <- "Cell Type"

#add new column: summary
ref <- ref %>%
  mutate(Summary = ifelse(`Cell Type` == 'Unknown', 'Unknown', 'Labeled'))

ref_summary <- ref %>%
  group_by(Source, Summary) %>%
  summarize(total_percentage = sum(Percentage))
```

    ## `summarise()` has grouped output by 'Source'. You can override using the
    ## `.groups` argument.

``` r
#make a bar plot
custom_colors <- c("Labeled"="blue",
                   "Unknown"="grey")
ggplot(ref_summary, aes(fill = Summary, y = total_percentage, x = Source)) + 
  geom_bar(stat = "identity", position = "stack") +
  labs(title="Cell Typist - Reference Annotation", y = "Percentage", x = "Dataset Sources", fill ="") + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(axis.text = element_text(size = 11),
        axis.text.y = element_text(size = 10)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) 
```

![](reference_plots_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Scatter plot of the top 5 labeled cell types

``` r
#select top 5 annotations
# Assuming df has columns: group, value
top_5_each_group <- ref %>%
  group_by(Source) %>%
  top_n(5, Percentage) %>%
  ungroup()


#Plot as bar and scatter
top_5_each_group$Source<- factor(top_5_each_group$Source, levels = c("Polioudakis_17-18GW","Herring_22GW", "Herring_24GW", "Herring_34GW", 
  "Velasco_12w","Velasco_24w", ordered = TRUE))


custom_colors <- c("Polioudakis_17-18GW"="#402C49","Herring_22GW"="#F4A460",                      "Herring_24GW"="#D2691E", "Herring_34GW"="#A0522D", 
                   "Velasco_12w"="#DEB887","Velasco_24w" = "#D2B48C")

top_5_each_group$`Cell Type` <- factor(top_5_each_group$`Cell Type`, 
                                      levels = c("Ventral midbrain radial glia",
                                                  "Striatum radial glia",
                                                  "Telencephalon radial glia", 
                                                  "Cortex neuronal IPC",
                                                  "Forebrain neuronal IPC",
                                                  "Cortex neuron",
                                                  "Forebrain neuron",
                                                  "Ventral midbrain neuron",
                                                  "Dorsal midbrain neuron",
                                                  "Hippocampus neuron",
                                                  "Cortex neuron|Subcortex neuron",
                                                  "Cortex neuron|Hippocampus neuron",
                                                  "Subcortex neuron",
                                                  "Pons glioblast",
                                                  "Unknown",
                                                  ordered = TRUE))

#set all cell types that are not present in each dataset to 0
df_all_combinations <- expand.grid(Source = unique(top_5_each_group$Source), `Cell Type` = unique(top_5_each_group$`Cell Type`))
# Merge the new data frame with the original data frame
combo_full <- merge(df_all_combinations, top_5_each_group, by = c("Source", "Cell Type"), all.x = TRUE)

combo_full <- replace(combo_full, is.na(combo_full), 0)

# scatter plot
ggplot(top_5_each_group, aes(x=`Cell Type`, y=Percentage, color = Source)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  scale_color_manual(values = custom_colors) +
  labs(title="Cell Typist - Reference Annotation",
       x="Cell Type", y = "Proportion (%)")+ 
  theme_classic() +
  scale_x_discrete(labels = function(x) str_replace_all(x, " ", "\n"))+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text = element_text(size = 8),
    axis.text.y = element_text(size = 10),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.text = element_text(size = 8),  # Smaller legend text
    legend.title = element_text(size = 7) )
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](reference_plots_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
