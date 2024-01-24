Visualizing Annotations of Fetal Brain & Brain Organoid scRNAseq
================
Fallon Ratner
2023-11-23

## Setting up

-   Load libraries
-   Load dataframes
-   Calculate the proportion as a percentage
-   Plot stacked area chart

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

setwd("C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/output/Annotation/Final_10_08")


#Make a dataframe for the in vitro datasets
# List the in vitro datasets
datasets <- c("han11.1", "han11.2", "han12","han13",
             "polio", "liu","herr22", "herr24", "herr34",
             "herr2d", "cou13", "cou17", "cou19") 

source_values <- c("Han_GW11.1" = "han11.1",
                   "Han_GW11.2" = "han11.2",
                   "Han_GW12" = "han12",
                   "Han_GW13"= "han13",
                   "Couturier_GW13"= "cou13",
                   "Couturier_GW17"= "cou17",
                   "Couturier_GW19"= "cou19",
                   "Polioudakis_GW17-18"= "polio",
                   "Liu_GW17-19"= "liu",
                   "Herring_GW22"= "herr22",
                   "Herring_GW24"= "herr24",
                   "Herring_GW34"= "herr34",
                   "Herring_Day2"= "herr2d"
)
# Initialize an empty dataframe for storing combined data
vivo <- data.frame()
# Loop through each datasets
for (dataset in datasets) {
  # Construct the filename for each datasets
  csv_file <- paste0(dataset, "_sctype25_7.csv")
  
  # Read the csv with stringsAsFactors set to FALSE
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  # Add source as a column in the df 
  source_name <- names(source_values)[source_values == dataset]
  df$source <- source_name
  
  # Combine into one dataframe
  if (nrow(vivo) == 0) {
    vivo <- df
  } else {
    vivo <- rbind(vivo, df)
  }
}


# Initialize an empty dataframe for storing combined data
vivo_gb <- data.frame()
# Loop through each datasets
for (dataset in datasets) {
  # Construct the filename for each datasets
  csv_file <- paste0(dataset, "_sctype_INs_25_7.csv")
  
  # If the file does not exist, print a warning and skip to the next dataset
  if (!file.exists(csv_file)) {
    print(paste("Warning: File", csv_file, "does not exist. Skipping..."))
    next
  }
  
  # Read the csv with stringsAsFactors set to FALSE
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  # Add source as a column in the df
  source_name <- names(source_values)[source_values == dataset]
  
  # If no match is found, it's a problem
  if(length(source_name) == 0) {
    stop(paste("No source value found for dataset:", dataset))
  }
  
  df$source <- source_name
  
  # Combine into one dataframe
  if (nrow(vivo_gb) == 0) {
    vivo_gb <- df
  } else {
    vivo_gb <- rbind(vivo_gb, df)
  }
}
```

    ## [1] "Warning: File han11.1_sctype_INs_25_7.csv does not exist. Skipping..."
    ## [1] "Warning: File liu_sctype_INs_25_7.csv does not exist. Skipping..."

``` r
#Read in the in vitro data
vitro <- read.csv("outputvitro_int_sctype1011.csv")
vitro_gb <- read.csv("outputvitro_gaba_int_sctype1011.csv")
```

## Make a Stacked bar chart for the Fetal Brain Datasets

``` r
#Add a column with the value ages to vivo data
#Polioudakis & Liu ages are the average value and Herring Day2 is converted to weeks
age_values <- c("Han_GW11.1" = 11,
                   "Han_GW11.2" = 11,
                   "Han_GW12" = 12,
                   "Han_GW13"= 13,
                   "Couturier_GW13"= 13,
                   "Couturier_GW17"= 17,
                   "Couturier_GW19"= 17,
                   "Polioudakis_GW17-18"= 17.5,
                   "Liu_GW17-19"= 18,
                   "Herring_GW22"= 22,
                   "Herring_GW24"= 24,
                   "Herring_GW34"= 34,
                   "Herring_Day2"= 40)

vivo$age <- age_values[vivo$source]

# Calculate total cells for each age
age_totals <- vivo %>%
  group_by(age) %>%
  summarise(total_cells_for_age = sum(ncells))

# Join the totals back to the original data and calculate the percentage
avg_percentages <- vivo %>%
  group_by(age, type) %>%
  summarise(total_cells = sum(ncells), .groups = "keep") %>%
  left_join(age_totals, by = "age") %>%
  mutate(percentage = (total_cells / total_cells_for_age) * 100) %>%
  select(-total_cells_for_age) 

#Set the levels for the data for the specified order in the stacked area plot, and give each cell type a specific color
# Give a specific order:
avg_percentages$type <- factor(avg_percentages$type, levels=c("outer radial glia", "ventral Radial Glia", "Intermediate Progenitors","Neuroblasts","Immature Excitatory", "Maturing Excitatory", "Migrating Excitatory","Deep Excitatory Layers", "Upper Excitatory Layers","Interneuron Precursors", "MGE INs", "CGE INs","Glioblasts", "Immature Astrocytes", "Mature Astrocytes", "OPCs","Oligodendrocytes", "Microglia", "Endothelial Cells", "Mural") )

colors <- c("CGE INs" = "#2E8B57", "MGE INs" = "#15693A", "Interneuron Precursors" = "#3CB371","Neuroblasts" = "#CAE7FA","Glioblasts" = "#DDA0DD", "Immature Astrocytes" = "#DB7093","Mature Astrocytes" = "#C71585","Immature Excitatory" = "#87CEFA",
"Maturing Excitatory" = "#00BFFF","Migrating Excitatory" = "#6495ED","Deep Excitatory Layers" = "#7B68EE", "Upper Excitatory Layers" = "#6A5ACD", "Microglia" = "#8A2BE2", "outer radial glia" = "#DEB887", "ventral Radial Glia" = "#D2B48C","Intermediate Progenitors" = "#F4A460","Endothelial Cells" = "#B22222", "Mural" = "#BDB76B","OPCs" ="#BA55D3", "Oligodendrocytes"="#9932CC" ) 

avg_percentages <- avg_percentages %>% arrange(age, type)

ggplot(avg_percentages, aes(x=age, y=percentage, fill=type)) + 
  geom_area(position = 'stack') +
  scale_fill_manual(values = colors) + 
  labs(title = "Neurodevelopment In Vivo",y="Average Proportion (%)", x= "Age (GW)", fill= "Cell Type") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white"), 
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"))+
   theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
```

![](visualization_fetal_organoids_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggplot(avg_percentages, aes(x=as.factor(age), y=percentage, fill=type)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = colors) + 
  labs(title = "Neurodevelopment in the Fetal Brain",y="Average Proportion (%)", x= "Age (GW)", fill= "Cell Type") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white"), 
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"))+
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
```

![](visualization_fetal_organoids_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

## Make a Stacked bar chart for the Fetal Brain Interneuron Sub-Types

``` r
vivo_gb$age <- age_values[vivo_gb$source]

# Calculate total cells for each age
age_totals_gb <- vivo_gb %>%
  group_by(age) %>%
  summarise(total_cells_for_age = sum(ncells))

# Join the totals back to the original data and calculate the percentage
avg_percentages_gb <- vivo_gb %>%
  group_by(age, type) %>%
  summarise(total_cells = sum(ncells), .groups = "keep") %>%
  left_join(age_totals_gb, by = "age") %>%
  mutate(percentage = (total_cells / total_cells_for_age) * 100) %>%
  select(-total_cells_for_age) 

avg_percentages_gb$type <- factor(avg_percentages_gb$type, levels=c("ID2","VIP", "NDNF", "SST", "PVALB", "Unknown") )


colors2 <- c("ID2" = "#7DECB4","VIP" ="#4DC086","NDNF" = "#359264","SST"= "#8E92F6", "PVALB" ="#3B40C9", "Unknown" ="#696969")

ggplot(avg_percentages_gb, aes(x=age, y=percentage, fill=type)) + 
  geom_area() +
  scale_fill_manual(values = colors2) + 
  labs(title = "Interneuron Development in the Fetal Brain",y="Average Proportion (%)", x= "Age (GW)", fill= "Cell Type") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white"),  # Ensure entire plot has a white background
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )
```

![](visualization_fetal_organoids_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggplot(avg_percentages_gb, aes(x=as.factor(age), y=percentage, fill=type)) +   geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = colors2) + 
  labs(title = "Interneuron Development in the Fetal Brain",y="Average Proportion (%)", x= "Age (GW)", fill= "Cell Type") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white"), 
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"))+
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
```

![](visualization_fetal_organoids_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

## Make a Stacked bar chart for Brain Organoids

``` r
#Add age column to vitro df
age_values2 <- c("Bhaduri_H1S_10w" = 10,
                "Bhaduri_H1X_10w" = 10,
                "Bhaduri_H1S_5w" = 5,
                "Bhaduri_H1X_5w"= 5,
                "Bhaduri_H1S_8w"= 8,
                "Bhaduri_H1X_8w"= 8,
                "Fair_12w"= 12,
                "Fair_20w"= 20,
                "Giandomenico_10w"= 10,
                "Madhavan_12w"= 12,
                "Popova_batch1_7w"= 7,
                "Popova_batch2_7w"= 7,
                "Trujillo_12w"= 12,
                "Trujillo_24w"= 24,
                "Trujillo_40w"= 40,
                "Trujillo_4w"= 4,
                "Velasco_12w"= 12,
                "Velasco_24w"= 24,
                "Xiang_10w"= 10,
                "Xiang_11w"= 11,
                "Xiangb1_4w"= 4,
                "Xiangb2_4w"= 4)

vitro$age <- age_values2[vitro$source]
# Calculate total cells for each age
age_totals2 <- vitro %>%
  group_by(age) %>%
  summarise(total_cells_for_age = sum(ncells))

# Join the totals back to the original data and calculate the percentage
avg_percentages2 <- vitro %>%
  group_by(age, cell_type) %>%
  summarise(total_cells = sum(ncells), .groups = "keep") %>%
  left_join(age_totals2, by = "age") %>%
  mutate(percentage = (total_cells / total_cells_for_age) * 100) %>%
  select(-total_cells_for_age) 

# Custom color palette
colors_org <- c("CGE INs" = "#2E8B57", "MGE INs" = "#15693A", "Interneuron Precursors" = "#3CB371",
            "Glioblasts" = "#DDA0DD", "Immature Astrocytes" = "#DB7093","Mature Astrocytes" = "#C71585",
            "Neuroblasts" = "#CAE7FA","Immature Excitatory" = "#87CEFA",
            "Maturing Excitatory" = "#00BFFF","Migrating Excitatory" = "#6495ED",
            "Deep Excitatory Layers" = "#7B68EE", "Upper Excitatory Layers" = "#6A5ACD", "Microglia" = "#8A2BE2", 
            "outer radial glia" = "#DEB887", "ventral Radial Glia" = "#D2B48C",
            "Intermediate Progenitors" = "#F4A460","Endothelial Cells" = "#B22222", "Mural" = "#BDB76B",
            "OPCs" ="#BA55D3", "Oligodendrocytes"="#9932CC", "Mixed" = "#696969") 

# Give a specific order:
avg_percentages2 $cell_type <- factor(avg_percentages2 $cell_type, levels=c("outer radial glia", "ventral Radial Glia", "Intermediate Progenitors",
                                                                        "Neuroblasts","Immature Excitatory", "Maturing Excitatory", "Migrating Excitatory",
                                                                        "Deep Excitatory Layers", "Upper Excitatory Layers", "Mixed",
                                                                        "Interneuron Precursors", "MGE INs", "CGE INs",
                                                                        "Glioblasts", "Immature Astrocytes", "Mature Astrocytes", "OPCs",
                                                                        "Oligodendrocytes", "Microglia", "Endothelial Cells", "Mural") )

ggplot(avg_percentages2 , aes(x=age, y=percentage, fill=cell_type)) + 
  geom_area() +
  scale_fill_manual(values = colors_org) + 
  labs(title = "Neurodevelopment in the Fetal Brain",y="Average Proportion (%)", x= "Age (Weeks in Culture)", fill= "Cell Type") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white"),  # Ensure entire plot has a white background
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )
```

![](visualization_fetal_organoids_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggplot(avg_percentages2, aes(x=as.factor(age), y=percentage, fill=cell_type)) + 
  geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values = colors_org) + 
  labs(title = "Neurodevelopment in Brain Organoids",y="Average Proportion (%)",  x= "Age (Weeks in Culture)", fill= "Cell Type") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white"), 
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"))+
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
```

![](visualization_fetal_organoids_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

## Make a Stacked bar chart for Brain Organoid Interneuron Sub-Types

``` r
vitro_gb$age <- age_values2[vitro_gb$source]

# Calculate total cells for each age
age_totals_gb2 <- vitro_gb %>%
  group_by(age) %>%
  summarise(total_cells_for_age = sum(ncells))

# Join the totals back to the original data and calculate the percentage
avg_percentages_gb2 <- vitro_gb %>%
  group_by(age, cell_type) %>%
  summarise(total_cells = sum(ncells), .groups = "keep") %>%
  left_join(age_totals_gb2, by = "age") %>%
  mutate(percentage = (total_cells / total_cells_for_age) * 100) %>%
  select(-total_cells_for_age) 

avg_percentages_gb2$cell_type <- factor(avg_percentages_gb2$cell_type, levels=c("ID2","VIP", "NDNF", "SST", "PVALB", "Unknown") )


colors2 <- c("ID2" = "#7DECB4","VIP" ="#4DC086","NDNF" = "#359264","SST"= "#8E92F6", "PVALB" ="#3B40C9", "Unknown" ="#696969")

ggplot(avg_percentages_gb2, aes(x=age, y=percentage, fill=cell_type)) + 
  geom_area() +
  scale_fill_manual(values = colors2) + 
  labs(title = "Interneuron Development in Brain Organoids",y="Average Proportion (%)", x= "Age (GW)", fill= "Cell Type") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white"),  # Ensure entire plot has a white background
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )
```

![](visualization_fetal_organoids_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggplot(avg_percentages_gb2, aes(x=as.factor(age), y=percentage, fill=cell_type)) +   geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = colors2) + 
  labs(title = "Interneuron Development in Brain Organoids",y="Average Proportion (%)", x= "Age (GW)", fill= "Cell Type") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white"), 
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"))+
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
```

![](visualization_fetal_organoids_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

## Make a Violin Plot

``` r
# Remove unnecessary columns
vivo <- vivo %>% select(-cluster, -scores)
vitro <- vitro %>% select(-seurat_cluster)

colnames(vivo)[colnames(vivo) == "type"] <- "cell_type"

#calculate the proportion as a percentage
vivo <- vivo %>%
  group_by(source, cell_type) %>%
  summarize(total_cells = sum(ncells), .groups = "drop") %>%
  group_by(source) %>%
  mutate(percentage = (total_cells / sum(total_cells)) * 100)

vitro <- vitro %>%
  group_by(source, cell_type) %>%
  summarize(total_cells = sum(ncells), .groups = "drop") %>%
  group_by(source) %>%
  mutate(percentage = (total_cells / sum(total_cells)) * 100)

# Add another column called condition and set custom colors & order
vivo$condition <- "Vivo"
vitro$condition <- "Vitro"
# Combine the data frames
vivo_vitro <- rbind(vivo, vitro)

custom_colors4 <- c("Vivo"="#808080","Vitro"="#0000FF")

vivo_vitro$condition <- factor(vivo_vitro$condition, levels = c("Vivo","Vitro"),ordered = TRUE)

#Order the Cells 
vivo_vitro$cell_type <- factor(vivo_vitro$cell_type, levels = c("outer radial glia", "ventral Radial Glia", "Intermediate Progenitors",
                                                                  "Neuroblasts","Immature Excitatory", "Maturing Excitatory", "Migrating Excitatory",
                                                                  "Deep Excitatory Layers", "Upper Excitatory Layers", "Mixed",
                                                                  "Interneuron Precursors", "MGE INs", "CGE INs",
                                                                  "Glioblasts", "Immature Astrocytes", "Mature Astrocytes", "OPCs",
                                                                  "Oligodendrocytes", "Microglia", "Endothelial Cells", "Mural"), ordered = TRUE)

ggplot(vivo_vitro, aes(x = condition, y = percentage, color = condition)) + 
  geom_violin(trim = FALSE) + 
  labs(title="Cell Types in Fetal Brain vs Brain Organoids", y = "Proportion (%)")+
  geom_jitter(width = 0.1, size = 1, aes(color = condition)) +  # Add individual data points
  facet_wrap(vars(cell_type)) + 
  theme_bw(base_size = 3) + 
  scale_color_manual(values = custom_colors4) +
  theme(legend.position = "top", 
        legend.justification = c(1, 1), 
        legend.box.just = "right",
        panel.background = element_blank(),  # Set the background to white
        panel.grid = element_blank(),  # Remove the grid
        axis.line = element_line(color = "black")) +
  stat_summary(fun.min = function(y) mean(y) - sd(y),
               fun = mean,
               fun.max = function(y) mean(y) + sd(y),
               geom = "pointrange", 
               color = "black")
```

    ## Warning: Groups with fewer than two data points have been dropped.

    ## Warning: Removed 1 rows containing missing values (`geom_segment()`).

![](visualization_fetal_organoids_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
