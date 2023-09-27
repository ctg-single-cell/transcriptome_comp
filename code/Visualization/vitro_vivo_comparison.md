Visualizing In Vivo vs In Vitro
================
18-09-2023

## Import Libraries

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
```

## Make a look-up map for the in vitro data

``` r
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
```

## Load each csv file, calculate the proportion, and combine into one dataframe for the in vivo data

``` r
vivo <- data.frame()

for (dataset in datasets) {
  csv_file <- paste0(dataset, "_sctype25_7.csv")
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  df <- df %>%
    group_by(type) %>%
    summarize(total_cells = sum(ncells), .groups = "drop") %>%
    mutate(percentage = (total_cells / sum(total_cells)) * 100)
  
  source_name <- names(source_values)[source_values == dataset]
  df$source <- source_name
  
  if (nrow(vivo) == 0) {
    vivo <- df
  } else {
    vivo <- rbind(vivo, df)
  }
}
```

## Load each csv file, calculate the proportion, and combine into one dataframe for the in vivo interneuron data

``` r
vivo_gb <- data.frame()

for (dataset in datasets) {
  csv_file <- paste0(dataset, "_sctype_INs_25_7.csv")
  
  if (!file.exists(csv_file)) {
    print(paste("Warning: File", csv_file, "does not exist. Skipping..."))
    next
  }
  
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  df <- df %>%
    group_by(type) %>%
    summarize(total_cells = sum(ncells), .groups = "drop") %>%
    mutate(percentage = (total_cells / sum(total_cells)) * 100)
  
  source_name <- names(source_values)[source_values == dataset]
  
  if(length(source_name) == 0) {
    stop(paste("No source value found for dataset:", dataset))
  }
  
  df$source <- source_name
  
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
head(vivo)
```

    ## # A tibble: 6 × 4
    ##   type                   total_cells percentage source    
    ##   <chr>                        <int>      <dbl> <chr>     
    ## 1 Glioblasts                     191       9.49 Han_GW11.1
    ## 2 Immature Astrocytes            349      17.3  Han_GW11.1
    ## 3 Interneuron Precursors         197       9.79 Han_GW11.1
    ## 4 Maturing Excitatory            773      38.4  Han_GW11.1
    ## 5 Microglia                       85       4.22 Han_GW11.1
    ## 6 Mural                           46       2.29 Han_GW11.1

## Load the in vitro csv files and calculate the proportions

``` r
vitro <- read.csv("outputvitro_int_sctype.csv")
vitro_gb <- read.csv("outputvitro_gaba_int_sctype.csv")

vitro <- vitro %>%
  group_by(source, cell_type) %>%
  summarize(total_cells = sum(ncells), .groups = "drop") %>%
  group_by(source) %>%
  mutate(percentage = (total_cells / sum(total_cells)) * 100)

vitro_gb <- vitro_gb %>%
  group_by(source, cell_type) %>%
  summarize(total_cells = sum(ncells), .groups = "drop") %>%
  group_by(source) %>%
  mutate(percentage = (total_cells / sum(total_cells)) * 100)
```

## Combine in vivo and in vitro data into one df

``` r
vivo <- vivo %>% rename(cell_type = type)
vivo_gb <- vivo_gb %>% rename(cell_type = type)

vivo_vitro <- rbind(vivo,vitro)
vivo_vitro_gb <- rbind(vivo_gb,vitro_gb)
```

## Wrangle the dataframe so cell types that are not present are given a value of 0

``` r
df_all_combinations <- expand.grid(source = unique(vivo_vitro$source), cell_type = unique(vivo_vitro$cell_type))

combo_full <- merge(df_all_combinations, vivo_vitro, by = c("source", "cell_type"), all.x = TRUE)
combo_full <- replace(combo_full, is.na(combo_full), 0)
```

## Set the levels for the data for the specified order in the bar plot, and give each dataset a specific color

``` r
combo_full$source <- factor(combo_full$source, levels = c("Han_GW11.1","Han_GW11.2", "Han_GW12", "Han_GW13", "Couturier_GW13", 
                                                         "Couturier_GW17", "Couturier_GW19","Polioudakis_GW17-18", "Liu_GW17-19",
                                                         "Herring_GW22", "Herring_GW24", "Herring_GW34", "Herring_Day2",
                                                         "Bhaduri_H1S_5w","Bhaduri_H1X_5w", "Bhaduri_H1S_8w",
                                                         "Bhaduri_H1X_8w","Bhaduri_H1S_10w", "Bhaduri_H1X_10w",
                                                         "Trujillo_4w", "Trujillo_12w", "Trujillo_24w", "Trujillo_40w",
                                                         "Xiangb1_4w", "Xiangb2_4w","Xiang_10w","Xiang_11w",
                                                         "Fair_12w", "Fair_20w","Velasco_12w","Velasco_24w",
                                                         "Popova_batch1_7w", "Popova_batch2_7w",
                                                         "Giandomenico_10w", "Madhavan_12w"), ordered = TRUE)

#Define custom color palette
custom_colors <- c( "Han_GW11.1"="#D3D3D3","Han_GW11.2"=    "#A9A9A9", "Han_GW12"="#808080", "Han_GW13"= "#696969", "Couturier_GW13"="#FFD700", 
                    "Couturier_GW17"="#DAA520", "Couturier_GW19"="#B8860B","Polioudakis_GW17-18"="#402C49", "Liu_GW17-19"="#BC8F8F",
                    "Herring_GW22"="#F4A460", "Herring_GW24"="#D2691E", "Herring_GW34"="#A0522D", "Herring_Day2"="#8B4513",
  "Bhaduri_H1S_5w"="#6495ED", "Bhaduri_H1X_5w"= "#1E90FF", "Bhaduri_H1S_8w"="#0000FF", 
  "Bhaduri_H1X_8w"="#0000CD", "Bhaduri_H1S_10w"="#00008B","Bhaduri_H1X_10w"="#000080", "Fair_12w"="#40E0D0",
  "Fair_20w"="#48D1CC", "Giandomenico_10w"="#008080", "Madhavan_12w"="#C44733", 
  "Popova_batch1_7w"="#9E7DBA", "Popova_batch2_7w"= "#896FA1", "Trujillo_4w"="#A9414B", 
  "Trujillo_12w"="#9B2A35", "Trujillo_24w"="#841E28","Trujillo_40w"="#730F19", "Velasco_12w"="#DEB887",
  "Velasco_24w" = "#D2B48C","Xiangb1_4w"="#599035", "Xiangb2_4w"="#598A38","Xiang_10w" = "#426828","Xiang_11w"="#2F5416")
```

\##Make a bar chart for all of the cell types

``` r
bar <- ggplot(combo_full, aes(fill = source, y = percentage, x = cell_type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

bar + labs(title="scType Annotation", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 7)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
```

![](vitro_vivo_comparison_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Make a Stacked bar chart

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

#Calculate the average proportion for the same ages
avg_percentages <- vivo %>%
  group_by(age,cell_type) %>%
  summarise(avg_percentage = mean(percentage, na.rm = TRUE))
```

    ## `summarise()` has grouped output by 'age'. You can override using the `.groups`
    ## argument.

``` r
#Set the levels for the data for the specified order in the stacked area plot, and give each cell type a specific color
# Give a specific order:
avg_percentages$cell_type <- factor(avg_percentages$cell_type, levels=c("outer radial glia", "ventral Radial Glia", "Intermediate Progenitors",
                                                                        "Immature Excitatory", "Maturing Excitatory", "Migrating Excitatory",
                                                                        "Deep Excitatory Layers", "Upper Excitatory Layers",
                                                                        "Interneuron Precursors", "MGE INs", "CGE INs",
                                                                        "Glioblasts", "Immature Astrocytes", "Mature Astrocytes", "OPCs",
                                                                        "Oligodendrocytes", "Microglia", "Endothelial Cells", "Mural") )

colors <- c("CGE INs" = "#2E8B57", "MGE INs" = "#15693A", "Interneuron Precursors" = "#3CB371",
            "Glioblasts" = "#DDA0DD", "Immature Astrocytes" = "#DB7093","Mature Astrocytes" = "#C71585","Immature Excitatory" = "#87CEFA",
            "Maturing Excitatory" = "#00BFFF","Migrating Excitatory" = "#6495ED",
             "Deep Excitatory Layers" = "#7B68EE", "Upper Excitatory Layers" = "#6A5ACD", "Microglia" = "#8A2BE2", 
             "outer radial glia" = "#DEB887", "ventral Radial Glia" = "#D2B48C",
            "Intermediate Progenitors" = "#F4A460","Endothelial Cells" = "#B22222", "Mural" = "#BDB76B",
            "OPCs" ="#BA55D3", "Oligodendrocytes"="#9932CC" ) 


ggplot(avg_percentages, aes(x=age, y=avg_percentage, fill=cell_type)) + 
  geom_area() +
  scale_fill_manual(values = colors) + 
  labs(title = "Neurodevelopment In Vivo",y="Average Proportion (%)", x= "Age (GW)", fill= "Cell Type") + 
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white"), 
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text.y = element_blank())+
   theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
```

![](vitro_vivo_comparison_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Make a Rose Plot: Combines Bar and Pie Plots

``` r
ggplot(data=combo_full,aes(x=cell_type,y=percentage,fill=source))+
  geom_bar(stat="identity") +
  coord_polar() + scale_fill_manual(values = custom_colors) +
  xlab("")+ylab("") + labs(y= "Proportion(%)",fill= "scRNAseq") +
  theme(axis.text.x = element_text(size=5, vjust=0.3))
```

![](vitro_vivo_comparison_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## Make a Scatter Plot

``` r
vivo_vitro$source<- factor(vivo_vitro$source, levels = c("Han_GW11.1","Han_GW11.2", "Han_GW12", "Han_GW13", "Couturier_GW13", 
                                                         "Couturier_GW17", "Couturier_GW19","Polioudakis_GW17-18", "Liu_GW17-19",
                                                         "Herring_GW22", "Herring_GW24", "Herring_GW34", "Herring_Day2",
                                                         "Bhaduri_H1S_5w","Bhaduri_H1X_5w", "Bhaduri_H1S_8w",
                                                         "Bhaduri_H1X_8w","Bhaduri_H1S_10w", "Bhaduri_H1X_10w",
                                                         "Trujillo_4w", "Trujillo_12w", "Trujillo_24w", "Trujillo_40w",
                                                         "Xiangb1_4w", "Xiangb2_4w","Xiang_10w","Xiang_11w",
                                                         "Fair_12w", "Fair_20w","Velasco_12w","Velasco_24w",
                                                         "Popova_batch1_7w", "Popova_batch2_7w",
                                                         "Giandomenico_10w", "Madhavan_12w"), ordered = TRUE)



custom_colors <- c( "Han_GW11.1"="#808080","Han_GW11.2"=    "#808080", "Han_GW12"="#808080", "Han_GW13"= "#808080", "Couturier_GW13"="#808080", 
                    "Couturier_GW17"="#808080", "Couturier_GW19"="#808080","Polioudakis_GW17-18"="#808080", "Liu_GW17-19"="#808080",
                    "Herring_GW22"="#808080", "Herring_GW24"="#808080", "Herring_GW34"="#808080", "Herring_Day2"="#808080",
                    "Bhaduri_H1S_5w"="#0000FF", "Bhaduri_H1X_5w"= "#0000FF", "Bhaduri_H1S_8w"="#0000FF", 
                    "Bhaduri_H1X_8w"="#0000FF", "Bhaduri_H1S_10w"="#0000FF","Bhaduri_H1X_10w"="#0000FF", "Fair_12w"="#0000FF",
                    "Fair_20w"="#0000FF", "Giandomenico_10w"="#0000FF", "Madhavan_12w"="#0000FF", 
                    "Popova_batch1_7w"="#0000FF", "Popova_batch2_7w"= "#0000FF", "Trujillo_4w"="#0000FF", 
                    "Trujillo_12w"="#0000FF", "Trujillo_24w"="#0000FF","Trujillo_40w"="#0000FF", "Velasco_12w"="#0000FF",
                    "Velasco_24w" = "#0000FF","Xiangb1_4w"="#0000FF", "Xiangb2_4w"="#0000FF","Xiang_10w" = "#0000FF","Xiang_11w"="#0000FF")

ggplot(vivo_vitro, aes(x=cell_type, y=percentage, color = source)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  scale_color_manual(values = custom_colors) +
  labs(title="Cell Types in Fetal Brain vs Brain Organoids",
       x="Cell Type", y = "Proportion (%)")+ 
  theme_classic() +
 scale_x_discrete(labels = function(x) str_replace_all(x, " ", "\n"))+
  theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  axis.text = element_text(size = 8),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.text = element_text(size = 6),  # Smaller legend text
  legend.title = element_text(size = 7) )
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](vitro_vivo_comparison_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## Make a Violin Plot

``` r
#Remove age column
vivo <-vivo %>% select(-age)
#Add another column called condition and set custom colors & order
vivo$condition <- "Vivo"
vitro$condition <- "Vitro"

vivo_vitro2 <- rbind(vivo,vitro)

custom_colors4 <- c("Vivo"="#808080","Vitro"="#0000FF")

vivo_vitro2$condition<- factor(vivo_vitro2$condition, levels = c("Vivo","Vitro"),ordered = TRUE)

#Order the Cells 
vivo_vitro2$cell_type <- factor(vivo_vitro2$cell_type, levels = c("outer radial glia", "ventral Radial Glia", "Intermediate Progenitors",
                                                                  "Neuroblasts","Immature Excitatory", "Maturing Excitatory", "Migrating Excitatory",
                                                                  "Deep Excitatory Layers", "Upper Excitatory Layers", "Mixed",
                                                                  "Interneuron Precursors", "MGE INs", "CGE INs",
                                                                  "Glioblasts", "Immature Astrocytes", "Mature Astrocytes", "OPCs",
                                                                  "Oligodendrocytes", "Microglia", "Endothelial Cells", "Mural"), ordered = TRUE)

 ggplot(vivo_vitro2, aes(x = condition, y = percentage, color = condition)) + 
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

![](vitro_vivo_comparison_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## Do Singificance Testing: Wilcoxon Rank Sum Test

``` r
results <- vivo_vitro2 %>%
  group_by(cell_type) %>%
  filter(n_distinct(condition) == 2) %>%
  summarize(wilcox_p_value = wilcox.test(percentage ~ condition)$p.value)

p_values <- results$wilcox_p_value

bonferroni_corrected_p <- p.adjust(p_values, method = "bonferroni")
results$bonferroni_p <- bonferroni_corrected_p
print(results)
```

    ## # A tibble: 14 × 3
    ##    cell_type                wilcox_p_value bonferroni_p
    ##    <ord>                             <dbl>        <dbl>
    ##  1 outer radial glia              0.166         1      
    ##  2 ventral Radial Glia            0.744         1      
    ##  3 Intermediate Progenitors       0.00404       0.0566 
    ##  4 Neuroblasts                    0.154         1      
    ##  5 Immature Excitatory            0.203         1      
    ##  6 Maturing Excitatory            0.969         1      
    ##  7 Upper Excitatory Layers        1             1      
    ##  8 Interneuron Precursors         0.262         1      
    ##  9 MGE INs                        0.000198      0.00277
    ## 10 CGE INs                        0.0832        1      
    ## 11 Immature Astrocytes            0.0233        0.326  
    ## 12 OPCs                           0.456         1      
    ## 13 Endothelial Cells              0.0121        0.170  
    ## 14 Mural                          0.521         1
