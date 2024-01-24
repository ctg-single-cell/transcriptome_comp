Visualizing the Author’s Annotations
================
Fallon Ratner
2023-11-22

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
library(VennDiagram)
```

    ## Loading required package: grid

    ## Loading required package: futile.logger

``` r
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
#select own wd
setwd("C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/output/Annotation/Meta_Ref_Sctype/Meta")


#Read in CSV files:
#velasco and polio
df1 <- read.csv("proportion_meta_velasco_polio.csv")
df1 <- df1 %>% rename("Source" = Group)
df1$Source[df1$Source == "3months"] <- "Velasco_12w"
df1$Source[df1$Source == "6months"] <- "Velasco_24w"
df1$Source[df1$Source == "17-18 GW"] <- "Polioudakis_17-18GW"

df2 <- read.csv("proportion_herring_meta_annot.csv")

df2 <- df2 %>% rename("CellType" = X)
df2 <- df2 %>% rename("Source" = Group)
df2$Source[df2$Source == "22 GA"] <- "Herring_22GW"
df2$Source[df2$Source == "24 GA"] <- "Herring_24GW"
df2$Source[df2$Source == "34 GA"] <- "Herring_34GW"

#Comdine df1 and df2
#remove columns that aren't the same
df1[c("X", "Count")] <- list(NULL)
df2[c("major_clust")] <- list(NULL)

df <- rbind(df1, df2)

#change name of df
Meta <- df
Meta <- Meta %>%
  mutate(group = case_when(
    Source == "Herring_22GW" ~ "Herring",
    Source == "Herring_24GW" ~ "Herring",
    Source == "Herring_34GW" ~ "Herring",
    Source == "Velasco_12w" ~ "Velasco",      
    Source == "Velasco_24w" ~ "Velasco",
    Source == "Polioudakis_17-18GW" ~ "Polioudakis",
    TRUE ~ NA_character_                     # Default case
  ))
```

## Make a VennDiagram to visualize the overlap in cell labels

``` r
# Filter by the group and then get unique CellType terms
herr_terms <- unique(Meta[Meta$group == "Herring",]$CellType)
pol_terms <- unique(Meta[Meta$group == "Polioudakis",]$CellType)
vel_terms <- unique(Meta[Meta$group == "Velasco",]$CellType)
# List of terms
listInput <- list(Herring = herr_terms, Polioudakis = pol_terms, Velasco = vel_terms)


venn.plot <- venn.diagram(
  x = listInput,
  category.names = c("Herring", "Polioudakis", "Velasco"),
  category.colours = c("red", "blue", "green"),
  filename = NULL,
    # Circles
  lwd = 1,
  lty = 'blank',
  fill = c(alpha("blue",0.3), alpha('red',0.3), alpha('green',0.3)),
  
  # Numbers
  cex = 1.0,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1.0,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-160, 0, 160),  # Adjusted angles
  cat.dist = c(0.07, 0.07, 0.07), # Slightly increased distances to ensure clarity
  cat.fontfamily = "sans",
  rotation = 1,
  
  # Adjust margin
  margin = 0.2  # Adjust as needed
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)
```

![](meta_plots_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Make a bar chart for author’s annotations

``` r
df$Source<- factor(Meta$Source, levels = c("Polioudakis_17-18GW", "Herring_22GW", "Herring_24GW", "Herring_34GW", 
                                         "Velasco_12w","Velasco_24w", ordered = TRUE))
                   


custom_colors <- c("Polioudakis_17-18GW"="#402C49",
                    "Herring_22GW"="#F4A460", "Herring_24GW"="#D2691E", "Herring_34GW"="#A0522D", 
                   "Velasco_12w"="#DEB887",
                    "Velasco_24w" = "#D2B48C")


#set all cell types that are not present in each dataset to 0
df_all_combinations <- expand.grid(Source = unique(Meta$Source), CellType = unique(Meta$CellType))
# Merge the new data frame with the original data frame
combo_full <- merge(df_all_combinations, Meta, by = c("Source", "CellType"), all.x = TRUE)

combo_full <- replace(combo_full, is.na(combo_full), 0)


ggplot(combo_full, aes(fill = Source, y = Percentage, x = CellType)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8) +
  labs(title="Author's Annotation", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 6)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) 
```

![](meta_plots_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Make an Upset plot for the author’s annotations

``` r
library(UpSetR)

# Create a matrix from distinct group-CellType combinations
mat <- Meta %>%
  distinct(group, CellType) %>%
  group_by(group, CellType) %>%
  tally() %>%
  # Pivoting data to wide format for the UpSetR matrix
  pivot_wider(names_from = group, values_from = n, values_fill = 0)

# Extract unique CellType values from the created matrix
row_names <- mat$CellType

# Assign row names for the matrix from the original dataframe
mat <- as.matrix(mat[, -1])
rownames(mat) <- row_names

# Convert matrix to dataframe for UpSetR 
df_for_upset <- as.data.frame(mat)

# Plot the UpSetR plot
upset(df_for_upset)
```

![](meta_plots_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
