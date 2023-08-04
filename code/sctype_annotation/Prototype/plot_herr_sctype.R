# Project: Transcriptome Comp
# Created by: Fallon Ratner on 01.05.23

library(ggplot2)
library(dplyr)
library(stringr)

#Read in the files
library(dplyr)

# Read the CSV files
ga22 <- read.csv("ga22_cluster_sctype.csv")
ga24 <- read.csv("ga24_cluster_sctype.csv")
ga34 <- read.csv("ga34_cluster_sctype.csv")

# Add age column to each data frame
ga22 <- ga22 %>% mutate(age = "ga22")
ga24 <- ga24 %>% mutate(age = "ga24")
ga34 <- ga34 %>% mutate(age = "ga34")

# Combine data frames into one
herr <- rbind(ga22, ga24, ga34)

# Calculate the percentage for each subset
herr <- herr %>%
  group_by(age, type) %>%
  summarize(total_cells = sum(ncells), .groups = "drop") %>%
  group_by(age) %>%
  mutate(percentage = (total_cells / sum(total_cells)) * 100)

#the following code is to make sure all of the bars are the same width
df_all_combinations <- expand.grid(age = unique(herr$age), type = unique(herr$type))
# Merge the new data frame with the original data frame
combo_full <- merge(df_all_combinations, herr, by = c("age", "type"), all.x = TRUE)

combo_full <- replace(combo_full, is.na(combo_full), 0)
# Convert the 'group' column into an ordered factor
combo_full$age <- factor(combo_full$age, levels = c("ga22","ga24", "ga34"), ordered = TRUE)

# Create the ggplot bar chart
bar <- ggplot(combo_full, aes(fill = age, y = percentage, x = type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

bar + labs(title="scType Annotation", subtitle="Herring et al., 2022", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values=c( 
    "ga22" = "#D3D3D3", 
    "ga24" = "#A9A9A9",
    "ga34"= "#696969"))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 9)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()