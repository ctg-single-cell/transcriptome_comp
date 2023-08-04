#Created by Fallon Ratner on02.08.23
#Bar Plots for Annotating Integrated Data
library(ggplot2)
library(dplyr)
library(stringr)

#In Vitro Integrated:  cells
# Read the CSV files
vitro <- read.csv("outputvitro_int_sctype.csv")
vitro_gb <- read.csv("outputvitro_gaba_int_sctype.csv")

#Calculate the proportion for each cell type per source as a percentage
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

#the following code is to make sure all of the bars are the same width
df_all_combinations <- expand.grid(source = unique(vitro$source), cell_type = unique(vitro$cell_type))
# Merge the new data frame with the original data frame
combo_full <- merge(df_all_combinations, vitro, by = c("source", "cell_type"), all.x = TRUE)

combo_full <- replace(combo_full, is.na(combo_full), 0)
# Convert the 'group' column into an ordered factor
combo_full$source <- factor(combo_full$source, levels = c("Bhaduri_H1S_3w","Bhaduri_H1X_3w", "Bhaduri_H1S_5w","Bhaduri_H1X_5w", "Bhaduri_H1S_8w",
                                                            "Bhaduri_H1X_8w","Bhaduri_H1S_10w", "Bhaduri_H1X_10w",
                                                            "Trujillo_4w", "Trujillo_12w", "Trujillo_24w", "Trujillo_40w",
                                                            "Xiangb1_4w", "Xiangb2_4w","Xiang_10w","Xiang_11w",
                                                            "Fair_12w", "Fair_20w","Velasco_12w","Velasco_24w", 
                                                            "Popova_batch1_7w", "Popova_batch2_7w",
                                                            "Giandomenico_10w", "Madhavan_12w"), ordered = TRUE)

# Define custom color palette
custom_colors <- c(
  "Bhaduri_H1S_3w"="#87CEFA","Bhaduri_H1X_3w"=	"#87CEEB", "Bhaduri_H1S_5w"="#6495ED", "Bhaduri_H1X_5w"= "#1E90FF", "Bhaduri_H1S_8w"="#0000FF", 
  "Bhaduri_H1X_8w"="#0000CD", "Bhaduri_H1S_10w"="#00008B","Bhaduri_H1X_10w"="#000080", "Fair_12w"="#40E0D0",
  "Fair_20w"="#48D1CC", "Giandomenico_10w"="#008080", "Madhavan_12w"="#4B0082",
  "Popova_batch1_7w"="#8B008B", "Popova_batch2_7w"= "#800080", "Trujillo_4w"="#FF6347", 
  "Trujillo_12w"="#DC143C", "Trujillo_24w"="#FF4500","Trujillo_40w"="#FF0000", "Velasco_12w"="#DAA520",
  "Velasco_24w" = "#B8860B","Xiangb1_4w"="#D2B48C", "Xiangb2_4w"="#DEB887","Xiang_10w" = "#A0522D","Xiang_11w"="#8B4513")
#Plot everything together
pdf("vitro_int.pdf")
barvit <- ggplot(combo_full, aes(fill = source, y = percentage, x = cell_type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

barvit + labs(title="scType Annotation of Integrated In Vitro Brain Organoids", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 9)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()


# Create separate bar plots for each cell type and save as PDF
cell_types <- unique(combo_full$cell_type)
for (cell_type in cell_types) {
  plot_data <- combo_full[combo_full$cell_type == cell_type, ]
  
  p <- ggplot(plot_data, aes(fill = source, y = percentage, x = cell_type)) +
    geom_bar(position = "dodge", stat = "identity", width = 0.8) +
    labs(title = paste("scType Annotation of Integrated Brain Organoids", cell_type),
         y = "Proportion (%)",
         fill = "scRNAseq Data") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text = element_text(size = 9),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = custom_colors)
  # Save the plot as PDF
  pdf_file <- paste0("barplot1_", cell_type, ".pdf")
  ggsave(filename = pdf_file, plot = p, width = 6, height = 4)
}

#########################################################################3
#the following code is to make sure all of the bars are the same width
df_all_combinations2 <- expand.grid(source = unique(vitro_gb$source), cell_type = unique(vitro_gb$cell_type))
# Merge the new data frame with the original data frame
combo_full2 <- merge(df_all_combinations2, vitro_gb, by = c("source", "cell_type"), all.x = TRUE)

combo_full2 <- replace(combo_full2, is.na(combo_full2), 0)
combo_full2$source <- factor(combo_full2$source, levels = c("Bhaduri_H1S_3w","Bhaduri_H1X_3w", "Bhaduri_H1S_5w","Bhaduri_H1X_5w", "Bhaduri_H1S_8w",
                                                          "Bhaduri_H1X_8w","Bhaduri_H1S_10w", "Bhaduri_H1X_10w",
                                                          "Trujillo_4w", "Trujillo_12w", "Trujillo_24w", "Trujillo_40w",
                                                          "Xiangb1_4w", "Xiangb2_4w","Xiang_10w","Xiang_11w",
                                                          "Fair_12w", "Fair_20w","Velasco_12w","Velasco_24w", 
                                                          "Popova_batch1_7w", "Popova_batch2_7w",
                                                          "Giandomenico_10w", "Madhavan_12w"), ordered = TRUE)

# Create the ggplot bar chart
pdf("vitro_gaba_int.pdf")
barvitgb <- ggplot(combo_full2, aes(fill = source, y = percentage, x = cell_type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

barvitgb + labs(title="scType Annotation of Integrated In Vitro Brain Organoids: Interneurons", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) +scale_fill_manual(values= custom_colors)+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 9)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()
