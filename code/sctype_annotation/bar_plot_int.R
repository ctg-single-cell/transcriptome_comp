#Created by Fallon Ratner on 13.07.23
#Bar Plots for Annotating Integrated Data
library(ggplot2)
library(dplyr)
library(stringr)
#In Vivo Integrated:76,212 cells
#In Vitro Integrated: 66,295 cells
# Read the CSV files
setwd("C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/snellius/output")
vivo <- read.csv("outputvivo_int_sctype.csv")
vitro <- read.csv("outputvitro_int_sctype.csv")

# Calculate the percentage for each subset
vivo <- vivo %>%
  group_by(source, cell_type) %>%
  summarize(total_cells = sum(ncells), .groups = "drop") %>%
  group_by(source) %>%
  mutate(percentage = (total_cells / sum(total_cells)) * 100)

#the following code is to make sure all of the bars are the same width
df_all_combinations <- expand.grid(source = unique(vivo$source), cell_type = unique(vivo$cell_type))
# Merge the new data frame with the original data frame
combo_full <- merge(df_all_combinations, vivo, by = c("source", "cell_type"), all.x = TRUE)

combo_full <- replace(combo_full, is.na(combo_full), 0)
# Convert the 'group' column into an ordered factor
combo_full$source<- factor(combo_full$source, levels = c("Han_GW11.1","Han_GW11.2", "Han_GW12", "Han_GW13", "Couturier_GW13", 
                                                       "Couturier_GW17", "Couturier_GW19","Polioudakis_GW17-18", "Liu_GW17-19",
                                                       "Herring_GA22", "Herring_GA24", "Herring_GA34"), ordered = TRUE)
# Create the ggplot bar chart
pdf("vivo_int_v1.pdf")
barviv <- ggplot(combo_full, aes(fill = source, y = percentage, x = cell_type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

barviv + labs(title="scType Annotation of Integrated In Vivo Brain", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) +scale_fill_manual(values=c( 
  "Han_GW11.1"="#D3D3D3","Han_GW11.2"=	"#A9A9A9", "Han_GW12"="#808080", "Han_GW13"= "#696969", "Couturier_GW13"="#FFEFD5", 
  "Couturier_GW17"="#FFE4B5", "Couturier_GW19"="#FFDEAD","Polioudakis_GW17-18"="#4B0082", "Liu_GW17-19"="#0000CD",
  "Herring_GA22"="#F4A460", "Herring_GA24"="#D2691E", "Herring_GA34"="#A0522D"))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 9)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

pdf("vivo_int_v2.pdf")
barviv2 <- ggplot(combo_full, aes(fill = source, y = total_cells, x = cell_type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

barviv2 + labs(title="scType Annotation of Integrated In Vivo Brain", y = "# of Cells", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) +scale_fill_manual(values=c( 
    "Han_GW11.1"="#D3D3D3","Han_GW11.2"=	"#A9A9A9", "Han_GW12"="#808080", "Han_GW13"= "#696969", "Couturier_GW13"="#FFEFD5", 
    "Couturier_GW17"="#FFE4B5", "Couturier_GW19"="#FFDEAD","Polioudakis_GW17-18"="#4B0082", "Liu_GW17-19"="#0000CD",
    "Herring_GA22"="#F4A460", "Herring_GA24"="#D2691E", "Herring_GA34"="#A0522D"))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 9)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()
#########################################################################
#########################################################################
vitro <- vitro %>%
  group_by(source, cell_type) %>%
  summarize(total_cells = sum(ncells), .groups = "drop") %>%
  group_by(source) %>%
  mutate(percentage = (total_cells / sum(total_cells)) * 100)

#the following code is to make sure all of the bars are the same width
df_all_combinations2 <- expand.grid(source = unique(vitro$source), cell_type = unique(vitro$cell_type))
# Merge the new data frame with the original data frame
combo_full2 <- merge(df_all_combinations2, vitro, by = c("source", "cell_type"), all.x = TRUE)

combo_full2 <- replace(combo_full2, is.na(combo_full2), 0)
# Convert the 'group' column into an ordered factor
combo_full2$source <- factor(combo_full2$source, levels = c("Trujillo_1mo","Meng_U1M_1.5mo", "Meng_U2F_1.5mo", "Giandomenico_2.5mo", "Trujillo_3mo", 
                                                       "Velasco_3mo", "Velasco_6mo","Trujillo_6mo", "Trujillo_10mo"), ordered = TRUE)
# Create the ggplot bar chart
pdf("vitro_int_v1.pdf")
barvit <- ggplot(combo_full2, aes(fill = source, y = percentage, x = cell_type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

barvit + labs(title="scType Annotation of Integrated In Vitro Brain Organoids", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) +scale_fill_manual(values=c( 
    "Trujillo_1mo"="#3CB371","Meng_U1M_1.5mo"=	"#DAA520", "Meng_U2F_1.5mo"="#B8860B", "Giandomenico_2.5mo"= "#DC143C", "Trujillo_3mo"="#2E8B57", 
    "Velasco_3mo"="#5F9EA0", "Velasco_6mo"="#008080","Trujillo_6mo"="#008000", "Trujillo_10mo"="#006400"))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 9)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

# Convert the 'group' column into an ordered factor
combo_full2$source <- factor(combo_full2$source, levels = c("Meng_U1M_1.5mo", "Meng_U2F_1.5mo", "Giandomenico_2.5mo",  
                                                            "Velasco_3mo", "Velasco_6mo","Trujillo_1mo","Trujillo_3mo","Trujillo_6mo", "Trujillo_10mo"), ordered = TRUE)
# Create the ggplot bar chart
pdf("vitro_int_v2.pdf")
barvit2 <- ggplot(combo_full2, aes(fill = source, y = percentage, x = cell_type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

barvit2 + labs(title="scType Annotation of Integrated In Vitro Brain Organoids", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) +scale_fill_manual(values=c( 
    "Trujillo_1mo"="#3CB371","Meng_U1M_1.5mo"=	"#DAA520", "Meng_U2F_1.5mo"="#B8860B", "Giandomenico_2.5mo"= "#DC143C", "Trujillo_3mo"="#2E8B57", 
    "Velasco_3mo"="#5F9EA0", "Velasco_6mo"="#008080","Trujillo_6mo"="#008000", "Trujillo_10mo"="#006400"))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 9)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

pdf("vitro_int_v3.pdf")
barvit3 <- ggplot(combo_full2, aes(fill = source, y = total_cells, x = cell_type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

barvit3 + labs(title="scType Annotation of Integrated In Vitro Brain Organoids", y = "# of Cells", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) +scale_fill_manual(values=c( 
    "Trujillo_1mo"="#3CB371","Meng_U1M_1.5mo"=	"#DAA520", "Meng_U2F_1.5mo"="#B8860B", "Giandomenico_2.5mo"= "#DC143C", "Trujillo_3mo"="#2E8B57", 
    "Velasco_3mo"="#5F9EA0", "Velasco_6mo"="#008080","Trujillo_6mo"="#008000", "Trujillo_10mo"="#006400"))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 9)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

#########################################################################
#########################################################################
# Combine data frames into one
combo <- rbind(vivo, vitro)
#the following code is to make sure all of the bars are the same width
df_all_combinations3 <- expand.grid(source = unique(combo$source), cell_type = unique(combo$cell_type))
# Merge the new data frame with the original data frame
combo_full3 <- merge(df_all_combinations3, combo, by = c("source", "cell_type"), all.x = TRUE)

combo_full3 <- replace(combo_full3, is.na(combo_full3), 0)
# Convert the 'group' column into an ordered factor
combo_full3$source <- factor(combo_full3$source, levels = c("Han_GW11.1","Han_GW11.2", "Han_GW12", "Han_GW13", "Couturier_GW13", 
                                                            "Couturier_GW17", "Couturier_GW19","Polioudakis_GW17-18", "Liu_GW17-19",
                                                            "Herring_GA22", "Herring_GA24", "Herring_GA34","Meng_U1M_1.5mo", "Meng_U2F_1.5mo", "Giandomenico_2.5mo",  
                                                            "Velasco_3mo", "Velasco_6mo","Trujillo_1mo","Trujillo_3mo","Trujillo_6mo", "Trujillo_10mo"), ordered = TRUE)


# Define custom color palette
custom_colors <- c(
  "Han_GW11.1"="#D3D3D3","Han_GW11.2"=	"#A9A9A9", "Han_GW12"="#808080", "Han_GW13"= "#696969", "Couturier_GW13"="#FFEFD5", 
  "Couturier_GW17"="#FFE4B5", "Couturier_GW19"="#FFDEAD","Polioudakis_GW17-18"="#4B0082", "Liu_GW17-19"="#0000CD",
  "Herring_GA22"="#F4A460", "Herring_GA24"="#D2691E", "Herring_GA34"="#A0522D",
  "Trujillo_1mo"="#3CB371","Meng_U1M_1.5mo"=	"#DAA520", "Meng_U2F_1.5mo"="#B8860B", "Giandomenico_2.5mo"= "#DC143C", "Trujillo_3mo"="#2E8B57", 
  "Velasco_3mo"="#5F9EA0", "Velasco_6mo"="#008080","Trujillo_6mo"="#008000", "Trujillo_10mo"="#006400"
)

# Create separate bar plots for each cell type and save as PDF
cell_types <- unique(combo_full3$cell_type)
for (cell_type in cell_types) {
  plot_data <- combo_full3[combo_full3$cell_type == cell_type, ]
  
  p <- ggplot(plot_data, aes(fill = source, y = percentage, x = cell_type,pattern = new_column)) +
    geom_bar(position = "dodge", stat = "identity", width = 0.8) +
    labs(title = paste("scType Annotation of Integrated In vivo vs In vitro Brain", cell_type),
         y = "Proportion (%)",
         fill = "scRNAseq Data",
         pattern = "Vivo/Vitro") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text = element_text(size = 9),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = custom_colors) + 
    geom_bar_pattern(position = position_dodge(preserve = "single"),
                                                                color = "black", 
                                                                pattern_fill = "black",
                                                                pattern_angle = 45,
                                                                pattern_density = 0.1,
                                                                pattern_spacing = 0.025,
                                                                pattern_key_scale_factor = 0.6) + 
    scale_pattern_manual(values = c(Vitro = "stripe", Vivo = "none")) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none")))
  # Save the plot as PDF
  pdf_file <- paste0("barplot1_", cell_type, ".pdf")
  ggsave(filename = pdf_file, plot = p, width = 6, height = 4)
}


  
  
  

