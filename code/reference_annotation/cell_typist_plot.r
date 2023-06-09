# Project: Transcriptome Comp
# Created by: Fallon Ratner 


library(ggplot2)
library(dplyr)

#########################################################################################################################
#Plot generation for Cell_Typist Output
df <- read.csv("proportion_17-18 GW.csv")

#Subset df
df1 <- df[1:10,]

ggplot(subset(df1, Group == "17-18 GW"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 17-18 GW Cortex", subtitle = "Polioudakis et al., 2019", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "blue") + 
  scale_color_manual(values = "blue")

#Subset df
df2 <- df[10:105,]

ggplot(subset(df2, Group == "17-18 GW"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 17-18 GW Cortex", subtitle = "Polioudakis et al., 2019", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(size = 8),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "blue") + 
  scale_color_manual(values = "blue")
#########################################################################################################################
#Plot generation for Cell_Typist Output
df <- read.csv("proportion_22_GA.csv")

#Subset df
df1 <- df[1:10,]

#Subset df
df2 <- df[10:105,]


ggplot(subset(df1, Group == "22 GA"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 22 GA Cortex", subtitle = "Herring et al., 2022", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "#D3D3D3") + 
  scale_color_manual(values = "#D3D3D3")

ggplot(subset(df2, Group == "22 GA"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 22 GA Cortex", subtitle = "Herring et al., 2022", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(size = 8),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "#D3D3D3") + 
  scale_color_manual(values = "#D3D3D3")
#########################################################################################################################
#Plot generation for Cell_Typist Output
df <- read.csv("proportion_24_GA.csv")

#Subset df
df1 <- df[1:10,]

#Subset df
df2 <- df[10:105,]


ggplot(subset(df1, Group == "24_GA"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 24 GA Cortex", subtitle = "Herring et al., 2022", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "dark gray") + 
  scale_color_manual(values = "dark gray")

ggplot(subset(df2, Group == "24_GA"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 24 GA Cortex", subtitle = "Herring et al., 2022", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(size = 8),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "#A9A9A9") + 
  scale_color_manual(values = "#A9A9A9")
#########################################################################################################################
#Plot generation for Cell_Typist Output
df <- read.csv("proportion_34_GA.csv")

#Subset df
df1 <- df[1:10,]

#Subset df
df2 <- df[10:105,]


ggplot(subset(df1, Group == "34_GA"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 24 GA Cortex", subtitle = "Herring et al., 2022", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "#696969") + 
  scale_color_manual(values = "#696969")

ggplot(subset(df2, Group == "34_GA"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 24 GA Cortex", subtitle = "Herring et al., 2022", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(size = 8),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "#696969") + 
  scale_color_manual(values = "#696969")

#########################################################################################################################
#Plot generation for Cell_Typist Output
df <- read.csv("proportion_3 months.csv")

#Subset df
df1 <- df[1:10,]

#Subset df
df2 <- df[10:105,]


ggplot(subset(df1, Group == "3 months"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 3 months Organoids", subtitle = "Velasco et al., 2019 ", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "lightgray") + 
  scale_color_manual(values = "lightgray")

ggplot(subset(df2, Group == "3 months"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 3 months Organoids", subtitle = "Velasco et al., 2019 ", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(size = 8),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "lightgray") + 
  scale_color_manual(values = "lightgray")
#########################################################################################################################
#Plot generation for Cell_Typist Output
df <- read.csv("proportion_6 months.csv")

#Subset df
df1 <- df[1:10,]

#Subset df
df2 <- df[10:105,]


ggplot(subset(df1, Group == "6 months"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 6 months Organoids", subtitle = "Velasco et al., 2019 ", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "darkgray") + 
  scale_color_manual(values = "darkgray")

ggplot(subset(df2, Group == "6 months"), aes(y = Percentage, x = reorder(predicted_labels, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 6 months Organoids", subtitle = "Velasco et al., 2019 ", y = "Proportion (%)", x = NULL, fill ="scRNAseq Data") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(size = 8),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "darkgray") + 
  scale_color_manual(values = "darkgray")
