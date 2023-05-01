# Project: Transcriptome Comp
# Created by: Fallon Ratner on 01.05.23


library(ggplot2)
library(dplyr)

#Read in the file with the metatdata for Velasco and Polioudakis
df <- read.csv("proportion_meta_velasco_polio.csv")
df <- df[, -1]
df_all_combinations <- expand.grid(Group = unique(df$Group), CellType = unique(df$CellType))

# Merge the new data frame with the original data frame
combo_full <- merge(df_all_combinations, df, by = c("Group", "CellType"), all.x = TRUE)

combo_full <- replace(combo_full, is.na(combo_full), 0)

# Create ggplot bar chart
pdf("vel_pol_meta_prop_unstackv3.pdf")
bar <- ggplot(combo_full, aes(fill = Group, y = Percentage, x = CellType)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8) +
  labs(title = "Percentage of Cells Meta Annotation",
       subtitle = "Velasco et al., 2019 & Polioudakis et al., 2019",
       y = "Proportion (%)", fill = "scRNAseq Data") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 8),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = c("17-18 GW" = "blue", 
                               "3months" = "lightgray", 
                               "6months" = "darkgray")) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15))
dev.off()
#######################################################################################
# Read the CSV file for Velasco & Poliodaskis cell counts & plot separately
combo1 <- read.csv("proportion_velasco_polio.csv") %>%
  select(-1)

ggplot(subset(combo1, Group == "3months"), aes(y = Count, x = reorder(CellType, -Count), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Count 3months Organoid", subtitle = "Velasco et al., 2019", y = "Count", x = NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "#3399FF") + 
  scale_color_manual(values = "#3399FF") +
  ggsave("count_velasco_3months_v2.pdf")

ggplot(subset(combo1, Group == "6months"), aes(y = Count, x = reorder(CellType, -Count), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Count 6months Organoid", subtitle = "Velasco et al., 2019", y = "Count", x = NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "#000099") + 
  scale_color_manual(values = "#000099") +
  ggsave("count_velasco_6months_v2.pdf")

ggplot(subset(combo1, Group == "17-18 GW"), aes(y = Count, x = reorder(CellType, -Count), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Count 17-18 GW Cortex", subtitle = "Polioudakis et al., 2019", y = "Count", x = NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "#FF3333") + 
  scale_color_manual(values = "#FF3333") +
  ggsave("count_polio_v2.pdf")

########################################################################
#Make plots for the cells present in Velasco & Polioudakis
#Subset df for oRG only
org <- combo1[grep("oRG", combo1$CellType), ]
head(org)
#Make a bar plot for org and save as pdf
pdf("prop_velasco_polio_org_v2.pdf")
bar2 <- ggplot(org, aes(y=Percentage, x=Group, fill=Group)) + 
    geom_col()
    
bar2 + labs(title="Percentage of Outer Radial Glia", subtitle="Velasco et al., 2019 & Polioudakis et al., 2019", y = "Proportion (%)") + theme(axis.title.x = element_blank()) + theme(legend.position = "none")  + scale_fill_manual(values=c("17-18 GW" ="#FF3333","3months" = "#3399FF","6months" = "#000099"))+ scale_color_manual(values = c("17-18 GW" = "#FF3333", "3months" = "#3399FF", "6months" = "#000099")) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))  
dev.off()
########################################################################
#Read in df for cycling only
cy <- read.csv("proportion_velasco_polio_cycling.csv") 
#remove first column
cy <- cy[, -1]
#Make a bar plot for cycling and save as pdf
pdf("prop_velasco_polio_cy_v2.pdf")
bar3 <- ggplot(cy, aes(y=Percentage, x=Group, fill=Group)) + 
    geom_col()
    
bar3 + labs(title="Percentage of Cycling Cells", subtitle="Velasco et al., 2019 & Polioudakis et al., 2019", y = "Proportion (%)") + theme(axis.title.x = element_blank()) + theme(legend.position = "none")  + scale_fill_manual(values=c("17-18 GW" ="#FF3333","3months" = "#3399FF","6months" = "#000099"))+ scale_color_manual(values = c("17-18 GW" = "#FF3333", "3months" = "#3399FF", "6months" = "#000099")) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))   
dev.off()
########################################################################
#Subset df for IPs only
ip <- combo1[grep("IP", combo1$CellType), ]
head(ip)
#Make a bar plot for org and save as pdf
pdf("prop_velasco_polio_IP_v2.pdf")
bar4 <- ggplot(ip, aes(y=Percentage, x=Group, fill=Group)) + 
    geom_col()
    
bar4 + labs(title="Percentage of Intermediate Progenitors", subtitle="Velasco et al., 2019 & Polioudakis et al., 2019", y = "Proportion (%)") + theme(axis.title.x = element_blank()) + theme(legend.position = "none")  + scale_fill_manual(values=c("17-18 GW" ="#FF3333","3months" = "#3399FF","6months" = "#000099"))+ scale_color_manual(values = c("17-18 GW" = "#FF3333", "3months" = "#3399FF", "6months" = "#000099")) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()
#########################################################################################################
#Make plots for the Velasco & Polioudakis Cell Annotations made w/ scType
combo <- read.csv("proportion_velasco_polio_sctype_clusterv2.csv")
ggplot(subset(combo1, Group == "3months"), aes(y = Percentage, x = reorder(CellType, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 3months Organoid", subtitle = "Velasco et al., 2019", y = "Percentage", x = NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "#3399FF") + 
  scale_color_manual(values = "#3399FF") +
  ggsave("Percentage_velasco_3months_v2.pdf")

ggplot(subset(combo1, Group == "6months"), aes(y = Percentage, x = reorder(CellType, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 6months Organoid", subtitle = "Velasco et al., 2019", y = "Percentage", x = NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "#000099") + 
  scale_color_manual(values = "#000099") +
  ggsave("Percentage_velasco_6months_v2.pdf")

ggplot(subset(combo1, Group == "17-18 GW"), aes(y = Percentage, x = reorder(CellType, -Percentage), fill = Group)) + 
  geom_col() +
  labs(title = "Cell Percentage 17-18 GW Cortex", subtitle = "Polioudakis et al., 2019", y = "Percentage", x = NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = "#FF3333") + 
  scale_color_manual(values = "#FF3333") +
  ggsave("Percentage_polio_v2.pdf")
#########################################################################################################################
#Plot generation for scType Output
combo <- read.csv("proportion_velasco_polio_sctype_cluster.csv")
combo <- combo[, -1]
#the following code is to make sure all of the bars are the same width
df_all_combinations <- expand.grid(Group = unique(combo$Group), type = unique(combo$type))
# Merge the new data frame with the original data frame
combo_full <- merge(df_all_combinations, combo, by = c("Group", "type"), all.x = TRUE)

combo_full <- replace(combo_full, is.na(combo_full), 0)


library(stringr)
pdf("sctype_cluster_prop_unstackv2.pdf")
# Create the ggplot bar chart
bar5 <- ggplot(combo_full, aes(fill = Group, y = Percentage, x = type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

bar5 + labs(title="Percentage of scType Annotation", subtitle="Velasco et al., 2019 & Polioudakis et al., 2019", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values=c("17-18 GW" = "blue", 
                                                                     "3months" = "lightgray", 
                                                                     "6months" = "darkgray"))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()
#################################################################################
pdf("sctype_cluster_prop_unstackv3.pdf")
# Create the ggplot bar chart
bar5 <- ggplot(combo_full, aes(fill = Group, y = Percentage, x = type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

bar5 + labs(title="Percentage of scType Annotation", subtitle="Velasco et al., 2019 & Polioudakis et al., 2019", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values=c("blue",
                                                                     "white",
                                                                     "white"))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()
##########################################
pdf("sctype_cluster_prop_unstackv4.pdf")
# Create the ggplot bar chart
bar5 <- ggplot(combo_full, aes(fill = Group, y = Percentage, x = type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

bar5 + labs(title="Percentage of scType Annotation", subtitle="Velasco et al., 2019 & Polioudakis et al., 2019", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values=c("blue",
                                                                     "lightgray",
                                                                     "white"))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()