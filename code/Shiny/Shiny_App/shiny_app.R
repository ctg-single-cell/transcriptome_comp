#R Shiny App
#Created by Fallon Ratner 14-09-23

library(ggplot2)
library(grid)
library(stringr)
library(dplyr)
library(shiny)
library(plotly)
library(Seurat)
library(DT)
setwd("C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/Shiny")

# Define UI ----
ui <- fluidPage(
  titlePanel("Comparing scRNAseq from Organoids and Fetal Brain"),
  sidebarLayout(
    sidebarPanel(
      selectInput("cell_type", "Select Cell Type:", 
                  choices = c("ventral Radial Glia",
                               "outer radial glia",
                               "Intermediate Progenitors",
                               "Neuroblasts",
                               "Glioblasts",
                               "Upper Excitatory Layers",
                               "Deep Excitatory Layers",
                               "Immature Excitatory",
                               "Maturing Excitatory",
                               "Migrating Excitatory",
                               "Interneuron Precursors",
                               "MGE INs",
                               "CGE INs",
                               "PVALB",
                               "SST",
                               "VIP",
                               "ID2",
                               "NDNF",
                               "OPCs",
                               "Oligodendrocytes",
                               "Immature Astrocytes",
                               "Mature Astrocytes",
                               "Microglia",
                               "Endothelial Cells",
                               "Mural")),
      selectInput("dev_stage", "Select Developmental Stage:", 
                  choices = c("Early Development (11-13GW)", 
                                   "Mid Development (17-19GW)",
                                   "Late Development (22GW-Day2)")),
      actionButton("goButton", "Go!")
    ),
    # Create a tabset panel
    mainPanel(
      tabsetPanel(
       # First tab
        tabPanel("Tab 1", 
               tags$p("This is a tool to compare cell types from various cortical organoid protocols to the fetal brain. 
                      Multiple public scRNAseq datasets from  embryonic human brain were processed and annotated separately.
                      Multiple public scRNAseq datsets from cortical organoids were integrated and annotated together.
                      Details about the datasets can be found in Tab3. Information
                      about the code can be found here: https://github.com/ctg-single-cell/transcriptome_comp"),
               h3("Organoid Protocol Selection Guide"),
               img(src = "shiny_description.jpg", height = 150, width = 600),
               # In Vivo Data section
               h3("In Vivo Data"),
               img(src = "vivo_umap_cc_2209.svg", height = 600, width = 600),
               img(src = "vivo_gaba_umap_cc_2209.svg",height = 600, width = 600),
               # In Vitro Data section
               h3("In Vitro Data"),
               img(src = "vitro_umap3_int_2109.svg", height = 600, width = 600),
               img(src="vitro_gaba_umap_int_2209.svg", height = 600, width = 600)
       ),
      
      # Second tab
      tabPanel("Tab 2", 
               h3("Correlation Between Cell Types: Fetal Brain vs Organoids"),
               DTOutput("correlationTable"),
               plotOutput("correlationHeatmap")
      ),
      
      # Add more tabs as needed
      tabPanel("Tab 3", 
               h3("Key Resources"),
               dataTableOutput("keyTable"),
               h3("In Vivo Datasets"),
               dataTableOutput("vivoTable"),
               h3("In Vitro Datasets"),
               dataTableOutput("vitroTable")
      )
    ) # Closing for tabsetPanel()
    ) # End of mainPanel
  ) # Closing for sidebarLayout()
) # Closing for fluidPage()

############ Define tables for tab 3
key <- data.frame(
  `Deposited Data` = c("scRNAseq of 13-19 GW Cortex", "scRNAseq of 11-13 GW brain","snRNAseq of 22GW-Day2 PFC",
                       "scRNAseq of 17-19 GW Cortex", "scRNAseq of 17-18 GW Cortex", "scRNAseq of hCOs",
                       "scRNAseq of hCOs", "scRNAseq of hCOs", "scRNAseq of hCOs", "scRNAseq of hCOs",
                       "scRNAseq of hCOs", "scRNAseq of hCOs"), 
  Source = c("Couturier et al., 2020", "Han et al., 2020", "Herring et al., 2022", "Liu et al., 2023",
             "Polioudakis et al., 2019", "Bhaduri et al., 2020", "Fair et al., 2020", "Giandomenico et al., 2019",
             "Madhavan et al., 2018", "Trujillo et al., 2019", "Velasco et al., 2019", "Xiang et al., 2017"), 
  Identifier = c("EGA: EGAS00001004422", "GEO: GSE134355","GEO: GSE168408", "BioProject: PRJNA798712",
                 "dbGaP: phs001836", "GEO: GSE132672", "GEO: GSE157019", "GEO: GSE124174", "GEO: GSE110006",
                 "GEO: GSE130238", "GEO: GSE129519", "GEO: GSE98201") 
)

vivo <- data.frame(
  `Dataset` = c("Couturier et al., 2020", "Han et al., 2020", "Herring et al., 2022", "Herring et al., 2022",
                "Herring et al., 2022", "Herring et al., 2022", "Liu et al., 2023", "Polioudakis et al., 2019"),
  Tissue = c("Telencephalon", "Brain", "PFC", "PFC", "PFC", "PFC", "Cortex", "Neocortex"),
  Region = c("N/A", "N/A", "BA9 & BA46", "BA9", "BA9", "BA8", "N/A", "N/A"),
  Age = c("13, 17, 19 GW", "11, 12, 13 GW", "22 GW", "24 GW", "34 GW", "Day 2", "17-19 GW", "17-18 GW"),
  Sequencing = c("Single Cell 3' (Droplet based)", "Single Cell 3' (smart-seq2)",
                 "Chromium Single Cell 3' (Chip B)", "Chromium Single Cell 3' (Chip B)",
                 "Chromium Single Cell 3' (Chip B)", "Chromium Single Cell 3' (Chip B)",
                 "Paired-end (Smart-seq2/3)", "Paired end(Drop-seq)"),
  Genome = c("hg38", "hg38", "hg19", "hg19", "hg19", "hg19", "hg38", "hg38"),
  Sex = c("N/A", "11: M&F;12:M; 13:F", "M", "M", "F", "F", "N/A", "17,17,18: F & 18:M")
)

vitro <- data.frame(
  `Dataset` = c("Bhaduri et al., 2020", "Fair et al., 2020", "Giandomenico et al., 2019",
                "Madhavan et al., 2018", "Popova et al., 2021", "Trujillo et al., 2019",
                "Velasco et al., 2019", "Xiang et al., 2017"),
  Type = c("Cortical: Forebrain", "Cerebral: whole-brain", "Cortical: Forebrain",
                    "Oligo-Cortical Spheroids", "Cortical: neuro-immune model", "Cortical",
                    "Cortical: Forebrain", "Cortical"),
  Genome = c("hg38", "hg38", "hg38", "hg19", "hg38", "hg38", "hg38", "hg19"),
  Sequencing = c("Single Cell 3' v2", "Single Cell 3' v2", "Chromium Single Cell 3' Chip",
                 "Single Cell 3'", "Single Cell 3'", "Single Cell 3’ v2", "Single Cell 3’ v2",
                 "Single Cell 3’ v1 & v2"),
  Protocol = c("Least & Most Directed", "Undirected", "Undirected", "Directed",
               "Directed", "Directed","Directed", "Directed"),
  Cell_Line = c("H1: hES (male)", "iPSC Line A18945", "H1: hES (male) and H9: hES (femael):pooled",
                "H7: hES (female)", "iPSC Lines 28126(male), 1323-4(female), YH10: mixed", 
                "iPSC", "PGP1: iPSC", "HES-3 NKX2-1GFP: hES"),
  Time_Sequenced = c("3,5,8,10 weeks", "93,140 days", "75 days", "12 weeks",
                     "10 weeks", "1,3,6,10 months", "3,6 months", "30,72,79 days")
  
)

# Define the mapping mechanism
dev_stage_map <- list(
  "Early Development (11-13GW)" = "early_devo",
  "Mid Development (17-19GW)" = "mid_devo",
  "Late Development (22GW-Day2)" = "late_devo"
)

# Define server logic
server <- function(input, output) {
  
  output$keyTable <- renderDataTable({
    key
  })
  
  output$vivoTable <- renderDataTable({
    vivo
  })
  
  output$vitroTable <- renderDataTable({
    vitro
  })
  
  # Default outputs
  output$correlationTable <- renderDT({
    datatable(data.frame())
  })
  
  output$correlationHeatmap <- renderPlot({})
  
  observeEvent(input$goButton, {
    special_cell_types <- c("Glioblasts", "Upper Excitatory Layers", "Migrating Excitatory",
                            "PVALB", "VIP", "ID2", "Oligodendrocytes", "Mature Astrocytes",
                            "Microglia")
    
    if (input$cell_type %in% special_cell_types) {
      showModal(modalDialog(
        title = "Information",
        paste("This cell type:", input$cell_type," was not annotated in any of the organoid datasets.")
      ))
    } else {
      input_dev_stage <- input$dev_stage 
      filename_dev_stage <- dev_stage_map[[input_dev_stage]]
      csv_filename <- paste0("pairwise_correlations_", filename_dev_stage, "_vs_datasets_", input$cell_type, ".csv")
      print(csv_filename)
      
      if(file.exists(csv_filename)) {
        correlation_df <- read.csv(csv_filename)
        # Update the data table
        output$correlationTable <- renderDT({
          datatable(correlation_df)
        })
        # Reorder the datasets based on the correlation
        correlation_df$Dataset <- with(correlation_df, reorder(Dataset, Correlation))
        # Update the heatmap
        output$correlationHeatmap <- renderPlot({
          ggplot(correlation_df, aes(x = Dataset, y = "", fill = Correlation)) +
            geom_tile() +
            scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
            theme_minimal() +
            ylab(NULL) +
            # Angle the x-axis labels
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            theme(axis.text.x = element_text(size = 9))
  
        })
      } else {
        showModal(modalDialog(
          title = "Warning",
          "The selected inputs are not available."
        ))
      }
    }
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)