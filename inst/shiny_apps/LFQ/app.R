source(paste(system.file("shiny_apps", "LFQ", package = "DEP2"), "/source.R", sep = ""))

library(shinyWidgets)
library(dplyr)
library(tibble)
library(SummarizedExperiment)
library(DEP2)
library(shiny)
library(shinydashboard)
library(shinyBS)

convertMenuItem <- function(mi,tabName) {
  mi$children[[1]]$attribs['data-toggle']="tab"
  mi$children[[1]]$attribs['data-value'] = tabName
  if(length(mi$attribs$class)>0 && mi$attribs$class=="treeview"){
    mi$attribs$class=NULL
  }
  mi
}

ui <- shinyUI(
  dashboardPage(skin = "blue",
    dashboardHeader(title = "DEP - LFQ"),
    dashboardSidebar(width = 335,
      sidebarMenu(#* sidebarMenu----
        convertMenuItem(#** Welcome menuItem ----
          menuItem("Welcome", selected = FALSE, tabName = "Welcome_tab", icon = icon("th")), "Welcome_tab"
          ),
        convertMenuItem(#** DEP menuItem ----
        menuItem("DEP-LFQ options", selected = TRUE, tabName = "DEP_tab", icon = icon("th"),
                 menuItem("Files", selected = TRUE,
                   fluidRow(
                    #column(width = 2,
                    #actionButton("help_format_DEP",width = 12,
                    #              label = "", icon = icon("bell")
                    #      )),
                   column(width = 9,
                   fileInput('file1',width = "300px",
                           'ProteinGroups.txt',
                           accept=c('text/csv',
                                    'text/comma-separated-values,text/plain',
                                    '.csv')),
                   fileInput('file2',width = "300px",
                             'ExperimentalDesign.txt',
                             accept=c('text/csv',
                                      'text/comma-separated-values,text/plain',
                                      '.csv')),
                   fileInput('resultRData','load saved result RData',width = "300px")
                   ),
                   column(
                     width = 1, h4(),
                     br(),
                     actionButton("help_format_DEP",
                                  label = "", icon = icon("question-circle"),
                                  style = "color: #f6f6f6; background-color: #2c3b41; border-color: #2c3b41"
                     ),
                     shinyBS::bsTooltip(
                       "help_format_DEP",
                       "How to provide your input data",
                       "top",
                       options = list(container = "body")
                     )
                   )),
                  radioButtons("anno",
                               "Sample annotation",
                               choices = list("Parse from columns" = "columns",
                                              "Use Experimental Design" = "expdesign"),
                               selected = "columns")
                  ),
                 shinyBS::bsTooltip("anno", "Sample annotation type", "top", options = list(container = "body")),
                 # menuItemOutput("columns"),
                 menuItem("Columns",
                  uiOutput("name"),
                  uiOutput("id"),
                  uiOutput("filt"),
                  numericInput("thr", "allowed max.miss.num at least one condition", min = 0, max = 20, value = 0),
                  uiOutput("control"),
                  uiOutput("test_manual")
                  ),
                 shinyBS::bsTooltip("name", "Name of the column containing feature names", "right", options = list(container = "body")),
                 shinyBS::bsTooltip("id", "Name of the column containing feature IDs", "right", options = list(container = "body")),
                 shinyBS::bsTooltip("filt", "The filtered columns based on", "right", options = list(container = "body")),
                 shinyBS::bsTooltip("thr", "Set the threshold for the allowed number of missing values in at least one condition", "right", options = list(container = "body")),
                 shinyBS::bsTooltip("control", "Choose your control condition", "right", options = list(container = "body")),
                 shinyBS::bsTooltip("test_manual", "Choose the contrasts that you want to test", "right", options = list(container = "body")),
                 menuItem("Imputation options",
                 radioButtons("imputation",
                              "Imputation type",
                              choices = c(c("man",c("bpca", "knn", "QRILC", "MLE", "MinDet",
                                                    "MinProb", "min", "zero", "mixed",
                                                    "nbavg", "with", "none"))[1:9], "mixed on proteins", "mixed on samples"),
                              selected = "MinProb"),
                 shinyBS::bsTooltip("imputation", "Choose an imputation method", "right", options = list(container = "body")),
                 actionButton("help_imputation",
                    label = "Detailed information", #icon = icon("question-circle"),
                    style = "color: #f6f6f6; background-color: #2c3b41; border-color: #2c3b41; text-align: left; margin: 0 0 0px 0px",
                    #rgba(255, 255, 255, 0.1) rgba(0, 0, 0, 0.1)
                  ),
                 shinyBS::bsTooltip(
                    "help_imputation",
                    "The detailed information of imputation methods",
                    "right",
                    options = list(container = "body")
                  )
                 # p(a("Detailed information link",
                 #     href = "https://www.rdocumentation.org/packages/MSnbase/versions/1.20.7/topics/impute-methods",
                 #     target="_blank"))
        ),
                 menuItem("Type of FDR correction",
                  radioButtons("fdr_correction",
                              "Type of FDR correction",
                              choices = c("BH", "fdrtool"),
                              selected = "fdrtool"),
                  shinyBS::bsTooltip("fdr_correction", "Choose the method of pvalue adjustment", "right", options = list(container = "body"))
                  ),
                 actionButton("analyze", "Analyze"),
                 shinyBS::bsTooltip("analyze", "Click on it to analyze your data", "right", options = list(container = "body")),
                 tags$hr(),
                 tags$style(type="text/css", "#downloadData {background-color:white;color: black;font-family: Source Sans Pro}"),
                 uiOutput("downloadTable"),
                 uiOutput("downloadButton"),
                 h6(),
                 tags$style(type="text/css", "#Save_RData {background-color:white;color: black;font-family: Source Sans Pro}"),
                 uiOutput("downloadButton_for_save_RData"),
                 shinyBS::bsTooltip("downloadTable", "Choose a dataset to save, and here we offer four forms of datasets for downloading including results.txt, significant_proteins.txt, displayed_subset.txt, full_dataset.txt", "right", options = list(container = "body")),
                 shinyBS::bsTooltip("downloadButton", "Click on it to download the selected dataset", "right", options = list(container = "body")),
                 shinyBS::bsTooltip("downloadButton_for_save_RData", "Click on it to Save the RData, and when you upload the results.RData in the menuItem [Files] next time, you can get the same result as this time", "right", options = list(container = "body"))
        ), "DEP_tab"
          ),
        convertMenuItem(#** Annotation menuItem ----
        menuItem("Annotation options", selected = TRUE, tabName = "Annotation_tab", icon = icon("th"),
          menuItem("Import from", selected = TRUE,
                      checkboxInput("import_from_for_annotation", "DEP-LFQ or DEG-RNAseq", value = FALSE)
            ),
          uiOutput("import_for_annotation"),
          uiOutput("import_contrast_for_annotation"),
          # fluidRow(column(width = 5, checkboxInput("From_LFQ", "DEP-LFQ", value = FALSE)),
          #          column(width = 6, checkboxInput("From_RNAseq", "DEG-RNAseq", value = FALSE))),
          uiOutput("text_input_for_annotation"),
          uiOutput("organism_for_annotation"),
          # textAreaInput(inputId = "text_input_for_annotation", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%"),
          # selectizeInput("organism_for_annotation",
          #            "Select organism",
          #            choices=c("human","mouse"),
          #            selected="human"),
          actionButton("analyze_for_annotation", "Analyze"),
          tags$hr(),
          tags$style(type="text/css", "#downloadannotation {background-color:white;color: black;font-family: Source Sans Pro}"),
          uiOutput("downloadannotation_for_output"),
          shinyBS::bsTooltip("import_from_for_annotation", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),
          shinyBS::bsTooltip("import_for_annotation", "Choose genes that you want to do the gene annotation. UPregu for DEP-LFQ: up-regulated genes for DEP-LFQ panel; DOWNregu for DEP-LFQ: down-regulated genes for DEP-LFQ panel; UPDOWN for DEP-LFQ: all regulated genes, up- and down- regulated for DEP-LFQ panel; UPregu for DEG-RNAseq: up-regulated genes for DEG-RNAseq panel; DOWNregu for DEG-RNAseq: down-regulated genes for DEG-RNAseq panel; UPDOWN for DEG-RNAseq: all regulated genes, up- and down- regulated for DEG-RNAseq panel", "right", options = list(container = "body")),
          shinyBS::bsTooltip("import_contrast_for_annotation", "The contrast that you want to do the gene annotation. Note that: [Any significant] represents the genes that are significant in at least one contrast", "right", options = list(container = "body")),
          shinyBS::bsTooltip("text_input_for_annotation", "Paste gene list here", "right", options = list(container = "body")),
          shinyBS::bsTooltip("organism_for_annotation", "Select the organism", "right", options = list(container = "body")),
          shinyBS::bsTooltip("analyze_for_annotation", "Click on it to analyze", "right", options = list(container = "body")),
          shinyBS::bsTooltip("downloadannotation_for_output", "Click on it to download the annotation table", "right", options = list(container = "body"))
          ), "Annotation_tab"
          ),
        menuItem("ORA options", selected = TRUE, icon = icon("th"),
          convertMenuItem(#** GO analysis ----
            menuItem("GO options", selected = TRUE, tabName = "GO_tab", icon = icon("angle-double-right"),#bookmark bars book-open
              menuItem("Import from", selected = TRUE,
                          checkboxInput("import_from_for_go", "DEP-LFQ or DEG-RNAseq", value = FALSE)
                ),
              uiOutput("import_for_go"),
              uiOutput("import_contrast_for_go"),
              uiOutput("text_input_for_go"),
              uiOutput("organism_for_go"),
              uiOutput("df_with_lg2fc"),
              # textAreaInput(inputId = "text_input_for_go", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%"),
              # selectizeInput("organism_for_go",
              #            "Select organism",
              #            choices=c("human","mouse"),
              #            selected="human"),
              # checkboxInput("df_with_lg2fc", "If with log2 fold change", value = FALSE),
              # radioButtons("go_color",
              #                "colorBy",
              #                c("pvalue", "p.adjust"),
              #                selected = "p.adjust"),
              actionButton("analyze_for_go", "Analyze"),
              tags$hr(),
              tags$style(type="text/css", "#downloadgo {background-color:white;color: black;font-family: Source Sans Pro}"),
              uiOutput("downloadTable_go"),
              uiOutput("downloadButton_go"),
              # downloadButton("downloadgo", "Save table")
              shinyBS::bsTooltip("import_from_for_go", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),
              shinyBS::bsTooltip("import_for_go", "Choose genes that you want to do the GO analysis. UPregu for DEP-LFQ: up-regulated genes for DEP-LFQ panel; DOWNregu for DEP-LFQ: down-regulated genes for DEP-LFQ panel; UPDOWN for DEP-LFQ: all regulated genes, up- and down- regulated for DEP-LFQ panel; UPregu for DEG-RNAseq: up-regulated genes for DEG-RNAseq panel; DOWNregu for DEG-RNAseq: down-regulated genes for DEG-RNAseq panel; UPDOWN for DEG-RNAseq: all regulated genes, up- and down- regulated for DEG-RNAseq panel", "right", options = list(container = "body")),
              shinyBS::bsTooltip("import_contrast_for_go", "The contrast that you want to do the GO analysis. Note that: [Any significant] represents the genes that are significant in at least one contrast, and when you select [Any significant], you should not choose [If with log2 fold change]", "right", options = list(container = "body")),
              shinyBS::bsTooltip("text_input_for_go", "Paste your gene list here", "right", options = list(container = "body")),
              shinyBS::bsTooltip("organism_for_go", "Select the organism", "right", options = list(container = "body")),
              shinyBS::bsTooltip("df_with_lg2fc", "Whether your gene list with log2 fold change", "right", options = list(container = "body")),
              shinyBS::bsTooltip("analyze_for_go", "Click on it to analyze", "right", options = list(container = "body")),
              shinyBS::bsTooltip("downloadTable_go", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
              shinyBS::bsTooltip("downloadButton_go", "Click on it to download the result table", "right", options = list(container = "body"))
              ), "GO_tab"
            ),
          convertMenuItem(#** KEGG analysis ----
            menuItem("KEGG options", selected = TRUE, tabName = "KEGG_tab", icon = icon("angle-double-right"),
              # textAreaInput(inputId = "text_input_for_kegg", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%"),
              # selectizeInput("organism_for_kegg",
                         # "Select organism",
                         # choices=c("hsa","mmu"),
                         # selected="hsa"),
              # checkboxInput("df_with_lg2fc_for_kegg", "If with log2 fold change", value = FALSE),
              # radioButtons("kegg_color",
              #                "colorBy",
              #                c("pvalue", "p.adjust"),
              #                selected = "p.adjust"),
              menuItem("Import from", selected = TRUE,
                          checkboxInput("import_from_for_kegg", "DEP-LFQ or DEG-RNAseq", value = FALSE)
                ),
              uiOutput("import_for_kegg"),
              uiOutput("import_contrast_for_kegg"),
              uiOutput("text_input_for_kegg"),
              uiOutput("organism_for_kegg"),
              uiOutput("df_with_lg2fc_for_kegg"),
              actionButton("analyze_for_kegg", "Analyze"),
              tags$hr(),
              tags$style(type="text/css", "#downloadkegg {background-color:white;color: black;font-family: Source Sans Pro}"),
              uiOutput("downloadTable_kegg"),
              uiOutput("downloadButton_kegg"),
              # downloadButton("downloadkegg", "Save table")
              shinyBS::bsTooltip("import_from_for_kegg", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),
              shinyBS::bsTooltip("import_for_kegg", "Choose genes that you want to do the KEGG analysis. UPregu for DEP-LFQ: up-regulated genes for DEP-LFQ panel; DOWNregu for DEP-LFQ: down-regulated genes for DEP-LFQ panel; UPDOWN for DEP-LFQ: all regulated genes, up- and down- regulated for DEP-LFQ panel; UPregu for DEG-RNAseq: up-regulated genes for DEG-RNAseq panel; DOWNregu for DEG-RNAseq: down-regulated genes for DEG-RNAseq panel; UPDOWN for DEG-RNAseq: all regulated genes, up- and down- regulated for DEG-RNAseq panel", "right", options = list(container = "body")),
              shinyBS::bsTooltip("import_contrast_for_kegg", "The contrast that you want to do the KEGG analysis. Note that: [Any significant] represents the genes that are significant in at least one contrast, and when you select [Any significant], you should not choose [If with log2 fold change]", "right", options = list(container = "body")),
              shinyBS::bsTooltip("text_input_for_kegg", "Paste your gene list here", "right", options = list(container = "body")),
              shinyBS::bsTooltip("organism_for_kegg", "Select the organism: hsa represents human, mmu represents mouse", "right", options = list(container = "body")),
              shinyBS::bsTooltip("df_with_lg2fc_for_kegg", "Whether your gene list with log2 fold change", "right", options = list(container = "body")),
              shinyBS::bsTooltip("analyze_for_kegg", "Click on it to analyze", "right", options = list(container = "body")),
              shinyBS::bsTooltip("downloadTable_kegg", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
              shinyBS::bsTooltip("downloadButton_kegg", "Click on it to download the result table", "right", options = list(container = "body"))
              ), "KEGG_tab"
            ),
          convertMenuItem(#** Reactome analysis ----
            menuItem("Reactome options", selected = TRUE, tabName = "Reactome_tab", icon = icon("angle-double-right"),
              # textAreaInput(inputId = "text_input_for_reactome", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%"),
              # selectizeInput("organism_for_reactome",
              #            "Select organism",
              #            choices=c("human","mouse"),
              #            selected="human"),
              # checkboxInput("df_with_lg2fc_for_reactome", "If with log2 fold change", value = FALSE),
              # radioButtons("reactome_color",
              #                "colorBy",
              #                c("pvalue", "p.adjust"),
              #                selected = "p.adjust"),
              menuItem("Import from", selected = TRUE,
                          checkboxInput("import_from_for_reactome", "DEP-LFQ or DEG-RNAseq", value = FALSE)
                ),
              uiOutput("import_for_reactome"),
              uiOutput("import_contrast_for_reactome"),
              uiOutput("text_input_for_reactome"),
              uiOutput("organism_for_reactome"),
              uiOutput("df_with_lg2fc_for_reactome"),
              actionButton("analyze_for_reactome", "Analyze"),
              tags$hr(),
              tags$style(type="text/css", "#downloadreactome {background-color:white;color: black;font-family: Source Sans Pro}"),
              uiOutput("downloadTable_reactome"),
              uiOutput("downloadButton_reactome"),
              # downloadButton("downloadreactome", "Save table")
              shinyBS::bsTooltip("import_from_for_reactome", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),
              shinyBS::bsTooltip("import_for_reactome", "Choose genes that you want to do the reactome analysis. UPregu for DEP-LFQ: up-regulated genes for DEP-LFQ panel; DOWNregu for DEP-LFQ: down-regulated genes for DEP-LFQ panel; UPDOWN for DEP-LFQ: all regulated genes, up- and down- regulated for DEP-LFQ panel; UPregu for DEG-RNAseq: up-regulated genes for DEG-RNAseq panel; DOWNregu for DEG-RNAseq: down-regulated genes for DEG-RNAseq panel; UPDOWN for DEG-RNAseq: all regulated genes, up- and down- regulated for DEG-RNAseq panel", "right", options = list(container = "body")),
              shinyBS::bsTooltip("import_contrast_for_reactome", "The contrast that you want to do the reactome analysis. Note that: [Any significant] represents the genes that are significant in at least one contrast, and when you select [Any significant], you should not choose [If with log2 fold change]", "right", options = list(container = "body")),
              shinyBS::bsTooltip("text_input_for_reactome", "Paste your gene list here", "right", options = list(container = "body")),
              shinyBS::bsTooltip("organism_for_reactome", "Select the organism", "right", options = list(container = "body")),
              shinyBS::bsTooltip("df_with_lg2fc_for_reactome", "Whether your gene list with log2 fold change", "right", options = list(container = "body")),
              shinyBS::bsTooltip("analyze_for_reactome", "Click on it to analyze", "right", options = list(container = "body")),
              shinyBS::bsTooltip("downloadTable_reactome", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
              shinyBS::bsTooltip("downloadButton_reactome", "Click on it to download the result table", "right", options = list(container = "body"))
              ), "Reactome_tab"
            )
          ),
        menuItem("GSEA options", selected = TRUE, icon = icon("th"),
          convertMenuItem(#** gseGO analysis----
            menuItem("gseGO options", selected = TRUE, tabName = "gseGO_tab", icon = icon("angle-double-right"),
            # textAreaInput(inputId = "text_input_for_gsego", label = "Please paste your order ranked geneList", placeholder = "TP53 2.1\nPTEN 1.2", rows = 8, width = "100%"),
            # selectizeInput("organism_for_gsego",
                       # "Select organism",
                       # choices=c("human","mouse"),
                       # selected="human"),
            # radioButtons("gsego_color",
            #                "colorBy",
            #                c("pvalue", "p.adjust"),
            #                selected = "p.adjust"),
            menuItem("Import from", selected = TRUE,
                        checkboxInput("import_from_for_gsego", "DEP-LFQ or DEG-RNAseq", value = FALSE)
              ),
            uiOutput("import_for_gsego"),
            uiOutput("import_contrast_for_gsego"),
            uiOutput("text_input_for_gsego"),
            uiOutput("organism_for_gsego"),
            actionButton("analyze_for_gsego", "Analyze"),
            tags$hr(),
            tags$style(type="text/css", "#downloadgsego {background-color:white;color: black;font-family: Source Sans Pro}"),
            uiOutput("downloadTable_gsego"),
            uiOutput("downloadButton_gsego"),
            # downloadButton("downloadgsego", "Save table")
            shinyBS::bsTooltip("import_from_for_gsego", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),
            shinyBS::bsTooltip("import_for_gsego", "Choose genes that you want to do the gseGO analysis. DEP-LFQ: import from DEP-LFQ panel; DEG-RNAseq: import from DEG-RNAseq panel", "right", options = list(container = "body")),
            shinyBS::bsTooltip("import_contrast_for_gsego", "The contrast that you want to do the gseGO analysis", "right", options = list(container = "body")),
            shinyBS::bsTooltip("text_input_for_gsego", "Paste your gene list with log2 fold change here", "right", options = list(container = "body")),
            shinyBS::bsTooltip("organism_for_gsego", "Select the organism", "right", options = list(container = "body")),
            shinyBS::bsTooltip("analyze_for_gsego", "Click on it to analyze", "right", options = list(container = "body")),
            shinyBS::bsTooltip("downloadTable_gsego", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
            shinyBS::bsTooltip("downloadButton_gsego", "Click on it to download the result table", "right", options = list(container = "body"))
            ), "gseGO_tab"
            ),
          convertMenuItem(#** gseKEGG analysis----
            menuItem("gseKEGG options", selected = TRUE, tabName = "gseKEGG_tab", icon = icon("angle-double-right"),
              # textAreaInput(inputId = "text_input_for_gsekegg", label = "Please paste your order ranked geneList", placeholder = "TP53 2.1\nPTEN 1.2", rows = 8, width = "100%"),
              # selectizeInput("organism_for_gsekegg",
              #            "Select organism",
              #            choices=c("hsa","mmu"),
              #            selected="hsa"),
              # radioButtons("gsekegg_color",
              #                "colorBy",
              #                c("pvalue", "p.adjust"),
              #                selected = "p.adjust"),
              menuItem("Import from", selected = TRUE,
                          checkboxInput("import_from_for_gsekegg", "DEP-LFQ or DEG-RNAseq", value = FALSE)
                ),
              uiOutput("import_for_gsekegg"),
              uiOutput("import_contrast_for_gsekegg"),
              uiOutput("text_input_for_gsekegg"),
              uiOutput("organism_for_gsekegg"),
              actionButton("analyze_for_gsekegg", "Analyze"),
              tags$hr(),
              tags$style(type="text/css", "#downloadgsekegg {background-color:white;color: black;font-family: Source Sans Pro}"),
              uiOutput("downloadTable_gsekegg"),
              uiOutput("downloadButton_gsekegg"),
              # downloadButton("downloadgsekegg", "Save table")
              shinyBS::bsTooltip("import_from_for_gsekegg", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),
              shinyBS::bsTooltip("import_for_gsekegg", "Choose genes that you want to do the gseKEGG analysis. DEP-LFQ: import from DEP-LFQ panel; DEG-RNAseq: import from DEG-RNAseq panel", "right", options = list(container = "body")),
              shinyBS::bsTooltip("import_contrast_for_gsekegg", "The contrast that you want to do the gseKEGG analysis", "right", options = list(container = "body")),
              shinyBS::bsTooltip("text_input_for_gsekegg", "Paste your gene list with log2 fold change here", "right", options = list(container = "body")),
              shinyBS::bsTooltip("organism_for_gsekegg", "Select the organism: hsa represents human, mmu represents mouse", "right", options = list(container = "body")),
              shinyBS::bsTooltip("analyze_for_gsekegg", "Click on it to analyze", "right", options = list(container = "body")),
              shinyBS::bsTooltip("downloadTable_gsekegg", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
              shinyBS::bsTooltip("downloadButton_gsekegg", "Click on it to download the result table", "right", options = list(container = "body"))
              ), "gseKEGG_tab"
            ),
          convertMenuItem(#** gseReactome analysis----
            menuItem("gseReactome options", selected = TRUE, tabName = "gseReactome_tab", icon = icon("angle-double-right"),
              # textAreaInput(inputId = "text_input_for_gsereactome", label = "Please paste your order ranked geneList", placeholder = "TP53 2.1\nPTEN 1.2", rows = 8, width = "100%"),
              # selectizeInput("organism_for_gsereactome",
              #            "Select organism",
              #            choices=c("human","mouse"),
              #            selected="human"),
              # radioButtons("gsereactome_color",
              #                "colorBy",
              #                c("pvalue", "p.adjust"),
              #                selected = "p.adjust"),
              menuItem("Import from", selected = TRUE,
                          checkboxInput("import_from_for_gsereactome", "DEP-LFQ or DEG-RNAseq", value = FALSE)
                ),
              uiOutput("import_for_gsereactome"),
              uiOutput("import_contrast_for_gsereactome"),
              uiOutput("text_input_for_gsereactome"),
              uiOutput("organism_for_gsereactome"),
              actionButton("analyze_for_gsereactome", "Analyze"),
              tags$hr(),
              tags$style(type="text/css", "#downloadgsereactome {background-color:white;color: black;font-family: Source Sans Pro}"),
              uiOutput("downloadTable_gsereactome"),
              uiOutput("downloadButton_gsereactome"),
              # downloadButton("downloadgsereactome", "Save table")
              shinyBS::bsTooltip("import_from_for_gsereactome", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),
              shinyBS::bsTooltip("import_for_gsereactome", "Choose genes that you want to do the gseReactome analysis. DEP-LFQ: import from DEP-LFQ panel; DEG-RNAseq: import from DEG-RNAseq panel", "right", options = list(container = "body")),
              shinyBS::bsTooltip("import_contrast_for_gsereactome", "The contrast that you want to do the gseReactome analysis", "right", options = list(container = "body")),
              shinyBS::bsTooltip("text_input_for_gsereactome", "Paste your gene list with log2 fold change here", "right", options = list(container = "body")),
              shinyBS::bsTooltip("organism_for_gsereactome", "Select the organism", "right", options = list(container = "body")),
              shinyBS::bsTooltip("analyze_for_gsereactome", "Click on it to analyze", "right", options = list(container = "body")),
              shinyBS::bsTooltip("downloadTable_gsereactome", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
              shinyBS::bsTooltip("downloadButton_gsereactome", "Click on it to download the result table", "right", options = list(container = "body"))
              ), "gseReactome_tab"
            )
          ),
        convertMenuItem(#** PPI menuItem----
        menuItem("PPITools options", selected = TRUE, tabName = "PPI_tab", icon = icon("th"),
                 # selectInput("organism",
                 #             "select organism",
                 #             choices=dir(system.file("PPIdata",package = "DEP2")),
                 #             selected="human"),
                 # # actionButton("load_PPIdata","Load Datas First!"),
                 # helpText("choose organism string data to load", br(),"If it doesn't exist app will report an error.", br(),"you need to download file from string following it"),
                 # h4("paste your gene list"),
          menuItem("Import from", selected = TRUE,
                      checkboxInput("import_from_for_ppi", "DEP-LFQ or DEG-RNAseq", value = FALSE)
            ),
          uiOutput("import_for_ppi"),
          uiOutput("import_contrast_for_ppi"),
                 # tags$style(type="text/css", "textarea {width:80%}"),
                 # textAreaInput(inputId = "input_text", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%"),
                 uiOutput("input_text"),
                 selectInput("organism",
                             "Select organism",
                             choices=dir(system.file("PPIdata",package = "DEP2")),
                             selected="human"),
                 # actionButton("load_PPIdata","Load Datas First!"),
                 helpText("choose organism string data to load", br(),"If it doesn't exist, app will report an error.", br(),"you need to download file from string first"),
                 selectizeInput("chooseScore",
                                "Select concerned scores" ,
                                choices=c("combined_score" ,"neighborhood","fusion","cooccurence","coexpression","experimental","database","textmining"),
                                multiple = TRUE,
                                selected = c("combined_score" ,"neighborhood","fusion","cooccurence","coexpression","experimental","database","textmining")),
                 fluidRow(
                   column(width = 10,
                          numericInput("scorecutoff",
                                       "scores cutoff",
                                       min = 400, max = 1000, value = 400))),
                 actionButton("String_annalysis","String It!"),
                 # radioButtons("reactome_color",
                 #                "colorBy",
                 #                c("pvalue", "p.adjust"),
                 #                selected = "p.adjust"),
                 # h4("\n"),
                 tags$hr(),
                 tags$style(type="text/css", "#downloadPPITable {background-color:white;color: black;font-family: Source Sans Pro}"),
                 tags$style(type="text/css", "#downloadPPInetwork {background-color:white;color: black;font-family: Source Sans Pro}"),
                 # h4("\n"),
                 uiOutput("download_PPITable"),
                 h6(),
                 uiOutput("download_PPInetwork"),
                 # downloadButton("downloadPPITable", "Save table"),
                 # downloadButton("downloadPPInetwork", "Save network")
          shinyBS::bsTooltip("import_from_for_ppi", "Import genes from the result of DEP-LFQ or DEG-RNAseq panel", "top", options = list(container = "body")),
          shinyBS::bsTooltip("import_for_ppi", "Choose genes that you want to do the PPI(Protein-Protein Interaction Networks). UPregu for DEP-LFQ: up-regulated genes for DEP-LFQ panel; DOWNregu for DEP-LFQ: down-regulated genes for DEP-LFQ panel; UPDOWN for DEP-LFQ: all regulated genes, up- and down- regulated for DEP-LFQ panel; UPregu for DEG-RNAseq: up-regulated genes for DEG-RNAseq panel; DOWNregu for DEG-RNAseq: down-regulated genes for DEG-RNAseq panel; UPDOWN for DEG-RNAseq: all regulated genes, up- and down- regulated for DEG-RNAseq panel", "right", options = list(container = "body")),
          shinyBS::bsTooltip("import_contrast_for_ppi", "The contrast that you want to do the PPI(Protein-Protein Interaction Networks). Note that: [Any significant] represents the genes that are significant in at least one contrast", "right", options = list(container = "body")),
          shinyBS::bsTooltip("input_text", "Paste gene list here", "top", options = list(container = "body")),
          shinyBS::bsTooltip("organism", "Select the organism", "right", options = list(container = "body")),
          shinyBS::bsTooltip("chooseScore", "Select which type of evidence will contribute to the prediction of the score under the active interaction sources", "right", options = list(container = "body")),
          shinyBS::bsTooltip("scorecutoff", "Minimum required interaction score. The minimum required interaction score puts a threshold on the confidence score, such that only interaction above this score are included in the predicted network. Confidence limits are as follows: low confidence - 0.15 (or better); medium confidence - 0.4; high confidence - 0.7; highest confidence - 0.9", "right", options = list(container = "body")),
          shinyBS::bsTooltip("String_annalysis", "Click on it to analyze", "right", options = list(container = "body")),
          shinyBS::bsTooltip("download_PPITable", "Click on it to download the Protein-Protein Interaction table", "right", options = list(container = "body")),
          shinyBS::bsTooltip("download_PPInetwork", "Click on it to download the PPI(Protein-Protein Interaction Network)", "right", options = list(container = "body"))
        ), "PPI_tab"
          ),
        convertMenuItem(#** DEG-RNAseq menuItem----
        menuItem("DEG-RNAseq options", selected = TRUE, tabName = "RNAseq_tab", icon = icon("th"),
                 menuItem("Files", selected = TRUE,
                   fluidRow(
                    column(width = 9,
                    fileInput('file1_for_DERNAseq',width = "300px",
                           'Countmatrix.txt',
                           accept=c('text/csv',
                                    'text/comma-separated-values,text/plain',
                                    '.csv')),
                    fileInput('file2_for_DERNAseq',width = "300px",
                           'ExperimentalDesign.txt',
                           accept=c('text/csv',
                                    'text/comma-separated-values,text/plain',
                                    '.csv'))),
                    column(width = 1,h4(),
                           br(),
                     actionButton("help_format_RNA",
                              label = "", icon = icon("question-circle"),
                              style = "color: #f6f6f6; background-color: #2c3b41; border-color: #2c3b41"
                 ),
                 shinyBS::bsTooltip(
                   "help_format_RNA",
                   "How to provide your input data",
                   "top",
                   options = list(container = "body")
                 ))
        )
        ),
                 menuItem("RNAseq_settings", selected = TRUE,
                  uiOutput("dds_design"),
                  uiOutput("choose_expfac"),
                  uiOutput("whether_LRT"),
                  uiOutput("reduced_model"),
                  # uiOutput("choose_expfac1"),
                  # uiOutput("choose_expfac2"),
                  uiOutput("control_for_RNAseq"),
                  uiOutput("test_manual_for_RNAseq"),
                  numericInput("filter_rowsum", "Threshold row sums", min = 0, max = 10000, value = 10)
                  ),
                 menuItem("ID Transformation", selected = TRUE,
                  checkboxInput("transid_for_RNAseq", "ID transformation", value = FALSE),
                  uiOutput("ui_set_no_map_to_rowname"),
                  uiOutput("ui_selectspecies"),
                  uiOutput("ui_idtype")
                  ),
                 menuItem("result options", selected = TRUE,
                  checkboxInput("independent_filtering", "Independent filtering", value = TRUE),
                  checkboxInput("Shrink_lfc", "Shrink lfc", value = FALSE),
                  checkboxInput("Use_IHW", "Use IHW", value = FALSE)
                  ),
                 uiOutput("ui_nrcores"),
                 actionButton("analyze_for_DERNAseq", "Analyze"),
                 tags$hr(),
                 tags$style(type="text/css", "#download_DERNAseq {background-color:white;color: black;font-family: Source Sans Pro}"),
                 uiOutput("downloadTable_for_DERNAseq"),
                 uiOutput("downloadButton_for_DERNAseq"),
          shinyBS::bsTooltip("dds_design", "Select your design, and it can be multi-factor", "right", options = list(container = "body")),
          shinyBS::bsTooltip("choose_expfac", "Choose the contrast to build upon", "right", options = list(container = "body")),
          shinyBS::bsTooltip("whether_LRT", "I can perform a LRT test on the chosen factor, choose whether perform a LRT test", "right", options = list(container = "body")),
          shinyBS::bsTooltip("reduced_model", "Select the reduced model. Note: 1) the full model is the same as your selected design; 2) if the number of your selected design is 1, select 1", "right", options = list(container = "body")),
          # shinyBS::bsTooltip("choose_expfac1", "Select the name of the numerator level for the fold change", "right", options = list(container = "body")),
          # shinyBS::bsTooltip("choose_expfac2", "Select the name of the denominator level for the fold change (must be different from the numerator)", "right", options = list(container = "body")),
          shinyBS::bsTooltip("control_for_RNAseq", "Choose your control condition", "right", options = list(container = "body")),
          shinyBS::bsTooltip("test_manual_for_RNAseq", "Choose the contrasts that you want to test", "right", options = list(container = "body")),
          shinyBS::bsTooltip("filter_rowsum", "Threshold on the row sums of the counts", "right", options = list(container = "body")),
          shinyBS::bsTooltip("transid_for_RNAseq", "Whether transform your id to gene symbol", "right", options = list(container = "body")),
          shinyBS::bsTooltip("ui_set_no_map_to_rowname", "Whether the ids that can not map to a gene symbol keep the orignal ids", "right", options = list(container = "body")),
          shinyBS::bsTooltip("ui_selectspecies", "Select the species", "right", options = list(container = "body")),
          shinyBS::bsTooltip("ui_idtype", "Select your id type", "right", options = list(container = "body")),
          shinyBS::bsTooltip("independent_filtering", "whether independent filtering should be applied automatically", "right", options = list(container = "body")),
          shinyBS::bsTooltip("Shrink_lfc", "whether shrink the log fold change for the contrast of interest", "right", options = list(container = "body")),
          shinyBS::bsTooltip("Use_IHW", "Whether use Independent Hypothesis Weighting (IHW) as a filtering function", "right", options = list(container = "body")),
          shinyBS::bsTooltip("ui_nrcores", "Choose how many cores to use for computing", "right", options = list(container = "body")),
          shinyBS::bsTooltip("analyze_for_DERNAseq", "Click on it to analyze", "right", options = list(container = "body")),
          shinyBS::bsTooltip("downloadTable_for_DERNAseq", "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
          shinyBS::bsTooltip("downloadButton_for_DERNAseq", "Click on it to download the result table", "right", options = list(container = "body"))
        ), "RNAseq_tab"
          ),
        convertMenuItem(#** PR-Heatmap menuItem----
        menuItem("PR-Heatmap options", selected = TRUE, tabName = "PRHeatmap_tab", icon = icon("th"),
                menuItem("Import from", selected = TRUE,
                  checkboxInput("import_from_for_PR_Heatmap", "DEP-LFQ and DEG-RNAseq", value = FALSE),
                  uiOutput("import_contrast_for_PR_Heatmap")),
                # menuItem(
                   fluidRow(
                  column(width = 9,
                   uiOutput("file1_for_pro_rna"),
                   uiOutput("file2_for_pro_rna"),
                   uiOutput("file3_for_pro_rna")),
                column(width = 1,h4(),
                       br(),
                       uiOutput("help_format_PRO_RNA"),
                       shinyBS::bsTooltip(
                         "help_format_PRO_RNA",
                         "How to provide your input data",
                         "top",
                         options = list(container = "body")
                       ))
                   ),
                  # ),
        #          menuItem("Files", selected = TRUE,
        #            fluidRow(
        #           column(width = 9,
        #            # fileInput('file1_for_pro_rna',width = "300px",
        #            #         'Proteinmatrix.txt',
        #            #         accept=c('text/csv',
        #            #                  'text/comma-separated-values,text/plain',
        #            #                  '.csv')),
        #          fileInput('file2_for_pro_rna',width = "300px",
        #                    'RNAseq_normalized_count\nmatrix.txt',
        #                    accept=c('text/csv',
        #                             'text/comma-separated-values,text/plain',
        #                             '.csv')),
        #          fileInput('file3_for_pro_rna',width = "300px",
        #                    'Genelist.txt',
        #                    accept=c('text/csv',
        #                             'text/comma-separated-values,text/plain',
        #                             '.csv'))),
        #         column(width = 1,h4(),
        #                br(),
        #                actionButton("help_format_PRO_RNA",
        #                             label = "", icon = icon("question-circle"),
        #                             style = "color: #f6f6f6; background-color: #2c3b41; border-color: #2c3b41"
        #                ),
        #                shinyBS::bsTooltip(
        #                  "help_format_PRO_RNA",
        #                  "How to provide your input data",
        #                  "top",
        #                  options = list(container = "body")
        #                ))
        # )
        # ),
                 actionButton("analyze_for_pro_rna", "Start"),
                 tags$hr(),
                 tags$style(type="text/css", "#download_table_for_pro_rna {background-color:white;color: black;font-family: Source Sans Pro}"),
                 uiOutput("downloadbutton_table_for_pro_rna"),
          shinyBS::bsTooltip("import_from_for_PR_Heatmap", "Import results from the result of DEP-LFQ and DEG-RNAseq panel", "top", options = list(container = "body")),
          shinyBS::bsTooltip("import_contrast_for_PR_Heatmap", "The contrast that you are interested in. Note that: 1). The contrast tested in both DEP-LFQ panel and DEG-RNAseq panel  must be the same; 2) It will extract the significant results of both DEP-LFQ and DEG-RNAseq upon your selected contrast; 3).[Any significant] represents the genes that are significant in at least one contrast", "right", options = list(container = "body")),
          shinyBS::bsTooltip("analyze_for_pro_rna", "Click on it to start", "right", options = list(container = "body")),
          shinyBS::bsTooltip("downloadbutton_table_for_pro_rna", "Click on it to download the table", "right", options = list(container = "body"))
        ), "PRHeatmap_tab"
          ),
        convertMenuItem(#** Extralinks menuItem ----
          menuItem("Extralinks", selected = FALSE, tabName = "Extralinks_tab", icon = icon("th")), "Extralinks_tab"
          )
      )
    ),

    dashboardBody(
  #     setBackgroundColor(
  #   color = "white",
  #   gradient = "linear",
  #   direction = "bottom",shinydashboard = TRUE
  # ),
  # tags$head(tags$style(HTML('/* body */.content-wrapper, .right-side {background-color: white;}'))),
  # tags$style(".content {background-color: white;}"),
      includeCSS(system.file("shiny_apps/LFQ/www", "style.css", package = "DEP2")),
      tabItems(#* tabItems----
      tabItem(#** Welcome tabItem----
        tabName = "Welcome_tab",
            fluidRow(
              column(style = "background-color:white;",
                width = 12,
                includeMarkdown(system.file("extdata", "Tutorial.md", package = "DEP2"))
              )
            )
        ),
        tabItem(#** DEP tabItem----
          tabName = "DEP_tab",
          helpText("Please cite: "),
      fluidRow(
        box(numericInput("p",
                         "adj. P value",
                         min = 0.0001, max = 0.1, value = 0.05),
            width = 2),
        box(numericInput("lfc",
                         "Log2 fold change",
                         min = 0, max = 10, value = 1),
            width = 2),
        infoBoxOutput("significantBox"),
        box(radioButtons("pres",
                         "Data presentation",
                         c("contrast", "centered"),
                         selected = "centered"),
            width = 2),
        box(radioButtons(inputId = "contrasts",
                         label = "Contrasts",
                         choices = c("control", "all", "manual"),
                         selected = "control", inline = TRUE),
            width = 2),
          shinyBS::bsTooltip("p", "Set the false discovery rate threshold", "top", options = list(container = "body")),
          shinyBS::bsTooltip("lfc", "Set the log2 fold change threshold", "top", options = list(container = "body")),
          shinyBS::bsTooltip("pres", "The type of data scaling used for heatmap plotting. Either the fold change (contrast) or the centered log2-intensity (centered)", "top", options = list(container = "body")),
          shinyBS::bsTooltip("contrasts", "The type of contrasts that will be tested. This can be all possible pairwise comparisons (all), limited to the comparisons versus the control (control), or manually defined contrasts (manual)", "top", options = list(container = "body"))
      ),
      fluidRow(
        column(width = 7,
               box(title = "Top Table",
                   box(uiOutput("select"), width = 6),
                   box(uiOutput("exclude"), width = 6),
                   DT::dataTableOutput("table"), width = 12),
          shinyBS::bsTooltip("select", "Select the rows that are significant of the selected contrasts", "top", options = list(container = "body")),
          shinyBS::bsTooltip("exclude", "Exclude the rows that are significant of the selected contrasts", "top", options = list(container = "body")),
          shinyBS::bsTooltip("table", "Table of the significant results", "top", options = list(container = "body"))
        ),
        column(width = 5,
               tabBox(title = "Result Plots", width = 12,
                      tabPanel(title = "Selected Protein",
                        fluidRow(
                          box(numericInput("selected_plot_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                          box(numericInput("selected_plot_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                          ),
                        fluidRow(
                          plotOutput("selected_plot"),
                               downloadButton('downloadPlot', 'Save plot')
                          ),
          shinyBS::bsTooltip("selected_plot_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("selected_plot_Height", "Height of the figure to export", "top", options = list(container = "body"))
                               ),
                      tabPanel(title = "Heatmap",
                              fluidRow(
                                column(width = 6,
                                    box(checkboxInput("manual_heatmap",
                                                   "Manual heatmap",
                                                   value = FALSE), width = 12)
                                    ),
                                column(width = 6,
                                    box(selectizeInput("colorbar",
                                       "colorbar",
                                       choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                                       selected = c("RdBu"), multiple = FALSE), width = 12)
                                    )
                               ),
                               fluidRow(
                                box(uiOutput("heatmap_cntrst"), width = 12)
                                ),
                               # fluidRow(
                               #   box(checkboxInput("manual_heatmap",
                               #                     "Manual heatmap",
                               #                     value = FALSE),
                               #   # checkboxInput("same_trend",
                               #   #                   "Same trend heatmap",
                               #   #                   value = FALSE),
                               #    uiOutput("heatmap_cntrst"), width = 12, collapsible = FALSE, collapsed = FALSE)
                               #   #box(checkboxInput("manual_heatmap",
                               #   #                  "Manual heatmap",
                               #   #                  value = FALSE)),
                               # ),
                               fluidRow(
                                box(numericInput("k",
                                                  "Kmeans",
                                                  min = 0, max = 15, value = 6),
                                     width = 4),
                                box(numericInput("limit",
                                                  "Color",
                                                  min = 0, max = 16, value = 4),
                                     width = 4),
                                box(numericInput("size",
                                                  "Heat size(4-30)",
                                                  min = 4, max = 30, value = 10),
                                     width = 4)

                                ),
                               fluidRow(
                                 box(
                                  # selectizeInput("colorbar",
                                  #      "colorbar",
                                  #      choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                                  #      selected = c("RdBu"), multiple = FALSE),
                                  checkboxInput("if_rowname_color",
                                                   "rowname color",
                                                   value = FALSE),

                                  numericInput("row_font_size",
                                                  "row font size",
                                                  min = 0, max = 15, value = 6),
                                 numericInput("col_font_size",
                                                  "col font size",
                                                  min = 0, max = 16, value = 10),
                                 numericInput("Width",
                                                  "heatmap width",
                                                  min = 1, max = 30, value = 7),
                                 numericInput("Height",
                                                  "heatmap height",
                                                  min = 1, max = 30, value = 10),
                                 checkboxInput("cluster_columns",
                                                   "Cluster columns",
                                                   value = TRUE),
                                 uiOutput("Custom_columns_order"),
                                 checkboxInput("if_mysplit",
                                                   "If mysplit",
                                                   value = FALSE),
                                  uiOutput("mysplit"),
                                     width = 12, collapsible = TRUE, collapsed = TRUE),
                               ),
                               fluidRow(
                                 uiOutput("plot"),
                                 downloadButton('downloadHeatmap', 'Save heatmap')),
          shinyBS::bsTooltip("colorbar", "Choose the color palette of the colorbar", "top", options = list(container = "body")),
          shinyBS::bsTooltip("manual_heatmap", "Whether plot the manual heatmap", "top", options = list(container = "body")),
          shinyBS::bsTooltip("heatmap_cntrst", "Choose the contrasts that you want to show, act when [Manual heatmap] is TRUE, can be one or more contrasts", "top", options = list(container = "body")),
          shinyBS::bsTooltip("k", "Set the number of k-means clusters", "top", options = list(container = "body")),
          shinyBS::bsTooltip("limit", "Set the outer limits of the color scale", "top", options = list(container = "body")),
          shinyBS::bsTooltip("size", "Set the height of the plot region", "top", options = list(container = "body")),
          shinyBS::bsTooltip("if_rowname_color", "Whether set the color of the row names,(purple: Peptides == 1, blue: Peptides == 2)", "top", options = list(container = "body")),
          shinyBS::bsTooltip("row_font_size", "Set the size of row labels", "top", options = list(container = "body")),
          shinyBS::bsTooltip("col_font_size", "Set the size of column labels", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Height", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("cluster_columns", "Whether make cluster on columns", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Custom_columns_order", "Order of column, act when [Cluster columns] is FALSE", "top", options = list(container = "body")),
          shinyBS::bsTooltip("if_mysplit", "Whether order the splices manually", "top", options = list(container = "body")),
          shinyBS::bsTooltip("mysplit", "The ordered splices of the heatmap (from top to bottom), only act when both [If mysplit] and [Kmeans] are TRUE", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Volcano plot",
                               fluidRow(
                                 box(uiOutput("volcano_cntrst"), width = 9),
                                 box(numericInput("fontsize",
                                                  "Font size",
                                                  min = 0, max = 8, value = 4),

                                     width = 3)
                               ),
                               fluidRow(
                                 box(checkboxInput("my_breaks",
                                                   "My breaks",
                                                   value = FALSE),
                                 multiInput("mybreaks","Mybreaks",choices = seq(-40, 40, 0.5)),
                                 numericInput("Volcano_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7),
                                 numericInput("Volcano_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7),
                                 #selectizeInput("mybreaks","Mybreaks",choices = seq(-40, 40, 1), multiple = TRUE, size = 10),
                                     width = 8, collapsible = TRUE, collapsed = TRUE),
                                 box(checkboxInput("check_names",
                                                   "Display names",
                                                   value = FALSE),
                                     checkboxInput("p_adj",
                                                   "Adjusted p values",
                                                   value = FALSE),
                                     checkboxInput("same_width",
                                                   "Same width",
                                                   value = FALSE),
                                     width = 4)),
                               fluidRow(
                                 plotOutput("volcano", height = 600),
                                 downloadButton('downloadVolcano', 'Save volcano')
                               ),
          shinyBS::bsTooltip("volcano_cntrst", "Choose the contrast to plot", "top", options = list(container = "body")),
          shinyBS::bsTooltip("fontsize", "Set the size of name labels", "top", options = list(container = "body")),
          shinyBS::bsTooltip("my_breaks", "Whether set the breaks of x axis by yourself", "top", options = list(container = "body")),
          shinyBS::bsTooltip("mybreaks", "Choose the breaks of x axis, act when [My breaks] is TRUE", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Volcano_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Volcano_Height", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("check_names", "Whether or not to plot names", "top", options = list(container = "body")),
          shinyBS::bsTooltip("p_adj", "Whether or not to use adjusted p values", "top", options = list(container = "body")),
          shinyBS::bsTooltip("same_width", "Whether the x axis to have the same width from 0", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Custom Volcano",
                        fluidRow(
                            column(width = 4,
                             colourInput(inputId = "up_color",
                                                      label = "up",
                                                      showColour = "both",
                                                      palette = "square",
                                                      value = "#CD8500")
                                ),
                            column(width = 4,
                                colourInput(inputId = "down_color",
                                                      label = "down",
                                                      showColour = "both",
                                                      palette = "square",
                                                      value = "#0099CC")
                                ),
                            column(width = 4,
                                colourInput(inputId = "stable_color",
                                                      label = "not significant",
                                                      showColour = "both",
                                                      palette = "square",
                                                      value = "#BEBEBE")
                                )
                            ),
                               fluidRow(
                                 box(uiOutput("Volcano_cntrst"), width = 8),
                                 box(selectizeInput("labelWay",
                                       "label way",
                                       choices = c("all significant", "up","down","selected proteins"),
                                       selected = c("all significant"), multiple = FALSE),width = 4)
                               ),
                               fluidRow(
                                 box(
                                  checkboxInput("if_label_rectangle",
                                                   "label with rectangle",
                                                   value = FALSE),
                                  checkboxInput("if_peptide_color",
                                                   "peptide color",
                                                   value = FALSE),
                                  # column(width = 6,
                                    uiOutput("Peptides1"),
                                    uiOutput("Peptides2"),
                                    # ),
                                  # column(width = 6, uiOutput("Peptides2")),
                                  numericInput("stroke",
                                                  "point outside width",
                                                  min = 0, max = 5, value = 0.6),
                                  uiOutput("selected_proteins"),
                                     numericInput("showNum",
                                                  "show number",
                                                  min = 0, max = 100000, value = 20),
                                     numericInput("fontSize",
                                                  "font size",
                                                  min = 0, max = 20, value = 9),
                                     numericInput("dotsize",
                                                  "dot size",
                                                  min = 0, max = 20, value = 2),
                                     numericInput("custom_volcano_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7),
                                 numericInput("custom_volcano_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 8, collapsible = TRUE, collapsed = TRUE),
                                 box(checkboxInput("P_adj",
                                                   "Adjusted p values",
                                                   value = FALSE),
                                     checkboxInput("Same_width",
                                                   "Same width",
                                                   value = FALSE),
                                     width = 4)),
                               fluidRow(
                                 plotOutput("custom_volcano", height = 600),
                                 downloadButton('download_custom_volcano', 'Save custom volcano')
                               ),
          shinyBS::bsTooltip("up_color", "Set the color of the points those are up-regulated", "top", options = list(container = "body")),
          shinyBS::bsTooltip("down_color", "Set the color of the points those are down-regulated", "top", options = list(container = "body")),
          shinyBS::bsTooltip("stable_color", "Set the color of the points those are not significant", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Volcano_cntrst", "Choose the contrast to plot", "top", options = list(container = "body")),
          shinyBS::bsTooltip("labelWay", "Choose the way of adding labels, one of all significant, up, down, and selected proteins", "top", options = list(container = "body")),
          shinyBS::bsTooltip("if_label_rectangle", "Whether add rectangle underneath the text, making it easier to read", "top", options = list(container = "body")),
          shinyBS::bsTooltip("if_peptide_color", "Whether set the color of the points,(default, purple: Peptides == 1, blue: Peptides == 2)", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Peptides1", "Set the color of the points which the number of Peptides == 1", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Peptides2", "Set the color of the points which the number of Peptides == 2", "top", options = list(container = "body")),
          shinyBS::bsTooltip("stroke", "Set the thickness of black line around the point", "top", options = list(container = "body")),
          shinyBS::bsTooltip("selected_proteins", "Choose the point labels to show, act when [label Way] is selected proteins", "top", options = list(container = "body")),
          shinyBS::bsTooltip("showNum", "Set the number of the labels to add, act when [label way] is significant, up or down. Note that, when it is larger than the number of significant result, it is seted to the number of significant result", "top", options = list(container = "body")),
          shinyBS::bsTooltip("fontSize", "Set the font size of the label", "top", options = list(container = "body")),
          shinyBS::bsTooltip("dotsize", "Set the dot size of the label", "top", options = list(container = "body")),
          shinyBS::bsTooltip("custom_volcano_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("custom_volcano_Height", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("P_adj", "Whether or not to use adjusted p values", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Same_width", "Whether the x axis to have the same width from 0", "top", options = list(container = "body"))
                      )

               ),
               tabBox(title = "QC Plots", width = 12,
                      tabPanel(title = "Pca plot",
                               fluidRow(
                                box(selectizeInput("Indicate",
                                       "Color and shape",
                                       choices = c("condition", "replicate", "Condition", "Replicate"),
                                       selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                box(checkboxInput("if_square",
                                                   "if square",
                                                   value = FALSE),
                                     width = 3)),
                               fluidRow(
                                box(numericInput("pca_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                                box(numericInput("pca_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                                ),
                               fluidRow(
                                 plotOutput("pca", height = 600),
                                 downloadButton('downloadPca', 'Save')
                               ),
          shinyBS::bsTooltip("Indicate", "Set the color, shape and facet_wrap of the plot", "top", options = list(container = "body")),
          shinyBS::bsTooltip("if_square", "Whether x limit is equal to y limit", "top", options = list(container = "body")),
          shinyBS::bsTooltip("pca_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("pca_Height", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      #tabPanel(title = "Pca plot",
                      #         plotOutput("pca", height = 600),
                      #         downloadButton('downloadPca', 'Save')
                      #)
                       tabPanel(title = "Pearson correlation",
                               fluidRow(
                                box(selectizeInput("Pearson_pal",
                                       "color panel",
                                       choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Blues",  "BuGn", "BuPu", "GnBu", "Greens" , "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
                                       selected = c("PRGn"), multiple = FALSE),
                                    checkboxInput("Pearson_pal_rev",
                                                   "pal rev",
                                                   value = FALSE), width = 4),
                                box(numericInput("Pearson_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7),
                                    numericInput("Pearson_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 4),
                                box(numericInput("Pearson_lower",
                                                  "lower",
                                                  min = -1, max = 1, value = -1),
                                    numericInput("Pearson_upper",
                                                  "upper",
                                                  min = -1, max = 1, value = 1), width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("add_values_for_DEP_person",
                                                   "Add values",
                                                   value = FALSE),
                                    numericInput("value_size_for_DEP_person",
                                                  "Value size",
                                                  min = 1, max = 30, value = 10),
                                    numericInput("value_digits_for_DEP_person",
                                                  "Value digits",
                                                  min = 1, max = 30, value = 2), width = 12, collapsible = TRUE, collapsed = TRUE)
                                ),
                               fluidRow(
                                 plotOutput("Pearson_correlation", height = 600),
                                 downloadButton('download_Pearson_correlation', 'Save')
                               ),
          shinyBS::bsTooltip("Pearson_pal", "Set the color panel (from RColorBrewer)", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Pearson_pal_rev", "Whether or not to invert the color palette", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Pearson_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Pearson_Height", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Pearson_lower", "Set the lower limit of the color scale", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Pearson_upper", "Set the upper limit of the color scale", "top", options = list(container = "body"))
                      ),
                        tabPanel(title = "Gower's distance",
                             fluidRow(
                                box(selectizeInput("Gower_pal",
                                       "color panel",
                                       choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Blues",  "BuGn", "BuPu", "GnBu", "Greens" , "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
                                       selected = c("YlOrRd"), multiple = FALSE),
                                   checkboxInput("Gower_pal_rev",
                                                   "pal rev",
                                                   value = TRUE), width = 6),
                                box(numericInput("Gower_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7),
                                    numericInput("Gower_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                                ),
                             fluidRow(
                                box(checkboxInput("add_values_for_DEP_gower",
                                                   "Add values",
                                                   value = FALSE),
                                    numericInput("value_size_for_DEP_gower",
                                                  "Value size",
                                                  min = 1, max = 30, value = 10),
                                    numericInput("value_digits_for_DEP_gower",
                                                  "Value digits",
                                                  min = 1, max = 30, value = 2), width = 12, collapsible = TRUE, collapsed = TRUE)
                                ),
                               fluidRow(
                                 plotOutput("Gowers_distance", height = 600),
                                 downloadButton('download_Gowers_distance', 'Save')
                               ),
          shinyBS::bsTooltip("Gower_pal", "Set the color panel (from RColorBrewer)", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Gower_pal_rev", "Whether or not to invert the color palette", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Gower_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Gower_Height", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                        tabPanel(title = "Sample CVs",
                             fluidRow(
                                box(numericInput("Sample_CVs_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7),
                                    width = 6),
                                box(
                                  numericInput("Sample_CVs_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6
                                  )
                                ),
                               fluidRow(
                                 plotOutput("Sample_CVs", height = 600),
                                 downloadButton('download_Sample_CVs', 'Save')
                               ),
          shinyBS::bsTooltip("Sample_CVs_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Sample_CVs_Height", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Protein Numbers",
                        fluidRow(
                          box(numericInput("numbers_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                          box(numericInput("numbers_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                          ),
                        fluidRow(
                          plotOutput("numbers", height = 600),
                          downloadButton('downloadNumbers', 'Save')
                          ),
          shinyBS::bsTooltip("numbers_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("numbers_Height", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Sample coverage",
                        fluidRow(
                          box(numericInput("coverage_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                          box(numericInput("coverage_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                          ),
                        fluidRow(
                            plotOutput("coverage", height = 600),
                            downloadButton('downloadCoverage', 'Save')

                          ),
          shinyBS::bsTooltip("coverage_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("coverage_Height", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Normalization",
                        fluidRow(
                          box(numericInput("norm_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                          box(numericInput("norm_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                          ),
                        fluidRow(
                            plotOutput("norm", height = 600),
                            downloadButton('downloadNorm', 'Save')

                          ),
          shinyBS::bsTooltip("norm_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("norm_Height", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Missing values - Quant",
                        fluidRow(
                          box(numericInput("detect_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                          box(numericInput("detect_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                          ),
                        fluidRow(
                            plotOutput("detect", height = 600),
                            downloadButton('downloadDetect', 'Save')
                          ),
          shinyBS::bsTooltip("detect_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("detect_Height", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Missing values - Heatmap",
                        fluidRow(
                                 box(
                                 checkboxInput("cluster_columns_for_missval_heatmap",
                                                   "Cluster columns",
                                                   value = TRUE), width = 4),
                                 box(uiOutput("Custom_columns_order_for_missval_heatmap"),
                                     width = 8)
                               ),
                        fluidRow(
                                box(numericInput("missval_heatmap_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                                 box(numericInput("missval_heatmap_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 10), width = 6)),
                        fluidRow(
                          plotOutput("missval", height = 600),
                               downloadButton('downloadMissval', 'Save')
                          ),
          shinyBS::bsTooltip("cluster_columns_for_missval_heatmap", "Whether make cluster on columns", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Custom_columns_order_for_missval_heatmap", "Order of column, act when [Cluster columns] is FALSE", "top", options = list(container = "body")),
          shinyBS::bsTooltip("missval_heatmap_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("missval_heatmap_Height", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Imputation",
                        fluidRow(
                          box(numericInput("imputation_Width",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                          box(numericInput("imputation_Height",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                          ),
                        fluidRow(
                          plotOutput("imputation", height = 600),
                          downloadButton('downloadImputation', 'Save')
                          ),
          shinyBS::bsTooltip("imputation_Width", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("imputation_Height", "Height of the figure to export", "top", options = list(container = "body"))
                      )
               )
        )
      )),
        tabItem(#** Annotation tabItem----
          tabName = "Annotation_tab",
          #         tags$head(tags$style(type="text/css", "
      #                  #loadmessage {
      #                  top: 0px; left: 0px;
      #                  width: 100%; padding: 5px 0px 5px 0px;
      #                  text-align: center; font-weight: bold;
      #                  font-size: 100%; color: #000000;
      #                  background-color: #FFC1C1; z-index: 105;}")), ## 提示条的样式
      # conditionalPanel(condition="$('html').hasClass('shiny-busy')",
      #                  tags$div("calculating...please wait...",id="loadmessage")),
          fluidRow(
              box(title = "Annotation Table",
                   DT::dataTableOutput("annotation_Table"), width = 12)
              # DT::dataTableOutput("annotation_Table")
                          ),
          shinyBS::bsTooltip("annotation_Table", "The table of annotation result", "top", options = list(container = "body"))
          ),
             tabItem(tabName = "GO_tab",
            fluidRow(
        box(selectizeInput("go_ont",
                          "Ontology",
                          choices = c("BP", "CC", "MF", "ALL"),
                          selected = c("BP"), multiple = FALSE),
            width = 2),
        box(numericInput("go_p",
                         "P value",
                         min = 0, max = 1, value = 0.05),
            width = 2),
        box(numericInput("go_padj",
                         "adj. P value",
                         min = 0, max = 1, value = 0.05),
            width = 2),
        box(numericInput("go_qvalue",
                         "q value",
                         min = 0, max = 1, value = 0.2),
            width = 2),
        box(radioButtons("go_color",
                         "colorBy",
                         c("pvalue", "p.adjust"),
                         selected = "p.adjust"), width = 4),
        # infoBoxOutput("significantBox_for_go", width = 4)
          shinyBS::bsTooltip("go_ont", "One of BP, MF, and CC subontologies, or ALL for all three", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_qvalue", "qvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
      ),
          fluidRow(
            column(width = 7,
              fluidRow(
                box(checkboxInput("go_simplify", "removed redundancy of enriched GO terms", value = FALSE), width = 5),
                # box(radioButtons("go_color",
                #          "colorBy",
                #          c("pvalue", "p.adjust"),
                #          selected = "p.adjust"), width = 4),
                infoBoxOutput("significantBox_for_go", width = 7)
                ),
              fluidRow(
                 box(title = "GO Top Table",
                   DT::dataTableOutput("go_Table"), width = 12)
                ),
          shinyBS::bsTooltip("go_simplify", "Whether simplify output by removing redundancy of enriched GO terms", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_Table", "The table of the significant result", "top", options = list(container = "body"))
              ),
            column(width = 5,
              tabBox(title = "Result Plots", width = 12,
                       tabPanel(title = "Bar plot",
                               fluidRow(
                                box(numericInput("go_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("go_wide_bar", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("go_high_bar", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("go_bar_if_Split_for_ont_ALL",
                                                   "If splited by ontology",
                                                   value = FALSE), width = 4)
                                ),
                               fluidRow(
                                 plotOutput("go_barplot", height = 500),
                                 downloadButton('download_go_barplot', 'Save')
                               ),
          shinyBS::bsTooltip("go_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_high_bar", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_bar_if_Split_for_ont_ALL", "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each BP, CC and MF", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Dot plot",
                               fluidRow(
                                box(numericInput("go_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("go_wide_dot", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("go_high_dot", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("go_dot_if_Split_for_ont_ALL",
                                                   "If splited by ontology",
                                                   value = FALSE), width = 4)
                                ),
                               fluidRow(
                                 plotOutput("go_dotplot", height = 500),
                                 downloadButton('download_go_dotplot', 'Save')
                               ),
          shinyBS::bsTooltip("go_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_high_dot", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_dot_if_Split_for_ont_ALL", "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each BP, CC and MF", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Dot opt",
                               fluidRow(
                                box(numericInput("go_ShowCategory_dot_opt", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("go_wide_dot_opt", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("go_high_dot_opt", "Height", min = 1, max = 50, value = 5.7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("go_dotplot_opt", height = 500),
                                 downloadButton('download_go_dotplot_opt', 'Save')
                               ),
          shinyBS::bsTooltip("go_ShowCategory_dot_opt", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_wide_dot_opt", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_high_dot_opt", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Heatplot",
                               fluidRow(
                                box(numericInput("go_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("go_wide_heat", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("go_high_heat", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("go_heatplot", height = 500),
                                 downloadButton('download_go_heatplot', 'Save')
                               ),
          shinyBS::bsTooltip("go_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Cnetplot",
                               fluidRow(
                                box(numericInput("go_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                     width = 4),
                                box(numericInput("go_wide_cnet", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("go_high_cnet", "Height", min = 1, max = 50, value = 6),
                                     width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("go_circular_cnet",
                                                   "Circular",
                                                   value = TRUE),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("go_cnetplot", height = 500),
                                 downloadButton('download_go_cnetplot', 'Save')
                               ),
          shinyBS::bsTooltip("go_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Emaplot",
                               fluidRow(
                                box(numericInput("go_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("go_wide_ema", "Width", min = 1, max = 50, value = 11),
                                     width = 4),
                                box(numericInput("go_high_ema", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("go_emaplot", height = 500),
                                 downloadButton('download_go_emaplot', 'Save')
                               ),
          shinyBS::bsTooltip("go_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                              tabPanel(title = "Goplot",
                               fluidRow(
                                box(numericInput("go_ShowCategory_go", "ShowCategory", min = 1, max = 100, value = 10),
                                     width = 4),
                                box(numericInput("go_wide_go", "Width", min = 1, max = 50, value = 11),
                                     width = 4),
                                box(numericInput("go_high_go", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("go_circular_go",
                                                   "Circular",
                                                   value = TRUE),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("go_goplot", height = 500),
                                 downloadButton('download_go_goplot', 'Save')
                               ),
          shinyBS::bsTooltip("go_ShowCategory_go", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_wide_go", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_high_go", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_circular_go", "whether using circular layout", "top", options = list(container = "body"))
                      ),
                              tabPanel(title = "GOgraph",
                               fluidRow(
                                box(numericInput("go_ShowCategory_GOgraph", "FirstSigNodes", min = 1, max = 100, value = 10),
                                     width = 4),
                                box(numericInput("go_wide_GOgraph", "Width", min = 1, max = 50, value = 7),
                                     width = 4),
                                box(numericInput("go_high_GOgraph", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(actionButton("plot_for_GOgraphplot", "Plot", width = "10%", icon = icon("caret-right"))),#step-forward
                               fluidRow(
                                 plotOutput("go_GOgraphplot", height = 500),
                                 downloadButton('download_go_GOgraphplot', 'Save')
                               ),
          shinyBS::bsTooltip("go_ShowCategory_GOgraph", "Number of significant nodes (retangle nodes in the graph)", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_wide_GOgraph", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("go_high_GOgraph", "Height of the figure to export", "top", options = list(container = "body"))
                      )
                )

              )
                          )
          ),
             tabItem(tabName = "KEGG_tab",
            fluidRow(
        box(numericInput("kegg_p",
                         "P value",
                         min = 0, max = 1, value = 0.05),
            width = 2),
        box(numericInput("kegg_padj",
                         "adj. P value",
                         min = 0, max = 1, value = 0.05),
            width = 2),
        box(numericInput("kegg_qvalue",
                         "q value",
                         min = 0, max = 1, value = 0.2),
            width = 2),
        infoBoxOutput("significantBox_for_kegg", width = 4),
        box(radioButtons("kegg_color",
                         "colorBy",
                         c("pvalue", "p.adjust"),
                         selected = "p.adjust"), width = 2),
          shinyBS::bsTooltip("kegg_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_qvalue", "qvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
      ),
          fluidRow(
            column(width = 7,
              fluidRow(
                 box(title = "KEGG Top Table",
                   DT::dataTableOutput("kegg_Table"), width = 12)
                ),
          shinyBS::bsTooltip("kegg_Table", "The table of the significant result", "top", options = list(container = "body"))
              ),
            column(width = 5,
              tabBox(title = "Result Plots", width = 12,
                       tabPanel(title = "Bar plot",
                               fluidRow(
                                box(numericInput("kegg_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("kegg_wide_bar", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("kegg_high_bar", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("kegg_barplot", height = 500),
                                 downloadButton('download_kegg_barplot', 'Save')
                               ),
          shinyBS::bsTooltip("kegg_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_high_bar", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Dot plot",
                               fluidRow(
                                box(numericInput("kegg_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("kegg_wide_dot", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("kegg_high_dot", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("kegg_dotplot", height = 500),
                                 downloadButton('download_kegg_dotplot', 'Save')
                               ),
          shinyBS::bsTooltip("kegg_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_high_dot", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Dot opt",
                               fluidRow(
                                box(numericInput("kegg_ShowCategory_dot_opt", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("kegg_wide_dot_opt", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("kegg_high_dot_opt", "Height", min = 1, max = 50, value = 5.7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("kegg_dotplot_opt", height = 500),
                                 downloadButton('download_kegg_dotplot_opt', 'Save')
                               ),
          shinyBS::bsTooltip("kegg_ShowCategory_dot_opt", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_wide_dot_opt", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_high_dot_opt", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Heatplot",
                               fluidRow(
                                box(numericInput("kegg_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("kegg_wide_heat", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("kegg_high_heat", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("kegg_heatplot", height = 500),
                                 downloadButton('download_kegg_heatplot', 'Save')
                               ),
          shinyBS::bsTooltip("kegg_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Cnetplot",
                               fluidRow(
                                box(numericInput("kegg_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                     width = 4),
                                box(numericInput("kegg_wide_cnet", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("kegg_high_cnet", "Height", min = 1, max = 50, value = 6),
                                     width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("kegg_circular_cnet",
                                                   "Circular",
                                                   value = TRUE),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("kegg_cnetplot", height = 500),
                                 downloadButton('download_kegg_cnetplot', 'Save')
                               ),
          shinyBS::bsTooltip("kegg_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Emaplot",
                               fluidRow(
                                box(numericInput("kegg_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("kegg_wide_ema", "Width", min = 1, max = 50, value = 11),
                                     width = 4),
                                box(numericInput("kegg_high_ema", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("kegg_emaplot", height = 500),
                                 downloadButton('download_kegg_emaplot', 'Save')
                               ),
          shinyBS::bsTooltip("kegg_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("kegg_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                      )
                )

              )
                          )
          ),
             tabItem(tabName = "Reactome_tab",
          #         tags$head(tags$style(type="text/css", "
      #                  #loadmessage {
      #                  top: 0px; left: 0px;
      #                  width: 100%; padding: 5px 0px 5px 0px;
      #                  text-align: center; font-weight: bold;
      #                  font-size: 100%; color: #000000;
      #                  background-color: #FFC1C1; z-index: 105;}")), ## 提示条的样式
      # conditionalPanel(condition="$('html').hasClass('shiny-busy')",
      #                  tags$div("calculating...please wait...",id="loadmessage")),
            fluidRow(
        box(numericInput("reactome_p",
                         "P value",
                         min = 0, max = 1, value = 0.05),
            width = 2),
        box(numericInput("reactome_padj",
                         "adj. P value",
                         min = 0, max = 1, value = 0.05),
            width = 2),
        box(numericInput("reactome_qvalue",
                         "q value",
                         min = 0, max = 1, value = 0.2),
            width = 2),
        infoBoxOutput("significantBox_for_reactome", width = 4),
        box(radioButtons("reactome_color",
                         "colorBy",
                         c("pvalue", "p.adjust"),
                         selected = "p.adjust"), width = 2),
          shinyBS::bsTooltip("reactome_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_qvalue", "qvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
      ),
          fluidRow(
            column(width = 7,
              fluidRow(
                 box(title = "Reactome Top Table",
                   DT::dataTableOutput("reactome_Table"), width = 12),
          shinyBS::bsTooltip("reactome_Table", "The table of the significant result", "top", options = list(container = "body"))
                )
              ),
            column(width = 5,
              tabBox(title = "Result Plots", width = 12,
                       tabPanel(title = "Bar plot",
                               fluidRow(
                                box(numericInput("reactome_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("reactome_wide_bar", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("reactome_high_bar", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("reactome_barplot", height = 500),
                                 downloadButton('download_reactome_barplot', 'Save')
                               ),
          shinyBS::bsTooltip("reactome_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_high_bar", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Dot plot",
                               fluidRow(
                                box(numericInput("reactome_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("reactome_wide_dot", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("reactome_high_dot", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("reactome_dotplot", height = 500),
                                 downloadButton('download_reactome_dotplot', 'Save')
                               ),
          shinyBS::bsTooltip("reactome_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_high_dot", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Dot opt",
                               fluidRow(
                                box(numericInput("reactome_ShowCategory_dot_opt", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("reactome_wide_dot_opt", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("reactome_high_dot_opt", "Height", min = 1, max = 50, value = 5.7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("reactome_dotplot_opt", height = 500),
                                 downloadButton('download_reactome_dotplot_opt', 'Save')
                               ),
          shinyBS::bsTooltip("reactome_ShowCategory_dot_opt", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_wide_dot_opt", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_high_dot_opt", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Heatplot",
                               fluidRow(
                                box(numericInput("reactome_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("reactome_wide_heat", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("reactome_high_heat", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("reactome_heatplot", height = 500),
                                 downloadButton('download_reactome_heatplot', 'Save')
                               ),
          shinyBS::bsTooltip("reactome_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Cnetplot",
                               fluidRow(
                                box(numericInput("reactome_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                     width = 4),
                                box(numericInput("reactome_wide_cnet", "Width", min = 1, max = 50, value = 10),
                                     width = 4),
                                box(numericInput("reactome_high_cnet", "Height", min = 1, max = 50, value = 6),
                                     width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("reactome_circular_cnet",
                                                   "Circular",
                                                   value = TRUE),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("reactome_cnetplot", height = 500),
                                 downloadButton('download_reactome_cnetplot', 'Save')
                               ),
          shinyBS::bsTooltip("reactome_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Emaplot",
                               fluidRow(
                                box(numericInput("reactome_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("reactome_wide_ema", "Width", min = 1, max = 50, value = 11),
                                     width = 4),
                                box(numericInput("reactome_high_ema", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("reactome_emaplot", height = 500),
                                 downloadButton('download_reactome_emaplot', 'Save')
                               ),
          shinyBS::bsTooltip("reactome_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("reactome_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                      )
                )

              )
                          )
          ),
      tabItem(tabName = "gseGO_tab",
   fluidRow(
        box(selectizeInput("gsego_ont",
                          "Ontology",
                          choices = c("BP", "CC", "MF", "ALL"),
                          selected = c("BP"), multiple = FALSE),
            width = 2),
        box(numericInput("gsego_p",
                         "P value",
                         min = 0, max = 1, value = 0.05),
            width = 2),
        box(numericInput("gsego_padj",
                         "adj. P value",
                         min = 0, max = 1, value = 0.25),
            width = 2),
        box(numericInput("gsego_NES",
                         "NES",
                         min = 0, max = 10, value = 1),
            width = 2),
        box(selectizeInput("gsego_Phenotype",
                          "Phenotype",
                          choices = c("activated", "suppressed"),
                          selected = c("activated", "suppressed"), multiple = TRUE), width = 2),
        box(radioButtons("gsego_color",
                         "colorBy",
                         c("pvalue", "p.adjust"),
                         selected = "p.adjust"), width = 2),
          shinyBS::bsTooltip("gsego_ont", "One of BP, MF, and CC subontologies, or ALL for all three", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_NES", "the |NES| (The absolute value of normalized enrichment score) cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_Phenotype", "The phenotype that you want to show, one or both of activated and suppressed. For example, if your log2 fold change is A_vs_B, then, term enriched in A is activated, and term enriched in B is suppressed", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
      ),
          fluidRow(
            column(width = 7,
              fluidRow(
                box(checkboxInput("gsego_simplify", "removed redundancy of enriched GO terms", value = FALSE), width = 5),
                # box(radioButtons("go_color",
                #          "colorBy",
                #          c("pvalue", "p.adjust"),
                #          selected = "p.adjust"), width = 4),
                infoBoxOutput("significantBox_for_gsego", width = 7)
                ),
              fluidRow(
                 box(title = "gseGO Top Table",
                   DT::dataTableOutput("gsego_Table"), width = 12)
                ),
          shinyBS::bsTooltip("gsego_simplify", "Whether simplify output by removing redundancy of enriched GO terms", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_Table", "The table of the significant result", "top", options = list(container = "body"))
              ),
            column(width = 5,
              tabBox(title = "Result Plots", width = 12,
                       tabPanel(title = "Bar plot",
                               fluidRow(
                                box(numericInput("gsego_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("gsego_wide_bar", "Width", min = 1, max = 50, value = 12),
                                     width = 4),
                                box(numericInput("gsego_high_bar", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("gsego_bar_if_Split_for_ont_ALL",
                                                   "If splited by ontology",
                                                   value = FALSE), width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsego_barplot", height = 500),
                                 downloadButton('download_gsego_barplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsego_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_high_bar", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_bar_if_Split_for_ont_ALL", "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of each phenotype that you selected of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each phenotype that you selected based on each BP, CC and MF", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Dot plot",
                               fluidRow(
                                box(numericInput("gsego_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("gsego_wide_dot", "Width", min = 1, max = 50, value = 12),
                                     width = 4),
                                box(numericInput("gsego_high_dot", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("gsego_dot_if_Split_for_ont_ALL",
                                                   "If splited by ontology",
                                                   value = FALSE), width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsego_dotplot", height = 500),
                                 downloadButton('download_gsego_dotplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsego_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_high_dot", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_dot_if_Split_for_ont_ALL", "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of each phenotype that you selected of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each phenotype that you selected based on each BP, CC and MF", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Heatplot",
                               fluidRow(
                                box(numericInput("gsego_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("gsego_wide_heat", "Width", min = 1, max = 200, value = 80),
                                     width = 4),
                                box(numericInput("gsego_high_heat", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsego_heatplot", height = 500),
                                 downloadButton('download_gsego_heatplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsego_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Cnetplot",
                               fluidRow(
                                box(numericInput("gsego_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                     width = 4),
                                box(numericInput("gsego_wide_cnet", "Width", min = 1, max = 200, value = 30),
                                     width = 4),
                                box(numericInput("gsego_high_cnet", "Height", min = 1, max = 200, value = 25),
                                     width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("gsego_circular_cnet",
                                                   "Circular",
                                                   value = TRUE),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsego_cnetplot", height = 500),
                                 downloadButton('download_gsego_cnetplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsego_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Emaplot",
                               fluidRow(
                                box(numericInput("gsego_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("gsego_wide_ema", "Width", min = 1, max = 50, value = 11),
                                     width = 4),
                                box(numericInput("gsego_high_ema", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsego_emaplot", height = 500),
                                 downloadButton('download_gsego_emaplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsego_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                              tabPanel(title = "Gseaplot",
                               fluidRow(box(uiOutput("gsego_term"), width = 12)),
                               fluidRow(
                                box(numericInput("gsego_wide_Gsea", "Width", min = 1, max = 50, value = 8),
                                     width = 4),
                                box(numericInput("gsego_high_Gsea", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsego_Gseaplot", height = 500),
                                 downloadButton('download_gsego_Gseaplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsego_term", "The term that you want to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_wide_Gsea", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsego_high_Gsea", "Height of the figure to export", "top", options = list(container = "body"))
                      )
                )

              )
                          )
          ),
      tabItem(tabName = "gseKEGG_tab",
   fluidRow(
        box(numericInput("gsekegg_p",
                         "P value",
                         min = 0, max = 1, value = 0.05),
            width = 2),
        box(numericInput("gsekegg_padj",
                         "adj. P value",
                         min = 0, max = 1, value = 0.25),
            width = 2),
        box(numericInput("gsekegg_NES",
                         "NES",
                         min = 0, max = 10, value = 1),
            width = 2),
        box(selectizeInput("gsekegg_Phenotype",
                          "Phenotype",
                          choices = c("activated", "suppressed"),
                          selected = c("activated", "suppressed"), multiple = TRUE), width = 4),
        box(radioButtons("gsekegg_color",
                         "colorBy",
                         c("pvalue", "p.adjust"),
                         selected = "p.adjust"), width = 2),
          shinyBS::bsTooltip("gsekegg_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_NES", "the |NES| (The absolute value of normalized enrichment score) cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_Phenotype", "The phenotype that you want to show, one or both of activated and suppressed. For example, if your log2 fold change is A_vs_B, then, term enriched in A is activated, and term enriched in B is suppressed", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
      ),
          fluidRow(
            column(width = 7,
              fluidRow(
                infoBoxOutput("significantBox_for_gsekegg", width = 7)
                ),
              fluidRow(
                 box(title = "gseKEGG Top Table",
                   DT::dataTableOutput("gsekegg_Table"), width = 12)
                ),
          shinyBS::bsTooltip("gsekegg_Table", "The table of the significant result", "top", options = list(container = "body"))
              ),
            column(width = 5,
              tabBox(title = "Result Plots", width = 12,
                       tabPanel(title = "Bar plot",
                               fluidRow(
                                box(numericInput("gsekegg_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("gsekegg_wide_bar", "Width", min = 1, max = 50, value = 12),
                                     width = 4),
                                box(numericInput("gsekegg_high_bar", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsekegg_barplot", height = 500),
                                 downloadButton('download_gsekegg_barplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsekegg_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_high_bar", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Dot plot",
                               fluidRow(
                                box(numericInput("gsekegg_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("gsekegg_wide_dot", "Width", min = 1, max = 50, value = 12),
                                     width = 4),
                                box(numericInput("gsekegg_high_dot", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsekegg_dotplot", height = 500),
                                 downloadButton('download_gsekegg_dotplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsekegg_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_high_dot", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Heatplot",
                               fluidRow(
                                box(numericInput("gsekegg_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("gsekegg_wide_heat", "Width", min = 1, max = 200, value = 80),
                                     width = 4),
                                box(numericInput("gsekegg_high_heat", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsekegg_heatplot", height = 500),
                                 downloadButton('download_gsekegg_heatplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsekegg_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Cnetplot",
                               fluidRow(
                                box(numericInput("gsekegg_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                     width = 4),
                                box(numericInput("gsekegg_wide_cnet", "Width", min = 1, max = 200, value = 30),
                                     width = 4),
                                box(numericInput("gsekegg_high_cnet", "Height", min = 1, max = 200, value = 25),
                                     width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("gsekegg_circular_cnet",
                                                   "Circular",
                                                   value = TRUE),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsekegg_cnetplot", height = 500),
                                 downloadButton('download_gsekegg_cnetplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsekegg_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Emaplot",
                               fluidRow(
                                box(numericInput("gsekegg_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("gsekegg_wide_ema", "Width", min = 1, max = 50, value = 11),
                                     width = 4),
                                box(numericInput("gsekegg_high_ema", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsekegg_emaplot", height = 500),
                                 downloadButton('download_gsekegg_emaplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsekegg_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                              tabPanel(title = "Gseaplot",
                               fluidRow(box(uiOutput("gsekegg_term"), width = 12)),
                               fluidRow(
                                box(numericInput("gsekegg_wide_Gsea", "Width", min = 1, max = 50, value = 8),
                                     width = 4),
                                box(numericInput("gsekegg_high_Gsea", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsekegg_Gseaplot", height = 500),
                                 downloadButton('download_gsekegg_Gseaplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsekegg_term", "The term that you want to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_wide_Gsea", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsekegg_high_Gsea", "Height of the figure to export", "top", options = list(container = "body"))
                      )
                )

              )
                          )
          ),
      tabItem(tabName = "gseReactome_tab",
   fluidRow(
        box(numericInput("gsereactome_p",
                         "P value",
                         min = 0, max = 1, value = 0.05),
            width = 2),
        box(numericInput("gsereactome_padj",
                         "adj. P value",
                         min = 0, max = 1, value = 0.25),
            width = 2),
        box(numericInput("gsereactome_NES",
                         "NES",
                         min = 0, max = 10, value = 1),
            width = 2),
        box(selectizeInput("gsereactome_Phenotype",
                          "Phenotype",
                          choices = c("activated", "suppressed"),
                          selected = c("activated", "suppressed"), multiple = TRUE), width = 4),
        box(radioButtons("gsereactome_color",
                         "colorBy",
                         c("pvalue", "p.adjust"),
                         selected = "p.adjust"), width = 2),
          shinyBS::bsTooltip("gsereactome_p", "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_padj", "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_NES", "the |NES| (The absolute value of normalized enrichment score) cutoff on enrichment tests to report", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_Phenotype", "The phenotype that you want to show, one or both of activated and suppressed. For example, if your log2 fold change is A_vs_B, then, term enriched in A is activated, and term enriched in B is suppressed", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_color", "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
      ),
          fluidRow(
            column(width = 7,
              fluidRow(
                infoBoxOutput("significantBox_for_gsereactome", width = 7)
                ),
              fluidRow(
                 box(title = "gseReactome Top Table",
                   DT::dataTableOutput("gsereactome_Table"), width = 12)
                ),
          shinyBS::bsTooltip("gsereactome_Table", "The table of the significant result", "top", options = list(container = "body"))
              ),
            column(width = 5,
              tabBox(title = "Result Plots", width = 12,
                       tabPanel(title = "Bar plot",
                               fluidRow(
                                box(numericInput("gsereactome_ShowCategory_bar", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("gsereactome_wide_bar", "Width", min = 1, max = 50, value = 12),
                                     width = 4),
                                box(numericInput("gsereactome_high_bar", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsereactome_barplot", height = 500),
                                 downloadButton('download_gsereactome_barplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsereactome_ShowCategory_bar", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_wide_bar", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_high_bar", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Dot plot",
                               fluidRow(
                                box(numericInput("gsereactome_ShowCategory_dot", "ShowCategory", min = 1, max = 100, value = 20),
                                     width = 4),
                                box(numericInput("gsereactome_wide_dot", "Width", min = 1, max = 50, value = 12),
                                     width = 4),
                                box(numericInput("gsereactome_high_dot", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsereactome_dotplot", height = 500),
                                 downloadButton('download_gsereactome_dotplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsereactome_ShowCategory_dot", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_wide_dot", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_high_dot", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Heatplot",
                               fluidRow(
                                box(numericInput("gsereactome_ShowCategory_heat", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("gsereactome_wide_heat", "Width", min = 1, max = 200, value = 80),
                                     width = 4),
                                box(numericInput("gsereactome_high_heat", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsereactome_heatplot", height = 500),
                                 downloadButton('download_gsereactome_heatplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsereactome_ShowCategory_heat", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_wide_heat", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_high_heat", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Cnetplot",
                               fluidRow(
                                box(numericInput("gsereactome_ShowCategory_cnet", "ShowCategory", min = 1, max = 100, value = 5),
                                     width = 4),
                                box(numericInput("gsereactome_wide_cnet", "Width", min = 1, max = 200, value = 30),
                                     width = 4),
                                box(numericInput("gsereactome_high_cnet", "Height", min = 1, max = 200, value = 25),
                                     width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("gsereactome_circular_cnet",
                                                   "Circular",
                                                   value = TRUE),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsereactome_cnetplot", height = 500),
                                 downloadButton('download_gsereactome_cnetplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsereactome_ShowCategory_cnet", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_wide_cnet", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_high_cnet", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_circular_cnet", "whether using circular layout", "top", options = list(container = "body"))
                      ),
                       tabPanel(title = "Emaplot",
                               fluidRow(
                                box(numericInput("gsereactome_ShowCategory_ema", "ShowCategory", min = 1, max = 100, value = 30),
                                     width = 4),
                                box(numericInput("gsereactome_wide_ema", "Width", min = 1, max = 50, value = 11),
                                     width = 4),
                                box(numericInput("gsereactome_high_ema", "Height", min = 1, max = 50, value = 10),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsereactome_emaplot", height = 500),
                                 downloadButton('download_gsereactome_emaplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsereactome_ShowCategory_ema", "Number of categories to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_wide_ema", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_high_ema", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                              tabPanel(title = "Gseaplot",
                               fluidRow(box(uiOutput("gsereactome_term"), width = 12)),
                               fluidRow(
                                box(numericInput("gsereactome_wide_Gsea", "Width", min = 1, max = 50, value = 8),
                                     width = 4),
                                box(numericInput("gsereactome_high_Gsea", "Height", min = 1, max = 50, value = 7),
                                     width = 4)
                                ),
                               fluidRow(
                                 plotOutput("gsereactome_Gseaplot", height = 500),
                                 downloadButton('download_gsereactome_Gseaplot', 'Save')
                               ),
          shinyBS::bsTooltip("gsereactome_term", "The term that you want to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_wide_Gsea", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("gsereactome_high_Gsea", "Height of the figure to export", "top", options = list(container = "body"))
                      )
                )

              )
                          )
          ),
   tabItem(#** PPI tabItem----
    tabName = "PPI_tab",
            # tags$head(tags$style(type="text/css", "
            #             #loadmessage {
            #             top: 0px; left: 0px;
            #             width: 100%; padding: 5px 0px 5px 0px;
            #             text-align: center; font-weight: bold;
            #             font-size: 100%; color: #000000;
            #             background-color: #FFC1C1; z-index: 105;}")), ## 提示条的样式
            # conditionalPanel(condition="$('html').hasClass('shiny-busy')",
            #                  tags$div("calculating...please wait...",id="loadmessage")),
            fluidRow(
              column(width = 12,
                     tabBox(width = 12,
                            tabPanel(
                              tags$head(tags$style(type="text/css", "
                        #loadmessage {
                        top: 0px; left: 0px;
                        width: 100%; padding: 5px 0px 5px 0px;
                        text-align: center; font-weight: bold;
                        font-size: 100%; color: #000000;
                        background-color: #FFC1C1; z-index: 105;}")), ## 提示条的样式
                            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                             tags$div("calculating...please wait...",id="loadmessage")),
                              fluidRow(height=30,
                                       column(width = 3,
                                        colourInput(inputId = "nodecolor",
                                                      label = "input node color(like:#2EA9DF)",
                                                      showColour = "both",
                                                      palette = "square",
                                                      value = "#2EA9DF")
                                              # textInput("nodecolor", "input node color(like:#2EA9DF)", "#2EA9DF")
                                              ),
                                       column(width = 3,
                                              checkboxInput("changesize",
                                                            "whether change nodesizes by connection",
                                                            value = FALSE)),
                                       column(width = 3,
                                              numericInput("nodesize",
                                                           "basic node size",
                                                           min = 10, max = 150, value = 20)),
                                       column(width = 3,
                                              numericInput("PPIfontsize",
                                                           "fonts size",
                                                           min = 10, max = 150, value = 25))
                              ),
                              fluidRow(height=30,
                                       column(width = 3,
                                        colourInput(inputId = "linecolor",
                                                      label = "input line color(like:  #ADD8E6)",
                                                      showColour = "both",
                                                      palette = "square",
                                                      value = "#ADD8E6")
                                              # textInput("linecolor", "input line color(like:  #ADD8E6)", "#ADD8E6")
                                              ),

                                       column(width = 3,
                                              checkboxInput("changewidth",
                                                            "whether change link width by connection",
                                                            value = FALSE)),
                                       column(width = 3,
                                              numericInput("linewidth",
                                                           "basic line width",
                                                           min = 3, max = 80, value = 5)),
                                       column(width = 3,
                                              selectizeInput("layoutway",
                                                             "choose layout style of network",
                                                             choices=c("layout_as_star","layout_components","layout_in_circle", "layout_nicely" ,"layout_on_grid","layout_on_sphere",
                                                                       "layout_randomly","layout_with_dh","layout_with_drl" ,"layout_with_fr","layout_with_gem","layout_with_graphopt",
                                                                       "layout_with_kk","layout_with_lgl"  ),
                                                             selected="layout_components"))
                              ),
                              fluidRow(height=30,
                                       column(width = 2,
                                              selectizeInput("nodeshape",
                                                             "choose node shape",
                                                             choices=c("square", "triangle", "box", "circle", "dot", "star",
                                                                       "ellipse", "database", "text", "diamond"),
                                                             selected="dot")),
                                       column(width = 2,
                                              checkboxInput("highlightkey",
                                                            "highlight key nodes?",
                                                            value = TRUE)),
                                       column(width = 2,
                                              checkboxInput("smoothline",
                                                            "smooth lines",
                                                            value = FALSE)),
                                       column(width = 2,
                                              selectizeInput("smoothtype",
                                                             "type to smooth the line",
                                                             choices=c('dynamic', 'continuous', 'discrete', 'diagonalCross', 'straightCross', 'horizontal',
                                                                       'vertical', 'curvedCW', 'curvedCCW', 'cubicBezier'),
                                                             selected="continuous")),
                                       column(width = 3,
                                              numericInput("plothigh",
                                                           "plot high",
                                                           min = 10, max = 200, value = 75))
                              ),
                              fluidRow(
                                # visNetworkOutput("network",height=input$plothigh*10)
                                uiOutput("reactome_plot")
                              ),
                              fluidRow(
                                DT::dataTableOutput("String_Table")
                              ),
          shinyBS::bsTooltip("nodecolor", "choose the color of input node", "top", options = list(container = "body")),
          shinyBS::bsTooltip("changesize", "whether change nodesizes by connection", "top", options = list(container = "body")),
          shinyBS::bsTooltip("nodesize", "basic node size", "top", options = list(container = "body")),
          shinyBS::bsTooltip("PPIfontsize", "fonts size", "top", options = list(container = "body")),
          shinyBS::bsTooltip("linecolor", "choose the color of input line", "top", options = list(container = "body")),
          shinyBS::bsTooltip("changewidth", "whether change link width by connection", "top", options = list(container = "body")),
          shinyBS::bsTooltip("linewidth", "basic line width", "top", options = list(container = "body")),
          shinyBS::bsTooltip("layoutway", "choose layout style of network", "top", options = list(container = "body")),
          shinyBS::bsTooltip("nodeshape", "choose node shape", "top", options = list(container = "body")),
          shinyBS::bsTooltip("highlightkey", "highlight key nodes?", "top", options = list(container = "body")),
          shinyBS::bsTooltip("smoothline", "smooth lines", "top", options = list(container = "body")),
          shinyBS::bsTooltip("smoothtype", "type to smooth the line", "top", options = list(container = "body")),
          shinyBS::bsTooltip("plothigh", "the height of the plot region", "top", options = list(container = "body")),
          shinyBS::bsTooltip("String_Table", "the table of PPI (Protein-Protein Interaction)", "top", options = list(container = "body"))
                            )
                     )

              )
            )
   ),
      tabItem(#** RNAseq tabItem----
        tabName = "RNAseq_tab",
      fluidRow(
        box(numericInput("p_for_DERNAseq",
                         "adj. P value",
                         min = 0.0001, max = 0.1, value = 0.05),
            width = 2),
        box(numericInput("lfc_for_DERNAseq",
                         "Log2 fold change",
                         min = 0, max = 10, value = 1),
            width = 2),
        infoBoxOutput("significantBox_for_DERNAseq"),
        box(radioButtons("pres_for_RNAseq",
                         "Data presentation",
                         c("raw", "centered"),
                         selected = "centered"),
            width = 2),
        box(radioButtons(inputId = "contrasts_for_RNAseq",
                         label = "Contrasts",
                         choices = c("control", "all", "manual"),
                         selected = "control", inline = TRUE),
            width = 2),
          shinyBS::bsTooltip("p_for_DERNAseq", "Set the false discovery rate threshold", "top", options = list(container = "body")),
          shinyBS::bsTooltip("lfc_for_DERNAseq", "Set the log2 fold change threshold", "top", options = list(container = "body")),
          shinyBS::bsTooltip("pres_for_RNAseq", "The type of data scaling used for heatmap plotting. Either the raw or the centered log2-intensity (centered)", "top", options = list(container = "body")),
          shinyBS::bsTooltip("contrasts_for_RNAseq", "The type of contrasts that will be tested. This can be all possible pairwise comparisons (all), limited to the comparisons versus the control (control), or manually defined contrasts (manual)", "top", options = list(container = "body"))
      ),
      fluidRow(
        column(width = 7,
               box(title = "Top Table",
                   DT::dataTableOutput("table_for_DERNAseq"), width = 12),
          shinyBS::bsTooltip("table_for_DERNAseq", "Table of the significant results", "top", options = list(container = "body"))
        ),
        column(width = 5,
               tabBox(title = "Result Plots", width = 12,
                    tabPanel(title = "Diagno dispests",
                         fluidRow(
                          box(numericInput("Diagno_Width_for_RNAseq",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                          box(numericInput("Diagno_Height_for_RNAseq",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                          ),
                               fluidRow(
                                 plotOutput("diagno_dispests_for_DERNAseq", height = 600),
                                 downloadButton('download_diagno_dispests_for_DERNAseq', 'Save')
                               ),
          shinyBS::bsTooltip("Diagno_Width_for_RNAseq", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Diagno_Height_for_RNAseq", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                                           tabPanel(title = "Pca plot",
                               fluidRow(
                                box(selectizeInput("Indicate_for_RNAseq",
                                       "Color and shape",
                                       choices = c("condition", "replicate", "Condition", "Replicate"),
                                       selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                box(checkboxInput("if_square_for_RNAseq",
                                                   "if square",
                                                   value = FALSE),
                                     width = 3)),
                               fluidRow(
                                box(numericInput("pca_Width_for_RNAsqe",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                                box(numericInput("pca_Height_for_RNAseq",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                                ),
                               fluidRow(
                                 plotOutput("pca_for_RNAseq", height = 600),
                                 downloadButton('downloadPca_for_RNAseq', 'Save')
                               ),
          shinyBS::bsTooltip("Indicate_for_RNAseq", "Set the color, shape and facet_wrap of the plot", "top", options = list(container = "body")),
          shinyBS::bsTooltip("if_square_for_RNAseq", "Whether x limit is equal to y limit", "top", options = list(container = "body")),
          shinyBS::bsTooltip("pca_Width_for_RNAsqe", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("pca_Height_for_RNAseq", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      #tabPanel(title = "Pca plot",
                      #         plotOutput("pca", height = 600),
                      #         downloadButton('downloadPca', 'Save')
                      #)
                       tabPanel(title = "Pearson correlation",
                               fluidRow(
                                box(selectizeInput("Pearson_pal_for_RNAseq",
                                       "color panel",
                                       choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Blues",  "BuGn", "BuPu", "GnBu", "Greens" , "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
                                       selected = c("PRGn"), multiple = FALSE),
                                    checkboxInput("Pearson_pal_rev_for_RNAseq",
                                                   "pal rev",
                                                   value = FALSE), width = 4),
                                box(numericInput("Pearson_Width_for_RNAseq",
                                                  "width",
                                                  min = 1, max = 30, value = 7),
                                    numericInput("Pearson_Height_for_RNAseq",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 4),
                                box(numericInput("Pearson_lower_for_RNAseq",
                                                  "lower",
                                                  min = -1, max = 1, value = -1),
                                    numericInput("Pearson_upper_for_RNAseq",
                                                  "upper",
                                                  min = -1, max = 1, value = 1), width = 4)
                                ),
                               fluidRow(
                                box(checkboxInput("add_values_for_RNAseq_person",
                                                   "Add values",
                                                   value = FALSE),
                                    numericInput("value_size_for_RNAseq_person",
                                                  "Value size",
                                                  min = 1, max = 30, value = 10),
                                    numericInput("value_digits_for_RNAseq_person",
                                                  "Value digits",
                                                  min = 1, max = 30, value = 2), width = 12, collapsible = TRUE, collapsed = TRUE)
                                ),
                               fluidRow(
                                 plotOutput("Pearson_correlation_for_RNAseq", height = 600),
                                 downloadButton('download_Pearson_correlation_for_RNAseq', 'Save')
                               ),
          shinyBS::bsTooltip("Pearson_pal_for_RNAseq", "Set the color panel (from RColorBrewer)", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Pearson_pal_rev_for_RNAseq", "Whether or not to invert the color palette", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Pearson_Width_for_RNAseq", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Pearson_Height_for_RNAseq", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Pearson_lower_for_RNAseq", "Set the lower limit of the color scale", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Pearson_upper_for_RNAseq", "Set the upper limit of the color scale", "top", options = list(container = "body"))
                      ),
                        tabPanel(title = "Gower's distance",
                             fluidRow(
                                box(selectizeInput("Gower_pal_for_RNAseq",
                                       "color panel",
                                       choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Blues",  "BuGn", "BuPu", "GnBu", "Greens" , "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
                                       selected = c("YlOrRd"), multiple = FALSE),
                                   checkboxInput("Gower_pal_rev_for_RNAseq",
                                                   "pal rev",
                                                   value = TRUE), width = 6),
                                box(numericInput("Gower_Width_for_RNAseq",
                                                  "width",
                                                  min = 1, max = 30, value = 7),
                                    numericInput("Gower_Height_for_RNAseq",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                                ),
                             fluidRow(
                                box(checkboxInput("add_values_for_RNAseq_gower",
                                                   "Add values",
                                                   value = FALSE),
                                    numericInput("value_size_for_RNAseq_gower",
                                                  "Value size",
                                                  min = 1, max = 30, value = 10),
                                    numericInput("value_digits_for_RNAseq_gower",
                                                  "Value digits",
                                                  min = 1, max = 30, value = 2), width = 12, collapsible = TRUE, collapsed = TRUE)
                                ),
                               fluidRow(
                                 plotOutput("Gowers_distance_for_RNAseq", height = 600),
                                 downloadButton('download_Gowers_distance_for_RNAseq', 'Save')
                               ),
          shinyBS::bsTooltip("Gower_pal_for_RNAseq", "Set the color panel (from RColorBrewer)", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Gower_pal_rev_for_RNAseq", "Whether or not to invert the color palette", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Gower_Width_for_RNAseq", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Gower_Height_for_RNAseq", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Heatmap",
                               fluidRow(
                                column(width = 6,
                                    box(checkboxInput("manual_heatmap_for_RNAseq",
                                                   "Manual heatmap",
                                                   value = FALSE), width = 12)
                                    ),
                                column(width = 6,
                                    box(selectizeInput("colorbar_for_RNAseq",
                                       "colorbar",
                                       choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                                       selected = c("RdBu"), multiple = FALSE), width = 12)
                                    )
                                 ),
                                fluidRow(
                                    box(uiOutput("heatmap_cntrst_for_RNAseq"), width = 12)
                                    ),
                               # fluidRow(
                               #   box(checkboxInput("manual_heatmap_for_RNAseq",
                               #                     "Manual heatmap",
                               #                     value = FALSE),
                               #   # checkboxInput("same_trend",
                               #   #                   "Same trend heatmap",
                               #   #                   value = FALSE),
                               #    uiOutput("heatmap_cntrst_for_RNAseq"), width = 12, collapsible = FALSE, collapsed = FALSE)
                               #   #box(checkboxInput("manual_heatmap",
                               #   #                  "Manual heatmap",
                               #   #                  value = FALSE)),
                               # ),
                               fluidRow(
                                box(numericInput("k_for_RNAseq",
                                                  "Kmeans",
                                                  min = 0, max = 15, value = 2),
                                     width = 4),
                                box(numericInput("limit_for_RNAseq",
                                                  "Color",
                                                  min = 0, max = 16, value = 4),
                                     width = 4),
                                box(numericInput("size_for_RNAseq",
                                                  "Heat size(4-30)",
                                                  min = 4, max = 30, value = 10),
                                     width = 4)

                                ),
                               fluidRow(
                                 box(
                                  numericInput("row_font_size_for_RNAseq",
                                                  "row font size",
                                                  min = 0, max = 15, value = 6),
                                 numericInput("col_font_size_for_RNAseq",
                                                  "col font size",
                                                  min = 0, max = 16, value = 10),
                                 numericInput("Width_for_RNAseq",
                                                  "heatmap width",
                                                  min = 1, max = 30, value = 7),
                                 numericInput("Height_for_RNAseq",
                                                  "heatmap height",
                                                  min = 1, max = 30, value = 10),
                                 checkboxInput("cluster_columns_for_RNAseq",
                                                   "Cluster columns",
                                                   value = TRUE),
                                 uiOutput("Custom_columns_order_for_RNAseq"),
                                 checkboxInput("if_mysplit_for_RNAseq",
                                                   "If mysplit",
                                                   value = FALSE),
                                  uiOutput("mysplit_for_RNAseq"),
                                     width = 12, collapsible = TRUE, collapsed = TRUE)
                               ),
                               fluidRow(
                                  uiOutput("plot_heatmap_for_RNAseq"),
                                 downloadButton('downloadHeatmap_for_DERNAseq', 'Save heatmap')),
          shinyBS::bsTooltip("colorbar_for_RNAseq", "Choose the color palette of the colorbar", "top", options = list(container = "body")),
          shinyBS::bsTooltip("manual_heatmap_for_RNAseq", "Whether plot the manual heatmap", "top", options = list(container = "body")),
          shinyBS::bsTooltip("heatmap_cntrst_for_RNAseq", "Choose the contrasts that you want to show, act when [Manual heatmap] is TRUE, can be one or more contrasts", "top", options = list(container = "body")),
          shinyBS::bsTooltip("k_for_RNAseq", "Set the number of k-means clusters", "top", options = list(container = "body")),
          shinyBS::bsTooltip("limit_for_RNAseq", "Set the outer limits of the color scale", "top", options = list(container = "body")),
          shinyBS::bsTooltip("size_for_RNAseq", "Set the height of the plot region", "top", options = list(container = "body")),
          shinyBS::bsTooltip("row_font_size_for_RNAseq", "Set the size of row labels", "top", options = list(container = "body")),
          shinyBS::bsTooltip("col_font_size_for_RNAseq", "Set the size of column labels", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Width_for_RNAseq", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Height_for_RNAseq", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("cluster_columns_for_RNAseq", "Whether make cluster on columns", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Custom_columns_order_for_RNAseq", "Order of column, act when [Cluster columns] is FALSE", "top", options = list(container = "body")),
          shinyBS::bsTooltip("if_mysplit_for_RNAseq", "Whether order the splices manually", "top", options = list(container = "body")),
          shinyBS::bsTooltip("mysplit_for_RNAseq", "The ordered splices of the heatmap (from top to bottom), only act when both [If mysplit] and [Kmeans] are TRUE", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Custom Volcano",
                        fluidRow(
                            column(width = 4,
                             colourInput(inputId = "up_color_for_RNAseq",
                                                      label = "up",
                                                      showColour = "both",
                                                      palette = "square",
                                                      value = "#CD8500")
                                ),
                            column(width = 4,
                                colourInput(inputId = "down_color_for_RNAseq",
                                                      label = "down",
                                                      showColour = "both",
                                                      palette = "square",
                                                      value = "#0099CC")
                                ),
                            column(width = 4,
                                colourInput(inputId = "stable_color_for_RNAseq",
                                                      label = "not significant",
                                                      showColour = "both",
                                                      palette = "square",
                                                      value = "#BEBEBE")
                                )
                            ),
                               fluidRow(
                                 box(selectizeInput("labelWay_for_RNAseq",
                                       "label way",
                                       choices = c("all significant", "up","down","selected proteins"),
                                       selected = c("all significant"), multiple = FALSE),width = 4),
                                 box(uiOutput("contrast_for_RNAseq_customvolcano"), width = 8)
                               ),
                               fluidRow(
                                 box(
                                  checkboxInput("if_label_rectangle_for_RNAseq",
                                                   "label with rectangle",
                                                   value = FALSE),
                                  numericInput("stroke_for_RNAseq",
                                                  "point outside width",
                                                  min = 0, max = 5, value = 0.6),
                                  checkboxInput("labeled_for_RNAseq_customvolcano",
                                                   "labeled",
                                                   value = FALSE),
                                  uiOutput("selected_proteins_for_RNAseq"),
                                     numericInput("showNum_for_RNAseq",
                                                  "show number",
                                                  min = 0, max = 100000, value = 20),
                                     numericInput("fontSize_for_RNAseq",
                                                  "font size",
                                                  min = 0, max = 20, value = 9),
                                     numericInput("dotsize_for_RNAseq",
                                                  "dot size",
                                                  min = 0, max = 20, value = 2),
                                     numericInput("custom_volcano_Width_for_RNAseq",
                                                  "width",
                                                  min = 1, max = 30, value = 7),
                                 numericInput("custom_volcano_Height_for_RNAseq",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 8, collapsible = TRUE, collapsed = TRUE),
                                 box(checkboxInput("P_adj_for_RNAseq",
                                                   "Adjusted p values",
                                                   value = FALSE),
                                     checkboxInput("Same_width_for_RNAseq",
                                                   "Same width",
                                                   value = FALSE),
                                     width = 4)),
                               fluidRow(
                                 plotOutput("custom_volcano_for_RNAseq", height = 600),
                                 downloadButton('download_custom_volcano_for_RNAseq', 'Save custom volcano')
                               ),
          shinyBS::bsTooltip("up_color_for_RNAseq", "Set the color of the points those are up-regulated", "top", options = list(container = "body")),
          shinyBS::bsTooltip("down_color_for_RNAseq", "Set the color of the points those are down-regulated", "top", options = list(container = "body")),
          shinyBS::bsTooltip("stable_color_for_RNAseq", "Set the color of the points those are not significant", "top", options = list(container = "body")),
          shinyBS::bsTooltip("labelWay_for_RNAseq", "Choose the way of adding labels, one of all significant, up, down, and selected proteins", "top", options = list(container = "body")),
          shinyBS::bsTooltip("contrast_for_RNAseq_customvolcano", "Choose the contrast that you want to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("if_label_rectangle_for_RNAseq", "Whether add rectangle underneath the text, making it easier to read", "top", options = list(container = "body")),
          shinyBS::bsTooltip("stroke_for_RNAseq", "Set the thickness of black line around the point", "top", options = list(container = "body")),
          shinyBS::bsTooltip("labeled_for_RNAseq_customvolcano", "Whether use the labelway : selected proteins", "top", options = list(container = "body")),
          shinyBS::bsTooltip("selected_proteins_for_RNAseq", "Choose the point labels to show, act when [label Way] is selected proteins", "top", options = list(container = "body")),
          shinyBS::bsTooltip("showNum_for_RNAseq", "Set the number of the labels to add, act when [label way] is significant, up or down. Note that, when it is larger than the number of significant result, it is seted to the number of significant result", "top", options = list(container = "body")),
          shinyBS::bsTooltip("fontSize_for_RNAseq", "Set the font size of the label", "top", options = list(container = "body")),
          shinyBS::bsTooltip("dotsize_for_RNAseq", "Set the dot size of the label", "top", options = list(container = "body")),
          shinyBS::bsTooltip("custom_volcano_Width_for_RNAseq", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("custom_volcano_Height_for_RNAseq", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("P_adj_for_RNAseq", "Whether or not to use adjusted p values", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Same_width_for_RNAseq", "Whether the x axis to have the same width from 0", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "MA plot",
                              fluidRow(
                                box(uiOutput("selected_proteins_for_RNAseq_MAplot"),
                                  uiOutput("contrast_for_RNAseq_MA"), width = 6),
                                box(checkboxInput("add_rug_for_RNAseq",
                                                   "Add rug",
                                                   value = FALSE),
                                    checkboxInput("selected_genes_for_RNAseq_MAplot",
                                                   "Selected genes",
                                                   value = FALSE),
                                    checkboxInput("labeled_for_RNAseq_MAplot",
                                                   "Labeled",
                                                   value = FALSE),
                                width = 6)
                                ),
                               fluidRow(
                                box(numericInput("MAplot_Width_for_RNAsqe",
                                                  "width",
                                                  min = 1, max = 30, value = 9), width = 6),
                                box(numericInput("MAplot_Height_for_RNAseq",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                                ),
                               fluidRow(
                                  plotOutput("plot_MA_for_RNAseq"),
                                 downloadButton('downloadMAplot_for_RNAseq', 'Save')),
          shinyBS::bsTooltip("selected_proteins_for_RNAseq_MAplot", "Choose the point labels to show, act when [Labeled] is selected", "top", options = list(container = "body")),
          shinyBS::bsTooltip("contrast_for_RNAseq_MA", "Choose the contrast that you want to show", "top", options = list(container = "body")),
          shinyBS::bsTooltip("add_rug_for_RNAseq", "whether add rug plots in the margins", "top", options = list(container = "body")),
          shinyBS::bsTooltip("selected_genes_for_RNAseq_MAplot", "whether show genes, so that you can choose genes and add them on the plotting by setting [Labeled] TRUE", "top", options = list(container = "body")),
          shinyBS::bsTooltip("labeled_for_RNAseq_MAplot", "whether add labels on the plotting, act when [Selected genes] is TRUE", "top", options = list(container = "body")),
          shinyBS::bsTooltip("MAplot_Width_for_RNAsqe", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("MAplot_Height_for_RNAseq", "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Selected Gene",
                               fluidRow(
                                 box(
                                  uiOutput("selected_genes_for_RNAseq"),
                                  uiOutput("selected_group_for_RNAseq"), width = 8),
                                 box(checkboxInput("ylimZero_for_RNAseq",
                                                   "ylimZero",
                                                   value = TRUE),
                                     checkboxInput("add_labels_for_RNAseq",
                                                   "Add labels",
                                                   value = FALSE),
                                     checkboxInput("labels_repel_for_RNAseq",
                                                   "Labels repel",
                                                   value = TRUE), width = 4)),
                               fluidRow(
                                box(numericInput("selected_genes_Width_for_RNAsqe",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                                box(numericInput("selected_genes_Height_for_RNAseq",
                                                  "height",
                                                  min = 1, max = 30, value = 6), width = 6)
                                ),
                               fluidRow(
                                 plotOutput("plot_selected_genes_for_RNAseq", height = 600),
                                 downloadButton('download_selected_genes_for_RNAseq', 'Save')
                               ),
          shinyBS::bsTooltip("selected_genes_for_RNAseq", "Specify the name of the gene to plot", "top", options = list(container = "body")),
          shinyBS::bsTooltip("selected_group_for_RNAseq", "Choose interesting groups", "top", options = list(container = "body")),
          shinyBS::bsTooltip("ylimZero_for_RNAseq", "Whether set the lower y-axis limit to 0", "top", options = list(container = "body")),
          shinyBS::bsTooltip("add_labels_for_RNAseq", "Whether add labels", "top", options = list(container = "body")),
          shinyBS::bsTooltip("labels_repel_for_RNAseq", "Whether use functions of R package ggrepel to place labels", "top", options = list(container = "body")),
          shinyBS::bsTooltip("selected_genes_Width_for_RNAsqe", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("selected_genes_Height_for_RNAseq", "Height of the figure to export", "top", options = list(container = "body"))
                      )

               )
               )
        )
      ),
      tabItem(#** PR-Heatmap tabItem----
        tabName = "PRHeatmap_tab",
        fluidRow(
          column(width = 7,
            box(
              fluidRow(
                box(numericInput("row_font_size_for_pro_rna", "row font size", min = 0, max = 100, value = 6), width = 4),
                box(numericInput("col_font_size_for_pro_rna", "col font size", min = 0, max = 100, value = 8), width = 4),
                box(numericInput("size_for_pro_rna", "Heat size(4-30)", min = 4, max = 30, value = 6), width = 4)
                ),
              fluidRow(
                box(numericInput("column_title_size_for_pro_rna", "title size", min = 0, max = 100, value = 10), width = 4),
                box(numericInput("Numcol_pro_for_pro_rna", "color for protein", min = 0, max = 100, value = 2), width = 4),
                box(numericInput("Numcol_rna_for_pro_rna", "color for rna", min = 0, max = 100, value = 2), width = 4)
                ),
              # fluidRow(
              #   box(checkboxInput("cluster_columns_for_pro_rna", "Cluster columns", value = TRUE), width = 4),
              #   box(numericInput("heatmap_Width_for_pro_rna", "width", min = 1, max = 30, value = 7), width = 4),
              #   box(numericInput("heatmap_Height_for_pro_rna", "height", min = 1, max = 30, value = 7), width = 4)
              #   ),
              fluidRow(
                box(title = "Table",
                   DT::dataTableOutput("table_for_pro_rna"), width = 12)
                ), width = 12
              ),
               # box(title = "Table",
               #     DT::dataTableOutput("table_for_pro_rna"), width = 12),
          shinyBS::bsTooltip("table_for_pro_rna", "The table that used to plot the heatmap", "top", options = list(container = "body"))
        ),
          column(width = 5,
               box(title = "Heatmap",
                # fluidRow(
                # box(radioButtons(inputId = "Type_for_PRheatmap", label = "Type", choices = c("all", "both", "at least one"), selected = "all", inline = TRUE), width = 12)
                # ),
                               # fluidRow(
                               #  box(numericInput("row_font_size_for_pro_rna",
                               #                    "row font size",
                               #                    min = 0, max = 100, value = 6),
                               #       width = 4),
                               #  box(numericInput("col_font_size_for_pro_rna",
                               #                    "col font size",
                               #                    min = 0, max = 100, value = 8),
                               #       width = 4),
                               #  box(numericInput("size_for_pro_rna",
                               #                    "Heat size(4-30)",
                               #                    min = 4, max = 30, value = 6),
                               #       width = 4)

                               #  ),
                               # fluidRow(
                               #  box(numericInput("column_title_size_for_pro_rna",
                               #                    "title size",
                               #                    min = 0, max = 100, value = 10),
                               #       width = 4),
                               #  box(numericInput("Numcol_pro_for_pro_rna",
                               #                    "color for protein",
                               #                    min = 0, max = 100, value = 2),
                               #       width = 4),
                               #  box(numericInput("Numcol_rna_for_pro_rna",
                               #                    "color for rna",
                               #                    min = 0, max = 100, value = 2),
                               #       width = 4)
                               #  ),
                               fluidRow(
                                box(checkboxInput("cluster_columns_for_pro_rna",
                                                   "Cluster columns",
                                                   value = TRUE), width = 4),

                                box(numericInput("heatmap_Width_for_pro_rna",
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 4),
                                box(numericInput("heatmap_Height_for_pro_rna",
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 4)
                                ),
                               fluidRow(
                                column(width = 6,
                                    box(uiOutput("Type_for_PRheatmap"), width = 12)
                                    ),
                                column(width = 6,
                                    box(selectizeInput("colorbar_for_PRheatmap",
                                       "colorbar",
                                       choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                                       selected = c("RdBu"), multiple = FALSE), width = 12)
                                    )

                                 ),
                               # fluidRow(
                               #   box(uiOutput("Type_for_PRheatmap"), width = 12)
                               #   ),
                               fluidRow(
                                  uiOutput("plot_heatmap_for_pro_rna"),
                                 downloadButton('downloadHeatmap_for_pro_rna', 'Save heatmap'))
                   , width = 12)
        )
          ),
          shinyBS::bsTooltip("colorbar_for_PRheatmap", "Choose the color palette of the colorbar", "top", options = list(container = "body")),
          shinyBS::bsTooltip("row_font_size_for_pro_rna", "Set the size of row labels", "top", options = list(container = "body")),
          shinyBS::bsTooltip("col_font_size_for_pro_rna", "Set the size of column labels", "top", options = list(container = "body")),
          shinyBS::bsTooltip("size_for_pro_rna", "Set the height of the plot region", "top", options = list(container = "body")),
          shinyBS::bsTooltip("column_title_size_for_pro_rna", "Set the size of the title", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Numcol_pro_for_pro_rna", "Set the outer limits of the color scale of protein data", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Numcol_rna_for_pro_rna", "Set the outer limits of the color scale of RNAseq data", "top", options = list(container = "body")),
          shinyBS::bsTooltip("cluster_columns_for_pro_rna", "Whether make cluster on columns", "top", options = list(container = "body")),
          shinyBS::bsTooltip("heatmap_Width_for_pro_rna", "Width of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("heatmap_Height_for_pro_rna", "Height of the figure to export", "top", options = list(container = "body")),
          shinyBS::bsTooltip("Type_for_PRheatmap", "When [import from DEP-LFQ and DEG-RNAseq] is FALSE: all represents all genes that you uploaded; both represents genes that both have RNAseq data and protein data of your uploaded; at least one represents genes that at least have one of RNAseq and protein data of your uploaded. When [import from DEP-LFQ and DEG-RNAseq] is TRUE: based on the union of your selected RNAseq and protein significantly differential genes, at least one represents genes that at least have one of RNAseq and protein data; both represents genes that both have RNAseq data and protein data; both significant represents that both RNAseq data and protein data are significant", "left", options = list(container = "body")),

        ),
      tabItem(#** Extralinks tabItem----
        tabName = "Extralinks_tab",
            fluidRow(
              column(
                width = 12,
                includeMarkdown(system.file("extdata", "extralinks.md", package = "DEP2"))

              )
            )

        )
)

        )

    )
  )


server <- shinyServer(function(input, output) {
  options(shiny.maxRequestSize=5000*1024^2)

  ### UI functions ### --------------------------------------------------------
  output$name <- renderUI({
    selectizeInput("name",
                            "Name column",
                            choices=colnames(data()),
                            selected = "Gene.names")
    })

  output$id <- renderUI({
    selectizeInput("id",
                            "ID column",
                            choices=colnames(data()),
                            selected = "Protein.IDs")
    })

  output$filt <- renderUI({
    selectizeInput("filt",
                            "Filter on columns" ,
                            colnames(data()),
                            multiple = TRUE,
                            selected = c("Reverse","Potential.contaminant"))
    })

  output$control <- renderUI({
    validate(need(!is.null(input$file1), ""))
    if (input$anno == "columns" & !is.null(my_data) & input$contrasts == "control") {
               cols <<- grep("^LFQ", colnames(my_data))
               prefix <<- get_prefix(my_data[,cols] %>% colnames())
               selectizeInput("control", "Control",
                              choices=make.names(unlist(lapply(colnames(my_data)[cols] %>% gsub(prefix,"",.) %>% strsplit(., split = "_"), function(x){x[1]}))), selected = NULL)
             } else {
              if (input$anno == "expdesign" & !is.null(expdesign()) & input$contrasts == "control") {
               selectizeInput("control",
                              "Control",
                              choices = make.names(expdesign()$condition))
             }
             }
    # if (input$anno == "expdesign" & !is.null(expdesign()) & input$contrasts == "control") {
    #            selectizeInput("control",
    #                           "Control",
    #                           choices = make.names(expdesign()$condition))
    #          }
    })

  output$test_manual <- renderUI({
    validate(need(!is.null(input$file1), ""))
    if(!is.null(my_data) & input$contrasts == "manual"){
              cols <<- grep("^LFQ", colnames(my_data))
              prefix <<- get_prefix(my_data[,cols] %>% colnames())
              test_manual_name <- unique(make.names(unlist(lapply(colnames(my_data)[cols] %>% gsub(prefix,"",.) %>% strsplit(., split = "_"), function(x){x[1]}))))
              test_manual_name <- cbind(combn(test_manual_name,2),combn(test_manual_name,2, FUN = rev))
              test_manual_name <- apply(test_manual_name, 2, function(i){paste(i[1], i[2], sep = "_vs_")})
              selectizeInput("test_manual", "Manual test",
                              choices=test_manual_name, selected = NULL, multiple = TRUE)
             }
    })

  output$Peptides1 <- renderUI({
    if(input$if_peptide_color) {
     colourInput(inputId = "Peptides_1",
                 label = "Peptides1",
                 showColour = "both",
                 palette = "square",
                 value = "#A020F0")
    }
    })

  output$Peptides2 <- renderUI({
    if(input$if_peptide_color) {
    colourpicker::colourInput(inputId = "Peptides_2",
                label = "Peptides2",
                showColour = "both",
                palette = "square",
                value = "#0000FF")#
    }
    })
  ### Reactive functions ### --------------------------------------------------
  expdesign <- reactive({
    inFile <- input$file2
    if (is.null(inFile))
      return(NULL)
     read.csv(inFile$datapath, header = TRUE,
             sep = "\t", stringsAsFactors = FALSE) %>%
      mutate(id = row_number())
  })

  data <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    my_data <<- read.csv(inFile$datapath, header = TRUE,
             sep = "\t", stringsAsFactors = FALSE) %>%
      mutate(id = row_number())
  })

  filt <- reactive({
    data <<- data()
    cols <- grep("^LFQ", colnames(data))

    filtered <- DEP2:::filter_MaxQuant(data, input$filt)
    unique_names <- make_unique(filtered, input$name, input$id)

    if (input$anno == "columns") {
      se <- make_se_parse(unique_names, cols, mode = "delim", sep = "_")
    }
    if (input$anno == "expdesign") {
      se <- make_se(unique_names, cols, expdesign())
    }
     my_filt <<- filter_missval(se, thr = input$thr)
  })

  norm <- reactive({
     my_norm <<- normalize_vsn(filt())
  })


  imp <-  reactive({
     inFile1 <- input$resultRData

     if(is.null(inFile1)){
      if(input$imputation == "mixed on proteins") {
        # Extract protein names with missing values
        # in all replicates of at least one condition
        proteins_MNAR <- get_df_long(norm()) %>%
          dplyr::group_by(name, condition) %>%
          dplyr::summarize(NAs = all(is.na(intensity))) %>%
          dplyr::filter(NAs) %>%
          dplyr::pull(name) %>%
          unique()

        # Get a logical vector
        MNAR <- names(norm()) %in% proteins_MNAR

        # Perform a mixed imputation
        my_imp <<- impute(
          norm(),
          fun = "mixed",
          randna = !MNAR, # we have to define MAR which is the opposite of MNAR
          mar = "knn", # imputation function for MAR
          mnar = "MinDet") # imputation function for MNAR
        validate(
      need(length(which(is.na((assay(my_imp)[which(MNAR),])))) == 0, "Please select another imputation method")
    )
        my_imp

        } else {
          if(input$imputation == "mixed on samples") {
            sample_specific_imputation <- norm()
            MSnSet <- as(sample_specific_imputation, "MSnSet")

            # Impute differently for two sets of samples
            MSnSet_imputed1 <- MSnbase::impute(MSnSet[, which(gsub("_\\d$", "", colnames(assay(sample_specific_imputation))) %in% input$control)], method = "MinProb") # control samples that need not random (MNAR) left-censored imputation
            MSnSet_imputed2 <- MSnbase::impute(MSnSet[, -which(gsub("_\\d$", "", colnames(assay(sample_specific_imputation))) %in% input$control)], method = "knn")# treatment samples that need random(MAR) imputation

            # Combine into the SummarizedExperiment object
            assay(sample_specific_imputation) <- cbind(
                MSnbase::exprs(MSnSet_imputed1),
                MSnbase::exprs(MSnSet_imputed2))
            my_imp <<- sample_specific_imputation
            validate(
      need(length(which(is.na((assay(my_imp))))) == 0, "Please select another imputation method")
    )
            my_imp
            } else {
              #set.seed(12345)
              my_imp <<- DEP2::impute(norm(), input$imputation)
            }
        }
    } else {
     return(NULL)
    }
    })


  df <- reactive({
    if(input$contrasts == "control") {
      validate(
      need(input$control != "", "Please select a control condition under menuItem Columns in the DEP-LFQ options of the sidebar")
    )
    }

    if(input$contrasts == "manual") {
      validate(
      need(input$test_manual != "", "Please select manual contrasts to test under menuItem Columns")
    )
    }

     inFile1 <- input$resultRData
     if(is.null(inFile1)){
      if(input$contrasts == "control"){
        if(input$fdr_correction == "fdrtool") {
          df <- test_diff(se = imp(), type = input$contrasts, control = input$control)
          } else {
            df <- test_limma(se = imp(), type = input$contrasts, control = input$control)
          }
      }

      if(input$contrasts == "all") {
        if(input$fdr_correction == "fdrtool") {
          df <- test_diff(se = imp(), type = input$contrasts)
          } else {
            df <- test_limma(se = imp(), type = input$contrasts)
          }
      }

      if(input$contrasts == "manual") {
        if(input$fdr_correction == "fdrtool") {
          df <- test_diff(se = imp(), type = input$contrasts, test = input$test_manual)
          } else {
            df <- test_limma(se = imp(), type = input$contrasts, test = input$test_manual)
          }
      }
      } else {
      #load your saved RData in order to get the same result (imp is the key of if result from two analysis being the same)
      load(file = inFile1$datapath)
      my_imp <<- my_imp
      if(input$contrasts == "control"){
        if(input$fdr_correction == "fdrtool") {
          df <- test_diff(se = my_imp, type = input$contrasts, control = input$control)
          } else {
            df <- test_limma(se = my_imp, type = input$contrasts, control = input$control)
          }
      }

      if(input$contrasts == "all") {
        if(input$fdr_correction == "fdrtool") {
          df <- test_diff(se = my_imp, type = input$contrasts)
          } else {
            df <- test_limma(se = my_imp, type = input$contrasts)
          }
      }

      if(input$contrasts == "manual") {
        if(input$fdr_correction == "fdrtool") {
          df <- test_diff(se = my_imp, type = input$contrasts, test = input$test_manual)
          } else {
            df <- test_limma(se = my_imp, type = input$contrasts, test = input$test_manual)
          }
      }
      }
      df
    })

  dep <- reactive({
    my_dep <<- add_rejections(df(), input$p, input$lfc)
  })


  output$Save_RData <- downloadHandler(
    filename = function() { paste("results", ".RData", sep = "") },
    content = function(file) {
    	withProgress(message = 'Please wait ...', value = 0.66, {
      save(my_data, my_filt, my_norm, my_imp, my_dep, file=file)})}
  )

  ## All object and functions upon 'Analyze' input  ### ----------------------

  observeEvent(input$analyze, {

    ### Interactive UI functions ### ------------------------------------------
    output$downloadTable <- renderUI({
      selectizeInput("dataset",
                     "Choose a dataset to save" ,
                     c("results","significant_proteins",
                       "displayed_subset","full_dataset"))
    })

    output$downloadButton <- renderUI({
      downloadButton('downloadData', 'Save', class = "downloadData")
    })

    output$downloadButton_for_save_RData <- renderUI({
      downloadButton("Save_RData", "save result RData", class = "Save_RData")
      })

    output$significantBox <- renderInfoBox({
      num_total <- dep() %>%
        nrow()
      num_signif <- dep() %>%
        .[rowData(.)$significant, ] %>%
        nrow()
      frac <- num_signif / num_total

      if(frac > 0.2) {
        info_box <- infoBox("Significant proteins",
                            paste0(num_signif,
                                   " out of ",
                                   num_total),
                            paste0("Too large fraction (",
                                   signif(frac * 100, digits = 3),
                                   "%) of proteins differentially expressed"),
                            icon = icon("minus", lib = "glyphicon"),
                            color = "orange",
                            width = 4)
      }
      if(frac == 0) {
        info_box <- infoBox("Significant proteins",
                            paste0(num_signif,
                                   " out of ",
                                   num_total),
                            "No proteins differentially expressed",
                            icon = icon("thumbs-down", lib = "glyphicon"),
                            color = "red",
                            width = 4)
      }
      if(frac > 0 & frac <= 0.2) {
        info_box <-     infoBox("Significant proteins",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              paste0(signif(frac * 100, digits = 3),
                                     "% of proteins differentially expressed"),
                              icon = icon("thumbs-up", lib = "glyphicon"),
                              color = "green",
                              width = 4)
      }
      info_box
    })

    output$select <- renderUI({
      row_data <- rowData(dep())
      cols <- grep("_significant", colnames(row_data))
      names <- colnames(row_data)[cols]
      names <- gsub("_significant", "", names)
      selectizeInput("select",
                     "Select direct comparisons",
                     choices=names,
                     multiple = TRUE)
    })

    output$exclude <- renderUI({
      row_data <- rowData(dep())
      cols <- grep("_significant", colnames(row_data))
      names <- colnames(row_data)[cols]
      names <- gsub("_significant","",names)
      selectizeInput("exclude",
                     "Exclude direct comparisons",
                     choices = names,
                     multiple = TRUE)
    })


    output$heatmap_cntrst <- renderUI({
      if (!is.null(selected())) {
        df <- rowData(selected())
        cols <- grep("_significant$",colnames(df))
        selectizeInput("heatmap_cntrst",
                       "Contrast",
                       choices = gsub("_significant", "", colnames(df)[cols]), multiple = TRUE)
      }
    })


    output$mysplit <- renderUI({
        selectizeInput("mysplit",
                       "my split",
                       choices = c(1 : input$k), multiple = TRUE)

    })

    output$Custom_columns_order <- renderUI({
        selectizeInput("Custom_columns_order",
                       "Custom columns order",
                       choices = colnames(assay(dep())), multiple = TRUE)

    })

    output$Custom_columns_order_for_missval_heatmap <- renderUI({
      selectizeInput("Custom_columns_order_for_missval_heatmap",
                       "Custom columns order",
                       choices = colnames(assay(dep())), multiple = TRUE)
      })


    output$volcano_cntrst <- renderUI({
      if (!is.null(selected())) {
        df <- rowData(selected())
        cols <- grep("_significant$",colnames(df))
        selectizeInput("volcano_cntrst",
                       "Contrast",
                       choices = gsub("_significant", "", colnames(df)[cols]))
      }
    })

    output$Volcano_cntrst <- renderUI({
      if (!is.null(selected())) {
        df <- rowData(selected())
        cols <- grep("_significant$",colnames(df))
        selectizeInput("Volcano_cntrst",
                       "Contrast",
                       choices = gsub("_significant", "", colnames(df)[cols]))
      }
    })


    output$selected_proteins <- renderUI({
        row_Data <- rowData(dep())
        selectizeInput("selected_proteins",
                       "selected proteins",
                       choices = row_Data$name, selected = NULL, multiple = TRUE)
    })


    ### Reactive functions ### ------------------------------------------------
    excluded <- reactive({
      DEP2:::exclude_deps(dep(), input$exclude)
    })

    selected <- reactive({
      DEP2:::select_deps(excluded(), input$select)
    })

    res <- reactive({
      get_results(selected())
    })

    table <- reactive({
      DEP2:::get_table(res(), input$pres)
    })

    selected_plot_input <- reactive ({
      if(!is.null(input$table_rows_selected)) {
        selected_id <- table()[input$table_rows_selected,1]
        plot_single(selected(), selected_id, input$pres)
      }
    })

    heatmap_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        plot_heatmap(dep = selected(),
                     type = input$pres,
                     manual = input$manual_heatmap,
                     manual_name = input$heatmap_cntrst,
                     kmeans = TRUE,
                     k = input$k,
                     color = input$colorbar,
                     col_limit = input$limit,
                     #same_trend = input$same_trend,
                     row_font_size = input$row_font_size,
                     col_font_size = input$col_font_size,
                     cluster_columns = input$cluster_columns,
                     if_mysplit = input$if_mysplit,
                     mysplit = input$mysplit,
                     if_rowname_color = input$if_rowname_color,
                     column_order = if(input$cluster_columns) {NULL} else {input$Custom_columns_order})
      })
    })

    volcano_input <- reactive({
      if(!is.null(input$volcano_cntrst)) {
        plot_volcano(selected(),
                     input$volcano_cntrst,
                     input$fontsize,
                     input$check_names,
                     input$p_adj,
                     same_width = input$same_width,
                     my_breaks = input$my_breaks,
                     mybreaks = as.numeric(input$mybreaks))
      }
    })
    custom_volcano_input <- reactive({
      if(!is.null(input$Volcano_cntrst)) {
      withProgress(message = 'Plotting', value = 0.66, {
         my_plot_volcano(dep = selected(),
                     contrast = input$Volcano_cntrst,
                     adjusted = input$P_adj,
                     labelWay = input$labelWay,
                     showNum = input$showNum,
                     chooseTolabel = input$selected_proteins,
                     fontSize = input$fontSize,
                     dotsize = input$dotsize,
                     same_width = input$Same_width,
                     fcCutoff = input$lfc,
                     adjpCutoff = input$p,
                     if_peptide_color = input$if_peptide_color,
                     label.rectangle = input$if_label_rectangle,
                     stroke = input$stroke,
                     down_color = input$down_color,
                     stable_color = input$stable_color,
                     up_color = input$up_color,
                     peptide_1_color = input$Peptides_1,
                     peptide_2_color = input$Peptides_2
                     )
        })
      }
    })

    pca_input <- reactive({
      plot_pca(dep(), indicate = input$Indicate, if_square = input$if_square)
    })

    Pearson_correlation_input <- reactive({
      plot_cor(dep(), pal = input$Pearson_pal, pal_rev = input$Pearson_pal_rev, lower = input$Pearson_lower, upper = input$Pearson_upper, add_values = input$add_values_for_DEP_person, value_size = input$value_size_for_DEP_person, digits = input$value_digits_for_DEP_person)
    })

    Gowers_distance_input <- reactive({
      plot_dist(dep(), pal = input$Gower_pal, pal_rev = input$Gower_pal_rev, add_values = input$add_values_for_DEP_gower, value_size = input$value_size_for_DEP_gower, digits = input$value_digits_for_DEP_gower)
    })

    Sample_CVs_input <- reactive({
      plot_cvs(dep())
    })


    norm_input <- reactive({
      plot_normalization(filt(),
                         norm())
    })

    missval_input <- reactive({
      plot_missval(norm(),
                  cluster_columns = input$cluster_columns_for_missval_heatmap,
                  column_order = if(input$cluster_columns_for_missval_heatmap) {NULL} else {input$Custom_columns_order_for_missval_heatmap}
        )
    })

    detect_input <- reactive({
      plot_detect(norm())
    })

    imputation_input <- reactive({
      plot_imputation(norm(),
                      df())
    })

    numbers_input <- reactive({
      plot_numbers(norm())
    })

    coverage_input <- reactive({
      plot_coverage(norm())
    })

    ### Output functions ### --------------------------------------------------
    output$table <- DT::renderDataTable({
      table()
    }, options = list(pageLength = 25, scrollX = T),
    selection = list(selected = c(1)))

    output$selected_plot <- renderPlot({
      selected_plot_input()
    })

    output$heatmap <- renderPlot({
      heatmap_input()
    })

    output$volcano <- renderPlot({
      volcano_input()
    })

    output$custom_volcano <- renderPlot({
      custom_volcano_input()
    })

    output$pca <- renderPlot({
      pca_input()
    })

    output$Pearson_correlation <- renderPlot({
      Pearson_correlation_input()
    })

    output$Gowers_distance <- renderPlot({
      Gowers_distance_input()
    })

     output$Sample_CVs <- renderPlot({
      Sample_CVs_input()
    })

    output$norm <- renderPlot({
      norm_input()
    })

    output$missval <- renderPlot({
      missval_input()
    })

    output$detect <- renderPlot({
      detect_input()
    })

    output$imputation <- renderPlot({
      imputation_input()
    })

    output$numbers <- renderPlot({
      numbers_input()
    })

    output$coverage <- renderPlot({
      coverage_input()
    })

    observe({
      output$plot <- renderUI({
        plotOutput("heatmap", height = (100 * as.numeric(input$size)))
      })
    })

    ### Download objects and functions ### ------------------------------------
    datasetInput <- reactive({
      switch(input$dataset,
             "results" = get_results(dep()),
             "significant_proteins" = get_results(
              dep()) %>%
               dplyr::filter(significant) %>%
               dplyr::select(-significant),
             "displayed_subset" = res() %>%
               dplyr::filter(significant) %>%
               dplyr::select(-significant),
             "full_dataset" = get_df_wide(dep()))
    })

    output$downloadData <- downloadHandler(
      filename = function() { paste(input$dataset, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput(),
        file,
        col.names = TRUE,
        row.names = FALSE,
        sep ="\t") }
    )

    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0("Barplot_", table()[input$table_rows_selected,1], ".pdf")
      },
      content = function(file) {
        pdf(file, width = input$selected_plot_Width, height = input$selected_plot_Height)
        print(selected_plot_input())
        dev.off()
      }
    )

    output$downloadHeatmap <- downloadHandler(
      filename = 'Heatmap.pdf',
      content = function(file) {
        pdf(file, width = input$Width, height = input$Height)
        print(heatmap_input())
        dev.off()
      }
    )

    output$downloadVolcano <- downloadHandler(
      filename = function() {
        paste0("Volcano_", input$volcano_cntrst, ".pdf")
      },
      content = function(file) {
        pdf(file, width = input$Volcano_Width, height = input$Volcano_Height)
        print(volcano_input())
        dev.off()
      }
    )

    output$download_custom_volcano <- downloadHandler(
      filename = function() {
        paste0("custom_Volcano_", input$volcano_cntrst, ".pdf")
      },
      content = function(file) {
        pdf(file, width = input$custom_volcano_Width, height = input$custom_volcano_Height)
        print(custom_volcano_input())
        dev.off()
      }
    )

    output$downloadPca <- downloadHandler(
      filename = 'Pca.pdf',
      content = function(file) {
        pdf(file, width = input$pca_Width, height = input$pca_Height)
        print(pca_input())
        dev.off()
      }
    )

    output$download_Pearson_correlation <- downloadHandler(
      filename = 'Pearson_correlation.pdf',
      content = function(file) {
        pdf(file, width = input$Pearson_Width, height = input$Pearson_Height)
        print(Pearson_correlation_input())
        dev.off()
      }
    )

    output$download_Gowers_distance <- downloadHandler(
      filename = 'Gowers_distance.pdf',
      content = function(file) {
        pdf(file, width = input$Gower_Width, height = input$Gower_Height)
        print(Gowers_distance_input())
        dev.off()
      }
    )

    output$download_Sample_CVs <- downloadHandler(
      filename = 'Sample_CVs.pdf',
      content = function(file) {
        pdf(file, width = input$Sample_CVs_Width, height = input$Sample_CVs_Height)
        print(Sample_CVs_input())
        dev.off()
      }
    )

    output$downloadNorm <- downloadHandler(
      filename = "normalization.pdf",
      content = function(file) {
        pdf(file, width = input$norm_Width, height = input$norm_Height)
        print(norm_input())
        dev.off()
      }
    )

    output$downloadMissval <- downloadHandler(
      filename = "missing_values_heatmap.pdf",
      content = function(file) {
        pdf(file, width = input$missval_heatmap_Width, height = input$missval_heatmap_Height)
        print(missval_input())
        dev.off()
      }
    )

    output$downloadDetect <- downloadHandler(
      filename = "missing_values_quant.pdf",
      content = function(file) {
        pdf(file, width = input$detect_Width, height = input$detect_Height)
        gridExtra::grid.arrange(detect_input())
        dev.off()
      }
    )

    output$downloadImputation <- downloadHandler(
      filename = "imputation.pdf",
      content = function(file) {
        pdf(file, width = input$imputation_Width, height = input$imputation_Height)
        print(imputation_input())
        dev.off()
      }
    )

    output$downloadNumbers <- downloadHandler(
      filename = "numbers.pdf",
      content = function(file) {
        pdf(file, width = input$numbers_Width, height = input$numbers_Height)
        print(numbers_input())
        dev.off()
      }
    )

    output$downloadCoverage <- downloadHandler(
      filename = "coverage.pdf",
      content = function(file) {
        pdf(file, width = input$coverage_Width, height = input$coverage_Height)
        print(coverage_input())
        dev.off()
      }

    )
  })

## for annotation
output$import_for_annotation <- renderUI({
    if(input$import_from_for_annotation){
                 radioButtons("import_for_annotation",
                              "Choose genes",
                              choices = c("UPregu for DEP-LFQ", "DOWNregu for DEP-LFQ", "UPDOWN for DEP-LFQ", "UPregu for DEG-RNAseq", "DOWNregu for DEG-RNAseq", "UPDOWN for DEG-RNAseq"),
                              selected = "UPregu for DEP-LFQ")
    }
    })

output$import_contrast_for_annotation <- renderUI({
  if(input$import_from_for_annotation) {
if(input$import_for_annotation == "UPregu for DEP-LFQ" | input$import_for_annotation == "DOWNregu for DEP-LFQ" | input$import_for_annotation == "UPDOWN for DEP-LFQ") {
  validate(need(!is.null(data()), "Please go to the DEP-LFQ options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_annotation", "Contrast", choices = c(get_results(dep()) %>%
  dplyr::filter(significant) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique(), "Any significant"), selected = NULL)
    } else {
      validate(need(!is.null(countData()), "Please go to the DEG-RNAseq options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_annotation", "Contrast", choices = c(as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique(), "Any significant"), selected = NULL)
    }
  } else {
    return(NULL)
  }

  })

output$text_input_for_annotation <- renderUI({
  if(!input$import_from_for_annotation){
    textAreaInput(inputId = "text_input_for_annotation", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%")
  }
  })

output$organism_for_annotation <- renderUI({
  # if(!input$import_from_for_annotation){
    selectizeInput("organism_for_annotation", "Select organism", choices=c("human","mouse"), selected="human")
  # }
  })

######### observeEvent for gene annotation input
observeEvent(input$analyze_for_annotation,{
    withProgress(message = 'Please wait', value = 0.66, {
      output$downloadannotation_for_output <- renderUI({
          #align left;width: 10em; vertical-align: center
          downloadButton("downloadannotation", "Save table", class = "downloadannotation")#btn btn-success
        })

    if(!input$import_from_for_annotation) {
    # genelist <- reactive({ strsplit(input$text_input_for_annotation,'\n')[[1]] })
    # gene_name <-reactive(unlist(strsplit(genelist(),";")[]))
    # gene_df <<- data.frame(name=gene_name())
    genelist <- strsplit(input$text_input_for_annotation,'\n')[[1]]
    gene_name <-unique(unlist(strsplit(genelist,";")[]))
    gene_df <<- data.frame(name=gene_name)
    }

    if(input$import_from_for_annotation) {
      if(input$import_for_annotation == "UPregu for DEP-LFQ" | input$import_for_annotation == "DOWNregu for DEP-LFQ" | input$import_for_annotation == "UPDOWN for DEP-LFQ") {
        if(input$import_contrast_for_annotation == "Any significant") {
          if(input$import_for_annotation == "UPDOWN for DEP-LFQ") {
            gene_df <<- get_results(dep()) %>% dplyr::filter(significant) %>% dplyr::select(name)
            } else {
              if(input$import_for_annotation == "UPregu for DEP-LFQ") {
                index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
                df <- get_results(dep()) %>% dplyr::filter(significant)
                cols_diff = df %>% dplyr::select(index)
                cols_diff_reject = cols_diff > 0
                cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                index1 = apply(cols_diff_reject, 1, all)
                gene_df <<- df[index1, ] %>% dplyr::select(name)
                } else {
                  if(input$import_for_annotation == "DOWNregu for DEP-LFQ") {
                    index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
                    df <- get_results(dep()) %>% dplyr::filter(significant)
                    cols_diff = df %>% dplyr::select(index)
                    cols_diff_reject = cols_diff < 0
                    cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                    index1 = apply(cols_diff_reject, 1, all)
                    gene_df <<- df[index1, ] %>% dplyr::select(name)
                  }
                }
            }
          } else {
            if(input$import_for_annotation == "UPDOWN for DEP-LFQ") {
              gene_df <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_annotation, "significant", sep = "_"))) %>% dplyr::select(name)
              } else {
                if(input$import_for_annotation == "UPregu for DEP-LFQ") {
                  gene_df <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_annotation, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_annotation, "ratio", sep = "_")) > 0) %>% dplyr::select(name)
                  } else {
                    if(input$import_for_annotation == "DOWNregu for DEP-LFQ") {
                      gene_df <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_annotation, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_annotation, "ratio", sep = "_")) < 0) %>% dplyr::select(name)
                    }
                  }
              }
          }
        } else {
          if(input$import_for_annotation == "UPregu for DEG-RNAseq" | input$import_for_annotation == "DOWNregu for DEG-RNAseq" | input$import_for_annotation == "UPDOWN for DEG-RNAseq") {
            if(input$import_contrast_for_annotation == "Any significant") {
              if(input$import_for_annotation == "UPDOWN for DEG-RNAseq") {
                gene_df <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant) %>% dplyr::select(name)
                } else {
                  if(input$import_for_annotation == "UPregu for DEG-RNAseq") {
                    index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
                    df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
                    cols_diff = df %>% dplyr::select(index)
                    cols_diff_reject = cols_diff > 0
                    cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                    index1 = apply(cols_diff_reject, 1, all)
                    gene_df <<- df[index1, ] %>% dplyr::select(name)
                    } else {
                      if(input$import_for_annotation == "DOWNregu for DEG-RNAseq") {
                        index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
                        df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
                        cols_diff = df %>% dplyr::select(index)
                        cols_diff_reject = cols_diff < 0
                        cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                        index1 = apply(cols_diff_reject, 1, all)
                        gene_df <<- df[index1, ] %>% dplyr::select(name)
                      }
                    }
                }
              } else {
                if(input$import_for_annotation == "UPDOWN for DEG-RNAseq") {
                  gene_df <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_annotation, "significant", sep = "_"))) %>% dplyr::select(name)
                  } else {
                    if(input$import_for_annotation == "UPregu for DEG-RNAseq") {
                      gene_df <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_annotation, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_annotation, "log2FoldChange", sep = "_")) > 0) %>% dplyr::select(name)
                      } else {
                        if(input$import_for_annotation == "DOWNregu for DEG-RNAseq") {
                          gene_df <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_annotation, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_annotation, "log2FoldChange", sep = "_")) < 0) %>% dplyr::select(name)
                        }
                      }
                  }
              }
          }

        }
    }

    gene_df$name = rm_digit_end(gene_df$name)
    gene_df = data.frame(name = gene_df[!duplicated(gene_df$name),], stringsAsFactors = F)

    if(!nrow(gene_df) == 0) {
      gene_df <<- gene_df %>% dplyr::filter(!is.na(name))
    }

    if(nrow(gene_df) == 0) {
          output$annotation_Table <- DT::renderDataTable({
          shiny::validate(need(nrow(gene_df) != 0, message = "No genes meet your requirements, and can not do the annotation analysis"))
          DT::datatable(data.frame(a = "I am for check"), filter = 'top', options = list( autoWidth = F,scrollX = TRUE))
    })
    }
    shiny::validate(need(nrow(gene_df) != 0, message = "No genes meet your requirements, and can not do the annotation analysis"))

    gene_df$name <-  sapply(gene_df$name,the1stname)
    if(input$organism_for_annotation == "human"){
      library(org.Hs.eg.db)
      orgDB <- org.Hs.eg.db
      kegg_organism <- "hsa"
    }else if(input$organism_for_annotation == "mouse"){
      library(org.Mm.eg.db)
      orgDB <- org.Mm.eg.db
      kegg_organism <- "mmu"
    }
    # ids <- gene_df$ENTREZID <-  mapIds(x = orgDB,
    #                              keys = gene_df$name,
    #                              keytype = "SYMBOL",
    #                              column = "ENTREZID")
    check_ids <<- my_to_entrezid(orgDB = orgDB, gene = as.character(gene_df$name))$id
    # shiny::validate(need(class(check_ids) != "logical", message = "No genes can be mapped to ENTREZID, and can not do the  annotation analysis"))
    if(class(check_ids) == "logical") {
          output$annotation_Table <- DT::renderDataTable({
            shiny::validate(need(class(check_ids) != "logical", message = "No genes can be mapped to ENTREZID, and can not do the  annotation analysis"))
            DT::datatable(data.frame(a = "I am for check"), filter = 'top', options = list( autoWidth = F,scrollX = TRUE))
    })
    }
    shiny::validate(need(class(check_ids) != "logical", message = "No genes can be mapped to ENTREZID, and can not do the  annotation analysis"))
    ids <- gene_df$ENTREZID <- check_ids

    gene_df$ENTREZID <- as.character(gene_df$ENTREZID)
    #PG2<-dplyr::select(PG,c(Gene.names,EntrezID))
    ids <- na.omit(unlist(ids))
    ##ids <- ids[1:10]
    ##uds2 <- as.data.frame(ids)
    ids <<- as.vector(ids)
    # if(input$organism_for_annotation == "human"){
    #   library(org.Hs.eg.db)
    #   orgDB <- org.Hs.eg.db
    #   kegg_organism <- "hsa"
    # }else if(input$organism_for_annotation == "mouse"){
    #   library(org.Mm.eg.db)
    #   orgDB <- org.Mm.eg.db
    #   kegg_organism <- "mmu"
    # }
    bgo <<- Geneannotate(ids,genedb=orgDB)
    IDS <- unique(bgo$ENTREZID)
    # print("antating")
    # library(future.apply)
    # plan(multiprocess)
    # system.time(annots <- future_sapply(IDS,mergego,bgo=bgo) )

    library(parallel)
    cpus = detectCores(logical = F)
    cl <- makeCluster(cpus)
    annots <<- parSapply(cl,IDS,FUN  = mergego,bgo=bgo)
    stopCluster(cl)
    #rownames(annots)<-c("Ontology.CellComponent","Definition.CC","Ontology.biologicalProcess","Definition.BP","Ontology.MolecularFunction","Definition.MF")
    # rownames(annots)<-c("Ontology.CellComponent","Ontology.biologicalProcess","Ontology.MolecularFunction")
    # annots2 <- annots %>% t() %>%as.data.frame() %>%mutate(ENTREZID=colnames(annots))
    # mergepg <<- merge(gene_df,annots2,by.x="EntrezID",by.y="ENTREZID")
    rownames(annots)<-c("GO.CC.name","GO.BP.name","GO.MF.name", "PFAM", "Gene.Description", "KEGG.name", "reactome.name")
    annots2 <<- annots %>% t() %>%as.data.frame() %>%mutate(ENTREZID=colnames(annots))
    mergepg <<- merge(gene_df,annots2,by.x="ENTREZID",by.y="ENTREZID")
    mergepg = mergepg[ , c(2, 7, 1, 3:6, 8:9)]


    output$annotation_Table <- DT::renderDataTable(
      mergepg, filter = 'top', options = list( autoWidth = F,scrollX = TRUE
      )
    )

    output$downloadannotation <- downloadHandler(
      filename = function() { paste("gene annotation", ".txt", sep = "") },
      content = function(file) {
        write.table(mergepg,
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") }
    )
      })
  })

## for go analysis
output$import_for_go <- renderUI({
    if(input$import_from_for_go){
                 radioButtons("import_for_go",
                              "Choose genes",
                              choices = c("UPregu for DEP-LFQ", "DOWNregu for DEP-LFQ", "UPDOWN for DEP-LFQ", "UPregu for DEG-RNAseq", "DOWNregu for DEG-RNAseq", "UPDOWN for DEG-RNAseq"),
                              selected = "UPregu for DEP-LFQ")
    }
    })

output$import_contrast_for_go <- renderUI({
  if(input$import_from_for_go) {
if(input$import_for_go == "UPregu for DEP-LFQ" | input$import_for_go == "DOWNregu for DEP-LFQ" | input$import_for_go == "UPDOWN for DEP-LFQ") {
  validate(need(!is.null(data()), "Please go to the DEP-LFQ options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_go", "Contrast", choices = c(get_results(dep()) %>%
  dplyr::filter(significant) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique(), "Any significant"), selected = NULL)
    } else {
      validate(need(!is.null(countData()), "Please go to the DEG-RNAseq options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_go", "Contrast", choices = c(as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique(), "Any significant"), selected = NULL)
    }
  } else {
    return(NULL)
  }

  })

output$text_input_for_go <- renderUI({
  if(!input$import_from_for_go){
    textAreaInput(inputId = "text_input_for_go", label = "Please paste your gene list", placeholder = "TP53\nPTEN\n\nor\n\nTP53  6.21\nPTEN  -1.53", rows = 8, width = "100%")
  }
  })

output$organism_for_go <- renderUI({
  # if(!input$import_from_for_go){
    selectizeInput("organism_for_go", "Select organism", choices=c("human","mouse"), selected="human")
  # }
  })

output$df_with_lg2fc <- renderUI({
  # if(!input$import_from_for_go) {
    checkboxInput("df_with_lg2fc", "If with log2 fold change", value = FALSE)
  # }
  })

  # observeEvent for go analysis input
  observeEvent(input$analyze_for_go,{
      ### Interactive UI functions ### ------------------------------------------
    output$downloadTable_go <- renderUI({
      selectizeInput("dataset_for_go",
                     "Choose a dataset to save" ,
                     c("full_results","significant_results"
                       ))
    })

    output$downloadButton_go <- renderUI({
      downloadButton('downloadgo', 'Save table', class = "downloadgo")
    })

    ## reactive functions
    if(!input$import_from_for_go) {
      genelist <- reactive({ strsplit(input$text_input_for_go,'\n')[[1]] })
      gene_name <-reactive(unlist(strsplit(genelist(),";")))

      gene_df <- reactive({
        gene_df <<- data.frame(name=gene_name())
        gene_df$name <-  sapply(gene_df$name,the1stname)
        gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
        rownames(gene_df) = NULL
        if(!input$df_with_lg2fc){
        colnames(gene_df) = "name"
        gene_df <- as.data.frame(gene_df)
        test <<- gene_df
      } else {
        colnames(gene_df) = c("name","fc")
        gene_df <- as.data.frame(gene_df)
        gene_df$fc = as.numeric(as.character(gene_df$fc))
        gene_df$name = as.character(gene_df$name)
        test <<- gene_df
      }
        gene_df
        })
    }



    if(input$import_from_for_go) {
      if(input$import_for_go == "UPregu for DEP-LFQ" | input$import_for_go == "DOWNregu for DEP-LFQ" | input$import_for_go == "UPDOWN for DEP-LFQ") {
        if(input$import_contrast_for_go == "Any significant") {
          if(input$import_for_go == "UPDOWN for DEP-LFQ") {
            gene_df <<- get_results(dep()) %>% dplyr::filter(significant) %>% dplyr::select(name)
            } else {
              if(input$import_for_go == "UPregu for DEP-LFQ") {
                index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
                df <- get_results(dep()) %>% dplyr::filter(significant)
                cols_diff = df %>% dplyr::select(index)
                cols_diff_reject = cols_diff > 0
                cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                index1 = apply(cols_diff_reject, 1, all)
                gene_df <<- df[index1, ] %>% dplyr::select(name)
                } else {
                  if(input$import_for_go == "DOWNregu for DEP-LFQ") {
                    index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
                    df <- get_results(dep()) %>% dplyr::filter(significant)
                    cols_diff = df %>% dplyr::select(index)
                    cols_diff_reject = cols_diff < 0
                    cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                    index1 = apply(cols_diff_reject, 1, all)
                    gene_df <<- df[index1, ] %>% dplyr::select(name)
                  }
                }
            }
          } else {
            if(input$import_for_go == "UPDOWN for DEP-LFQ") {
              if(!input$df_with_lg2fc) {
                gene_df <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% dplyr::select(name)
                } else {
                  gene_df <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% dplyr::select(name, paste(input$import_contrast_for_go, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_go, "ratio", sep = "_"))
                }
              } else {
                if(input$import_for_go == "UPregu for DEP-LFQ") {
                  if(!input$df_with_lg2fc) {
                    gene_df <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_go, "ratio", sep = "_")) > 0) %>% dplyr::select(name)
                    } else {
                      gene_df <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_go, "ratio", sep = "_")) > 0) %>% dplyr::select(name, paste(input$import_contrast_for_go, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_go, "ratio", sep = "_"))
                    }
                  } else {
                    if(input$import_for_go == "DOWNregu for DEP-LFQ") {
                      if(!input$df_with_lg2fc) {
                        gene_df <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_go, "ratio", sep = "_")) < 0) %>% dplyr::select(name)
                        } else {
                          gene_df <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_go, "ratio", sep = "_")) < 0) %>% dplyr::select(name, paste(input$import_contrast_for_go, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_go, "ratio", sep = "_"))
                        }
                    }
                  }
              }
          }
        } else {
          if(input$import_for_go == "UPregu for DEG-RNAseq" | input$import_for_go == "DOWNregu for DEG-RNAseq" | input$import_for_go == "UPDOWN for DEG-RNAseq") {
            if(input$import_contrast_for_go == "Any significant") {
              if(input$import_for_go == "UPDOWN for DEG-RNAseq") {
                gene_df <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant) %>% dplyr::select(name)
                } else {
                  if(input$import_for_go == "UPregu for DEG-RNAseq") {
                    index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
                    df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
                    cols_diff = df %>% dplyr::select(index)
                    cols_diff_reject = cols_diff > 0
                    cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                    index1 = apply(cols_diff_reject, 1, all)
                    gene_df <<- df[index1, ] %>% dplyr::select(name)
                    } else {
                      if(input$import_for_go == "DOWNregu for DEG-RNAseq") {
                        index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
                        df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
                        cols_diff = df %>% dplyr::select(index)
                        cols_diff_reject = cols_diff < 0
                        cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                        index1 = apply(cols_diff_reject, 1, all)
                        gene_df <<- df[index1, ] %>% dplyr::select(name)
                      }
                    }
                }
              } else {
                if(input$import_for_go == "UPDOWN for DEG-RNAseq") {
                  if(!input$df_with_lg2fc) {
                    gene_df <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% dplyr::select(name)
                    } else {
                      gene_df <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% dplyr::select(name, paste(input$import_contrast_for_go, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_go, "log2FoldChange", sep = "_"))
                    }
                  } else {
                    if(input$import_for_go == "UPregu for DEG-RNAseq") {
                      if(!input$df_with_lg2fc) {
                        gene_df <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "log2FoldChange", sep = "_")) > 0) %>% dplyr::select(name)
                        } else {
                          gene_df <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "log2FoldChange", sep = "_")) > 0) %>% dplyr::select(name, paste(input$import_contrast_for_go, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_go, "log2FoldChange", sep = "_"))
                        }
                      } else {
                        if(input$import_for_go == "DOWNregu for DEG-RNAseq") {
                          if(!input$df_with_lg2fc) {
                            gene_df <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "log2FoldChange", sep = "_")) < 0) %>% dplyr::select(name)
                            } else {
                              gene_df <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_go, "log2FoldChange", sep = "_")) < 0) %>% dplyr::select(name, paste(input$import_contrast_for_go, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_go, "log2FoldChange", sep = "_"))
                            }
                        }
                      }
                  }
              }
          }

        }
    }


      if(!input$import_from_for_go) {
        gene_df_check <- try(gene_df())
        if(!class(gene_df_check) == "try-error") {
           gene_df_check$name = rm_digit_end(gene_df_check$name)
           if(ncol(gene_df_check) == 1) {
             gene_df <- data.frame(name = gene_df_check[!duplicated(gene_df_check$name),], stringsAsFactors = F)
           }
           if(ncol(gene_df_check) == 2) {
             gene_df <- gene_df_check[!duplicated(gene_df_check$name),]
           }

        }
  }

      if(input$import_from_for_go) {
           gene_df$name = rm_digit_end(gene_df$name)
           if(ncol(gene_df) == 1) {
             gene_df <- data.frame(name = gene_df[!duplicated(gene_df$name),], stringsAsFactors = F)
           }
           if(ncol(gene_df) == 2) {
             gene_df <- gene_df[!duplicated(gene_df$name),]
           }
  }
      # reat <<- reactive({

        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Please wait ...", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        reat <<- try(goAnalysis(df = gene_df, df_with_lg2fc = input$df_with_lg2fc, organism = input$organism_for_go), silent = TRUE)
        # })

     # test1 <<- reat()
    res <<- reactive({
      try(giveGO_res_and_table(reat = reat, ont = input$go_ont, pCutoff = input$go_p, p.adj.cutoff = input$go_padj, q.cutoff = input$go_qvalue, simplify = input$go_simplify), silent = TRUE)
      })

    test2 <<- res()

    output$go_Table <- DT::renderDataTable({
    if(!input$import_from_for_go) {
      shiny::validate(
      need(!class(gene_df_check) == "try-error", message = "No genes meet your requirements, and can not do the GO analysis")
    )
    }

    if(input$import_from_for_go) {
      shiny::validate(
      need(nrow(gene_df) != 0, message = "No genes meet your requirements, and can not do the GO analysis")
    )
    }

    validate( need(reat, "No genes can be mapped to ENTREZID, and can not do the GO analysis"))
      DT::datatable(res()$sig_table, filter = 'top', options = list( autoWidth = F,scrollX = TRUE))
    })

    # output$table <- DT::renderDataTable({
    #   table()
    # }, options = list(pageLength = 25, scrollX = T),
    # selection = list(selected = c(1)))


    output$significantBox_for_go <- renderInfoBox({
    if(!input$import_from_for_go) {
      shiny::validate(
      need(!class(gene_df_check) == "try-error", message = "No genes meet your requirements, and can not do the GO analysis")
    )
    }

    if(input$import_from_for_go) {
      shiny::validate(
      need(nrow(gene_df) != 0, message = "No genes meet your requirements, and can not do the GO analysis")
    )
    }
    validate( need(reat, "No genes can be mapped to ENTREZID, and can not do the GO analysis"))

      num_total_go <- res()$all_table %>%
        nrow()
      num_signif_go <- res()$sig_table %>%
        nrow()
      frac_go <- num_signif_go / num_total_go

        if(frac_go == 0) {
        info_box_go <- infoBox("Significant terms",
                            paste0(num_signif_go,
                                   " out of ",
                                   num_total_go),
                            "No terms enriched",
                            icon = icon("thumbs-down", lib = "glyphicon"),
                            color = "red",
                            width = 4)
       }
      if(!frac_go == 0) {
        info_box_go <-     infoBox("Significant terms",
                              paste0(num_signif_go,
                                     " out of ",
                                     num_total_go),
                              paste0(signif(frac_go * 100, digits = 3),
                                     "% of terms enriched"),
                              icon = icon("thumbs-up", lib = "glyphicon"),
                              color = "green",
                              width = 4)
      }
      info_box_go
    })

    go_barplot_input <- reactive({
      # # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        my_barplot(res = res(), ShowCategory = input$go_ShowCategory_bar, color = input$go_color, ont = input$go_ont, Split = input$go_bar_if_Split_for_ont_ALL)
    })

    go_dotplot_input <- reactive({
      progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        my_dotplot(res = res(), ShowCategory = input$go_ShowCategory_dot, color = input$go_color, ont = input$go_ont, Split = input$go_dot_if_Split_for_ont_ALL)
      })

    go_dotplot_opt_input <- reactive({
      validate(
      need(input$go_ont != "ALL", "Please go to the panel : Dot plot ")
    )
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        my_dotplot_opt(res = res(), color = input$go_color, size = "Count", title = "", decreasing = TRUE, ShowCategory = input$go_ShowCategory_dot_opt)
      })

    go_heatplot_input <- reactive({
      # progress <- shiny::Progress$new()
      #   progress$set(message = "Plotting", value = 0.66)
      #   # Close the progress when this reactive exits (even if there's an error)
      #   on.exit(progress$close())
      withProgress(message = 'Plotting', value = 0.66, {
      my_heatplot(res = res(), ShowCategory = input$go_ShowCategory_heat, df_with_lg2fc = input$df_with_lg2fc, ont = input$go_ont)
        })
      })

    go_cnetplot_input <- reactive({
      # progress <- shiny::Progress$new()
      #   progress$set(message = "Plotting", value = 0.66)
      #   # Close the progress when this reactive exits (even if there's an error)
      #   on.exit(progress$close())
      withProgress(message = 'Plotting', value = 0.66, {
      my_cnetplot(res = res(), ShowCategory = input$go_ShowCategory_cnet, circular = input$go_circular_cnet, colorEdge = TRUE, df_with_lg2fc = input$df_with_lg2fc, ont = input$go_ont)
        })
      })

    go_emaplot_input <- reactive({
      # progress <- shiny::Progress$new()
      #   progress$set(message = "Plotting", value = 0.66)
      #   # Close the progress when this reactive exits (even if there's an error)
      #   on.exit(progress$close())
      withProgress(message = 'Plotting', value = 0.66, {
      my_emaplot(res = res(), ShowCategory = input$go_ShowCategory_ema, color = input$go_color, layout = "kk", ont = input$go_ont)
        })
      })

    go_goplot_input <- reactive({
      validate(
      need(input$go_ont != "ALL", "Ontology ALL: can not plot for goplot")
    )
      # progress <- shiny::Progress$new()
      #   progress$set(message = "Plotting", value = 0.66)
      #   # Close the progress when this reactive exits (even if there's an error)
      #   on.exit(progress$close())
      withProgress(message = 'Plotting', value = 0.66, {
      my_goplot(res = res(), ShowCategory = input$go_ShowCategory_go, color = input$go_color, ont = input$go_ont, Layout = "kk",circular = input$go_circular_go)
        })
      })

    # go_GOgraphplot_input <- reactive({
    #   validate(
    #   need(input$go_ont != "ALL", "Ontology ALL: can not plot for plotGOgraph")
    # )
    #   # progress <- shiny::Progress$new()
    #   #   progress$set(message = "Plotting", value = 0.66)
    #   #   # Close the progress when this reactive exits (even if there's an error)
    #   #   on.exit(progress$close())
    #   # withProgress(message = 'Plotting', value = 0.66, {
    #   my_plotGOgraph(res = res(), firstSigNodes = input$go_ShowCategory_GOgraph, ont = input$go_ont)
    #     # })
    #   })

    output$go_barplot <- renderPlot({
      go_barplot_input()
    })

    output$go_dotplot <- renderPlot({
      go_dotplot_input()
    })

    output$go_dotplot_opt <- renderPlot({
      go_dotplot_opt_input()
    })

    output$go_heatplot <- renderPlot({
      go_heatplot_input()
    })

     output$go_cnetplot <- renderPlot({
      go_cnetplot_input()
    })

    output$go_emaplot <- renderPlot({
      go_emaplot_input()
    })

    output$go_goplot <- renderPlot({
      go_goplot_input()
    })

    #  output$go_GOgraphplot <- renderPlot({
    #   go_GOgraphplot_input()
    # })


    ### Download objects and functions ### ------------------------------------
    datasetInput_for_go <- reactive({
      table_for_go = res()
      if(input$go_ont == "ALL") {
        if(input$go_bar_if_Split_for_ont_ALL | input$go_dot_if_Split_for_ont_ALL){
          switch(input$dataset_for_go,
             "full_results" = res()$all_table,
             "significant_results" = res()$sig_table)

          } else {
            switch(input$dataset_for_go,
             "full_results" = res()$all_table %>% dplyr::arrange(p.adjust),
             "significant_results" = res()$sig_table %>% dplyr::arrange(p.adjust))

          }

        } else {
          switch(input$dataset_for_go,
             "full_results" = res()$all_table,
             "significant_results" = res()$sig_table)
        }
      # switch(input$dataset_for_go,
      #        "full_results" = res()$all_table,
      #        "significant_results" = res()$sig_table)
    })

    output$downloadgo <- downloadHandler(
      filename = function() { paste(input$dataset_for_go, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput_for_go(),
        file,
        col.names = TRUE,
        row.names = FALSE,
        sep ="\t") }
    )

    output$download_go_barplot <- downloadHandler(
      filename = 'barplot.pdf',
      content = function(file) {
        pdf(file, width = input$go_wide_bar, height = input$go_high_bar)
        print(go_barplot_input())
        dev.off()
      }
    )

    output$download_go_dotplot <- downloadHandler(
      filename = 'dotplot.pdf',
      content = function(file) {
        pdf(file, width = input$go_wide_dot, height = input$go_high_dot)
        print(go_dotplot_input())
        dev.off()
      }
    )

    output$download_go_dotplot_opt <- downloadHandler(
      filename = 'dotplot_opt.pdf',
      content = function(file) {
        pdf(file, width = input$go_wide_dot_opt, height = input$go_high_dot_opt)
        print(go_dotplot_opt_input())
        dev.off()
      }
    )

    output$download_go_heatplot <- downloadHandler(
      filename = 'heatplot.pdf',
      content = function(file) {
        pdf(file, width = input$go_wide_heat, height = input$go_high_heat)
        print(go_heatplot_input())
        dev.off()
      }
    )

    output$download_go_cnetplot <- downloadHandler(
      filename = 'cnetplot.pdf',
      content = function(file) {
        pdf(file, width = input$go_wide_cnet, height = input$go_high_cnet)
        print(go_cnetplot_input())
        dev.off()
      }
    )

    output$download_go_emaplot <- downloadHandler(
      filename = 'emaplot.pdf',
      content = function(file) {
        pdf(file, width = input$go_wide_ema, height = input$go_high_ema)
        print(go_emaplot_input())
        dev.off()
      }
    )

    output$download_go_goplot <- downloadHandler(
      filename = 'goplot.pdf',
      content = function(file) {
        pdf(file, width = input$go_wide_go, height = input$go_high_go)
        print(go_goplot_input())
        dev.off()
      }
    )

    # output$download_go_GOgraphplot <- downloadHandler(
    #   filename = 'GOgraphplot.pdf',
    #   content = function(file) {
    #     pdf(file, width = input$go_wide_GOgraph, height = input$go_high_GOgraph)
    #     print(go_GOgraphplot_input())
    #     dev.off()
    #   }
    # )

    observeEvent(input$plot_for_GOgraphplot, {
      go_GOgraphplot_input <- reactive({
      validate(
      need(input$go_ont != "ALL", "Ontology ALL: can not plot for plotGOgraph")
    )
      # progress <- shiny::Progress$new()
      #   progress$set(message = "Plotting", value = 0.66)
      #   # Close the progress when this reactive exits (even if there's an error)
      #   on.exit(progress$close())
      withProgress(message = 'Plotting', value = 0.66, {
      my_plotGOgraph(res = res(), firstSigNodes = input$go_ShowCategory_GOgraph, ont = input$go_ont)
        })
      })

       output$go_GOgraphplot <- renderPlot({
      go_GOgraphplot_input()
    })

       output$download_go_GOgraphplot <- downloadHandler(
      filename = 'GOgraphplot.pdf',
      content = function(file) {
        pdf(file, width = input$go_wide_GOgraph, height = input$go_high_GOgraph)
        print(go_GOgraphplot_input())
        dev.off()
      }
    )
      })
  })

  ## for kegg analysis
output$import_for_kegg <- renderUI({
    if(input$import_from_for_kegg){
                 radioButtons("import_for_kegg",
                              "Choose genes",
                              choices = c("UPregu for DEP-LFQ", "DOWNregu for DEP-LFQ", "UPDOWN for DEP-LFQ", "UPregu for DEG-RNAseq", "DOWNregu for DEG-RNAseq", "UPDOWN for DEG-RNAseq"),
                              selected = "UPregu for DEP-LFQ")
    }
    })

output$import_contrast_for_kegg <- renderUI({
  if(input$import_from_for_kegg) {
if(input$import_for_kegg == "UPregu for DEP-LFQ" | input$import_for_kegg == "DOWNregu for DEP-LFQ" | input$import_for_kegg == "UPDOWN for DEP-LFQ") {
  validate(need(!is.null(data()), "Please go to the DEP-LFQ options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_kegg", "Contrast", choices = c(get_results(dep()) %>%
  dplyr::filter(significant) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique(), "Any significant"), selected = NULL)
    } else {
      validate(need(!is.null(countData()), "Please go to the DEG-RNAseq options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_kegg", "Contrast", choices = c(as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique(), "Any significant"), selected = NULL)
    }
  } else {
    return(NULL)
  }

  })

output$text_input_for_kegg <- renderUI({
  if(!input$import_from_for_kegg){
    textAreaInput(inputId = "text_input_for_kegg", label = "Please paste your gene list", placeholder = "TP53\nPTEN\n\nor\n\nTP53  6.21\nPTEN  -1.53", rows = 8, width = "100%")
  }
  })

output$organism_for_kegg <- renderUI({
  # if(!input$import_from_for_kegg){
    selectizeInput("organism_for_kegg", "Select organism", choices=c("hsa","mmu"), selected="hsa")
  # }
  })

output$df_with_lg2fc_for_kegg <- renderUI({
  # if(!input$import_from_for_kegg) {
    checkboxInput("df_with_lg2fc_for_kegg", "If with log2 fold change", value = FALSE)
  # }
  })

  # observeEvent for KEGG analysis input
  observeEvent(input$analyze_for_kegg,{
      ### Interactive UI functions ### ------------------------------------------
    output$downloadTable_kegg <- renderUI({
      selectizeInput("dataset_for_kegg",
                     "Choose a dataset to save" ,
                     c("full_results","significant_results"
                       ))
    })

    output$downloadButton_kegg <- renderUI({
      downloadButton('downloadkegg', 'Save table', class = "downloadkegg")
    })

    ## reactive functions
    if(!input$import_from_for_kegg) {
        genelist_kegg <- reactive({ strsplit(input$text_input_for_kegg,'\n')[[1]] })
        gene_name_kegg <-reactive(unlist(strsplit(genelist_kegg(),";")))

        gene_df_kegg <- reactive({
          gene_df <<- data.frame(name=gene_name_kegg())
          gene_df$name <-  sapply(gene_df$name,the1stname)
          gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
          rownames(gene_df) = NULL
          if(!input$df_with_lg2fc_for_kegg){
          colnames(gene_df) = "name"
          gene_df = as.data.frame(gene_df)
        } else {
          colnames(gene_df) = c("name","fc")
          gene_df = as.data.frame(gene_df)
          gene_df$fc = as.numeric(as.character(gene_df$fc))
          gene_df$name = as.character(gene_df$name)
          # test <<- gene_df
        }
          gene_df
          })
    }

    if(input$import_from_for_kegg) {
      if(input$import_for_kegg == "UPregu for DEP-LFQ" | input$import_for_kegg == "DOWNregu for DEP-LFQ" | input$import_for_kegg == "UPDOWN for DEP-LFQ") {
        if(input$import_contrast_for_kegg == "Any significant") {
          if(input$import_for_kegg == "UPDOWN for DEP-LFQ") {
            gene_df_kegg <<- get_results(dep()) %>% dplyr::filter(significant) %>% dplyr::select(name)
            } else {
              if(input$import_for_kegg == "UPregu for DEP-LFQ") {
                index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
                df <- get_results(dep()) %>% dplyr::filter(significant)
                cols_diff = df %>% dplyr::select(index)
                cols_diff_reject = cols_diff > 0
                cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                index1 = apply(cols_diff_reject, 1, all)
                gene_df_kegg <<- df[index1, ] %>% dplyr::select(name)
                } else {
                  if(input$import_for_kegg == "DOWNregu for DEP-LFQ") {
                    index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
                    df <- get_results(dep()) %>% dplyr::filter(significant)
                    cols_diff = df %>% dplyr::select(index)
                    cols_diff_reject = cols_diff < 0
                    cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                    index1 = apply(cols_diff_reject, 1, all)
                    gene_df_kegg <<- df[index1, ] %>% dplyr::select(name)
                  }
                }
            }
          } else {
            if(input$import_for_kegg == "UPDOWN for DEP-LFQ") {
              if(!input$df_with_lg2fc_for_kegg) {
                gene_df_kegg <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% dplyr::select(name)
                } else {
                  gene_df_kegg <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% dplyr::select(name, paste(input$import_contrast_for_kegg, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_kegg, "ratio", sep = "_"))
                }
              } else {
                if(input$import_for_kegg == "UPregu for DEP-LFQ") {
                  if(!input$df_with_lg2fc_for_kegg) {
                    gene_df_kegg <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_kegg, "ratio", sep = "_")) > 0) %>% dplyr::select(name)
                    } else {
                      gene_df_kegg <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_kegg, "ratio", sep = "_")) > 0) %>% dplyr::select(name, paste(input$import_contrast_for_kegg, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_kegg, "ratio", sep = "_"))
                    }
                  } else {
                    if(input$import_for_kegg == "DOWNregu for DEP-LFQ") {
                      if(!input$df_with_lg2fc_for_kegg) {
                        gene_df_kegg <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_kegg, "ratio", sep = "_")) < 0) %>% dplyr::select(name)
                        } else {
                          gene_df_kegg <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_kegg, "ratio", sep = "_")) < 0) %>% dplyr::select(name, paste(input$import_contrast_for_kegg, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_kegg, "ratio", sep = "_"))
                        }
                    }
                  }
              }
          }
        } else {
          if(input$import_for_kegg == "UPregu for DEG-RNAseq" | input$import_for_kegg == "DOWNregu for DEG-RNAseq" | input$import_for_kegg == "UPDOWN for DEG-RNAseq") {
            if(input$import_contrast_for_kegg == "Any significant") {
              if(input$import_for_kegg == "UPDOWN for DEG-RNAseq") {
                gene_df_kegg <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant) %>% dplyr::select(name)
                } else {
                  if(input$import_for_kegg == "UPregu for DEG-RNAseq") {
                    index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
                    df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
                    cols_diff = df %>% dplyr::select(index)
                    cols_diff_reject = cols_diff > 0
                    cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                    index1 = apply(cols_diff_reject, 1, all)
                    gene_df_kegg <<- df[index1, ] %>% dplyr::select(name)
                    } else {
                      if(input$import_for_kegg == "DOWNregu for DEG-RNAseq") {
                        index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
                        df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
                        cols_diff = df %>% dplyr::select(index)
                        cols_diff_reject = cols_diff < 0
                        cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                        index1 = apply(cols_diff_reject, 1, all)
                        gene_df_kegg <<- df[index1, ] %>% dplyr::select(name)
                      }
                    }
                }
              } else {
                if(input$import_for_kegg == "UPDOWN for DEG-RNAseq") {
                  if(!input$df_with_lg2fc_for_kegg) {
                    gene_df_kegg <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% dplyr::select(name)
                    } else {
                      gene_df_kegg <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% dplyr::select(name, paste(input$import_contrast_for_kegg, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_kegg, "log2FoldChange", sep = "_"))
                    }
                  } else {
                    if(input$import_for_kegg == "UPregu for DEG-RNAseq") {
                      if(!input$df_with_lg2fc_for_kegg) {
                        gene_df_kegg <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "log2FoldChange", sep = "_")) > 0) %>% dplyr::select(name)
                        } else {
                          gene_df_kegg <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "log2FoldChange", sep = "_")) > 0) %>% dplyr::select(name, paste(input$import_contrast_for_kegg, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_kegg, "log2FoldChange", sep = "_"))
                        }
                      } else {
                        if(input$import_for_kegg == "DOWNregu for DEG-RNAseq") {
                          if(!input$df_with_lg2fc_for_kegg) {
                            gene_df_kegg <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "log2FoldChange", sep = "_")) < 0) %>% dplyr::select(name)
                            } else {
                              gene_df_kegg <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_kegg, "log2FoldChange", sep = "_")) < 0) %>% dplyr::select(name, paste(input$import_contrast_for_kegg, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_kegg, "log2FoldChange", sep = "_"))
                            }
                        }
                      }
                  }
              }
          }

        }
    }

          if(!input$import_from_for_kegg) {
        gene_df_kegg_check <- try(gene_df_kegg())
        if(!class(gene_df_kegg_check) == "try-error") {
           gene_df_kegg_check$name = rm_digit_end(gene_df_kegg_check$name)
           if(ncol(gene_df_kegg_check) == 1) {
             gene_df_kegg <- data.frame(name = gene_df_kegg_check[!duplicated(gene_df_kegg_check$name),], stringsAsFactors = F)
           }
           if(ncol(gene_df_kegg_check) == 2) {
             gene_df_kegg <- gene_df_kegg_check[!duplicated(gene_df_kegg_check$name),]
           }

        }
  }

      if(input$import_from_for_kegg) {
           gene_df_kegg$name = rm_digit_end(gene_df_kegg$name)
           if(ncol(gene_df_kegg) == 1) {
             gene_df_kegg <- data.frame(name = gene_df_kegg[!duplicated(gene_df_kegg$name),], stringsAsFactors = F)
           }
           if(ncol(gene_df_kegg) == 2) {
             gene_df_kegg <- gene_df_kegg[!duplicated(gene_df_kegg$name),]
           }
  }

      # reat_kegg <<- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Please wait ...", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        reat_kegg <- try(keggAnalysis(df = gene_df_kegg, organism = input$organism_for_kegg, df_with_lg2fc = input$df_with_lg2fc_for_kegg), silent = TRUE)
        # })

     # test1 <<- reat_kegg()
    res_kegg <<- reactive({
    # if(input$import_from_for_kegg) {
    #   shiny::validate(
    #   need(nrow(gene_df_kegg) != 0, message = "No genes meet your requirements, and can not do the KEGG analysis")
    # )
    # }
    #    # validate( need(reat_kegg, "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
    #    shiny::validate( need(class(reat_kegg) != "try-error", message = "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
      try(givekegg_reat_res_and_table(reat = reat_kegg, pCutoff = input$kegg_p, p.adj.cutoff = input$kegg_padj, q.cutoff = input$kegg_qvalue), silent = TRUE)
      })

    test2 <<- res_kegg()


    # output$kegg_Table <- DT::renderDataTable(
    #   res_kegg()$sig_table, filter = 'top', options = list( autoWidth = F,scrollX = TRUE
    #   )
    # )
   output$kegg_Table <- DT::renderDataTable({
      if(!input$import_from_for_kegg) {
      shiny::validate(
      need(!class(gene_df_kegg_check) == "try-error", message = "No genes meet your requirements, and can not do the KEGG analysis")
    )
    }

    if(input$import_from_for_kegg) {
      shiny::validate(
      need(nrow(gene_df_kegg) != 0, message = "No genes meet your requirements, and can not do the KEGG analysis")
    )
    }
       validate( need(reat_kegg, "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
       # shiny::validate( need(class(reat_kegg) != "try-error", message = "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
       DT::datatable(res_kegg()$sig_table, filter = 'top', options = list( autoWidth = F,scrollX = TRUE))
    })

    output$significantBox_for_kegg <- renderInfoBox({
      if(!input$import_from_for_kegg) {
      shiny::validate(
      need(!class(gene_df_kegg_check) == "try-error", message = "No genes meet your requirements, and can not do the KEGG analysis")
    )
    }

    if(input$import_from_for_kegg) {
      shiny::validate(
      need(nrow(gene_df_kegg) != 0, message = "No genes meet your requirements, and can not do the KEGG analysis")
    )
    }
    validate( need(reat_kegg, "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
      num_total_kegg <- res_kegg()$all_table %>%
        nrow()
      num_signif_kegg <- res_kegg()$sig_table %>%
        nrow()
      frac_kegg <- num_signif_kegg / num_total_kegg

        if(frac_kegg == 0) {
        info_box_kegg <- infoBox("Significant terms",
                            paste0(num_signif_kegg,
                                   " out of ",
                                   num_total_kegg),
                            "No terms enriched",
                            icon = icon("thumbs-down", lib = "glyphicon"),
                            color = "red",
                            width = 4)
       }
      if(!frac_kegg == 0) {
        info_box_kegg <-     infoBox("Significant terms",
                              paste0(num_signif_kegg,
                                     " out of ",
                                     num_total_kegg),
                              paste0(signif(frac_kegg * 100, digits = 3),
                                     "% of terms enriched"),
                              icon = icon("thumbs-up", lib = "glyphicon"),
                              color = "green",
                              width = 4)
      }
      info_box_kegg
    })

    kegg_barplot_input <- reactive({
      # # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        barplot(res_kegg()$sig_res, showCategory = input$kegg_ShowCategory_bar, color = input$kegg_color)
    })

    kegg_dotplot_input <- reactive({
      progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        dotplot(res_kegg()$sig_res, showCategory = input$kegg_ShowCategory_dot, color = input$kegg_color)
      })

    kegg_dotplot_opt_input <- reactive({
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        my_dotplot_opt(res = res_kegg(), color = input$kegg_color, size = "Count", title = "", decreasing = TRUE, ShowCategory = input$kegg_ShowCategory_dot_opt)
      })

    kegg_heatplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      heatplot_for_react_kegg(res = res_kegg(), ShowCategory = input$kegg_ShowCategory_heat, df_with_lg2fc = input$df_with_lg2fc_for_kegg)
        })
      })

    kegg_cnetplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      cnetplot_for_react_kegg(res = res_kegg(), ShowCategory = input$kegg_ShowCategory_cnet, circular = input$kegg_circular_cnet, colorEdge = TRUE, df_with_lg2fc = input$df_with_lg2fc_for_kegg)
        })
      })

    kegg_emaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      emaplot_for_react_kegg(res = res_kegg(), ShowCategory = input$kegg_ShowCategory_ema, color = input$kegg_color, layout = "kk")
        })
      })

    output$kegg_barplot <- renderPlot({
      kegg_barplot_input()
    })

    output$kegg_dotplot <- renderPlot({
      kegg_dotplot_input()
    })

    output$kegg_dotplot_opt <- renderPlot({
      kegg_dotplot_opt_input()
    })

    output$kegg_heatplot <- renderPlot({
      kegg_heatplot_input()
    })

     output$kegg_cnetplot <- renderPlot({
      kegg_cnetplot_input()
    })

    output$kegg_emaplot <- renderPlot({
      kegg_emaplot_input()
    })

    ### Download objects and functions ### ------------------------------------
    datasetInput_for_kegg <- reactive({
      table_for_kegg = res_kegg()
          switch(input$dataset_for_kegg,
             "full_results" = res_kegg()$all_table,
             "significant_results" = res_kegg()$sig_table)
    })

    output$downloadkegg <- downloadHandler(
      filename = function() { paste(input$dataset_for_kegg, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput_for_kegg(),
        file,
        col.names = TRUE,
        row.names = FALSE,
        sep ="\t") }
    )

    output$download_kegg_barplot <- downloadHandler(
      filename = 'barplot.pdf',
      content = function(file) {
        pdf(file, width = input$kegg_wide_bar, height = input$kegg_high_bar)
        print(kegg_barplot_input())
        dev.off()
      }
    )

    output$download_kegg_dotplot <- downloadHandler(
      filename = 'dotplot.pdf',
      content = function(file) {
        pdf(file, width = input$kegg_wide_dot, height = input$kegg_high_dot)
        print(kegg_dotplot_input())
        dev.off()
      }
    )

    output$download_kegg_dotplot_opt <- downloadHandler(
      filename = 'dotplot_opt.pdf',
      content = function(file) {
        pdf(file, width = input$kegg_wide_dot_opt, height = input$kegg_high_dot_opt)
        print(kegg_dotplot_opt_input())
        dev.off()
      }
    )

    output$download_kegg_heatplot <- downloadHandler(
      filename = 'heatplot.pdf',
      content = function(file) {
        pdf(file, width = input$kegg_wide_heat, height = input$kegg_high_heat)
        print(kegg_heatplot_input())
        dev.off()
      }
    )

    output$download_kegg_cnetplot <- downloadHandler(
      filename = 'cnetplot.pdf',
      content = function(file) {
        pdf(file, width = input$kegg_wide_cnet, height = input$kegg_high_cnet)
        print(kegg_cnetplot_input())
        dev.off()
      }
    )

    output$download_kegg_emaplot <- downloadHandler(
      filename = 'emaplot.pdf',
      content = function(file) {
        pdf(file, width = input$kegg_wide_ema, height = input$kegg_high_ema)
        print(kegg_emaplot_input())
        dev.off()
      }
    )
  })

  ## for reactome analysis
output$import_for_reactome <- renderUI({
    if(input$import_from_for_reactome){
                 radioButtons("import_for_reactome",
                              "Choose genes",
                              choices = c("UPregu for DEP-LFQ", "DOWNregu for DEP-LFQ", "UPDOWN for DEP-LFQ", "UPregu for DEG-RNAseq", "DOWNregu for DEG-RNAseq", "UPDOWN for DEG-RNAseq"),
                              selected = "UPregu for DEP-LFQ")
    }
    })

output$import_contrast_for_reactome <- renderUI({
  if(input$import_from_for_reactome) {
if(input$import_for_reactome == "UPregu for DEP-LFQ" | input$import_for_reactome == "DOWNregu for DEP-LFQ" | input$import_for_reactome == "UPDOWN for DEP-LFQ") {
  validate(need(!is.null(data()), "Please go to the DEP-LFQ options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_reactome", "Contrast", choices = c(get_results(dep()) %>%
  dplyr::filter(significant) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique(), "Any significant"), selected = NULL)
    } else {
      validate(need(!is.null(countData()), "Please go to the DEG-RNAseq options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_reactome", "Contrast", choices = c(as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique(), "Any significant"), selected = NULL)
    }
  } else {
    return(NULL)
  }

  })

output$text_input_for_reactome <- renderUI({
  if(!input$import_from_for_reactome){
    textAreaInput(inputId = "text_input_for_reactome", label = "Please paste your gene list", placeholder = "TP53\nPTEN\n\nor\n\nTP53  6.21\nPTEN  -1.53", rows = 8, width = "100%")
  }
  })

output$organism_for_reactome <- renderUI({
  # if(!input$import_from_for_reactome){
    selectizeInput("organism_for_reactome", "Select organism", choices=c("human","mouse"), selected="human")
  # }
  })

output$df_with_lg2fc_for_reactome <- renderUI({
  # if(!input$import_from_for_reactome) {
    checkboxInput("df_with_lg2fc_for_reactome", "If with log2 fold change", value = FALSE)
  # }
  })

  # observeEvent for reactome analysis input
  observeEvent(input$analyze_for_reactome,{
      ### Interactive UI functions ### ------------------------------------------
    output$downloadTable_reactome <- renderUI({
      selectizeInput("dataset_for_reactome",
                     "Choose a dataset to save" ,
                     c("full_results","significant_results"
                       ))
    })

    output$downloadButton_reactome <- renderUI({
      downloadButton('downloadreactome', 'Save table', class = "downloadreactome")
    })

    ## reactive functions
    if(!input$import_from_for_reactome) {
      genelist_reactome <- reactive({ strsplit(input$text_input_for_reactome,'\n')[[1]] })
      gene_name_reactome <-reactive(unlist(strsplit(genelist_reactome(),";")))

      gene_df_reactome <- reactive({
        gene_df <<- data.frame(name=gene_name_reactome())
        gene_df$name <-  sapply(gene_df$name,the1stname)
        gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
        rownames(gene_df) = NULL
        if(!input$df_with_lg2fc_for_reactome){
        colnames(gene_df) = "name"
        gene_df = as.data.frame(gene_df)
      } else {
        colnames(gene_df) = c("name","fc")
        gene_df = as.data.frame(gene_df)
        gene_df$fc = as.numeric(as.character(gene_df$fc))
        gene_df$name = as.character(gene_df$name)
        # test <<- gene_df
      }
        gene_df
        })
    }

    if(input$import_from_for_reactome) {
      if(input$import_for_reactome == "UPregu for DEP-LFQ" | input$import_for_reactome == "DOWNregu for DEP-LFQ" | input$import_for_reactome == "UPDOWN for DEP-LFQ") {
        if(input$import_contrast_for_reactome == "Any significant") {
          if(input$import_for_reactome == "UPDOWN for DEP-LFQ") {
            gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(significant) %>% dplyr::select(name)
            } else {
              if(input$import_for_reactome == "UPregu for DEP-LFQ") {
                index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
                df <- get_results(dep()) %>% dplyr::filter(significant)
                cols_diff = df %>% dplyr::select(index)
                cols_diff_reject = cols_diff > 0
                cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                index1 = apply(cols_diff_reject, 1, all)
                gene_df_reactome <<- df[index1, ] %>% dplyr::select(name)
                } else {
                  if(input$import_for_reactome == "DOWNregu for DEP-LFQ") {
                    index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
                    df <- get_results(dep()) %>% dplyr::filter(significant)
                    cols_diff = df %>% dplyr::select(index)
                    cols_diff_reject = cols_diff < 0
                    cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                    index1 = apply(cols_diff_reject, 1, all)
                    gene_df_reactome <<- df[index1, ] %>% dplyr::select(name)
                  }
                }
            }
          } else {
            if(input$import_for_reactome == "UPDOWN for DEP-LFQ") {
              if(!input$df_with_lg2fc_for_reactome) {
                gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::select(name)
                } else {
                  gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "ratio", sep = "_"))
                }
              } else {
                if(input$import_for_reactome == "UPregu for DEP-LFQ") {
                  if(!input$df_with_lg2fc_for_reactome) {
                    gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_reactome, "ratio", sep = "_")) > 0) %>% dplyr::select(name)
                    } else {
                      gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_reactome, "ratio", sep = "_")) > 0) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "ratio", sep = "_"))
                    }
                  } else {
                    if(input$import_for_reactome == "DOWNregu for DEP-LFQ") {
                      if(!input$df_with_lg2fc_for_reactome) {
                        gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_reactome, "ratio", sep = "_")) < 0) %>% dplyr::select(name)
                        } else {
                          gene_df_reactome <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_reactome, "ratio", sep = "_")) < 0) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "ratio", sep = "_"))
                        }
                    }
                  }
              }
          }
        } else {
          if(input$import_for_reactome == "UPregu for DEG-RNAseq" | input$import_for_reactome == "DOWNregu for DEG-RNAseq" | input$import_for_reactome == "UPDOWN for DEG-RNAseq") {
            if(input$import_contrast_for_reactome == "Any significant") {
              if(input$import_for_reactome == "UPDOWN for DEG-RNAseq") {
                gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant) %>% dplyr::select(name)
                } else {
                  if(input$import_for_reactome == "UPregu for DEG-RNAseq") {
                    index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
                    df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
                    cols_diff = df %>% dplyr::select(index)
                    cols_diff_reject = cols_diff > 0
                    cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                    index1 = apply(cols_diff_reject, 1, all)
                    gene_df_reactome <<- df[index1, ] %>% dplyr::select(name)
                    } else {
                      if(input$import_for_reactome == "DOWNregu for DEG-RNAseq") {
                        index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
                        df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
                        cols_diff = df %>% dplyr::select(index)
                        cols_diff_reject = cols_diff < 0
                        cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                        index1 = apply(cols_diff_reject, 1, all)
                        gene_df_reactome <<- df[index1, ] %>% dplyr::select(name)
                      }
                    }
                }
              } else {
                if(input$import_for_reactome == "UPDOWN for DEG-RNAseq") {
                  if(!input$df_with_lg2fc_for_reactome) {
                    gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::select(name)
                    } else {
                      gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_"))
                    }
                  } else {
                    if(input$import_for_reactome == "UPregu for DEG-RNAseq") {
                      if(!input$df_with_lg2fc_for_reactome) {
                        gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) > 0) %>% dplyr::select(name)
                        } else {
                          gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) > 0) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_"))
                        }
                      } else {
                        if(input$import_for_reactome == "DOWNregu for DEG-RNAseq") {
                          if(!input$df_with_lg2fc_for_reactome) {
                            gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) < 0) %>% dplyr::select(name)
                            } else {
                              gene_df_reactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) < 0) %>% dplyr::select(name, paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_reactome, "log2FoldChange", sep = "_"))
                            }
                        }
                      }
                  }
              }
          }

        }
    }

          if(!input$import_from_for_reactome) {
        gene_df_reactome_check <<- try(gene_df_reactome())
        if(!class(gene_df_reactome_check) == "try-error") {
           gene_df_reactome_check$name = rm_digit_end(gene_df_reactome_check$name)
           if(ncol(gene_df_reactome_check) == 1) {
             gene_df_reactome <- data.frame(name = gene_df_reactome_check[!duplicated(gene_df_reactome_check$name),], stringsAsFactors = F)
           }
           if(ncol(gene_df_reactome_check) == 2) {
             gene_df_reactome <- gene_df_reactome_check[!duplicated(gene_df_reactome_check$name),]
           }

        }
  }

      if(input$import_from_for_reactome) {
           gene_df_reactome$name = rm_digit_end(gene_df_reactome$name)
           if(ncol(gene_df_reactome) == 1) {
             gene_df_reactome <- data.frame(name = gene_df_reactome[!duplicated(gene_df_reactome$name),], stringsAsFactors = F)
           }
           if(ncol(gene_df_reactome) == 2) {
             gene_df_reactome <- gene_df_reactome[!duplicated(gene_df_reactome$name),]
           }
  }

      # reat_reactome <<- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Please wait ...", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        reat_reactome <<- try(reactAnalysis(df = gene_df_reactome, organism = input$organism_for_reactome, df_with_lg2fc = input$df_with_lg2fc_for_reactome), silent = TRUE)
        # })

     # test1 <<- reat_reactome()
    res_reactome <<- reactive({
      try(givekegg_reat_res_and_table(reat = reat_reactome, pCutoff = input$reactome_p, p.adj.cutoff = input$reactome_padj, q.cutoff = input$reactome_qvalue), silent = TRUE)
      })

    test2 <<- res_reactome()

    output$reactome_Table <- DT::renderDataTable({
      if(!input$import_from_for_reactome) {
      shiny::validate(
      need(!class(gene_df_reactome_check) == "try-error", message = "No genes meet your requirements, and can not do the reactome analysis")
    )
    }

    if(input$import_from_for_reactome) {
      shiny::validate(
      need(nrow(gene_df_reactome) != 0, message = "No genes meet your requirements, and can not do the reactome analysis")
    )
    }
    validate( need(res_reactome(), "No genes can be mapped to ENTREZID, and can not do the reactome analysis"))
      DT::datatable(res_reactome()$sig_table, filter = 'top', options = list( autoWidth = F,scrollX = TRUE))
    })

    output$significantBox_for_reactome <- renderInfoBox({
      if(!input$import_from_for_reactome) {
      shiny::validate(
      need(!class(gene_df_reactome_check) == "try-error", message = "No genes meet your requirements, and can not do the reactome analysis")
    )
    }

    if(input$import_from_for_reactome) {
      shiny::validate(
      need(nrow(gene_df_reactome) != 0, message = "No genes meet your requirements, and can not do the reactome analysis")
    )
    }
    validate( need(res_reactome(), "No genes can be mapped to ENTREZID, and can not do the reactome analysis"))

      num_total_reactome <- res_reactome()$all_table %>%
        nrow()
      num_signif_reactome <- res_reactome()$sig_table %>%
        nrow()
      frac_reactome <- num_signif_reactome / num_total_reactome

        if(frac_reactome == 0) {
        info_box_reactome <- infoBox("Significant terms",
                            paste0(num_signif_reactome,
                                   " out of ",
                                   num_total_reactome),
                            "No terms enriched",
                            icon = icon("thumbs-down", lib = "glyphicon"),
                            color = "red",
                            width = 4)
       }
      if(!frac_reactome == 0) {
        info_box_reactome <-     infoBox("Significant terms",
                              paste0(num_signif_reactome,
                                     " out of ",
                                     num_total_reactome),
                              paste0(signif(frac_reactome * 100, digits = 3),
                                     "% of terms enriched"),
                              icon = icon("thumbs-up", lib = "glyphicon"),
                              color = "green",
                              width = 4)
      }
      info_box_reactome
    })

    reactome_barplot_input <- reactive({
      # # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        barplot(res_reactome()$sig_res, showCategory = input$reactome_ShowCategory_bar, color = input$reactome_color)
    })

    reactome_dotplot_input <- reactive({
      progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        dotplot(res_reactome()$sig_res, showCategory = input$reactome_ShowCategory_dot, color = input$reactome_color)
      })

    reactome_dotplot_opt_input <- reactive({
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        my_dotplot_opt(res = res_reactome(), color = input$reactome_color, size = "Count", title = "", decreasing = TRUE, ShowCategory = input$reactome_ShowCategory_dot_opt)
      })

    reactome_heatplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      heatplot_for_react_kegg(res = res_reactome(), ShowCategory = input$reactome_ShowCategory_heat, df_with_lg2fc = input$df_with_lg2fc_for_reactome)
        })
      })

    reactome_cnetplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      cnetplot_for_react_kegg(res = res_reactome(), ShowCategory = input$reactome_ShowCategory_cnet, circular = input$reactome_circular_cnet, colorEdge = TRUE, df_with_lg2fc = input$df_with_lg2fc_for_reactome)
        })
      })

    reactome_emaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      emaplot_for_react_kegg(res = res_reactome(), ShowCategory = input$reactome_ShowCategory_ema, color = input$reactome_color, layout = "kk")
        })
      })

    output$reactome_barplot <- renderPlot({
      reactome_barplot_input()
    })

    output$reactome_dotplot <- renderPlot({
      reactome_dotplot_input()
    })

    output$reactome_dotplot_opt <- renderPlot({
      reactome_dotplot_opt_input()
    })

    output$reactome_heatplot <- renderPlot({
      reactome_heatplot_input()
    })

     output$reactome_cnetplot <- renderPlot({
      reactome_cnetplot_input()
    })

    output$reactome_emaplot <- renderPlot({
      reactome_emaplot_input()
    })

    ### Download objects and functions ### ------------------------------------
    datasetInput_for_reactome <- reactive({
      table_for_reactome = res_reactome()
          switch(input$dataset_for_reactome,
             "full_results" = res_reactome()$all_table,
             "significant_results" = res_reactome()$sig_table)
    })

    output$downloadreactome <- downloadHandler(
      filename = function() { paste(input$dataset_for_reactome, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput_for_reactome(),
        file,
        col.names = TRUE,
        row.names = FALSE,
        sep ="\t") }
    )

    output$download_reactome_barplot <- downloadHandler(
      filename = 'barplot.pdf',
      content = function(file) {
        pdf(file, width = input$reactome_wide_bar, height = input$reactome_high_bar)
        print(reactome_barplot_input())
        dev.off()
      }
    )

    output$download_reactome_dotplot <- downloadHandler(
      filename = 'dotplot.pdf',
      content = function(file) {
        pdf(file, width = input$reactome_wide_dot, height = input$reactome_high_dot)
        print(reactome_dotplot_input())
        dev.off()
      }
    )

    output$download_reactome_dotplot_opt <- downloadHandler(
      filename = 'dotplot_opt.pdf',
      content = function(file) {
        pdf(file, width = input$reactome_wide_dot_opt, height = input$reactome_high_dot_opt)
        print(reactome_dotplot_opt_input())
        dev.off()
      }
    )

    output$download_reactome_heatplot <- downloadHandler(
      filename = 'heatplot.pdf',
      content = function(file) {
        pdf(file, width = input$reactome_wide_heat, height = input$reactome_high_heat)
        print(reactome_heatplot_input())
        dev.off()
      }
    )

    output$download_reactome_cnetplot <- downloadHandler(
      filename = 'cnetplot.pdf',
      content = function(file) {
        pdf(file, width = input$reactome_wide_cnet, height = input$reactome_high_cnet)
        print(reactome_cnetplot_input())
        dev.off()
      }
    )

    output$download_reactome_emaplot <- downloadHandler(
      filename = 'emaplot.pdf',
      content = function(file) {
        pdf(file, width = input$reactome_wide_ema, height = input$reactome_high_ema)
        print(reactome_emaplot_input())
        dev.off()
      }
    )
  })

  ## for gseGO analysis
output$import_for_gsego <- renderUI({
    if(input$import_from_for_gsego){
                 radioButtons("import_for_gsego",
                              "Choose genes",
                              choices = c("DEP-LFQ", "DEG-RNAseq"),
                              selected = "DEP-LFQ")
    }
    })

output$import_contrast_for_gsego <- renderUI({
  if(input$import_from_for_gsego) {
if(input$import_for_gsego == "DEP-LFQ") {
  validate(need(!is.null(data()), "Please go to the DEP-LFQ options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_gsego", "Contrast", choices = c(get_results(dep()) %>%
  dplyr::filter(significant) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique()), selected = NULL)
    } else {
      validate(need(!is.null(countData()), "Please go to the DEG-RNAseq options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_gsego", "Contrast", choices = c(as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique()), selected = NULL)
    }
  } else {
    return(NULL)
  }

  })

output$text_input_for_gsego <- renderUI({
  if(!input$import_from_for_gsego){
    textAreaInput(inputId = "text_input_for_gsego", label = "Please paste your gene list", placeholder = "TP53  6.21\nPTEN  -1.53", rows = 8, width = "100%")
  }
  })

output$organism_for_gsego <- renderUI({
  # if(!input$import_from_for_gsego){
    selectizeInput("organism_for_gsego", "Select organism", choices=c("human","mouse"), selected="human")
  # }
  })

  # observeEvent for gsego analysis input
  observeEvent(input$analyze_for_gsego,{
      ### Interactive UI functions ### ------------------------------------------
    output$downloadTable_gsego <- renderUI({
      selectizeInput("dataset_for_gsego",
                     "Choose a dataset to save" ,
                     c("full_results","significant_results"
                       ))
    })

    output$downloadButton_gsego <- renderUI({
      downloadButton('downloadgsego', 'Save table', class = "downloadgsego")
    })

    ## reactive functions
    if(!input$import_from_for_gsego) {
      genelist <- reactive({ strsplit(input$text_input_for_gsego,'\n')[[1]] })
      gene_name <-reactive(unlist(strsplit(genelist(),";")))

      gene_df_gsego <- reactive({
        gene_df <<- data.frame(name=gene_name())
        gene_df$name <-  sapply(gene_df$name,the1stname)
        gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
        rownames(gene_df) = NULL
        colnames(gene_df) = c("name","fc")
        gene_df = as.data.frame(gene_df)
        gene_df$fc = as.numeric(as.character(gene_df$fc))
        gene_df$name = as.character(gene_df$name)
        gene_df
        })
    }

    if(input$import_from_for_gsego) {
      if(input$import_for_gsego == "DEP-LFQ") {
        gene_df_gsego <<- get_results(dep()) %>% dplyr::select(name, paste(input$import_contrast_for_gsego, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsego, "ratio", sep = "_"))
        } else {
          if(input$import_for_gsego == "DEG-RNAseq") {
            gene_df_gsego <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::select(name, paste(input$import_contrast_for_gsego, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsego, "log2FoldChange", sep = "_"))
          }
        }
    }
      # test <<- gene_df_gsego()

      # reat_gsego <<- reactive({
        #Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Please wait ...", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        reat_gsego <<- try(gsegoAnalysis(df = if(!input$import_from_for_gsego) {gene_df_gsego()} else {gene_df_gsego}, organism = input$organism_for_gsego), silent = TRUE)
        # })

     # test1 <<- reat_gsego()
    res_gsego <- reactive({
    #   validate(
    #   need(length(input$gsego_Phenotype) != 0, "Please choose at least one Phenotype"),
    #   need(class(reat_gsego()) != "try-error", "Your order ranked geneList can not do gsea analysis, Please go to GO panel if desired")
    # )
      try(givegseGO_res_and_table(reat = reat_gsego, ont = input$gsego_ont, pCutoff = input$gsego_p, p.adj.cutoff = input$gsego_padj, NES.cutoff = input$gsego_NES, simplify = input$gsego_simplify, Phenotype = input$gsego_Phenotype), silent = TRUE)
      })

    test2 <<- res_gsego()

    output$gsego_Table <- DT::renderDataTable({
    if(input$import_from_for_gsego) {
      shiny::validate(
      need(nrow(gene_df_gsego) != 0, message = "No genes meet your requirements, and can not do the gseGO analysis")
    )
    }
      validate(
      need(length(input$gsego_Phenotype) != 0, "Please choose at least one Phenotype"),
      need(reat_gsego, "Sorry, I can not do the gseGO analysis")
        )
      DT::datatable(res_gsego()$sig_table, filter = 'top', options = list(autoWidth = F,scrollX = TRUE))
      })#selection = list(selected = c(1))


    output$significantBox_for_gsego <- renderInfoBox({
    if(input$import_from_for_gsego) {
      shiny::validate(
      need(nrow(gene_df_gsego) != 0, message = "No genes meet your requirements, and can not do the gseGO analysis")
    )
    }
      validate(
      need(length(input$gsego_Phenotype) != 0, "Please choose at least one Phenotype"),
      need(reat_gsego, "Sorry, I can not do the gseGO analysis")
        )

      num_total_gsego <- res_gsego()$all_table %>%
        nrow()
      num_signif_gsego <- res_gsego()$sig_table %>%
        nrow()
      frac_gsego <- num_signif_gsego / num_total_gsego

        if(frac_gsego == 0) {
        info_box_gsego <- infoBox("Significant terms",
                            paste0(num_signif_gsego,
                                   " out of ",
                                   num_total_gsego),
                            "No terms enriched",
                            icon = icon("thumbs-down", lib = "glyphicon"),
                            color = "red",
                            width = 4)
       }
      if(!frac_gsego == 0) {
        info_box_gsego <-     infoBox("Significant terms",
                              paste0(num_signif_gsego,
                                     " out of ",
                                     num_total_gsego),
                              paste0(signif(frac_gsego * 100, digits = 3),
                                     "% of terms enriched"),
                              icon = icon("thumbs-up", lib = "glyphicon"),
                              color = "green",
                              width = 4)
      }
      info_box_gsego
    })

    output$gsego_term <- renderUI({
      if (!nrow(res_gsego()$sig_table) == 0) {
        selectizeInput("gsego_term",
                       "Term",
                       choices = res_gsego()$sig_table$Description, selected = res_gsego()$sig_table$Description[1], multiple = FALSE)
      }
    })

    gsego_barplot_input <- reactive({
      # # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        gse_barplot(res = res_gsego(), ShowCategory = input$gsego_ShowCategory_bar, color = input$gsego_color, ont = input$gsego_ont, Split = input$gsego_bar_if_Split_for_ont_ALL)
    })

    gsego_dotplot_input <- reactive({
      progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        gse_dotplot(res = res_gsego(), ShowCategory = input$gsego_ShowCategory_dot, color = input$gsego_color, ont = input$gsego_ont, Split = input$gsego_dot_if_Split_for_ont_ALL)
      })

    gsego_heatplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      try(Heatplot(res_gsego()$sig_res, showCategory = input$gsego_ShowCategory_heat, foldChange = res_gsego()$de), silent = T)
        })
      })

    gsego_cnetplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      try(Cnetplot(x = res_gsego()$sig_res, showCategory = input$gsego_ShowCategory_cnet, foldChange = res_gsego()$de, circular = input$gsego_circular_cnet, colorEdge = TRUE), silent = T)
        })
      })

    gsego_emaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      try(Emapplot(res_gsego()$sig_res, showCategory = input$gsego_ShowCategory_ema, color = input$gsego_color, layout = "kk"), silent = T)
        })
      })

    gsego_Gseaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      my_gseaplot2(res_gsego()$sig_res, geneSetID = match(input$gsego_term, res_gsego()$sig_res$Description), title = input$gsego_term,# geneSetID = input$gsego_Table_rows_selected, title = res_gsego()$sig_table[input$gsego_Table_rows_selected,3]
                   color = "green",
                   base_size = 11,
                   rel_heights = c(1.5, 0.5, 1),
                   subplots = 1:3,
                   pvalue_table = TRUE,
                   ES_geom = "line")
        })
      })

    output$gsego_barplot <- renderPlot({
      gsego_barplot_input()
    })

    output$gsego_dotplot <- renderPlot({
      gsego_dotplot_input()
    })

    output$gsego_heatplot <- renderPlot({
      gsego_heatplot_input()
    })

     output$gsego_cnetplot <- renderPlot({
      gsego_cnetplot_input()
    })

    output$gsego_emaplot <- renderPlot({
      gsego_emaplot_input()
    })

    output$gsego_Gseaplot <- renderPlot({
      gsego_Gseaplot_input()
    })

    ### Download objects and functions ### ------------------------------------
    datasetInput_for_gsego <- reactive({
      table_for_gsego = res_gsego()
      if(input$gsego_ont == "ALL") {
        if(input$gsego_bar_if_Split_for_ont_ALL | input$gsego_dot_if_Split_for_ont_ALL){
          switch(input$dataset_for_gsego,
             "full_results" = res_gsego()$all_table %>% dplyr::arrange(ONTOLOGY),
             "significant_results" = res_gsego()$sig_table %>% dplyr::arrange(ONTOLOGY))

          } else {
            switch(input$dataset_for_gsego,
             "full_results" = res_gsego()$all_table,
             "significant_results" = res_gsego()$sig_table)

          }

        } else {
          switch(input$dataset_for_gsego,
             "full_results" = res_gsego()$all_table,
             "significant_results" = res_gsego()$sig_table)
        }
      # switch(input$dataset_for_gsego,
      #        "full_results" = res_gsego()$all_table,
      #        "significant_results" = res_gsego()$sig_table)
    })

    output$downloadgsego <- downloadHandler(
      filename = function() { paste(input$dataset_for_gsego, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput_for_gsego(),
        file,
        col.names = TRUE,
        row.names = FALSE,
        sep ="\t") }
    )

    output$download_gsego_barplot <- downloadHandler(
      filename = 'barplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsego_wide_bar, height = input$gsego_high_bar)
        print(gsego_barplot_input())
        dev.off()
      }
    )

    output$download_gsego_dotplot <- downloadHandler(
      filename = 'dotplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsego_wide_dot, height = input$gsego_high_dot)
        print(gsego_dotplot_input())
        dev.off()
      }
    )

    output$download_gsego_heatplot <- downloadHandler(
      filename = 'heatplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsego_wide_heat, height = input$gsego_high_heat)
        print(gsego_heatplot_input())
        dev.off()
      }
    )

    output$download_gsego_cnetplot <- downloadHandler(
      filename = 'cnetplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsego_wide_cnet, height = input$gsego_high_cnet)
        print(gsego_cnetplot_input())
        dev.off()
      }
    )

    output$download_gsego_emaplot <- downloadHandler(
      filename = 'emaplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsego_wide_ema, height = input$gsego_high_ema)
        print(gsego_emaplot_input())
        dev.off()
      }
    )

    output$download_gsego_Gseaplot <- downloadHandler(
      filename = 'Gseaplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsego_wide_Gsea, height = input$gsego_high_Gsea)
        print(gsego_Gseaplot_input())
        dev.off()
      }
    )
  })

  ## for gseKEGG analysis
output$import_for_gsekegg <- renderUI({
    if(input$import_from_for_gsekegg){
                 radioButtons("import_for_gsekegg",
                              "Choose genes",
                              choices = c("DEP-LFQ", "DEG-RNAseq"),
                              selected = "DEP-LFQ")
    }
    })

output$import_contrast_for_gsekegg <- renderUI({
  if(input$import_from_for_gsekegg) {
if(input$import_for_gsekegg == "DEP-LFQ") {
  validate(need(!is.null(data()), "Please go to the DEP-LFQ options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_gsekegg", "Contrast", choices = c(get_results(dep()) %>%
  dplyr::filter(significant) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique()), selected = NULL)
    } else {
      validate(need(!is.null(countData()), "Please go to the DEG-RNAseq options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_gsekegg", "Contrast", choices = c(as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique()), selected = NULL)
    }
  } else {
    return(NULL)
  }

  })

output$text_input_for_gsekegg <- renderUI({
  if(!input$import_from_for_gsekegg){
    textAreaInput(inputId = "text_input_for_gsekegg", label = "Please paste your gene list", placeholder = "TP53  6.21\nPTEN  -1.53", rows = 8, width = "100%")
  }
  })

output$organism_for_gsekegg <- renderUI({
  # if(!input$import_from_for_gsekegg){
    selectizeInput("organism_for_gsekegg", "Select organism", choices=c("hsa","mmu"), selected="hsa")
  # }
  })

  # observeEvent for gsekegg analysis input
  observeEvent(input$analyze_for_gsekegg,{
      ### Interactive UI functions ### ------------------------------------------
    output$downloadTable_gsekegg <- renderUI({
      selectizeInput("dataset_for_gsekegg",
                     "Choose a dataset to save" ,
                     c("full_results","significant_results"
                       ))
    })

    output$downloadButton_gsekegg <- renderUI({
      downloadButton('downloadgsekegg', 'Save table', class = "downloadgsekegg")
    })

    ## reactive functions
    if(!input$import_from_for_gsekegg) {
      genelist <- reactive({ strsplit(input$text_input_for_gsekegg,'\n')[[1]] })
      gene_name <-reactive(unlist(strsplit(genelist(),";")))

      gene_df_gsekegg <- reactive({
        gene_df <<- data.frame(name=gene_name())
        gene_df$name <-  sapply(gene_df$name,the1stname)
        gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
        rownames(gene_df) = NULL
        colnames(gene_df) = c("name","fc")
        gene_df = as.data.frame(gene_df)
        gene_df$fc = as.numeric(as.character(gene_df$fc))
        gene_df$name = as.character(gene_df$name)
        gene_df
        })
    }

    if(input$import_from_for_gsekegg) {
      if(input$import_for_gsekegg == "DEP-LFQ") {
        gene_df_gsekegg <<- get_results(dep()) %>% dplyr::select(name, paste(input$import_contrast_for_gsekegg, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsekegg, "ratio", sep = "_"))
        } else {
          if(input$import_for_gsekegg == "DEG-RNAseq") {
            gene_df_gsekegg <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::select(name, paste(input$import_contrast_for_gsekegg, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsekegg, "log2FoldChange", sep = "_"))
          }
        }
    }

      # reat_gsekegg <<- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Please wait ...", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        reat_gsekegg <- try(gsekeggAnalysis(df = if(!input$import_from_for_gsekegg) {gene_df_gsekegg()} else {gene_df_gsekegg}, organism = input$organism_for_gsekegg), silent = TRUE)
        # })

    res_gsekegg <<- reactive({
      try(give_gsekegg_gsereat_res_and_table(reat = reat_gsekegg, pCutoff = input$gsekegg_p, p.adj.cutoff = input$gsekegg_padj, NES.cutoff = input$gsekegg_NES, Phenotype = input$gsekegg_Phenotype), silent = TRUE)
      })

    test2 <<- res_gsekegg()

    output$gsekegg_Table <- DT::renderDataTable({
    if(input$import_from_for_gsekegg) {
      shiny::validate(
      need(nrow(gene_df_gsekegg) != 0, message = "No genes meet your requirements, and can not do the gsekegg analysis")
    )
    }
      validate(
      need(length(input$gsekegg_Phenotype) != 0, "Please choose at least one Phenotype"),
      need(reat_gsekegg, "Sorry, I can not do the gsekegg analysis")
        )
      DT::datatable(res_gsekegg()$sig_table, filter = 'top', options = list(autoWidth = F,scrollX = TRUE))
      })#selection = list(selected = c(1))

    output$significantBox_for_gsekegg <- renderInfoBox({
    if(input$import_from_for_gsekegg) {
      shiny::validate(
      need(nrow(gene_df_gsekegg) != 0, message = "No genes meet your requirements, and can not do the gsekegg analysis")
    )
    }
      validate(
      need(length(input$gsekegg_Phenotype) != 0, "Please choose at least one Phenotype"),
      need(reat_gsekegg, "Sorry, I can not do the gsekegg analysis")
        )
      num_total_gsekegg <- res_gsekegg()$all_table %>%
        nrow()
      num_signif_gsekegg <- res_gsekegg()$sig_table %>%
        nrow()
      frac_gsekegg <- num_signif_gsekegg / num_total_gsekegg

        if(frac_gsekegg == 0) {
        info_box_gsekegg <- infoBox("Significant terms",
                            paste0(num_signif_gsekegg,
                                   " out of ",
                                   num_total_gsekegg),
                            "No terms enriched",
                            icon = icon("thumbs-down", lib = "glyphicon"),
                            color = "red",
                            width = 4)
       }
      if(!frac_gsekegg == 0) {
        info_box_gsekegg <-     infoBox("Significant terms",
                              paste0(num_signif_gsekegg,
                                     " out of ",
                                     num_total_gsekegg),
                              paste0(signif(frac_gsekegg * 100, digits = 3),
                                     "% of terms enriched"),
                              icon = icon("thumbs-up", lib = "glyphicon"),
                              color = "green",
                              width = 4)
      }
      info_box_gsekegg
    })

    output$gsekegg_term <- renderUI({
      if (!nrow(res_gsekegg()$sig_table) == 0) {
        selectizeInput("gsekegg_term",
                       "Term",
                       choices = res_gsekegg()$sig_table$Description, selected = res_gsekegg()$sig_table$Description[1], multiple = FALSE)
      }
    })

    gsekegg_barplot_input <- reactive({
      # # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        enrichplot:::barplot.enrichResult(res_gsekegg()$sig_res, showCategory = input$gsekegg_ShowCategory_bar, color = input$gsekegg_color, x = "NES", split="phenotype") + labs(y = "NES")
    })

    gsekegg_dotplot_input <- reactive({
      progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        enrichplot:::dotplot(res_gsekegg()$sig_res, showCategory = input$gsekegg_ShowCategory_dot, color = input$gsekegg_color, x = "NES", split="phenotype")
      })

    gsekegg_heatplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      try(Heatplot(res_gsekegg()$sig_res, showCategory = input$gsekegg_ShowCategory_heat, foldChange = res_gsekegg()$de), silent = T)
        })
      })

    gsekegg_cnetplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      try(Cnetplot(x = res_gsekegg()$sig_res, showCategory = input$gsekegg_ShowCategory_cnet, foldChange = res_gsekegg()$de, circular = input$gsekegg_circular_cnet, colorEdge = TRUE), silent = T)
        })
      })

    gsekegg_emaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      try(Emapplot(res_gsekegg()$sig_res, showCategory = input$gsekegg_ShowCategory_ema, color = input$gsekegg_color, layout = "kk"), silent = T)
        })
      })

    gsekegg_Gseaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      my_gseaplot2(res_gsekegg()$sig_res, geneSetID = match(input$gsekegg_term, res_gsekegg()$sig_res$Description), title = input$gsekegg_term,# geneSetID = input$gsekegg_Table_rows_selected, title = res_gsekegg()$sig_table[input$gsekegg_Table_rows_selected,3]
                   color = "green",
                   base_size = 11,
                   rel_heights = c(1.5, 0.5, 1),
                   subplots = 1:3,
                   pvalue_table = TRUE,
                   ES_geom = "line")
        })
      })

    output$gsekegg_barplot <- renderPlot({
      gsekegg_barplot_input()
    })

    output$gsekegg_dotplot <- renderPlot({
      gsekegg_dotplot_input()
    })

    output$gsekegg_heatplot <- renderPlot({
      gsekegg_heatplot_input()
    })

     output$gsekegg_cnetplot <- renderPlot({
      gsekegg_cnetplot_input()
    })

    output$gsekegg_emaplot <- renderPlot({
      gsekegg_emaplot_input()
    })

    output$gsekegg_Gseaplot <- renderPlot({
      gsekegg_Gseaplot_input()
    })

    ### Download objects and functions ### ------------------------------------
    datasetInput_for_gsekegg <- reactive({
      table_for_gsekegg = res_gsekegg()
          switch(input$dataset_for_gsekegg,
             "full_results" = res_gsekegg()$all_table,
             "significant_results" = res_gsekegg()$sig_table)
    })

    output$downloadgsekegg <- downloadHandler(
      filename = function() { paste(input$dataset_for_gsekegg, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput_for_gsekegg(),
        file,
        col.names = TRUE,
        row.names = FALSE,
        sep ="\t") }
    )

    output$download_gsekegg_barplot <- downloadHandler(
      filename = 'barplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_bar, height = input$gsekegg_high_bar)
        print(gsekegg_barplot_input())
        dev.off()
      }
    )

    output$download_gsekegg_dotplot <- downloadHandler(
      filename = 'dotplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_dot, height = input$gsekegg_high_dot)
        print(gsekegg_dotplot_input())
        dev.off()
      }
    )

    output$download_gsekegg_heatplot <- downloadHandler(
      filename = 'heatplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_heat, height = input$gsekegg_high_heat)
        print(gsekegg_heatplot_input())
        dev.off()
      }
    )

    output$download_gsekegg_cnetplot <- downloadHandler(
      filename = 'cnetplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_cnet, height = input$gsekegg_high_cnet)
        print(gsekegg_cnetplot_input())
        dev.off()
      }
    )

    output$download_gsekegg_emaplot <- downloadHandler(
      filename = 'emaplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_ema, height = input$gsekegg_high_ema)
        print(gsekegg_emaplot_input())
        dev.off()
      }
    )

    output$download_gsekegg_Gseaplot <- downloadHandler(
      filename = 'Gseaplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsekegg_wide_Gsea, height = input$gsekegg_high_Gsea)
        print(gsekegg_Gseaplot_input())
        dev.off()
      }
    )
  })

  ## for gseReactome analysis
output$import_for_gsereactome <- renderUI({
    if(input$import_from_for_gsereactome){
                 radioButtons("import_for_gsereactome",
                              "Choose genes",
                              choices = c("DEP-LFQ", "DEG-RNAseq"),
                              selected = "DEP-LFQ")
    }
    })

output$import_contrast_for_gsereactome <- renderUI({
  if(input$import_from_for_gsereactome) {
if(input$import_for_gsereactome == "DEP-LFQ") {
  validate(need(!is.null(data()), "Please go to the DEP-LFQ options, and the differential analysis first"))
      selectizeInput("import_contrast_for_gsereactome", "Contrast", choices = c(get_results(dep()) %>%
  dplyr::filter(significant) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique()), selected = NULL)
    } else {
      validate(need(!is.null(countData()), "Please go to the DEG-RNAseq options, and the differential analysis first"))
      selectizeInput("import_contrast_for_gsereactome", "Contrast", choices = c(as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique()), selected = NULL)
    }
  } else {
    return(NULL)
  }

  })

output$text_input_for_gsereactome <- renderUI({
  if(!input$import_from_for_gsereactome){
    textAreaInput(inputId = "text_input_for_gsereactome", label = "Please paste your gene list", placeholder = "TP53  6.21\nPTEN  -1.53", rows = 8, width = "100%")
  }
  })

output$organism_for_gsereactome <- renderUI({
  # if(!input$import_from_for_gsereactome){
    selectizeInput("organism_for_gsereactome", "Select organism", choices=c("human","mouse"), selected="human")
  # }
  })

  # observeEvent for gsereactome analysis input
  observeEvent(input$analyze_for_gsereactome,{
      ### Interactive UI functions ### ------------------------------------------
    output$downloadTable_gsereactome <- renderUI({
      selectizeInput("dataset_for_gsereactome",
                     "Choose a dataset to save" ,
                     c("full_results","significant_results"
                       ))
    })

    output$downloadButton_gsereactome <- renderUI({
      downloadButton('downloadgsereactome', 'Save table', class = "downloadgsereactome")
    })

    ## reactive functions
    if(!input$import_from_for_gsereactome) {
      genelist <- reactive({ strsplit(input$text_input_for_gsereactome,'\n')[[1]] })
      gene_name <-reactive(unlist(strsplit(genelist(),";")))

      gene_df_gsereactome <- reactive({
        gene_df <<- data.frame(name=gene_name())
        gene_df$name <-  sapply(gene_df$name,the1stname)
        gene_df = t(as.data.frame(sapply(gene_df$name, function(i){strsplit(i, split = "\t")})))
        rownames(gene_df) = NULL
        colnames(gene_df) = c("name","fc")
        gene_df = as.data.frame(gene_df)
        gene_df$fc = as.numeric(as.character(gene_df$fc))
        gene_df$name = as.character(gene_df$name)
        gene_df
        })
    }

    if(input$import_from_for_gsereactome) {
      if(input$import_for_gsereactome == "DEP-LFQ") {
        gene_df_gsereactome <<- get_results(dep()) %>% dplyr::select(name, paste(input$import_contrast_for_gsereactome, "ratio", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsereactome, "ratio", sep = "_"))
        } else {
          if(input$import_for_gsereactome == "DEG-RNAseq") {
            gene_df_gsereactome <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::select(name, paste(input$import_contrast_for_gsereactome, "log2FoldChange", sep = "_")) %>% dplyr::rename(fc = paste(input$import_contrast_for_gsereactome, "log2FoldChange", sep = "_"))
          }
        }
    }

      # reat_gsereactome <<- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Please wait ...", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        reat_gsereactome <- try(gsereactAnalysis(df = if(!input$import_from_for_gsereactome) {gene_df_gsereactome()} else {gene_df_gsereactome}, organism = input$organism_for_gsereactome), silent = TRUE)
        # })

    res_gsereactome <<- reactive({
      try(give_gsekegg_gsereat_res_and_table(reat = reat_gsereactome, pCutoff = input$gsereactome_p, p.adj.cutoff = input$gsereactome_padj, NES.cutoff = input$gsereactome_NES, Phenotype = input$gsereactome_Phenotype), silent = TRUE)
      })

    test2 <<- res_gsereactome()

    output$gsereactome_Table <- DT::renderDataTable({
    if(input$import_from_for_gsereactome) {
      shiny::validate(
      need(nrow(gene_df_gsereactome) != 0, message = "No genes meet your requirements, and can not do the gseReactome analysis")
    )
    }
      validate(
      need(length(input$gsereactome_Phenotype) != 0, "Please choose at least one Phenotype"),
      need(reat_gsereactome, "Sorry, I can not do the gsereactome analysis")
        )
      DT::datatable(res_gsereactome()$sig_table, filter = 'top', options = list(autoWidth = F,scrollX = TRUE))
      })#selection = list(selected = c(1))

    output$significantBox_for_gsereactome <- renderInfoBox({
    if(input$import_from_for_gsereactome) {
      shiny::validate(
      need(nrow(gene_df_gsereactome) != 0, message = "No genes meet your requirements, and can not do the gseReactome analysis")
    )
    }
      validate(
      need(length(input$gsereactome_Phenotype) != 0, "Please choose at least one Phenotype"),
      need(reat_gsereactome, "Sorry, I can not do the gsereactome analysis")
        )
      num_total_gsereactome <- res_gsereactome()$all_table %>%
        nrow()
      num_signif_gsereactome <- res_gsereactome()$sig_table %>%
        nrow()
      frac_gsereactome <- num_signif_gsereactome / num_total_gsereactome

        if(frac_gsereactome == 0) {
        info_box_gsereactome <- infoBox("Significant terms",
                            paste0(num_signif_gsereactome,
                                   " out of ",
                                   num_total_gsereactome),
                            "No terms enriched",
                            icon = icon("thumbs-down", lib = "glyphicon"),
                            color = "red",
                            width = 4)
       }
      if(!frac_gsereactome == 0) {
        info_box_gsereactome <-     infoBox("Significant terms",
                              paste0(num_signif_gsereactome,
                                     " out of ",
                                     num_total_gsereactome),
                              paste0(signif(frac_gsereactome * 100, digits = 3),
                                     "% of terms enriched"),
                              icon = icon("thumbs-up", lib = "glyphicon"),
                              color = "green",
                              width = 4)
      }
      info_box_gsereactome
    })

    output$gsereactome_term <- renderUI({
      if (!nrow(res_gsereactome()$sig_table) == 0) {
        selectizeInput("gsereactome_term",
                       "Term",
                       choices = res_gsereactome()$sig_table$Description, selected = res_gsereactome()$sig_table$Description[1], multiple = FALSE)
      }
    })

    gsereactome_barplot_input <- reactive({
      # # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        enrichplot:::barplot.enrichResult(res_gsereactome()$sig_res, showCategory = input$gsereactome_ShowCategory_bar, color = input$gsereactome_color, x = "NES", split="phenotype") + labs(y = "NES")
    })

    gsereactome_dotplot_input <- reactive({
      progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        enrichplot:::dotplot(res_gsereactome()$sig_res, showCategory = input$gsereactome_ShowCategory_dot, color = input$gsereactome_color, x = "NES", split="phenotype")
      })

    gsereactome_heatplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      try(Heatplot(res_gsereactome()$sig_res, showCategory = input$gsereactome_ShowCategory_heat, foldChange = res_gsereactome()$de), silent = T)
        })
      })

    gsereactome_cnetplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      try(Cnetplot(x = res_gsereactome()$sig_res, showCategory = input$gsereactome_ShowCategory_cnet, foldChange = res_gsereactome()$de, circular = input$gsereactome_circular_cnet, colorEdge = TRUE), silent = T)
        })
      })

    gsereactome_emaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      try(Emapplot(res_gsereactome()$sig_res, showCategory = input$gsereactome_ShowCategory_ema, color = input$gsereactome_color, layout = "kk"), silent = T)
        })
      })

    gsereactome_Gseaplot_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
      my_gseaplot2(res_gsereactome()$sig_res, geneSetID = match(input$gsereactome_term, res_gsereactome()$sig_res$Description), title = input$gsereactome_term,# geneSetID = input$gsereactome_Table_rows_selected, title = res_gsereactome()$sig_table[input$gsereactome_Table_rows_selected,3]
                   color = "green",
                   base_size = 11,
                   rel_heights = c(1.5, 0.5, 1),
                   subplots = 1:3,
                   pvalue_table = TRUE,
                   ES_geom = "line")
        })
      })

    output$gsereactome_barplot <- renderPlot({
      gsereactome_barplot_input()
    })

    output$gsereactome_dotplot <- renderPlot({
      gsereactome_dotplot_input()
    })

    output$gsereactome_heatplot <- renderPlot({
      gsereactome_heatplot_input()
    })

     output$gsereactome_cnetplot <- renderPlot({
      gsereactome_cnetplot_input()
    })

    output$gsereactome_emaplot <- renderPlot({
      gsereactome_emaplot_input()
    })

    output$gsereactome_Gseaplot <- renderPlot({
      gsereactome_Gseaplot_input()
    })

    ### Download objects and functions ### ------------------------------------
    datasetInput_for_gsereactome <- reactive({
      table_for_gsereactome = res_gsereactome()
          switch(input$dataset_for_gsereactome,
             "full_results" = res_gsereactome()$all_table,
             "significant_results" = res_gsereactome()$sig_table)
    })

    output$downloadgsereactome <- downloadHandler(
      filename = function() { paste(input$dataset_for_gsereactome, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput_for_gsereactome(),
        file,
        col.names = TRUE,
        row.names = FALSE,
        sep ="\t") }
    )

    output$download_gsereactome_barplot <- downloadHandler(
      filename = 'barplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_bar, height = input$gsereactome_high_bar)
        print(gsereactome_barplot_input())
        dev.off()
      }
    )

    output$download_gsereactome_dotplot <- downloadHandler(
      filename = 'dotplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_dot, height = input$gsereactome_high_dot)
        print(gsereactome_dotplot_input())
        dev.off()
      }
    )

    output$download_gsereactome_heatplot <- downloadHandler(
      filename = 'heatplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_heat, height = input$gsereactome_high_heat)
        print(gsereactome_heatplot_input())
        dev.off()
      }
    )

    output$download_gsereactome_cnetplot <- downloadHandler(
      filename = 'cnetplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_cnet, height = input$gsereactome_high_cnet)
        print(gsereactome_cnetplot_input())
        dev.off()
      }
    )

    output$download_gsereactome_emaplot <- downloadHandler(
      filename = 'emaplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_ema, height = input$gsereactome_high_ema)
        print(gsereactome_emaplot_input())
        dev.off()
      }
    )

    output$download_gsereactome_Gseaplot <- downloadHandler(
      filename = 'Gseaplot.pdf',
      content = function(file) {
        pdf(file, width = input$gsereactome_wide_Gsea, height = input$gsereactome_high_Gsea)
        print(gsereactome_Gseaplot_input())
        dev.off()
      }
    )
  })


    ## for PPITools
output$import_for_ppi <- renderUI({
    if(input$import_from_for_ppi){
                 radioButtons("import_for_ppi",
                              "Choose genes",
                              choices = c("UPregu for DEP-LFQ", "DOWNregu for DEP-LFQ", "UPDOWN for DEP-LFQ", "UPregu for DEG-RNAseq", "DOWNregu for DEG-RNAseq", "UPDOWN for DEG-RNAseq"),
                              selected = "UPregu for DEP-LFQ")
    }
    })

output$import_contrast_for_ppi <- renderUI({
  if(input$import_from_for_ppi) {
if(input$import_for_ppi == "UPregu for DEP-LFQ" | input$import_for_ppi == "DOWNregu for DEP-LFQ" | input$import_for_ppi == "UPDOWN for DEP-LFQ") {
  validate(need(!is.null(data()), "Please go to the DEP-LFQ options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_ppi", "Contrast", choices = c(get_results(dep()) %>%
  dplyr::filter(significant) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique(), "Any significant"), selected = NULL)
    } else {
      validate(need(!is.null(countData()), "Please go to the DEG-RNAseq options, and do the differential analysis first"))
      selectizeInput("import_contrast_for_ppi", "Contrast", choices = c(as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique(), "Any significant"), selected = NULL)
    }
  } else {
    return(NULL)
  }

  })

output$input_text <- renderUI({
  if(!input$import_from_for_ppi){
    textAreaInput(inputId = "input_text", label = "Please paste your gene list", placeholder = "TP53\nPTEN", rows = 8, width = "100%")
  }
  })
  #### observeRvent for PPITools ####
  observeEvent(input$String_annalysis,{
    withProgress(message = 'Please wait', value = 0, {
      memory.limit(2000000)#2GB
      output$download_PPITable <- renderUI({
      downloadButton("downloadPPITable", "Save table", class = "downloadPPITable")
        })

      output$download_PPInetwork <- renderUI({
      downloadButton("downloadPPInetwork", "Save network", class = "downloadPPInetwork")
        })

      if(!input$import_from_for_ppi) {
        genelist <- strsplit(input$input_text,'\n')[[1]]
        gene_name <<-data.frame(name = unlist(strsplit(genelist,";")[]))
      }
      # genelist1 <<- genelist()

    if(input$import_from_for_ppi) {
      if(input$import_for_ppi == "UPregu for DEP-LFQ" | input$import_for_ppi == "DOWNregu for DEP-LFQ" | input$import_for_ppi == "UPDOWN for DEP-LFQ") {
        if(input$import_contrast_for_ppi == "Any significant") {
          if(input$import_for_ppi == "UPDOWN for DEP-LFQ") {
            gene_name <<- get_results(dep()) %>% dplyr::filter(significant) %>% dplyr::select(name)
            } else {
              if(input$import_for_ppi == "UPregu for DEP-LFQ") {
                index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
                df <- get_results(dep()) %>% dplyr::filter(significant)
                cols_diff = df %>% dplyr::select(index)
                cols_diff_reject = cols_diff > 0
                cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                index1 = apply(cols_diff_reject, 1, all)
                gene_name <<- df[index1, ] %>% dplyr::select(name)
                } else {
                  if(input$import_for_ppi == "DOWNregu for DEP-LFQ") {
                    index <- get_results(dep()) %>% colnames(.) %>% grep("_ratio", ., value = T)
                    df <- get_results(dep()) %>% dplyr::filter(significant)
                    cols_diff = df %>% dplyr::select(index)
                    cols_diff_reject = cols_diff < 0
                    cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                    index1 = apply(cols_diff_reject, 1, all)
                    gene_name <<- df[index1, ] %>% dplyr::select(name)
                  }
                }
            }
          } else {
            if(input$import_for_ppi == "UPDOWN for DEP-LFQ") {
              gene_name <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_ppi, "significant", sep = "_"))) %>% dplyr::select(name)
              } else {
                if(input$import_for_ppi == "UPregu for DEP-LFQ") {
                  gene_name <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_ppi, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_ppi, "ratio", sep = "_")) > 0) %>% dplyr::select(name)
                  } else {
                    if(input$import_for_ppi == "DOWNregu for DEP-LFQ") {
                      gene_name <<- get_results(dep()) %>% dplyr::filter(get(paste(input$import_contrast_for_ppi, "significant", sep = "_"))) %>% filter(get(paste(input$import_contrast_for_ppi, "ratio", sep = "_")) < 0) %>% dplyr::select(name)
                    }
                  }
              }
          }
        } else {
          if(input$import_for_ppi == "UPregu for DEG-RNAseq" | input$import_for_ppi == "DOWNregu for DEG-RNAseq" | input$import_for_ppi == "UPDOWN for DEG-RNAseq") {
            if(input$import_contrast_for_ppi == "Any significant") {
              if(input$import_for_ppi == "UPDOWN for DEG-RNAseq") {
                gene_name <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant) %>% dplyr::select(name)
                } else {
                  if(input$import_for_ppi == "UPregu for DEG-RNAseq") {
                    index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
                    df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
                    cols_diff = df %>% dplyr::select(index)
                    cols_diff_reject = cols_diff > 0
                    cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                    index1 = apply(cols_diff_reject, 1, all)
                    gene_name <<- df[index1, ] %>% dplyr::select(name)
                    } else {
                      if(input$import_for_ppi == "DOWNregu for DEG-RNAseq") {
                        index <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% colnames(.) %>% grep("_log2FoldChange", ., value = T)
                        df <- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(significant)
                        cols_diff = df %>% dplyr::select(index)
                        cols_diff_reject = cols_diff < 0
                        cols_diff_reject[is.na(cols_diff_reject)] <- FALSE
                        index1 = apply(cols_diff_reject, 1, all)
                        gene_name <<- df[index1, ] %>% dplyr::select(name)
                      }
                    }
                }
              } else {
                if(input$import_for_ppi == "UPDOWN for DEG-RNAseq") {
                  gene_name <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_ppi, "significant", sep = "_"))) %>% dplyr::select(name)
                  } else {
                    if(input$import_for_ppi == "UPregu for DEG-RNAseq") {
                      gene_name <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_ppi, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_ppi, "log2FoldChange", sep = "_")) > 0) %>% dplyr::select(name)
                      } else {
                        if(input$import_for_ppi == "DOWNregu for DEG-RNAseq") {
                          gene_name <<- as.data.frame(my_res_for_RNAseq_4$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., name = if(input$transid_for_RNAseq) {"symbol"} else {"rowname"}) %>% dplyr::filter(get(paste(input$import_contrast_for_ppi, "significant", sep = "_"))) %>% dplyr::filter(get(paste(input$import_contrast_for_ppi, "log2FoldChange", sep = "_")) < 0) %>% dplyr::select(name)
                        }
                      }
                  }
              }
          }

        }
    }

    if(!nrow(gene_name) == 0) {
      gene_name <<- gene_name %>% dplyr::filter(!is.na(name))
    }

    if(nrow(gene_name) == 0) {
         output$reactome_plot <- renderVisNetwork({
          validate(need(nrow(gene_name) != 0, "No genes meet your requirements, and can not do the PPI(Protein-Protein Interaction Network) analysis"))
          a <- reactive({
            stringNetwork(linksTable = links2() , layoutway=input$layoutway , nodecolor= input$nodecolor , nodeshape=input$nodeshape, linecolor=input$linecolor , nodesize=input$nodesize , changesize=input$changesize , fontsize=input$PPIfontsize , changewidth=input$changewidth , linewidth=input$linewidth,smoothline = input$smoothline, smoothtype = input$smoothtype,highlightkey=input$highlightkey)
            })
          assign("aplot",a(),envir =  strEnv)
          a()
        })
     }
      validate(need(nrow(gene_name) != 0, "No genes meet your requirements, and can not do the PPI(Protein-Protein Interaction Network) analysis"))
      strEnv = get_string_Env()
      speciesname <- input$organism

      use_cache = F
      if(exists(c("load succed"),envir = strEnv) && get("speciesname",envir = strEnv) == speciesname && get("load succed",envir = strEnv)){
        use_cache = T
      }
      if(!use_cache) load_PPIdata(speciesname)


      orgDbname = switch(speciesname,
                         "human" = "org.Hs.eg.db",
                         "mouse" = "org.Mm.eg.db",
                         "Arabidopsis thaliana" = "org.At.tair.db",
                         "rat" = "org.Rn.eg.db",
                         "Saccharomyces cerevisiae" = "org.Sc.sgd.db",
                         "zebrafish" = "org.Dr.eg.db")
      ## checkinstall
      if (!requireNamespace(orgDbname, quietly = TRUE)) BiocManager::install(orgDbname)
      ## librar orgDb
      eval(parse(text = paste("library(",orgDbname,")",sep = "")))
      orgDb <- get(orgDbname)
      ## transform geneID

      # geneID <- gene_name$name %>% bitr(fromType = "SYMBOL",
      #                                   toType = "ENTREZID",
      #                                   OrgDb = orgDb,
      #                                   drop = T)
      geneID <- my_to_entrezid(orgDB = orgDb, gene = as.character(gene_name$name)) %>% tibble::rownames_to_column() %>% dplyr::rename(SYMBOL = rowname, ENTREZID = id) %>% dplyr::select(SYMBOL, ENTREZID) %>% dplyr::filter(!is.na(ENTREZID))
      if(nrow(geneID) == 0) {
         output$reactome_plot <- renderVisNetwork({
          validate(need(nrow(geneID) != 0, "No genes can be mapped to ENTREZID, and can not do the PPI(Protein-Protein Interaction Network) analysis"))
          a <- reactive({
            stringNetwork(linksTable = links2() , layoutway=input$layoutway , nodecolor= input$nodecolor , nodeshape=input$nodeshape, linecolor=input$linecolor , nodesize=input$nodesize , changesize=input$changesize , fontsize=input$PPIfontsize , changewidth=input$changewidth , linewidth=input$linewidth,smoothline = input$smoothline, smoothtype = input$smoothtype,highlightkey=input$highlightkey)
            })
          assign("aplot",a(),envir =  strEnv)
          a()
        })
      }
      validate(need(nrow(geneID) != 0, "No genes meet your requirements, and can not do the PPI(Protein-Protein Interaction Network) analysis"))

      aliasDf <- get("aliasDf",envir = strEnv)
      proteinsDf <- get("proteinsDf",envir = strEnv)
      protein_links_detail <- get("protein_links_detail",envir = strEnv)
      data_mapped <<- geneID %>% mymap(my_data_frame_id_col_names = "ENTREZID",
                                       aliasDf=get("aliasDf",envir = strEnv),
                                       proteinsDf=get("proteinsDf",envir = strEnv), ## data_mapped 3cols ENTREZID SYMBOL STRING_id
                                       removeUnmappedRows = TRUE,reload = F)
      incProgress(0.3)
      link_table <- reactive({myget_interactions(protein_detail = get("protein_links_detail",envir = strEnv),mapped_data = data_mapped,choose_scores = input$chooseScore, score_cutoff = 400)})
      links <- reactive({
        link_table() %>%                                                ## links 3cols from to weight from&to is SYMBOL
          mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>%
          mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"])
      })
      links2 <- reactive({filterlink(links(), scorecutoff=input$scorecutoff) })

      output$String_Table <- DT::renderDataTable(
        links2(), filter = 'top', options = list( autoWidth = F,scrollX = TRUE
        )
      )
      incProgress(0.5)
      # link3 <<- links2()
      observe({
        output$reactome_plot <- renderUI({
          visNetworkOutput("network",height=input$plothigh*10)
        })
        output$network <- renderVisNetwork({
          validate(need(nrow(links2()) != 0, "There is no Protein-Protein Interaction Network"))
          a <- reactive({
            stringNetwork(linksTable = links2() , layoutway=input$layoutway , nodecolor= input$nodecolor , nodeshape=input$nodeshape, linecolor=input$linecolor , nodesize=input$nodesize , changesize=input$changesize , fontsize=input$PPIfontsize , changewidth=input$changewidth , linewidth=input$linewidth,smoothline = input$smoothline, smoothtype = input$smoothtype,highlightkey=input$highlightkey)
            })
          assign("aplot",a(),envir =  strEnv)
          a()
        })
        ###output and download
        output$downloadPPITable <- downloadHandler(
          filename = "stringResult.txt",
          content = function(file) {
            write.table(links2(),
                        file,
                        col.names = TRUE,
                        row.names = FALSE,
                        sep ="\t") }
        )

        output$downloadPPInetwork <- downloadHandler(
          filename = "stringNetwork.html",
          content = function(file) {
            aplot <- get("aplot", envir = strEnv)
            htmlwidgets::saveWidget(aplot, file)}
        )
      })
      incProgress(0.2)
    })

  })
# # analysis for DEG-RNAseq
#
### UI functions output ui### --------------------------------------------------------
    output$dds_design <- renderUI({
      if(!is.null(input$file2_for_DERNAseq)) {
              selectizeInput("dds_design",
              label = "Select design",
              choices = c("", colnames(exdesign())), selected = "", multiple = TRUE)
             }
      })

    output$choose_expfac <- renderUI({
      if(!is.null(input$file2_for_DERNAseq)) {
                selectizeInput("choose_expfac",
              label = "Build the contrast upon",
              choices = c("", colnames(exdesign())), selected = "", multiple = FALSE)
              }
      })

    output$whether_LRT <- renderUI({
      shiny::validate(
        need(
          !is.null(input$file2_for_DERNAseq),
          "Please upload your countmatrix and exdesign files first"
        )
      )
      if(!is.null(input$file2_for_DERNAseq)){
      choose_factor_num <- exdesign() %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character() %>% unique() %>% length()
      }
      if(!is.null(input$file2_for_DERNAseq) & choose_factor_num > 2) {
        # radioButtons("Test_for_RNAseq",
        #                       "Test",
        #                       choices = c("Wald Test", "LRT Test"),
        #                       selected = "LRT Test")
        checkboxInput("perform_LRT_test", "perform a LRT test", value = TRUE)
      }
      })

    output$reduced_model <- renderUI({
      shiny::validate(
        need(
          !is.null(input$file2_for_DERNAseq),
          ""
        )
      )
      if(!is.null(input$file2_for_DERNAseq)){
      choose_factor_num <- exdesign() %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character() %>% unique() %>% length()
      }
      shiny::validate(
        need(
          input$perform_LRT_test,
          ""
        )
      )
      shiny::validate(
        need(
          input$dds_design != "",
          ""
        )
      )
      if(!is.null(input$file2_for_DERNAseq) & choose_factor_num > 2 & input$perform_LRT_test & input$dds_design != "") {
        selectizeInput("reduced_model",
              label = "Select the reduced model",
              choices = c("", "1", input$dds_design), selected = "", multiple = TRUE)
      }
      })

    # output$choose_expfac1 <- renderUI({
    #   shiny::validate(
    #     need(
    #       !is.null(input$file2_for_DERNAseq),
    #       ""#Please upload your countmatrix and exdesign files first
    #     )
    #   )
    #   shiny::validate(
    #     need(
    #       input$choose_expfac != "",
    #       "Please select the design and an experimental factor to build the contrast upon"
    #     )
    #   )
    #    if(!is.null(input$file2_for_DERNAseq)) {
    #             selectizeInput("choose_expfac1",
    #           label = "Numerator for the fold change",
    #           choices = exdesign() %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character(), selected = NULL, multiple = FALSE)
    #           }
    #   })

    # output$choose_expfac2 <- renderUI({
    #   if(!is.null(input$file2_for_DERNAseq)) {
    #             selectizeInput("choose_expfac2",
    #           label = "Denominator for the fold change",
    #           choices = exdesign() %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character(), selected = NULL, multiple = FALSE)
    #           }
    #   })

  output$control_for_RNAseq <- renderUI({
      shiny::validate(
        need(
          !is.null(input$file2_for_DERNAseq),
          ""#Please upload your countmatrix and exdesign files first
        )
      )
      shiny::validate(
        need(
          input$choose_expfac != "",
          "Please select the design and an experimental factor to build the contrast upon"
        )
      )
       if(!is.null(input$file2_for_DERNAseq)  & input$contrasts_for_RNAseq == "control") {
                selectizeInput("control_for_RNAseq",
              label = "Control",
              choices = exdesign() %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character(), selected = NULL, multiple = FALSE)
              }
    })

  output$test_manual_for_RNAseq <- renderUI({
      if(!is.null(input$file2_for_DERNAseq)  & input$contrasts_for_RNAseq == "manual") {
              test_manual_name <- exdesign() %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character() %>% unique()
              test_manual_name <- cbind(combn(test_manual_name,2),combn(test_manual_name,2, FUN = rev))
              selectizeInput("test_manual_for_RNAseq",
              label = "Manual test",
              choices = apply(test_manual_name, 2, function(i){paste(i[1], i[2], sep = "_vs_")}), selected = NULL, multiple = TRUE)
              }
    })

    output$ui_nrcores <- renderUI({
      mincores <- 1
      maxcores <- BiocParallel::multicoreWorkers()
      sliderInput("nrcores",
        label = "Cores to use",
        min = mincores, max = maxcores, value = 1, step = 1
      )
    })

  annoSpecies_df <-
      data.frame(
        species = c(
          "", "Anopheles", "Arabidopsis", "Bovine", "Worm",
          "Canine", "Fly", "Zebrafish", "E coli strain K12",
          "E coli strain Sakai", "Chicken", "Human", "Mouse",
          "Rhesus", "Malaria", "Chimp", "Rat",
          "Yeast", "Streptomyces coelicolor", "Pig", "Toxoplasma gondii",
          "Xenopus"
        ),
        pkg = c(
          "", "org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db",
          "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db",
          "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db",
          "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db", "org.Rn.eg.db",
          "org.Sc.sgd.db", "org.Sco.eg.db", "org.Ss.eg.db", "org.Tgondii.eg.db",
          "org.Xl.eg.db"
        ),
        stringsAsFactors = FALSE
      )

  annoSpecies_df <- annoSpecies_df[order(annoSpecies_df$species), ]
  rownames(annoSpecies_df) <- annoSpecies_df$species # easier to access afterwards

  output$ui_selectspecies <- renderUI({
      if (input$transid_for_RNAseq) {
        selectInput("speciesSelect",
        label = "Select species",
        choices = annoSpecies_df$species, selected = ""
      )
      } else {
        return(NULL)
      }
    })

  output$ui_set_no_map_to_rowname <- renderUI({
      if (input$transid_for_RNAseq) {
        checkboxInput("set_no_map_to_rowname", "No-mapping to rowname", value = TRUE)
      } else {
        return(NULL)
      }
    })

   output$ui_idtype <- renderUI({
      if (input$transid_for_RNAseq) {
        shiny::validate(
        need(
          input$speciesSelect != "",
          "Please select the species"
        )
      )
        std_choices <- c("ENSEMBL", "ENTREZID", "REFSEQ", "SYMBOL")
        annopkg <- annoSpecies_df$pkg[annoSpecies_df$species == input$speciesSelect]
        shiny::validate(
        need(
          require(annopkg, character.only = TRUE),
          paste0("The package ", annopkg, " is not installed/available. Try installing it with BiocManager::install('", annopkg, "')")
        )
      )
        require(annopkg, character.only = TRUE)
        pkg_choices <- keytypes(get(annopkg))
        std_choices <- union(std_choices, pkg_choices)
        selectInput("idtype", "select your id type", choices = std_choices)
      } else {
        return(NULL)
      }
    })

### main reactive functions
  countData <- reactive({
    inFile <- input$file1_for_DERNAseq
    if (is.null(inFile))
      return(NULL)
     countData <- read.csv(inFile$datapath, sep = "\t", stringsAsFactors = F, header = T)
     row_name = countData[,1]
     countData = countData[ , -1]
     countData = as.matrix(countData)
     rownames(countData) = row_name
     return(countData)
  })

  exdesign <- reactive({
    inFile <- input$file2_for_DERNAseq
    if (is.null(inFile))
      return(NULL)
     exdesign = read.csv(inFile$datapath, sep = "\t", header = T)
     return(exdesign)
  })

  observeEvent(input$analyze_for_DERNAseq, {
    output$downloadTable_for_DERNAseq <- renderUI({
      selectizeInput("dataset_for_DERNAseq",
                     "Choose a dataset to save" ,
                     c("full_results","significant_results"
                       ))
    })

    output$downloadButton_for_DERNAseq <- renderUI({
      downloadButton('download_DERNAseq', 'Save table', class = "download_DERNAseq")
    })

    output$significantBox_for_DERNAseq <- renderInfoBox({
      num_total_RNAseq <- as.data.frame(res_for_RNAseq_1()) %>%
        nrow()
      num_signif_RNAseq <- as.data.frame(res_for_RNAseq_4()$res_for_RNAseq_Filter) %>%
        nrow()
      frac_RNAseq <- num_signif_RNAseq / num_total_RNAseq

        if(frac_RNAseq == 0) {
        info_box_RNAseq <- infoBox("Significant genes",
                            paste0(num_signif_RNAseq,
                                   " out of ",
                                   num_total_RNAseq),
                            "No genes differentially expressed",
                            icon = icon("thumbs-down", lib = "glyphicon"),
                            color = "red",
                            width = 4)
       }
      if(!frac_RNAseq == 0) {
        info_box_RNAseq <-     infoBox("Significant genes",
                              paste0(num_signif_RNAseq,
                                     " out of ",
                                     num_total_RNAseq),
                              paste0(signif(frac_RNAseq * 100, digits = 3),
                                     "% of genes differentially expressed"),
                              icon = icon("thumbs-up", lib = "glyphicon"),
                              color = "green",
                              width = 4)
      }
      info_box_RNAseq
    })

     observe({
      output$plot_heatmap_for_RNAseq <- renderUI({
        plotOutput("heatmap_for_DERNAseq", height = (100 * as.numeric(input$size_for_RNAseq)))
      })
    })


    output$heatmap_cntrst_for_RNAseq <- renderUI({
        df <- res_for_RNAseq_4()$res_for_RNAseq
        cols <- grep("_significant$",colnames(df))
        selectizeInput("heatmap_cntrst_for_RNAseq",
                       "Contrast",
                       choices = gsub("_significant", "", colnames(df)[cols]), multiple = TRUE)
    })


    output$mysplit_for_RNAseq <- renderUI({
        selectizeInput("mysplit_for_RNAseq",
                       "my split",
                       choices = c(1 : input$k_for_RNAseq), multiple = TRUE)

    })

    output$Custom_columns_order_for_RNAseq <- renderUI({
        selectizeInput("Custom_columns_order_for_RNAseq",
                       "Custom columns order",
                       choices = colnames(ntd()), multiple = TRUE)

    })

    output$selected_proteins_for_RNAseq <- renderUI({
      withProgress(message = 'Please Wait',value = 0.66, {
        if(input$labeled_for_RNAseq_customvolcano){
          res_for_RNAseq = res_for_RNAseq_1() %>% as.data.frame(.)
        selectizeInput("selected_proteins_for_RNAseq",
                       "selected genes",
                       choices = if("symbol" %in% colnames(res_for_RNAseq)) {c(NULL,res_for_RNAseq$symbol)} else {c(NULL, rownames(res_for_RNAseq))}, selected = NULL, multiple = TRUE)
        } else {
          return(NULL)
        }
        })
    })

  output$contrast_for_RNAseq_customvolcano <- renderUI({
        df <- res_for_RNAseq_4()$res_for_RNAseq
        cols <- grep("_significant$",colnames(df))
        selectizeInput("contrast_for_RNAseq_customvolcano",
                       "Contrast",
                       choices = gsub("_significant", "", colnames(df)[cols]), selected =  gsub("_significant", "", colnames(df)[cols])[1], multiple = FALSE)
    })

    output$selected_proteins_for_RNAseq_MAplot <- renderUI({
      if(input$selected_genes_for_RNAseq_MAplot) {
        # # Create a Progress object
      #   progress <- shiny::Progress$new()
      #   progress$set(message = "Please wait", value = 0.66)
      #   # Close the progress when this reactive exits (even if there's an error)
        # on.exit(progress$close())
      res_for_RNAseq = res_for_RNAseq_1() %>% as.data.frame(.)
        selectizeInput("selected_proteins_for_RNAseq_MAplot",
                       "selected genes",
                       choices = if("symbol" %in% colnames(res_for_RNAseq)) {c(NULL,res_for_RNAseq$symbol)} else {c(NULL, rownames(res_for_RNAseq))}, selected = NULL, multiple = TRUE)

      }
    })

  output$contrast_for_RNAseq_MA <- renderUI({
        df <- res_for_RNAseq_4()$res_for_RNAseq
        cols <- grep("_significant$",colnames(df))
        selectizeInput("contrast_for_RNAseq_MA",
                       "Contrast",
                       choices = gsub("_significant", "", colnames(df)[cols]), selected = gsub("_significant", "", colnames(df)[cols])[1], multiple = FALSE)
    })

    output$selected_genes_for_RNAseq <- renderUI({
      withProgress(message = 'Please Wait',value = 0.66, {
      res_for_RNAseq = res_for_RNAseq_1() %>% as.data.frame(.)
        selectizeInput("selected_genes_for_RNAseq",
                       "selected genes",
                       choices = if("symbol" %in% colnames(res_for_RNAseq)) {res_for_RNAseq$symbol} else {rownames(res_for_RNAseq)}, selected = if("symbol" %in% colnames(res_for_RNAseq)) {res_for_RNAseq$symbol[1]} else {rownames(res_for_RNAseq)[1]}, multiple = FALSE)

        })
    })

    output$selected_group_for_RNAseq <- renderUI({
        selectizeInput("selected_group_for_RNAseq",
                       "Color by",
                       choices = input$dds_design, selected = input$choose_expfac, multiple = FALSE)

    })



    dds_1 <- reactive({
          DESeqDataSetFromMatrix(countData = countData(),colData = exdesign(), design = as.formula(paste0("~", paste(input$dds_design, collapse = " + "))))
    })

    dds_2 <- reactive({
      dds_1()[rowSums(counts(dds_1())) >= input$filter_rowsum, ]
      })

    dds <- reactive({
      withProgress(message = 'Running DEseq', value = 0.66, {
        choose_factor_num <- exdesign() %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character() %>% unique() %>% length()
      if (input$nrcores == 1) {
        if(choose_factor_num > 2) {
          if(input$perform_LRT_test) {
            DESeq(dds_2(), test = "LRT", full = as.formula(paste0("~", paste(input$dds_design, collapse = " + "))), reduced = as.formula(paste0("~", paste(input$reduced_model, collapse = " + "))))
            } else {
              DESeq(dds_2())
            }
          } else {
            DESeq(dds_2())
          }
        } else {
          if(choose_factor_num > 2) {
            if(input$perform_LRT_test) {
              DESeq(dds_2(), test = "LRT", full = as.formula(paste0("~", paste(input$dds_design, collapse = " + "))), reduced = as.formula(paste0("~", paste(input$reduced_model, collapse = " + "))), parallel = TRUE, BPPARAM = MulticoreParam(workers = input$nrcores))
              } else {
                DESeq(dds_2(), parallel = TRUE, BPPARAM = MulticoreParam(workers = input$nrcores))
              }
            } else {
              DESeq(dds_2(), parallel = TRUE, BPPARAM = MulticoreParam(workers = input$nrcores))
            }
        }
        })
      })

    #res and sig_res
    res_for_RNAseq_1 <- reactive({
      withProgress(message = 'Get results', value = 0.66, {
      if(input$contrasts_for_RNAseq == "control"){
        res_for_RNAseq <<- get_res(dds = dds(), coldata = exdesign(), type = input$contrasts_for_RNAseq,  control = input$control_for_RNAseq, contrast_upon = input$choose_expfac, filter_ihw = input$Use_IHW, independentFiltering = input$independent_filtering, lfcshark = input$Shrink_lfc)
      }
      if(input$contrasts_for_RNAseq == "all") {
        res_for_RNAseq <<- get_res(dds = dds(), coldata = exdesign(), type = input$contrasts_for_RNAseq, contrast_upon = input$choose_expfac, filter_ihw = input$Use_IHW, independentFiltering = input$independent_filtering, lfcshark = input$Shrink_lfc)
      }
      if(input$contrasts_for_RNAseq == "manual") {
        res_for_RNAseq <<- get_res(dds = dds(), coldata = exdesign(), type = input$contrasts_for_RNAseq, test = input$test_manual_for_RNAseq, contrast_upon = input$choose_expfac, filter_ihw = input$Use_IHW, independentFiltering = input$independent_filtering, lfcshark = input$Shrink_lfc)
      }
      if(input$transid_for_RNAseq) {
        ann <<- mapIds(get(annoSpecies_df[input$speciesSelect, ]$pkg), keys = rownames(res_for_RNAseq), column = "SYMBOL", keytype = input$idtype, multiVals = "first")
        res_for_RNAseq$symbol = ann
        # test2 <<- get(annoSpecies_df[input$speciesSelect, ]$pkg)
      if(input$set_no_map_to_rowname) {
          res_for_RNAseq$symbol[which(is.na(res_for_RNAseq$symbol))] = rownames(res_for_RNAseq)[which(is.na(res_for_RNAseq$symbol))]
        }
      }

      return(res_for_RNAseq)
      # res_for_RNAseq$Status = ifelse(res_for_RNAseq$padj < input$p_for_DERNAseq & abs(res_for_RNAseq$log2FoldChange) >= input$lfc_for_DERNAseq, ifelse(res_for_RNAseq$log2FoldChange >= input$lfc_for_DERNAseq, "Up-regulated", "Down-regulated"), "Non-significant")
      # res_for_RNAseq_Filter <- res_for_RNAseq[which(res_for_RNAseq$padj < input$p_for_DERNAseq & abs(res_for_RNAseq$log2FoldChange) >= input$lfc_for_DERNAseq),]
      # res_for_RNAseq_Filter <- res_for_RNAseq_Filter[order(res_for_RNAseq_Filter$padj),]

      # return(list(res_for_RNAseq = res_for_RNAseq, res_for_RNAseq_Filter = res_for_RNAseq_Filter))
        })
      })

    res_for_RNAseq_4 <- reactive({
      res_for_RNAseq_3 <<- add_rejection(diff = res_for_RNAseq_1(), alpha = input$p_for_DERNAseq, lfc = input$lfc_for_DERNAseq)
      res_for_RNAseq_Filter <<- res_for_RNAseq_3[res_for_RNAseq_3$significant, ]
      my_res_for_RNAseq_4 <<- list(res_for_RNAseq = res_for_RNAseq_3, res_for_RNAseq_Filter = res_for_RNAseq_Filter)

      return(my_res_for_RNAseq_4)
      })

    output$table_for_DERNAseq <- DT::renderDataTable(
      as.data.frame(res_for_RNAseq_4()$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname), filter = 'top', options = list( autoWidth = F,scrollX = TRUE
      )
    )

    # data for plot pca, cor, dist
    rld <- reactive({
      rld <- rlog(dds(),blind=FALSE)
      as.data.frame(assay(rld))
      })

    ntd <- reactive({
      ntd = normTransform(dds())
      if(all(rownames(assay(dds())) == rownames(assay(ntd))) & all(rownames(assay(dds())) == rownames(res_for_RNAseq_1()))) {
        assay_ntd <- assay(ntd)
        return(assay_ntd)
      }
      })
    # res merge with ntd, only for download
    res_for_RNAseq_merge <- reactive({
      res_for_RNAseq_merge = as.data.frame(res_for_RNAseq_4()$res_for_RNAseq)
      res_for_RNAseq_merge$gene = rownames(res_for_RNAseq_merge)
      ntd_res <- merge(res_for_RNAseq_merge, as.data.frame(ntd()), by.x = "gene", by.y = "row.names")
      my_res_for_RNAseq_merge <<- ntd_res
      return(ntd_res)
      })

    ### Download objects and functions ### ------------------------------------
    datasetInput_for_RNAseq <- reactive({
          switch(input$dataset_for_DERNAseq,
             "full_results" = res_for_RNAseq_merge() %>% dplyr::arrange(desc(significant)),
             "significant_results" = as.data.frame(res_for_RNAseq_4()$res_for_RNAseq_Filter) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname))
    })

    output$download_DERNAseq <- downloadHandler(
      filename = function() { paste(input$dataset_for_DERNAseq, ".csv", sep = "") },
      content = function(file) {
        write.csv(datasetInput_for_RNAseq(),
        file, row.names = FALSE
        ) }
    )

    data_forRNAseq_heatmap <- reactive({
      if(input$transid_for_RNAseq) {
        data_1 <<- ntd()
        rownames(data_1) <<- res_for_RNAseq_1()$symbol
        } else {
          data_1 <<- ntd()
        }
        return(data_1)
      })
    ### plot
    heatmap_input_for_RNAseq <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        plot_heatmap_rnaseq(data = data_forRNAseq_heatmap(),
                     row_data = res_for_RNAseq_4()$res_for_RNAseq,
                     type = input$pres_for_RNAseq,
                     manual = input$manual_heatmap_for_RNAseq,
                     manual_name = input$heatmap_cntrst_for_RNAseq,
                     kmeans = TRUE,
                     k = input$k_for_RNAseq,
                     color = input$colorbar_for_RNAseq,
                     col_limit = input$limit_for_RNAseq,
                     row_font_size = input$row_font_size_for_RNAseq,
                     col_font_size = input$col_font_size_for_RNAseq,
                     cluster_columns = input$cluster_columns_for_RNAseq,
                     if_mysplit = input$if_mysplit_for_RNAseq,
                     mysplit = input$mysplit_for_RNAseq,
                     column_order = if(input$cluster_columns_for_RNAseq) {NULL} else {input$Custom_columns_order_for_RNAseq})
      })
    })

    custom_volcano_input_for_RNAseq <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
         plot_volcano_rnaseq(res = res_for_RNAseq_4()$res_for_RNAseq,
                     contrast = input$contrast_for_RNAseq_customvolcano,
                     adjusted = input$P_adj_for_RNAseq,
                     labelWay = input$labelWay_for_RNAseq,
                     showNum = input$showNum_for_RNAseq,
                     chooseTolabel = input$selected_proteins_for_RNAseq,
                     fontSize = input$fontSize_for_RNAseq,
                     dotsize = input$dotsize_for_RNAseq,
                     same_width = input$Same_width_for_RNAseq,
                     fcCutoff = input$lfc_for_DERNAseq,
                     adjpCutoff = input$p_for_DERNAseq,
                     label.rectangle = input$if_label_rectangle_for_RNAseq,
                     stroke = input$stroke_for_RNAseq,
                     down_color = input$down_color_for_RNAseq,
                     stable_color = input$stable_color_for_RNAseq,
                     up_color = input$up_color_for_RNAseq
                     )
        })
    })

    pca_input_for_RNAseq <- reactive({
      plot_my_pca(data = rld(), coldata = exdesign(), indicate = input$Indicate_for_RNAseq, if_square = input$if_square_for_RNAseq)
    })

    Pearson_correlation_input_for_RNAseq <- reactive({
      plot_my_cor(data = rld(), coldata = exdesign(), pal = input$Pearson_pal_for_RNAseq, pal_rev = input$Pearson_pal_rev_for_RNAseq, lower = input$Pearson_lower_for_RNAseq, upper = input$Pearson_upper_for_RNAseq, add_values = input$add_values_for_RNAseq_person, value_size = input$value_size_for_RNAseq_person, digits = input$value_digits_for_RNAseq_person)
    })

    Gowers_distance_input_for_RNAseq <- reactive({
      plot_my_dist(data = rld(), coldata = exdesign(), pal = input$Gower_pal_for_RNAseq, pal_rev = input$Gower_pal_rev_for_RNAseq, add_values = input$add_values_for_RNAseq_gower, value_size = input$value_size_for_RNAseq_gower, digits = input$value_digits_for_RNAseq_gower)
    })

    diagno_dispests_input_for_RNAseq <- reactive({
      plotDispEsts(dds())
      })

    MA_plot_input_for_RNAseq <- reactive({
      # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
      plot_ma(res = res_for_RNAseq_4()$res_for_RNAseq, contrast = input$contrast_for_RNAseq_MA, FDR = input$p_for_DERNAseq, log2fc = input$lfc_for_DERNAseq, hlines = input$lfc_for_DERNAseq, intgenes = if(input$labeled_for_RNAseq_MAplot & input$selected_genes_for_RNAseq_MAplot) {input$selected_proteins_for_RNAseq_MAplot} else {NULL}, add_rug = input$add_rug_for_RNAseq)
      })

    selected_genes_input_for_RNAseq <- reactive({
      # Create a Progress object
        # progress <- shiny::Progress$new()
        # progress$set(message = "Plotting", value = 0.66)
        # # Close the progress when this reactive exits (even if there's an error)
        # on.exit(progress$close())
        if("symbol" %in% colnames(as.data.frame(res_for_RNAseq_1()))) {
          Annotation_obj <<- as.data.frame(res_for_RNAseq_1())
          Gene <<- rownames(as.data.frame(res_for_RNAseq_1()))[match(input$selected_genes_for_RNAseq, as.data.frame(res_for_RNAseq_1())$symbol)]
          } else {
            Annotation_obj <<- NULL
            Gene <<- input$selected_genes_for_RNAseq
          }
       print(ggplot_Counts(dds = dds(), gene = Gene, annotation_obj = Annotation_obj, intgroup = input$selected_group_for_RNAseq, ylimZero = input$ylimZero_for_RNAseq, add_labels = input$add_labels_for_RNAseq, labels_repel = input$labels_repel_for_RNAseq))
      })

    output$heatmap_for_DERNAseq <- renderPlot({
      print(heatmap_input_for_RNAseq())
    })

    output$custom_volcano_for_RNAseq <- renderPlot({
      print(custom_volcano_input_for_RNAseq())
    })

    output$pca_for_RNAseq <- renderPlot({
      pca_input_for_RNAseq()
    })

    output$Pearson_correlation_for_RNAseq <- renderPlot({
      Pearson_correlation_input_for_RNAseq()
    })

    output$Gowers_distance_for_RNAseq <- renderPlot({
      Gowers_distance_input_for_RNAseq()
    })

    output$diagno_dispests_for_DERNAseq <- renderPlot({
      diagno_dispests_input_for_RNAseq()
      })

    output$plot_MA_for_RNAseq <- renderPlot({
      MA_plot_input_for_RNAseq()
      })

    output$plot_selected_genes_for_RNAseq <- renderPlot({
      print(selected_genes_input_for_RNAseq())
      })



    output$downloadHeatmap_for_DERNAseq <- downloadHandler(
      filename = 'Heatmap.pdf',
      content = function(file) {
        pdf(file, width = input$Width_for_RNAseq, height = input$Height_for_RNAseq)
        print(heatmap_input_for_RNAseq())
        dev.off()
      }
    )

    output$download_custom_volcano_for_RNAseq <- downloadHandler(
      filename = function() {
        paste0("custom_Volcano_", input$choose_expfac1, "_vs_", input$choose_expfac2, ".pdf")
      },
      content = function(file) {
        pdf(file, width = input$custom_volcano_Width_for_RNAseq, height = input$custom_volcano_Height_for_RNAseq)
        print(custom_volcano_input_for_RNAseq())
        dev.off()
      }
    )

    output$downloadPca_for_RNAseq <- downloadHandler(
      filename = 'Pca.pdf',
      content = function(file) {
        pdf(file, width = input$pca_Width_for_RNAsqe, height = input$pca_Height_for_RNAseq)
        print(pca_input_for_RNAseq())
        dev.off()
      }
    )

    output$download_Pearson_correlation_for_RNAseq <- downloadHandler(
      filename = 'Pearson_correlation.pdf',
      content = function(file) {
        pdf(file, width = input$Pearson_Width_for_RNAseq, height = input$Pearson_Height_for_RNAseq)
        print(Pearson_correlation_input_for_RNAseq())
        dev.off()
      }
    )

    output$download_Gowers_distance_for_RNAseq <- downloadHandler(
      filename = 'Gowers_distance.pdf',
      content = function(file) {
        pdf(file, width = input$Gower_Width_for_RNAseq, height = input$Gower_Height_for_RNAseq)
        print(Gowers_distance_input_for_RNAseq())
        dev.off()
      }
    )

    output$download_diagno_dispests_for_DERNAseq <- downloadHandler(
      filename = 'diagno_dispests.pdf',
      content = function(file) {
        pdf(file, width = input$Diagno_Width_for_RNAseq, height = input$Diagno_Height_for_RNAseq)
        plotDispEsts(dds())
        dev.off()
      }
    )

     output$downloadMAplot_for_RNAseq <- downloadHandler(
      filename = 'MA_plot.pdf',
      content = function(file) {
        pdf(file, width = input$MAplot_Width_for_RNAsqe, height = input$MAplot_Height_for_RNAseq)
        print(MA_plot_input_for_RNAseq())
        dev.off()
      }
    )
     output$download_selected_genes_for_RNAseq <- downloadHandler(
      filename = 'selected_genes.pdf',
      content = function(file) {
        pdf(file, width = input$selected_genes_Width_for_RNAsqe, height = input$selected_genes_Height_for_RNAseq)
        print(selected_genes_input_for_RNAseq())
        dev.off()
      }
    )
     test10 <<- exdesign()
     test11 <<- countData()
    test9 <<- res_for_RNAseq_1()
    test1 <<- res_for_RNAseq_4()
    # test8 <<- selected_genes_input_for_RNAseq()
    test7 <<- dds()
    # test6 <<- heatmap_input_for_RNAseq()
    test12 <<- data_forRNAseq_heatmap()
    test4 <<- res_for_RNAseq_merge()
    test3 <<- ntd()
    test2 <<- rld()
    })

# for PR-Heatmap
    output$file1_for_pro_rna <- renderUI({
      if(!input$import_from_for_PR_Heatmap) {
          fileInput('file1_for_pro_rna',width = "300px",
                    'Proteinmatrix.txt',
                    accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
      }
          })

    output$file2_for_pro_rna <- renderUI({
      if(!input$import_from_for_PR_Heatmap) {
      fileInput('file2_for_pro_rna',width = "300px",
                'RNAseq_normalized_count\nmatrix.txt',
                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
      }
      })

    output$file3_for_pro_rna <- renderUI({
      if(!input$import_from_for_PR_Heatmap) {
      fileInput('file3_for_pro_rna',width = "300px",
                'Genelist.txt',
                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
      }
      })

    output$help_format_PRO_RNA <- renderUI({
      if(!input$import_from_for_PR_Heatmap) {
      actionButton("help_format_PRO_RNA",
                  label = "", icon = icon("question-circle"),
                  style = "color: #f6f6f6; background-color: #2c3b41; border-color: #2c3b41")
      }
      })

    output$import_contrast_for_PR_Heatmap <- renderUI({
      if(input$import_from_for_PR_Heatmap) {
        validate(need(!is.null(data()), "Please go to the DEP-LFQ options, and do the differential analysis first"),
                 need(!is.null(countData()), "Please go to the DEG-RNAseq options, and do the differential analysis first")
          )
            selectizeInput("import_contrast_for_PR_Heatmap", "Contrast", choices = c(get_results(dep()) %>%
        dplyr::filter(significant) %>% colnames() %>% grep("_significant", ., value = T) %>% gsub("_significant", "", .) %>% unique(), "Any significant"), selected = NULL)
      }
      })

    output$Type_for_PRheatmap <- renderUI({
      if(!input$import_from_for_PR_Heatmap) {
        radioButtons(inputId = "Type_for_PRheatmap", label = "Choose the genes to show", choices = c("all", "both", "at least one"), selected = "all", inline = TRUE)
        } else {
          radioButtons(inputId = "Type_for_PRheatmap", label = "Choose the genes to show", choices = c("at least one", "both", "both significant"), selected = "both", inline = TRUE)
        }
      })

  observeEvent(input$analyze_for_pro_rna, {
    withProgress(message = "Please wait...", {
      output$downloadbutton_table_for_pro_rna <- renderUI({
        downloadButton("download_table_for_pro_rna", "Save table", class = "download_table_for_pro_rna")
        })

        ### Reactive functions ### --------------------------------------------------
  if(!input$import_from_for_PR_Heatmap) {
    pro_for_pro_rna <- reactive({
      inFile <- input$file1_for_pro_rna
      if (is.null(inFile))
        return(NULL)
       read.csv(inFile$datapath, header = TRUE,
               sep = "\t", stringsAsFactors = FALSE)
    })

    rna_for_pro_rna <- reactive({
      inFile <- input$file2_for_pro_rna
      if (is.null(inFile))
        return(NULL)
       read.csv(inFile$datapath, header = TRUE,
               sep = "\t", stringsAsFactors = FALSE)
    })

    gene_for_pro_rna <- reactive({
      inFile <- input$file3_for_pro_rna
      if (is.null(inFile)) {
        return(NULL)
        } else {
          file = read.csv(inFile$datapath, header = TRUE,
               sep = "\t", stringsAsFactors = FALSE)
          return(file)
        }
    })
  } else {
    pro_for_pro_rna <- reactive({
      pro = get_df_wide(dep())
      prepare_pro(pro)
      })

    rna_for_pro_rna <- reactive({
      rna = my_res_for_RNAseq_merge %>% dplyr::arrange(desc(significant))
      prepare_rna(rna)

      })
  }

  # test1 <<- pro_for_pro_rna()
  # b <<- rna_for_pro_rna()
  # test4 <<- gene_for_pro_rna()

   observe({
      output$plot_heatmap_for_pro_rna <- renderUI({
        plotOutput("heatmap_for_pro_rna", height = (100 * as.numeric(input$size_for_pro_rna)))
      })
    })



  ht_lis <- reactive({
    if(!input$import_from_for_PR_Heatmap) {
      plot_pro_rna_heatmap(
        rna = rna_for_pro_rna(),
        pro = pro_for_pro_rna(),
        gene = gene_for_pro_rna(),
        Type = input$Type_for_PRheatmap,
        row_font_size = input$row_font_size_for_pro_rna,
        col_font_size = input$col_font_size_for_pro_rna,
        column_title_size = input$column_title_size_for_pro_rna,
        Numcol_pro = input$Numcol_pro_for_pro_rna,
        Numcol_rna = input$Numcol_rna_for_pro_rna,
        cluster_columns = input$cluster_columns_for_pro_rna,
        color = input$colorbar_for_PRheatmap)
      } else {
        if(input$import_contrast_for_PR_Heatmap == "Any significant") {
          manual = FALSE
          contrast = NULL
        } else {
          manual = TRUE
          contrast = input$import_contrast_for_PR_Heatmap
        }
       plot_pro_rna_heatmap_for_whole(
        rna = rna_for_pro_rna(),
        pro = pro_for_pro_rna(),
        manual = manual,
        contrast = contrast,
        Type = input$Type_for_PRheatmap,
        row_font_size = input$row_font_size_for_pro_rna,
        col_font_size = input$col_font_size_for_pro_rna,
        column_title_size = input$column_title_size_for_pro_rna,
        Numcol_pro = input$Numcol_pro_for_pro_rna,
        Numcol_rna = input$Numcol_rna_for_pro_rna,
        cluster_columns = input$cluster_columns_for_pro_rna,
        color = input$colorbar_for_PRheatmap)
      }
    })

     output$table_for_pro_rna <- DT::renderDataTable(
      ht_lis()$dat, filter = 'top', options = list( autoWidth = F,scrollX = TRUE
      )
    )

    output$download_table_for_pro_rna <- downloadHandler(
      filename = function() { paste("data", ".txt", sep = "") },
      content = function(file) {
        write.table(ht_lis()$dat,
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") }
    )

    output$heatmap_for_pro_rna <-  renderPlot({
      print(ht_lis()$ht_merge)
    })

    output$downloadHeatmap_for_pro_rna <- downloadHandler(
      filename = 'Heatmap_merge.pdf',
      content = function(file) {
        pdf(file, width = input$heatmap_Width_for_pro_rna, height = input$heatmap_Height_for_pro_rna)
        print(ht_lis()$ht_merge)
        dev.off()
      }
    )

    test1 <<- ht_lis()
    test2 <<- pro_for_pro_rna()
    test3 <<- rna_for_pro_rna()
    # test4 <<- gene_for_pro_rna()

      })
    })

  ###"question for DEP"
  observeEvent(input$help_format_DEP, {
    showModal(modalDialog(
      title = "Format specifications for DEP-LFQ",
      includeMarkdown(system.file("extdata", "DEP_LFQ.md", package = "DEP2")),
      h4("Example:"),
      tags$img(
        src = base64enc::dataURI(file = system.file("extdata", "DEP_LFQ.png", package = "DEP2"), mime = "image/png"),
        width = 750
      ),
      easyClose = TRUE,
      footer = NULL,
      size = "l"
    ))
  })
  ### help imputation
  observeEvent(input$help_imputation, {
          showModal(modalDialog(
            title = "The detailed information of imputation methods",
            includeMarkdown(system.file("extdata", "impute.md", package = "DEP2")),
            easyClose = TRUE,
            footer = NULL,
            size = "l"
          ))
        })

   ###"question for RNA"
  observeEvent(input$help_format_RNA, {
    showModal(modalDialog(
      title = "Format specifications for DEG-RNAseq",
      includeMarkdown(system.file("extdata", "DEG_RNAseq.md", package = "DEP2")),
      h4("Example:"),
      tags$img(
        src = base64enc::dataURI(file = system.file("extdata", "DEG_RNAseq.png", package = "DEP2"), mime = "image/png"),
        width = 750
      ),
      easyClose = TRUE,
      footer = NULL,
      size = "l"
    ))
  })

  ###"question for pro-RNA"
  observeEvent(input$help_format_PRO_RNA, {
    showModal(modalDialog(
      title = "Format specifications for PR-Heatmap",
      includeMarkdown(system.file("extdata", "PR_Heatmap.md", package = "DEP2")),
      h4("Example:"),
      tags$img(
        src = base64enc::dataURI(file = system.file("extdata", "PR_Heatmap.png", package = "DEP2"), mime = "image/png"),
        width = 750
      ),
      easyClose = TRUE,
      footer = NULL,
      size = "l"
    ))
  })






})


shinyApp(ui = ui, server = server)
