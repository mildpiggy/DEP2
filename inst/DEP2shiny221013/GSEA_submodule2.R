
#* GO submodule UI/Server in GSEA module ----
GSEA_GO_UI <- function(id,tab_name){
  ns = NS(id)
  tabPanel(tab_name,
           fluidRow(
             box(selectizeInput(ns("gsego_ont"),
                                "Ontology",
                                choices = c("BP", "CC", "MF", "ALL"),
                                selected = c("BP"), multiple = FALSE),
                 width = 2),
             box(numericInput(ns("gsego_p"),
                              "P value",
                              min = 0, max = 1, value = 0.05),
                 width = 2),
             box(numericInput(ns("gsego_padj"),
                              "adj. P value",
                              min = 0, max = 1, value = 0.25),
                 width = 2),
             box(numericInput(ns("gsego_NES"),
                              "NES",
                              min = 0, max = 10, value = 1),
                 width = 2),
             box(selectizeInput(ns("gsego_Phenotype"),
                                "Phenotype",
                                choices = c("activated", "suppressed"),
                                selected = c("activated", "suppressed"), multiple = TRUE), width = 2),
             box(radioButtons(ns("gsego_color"),
                              "colorBy",
                              c("pvalue", "p.adjust"),
                              selected = "p.adjust"), width = 2),
             shinyBS::bsTooltip(ns("gsego_ont"), "One of BP, MF, and CC subontologies, or ALL for all three", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("gsego_p"), "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("gsego_padj"), "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("gsego_NES"), "the |NES| (The absolute value of normalized enrichment score) cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("gsego_Phenotype"), "The phenotype that you want to show, one or both of activated and suppressed. For example, if your log2 fold change is A_vs_B, then, term enriched in A is activated, and term enriched in B is suppressed", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("gsego_color"), "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
           ),
           fluidRow(
             column(width = 7,
                    fluidRow(
                      box(checkboxInput(ns("gsego_simplify"), "removed redundancy of enriched GO terms", value = FALSE), width = 5),
                      # box(radioButtons("go_color",
                      #          "colorBy",
                      #          c("pvalue", "p.adjust"),
                      #          selected = "p.adjust"), width = 4),
                      infoBoxOutput(ns("significantBox_for_gsego"), width = 7)
                    ),
                    fluidRow(
                      box(title = "Result Table",
                          fluidRow(
                            column(width = 5,
                                   uiOutput(ns("downloadTable_gsego"))
                            )
                          ),
                          fluidRow(
                            column(width = 4,
                                   uiOutput(ns("downloadButton_gsego"))
                            )
                          ),  
                          tags$br(),                    
                          DT::dataTableOutput(ns("gsego_Table")), width = 12)
                    ),
                    shinyBS::bsTooltip(ns("gsego_simplify"), "Whether simplify output by removing redundancy of enriched GO terms", "top", options = list(container = "body")),
                    shinyBS::bsTooltip(ns("gsego_Table"), "The table of the significant result", "top", options = list(container = "body")),
                    shinyBS::bsTooltip(ns("downloadTable_gsego"), "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
                    shinyBS::bsTooltip(ns("downloadButton_gsego"), "Click on it to download the result table", "right", options = list(container = "body"))
             ),
             column(width = 5,
                    tabBox(title = "Result Plots", width = 12,
                           tabPanel(title = "Bar plot",
                                    fluidRow(box(checkboxInput(ns("gsego_Select_bar"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("gsego_ShowCategory_bar"), "ShowCategory", min = 1, max = 100, value = 20),
                                          width = 4),
                                      box(numericInput(ns("gsego_wide_bar"), "Width", min = 1, max = 50, value = 12),
                                          width = 4),
                                      box(numericInput(ns("gsego_high_bar"), "Height", min = 1, max = 50, value = 10),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      box(checkboxInput(ns("gsego_bar_if_Split_for_ont_ALL"),
                                                        "If splited by ontology",
                                                        value = FALSE), width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gsego_barplot")),
                                      downloadButton(ns('download_gsego_barplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gsego_Select_bar"), "Display the selected genesets in res table or top genesets(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_ShowCategory_bar"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_wide_bar"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_high_bar"), "Height of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_bar_if_Split_for_ont_ALL"), "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of each phenotype that you selected of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each phenotype that you selected based on each BP, CC and MF", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Dot plot",
                                    fluidRow(box(checkboxInput(ns("gsego_Select_dot"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("gsego_ShowCategory_dot"), "ShowCategory", min = 1, max = 100, value = 20),
                                          width = 4),
                                      box(numericInput(ns("gsego_wide_dot"), "Width", min = 1, max = 50, value = 12),
                                          width = 4),
                                      box(numericInput(ns("gsego_high_dot"), "Height", min = 1, max = 50, value = 10),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      box(checkboxInput(ns("gsego_dot_if_Split_for_ont_ALL"),
                                                        "If splited by ontology",
                                                        value = FALSE), width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gsego_dotplot"), height = 500),
                                      downloadButton(ns('download_gsego_dotplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gsego_Select_dot"), "Display the selected genesets in res table or top genesets(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_ShowCategory_dot"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_wide_dot"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_high_dot"), "Height of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_dot_if_Split_for_ont_ALL"), "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of each phenotype that you selected of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each phenotype that you selected based on each BP, CC and MF", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Heatplot",
                                    fluidRow(box(checkboxInput(ns("gsego_Select_heat"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("gsego_ShowCategory_heat"), "ShowCategory", min = 1, max = 100, value = 30),
                                          width = 4),
                                      box(numericInput(ns("gsego_wide_heat"), "Width", min = 1, max = 200, value = 80),
                                          width = 4),
                                      box(numericInput(ns("gsego_high_heat"), "Height", min = 1, max = 50, value = 7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gsego_heatplot"), height = 500),
                                      downloadButton(ns('download_gsego_heatplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gsego_Select_heat"), "Display the selected genesets in res table or top genesets(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_ShowCategory_heat"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_wide_heat"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_high_heat"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Cnetplot",
                                    fluidRow(box(checkboxInput(ns("gsego_Select_cnet"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("gsego_ShowCategory_cnet"), "ShowCategory", min = 1, max = 100, value = 5),
                                          width = 4),
                                      box(numericInput(ns("gsego_wide_cnet"), "Width", min = 1, max = 200, value = 30),
                                          width = 4),
                                      box(numericInput(ns("gsego_high_cnet"), "Height", min = 1, max = 200, value = 25),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      box(checkboxInput(ns("gsego_circular_cnet"),
                                                        "Circular",
                                                        value = TRUE),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gsego_cnetplot"), height = 500),
                                      downloadButton(ns('download_gsego_cnetplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gsego_Select_cnet"), "Display the selected genesets in res table or top genesets(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_ShowCategory_cnet"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_wide_cnet"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_high_cnet"), "Height of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_circular_cnet"), "whether using circular layout", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Emaplot",
                                    fluidRow(box(checkboxInput(ns("gsego_Select_ema"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("gsego_ShowCategory_ema"), "ShowCategory", min = 1, max = 100, value = 30),
                                          width = 4),
                                      box(numericInput(ns("gsego_wide_ema"), "Width", min = 1, max = 50, value = 11),
                                          width = 4),
                                      box(numericInput(ns("gsego_high_ema"), "Height", min = 1, max = 50, value = 10),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gsego_emaplot"), height = 500),
                                      downloadButton(ns('download_gsego_emaplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gsego_Select_ema"), "Display the selected genesets in res table or top genesets(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_ShowCategory_ema"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_wide_ema"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_high_ema"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Gseaplot",
                                    fluidRow(box(uiOutput(ns("gsego_term")), width = 12)),
                                    fluidRow(
                                      box(numericInput(ns("gsego_wide_Gsea"), "Width", min = 1, max = 50, value = 8),
                                          width = 4),
                                      box(numericInput(ns("gsego_high_Gsea"), "Height", min = 1, max = 50, value = 7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gsego_Gseaplot"), height = 500),
                                      downloadButton(ns('download_gsego_Gseaplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gsego_term"), "The term that you want to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_wide_Gsea"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gsego_high_Gsea"), "Height of the figure to export", "top", options = list(container = "body"))
                           )
                    )
                    
             )
           )
  )
}

GSEA_GO_server_module <- function(id, gene_list, organism_for_GSEA, annoSpecies_df) {
  moduleServer(
    id = id,
    function(input, output, session) {
      print(id)
      print(organism_for_GSEA)
      # print(annoSpecies_df)
      print(gene_list)
      
      pkg_for_GSEA <- annoSpecies_df$pkg[annoSpecies_df$species == organism_for_GSEA]
      print(pkg_for_GSEA)
      output$downloadTable_gsego <- renderUI({
        selectizeInput(session$ns("dataset_for_gsego"),
                       "Choose a dataset to save" ,
                       c("full_results","significant_results"
                       ))
      })
      
      output$downloadButton_gsego <- renderUI({
        downloadButton(session$ns('downloadgsego'), 'Save table', class = "downloadgsego")
      })
     
      reat_gsego <- reactive({
        # cat(paste0("theOrg1:", organism_for_GSEA))
        try(gsegoAnalysis(gene_list = gene_list, organism = organism_for_GSEA, species_df = annoSpecies_df), silent = TRUE)
      })
      
      res_gsego <- reactive({
        #   validate(
        #   need(length(input$gsego_Phenotype) != 0, "Please choose at least one Phenotype"),
        #   need(class(reat_gsego()) != "try-error", "Your order ranked geneList can not do gsea analysis, Please go to GO panel if desired")
        # )
        reat_gsego_save <<- reat_gsego()
        try(give_gseGO_res_and_table(reat = reat_gsego(), ont = input$gsego_ont, pCutoff = input$gsego_p, p.adj.cutoff = input$gsego_padj, NES.cutoff = input$gsego_NES, simplify = input$gsego_simplify, Phenotype = input$gsego_Phenotype), silent = TRUE)
      })
      
      test22222 <<- res_gsego()
      
      output$gsego_Table <- DT::renderDataTable({
        shiny::validate(need(require(pkg_for_GSEA, character.only = TRUE), message = paste0("The package ", pkg_for_GSEA, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_GSEA, "')")))
        shiny::validate(
          need(!is.null(gene_list) && length(gene_list) > 0, message = "No genes meet your requirements, and can not do the gseGO analysis")
        )
        shiny::validate(
          need(length(input$gsego_Phenotype) != 0, "Please choose at least one Phenotype"),
          need(reat_gsego(), "Sorry, I can not do the gseGO analysis")
        )
        
        validate(
          need(nrow(res_gsego()$all_table) != 0, "Sorry, I can not do the gseGO analysis")
        )
        
        res_gsego()$sig_table
        # DT::datatable(res_gsego()$sig_table, filter = 'top', options = list(autoWidth = F,scrollX = TRUE,pageLength = 25), selection = list(selected = c(1:10)) )
      }, 
      options = list(pageLength = 25, scrollX = T),
      filter = 'top',
      selection = list(selected = c(1))
      ) 
      
      
      output$significantBox_for_gsego <- renderInfoBox({
        shiny::validate(need(require(pkg_for_GSEA, character.only = TRUE), message = paste0("The package ", pkg_for_GSEA, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_GSEA, "')")))
        shiny::validate(
          need(!is.null(gene_list) && length(gene_list) > 0, message = "No genes meet your requirements, and can not do the gseGO analysis")
        )
        shiny::validate(
          need(length(input$gsego_Phenotype) != 0, "Please choose at least one Phenotype"),
          need(reat_gsego(), "Sorry, I can not do the gseGO analysis")
        )
        
        shiny::validate(
          need(nrow(res_gsego()$all_table) != 0, "No significantly changed geneset")
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
          info_box_gsego <- infoBox("Significant terms",
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
          selectizeInput(session$ns("gsego_term"),
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
        if(input$gsego_Select_bar){
          if(is.null(input$gsego_Table_rows_selected)) return(NULL)
          res <- res_gsego()$sig_res %>% filter(.,1:nrow(.) %in% input$gsego_Table_rows_selected)
        }else{
          res <- res_gsego()$sig_res
        }
        gse_barplot(res = res, ShowCategory = input$gsego_ShowCategory_bar, color = input$gsego_color, ont = input$gsego_ont, Split = input$gsego_bar_if_Split_for_ont_ALL)
      })
      
      gsego_dotplot_input <- reactive({
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        if(input$gsego_Select_dot){
          if(is.null(input$gsego_Table_rows_selected)) return(NULL)
          res <- res_gsego()$sig_res %>% filter(.,1:nrow(.) %in% input$gsego_Table_rows_selected)
        }else{
          res <- res_gsego()$sig_res
        }
        on.exit(progress$close())
        gse_dotplot(res = res, ShowCategory = input$gsego_ShowCategory_dot, color = input$gsego_color, ont = input$gsego_ont, Split = input$gsego_dot_if_Split_for_ont_ALL)
      })
      
      gsego_heatplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          print(input$gsego_Table_rows_selected)
          showCategory = input$gsego_ShowCategory_heat
          if(input$gsego_Select_heat){
            if(is.null(input$gsego_Table_rows_selected)) return(NULL)
            res <- res_gsego()$sig_res %>% filter(.,1:nrow(.) %in% input$gsego_Table_rows_selected)
            # showCategory = nrow(res)
          }else{
            res <- res_gsego()$sig_res
          }
          try(Heatplot(res, showCategory = showCategory, foldChange = res_gsego()$de), silent = T)
        })
      })
      
      gsego_cnetplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          if(input$gsego_Select_cnet){
            if(is.null(input$gsego_Table_rows_selected)) return(NULL)
            res <- res_gsego()$sig_res %>% filter(.,1:nrow(.) %in% input$gsego_Table_rows_selected)
            # showCategory = nrow(res)
          }else{
            res <- res_gsego()$sig_res
          }
          try(Cnetplot(x = res, showCategory = input$gsego_ShowCategory_cnet, foldChange = res_gsego()$de, circular = input$gsego_circular_cnet, colorEdge = TRUE), silent = T)
        })
      })
      
      gsego_emaplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          if(input$gsego_Select_ema){
            if(is.null(input$gsego_Table_rows_selected)) return(NULL)
            res <- res_gsego()$sig_res %>% filter(.,1:nrow(.) %in% input$gsego_Table_rows_selected)
            # showCategory = nrow(res)
          }else{
            res <- res_gsego()$sig_res
          }
          try(Emapplot(res, showCategory = input$gsego_ShowCategory_ema, color = input$gsego_color, layout = "kk"), silent = T)
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
      }, height = input$gsego_high_bar * 50)
      
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
      
    }
  )}



#* other submodule UI/Server in GSEA module ----
GSEA_other_UI <- function(id, tab_name = id){
  ns = NS(id)
  tabPanel(tab_name,
           fluidRow(
             box(numericInput(ns("gseOther_p"),
                              "P value",
                              min = 0, max = 1, value = 0.05),
                 width = 2),
             box(numericInput(ns("gseOther_padj"),
                              "adj. P value",
                              min = 0, max = 1, value = 0.25),
                 width = 2),
             box(numericInput(ns("gseOther_NES"),
                              "NES",
                              min = 0, max = 10, value = 1),
                 width = 2),
             box(selectizeInput(ns("gseOther_Phenotype"),
                                "Phenotype",
                                choices = c("activated", "suppressed"),
                                selected = c("activated", "suppressed"), multiple = TRUE), width = 4),
             box(radioButtons(ns("gseOther_color"),
                              "colorBy",
                              c("pvalue", "p.adjust"),
                              selected = "p.adjust"), width = 2),
             shinyBS::bsTooltip(ns("gseOther_p"), "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("gseOther_padj"), "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("gseOther_NES"), "the |NES| (The absolute value of normalized enrichment score) cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("gseOther_Phenotype"), "The phenotype that you want to show, one or both of activated and suppressed. For example, if your log2 fold change is A_vs_B, then, term enriched in A is activated, and term enriched in B is suppressed", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("gseOther_color"), "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
           ),
           fluidRow(
             column(width = 7,
                    fluidRow(
                      infoBoxOutput(ns("significantBox_for_gseOther"), width = 7)
                    ),
                    fluidRow(
                      box(title = "Result Table",
                          fluidRow(
                            column(width = 5,
                                   uiOutput(ns("downloadTable_gseOther"))
                            )
                          ),
                          fluidRow(
                            column(width = 4,
                                   uiOutput(ns("downloadButton_gseOther"))
                            )
                          ),  
                          tags$br(),                    
                          DT::dataTableOutput(ns("gseOther_Table")), width = 12)
                    ),
                    shinyBS::bsTooltip(ns("gseOther_Table"), "The table of the significant result", "top", options = list(container = "body")),
                    shinyBS::bsTooltip(ns("downloadTable_gseOther"), "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
                    shinyBS::bsTooltip(ns("downloadButton_gseOther"), "Click on it to download the result table", "right", options = list(container = "body"))
             ),
             column(width = 5,
                    tabBox(title = "Result Plots", width = 12,
                           tabPanel(title = "Bar plot",
                                    fluidRow(
                                      box(numericInput(ns("gseOther_ShowCategory_bar"), "ShowCategory", min = 1, max = 100, value = 20),
                                          width = 4),
                                      box(numericInput(ns("gseOther_wide_bar"), "Width", min = 1, max = 50, value = 12),
                                          width = 4),
                                      box(numericInput(ns("gseOther_high_bar"), "Height", min = 1, max = 50, value = 10),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gseOther_barplot")),
                                      downloadButton(ns('download_gseOther_barplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gseOther_ShowCategory_bar"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_wide_bar"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_high_bar"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Dot plot",
                                    fluidRow(
                                      box(numericInput(ns("gseOther_ShowCategory_dot"), "ShowCategory", min = 1, max = 100, value = 20),
                                          width = 4),
                                      box(numericInput(ns("gseOther_wide_dot"), "Width", min = 1, max = 50, value = 12),
                                          width = 4),
                                      box(numericInput(ns("gseOther_high_dot"), "Height", min = 1, max = 50, value = 10),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gseOther_dotplot"), height = 500),
                                      downloadButton(ns('download_gseOther_dotplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gseOther_ShowCategory_dot"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_wide_dot"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_high_dot"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Heatplot",
                                    fluidRow(
                                      box(numericInput(ns("gseOther_ShowCategory_heat"), "ShowCategory", min = 1, max = 100, value = 30),
                                          width = 4),
                                      box(numericInput(ns("gseOther_wide_heat"), "Width", min = 1, max = 200, value = 80),
                                          width = 4),
                                      box(numericInput(ns("gseOther_high_heat"), "Height", min = 1, max = 50, value = 7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gseOther_heatplot"), height = 500),
                                      downloadButton(ns('download_gseOther_heatplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gseOther_ShowCategory_heat"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_wide_heat"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_high_heat"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Cnetplot",
                                    fluidRow(
                                      box(numericInput(ns("gseOther_ShowCategory_cnet"), "ShowCategory", min = 1, max = 100, value = 5),
                                          width = 4),
                                      box(numericInput(ns("gseOther_wide_cnet"), "Width", min = 1, max = 200, value = 30),
                                          width = 4),
                                      box(numericInput(ns("gseOther_high_cnet"), "Height", min = 1, max = 200, value = 25),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      box(checkboxInput(ns("gseOther_circular_cnet"),
                                                        "Circular",
                                                        value = TRUE),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gseOther_cnetplot"), height = 500),
                                      downloadButton(ns('download_gseOther_cnetplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gseOther_ShowCategory_cnet"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_wide_cnet"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_high_cnet"), "Height of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_circular_cnet"), "whether using circular layout", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Emaplot",
                                    fluidRow(
                                      box(numericInput(ns("gseOther_ShowCategory_ema"), "ShowCategory", min = 1, max = 100, value = 30),
                                          width = 4),
                                      box(numericInput(ns("gseOther_wide_ema"), "Width", min = 1, max = 50, value = 11),
                                          width = 4),
                                      box(numericInput(ns("gseOther_high_ema"), "Height", min = 1, max = 50, value = 10),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gseOther_emaplot"), height = 500),
                                      downloadButton(ns('download_gseOther_emaplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gseOther_ShowCategory_ema"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_wide_ema"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_high_ema"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Gseaplot",
                                    fluidRow(box(uiOutput(ns("gseOther_term")), width = 12)),
                                    fluidRow(
                                      box(numericInput(ns("gseOther_wide_Gsea"), "Width", min = 1, max = 50, value = 8),
                                          width = 4),
                                      box(numericInput(ns("gseOther_high_Gsea"), "Height", min = 1, max = 50, value = 7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("gseOther_Gseaplot"), height = 500),
                                      downloadButton(ns('download_gseOther_Gseaplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("gseOther_term"), "The term that you want to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_wide_Gsea"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("gseOther_high_Gsea"), "Height of the figure to export", "top", options = list(container = "body"))
                           )
                    )
                    
             )
           )
  )
}

GSEA_other_server_module <- function(id, gene_list, organism_for_GSEA, annoSpecies_df, 
                                     type, Msigdb_category = NULL , Msigdb_subcategory =NULL) {
  moduleServer(
    id = id,
    function(input, output, session) {
      print(id)
      print(organism_for_GSEA)
      # print(annoSpecies_df)
      print(gene_list)
      
      pkg_for_GSEA <- annoSpecies_df$pkg[annoSpecies_df$species == organism_for_GSEA]
      print(pkg_for_GSEA)
      
      if(type == "kegg"){
        # reat_gseOther <<- reactive({
        # Create a Progress object
        # progress <- shiny::Progress$new()
        # progress$set(message = "I am doing gseOther analysis, please wait ...", value = 0.66)
        # # Close the progress when this reactive exits (even if there's an error)
        # on.exit(progress$close())
        reat_gseOther <- reactive(try(gsekeggAnalysis(gene_list = gene_list, organism = organism_for_GSEA, species_df = annoSpecies_df), silent = TRUE))
        
        res_gseOther <- reactive({
          try(give_gseOther_res_and_table(reat = reat_gseOther(), pCutoff = input$gseOther_p, p.adj.cutoff = input$gseOther_padj, NES.cutoff = input$gseOther_NES, Phenotype = input$gseOther_Phenotype), silent = TRUE)
        })
      }else if(type == "reactome"){
        reat_gseOther <- reactive(try(gsereactAnalysis(gene_list = gene_list, organism = organism_for_GSEA, species_df = annoSpecies_df), silent = TRUE))
        
        res_gseOther <- reactive({
          try(give_gseOther_res_and_table(reat = reat_gseOther(), pCutoff = input$gseOther_p, p.adj.cutoff = input$gseOther_padj, NES.cutoff = input$gseOther_NES, Phenotype = input$gseOther_Phenotype), silent = TRUE)
        })
      }else if(type == "Msigdb"){
        reat_gseOther <- reactive(try(
          gsemsigdb_Analysis(gene_list = gene_list, organism = organism_for_GSEA, species_df = annoSpecies_df, 
                             category = Msigdb_category , subcategory =Msigdb_subcategory
                             # msigdb = "c2.cgp (chemical and genetic perturbations)"
                             ), silent = TRUE
        ))
        reat_gseOther_save <<- reat_gseOther
        res_gseOther <- reactive({
          try(give_gseOther_res_and_table(reat = reat_gseOther(), pCutoff = input$gseOther_p, p.adj.cutoff = input$gseOther_padj, NES.cutoff = input$gseOther_NES, Phenotype = input$gseOther_Phenotype), silent = TRUE)
        })
        res_gseOther_save <<- res_gseOther
      }
      
      ### Interactive UI functions ### ------------------------------------------
      output$downloadTable_gseOther <- renderUI({
        selectizeInput(session$ns("dataset_for_gseOther"),
                       "Choose a dataset to save" ,
                       c("full_results","significant_results"
                       ))
      })
      
      output$downloadButton_gseOther <- renderUI({
        downloadButton(session$ns('downloadgseOther'), 'Save table', class = "downloadgseOther")
      })
      
      
      
      test2 <<- res_gseOther()
      
      output$gseOther_Table <- DT::renderDataTable({
        shiny::validate(
          need(require(pkg_for_GSEA, character.only = TRUE), 
               message = paste0("The package ", pkg_for_GSEA, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_GSEA, "')")
               ))
        shiny::validate(
          need(!is.null(gene_list) && length(gene_list) > 0, message = "No genes meet your requirements, and can not do the gseGO analysis")
        )
        
        shiny::validate(
          need(length(input$gseOther_Phenotype) != 0, "Please choose at least one Phenotype"),
          need(reat_gseOther(), "Sorry, I can not do the gseOther analysis")
        )
        shiny::validate(
          need(nrow(res_gseOther()$all_table) != 0, "Sorry, I can not do the gseOther analysis")
        )
        # DT::datatable(res_gseOther()$sig_table, filter = 'top', options = list(autoWidth = F,scrollX = TRUE))
        res_gseOther()$sig_table
        # DT::datatable(res_gsego()$sig_table, filter = 'top', options = list(autoWidth = F,scrollX = TRUE,pageLength = 25), selection = list(selected = c(1:10)) )
      }, 
      options = list(pageLength = 25, scrollX = T),
      filter = 'top',
      selection = list(selected = c(1))
      ) 
      
      output$significantBox_for_gseOther <- renderInfoBox({
        shiny::validate(
          need(require(pkg_for_GSEA, character.only = TRUE), 
               message = paste0("The package ", pkg_for_GSEA, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_GSEA, "')")
          ))
        shiny::validate(
          need(!is.null(gene_list) && length(gene_list) > 0, message = "No genes meet your requirements, and can not do the gseGO analysis")
        )
        
        
        validate(
          need(length(input$gseOther_Phenotype) != 0, "Please choose at least one Phenotype"),
          need(reat_gseOther(), "Sorry, I can not do the gseOther analysis")
        )
        validate(
          need(nrow(res_gseOther()$all_table) != 0, "Sorry, I can not do the gseOther analysis")
        )
        num_total_gseOther <- res_gseOther()$all_table %>%
          nrow()
        num_signif_gseOther <- res_gseOther()$sig_table %>%
          nrow()
        frac_gseOther <- num_signif_gseOther / num_total_gseOther
        
        if(frac_gseOther == 0) {
          info_box_gseOther <- infoBox("Significant terms",
                                      paste0(num_signif_gseOther,
                                             " out of ",
                                             num_total_gseOther),
                                      "No terms enriched",
                                      icon = icon("thumbs-down", lib = "glyphicon"),
                                      color = "red",
                                      width = 4)
        }
        if(!frac_gseOther == 0) {
          info_box_gseOther <- infoBox("Significant terms",
                                       paste0(num_signif_gseOther,
                                              " out of ",
                                              num_total_gseOther),
                                       paste0(signif(frac_gseOther * 100, digits = 3),
                                              "% of terms enriched"),
                                       icon = icon("thumbs-up", lib = "glyphicon"),
                                       color = "green",
                                       width = 4)
        }
        info_box_gseOther
      })
      
      output$gseOther_term <- renderUI({
        if (!nrow(res_gseOther()$sig_table) == 0) {
          selectizeInput(session$ns("gseOther_term"),
                         "Term",
                         choices = res_gseOther()$sig_table$Description, selected = res_gseOther()$sig_table$Description[1], multiple = FALSE)
        }
      })
      
      gseOther_barplot_input <- reactive({
        # # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        enrichplot:::barplot.enrichResult(res_gseOther()$sig_res, showCategory = input$gseOther_ShowCategory_bar, color = input$gseOther_color, x = "NES", split="phenotype") + labs(y = "NES")
      })
      
      gseOther_dotplot_input <- reactive({
        progress <- shiny::Progress$new()
        progress$set(message = "Plotting", value = 0.66)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        enrichplot:::dotplot(res_gseOther()$sig_res, showCategory = input$gseOther_ShowCategory_dot, color = input$gseOther_color, x = "NES", split="phenotype")
      })
      
      gseOther_heatplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          try(Heatplot(res_gseOther()$sig_res, showCategory = input$gseOther_ShowCategory_heat, foldChange = res_gseOther()$de), silent = T)
        })
      })
      
      gseOther_cnetplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          try(Cnetplot(x = res_gseOther()$sig_res, showCategory = input$gseOther_ShowCategory_cnet, foldChange = res_gseOther()$de, circular = input$gseOther_circular_cnet, colorEdge = TRUE), silent = T)
        })
      })
      
      gseOther_emaplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          try(Emapplot(res_gseOther()$sig_res, showCategory = input$gseOther_ShowCategory_ema, color = input$gseOther_color, layout = "kk"), silent = T)
        })
      })
      
      gseOther_Gseaplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          my_gseaplot2(res_gseOther()$sig_res, geneSetID = match(input$gseOther_term, res_gseOther()$sig_res$Description), title = input$gseOther_term,# geneSetID = input$gseOther_Table_rows_selected, title = res_gseOther()$sig_table[input$gseOther_Table_rows_selected,3]
                       color = "green",
                       base_size = 11,
                       rel_heights = c(1.5, 0.5, 1),
                       subplots = 1:3,
                       pvalue_table = TRUE,
                       ES_geom = "line")
        })
      })
      
      output$gseOther_barplot <- renderPlot({
        gseOther_barplot_input()
      }, height = input$gseOther_high_bar*50 )
      
      output$gseOther_dotplot <- renderPlot({
        gseOther_dotplot_input()
      })
      
      output$gseOther_heatplot <- renderPlot({
        gseOther_heatplot_input()
      })
      
      output$gseOther_cnetplot <- renderPlot({
        gseOther_cnetplot_input()
      })
      
      output$gseOther_emaplot <- renderPlot({
        gseOther_emaplot_input()
      })
      
      output$gseOther_Gseaplot <- renderPlot({
        gseOther_Gseaplot_input()
      })
      
      ### Download objects and functions ### ------------------------------------
      datasetInput_for_gseOther <- reactive({
        table_for_gseOther = res_gseOther()
        switch(input$dataset_for_gseOther,
               "full_results" = res_gseOther()$all_table,
               "significant_results" = res_gseOther()$sig_table)
      })
      
      output$downloadgseOther <- downloadHandler(
        filename = function() { paste(input$dataset_for_gseOther, ".txt", sep = "") },
        content = function(file) {
          write.table(datasetInput_for_gseOther(),
                      file,
                      col.names = TRUE,
                      row.names = FALSE,
                      sep ="\t") }
      )
      
      output$download_gseOther_barplot <- downloadHandler(
        filename = 'barplot.pdf',
        content = function(file) {
          pdf(file, width = input$gseOther_wide_bar, height = input$gseOther_high_bar)
          print(gseOther_barplot_input())
          dev.off()
        }
      )
      
      output$download_gseOther_dotplot <- downloadHandler(
        filename = 'dotplot.pdf',
        content = function(file) {
          pdf(file, width = input$gseOther_wide_dot, height = input$gseOther_high_dot)
          print(gseOther_dotplot_input())
          dev.off()
        }
      )
      
      output$download_gseOther_heatplot <- downloadHandler(
        filename = 'heatplot.pdf',
        content = function(file) {
          pdf(file, width = input$gseOther_wide_heat, height = input$gseOther_high_heat)
          print(gseOther_heatplot_input())
          dev.off()
        }
      )
      
      output$download_gseOther_cnetplot <- downloadHandler(
        filename = 'cnetplot.pdf',
        content = function(file) {
          pdf(file, width = input$gseOther_wide_cnet, height = input$gseOther_high_cnet)
          print(gseOther_cnetplot_input())
          dev.off()
        }
      )
      
      output$download_gseOther_emaplot <- downloadHandler(
        filename = 'emaplot.pdf',
        content = function(file) {
          pdf(file, width = input$gseOther_wide_ema, height = input$gseOther_high_ema)
          print(gseOther_emaplot_input())
          dev.off()
        }
      )
      
      output$download_gseOther_Gseaplot <- downloadHandler(
        filename = 'Gseaplot.pdf',
        content = function(file) {
          pdf(file, width = input$gseOther_wide_Gsea, height = input$gseOther_high_Gsea)
          print(gseOther_Gseaplot_input())
          dev.off()
        }
      )
      
      
    }
  )
}



#* Msigdb submodule UI/Server in GSEA module ----
GSEA_Msigdb_UI <- function(id, tab_name = "Msigdb"){
  ns = NS(id)
  # id_save2 <<- id
  tabPanel(
    tab_name,
    uiOutput(ns("Msigdb_tabs"))
  )
}

GSEA_Msigdb_server_module <- function(id, gene_df, gene_list, organism_for_GSEA, annoSpecies_df, Msigdb_selection = NULL, Msigdb_selection2 = NULL) {
  # id = paste0(id,"_",suffix)
  # id_save <<- id
  
  moduleServer(
    id = id,
    function(input, output, session) {
      ns = session$ns
      print("000")
      Msigdb_selection_save <<- Msigdb_selection
      Msigdb_selection2_save <<- Msigdb_selection2
      print("111")
      #** Msigdb_tabs render in msigdb ----
      output$Msigdb_tabs <- renderUI({
        # tabs_save <<- tabs
        print("aa")
        tabsetPanel(
          id = session$ns("Msigdb_Tabs")
        )
      })
      observeEvent(input$Msigdb_Tabs,{
        
        if(length(Msigdb_selection) > 0){
          for (i in 1:length(Msigdb_selection)) { ## append subtabs UI
            print("bb")
            # tabs[[i]] <- GSEA_Msigdb_UI(ns(Msigdb_selection[i]), tab_name = Msigdb_selection[i])
            appendTab("Msigdb_Tabs",GSEA_other_UI(ns(Msigdb_selection[i]), tab_name = Msigdb_selection[i]))
            print("BB")
          }
          
          for (i in 1:length(Msigdb_selection)) { ## import subtabs server
            if(Msigdb_selection[i] != "other"){
              print("CC")
              category = strsplit(Msigdb_selection[i] ,"\\.")[[1]][1]
              subcategory = strsplit(Msigdb_selection[i] ,"\\.")[[1]][2]
              # subcategory = ifelse(is.na(subcategory), NULL, ifelse(subcategory == "TFT","TFT:GTRD", subcategory))
              if(is.na(subcategory)) {subcategory =NULL
              }else if(subcategory == "TFT") { subcategory == "TFT:GTRD"}
              GSEA_other_server_module(id = Msigdb_selection[i], gene_list = gene_list, organism_for_GSEA = organism_for_GSEA, annoSpecies_df = annoSpecies_df, 
                                      type = "Msigdb", Msigdb_category = category , Msigdb_subcategory = subcategory )
              print("ee")
            }else if(!is.null(Msigdb_selection2)){
              Msigdb_selection2_saved2 <<- Msigdb_selection2
              Msigdb_select2 <- Msigdb_selection2 %>% sapply(.,function(x) strsplit(x," \\(")[[1]][1]) %>% unname()
              GSEA_other_server_module(id = Msigdb_selection[i], gene_list = gene_list, organism_for_GSEA = organism_for_GSEA, annoSpecies_df = annoSpecies_df, 
                                      type = "Msigdb", Msigdb_category = NULL , Msigdb_subcategory = Msigdb_select2 )
              print("FF")
            }
          }
          
        }
      },
      once = T,ignoreInit=T , ignoreNULL=F)
    }
    
  )
}








