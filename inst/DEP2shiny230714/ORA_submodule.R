#* GO submodule UI/Server in ORA module ----
ORA_GO_UI <- function(id){
  ns = NS(id)
  tabPanel("GO",
           fluidRow(
             box(selectizeInput(ns("go_ont"),
                                "Ontology",
                                choices = c("BP", "CC", "MF", "ALL"),
                                selected = c("BP"), multiple = FALSE),
                 width = 2),
             box(numericInput(ns("go_p"),
                              "P value",
                              min = 0, max = 1, value = 0.05),
                 width = 2),
             box(numericInput(ns("go_padj"),
                              "adj. P value",
                              min = 0, max = 1, value = 0.05),
                 width = 2),
             box(numericInput(ns("go_qvalue"),
                              "q value",
                              min = 0, max = 1, value = 0.2),
                 width = 2),
             box(radioButtons(ns("go_color"),
                              "colorBy",
                              c("pvalue", "p.adjust"),
                              selected = "p.adjust"), width = 4),
             # infoBoxOutput("significantBox_for_go", width = 4)
             shinyBS::bsTooltip(ns("go_ont"), "One of BP, MF, and CC subontologies, or ALL for all three", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("go_p"), "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("go_padj"), "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("go_qvalue"), "qvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("go_color"), "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
           ),
           fluidRow(
             column(width = 7,
                    fluidRow(
                      box(checkboxInput(ns("go_simplify"), "removed redundancy of enriched GO terms", value = FALSE), width = 5),
                      # box(radioButtons("go_color",
                      #          "colorBy",
                      #          c("pvalue", "p.adjust"),
                      #          selected = "p.adjust"), width = 4),
                      infoBoxOutput(ns("significantBox_for_go"), width = 7)
                    ),
                    fluidRow(
                      box(title = "Result Table",
                          fluidRow(
                            column(width = 5,
                                   uiOutput(ns("downloadTable_go"))
                            )
                          ),
                          fluidRow(
                            column(width = 4,
                                   uiOutput(ns("downloadButton_go"))
                            )
                          ),
                          tags$br(),
                          # uiOutput(ns("downloadTable_go")),
                          # uiOutput(ns("downloadButton_go")),
                          DT::dataTableOutput(ns("Table")), width = 12)
                    ),
                    shinyBS::bsTooltip(ns("go_simplify"), "Whether simplify output by removing redundancy of enriched GO terms", "top", options = list(container = "body")),
                    shinyBS::bsTooltip(ns("Table"), "The table of the significant result", "top", options = list(container = "body")),
                    shinyBS::bsTooltip(ns("downloadTable_go"), "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
                    shinyBS::bsTooltip(ns("downloadButton_go"), "Click on it to download the result table", "right", options = list(container = "body"))
             ),
             column(width = 5,
                    tabBox(title = "Result Plots", width = 12,
                           tabPanel(title = "Bar plot",
                                    fluidRow(box(checkboxInput(ns("go_Select_bar"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("go_ShowCategory_bar"), "ShowCategory", min = 1, max = 100, value = 20),
                                          width = 4),
                                      box(numericInput(ns("go_wide_bar"), "Width", min = 1, max = 50, value = 10),
                                          width = 4),
                                      box(numericInput(ns("go_high_bar"), "Height", min = 1, max = 50, value = 7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      box(checkboxInput(ns("go_bar_if_Split_for_ont_ALL"),
                                                        "If splited by ontology",
                                                        value = FALSE), width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("go_barplot"), height = 500),
                                      downloadButton(ns('download_go_barplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("go_Select_bar"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_ShowCategory_bar"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_wide_bar"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_high_bar"), "Height of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_bar_if_Split_for_ont_ALL"), "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each BP, CC and MF", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Dot plot",
                                    fluidRow(box(checkboxInput(ns("go_Select_dot"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("go_ShowCategory_dot"), "ShowCategory", min = 1, max = 100, value = 20),
                                          width = 4),
                                      box(numericInput(ns("go_wide_dot"), "Width", min = 1, max = 50, value = 10),
                                          width = 4),
                                      box(numericInput(ns("go_high_dot"), "Height", min = 1, max = 50, value = 7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      box(checkboxInput(ns("go_dot_if_Split_for_ont_ALL"),
                                                        "If splited by ontology",
                                                        value = FALSE), width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("go_dotplot"), height = 500),
                                      downloadButton(ns('download_go_dotplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("go_Select_dot"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_ShowCategory_dot"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_wide_dot"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_high_dot"), "Height of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_dot_if_Split_for_ont_ALL"), "Whether split by ontology, act when [Ontology] is ALL. Note that: by default, the plot shows the top ShowCategory number of the significant table including BP,CC,MF result sorted by p.adjust values in ascending order together, and if it is selected , the plot shows the top ShowCategory number of each BP, CC and MF", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Dot opt",
                                    fluidRow(box(checkboxInput(ns("go_Select_dot_opt"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("go_ShowCategory_dot_opt"), "ShowCategory", min = 1, max = 100, value = 20),
                                          width = 4),
                                      box(numericInput(ns("go_wide_dot_opt"), "Width", min = 1, max = 50, value = 10),
                                          width = 4),
                                      box(numericInput(ns("go_high_dot_opt"), "Height", min = 1, max = 50, value = 5.7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("go_dotplot_opt"), height = 500),
                                      downloadButton(ns('download_go_dotplot_opt'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("go_Select_dot_opt"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_ShowCategory_dot_opt"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_wide_dot_opt"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_high_dot_opt"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Heatplot",
                                    fluidRow(box(checkboxInput(ns("go_Select_heat"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("go_ShowCategory_heat"), "ShowCategory", min = 1, max = 100, value = 30),
                                          width = 4),
                                      box(numericInput(ns("go_wide_heat"), "Width", min = 1, max = 50, value = 10),
                                          width = 4),
                                      box(numericInput(ns("go_high_heat"), "Height", min = 1, max = 50, value = 7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("go_heatplot"), height = 500),
                                      downloadButton(ns('download_go_heatplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("go_Select_heat"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_ShowCategory_heat"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_wide_heat"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_high_heat"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Cnetplot",
                                    fluidRow(box(checkboxInput(ns("go_Select_cnet"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("go_ShowCategory_cnet"), "ShowCategory", min = 1, max = 100, value = 5),
                                          width = 4),
                                      box(numericInput(ns("go_wide_cnet"), "Width", min = 1, max = 50, value = 10),
                                          width = 4),
                                      box(numericInput(ns("go_high_cnet"), "Height", min = 1, max = 50, value = 6),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      box(checkboxInput(ns("go_circular_cnet"),
                                                        "Circular",
                                                        value = TRUE),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("go_cnetplot"), height = 500),
                                      downloadButton(ns('download_go_cnetplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("go_Select_cnet"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_ShowCategory_cnet"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_wide_cnet"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_high_cnet"), "Height of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_circular_cnet"), "whether using circular layout", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Emaplot",
                                    fluidRow(box(checkboxInput(ns("go_Select_ema"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("go_ShowCategory_ema"), "ShowCategory", min = 1, max = 100, value = 30),
                                          width = 4),
                                      box(numericInput(ns("go_wide_ema"), "Width", min = 1, max = 50, value = 11),
                                          width = 4),
                                      box(numericInput(ns("go_high_ema"), "Height", min = 1, max = 50, value = 10),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("go_emaplot"), height = 500),
                                      downloadButton(ns('download_go_emaplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("go_Select_ema"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_ShowCategory_ema"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_wide_ema"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("go_high_ema"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           ## GOplot and GOgraph are not usable for all omics analysis result, hide it.
                           # tabPanel(title = "Goplot",
                           #          fluidRow(
                           #            box(numericInput(ns("go_ShowCategory_go"), "ShowCategory", min = 1, max = 100, value = 10),
                           #                width = 4),
                           #            box(numericInput(ns("go_wide_go"), "Width", min = 1, max = 50, value = 11),
                           #                width = 4),
                           #            box(numericInput(ns("go_high_go"), "Height", min = 1, max = 50, value = 10),
                           #                width = 4)
                           #          ),
                           #          fluidRow(
                           #            box(checkboxInput(ns("go_circular_go"),
                           #                              "Circular",
                           #                              value = TRUE),
                           #                width = 4)
                           #          ),
                           #          fluidRow(
                           #            plotOutput(ns("go_goplot"), height = 500),
                           #            downloadButton(ns('download_go_goplot'), 'Save')
                           #          ),
                           #          shinyBS::bsTooltip(ns("go_ShowCategory_go"), "Number of categories to show", "top", options = list(container = "body")),
                           #          shinyBS::bsTooltip(ns("go_wide_go"), "Width of the figure to export", "top", options = list(container = "body")),
                           #          shinyBS::bsTooltip(ns("go_high_go"), "Height of the figure to export", "top", options = list(container = "body")),
                           #          shinyBS::bsTooltip(ns("go_circular_go"), "whether using circular layout", "top", options = list(container = "body"))
                           # ),
                           # tabPanel(title = "GOgraph",
                           #          fluidRow(
                           #            box(numericInput(ns("go_ShowCategory_GOgraph"), "FirstSigNodes", min = 1, max = 100, value = 10),
                           #                width = 4),
                           #            box(numericInput(ns("go_wide_GOgraph"), "Width", min = 1, max = 50, value = 7),
                           #                width = 4),
                           #            box(numericInput(ns("go_high_GOgraph"), "Height", min = 1, max = 50, value = 7),
                           #                width = 4)
                           #          ),
                           #          fluidRow(actionButton(ns("plot_for_GOgraphplot"), "Plot", width = "10%", icon = icon("caret-right"))),#step-forward
                           #          fluidRow(
                           #            plotOutput(ns("go_GOgraphplot"), height = 500),
                           #            downloadButton(ns('download_go_GOgraphplot'), 'Save')
                           #          ),
                           #          shinyBS::bsTooltip(ns("go_ShowCategory_GOgraph"), "Number of significant nodes (retangle nodes in the graph)", "top", options = list(container = "body")),
                           #          shinyBS::bsTooltip(ns("go_wide_GOgraph"), "Width of the figure to export", "top", options = list(container = "body")),
                           #          shinyBS::bsTooltip(ns("go_high_GOgraph"), "Height of the figure to export", "top", options = list(container = "body"))
                           # )
                    )

             )
           )

  )
}

ORA_GO_server_module <- function(id, gene_id, organism_for_ORA, annoSpecies_df) {
  moduleServer(
    id = id,
    function(input, output, session) {

      # print(head(gene_df))
      print(id)
      print(organism_for_ORA)
      # print(annoSpecies_df)
      print(gene_id)
      print("b1")

      pkg_for_ORA <- annoSpecies_df$pkg[annoSpecies_df$species == organism_for_ORA]

      output$downloadTable_go <- renderUI({
        selectizeInput(session$ns("dataset_for_go"),
                       "Choose a dataset to save" ,
                       c("full_results","significant_results"
                       ))
      })

      output$downloadButton_go <- renderUI({
        downloadButton(session$ns('downloadgo'), 'Save table')
      })

      reat <- reactive({
        try(goAnalysis(gene_id = gene_id,
                       # df_with_lg2fc = F,
                       organism = organism_for_ORA, species_df = annoSpecies_df),
            silent = TRUE)
      })

      res <- reactive({
        try(giveGO_res_and_table(reat = reat(), ont = input$go_ont, pCutoff = input$go_p, p.adj.cutoff = input$go_padj, q.cutoff = input$go_qvalue, simplify = input$go_simplify), silent = TRUE)
      })

      res_select <- reactive({
        res <- res()
        if(is.null(input$Table_rows_selected)){
          res$sig_res <- NULL
        }else{
          res$sig_res <- res$sig_res %>% filter(.,1:nrow(.) %in% input$Table_rows_selected)
        }
        return(res)
      })

      output$Table <- DT::renderDataTable({

        shiny::validate(need(require(pkg_for_ORA, character.only = TRUE), message = paste0("The package ", pkg_for_ORA, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_ORA, "')")))
        shiny::validate( need( class(reat()) != "try-error", "No genes can be mapped to ENTREZID, and can not do the GO analysis"))
        res <- res()
        validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
        res()$sig_table
        # DT::datatable(res()$sig_table, filter = 'top', options = list( autoWidth = F,scrollX = TRUE, selection = list(selected = c(1)) ) )
      },
      options = list(pageLength = 25, scrollX = T),
      filter = 'top',
      selection = list(selected = c(1))
      )

      # output$table <- DT::renderDataTable({
      #   table()
      # }, options = list(pageLength = 25, scrollX = T),
      # selection = list(selected = c(1)))


      output$significantBox_for_go <- renderInfoBox({
        shiny::validate(need(require(pkg_for_ORA, character.only = TRUE), message = paste0("The package ", pkg_for_ORA, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_ORA, "')")))

        validate( need(class(reat()) != "try-error", "No genes can be mapped to ENTREZID, and can not do the GO analysis"))

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
          info_box_go <- infoBox("Significant terms",
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
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$go_Select_bar){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }
          my_barplot(res = res, ShowCategory = input$go_ShowCategory_bar, color = input$go_color, ont = input$go_ont, Split = input$go_bar_if_Split_for_ont_ALL)
        })
      })

      go_dotplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$go_Select_dot){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }
          my_dotplot(res = res, ShowCategory = input$go_ShowCategory_dot, color = input$go_color, ont = input$go_ont, Split = input$go_dot_if_Split_for_ont_ALL)
        })
      })

      go_dotplot_opt_input <- reactive({
        validate(
          need(input$go_ont != "ALL", "Please go to the panel : Dot plot ")
        )
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$go_Select_dot_opt){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }
          my_dotplot_opt(res = res, color = input$go_color, size = "Count", title = "", decreasing = TRUE, ShowCategory = input$go_ShowCategory_dot_opt)
        })
      })

      go_heatplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$go_Select_heat){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }
          my_heatplot(res = res, ShowCategory = input$go_ShowCategory_heat, df_with_lg2fc = F, ont = input$go_ont)
        })
      })

      go_cnetplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$go_Select_cnet){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }
          my_cnetplot(res = res, ShowCategory = input$go_ShowCategory_cnet, circular = input$go_circular_cnet, colorEdge = TRUE, df_with_lg2fc = F, ont = input$go_ont)
        })
      })

      go_emaplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$go_Select_ema){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }
          my_emaplot(res = res, ShowCategory = input$go_ShowCategory_ema, color = input$go_color, layout = "kk", ont = input$go_ont)
        })
      })

      go_goplot_input <- reactive({
        validate(
          need(input$go_ont != "ALL", "Ontology ALL: can not plot for goplot")
        )
        withProgress(message = 'Plotting', value = 0.66, {
          my_goplot(res = res(), ShowCategory = input$go_ShowCategory_go, color = input$go_color, ont = input$go_ont, Layout = "kk",circular = input$go_circular_go)
        })
      })

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


      observeEvent(input$plot_for_GOgraphplot, {
        go_GOgraphplot_input <- reactive({
          validate(
            need(input$go_ont != "ALL", "Ontology ALL: can not plot for plotGOgraph")
          )
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


    }
  )
}


#* other submodule UI/Server in ORA module, for KEGG and Reactome. ----
ORA_other_UI <- function(id, tab_name = id){
  ns = NS(id)
  tabPanel(tab_name,
           fluidRow(
             box(numericInput(ns("p"),
                              "P value",
                              min = 0, max = 1, value = 0.05),
                 width = 2),
             box(numericInput(ns("padj"),
                              "adj. P value",
                              min = 0, max = 1, value = 0.05),
                 width = 2),
             box(numericInput(ns("qvalue"),
                              "q value",
                              min = 0, max = 1, value = 0.2),
                 width = 2),
             infoBoxOutput(ns("significantBox_for"), width = 4),
             box(radioButtons(ns("color"),
                              "colorBy",
                              c("pvalue", "p.adjust"),
                              selected = "p.adjust"), width = 2),
             shinyBS::bsTooltip(ns("p"), "Pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("padj"), "Adjusted pvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("qvalue"), "qvalue cutoff on enrichment tests to report", "top", options = list(container = "body")),
             shinyBS::bsTooltip(ns("color"), "Color filled based on, one of pvalue, and p.adjust, and it is used in all output graphics that color is filled based on pvalue or p.adjust", "top", options = list(container = "body"))
           ),
           fluidRow(
             column(width = 7,
                    fluidRow(
                      box(title = "Top Table",
                          fluidRow(
                            column(width = 5,
                                   uiOutput(ns("downloadTable"))
                            )
                          ),
                          fluidRow(
                            column(width = 4,
                                   uiOutput(ns("downloadButton"))
                            )
                          ),
                          tags$br(),
                          # uiOutput(ns("downloadTable")),
                          # uiOutput(ns("downloadButton")),
                          DT::dataTableOutput(ns("Table")), width = 12)
                    ),
                    shinyBS::bsTooltip(ns("Table"), "The table of the significant result", "top", options = list(container = "body")),
                    shinyBS::bsTooltip(ns("downloadTable"), "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
                    shinyBS::bsTooltip(ns("downloadButton"), "Click on it to download the result table", "right", options = list(container = "body"))
             ),
             column(width = 5,
                    tabBox(title = "Result Plots", width = 12,
                           tabPanel(title = "Bar plot",
                                    fluidRow(box(checkboxInput(ns("Select_bar"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("ShowCategory_bar"), "ShowCategory", min = 1, max = 100, value = 20),
                                          width = 4),
                                      box(numericInput(ns("wide_bar"), "Width", min = 1, max = 50, value = 10),
                                          width = 4),
                                      box(numericInput(ns("high_bar"), "Height", min = 1, max = 50, value = 7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("barplot"), height = 500),
                                      downloadButton(ns('download_barplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("Select_bar"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("ShowCategory_bar"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("wide_bar"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("high_bar"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Dot plot",
                                    fluidRow(box(checkboxInput(ns("Select_dot"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("ShowCategory_dot"), "ShowCategory", min = 1, max = 100, value = 20),
                                          width = 4),
                                      box(numericInput(ns("wide_dot"), "Width", min = 1, max = 50, value = 10),
                                          width = 4),
                                      box(numericInput(ns("high_dot"), "Height", min = 1, max = 50, value = 7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("dotplot"), height = 500),
                                      downloadButton(ns('download_dotplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("Select_dot"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("ShowCategory_dot"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("wide_dot"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("high_dot"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Dot opt",
                                    fluidRow(box(checkboxInput(ns("Select_dot_opt"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("ShowCategory_dot_opt"), "ShowCategory", min = 1, max = 100, value = 20),
                                          width = 4),
                                      box(numericInput(ns("wide_dot_opt"), "Width", min = 1, max = 50, value = 10),
                                          width = 4),
                                      box(numericInput(ns("high_dot_opt"), "Height", min = 1, max = 50, value = 5.7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("dotplot_opt"), height = 500),
                                      downloadButton(ns('download_dotplot_opt'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("Select_dot_opt"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("ShowCategory_dot_opt"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("wide_dot_opt"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("high_dot_opt"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Heatplot",
                                    fluidRow(box(checkboxInput(ns("Select_heat"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("ShowCategory_heat"), "ShowCategory", min = 1, max = 100, value = 30),
                                          width = 4),
                                      box(numericInput(ns("wide_heat"), "Width", min = 1, max = 50, value = 10),
                                          width = 4),
                                      box(numericInput(ns("high_heat"), "Height", min = 1, max = 50, value = 7),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("heatplot"), height = 500),
                                      downloadButton(ns('download_heatplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("Select_heat"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("ShowCategory_heat"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("wide_heat"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("high_heat"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Cnetplot",
                                    fluidRow(box(checkboxInput(ns("Select_cnet"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("ShowCategory_cnet"), "ShowCategory", min = 1, max = 100, value = 5),
                                          width = 4),
                                      box(numericInput(ns("wide_cnet"), "Width", min = 1, max = 50, value = 10),
                                          width = 4),
                                      box(numericInput(ns("high_cnet"), "Height", min = 1, max = 50, value = 6),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      box(checkboxInput(ns("circular_cnet"),
                                                        "Circular",
                                                        value = TRUE),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("cnetplot"), height = 500),
                                      downloadButton(ns('download_cnetplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("Select_cnet"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("ShowCategory_cnet"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("wide_cnet"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("high_cnet"), "Height of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("circular_cnet"), "whether using circular layout", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Emaplot",
                                    fluidRow(box(checkboxInput(ns("Select_ema"), label = "Selected genesets"), width = 6)),
                                    fluidRow(
                                      box(numericInput(ns("ShowCategory_ema"), "ShowCategory", min = 1, max = 100, value = 30),
                                          width = 4),
                                      box(numericInput(ns("wide_ema"), "Width", min = 1, max = 50, value = 11),
                                          width = 4),
                                      box(numericInput(ns("high_ema"), "Height", min = 1, max = 50, value = 10),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("emaplot"), height = 500),
                                      downloadButton(ns('download_emaplot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("Select_ema"), "Display the selected genesets in Result Table or top enriched terms(ordered by p)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("ShowCategory_ema"), "Number of categories to show", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("wide_ema"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("high_ema"), "Height of the figure to export", "top", options = list(container = "body"))
                           )
                    )

             )
           )

  )
}

ORA_other_server_module <- function(id,
                                    # gene_df,
                                    gene_id,
                                    organism_for_ORA, annoSpecies_df,
                                    type, Msigdb_category = NULL , Msigdb_subcategory =NULL) {
  moduleServer(
    id,
    function(input, output, session,id=id) {

      print(organism_for_ORA)
      pkg_for_ORA <- annoSpecies_df$pkg[annoSpecies_df$species == organism_for_ORA]

      print(type)
      print(pkg_for_ORA)

      ### Interactive UI functions ### ------------------------------------------
      output$downloadTable <- renderUI({
        selectizeInput(session$ns("dataset_for"),
                       "Choose a dataset to save" ,
                       c("full_results","significant_results"
                       ))
      })

      output$downloadButton <- renderUI({
        downloadButton(session$ns('downloadkegg'), 'Save table')
      })


      if(type == "kegg"){
        reat <- reactive({try(keggAnalysis(gene_id = gene_id,
                                           # df = gene_df,
                                           organism = organism_for_ORA, df_with_lg2fc = F, species_df = annoSpecies_df),
                              silent = TRUE)})

      }else if(type == "reactome"){
        reat <- reactive({
          try(reactAnalysis(gene_id = gene_id,
                            # df = gene_df,
                            organism = organism_for_ORA, df_with_lg2fc = F, species_df = annoSpecies_df),
              silent = TRUE)
        })

      }else if(type == "Msigdb"){
        reat <- reactive({
          try(msigdb_Analysis(gene_id = gene_id,
                              # df = gene_df,
                              organism = organism_for_ORA, df_with_lg2fc = F, species_df = annoSpecies_df,
                              category = Msigdb_category, subcategory = Msigdb_subcategory), silent = TRUE)
        })
      }

      res <- reactive({
        try(give_res_and_table(reat = reat(), pCutoff = input$p, p.adj.cutoff = input$padj, q.cutoff = input$qvalue),
            silent = TRUE)
      })

      res_select <- reactive({
        res <- res()
        if(is.null(input$Table_rows_selected)){
          res$sig_res <- NULL
        }else{
          res$sig_res <- res$sig_res %>% filter(.,1:nrow(.) %in% input$Table_rows_selected)
        }
        return(res)
      })

      # output$Table <- DT::renderDataTable({
      #   shiny::validate(
      #     need(length(gene_id) != 0, message = "No genes meet your requirements, and can not do the KEGG analysis")
      #   )
      #
      #   shiny::validate(need(require(pkg_for_ORA, character.only = TRUE), message = paste0("The package ", pkg_for_ORA, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_ORA, "')")))
      #
      #   validate( need(class(reat()) != "try-error", "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
      #   # shiny::validate( need(class(reat) != "try-error", message = "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
      #   DT::datatable(res()$sig_table, filter = 'top', options = list( autoWidth = F,scrollX = TRUE, selection = list(selected = c(1)) ) )
      # })
      output$Table <- DT::renderDataTable({

        shiny::validate(need(require(pkg_for_ORA, character.only = TRUE), message = paste0("The package ", pkg_for_ORA, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_ORA, "')")))
        shiny::validate( need( class(reat()) != "try-error", "No genes can be mapped to ENTREZID, and can not do the GO analysis"))
        res <- res()
        validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
        res()$sig_table
        # DT::datatable(res()$sig_table, filter = 'top', options = list( autoWidth = F,scrollX = TRUE, selection = list(selected = c(1)) ) )
      },
      options = list(pageLength = 25, scrollX = T),
      filter = 'top',
      selection = list(selected = c(1))
      )

      output$significantBox_for <- renderInfoBox({
        shiny::validate(
          need(length(gene_id) != 0, message = "No genes meet your requirements, and can not do the KEGG analysis")
        )

        shiny::validate(need(require(pkg_for_ORA, character.only = TRUE), message = paste0("The package ", pkg_for_ORA, " is not installed/available. Try installing it with BiocManager::install('", pkg_for_ORA, "')")))

        validate( need(class(reat()) != "try-error", "No genes can be mapped to ENTREZID, and can not do the KEGG analysis"))
        num_total <- res()$all_table %>%
          nrow()
        num_signif <- res()$sig_table %>%
          nrow()
        frac <- num_signif / num_total

        if(frac == 0) {
          info_box <- infoBox("Significant terms",
                                   paste0(num_signif,
                                          " out of ",
                                          num_total),
                                   "No terms enriched",
                                   icon = icon("thumbs-down", lib = "glyphicon"),
                                   color = "red",
                                   width = 4)
        }
        if(!frac == 0) {
          info_box <- infoBox("Significant terms",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              paste0(signif(frac * 100, digits = 3),
                                     "% of terms enriched"),
                              icon = icon("thumbs-up", lib = "glyphicon"),
                              color = "green",
                              width = 4)
        }
        info_box
      })

      barplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$Select_bar){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }else{
            res <- res()
          }
          barplot(res$sig_res, showCategory = input$ShowCategory_bar, color = input$color)
        })
      })

      dotplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$Select_dot){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }
          dotplot(res$sig_res, showCategory = input$ShowCategory_dot, color = input$color)
        })
      })

      dotplot_opt_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$Select_dot_opt){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }
          my_dotplot_opt(res = res, color = input$color, size = "Count", title = "", decreasing = TRUE, ShowCategory = input$ShowCategory_dot_opt)
        })
      })

      heatplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$Select_heat){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }
          heatplot_for(res = res, ShowCategory = input$ShowCategory_heat, df_with_lg2fc = F)
        })
      })

      cnetplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$Select_cnet){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }
          cnetplot_for(res = res, ShowCategory = input$ShowCategory_cnet, circular = input$circular_cnet, colorEdge = TRUE, df_with_lg2fc = F)
        })
      })

      emaplot_input <- reactive({
        withProgress(message = 'Plotting', value = 0.66, {
          res <- res()
          validate(need(nrow(res$sig_res) > 0, "no term is enriched"))
          if(input$Select_ema){
            validate(need(!is.null(input$Table_rows_selected), "please select at least one enriched term in result table"))
            res <- res_select()
          }
          emaplot_for(res = res, ShowCategory = input$ShowCategory_ema, color = input$color, layout = "kk")
        })
      })

      output$barplot <- renderPlot({
        barplot_input()
      })

      output$dotplot <- renderPlot({
        dotplot_input()
      })

      output$dotplot_opt <- renderPlot({
        dotplot_opt_input()
      })

      output$heatplot <- renderPlot({
        heatplot_input()
      })

      output$cnetplot <- renderPlot({
        cnetplot_input()
      })

      output$emaplot <- renderPlot({
        emaplot_input()
      })

      ### Download objects and functions ### ------------------------------------
      datasetInput_for <- reactive({
        table_for = res()
        switch(input$dataset_for,
               "full_results" = res()$all_table,
               "significant_results" = res()$sig_table)
      })

      output$downloadkegg <- downloadHandler(
        filename = function() { paste(input$dataset_for, ".txt", sep = "") },
        content = function(file) {
          write.table(datasetInput_for(),
                      file,
                      col.names = TRUE,
                      row.names = FALSE,
                      sep ="\t") }
      )

      output$download_barplot <- downloadHandler(
        filename = 'barplot.pdf',
        content = function(file) {
          pdf(file, width = input$wide_bar, height = input$high_bar)
          print(barplot_input())
          dev.off()
        }
      )

      output$download_dotplot <- downloadHandler(
        filename = 'dotplot.pdf',
        content = function(file) {
          pdf(file, width = input$wide_dot, height = input$high_dot)
          print(dotplot_input())
          dev.off()
        }
      )

      output$download_dotplot_opt <- downloadHandler(
        filename = 'dotplot_opt.pdf',
        content = function(file) {
          pdf(file, width = input$wide_dot_opt, height = input$high_dot_opt)
          print(dotplot_opt_input())
          dev.off()
        }
      )

      output$download_heatplot <- downloadHandler(
        filename = 'heatplot.pdf',
        content = function(file) {
          pdf(file, width = input$wide_heat, height = input$high_heat)
          print(heatplot_input())
          dev.off()
        }
      )

      output$download_cnetplot <- downloadHandler(
        filename = 'cnetplot.pdf',
        content = function(file) {
          pdf(file, width = input$wide_cnet, height = input$high_cnet)
          print(cnetplot_input())
          dev.off()
        }
      )

      output$download_emaplot <- downloadHandler(
        filename = 'emaplot.pdf',
        content = function(file) {
          pdf(file, width = input$wide_ema, height = input$high_ema)
          print(emaplot_input())
          dev.off()
        }
      )



    }
  )}


#* Msigdb submodule UI/Server in ORA module ----
ORA_Msigdb_UI <- function(id, tab_name = "Msigdb"){
  ns = NS(id)

  tabPanel(
    tab_name,
    uiOutput(ns("Msigdb_tabs"))
  )
}

ORA_Msigdb_server_module <- function(id, gene_df, gene_id, organism_for_ORA, annoSpecies_df, Msigdb_selection = NULL, Msigdb_selection2 = NULL) {
  # id = paste0(id,"_",suffix)

  moduleServer(
    id = id,
    function(input, output, session) {
      ns = session$ns

      #** Msigdb_tabs render in msigdb ----
      output$Msigdb_tabs <- renderUI({

        tabsetPanel(
          id = session$ns("Msigdb_Tabs")
        )
      })
      observeEvent(input$Msigdb_Tabs,{
        if(length(Msigdb_selection) > 0){
          for (i in 1:length(Msigdb_selection)) { ## subtabs UI

            appendTab("Msigdb_Tabs",ORA_other_UI(ns(Msigdb_selection[i]), tab_name = Msigdb_selection[i]))
          }
          for (i in 1:length(Msigdb_selection)) { ## subtabs server
            if(Msigdb_selection[i] != "other"){

              category = strsplit(Msigdb_selection[i] ,"\\.")[[1]][1]
              subcategory = strsplit(Msigdb_selection[i] ,"\\.")[[1]][2]
              # subcategory = ifelse(is.na(subcategory), NULL, ifelse(subcategory == "TFT","TFT:GTRD", subcategory))
              if(is.na(subcategory)) {subcategory =NULL
              }else if(subcategory == "TFT") { subcategory == "TFT:GTRD"}
              ORA_other_server_module(id = Msigdb_selection[i], gene_id = gene_id, organism_for_ORA = organism_for_ORA, annoSpecies_df = DEP2:::annoSpecies_df(),
                                      type = "Msigdb", Msigdb_category = category , Msigdb_subcategory =subcategory )

            }else if(!is.null(Msigdb_selection2)){
              Msigdb_select2 <- Msigdb_selection2 %>% sapply(.,function(x) strsplit(x," \\(")[[1]][1]) %>% unname()
              ORA_other_server_module(id = Msigdb_selection[i], gene_id = gene_id, organism_for_ORA = organism_for_ORA, annoSpecies_df = DEP2:::annoSpecies_df(),
                                      type = "Msigdb", Msigdb_category = NULL , Msigdb_subcategory = Msigdb_select2 )

            }
          }
        }
      },
      once = T,ignoreInit=T , ignoreNULL=F)



      # print("GG")
      # for (i in 1:length(Msigdb_selection)) {
      #   ORA_other_server_module(Msigdb_selection[i])
      # }

    }
  )}




