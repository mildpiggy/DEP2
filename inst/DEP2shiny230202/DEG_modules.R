DEG_sidebar_mod <- function(id,labelname = "DEP-pg_sidabar"){
  ns = NS(id)
  tagList(
    sidebarPanel(
      width = 3,
      h3(labelname),
      tags$hr(),
      bsCollapse(
        bsCollapsePanel("Files",
                        style = "primary",
                        fluidRow(
                          column(width = 12,
                                 radioButtons(ns("uploadmode"),"Upload mode",
                                              choices = c("From result table" = "fromTable",
                                                          "From saved log file" = "fromLog"),
                                              selected = "fromTable")
                          )
                        ),
                        hr(),

                        # file input options, from Count table
                        conditionalPanel(
                          condition = paste0("input['",ns("uploadmode"),"']","== 'fromTable'"), ## condition in mod
                          fluidRow(
                            column(width = 9,
                                   fileInput(ns('file1'),width = "300px",
                                             'Countmatrix.txt',
                                             accept=c('text/csv',
                                                      'text/comma-separated-values,text/plain',
                                                      '.csv',
                                                      '.gz',
                                                      '.tsv')),
                                   conditionalPanel(
                                     condition = paste0("input['",ns("anno"),"']","== 'expdesign'"),
                                     fileInput(ns('file2'),width = "300px",
                                               'ExperimentalDesign.txt',
                                               accept=c('text/csv',
                                                        'text/comma-separated-values,text/plain',
                                                        '.csv','.gz','.tsv')),
                                     uiOutput(ns("for_loadcount")) ## only to force counts file is read in immediately
                                   )

                            ),
                            column(width = 1,h4(),
                                   br(),
                                   actionButton(ns("help_format_RNA"),
                                                label = "", icon = icon("question-circle"),
                                                style = "color: #f6f6f6; background-color: #2c3b41; border-color: #2c3b41"
                                   ),
                                   shinyBS::bsTooltip(
                                     ns("help_format_RNA"),
                                     "How to provide your input data",
                                     "top",
                                     options = list(container = "body")
                                   ))
                          )
                        ),

                        # file input options, from log file
                        conditionalPanel(
                          condition = paste0("input['",ns("uploadmode"),"']","== 'fromLog'"),
                          fluidRow(
                            column(width = 9,
                                   fileInput(ns('file_log'),width = "300px",
                                             'The saved ".logrds" file',
                                             accept=c('.logrds')),
                            )
                          )
                        ),
                        hr(),

                        ## how to generate experiment design
                        radioButtons(ns("anno"),
                                     "Sample annotation",
                                     choices = list("Parse from columns" = "columns",
                                                    "Use Experimental Design" = "expdesign"),
                                     selected = "columns")
        ),
        bsCollapsePanel("ID Transformation", style = "primary",
                        div(id = ns("transid_opts"),
                            checkboxInput(ns("transid_for_RNAseq"), "ID transformation", value = FALSE),
                            conditionalPanel(paste0("input['",ns("transid_for_RNAseq"),"']"),
                                             selectInput(ns("species_select"),
                                                         label = "Select species",
                                                         choices = DEP2:::annoSpecies_df()$species, selected = ""
                                             )),
                            conditionalPanel(paste0("input['",ns("transid_for_RNAseq"),"']"),
                                             checkboxInput(ns("set_no_map_to_rowname"), "Keep no-mapping to rowname", value = TRUE)),
                            conditionalPanel(paste0("input['",ns("transid_for_RNAseq"),"']"),
                                             uiOutput(ns("ui_idtype")))
                        )

                        # uiOutput(ns("ui_set_no_map_to_rowname")),
                        # uiOutput(ns("ui_selectspecies")),
                        # uiOutput(ns("ui_idtype"))
        ),
        # uiOutput(ns("ui_ID_transformation")),
        bsCollapsePanel("RNAseq_settings", style = "primary",
                        uiOutput(ns("dds_design")),
                        uiOutput(ns("choose_expfac")),
                        uiOutput(ns("whether_LRT")),
                        uiOutput(ns("reduced_model")),
                        # uiOutput("choose_expfac1"),
                        # uiOutput("choose_expfac2"),
                        uiOutput(ns("control")),
                        uiOutput(ns("test_manual_for_RNAseq")),
                        numericInput(ns("filter_rowsum"), "Threshold row sums", min = 0, max = 10000, value = 20),
                        numericInput(ns("filter_missnum"), "Threshold of missing number", min = 0, max = 10000, value = 2)
                        # numericInput(ns("thr"), "Allowed max.miss.num at least one condition", min = 0, max = 20, value = 3)
        ),
        bsCollapsePanel("Result options", style = "primary",
                        checkboxInput(ns("independent_filtering"), "Independent filtering", value = TRUE),
                        checkboxInput(ns("Shrink_lfc"), "Shrink lfc", value = FALSE),
                        checkboxInput(ns("Use_IHW"), "Use IHW", value = FALSE)
        ),
        open = "Files"
      ),

      uiOutput(ns("ui_nrcores")),
      # conditionalPanel(
      #   paste0(
      #     "typeof ",paste0("!input['",ns("choose_expfac"),"']")," !== 'undefined' && ",
      #     paste0("!input['",ns("choose_expfac"),"']"),".length > 0"),
      #   actionButton(ns("analyze_for_DERNAseq"), "Analyze")),
      actionButton(ns("analyze_for_DERNAseq"), "Analyze"),
      # tags$hr(),

      tags$hr(),
      tags$style(type="text/css", "#download_DERNAseq {background-color:white;color: black;font-family: Source Sans Pro}"),
      uiOutput(ns("downloadTable_for_DERNAseq")),
      uiOutput(ns("downloadButton_for_DERNAseq")),
      uiOutput(ns("downloadButton_for_log")),
      shinyBS::bsTooltip(ns("dds_design"), "Select your design, and it can be multi-factor", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("choose_expfac"), "Choose the contrast to build upon", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("whether_LRT"), "I can perform a LRT test on the chosen factor, choose whether perform a LRT test", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("reduced_model"), "Select the reduced model. Note: 1) the full model is the same as your selected design; 2) if the number of your selected design is 1, select 1", "right", options = list(container = "body")),
      # shinyBS::bsTooltip("choose_expfac1", "Select the name of the numerator level for the fold change", "right", options = list(container = "body")),
      # shinyBS::bsTooltip("choose_expfac2", "Select the name of the denominator level for the fold change (must be different from the numerator)", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("control"), "Choose your control condition", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("test_manual_for_RNAseq"), "Choose the contrasts that you want to test", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("filter_rowsum"), "Threshold on the row sums of the counts", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("filter_missnum"), "Threshold of row missing number (0 number)", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("transid_for_RNAseq"), "Whether transform your id to gene symbol", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("set_no_map_to_rowname"), "Whether keep the rows that map to a gene symbol", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("ui_selectspecies"), "Select the species", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("ui_idtype"), "Select your id type", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("independent_filtering"), "whether independent filtering should be applied automatically", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("Shrink_lfc"), "whether shrink the log fold change for the contrast of interest", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("Use_IHW"), "Whether use Independent Hypothesis Weighting (IHW) as a filtering function", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("ui_nrcores"), "Choose how many cores to use for computing", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("analyze_for_DERNAseq"), "Click on it to analyze", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("downloadTable_for_DERNAseq"), "Choose a dataset to save, and here we offer two forms of datasets for downloading including full_results.txt, significant_results.txt", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("downloadButton_for_DERNAseq"), "Click on it to download the result table", "right", options = list(container = "body"))

    )
  )
}

DEG_body_mod <- function(id){
  ns = NS(id)
  tagList(
    mainPanel(
      fluidRow(
        box(numericInput(ns("p"),
                         "adj. P value",
                         min = 0.0001, max = 0.1, value = 0.05),
            width = 2),
        box(numericInput(ns("lfc"),
                         "Log2 fold change",
                         min = 0, max = 10, value = 1),
            width = 2),
        infoBoxOutput(ns("significantBox_for_DERNAseq")),
        box(radioButtons(ns("pres_for_RNAseq"),
                         "Data presentation",
                         c("contrast", "centered"),
                         selected = "centered"),
            width = 2),
        box(radioButtons(inputId = ns("contrasts"),
                         label = "Contrasts",
                         choices = c("control", "all", "manual"),
                         selected = "control", inline = TRUE),
            width = 2),
        shinyBS::bsTooltip(ns("p"), "Set the false discovery rate threshold", "top", options = list(container = "body")),
        shinyBS::bsTooltip(ns("lfc"), "Set the log2 fold change threshold", "top", options = list(container = "body")),
        shinyBS::bsTooltip(ns("pres_for_RNAseq"), "The type of data scaling used for heatmap plotting. Either the raw or the centered log2-intensity (centered)", "top", options = list(container = "body")),
        shinyBS::bsTooltip(ns("contrasts"), "The type of contrasts that will be tested. This can be all possible pairwise comparisons (all), limited to the comparisons versus the control (control), or manually defined contrasts (manual)", "top", options = list(container = "body"))
      ),
      fluidRow(
        column(width = 7,
               box(title = "Top Table",
                   DT::dataTableOutput(ns("table_for_DERNAseq")), width = 12),
               shinyBS::bsTooltip(ns("table_for_DERNAseq"), "Table of the significant results", "top", options = list(container = "body"))
        ),
        column(width = 5,
               tabBox(title = "Result Plots", width = 12,
                      tabPanel(title = "Diagno dispests",
                               fluidRow(
                                 box(numericInput(ns("Diagno_Width_for_RNAseq"),
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                                 box(numericInput(ns("Diagno_Height_for_RNAseq"),
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                               ),
                               fluidRow(
                                 plotOutput(ns("diagno_dispests_for_DERNAseq"), height = 600),
                                 downloadButton(ns('download_diagno_dispests_for_DERNAseq'), 'Save')
                               ),
                               shinyBS::bsTooltip(ns("Diagno_Width_for_RNAseq"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Diagno_Height_for_RNAseq"), "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Pca plot",
                               fluidRow(
                                 box(selectizeInput(ns("Indicate_for_RNAseq"),
                                                    "Color and shape",
                                                    choices = c("condition", "replicate", "Condition", "Replicate"),
                                                    selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                 box(checkboxInput(ns("if_square_for_RNAseq"),
                                                   "if square",
                                                   value = FALSE),
                                     width = 3)),
                               fluidRow(
                                 box(uiOutput(ns("pca_top_n_for_RNAseq")),width = 9)
                               ),

                               fluidRow(
                                 box(numericInput(ns("pca_Width_for_RNAsqe"),
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                                 box(numericInput(ns("pca_Height_for_RNAseq"),
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                               ),
                               fluidRow(
                                 plotOutput(ns("pca_for_RNAseq"), height = 600),
                                 downloadButton(ns('downloadPca_for_RNAseq'), 'Save')
                               ),
                               shinyBS::bsTooltip(ns("Indicate_for_RNAseq"), "Set the color, shape and facet_wrap of the plot", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("if_square_for_RNAseq"), "Whether x limit is equal to y limit", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("pca_top_n_for_RNAseq"), "Set the number of top variable genes to consider", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("pca_Width_for_RNAsqe"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("pca_Height_for_RNAseq"), "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "UMAP plot",
                               fluidRow(
                                 box(selectizeInput(ns("umap_Indicate_for_RNAseq"),
                                                    "Color and shape",
                                                    choices = c("condition", "replicate", "Condition", "Replicate"),
                                                    selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                 box(checkboxInput(ns("umap_if_square_for_RNAseq"),
                                                   "if square",
                                                   value = FALSE),
                                     width = 6)),
                               fluidRow(
                                 box(uiOutput(ns("umap_n_neighbors_for_RNAseq")),width = 10)
                               ),
                               fluidRow(
                                 box(numericInput(ns("umap_Width_for_RNAseq"),
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                                 box(numericInput(ns("umap_Height_for_RNAseq"),
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                               ),
                               fluidRow(
                                 plotOutput(ns("umap_for_RNAseq"), height = 600),
                                 downloadButton(ns('downloadUMAP_for_RNAseq'), 'Save')
                               ),
                               shinyBS::bsTooltip(ns("umap_Indicate_for_RNAseq"), "Set the color, shape and facet_wrap of the plot", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("umap_if_square_for_RNAseq"), "Whether x limit is equal to y limit", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("umap_Width_for_RNAseq"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("umap_Height_for_RNAseq"), "Height of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("umap_n_neighbors_for_RNAseq"), "number of nearest neighbors for umap", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "t-SNE plot",
                               fluidRow(
                                 box(selectizeInput(ns("Tsne_Indicate_for_RNAseq"),
                                                    "Color and shape",
                                                    choices = c("condition", "replicate", "Condition", "Replicate"),
                                                    selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                 box(checkboxInput(ns("Tsne_if_square_for_RNAseq"),
                                                   "if square",
                                                   value = FALSE),
                                     width = 6)
                               ),
                               fluidRow(
                                 box(uiOutput(ns("Tsne_perplexity_for_RNAseq")),width = 8),
                                 box(numericInput(ns("Tsne_theseed_for_RNAseq"),
                                                  "set seed for t-SNE",
                                                  value = 42),
                                     width = 4)
                               ),
                               fluidRow(
                                 box(numericInput(ns("Tsne_Width_for_RNAseq"),
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                                 box(numericInput(ns("Tsne_Height_for_RNAseq"),
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                               ),
                               fluidRow(
                                 plotOutput(ns("Tsne_for_RNAseq"), height = 600),
                                 downloadButton(ns('downloadTSNE_for_RNAseq'), 'Save')
                               ),
                               shinyBS::bsTooltip(ns("Tsne_Indicate_for_RNAseq"), "Set the color, shape and facet_wrap of the plot", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Tsne_if_square_for_RNAseq"), "Whether x limit is equal to y limit", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Tsne_Width_for_RNAseq"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Tsne_Height_for_RNAseq"), "Height of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Tsne_perplexity_for_RNAseq"), "perplexity for t-SNE", "top", options = list(container = "body"))
                      ),
                      #tabPanel(title = "Pca plot",
                      #         plotOutput("pca", height = 600),
                      #         downloadButton('downloadPca', 'Save')
                      #)
                      tabPanel(title = "Pearson correlation",
                               fluidRow(
                                 box(selectizeInput(ns("Pearson_pal_for_RNAseq"),
                                                    "color panel",
                                                    choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Blues",  "BuGn", "BuPu", "GnBu", "Greens" , "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
                                                    selected = c("PRGn"), multiple = FALSE),
                                     checkboxInput(ns("Pearson_pal_rev_for_RNAseq"),
                                                   "pal rev",
                                                   value = FALSE), width = 4),
                                 box(numericInput(ns("Pearson_Width_for_RNAseq"),
                                                  "width",
                                                  min = 1, max = 30, value = 7),
                                     numericInput(ns("Pearson_Height_for_RNAseq"),
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 4),
                                 box(numericInput(ns("Pearson_lower_for_RNAseq"),
                                                  "lower",
                                                  min = -1, max = 1, value = -1),
                                     numericInput(ns("Pearson_upper_for_RNAseq"),
                                                  "upper",
                                                  min = -1, max = 1, value = 1), width = 4)
                               ),
                               fluidRow(
                                 box(checkboxInput(ns("add_values_for_RNAseq_person"),
                                                   "Add values",
                                                   value = FALSE),
                                     numericInput(ns("value_size_for_RNAseq_person"),
                                                  "Value size",
                                                  min = 1, max = 30, value = 10),
                                     numericInput(ns("value_digits_for_RNAseq_person"),
                                                  "Value digits",
                                                  min = 1, max = 30, value = 3), width = 12, collapsible = TRUE, collapsed = TRUE)
                               ),
                               fluidRow(
                                 plotOutput(ns("Pearson_correlation_for_RNAseq"), height = 600),
                                 downloadButton(ns('download_Pearson_correlation_for_RNAseq'), 'Save')
                               ),
                               shinyBS::bsTooltip(ns("Pearson_pal_for_RNAseq"), "Set the color panel (from RColorBrewer)", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Pearson_pal_rev_for_RNAseq"), "Whether or not to invert the color palette", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Pearson_Width_for_RNAseq"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Pearson_Height_for_RNAseq"), "Height of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Pearson_lower_for_RNAseq"), "Set the lower limit of the color scale", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Pearson_upper_for_RNAseq"), "Set the upper limit of the color scale", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Gower's distance",
                               fluidRow(
                                 box(selectizeInput(ns("Gower_pal_for_RNAseq"),
                                                    "color panel",
                                                    choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Blues",  "BuGn", "BuPu", "GnBu", "Greens" , "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
                                                    selected = c("YlOrRd"), multiple = FALSE),
                                     checkboxInput(ns("Gower_pal_rev_for_RNAseq"),
                                                   "pal rev",
                                                   value = TRUE), width = 6),
                                 box(numericInput(ns("Gower_Width_for_RNAseq"),
                                                  "width",
                                                  min = 1, max = 30, value = 7),
                                     numericInput(ns("Gower_Height_for_RNAseq"),
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                               ),
                               fluidRow(
                                 box(checkboxInput(ns("add_values_for_RNAseq_gower"),
                                                   "Add values",
                                                   value = FALSE),
                                     numericInput(ns("value_size_for_RNAseq_gower"),
                                                  "Value size",
                                                  min = 1, max = 30, value = 10),
                                     numericInput(ns("value_digits_for_RNAseq_gower"),
                                                  "Value digits",
                                                  min = 1, max = 30, value = 3), width = 12, collapsible = TRUE, collapsed = TRUE)
                               ),
                               fluidRow(
                                 plotOutput(ns("Gowers_distance_for_RNAseq"), height = 600),
                                 downloadButton(ns('download_Gowers_distance_for_RNAseq'), 'Save')
                               ),
                               shinyBS::bsTooltip(ns("Gower_pal_for_RNAseq"), "Set the color panel (from RColorBrewer)", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Gower_pal_rev_for_RNAseq"), "Whether or not to invert the color palette", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Gower_Width_for_RNAseq"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Gower_Height_for_RNAseq"), "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Heatmap",
                               fluidRow(
                               #   column(width = 6,
                               #          box(checkboxInput(ns("manual_heatmap_for_RNAseq"),
                               #                            "Manual heatmap",
                               #                            value = FALSE), width = 12)
                               #   ),
                                 column(width = 6,
                                        box(selectizeInput(ns("colorbar_for_RNAseq"),
                                                           "colorbar",
                                                           choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                                                           selected = c("RdBu"), multiple = FALSE), width = 12)
                                 )
                               ),
                               fluidRow(
                                 box(uiOutput(ns("heatmap_cntrst_for_RNAseq")), width = 12)
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
                                 conditionalPanel(condition = paste0("!input['",ns("if_chooseToshow_for_RNAseq"),"']"),
                                                    # "!input.if_chooseToshow_for_RNAseq",
                                                  box(numericInput(ns("k_for_RNAseq"),
                                                                   "Kmeans",
                                                                   min = 0, max = 15, value = 2),
                                                      width = 4)
                                 ),
                                 # box(numericInput("k_for_RNAseq",
                                 #                  "Kmeans",
                                 #                  min = 0, max = 15, value = 2),
                                 #     width = 4),
                                 box(numericInput(ns("limit_for_RNAseq"),
                                                  "Color",
                                                  min = 0, max = 16, value = 4),
                                     width = 4),
                                 box(numericInput(ns("size_for_RNAseq"),
                                                  "Heat size(4-30)",
                                                  min = 4, max = 30, value = 10),
                                     width = 4)

                               ),
                               fluidRow(
                                 dropdownButton(
                                   checkboxInput(ns("if_chooseToshow_for_RNAseq"),
                                                 "row selected",
                                                 value = FALSE),
                                   uiOutput(ns("chooseToshow_for_RNAseq")),
                                   numericInput(ns("row_font_size_for_RNAseq"),
                                                "row font size",
                                                min = 0, max = 15, value = 6),
                                   numericInput(ns("col_font_size_for_RNAseq"),
                                                "col font size",
                                                min = 0, max = 16, value = 10),
                                   numericInput(ns("Width_for_RNAseq"),
                                                "heatmap width",
                                                min = 1, max = 30, value = 7),
                                   numericInput(ns("Height_for_RNAseq"),
                                                "heatmap height",
                                                min = 1, max = 30, value = 10),
                                   checkboxInput(ns("cluster_columns_for_RNAseq"),
                                                 "Cluster columns",
                                                 value = TRUE),
                                   uiOutput(ns("Custom_columns_order_for_RNAseq")), ## order the cols, may some bugs here
                                   conditionalPanel(condition = paste0("!input['",ns("if_chooseToshow_for_RNAseq"),"']"),
                                                      # "!input.if_chooseToshow_for_RNAseq",
                                                    # checkboxInput(ns("if_mysplit_for_RNAseq"),
                                                    #               "If mysplit",
                                                    #               value = FALSE),
                                                    uiOutput(ns("mysplit_for_RNAseq"))
                                   ),
                                   # checkboxInput("if_mysplit_for_RNAseq",
                                   #               "If mysplit",
                                   #               value = FALSE),
                                   # uiOutput("mysplit_for_RNAseq"),
                                   tooltip = tooltipOptions(title = "More parameters"),
                                   icon = icon("gear"),
                                   # status = "info",
                                   width = "300px"),
                                 br("")
                               ),
                               fluidRow(
                                 uiOutput(ns("plot_heatmap_for_RNAseq")),
                                 downloadButton(ns('downloadHeatmap_for_DERNAseq'), 'Save heatmap')),
                               shinyBS::bsTooltip(ns("colorbar_for_RNAseq"), "Choose the color palette of the colorbar", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("manual_heatmap_for_RNAseq"), "Whether plot the manual heatmap", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("heatmap_cntrst_for_RNAseq"), "Choose the contrasts that you want to show, act when [Manual heatmap] is TRUE, can be one or more contrasts", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("k_for_RNAseq"), "Set the number of k-means clusters", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("limit_for_RNAseq"), "Set the outer limits of the color scale", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("size_for_RNAseq"), "Set the height of the plot region", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("if_chooseToshow_for_RNAseq"), "Whether only show the rows that you want to show", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("chooseToshow_for_RNAseq"), "Choose the rows that you want to show", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("row_font_size_for_RNAseq"), "Set the size of row labels", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("col_font_size_for_RNAseq"), "Set the size of column labels", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Width_for_RNAseq"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Height_for_RNAseq"), "Height of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("cluster_columns_for_RNAseq"), "Whether make cluster on columns", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Custom_columns_order_for_RNAseq"), "Order of column, act when [Cluster columns] is FALSE", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("if_mysplit_for_RNAseq"), "Whether order the splices manually", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("mysplit_for_RNAseq"), "The ordered splices of the heatmap (from top to bottom), only act when both [If mysplit] and [Kmeans] are TRUE", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Custom Volcano",
                               fluidRow(
                                 column(width = 4,
                                        colourpicker::colourInput(inputId = ns("up_color_for_RNAseq"),
                                                    label = "up",
                                                    showColour = "both",
                                                    palette = "square",
                                                    value = "#B2182B")
                                 ),
                                 column(width = 4,
                                        colourpicker::colourInput(inputId = ns("down_color_for_RNAseq"),
                                                    label = "down",
                                                    showColour = "both",
                                                    palette = "square",
                                                    value = "#2166AC")
                                 ),
                                 column(width = 4,
                                        colourpicker::colourInput(inputId = ns("stable_color_for_RNAseq"),
                                                    label = "not significant",
                                                    showColour = "both",
                                                    palette = "square",
                                                    value = "#BEBEBE")
                                 )
                               ),
                               fluidRow(
                                 box(selectizeInput(ns("labelWay_for_RNAseq"),
                                                    "label way",
                                                    choices = c("all significant", "up","down"),
                                                    selected = c("all significant"), multiple = FALSE),width = 4),
                                 box(uiOutput(ns("contrast_for_RNAseq_customvolcano")), width = 8)
                               ),
                               fluidRow(
                                 box(
                                   checkboxInput(ns("if_label_rectangle_for_RNAseq"),
                                                 "label with rectangle",
                                                 value = FALSE),
                                   # numericInput(ns("stroke_for_RNAseq"),
                                   #              "point outside width",
                                   #              min = 0, max = 5, value = 0.6),
                                   checkboxInput(ns("labeled_for_RNAseq_customvolcano"),
                                                 "labeled",
                                                 value = FALSE),
                                   uiOutput(ns("selected_proteins_for_RNAseq")),
                                   numericInput(ns("showNum_for_RNAseq"),
                                                "show number",
                                                min = 0, max = 100000, value = 20),
                                   numericInput(ns("fontSize_for_RNAseq"),
                                                "font size",
                                                min = 0, max = 20, value = 4),
                                   numericInput(ns("dotsize_for_RNAseq"),
                                                "dot size",
                                                min = 0, max = 20, value = 2),
                                   numericInput(ns("custom_volcano_Width_for_RNAseq"),
                                                "width",
                                                min = 1, max = 30, value = 7),
                                   numericInput(ns("custom_volcano_Height_for_RNAseq"),
                                                "height",
                                                min = 1, max = 30, value = 7), width = 8, collapsible = TRUE, collapsed = TRUE),
                                 box(checkboxInput(ns("P_adj_for_RNAseq"),
                                                   "Adjusted p values",
                                                   value = FALSE),
                                     checkboxInput(ns("Same_width_for_RNAseq"),
                                                   "Same width",
                                                   value = FALSE),
                                     width = 4)),
                               fluidRow(
                                 plotOutput(ns("custom_volcano_for_RNAseq"), height = 600),
                                 downloadButton(ns('download_custom_volcano_for_RNAseq'), 'Save custom volcano')
                               ),
                               shinyBS::bsTooltip(ns("up_color_for_RNAseq"), "Set the color of the points those are up-regulated", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("down_color_for_RNAseq"), "Set the color of the points those are down-regulated", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("stable_color_for_RNAseq"), "Set the color of the points those are not significant", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("labelWay_for_RNAseq"), "Choose the way of adding labels, one of all significant, up, down, and selected genes", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("contrast_for_RNAseq_customvolcano"), "Choose the contrast that you want to show", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("if_label_rectangle_for_RNAseq"), "Whether add rectangle underneath the text, making it easier to read", "top", options = list(container = "body")),
                               # shinyBS::bsTooltip(ns("stroke_for_RNAseq"), "Set the thickness of black line around the point", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("labeled_for_RNAseq_customvolcano"), "Whether use the labelway : selected genes", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("selected_proteins_for_RNAseq"), "Choose the point labels to show, act when [label Way] is selected genes", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("showNum_for_RNAseq"), "Set the number of the labels to add, act when [label way] is significant, up or down. Note that, when it is larger than the number of significant result, it is seted to the number of significant result", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("fontSize_for_RNAseq"), "Set the font size of the label", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("dotsize_for_RNAseq"), "Set the dot size of the label", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("custom_volcano_Width_for_RNAseq"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("custom_volcano_Height_for_RNAseq"), "Height of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("P_adj_for_RNAseq"), "Whether or not to use adjusted p values", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Same_width_for_RNAseq"), "Whether the x axis to have the same width from 0", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "MA plot",
                               fluidRow(
                                 box(uiOutput(ns("selected_proteins_for_RNAseq_MAplot")),
                                     uiOutput(ns("contrast_for_RNAseq_MA")), width = 6),
                                 box(checkboxInput(ns("add_rug_for_RNAseq"),
                                                   "Add rug",
                                                   value = FALSE),
                                     checkboxInput(ns("selected_genes_for_RNAseq_MAplot"),
                                                   "Selected genes",
                                                   value = FALSE),
                                     checkboxInput(ns("labeled_for_RNAseq_MAplot"),
                                                   "Labeled",
                                                   value = FALSE),
                                     width = 6)
                               ),
                               fluidRow(
                                 box(numericInput(ns("MAplot_Width_for_RNAsqe"),
                                                  "width",
                                                  min = 1, max = 30, value = 9), width = 6),
                                 box(numericInput(ns("MAplot_Height_for_RNAseq"),
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                               ),
                               fluidRow(
                                 plotOutput(ns("plot_MA_for_RNAseq")),
                                 downloadButton(ns('downloadMAplot_for_RNAseq'), 'Save')),
                               shinyBS::bsTooltip(ns("selected_proteins_for_RNAseq_MAplot"), "Choose the point labels to show, act when [Labeled] is selected", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("contrast_for_RNAseq_MA"), "Choose the contrast that you want to show", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("add_rug_for_RNAseq"), "whether add rug plots in the margins", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("selected_genes_for_RNAseq_MAplot"), "whether show genes, so that you can choose genes and add them on the plotting by setting [Labeled] TRUE", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("labeled_for_RNAseq_MAplot"), "whether add labels on the plotting, act when [Selected genes] is TRUE", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("MAplot_Width_for_RNAsqe"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("MAplot_Height_for_RNAseq"), "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Selected Gene",
                               fluidRow(
                                 box(
                                   uiOutput(ns("selected_genes_for_RNAseq")),
                                   uiOutput(ns("selected_group_for_RNAseq")), width = 8),
                                 box(checkboxInput(ns("ylimZero_for_RNAseq"),
                                                   "ylimZero",
                                                   value = TRUE),
                                     checkboxInput(ns("add_labels_for_RNAseq"),
                                                   "Add labels",
                                                   value = FALSE),
                                     checkboxInput(ns("labels_repel_for_RNAseq"),
                                                   "Labels repel",
                                                   value = TRUE), width = 4)),
                               fluidRow(
                                 box(numericInput(ns("selected_genes_Width_for_RNAsqe"),
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                                 box(numericInput(ns("selected_genes_Height_for_RNAseq"),
                                                  "height",
                                                  min = 1, max = 30, value = 6), width = 6)
                               ),
                               fluidRow(
                                 plotOutput(ns("plot_selected_genes_for_RNAseq"), height = 600),
                                 downloadButton(ns('download_selected_genes_for_RNAseq'), 'Save')
                               ),
                               shinyBS::bsTooltip(ns("selected_genes_for_RNAseq"), "Specify the name of the gene to plot", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("selected_group_for_RNAseq"), "Choose interesting groups", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("ylimZero_for_RNAseq"), "Whether set the lower y-axis limit to 0", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("add_labels_for_RNAseq"), "Whether add labels", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("labels_repel_for_RNAseq"), "Whether use functions of R package ggrepel to place labels", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("selected_genes_Width_for_RNAsqe"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("selected_genes_Height_for_RNAseq"), "Height of the figure to export", "top", options = list(container = "body"))
                      )

               )
        )
      )

    )
  )
}

DEG_server_module <- function(id){
  moduleServer(
    id,
    function(input,output,session){
      ns = session$ns

      annoSpecies_df <- DEP2:::annoSpecies_df()

      ## analysis for DEG-RNAseq
      #### UI functions output ui### --------------------------------------------------------
      output$for_loadcount <- renderUI({
        countData() ## to force countData() will execute immediately after counts file upload
        return(NULL)
      })

      output$dds_design <- renderUI({
        expdesign_1 <- expdesign()
        if(!is.null(expdesign())) {
          expdesign <- expdesign()
          expdesign111 <<- expdesign()
          ddschoices = which( (expdesign %>% apply(., 2, function(x){length(unique(x))})) < nrow(expdesign) )
          shiny::validate(need(
            length(ddschoices) > 0,
            "No column match requirement, check the experiement design"
          ))
          selectizeInput(ns("dds_design"),
                         label = "Select design",
                         choices = c(colnames(expdesign)[ddschoices]), selected = "condition", multiple = TRUE)
        }
      })

      output$choose_expfac <- renderUI({
        if(!is.null(expdesign()) && !is.null(input$dds_design) && input$dds_design != "") {
          selectizeInput(ns("choose_expfac"),
                         label = "Build the contrast upon",
                         choices = c(input$dds_design), selected = "condition", multiple = FALSE)
        }
      })

      output$whether_LRT <- renderUI({
        expdesign11 <<- expdesign()
        shiny::validate(
          need(
            !is.null(expdesign()),
            "Please upload your countmatrix (or log file) first"
          )
        )
        if(!is.null(expdesign())&& !is.null(input$choose_expfac) && input$choose_expfac != ""){
          choose_factor_num <- expdesign() %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character() %>% unique() %>% length()
          if(choose_factor_num > 2)
            return(checkboxInput("perform_LRT_test", "perform a LRT test", value = TRUE))
        }
        # if(!is.null(input$file2) & choose_factor_num > 2) {
        #   # radioButtons("Test_for_RNAseq",
        #   #                       "Test",
        #   #                       choices = c("Wald Test", "LRT Test"),
        #   #                       selected = "LRT Test")
        #   checkboxInput("perform_LRT_test", "perform a LRT test", value = TRUE)
        # }
      })

      output$reduced_model <- renderUI({
        shiny::validate(
          need(
            !is.null(expdesign()),
            ""
          )
        )
        if(!is.null(expdesign())&& !is.null(input$choose_expfac) && input$choose_expfac != ""){
          choose_factor_num <- expdesign() %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character() %>% unique() %>% length()
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
        if(!is.null(expdesign()) & choose_factor_num > 2 & input$perform_LRT_test & input$dds_design != "") {
          selectizeInput(ns("reduced_model"),
                         label = "Select the reduced model",
                         choices = c("", "1", input$dds_design), selected = "", multiple = TRUE)
        }
      })

      output$control <- renderUI({
        shiny::validate(
          need(
            !is.null(expdesign()),
            ""#Please upload your countmatrix and expdesign files first
          )
        )
        shiny::validate(
          need(
            input$choose_expfac != "",
            "Please select the design and an experimental factor to build the contrast upon"
          )
        )
        shiny::validate(
          need(
            !is.null(input$choose_expfac),
            "Please select the design and an experimental factor to build the contrast upon"
          )
        )

        if(!is.null(expdesign()) && input$contrasts == "control" && !is.null(input$choose_expfac) &&  input$choose_expfac != "") {
          expdesign = try(expdesign(),silent = T)
          # conditionalPanel(condition = paste0("input['",ns("choose_expfac"),"']"," != ''"),
          # if( (!is.null(upload_log())) & input$uploadmode == "fromLog"){
          #   selected = upload_log()$inputVals()$control
          #   selected11 <<- selected
          #   choices11 <<- expdesign %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character()
          # }else{
          #   selected = NULL
          # }

          selectizeInput(ns("control"),
                         label = "Control",
                         choices = expdesign %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character() %>% unique,
                         selected = NULL, multiple = FALSE)
          # )
        }
      })

      output$test_manual_for_RNAseq <- renderUI({
        if(!is.null(expdesign())  & input$contrasts == "manual") {
          test_manual_name <- expdesign() %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character() %>% unique()
          test_manual_name <- cbind(combn(test_manual_name,2),combn(test_manual_name,2, FUN = rev))
          selectizeInput(ns("test_manual_for_RNAseq"),
                         label = "Manual test",
                         choices = apply(test_manual_name, 2, function(i){paste(i[1], i[2], sep = "_vs_")}), selected = NULL, multiple = TRUE)
        }
      })

      output$ui_nrcores <- renderUI({
        mincores <- 1
        maxcores <- BiocParallel::multicoreWorkers()
        sliderInput(ns("nrcores"),
                    label = "Cores to use",
                    min = mincores, max = maxcores, value = 1, step = 1
        )
      })

      output$ui_idtype <- renderUI({
        if ((!is.null(input$transid_for_RNAseq)) && input$transid_for_RNAseq) {
          shiny::validate(
            need(
              input$species_select != "",
              "Please select the species"
            )
          )
          # std_choices <- c("ENSEMBL", "ENTREZID", "REFSEQ", "SYMBOL")
          std_choices <- c("ENSEMBL", "ENTREZID", "REFSEQ")
          annopkg <- annoSpecies_df$pkg[annoSpecies_df$species == input$species_select]
          shiny::validate(
            need(
              require(annopkg, character.only = TRUE),
              paste0("The package ", annopkg, " is not installed/available. Try installing it with BiocManager::install('", annopkg, "')")
            )
          )
          require(annopkg, character.only = TRUE)
          pkg_choices <- keytypes(get(annopkg))
          std_choices <- intersect(std_choices, pkg_choices)

          tagList(
            selectInput(ns("idtype"), "select your id type", choices = std_choices)
            # actionButton(ns("transformeID_butt"), "Transform ID")
          )
        } else {
          return(NULL)
        }
      })

      output$umap_n_neighbors_for_RNAseq <- renderUI(
        sliderInput(inputId = ns("umap_n_neighbors_for_RNAseq"),"number of nearest neighbors",min=2,max=min(15,ncol(countData())),value=4)
      )

      output$Tsne_perplexity_for_RNAseq <- renderUI(
        sliderInput(inputId = ns("Tsne_perplexity_for_RNAseq"),"perplexity",min=1,max= (ncol(countData())-1)/3 - 2^(-31) ,value=2)
      )


      # observeEvents upon uploadmode and log_file. fromLog -----------------------------------------------------------------

      ## updata options according uploadmode and upload_log
      Columns_to_updata = reactiveVal(value = 0)
      Columns_rendered = reactiveVal(value = 0)

      observeEvent(upload_log(),{

        cat("Update options 1 \n")
        Option_to_updata <- c("transid_for_RNAseq",
                              "independent_filtering", "Shrink_lfc","Use_IHW") # CheckboxInput
        for(i in Option_to_updata){
          updateCheckboxInput(session = session, inputId = i,
                              value = upload_log()$inputVals()[[i]] )
        }

        Option_to_updata2 <- c("filter_rowsum","filter_missnum") # NumericInput
        for(i in Option_to_updata2){
          updateNumericInput(session = session, inputId = i,
                             value = upload_log()$inputVals()[[i]] )
        }

        Option_to_updata3 <- c("contrasts")
        for(i in Option_to_updata3){
          updateRadioButtons(session = session, inputId = i,
                             selected = upload_log()$inputVals()[[i]] )
        }

        Option_to_updata4 <- c("species_select")
        for(i in Option_to_updata4){
          updateSelectInput(session = session, inputId = i,
                               selected = upload_log()$inputVals()[[i]] )
        }

        Columns_to_updata = Columns_to_updata(1)
        message("Update options 1 finished \n")
      },ignoreNULL = T,ignoreInit = T)

      observeEvent(
        input$dds_design,
        {
          Columns_rendered = Columns_rendered(1)
        },ignoreNULL = T,ignoreInit = T
      )

      observeEvent(
        c(Columns_to_updata(),Columns_rendered()),
        {
          if(Columns_to_updata() == 1 & Columns_rendered() == 1){
            cat("Update options in 'RNAseq_settings' \n")

            column_opt = c("dds_design", "choose_expfac","control")
            for(i in column_opt){
              # upload_log22 <<- upload_log()
              updateSelectizeInput(session = session, inputId = i,
                                   selected = upload_log()$inputVals()[[i]] )

            }

            # cat("Update control \n")
            # Sys.sleep(3)
            # updateSelectizeInput(session = session, inputId = "control",
            #                      selected = upload_log()$inputVals()[["control"]] )

            updateCheckboxInput(session = session, inputId = "perform_LRT_test",
                                value = upload_log()$inputVals()[["perform_LRT_test"]] )

            Threshold_opt = c("lfc","p","curvature","x0_fold")
            for(i in column_opt){
              # message("update ",i," to '", upload_log()$inputVals()[[i]],"'\n", sep = "")
              updateNumericInput(session = session, inputId = i,
                                 value = upload_log()$inputVals()[[i]] )
            }

            Columns_rendered = Columns_rendered(0) # back to 0 once restore the options
            Columns_to_updata = Columns_to_updata(0) # back to 0 once restore the options

            message("Update 'RNAseq_settings' finnished \n")
          }

        },
        ignoreNULL = T,ignoreInit = T
      )

      ## updata anno option according uploadmode
      observe({
        # uploadmode11 <<- input$uploadmode
        # cat("uploadmode shift")
        if(input$uploadmode == "fromTable"){
          updateRadioButtons(session = session,
                             inputId = "anno",
                             label = "Sample annotation",
                             choices = list("Parse from columns" = "columns",
                                            "Use Experimental Design" = "expdesign"),
                             selected = "columns"
          )
        }else if(input$uploadmode == "fromLog"){
          updateRadioButtons(session = session,
                             inputId = "anno",
                             label = "Sample annotation",
                             choices = list("Parse from columns" = "columns",
                                            "Use ExpDesign in log" = "logexpdesign"),
                             selected = "logexpdesign"
          )
        }

      })


      ### check input ### --------------------------------------------------------
      iv <- InputValidator$new()
      iv$add_rule("p", sv_between(0.0001, 0.1))
      iv$add_rule("lfc", sv_between(0, 10))
      iv$enable()

      #### main reactive functions ----

      # transformit <- reactiveVal(FALSE)

      upload_log <- reactive({
        logFile <- input$file_log

        if( !is.null(logFile) ){

          cat("Read infile_log\n")
          upload_log <- readRDS(logFile$datapath)

          upload_log111 <<- upload_log
          if( upload_log$appVersion == DEP2:::app_version && inherits(upload_log,"transcriptome_log") ){
            sendSweetAlert(
              session = shiny::getDefaultReactiveDomain(),
              title = "Load log successfully!",
              text = "Pass version check and load.",
              type = "success"
            )
            message("Read successfully\n")
            return(upload_log)

          }else if(!inherits(upload_log,"transcriptome_log")){
            sendSweetAlert(
              session = shiny::getDefaultReactiveDomain(),
              title = "Error !",
              text = "The uploaded file is not a log for transcriptome analysis, unload the log.\n
              Please ensure that: The log file is a log from RNA-seq analysis. And the app version is correct,
              DEP2 app is only compatible the log from same version. Old log may not be used on a updated app.
              ",
              type = "warning"
            )
            return(NULL)

          }else if(upload_log$appVersion != DEP2:::app_version){
            sendSweetAlert(
              session = shiny::getDefaultReactiveDomain(),
              title = "Error !",
              text = "The version of uploaded file don't match current app version now.\n
              Please ensure that: The log file is a log from RNA-seq analysis. And the app version is correct,
              DEP2 app is only compatible the log from same version. Old log may not be used on a updated app.
              ",
              type = "warning"
            )
            return(NULL)
          }

        }else{
          # cat("log is null \n")
          return(NULL)
        }

      })

      countData <- reactive({

        # from log
        upload_log = upload_log()
        if(!is.null(upload_log) & input$uploadmode == "fromLog"){
          # cat("a1")
          my_countData <<- upload_log$resultVals()$countData
          # reset("transid_opts")
          return(my_countData)
        }else if( is.null(upload_log) & input$uploadmode == "fromLog"){
          # cat("b1")
          reset("transid_opts")
          return(NULL)
        }

        # from table
        inFile <- input$file1
        # cat("c1")
        if (is.null(inFile))
          return(NULL)
        countData <- read.csv(inFile$datapath, sep = "\t", stringsAsFactors = F, header = T)
        row_name = countData[,1]
        countData = countData[ , -1]
        countData = as.matrix(countData)
        rownames(countData) = row_name

        # transformit(FALSE)
        reset("transid_opts") ## reset the ID transformation options after counts file change
        return(countData)
      })

      expdesign <- reactive({
        # from log
        if(input$uploadmode == "fromLog" && !is.null(countData()) ){
          # cat("aa \n")
          upload_log = upload_log()

          if(input$anno == "logexpdesign"){
            expdesign <- upload_log$resultVals()$expdesign

          }else if(input$anno == "columns"){

            countData <- countData()
            label = colnames(countData)
            expdesign <- get_exdesign_parse(label, mode = "delim", sep = "_",
                                            remove_prefix = F, remove_suffix = F)

          }
        }else if(input$uploadmode == "fromLog" && is.null(countData())){
          # cat("bb")
          return(NULL)
        }

        # from table
        if(input$uploadmode == "fromTable"){
          # cat("cc")
          inFile <- input$file2
          if ( input$anno == "columns" & (!is.null(countData())) ){

            countData <- countData()
            label = colnames(countData)
            expdesign <- get_exdesign_parse(label, mode = "delim", sep = "_",
                                            remove_prefix = F, remove_suffix = F)

          }else if(input$anno == "expdesign" & (!is.null(inFile)) & (!is.null(countData())) ){
            expdesign <- read.csv(inFile$datapath, header = TRUE,
                                  sep = "\t", stringsAsFactors = FALSE) %>%
              mutate(id = row_number())
          }else{
            return(NULL)
          }
        }

        # if(input$uploadmode == "fromTable"){
        #   inFile2 <- input$file2
        #   inFile1 <- input$file1
        #   countData <- countData()
        #   if(is.null(inFile2)){
        #     if(is.null(countData())){
        #       return(NULL)
        #     }else{
        #       countData <- countData()
        #       label_names = colnames(countData)
        #       expdesign = DEP2::get_exdesign_parse(label_names)
        #     }
        #   }else{
        #     expdesign = read.csv(inFile$datapath, sep = "\t", header = T)
        #   }
        # }

        my_expdesign <<- expdesign
        return(expdesign)




      })




      returnval <- reactiveVal()
      #### analyze event ----
      observeEvent(input$analyze_for_DERNAseq, {
        ## check inputs
        if(is.null(input$file1)) {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = "Please upload Countmatrix first.",
            type = "warning"
          )
        }
        req(!is.null(input$file1))

        if(is.null(input$dds_design) || input$dds_design == "") {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = "Please ensure that: select design columns in the RNAseq_settings opt",
            type = "warning"
          )
        }
        req(!is.null(input$dds_design))

        if(is.null(input$choose_expfac)) {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = "Please ensure that: select design columns in the RNAseq_settings opt",
            type = "warning"
          )
        }
        req(!is.null(input$choose_expfac))

        if(input$contrasts == "control" && (is.null(input$control) || input$control == "")) {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = "Please ensure that: Control can not be empty when Contrasts is control",
            type = "warning"
          )
          req(!is.null(input$control) && input$control != "")
        }

        output$downloadTable_for_DERNAseq <- renderUI({
          selectizeInput(ns("dataset_for_DERNAseq"),
                         "Choose a dataset to save" ,
                         c("full_results","significant_results"
                         ))
        })

        output$downloadButton_for_DERNAseq <- renderUI({
          downloadButton(ns('download_DERNAseq'), 'Save table', class = "download_DERNAseq")
        })

        output$downloadButton_for_log <- renderUI({ ## download log file
          downloadButton(session$ns("Save_log"), "Export log file")
        })

        output$pca_top_n_for_RNAseq  <- renderUI(
          sliderInput(inputId = ns("pca_top_n_for_RNAseq"),"number of top variable genes to consider",min = 2, max = nrow(deg()@rlg), value=500)
        )


        output$significantBox_for_DERNAseq <- renderInfoBox({
          num_total_RNAseq <- as.data.frame(diff()@test_result) %>%
            nrow()
          num_signif_RNAseq <- as.data.frame(deg_sig()@test_result) %>%
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
            plotOutput(ns("heatmap_for_DERNAseq"), height = (100 * as.numeric(input$size_for_RNAseq)))
          })
        })


        output$heatmap_cntrst_for_RNAseq <- renderUI({
          df <- deg_df()
          cols <- grep("_significant$",colnames(df))
          selectizeInput(ns("heatmap_cntrst_for_RNAseq"),
                         "Contrast",
                         choices = gsub("_significant", "", colnames(df)[cols]),
                         selected = NULL,
                         multiple = TRUE)
        })


        output$mysplit_for_RNAseq <- renderUI({
          selectizeInput(ns("mysplit_for_RNAseq"),
                         "my split",
                         choices = c(1 : input$k_for_RNAseq), multiple = TRUE)

        })

        output$Custom_columns_order_for_RNAseq <- renderUI({
          selectizeInput(ns("Custom_columns_order_for_RNAseq"),
                         "Custom columns order",
                         choices = colnames(deg()), multiple = TRUE)

        })

        output$selected_proteins_for_RNAseq <- renderUI({
          withProgress(message = 'Please Wait',value = 0.66, {
            if(input$labeled_for_RNAseq_customvolcano){
              res_for_RNAseq = diff()@test_result %>% as.data.frame(.)
              selectizeInput(ns("selected_proteins_for_RNAseq"),
                             "selected genes",
                             choices = if("symbol" %in% colnames(res_for_RNAseq)) {c(NULL,res_for_RNAseq$symbol)} else {c(NULL, rownames(res_for_RNAseq))}, selected = NULL, multiple = TRUE)
            } else {
              return(NULL)
            }
          })
        })

        output$contrast_for_RNAseq_customvolcano <- renderUI({
          df <- deg_df()
          cols <- grep("_significant$",colnames(df))
          selectizeInput(ns("contrast_for_RNAseq_customvolcano"),
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
            res_for_RNAseq = diff()@test_result %>% as.data.frame(.)
            selectizeInput(ns("selected_proteins_for_RNAseq_MAplot"),
                           "selected genes",
                           choices = if("symbol" %in% colnames(res_for_RNAseq)) {c(NULL,res_for_RNAseq$symbol)} else {c(NULL, rownames(res_for_RNAseq))}, selected = NULL, multiple = TRUE)

          }
        })

        output$contrast_for_RNAseq_MA <- renderUI({
          df <- deg_df()
          cols <- grep("_significant$",colnames(df))
          selectizeInput(ns("contrast_for_RNAseq_MA"),
                         "Contrast",
                         choices = gsub("_significant", "", colnames(df)[cols]), selected = gsub("_significant", "", colnames(df)[cols])[1], multiple = FALSE)
        })

        output$selected_genes_for_RNAseq <- renderUI({
          withProgress(message = 'Please Wait',value = 0.66, {
            res_for_RNAseq = diff()@test_result %>% as.data.frame(.)
            selectizeInput(ns("selected_genes_for_RNAseq"),
                           "selected genes",
                           choices = if("symbol" %in% colnames(res_for_RNAseq)) {res_for_RNAseq$symbol} else {rownames(res_for_RNAseq)}, selected = if("symbol" %in% colnames(res_for_RNAseq)) {res_for_RNAseq$symbol[1]} else {rownames(res_for_RNAseq)[1]}, multiple = FALSE)

          })
        })

        output$selected_group_for_RNAseq <- renderUI({
          selectizeInput(ns("selected_group_for_RNAseq"),
                         "Color by",
                         choices = input$dds_design, selected = input$choose_expfac, multiple = FALSE)

        })


        dds_1 <- reactive({
          dds_1 <- DESeqDataSetFromMatrix(countData = countData(),colData = expdesign(), design = as.formula(paste0("~", paste(input$dds_design, collapse = " + "))))
        })

        dds_filter <- reactive({
          dds_1 = dds_1()
          if(input$filter_missnum > ncol(dds_1)){
            filter_missnum =ncol(dds_1)
            warning("input filter missing number is to large, and reset to ",filter_missnum)
          }
          dds_filter <- filter_se(dds_1, rowsum_threshold = input$filter_rowsum, missnum = input$filter_missnum)
        })

        dds <- reactive({
          withProgress(message = 'Running DEseq', value = 0.66, {
            choose_factor_num <- expdesign() %>% dplyr::select(input$choose_expfac) %>% unlist() %>% as.character() %>% unique() %>% length()
            if (input$nrcores == 1) {
              if(choose_factor_num > 2) {
                if(input$perform_LRT_test) {
                  DESeq(dds_filter(), test = "LRT", full = as.formula(paste0("~", paste(input$dds_design, collapse = " + "))), reduced = as.formula(paste0("~", paste(input$reduced_model, collapse = " + "))))
                } else {
                  DESeq(dds_filter())
                }
              } else {
                DESeq(dds_filter())
              }
            } else {
              if(choose_factor_num > 2) {
                if(input$perform_LRT_test) {
                  DESeq(dds_filter(), test = "LRT", full = as.formula(paste0("~", paste(input$dds_design, collapse = " + "))), reduced = as.formula(paste0("~", paste(input$reduced_model, collapse = " + "))), parallel = TRUE, BPPARAM = MulticoreParam(workers = input$nrcores))
                } else {
                  DESeq(dds_filter(), parallel = TRUE, BPPARAM = MulticoreParam(workers = input$nrcores))
                }
              } else {
                DESeq(dds_filter(), parallel = TRUE, BPPARAM = MulticoreParam(workers = input$nrcores))
              }
            }

          })
        })

        transformit <- reactive({
          cat("check transform options...  ")
          transformit <- (input$transid_for_RNAseq && !is.null(input$species_select) && input$species_select != "" && !is.null(input$idtype))
          if(transformit) cat("transform ID\n")
          return(transformit)
        })

        ### differential test in DEseq2, generate a DEGdata
        diff <-  reactive({
          withProgress(message = 'Get results', value = 0.33, {
            diff <- DEP2:::get_res(dds = dds(),
                                   type = input$contrasts,
                                   control = input$control,
                                   test_contrasts = input$test_manual_for_RNAseq,
                                   contrast_upon = input$choose_expfac,
                                   filter_ihw = input$Use_IHW,
                                   independentFiltering = input$independent_filtering,
                                   lfcshark = input$Shrink_lfc)
            message("Get results finished.")

            diff <- DEP2::rlg_deg(diff)
            diff <- DEP2::ntf_deg(diff)
            message("rlg & ntf finished.")

            isolate({transformit_it <- transformit()})

            input$analyze_for_DERNAseq
            if(transformit_it){
              # incProgress(2/3,message = "Transforming ID")
              message("Transforming ID")
              withProgress(message = "Transforming ID, please wait!",value = 0.5,{
                diff <- DEP2::ID_transform(diff, from_columns = "rownames", fromtype = input$idtype, species = "Human", replace_rowname = "SYMBOL")
              })
              # cat("finished \n")
            }

            message("diff finished.")
            return(diff)

          })
        })

        ### add_rejection to exctract significant genes
        deg <- reactive({
          deg <- DEP2::add_rejections(diff = diff(), alpha = input$p, lfc = input$lfc)
          returnval <- returnval(deg)
          cat("deg")

          return(deg)
        })
        deg_df <- reactive(deg()@test_result)

        deg_sig <- reactive({
          deg <- deg()
          test_res <- deg@test_result
          deg_sig <- deg[test_res$significant, ]
          return(deg_sig)
        })


        output$table_for_DERNAseq <- DT::renderDataTable(
          as.data.frame(deg_sig()@test_result) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname),
          filter = 'top',
          options = list( autoWidth = F,scrollX = TRUE )
        )

        # res merge with ntd, only for download
        res_for_RNAseq_merge <- reactive({
          res_for_RNAseq_merge = as.data.frame(deg_df())
          res_for_RNAseq_merge$gene = rownames(res_for_RNAseq_merge)
          ntd_res <- merge(res_for_RNAseq_merge, as.data.frame(deg()@ntf), by.x = "gene", by.y = "row.names")
          return(ntd_res)
        })

        ### Download objects and functions ### ------------------------------------
        datasetInput_for_RNAseq <- reactive({
          switch(input$dataset_for_DERNAseq,
                 "full_results" = res_for_RNAseq_merge() %>% dplyr::arrange(desc(significant)),
                 "significant_results" = as.data.frame(deg_sig()@test_result) %>% tibble::rownames_to_column(.) %>% dplyr::rename(., Gene = rowname))
        })

        output$download_DERNAseq <- downloadHandler(
          filename = function() { paste(input$dataset_for_DERNAseq, ".csv", sep = "") },
          content = function(file) {
            write.csv(datasetInput_for_RNAseq(),
                      file, row.names = FALSE
            ) }
        )

        output$Save_log <- downloadHandler(
          filename = function() { paste("Export_log_rna", format(Sys.time(), "%m%d_%H%M"),
                                        ".logrds", sep = "") },
          content = function(file) {
            withProgress(message = 'Compressing result. Please wait ...', value = 0.66, {
              thelog = transcriptome_log$new(input = input,
                                             reactive_vals = list(
                                               countData = countData(),
                                               expdesign = expdesign(),
                                               returnval = returnval(),
                                               dds_1 = dds_1(),
                                               dds = dds(),
                                               diff = diff(),
                                               deg = deg(),
                                               deg_df = deg_df(),
                                               deg_sig = deg_sig()
                                             ),
                                             app_version = DEP2:::app_version
              )
              thelog11 <<- thelog
              saveRDS(thelog, file = file)
            })
          }
        )

        data_forRNAseq_heatmap <- reactive({
          if(input$transid_for_RNAseq) {
            data_1 <- deg()@ntf
            rownames(data_1) <- diff()@test_result$symbol
          } else {
            data_1 <- deg()@ntf
          }
          return(data_1)
        })

        #### plots ----
        heatmap_input_for_RNAseq <- reactive({
          # withProgress(message = 'Plotting', value = 0.66, {
            thedeg <- deg()
            if(input$cluster_columns_for_RNAseq || is.null(input$Custom_columns_order_for_RNAseq) || length(input$Custom_columns_order_for_RNAseq) < 2) {
              thedeg <- thedeg
            } else {
              cols <- input$Custom_columns_order_for_RNAseq
              thedeg <- thedeg[,cols]
              # DEP2::Order_cols(deg(), input$Custom_columns_order_for_RNAseq)
            }

            DEP2::plot_heatmap(object = thedeg,
                               type = input$pres_for_RNAseq,
                               manual_contrast = input$heatmap_cntrst_for_RNAseq,
                               kmeans = ifelse(input$k_for_RNAseq > 1, T , F),
                               k = input$k_for_RNAseq,
                               color = input$colorbar_for_RNAseq,
                               col_limit = input$limit_for_RNAseq,
                               row_font_size = input$row_font_size_for_RNAseq,
                               col_font_size = input$col_font_size_for_RNAseq,
                               cluster_columns = input$cluster_columns_for_RNAseq,
                               split_order = input$mysplit_for_RNAseq,
                               chooseToshow = input$chooseToshow_for_RNAseq
                               )

          # })
        })

        name_for_ht_choose <- reactive({
          ht_data <- DEP2::get_signicant(deg(),
                                         # contrasts = input$manual_name
                                         contrasts = input$heatmap_cntrst)
          heatmap_rownames = rownames(ht_data)
          return(heatmap_rownames)
        })

        output$chooseToshow_for_RNAseq <- renderUI({
          selectizeInput(ns("chooseToshow_for_RNAseq"),
                         "selected genes",
                         choices = name_for_ht_choose(), selected = NULL, multiple = TRUE)
        })

        custom_volcano_input_for_RNAseq <- reactive({
          withProgress(message = 'Plotting', value = 0.66, {

            DEP2::plot_volcano(object = deg(),
                               contrast = input$contrast_for_RNAseq_customvolcano,
                               adjusted = input$P_adj_for_RNAseq,
                               label_trend = switch (input$labelWay_for_RNAseq,
                                                     "all significant" = "all",
                                                     "up" = "up",
                                                     "down" = "down",
                                                     "selected genes" = "none"),
                               label_number = input$showNum_for_RNAseq,
                               chooseTolabel = input$selected_proteins_for_RNAseq,
                               label_size = input$fontSize_for_RNAseq,
                               dot_size = input$dotsize_for_RNAseq,
                               x_symmetry = input$Same_width_for_RNAseq,
                               fcCutoff = input$lfc,
                               pCutoff = input$p,
                               down_color = input$down_color_for_RNAseq,
                               stable_color = input$stable_color_for_RNAseq,
                               up_color = input$up_color_for_RNAseq
                               )
          })
        })

        pca_input_for_RNAseq <- reactive({
          # plot_my_pca(data = deg(), coldata = expdesign(), indicate = input$Indicate_for_RNAseq, if_square = input$if_square_for_RNAseq, n = input$pca_top_n_for_RNAseq)
          DEP2::plot_pca(object = deg(), features = "genes", indicate = input$Indicate_for_RNAseq, if_square = input$if_square_for_RNAseq, n = input$pca_top_n_for_RNAseq)
        })

        umap_input_for_RNAseq <- reactive({
          # plot_umap_for_rnaseq(data = deg()@rlg, coldata = expdesign(), indicate = input$umap_Indicate_for_RNAseq, if_square = input$umap_if_square_for_RNAseq, n_neighbors = input$umap_n_neighbors_for_RNAseq)
          DEP2::plot_umap(object = deg(), features = "genes", indicate = input$umap_Indicate_for_RNAseq, if_square = input$umap_if_square_for_RNAseq, n_neighbors = input$umap_n_neighbors_for_RNAseq)
        })

        Tsne_input_for_RNAseq <- reactive({
          DEP2::plot_Tsne(object = deg(), features = "genes", indicate = input$Tsne_Indicate_for_RNAseq, if_square = input$Tsne_if_square_for_RNAseq,
                          perplexity = input$Tsne_perplexity_for_RNAseq, theseed = input$Tsne_theseed_for_RNAseq)
          # plot_Tsne_for_rnaseq(data = deg()@rlg, coldata = expdesign(), indicate = input$Tsne_Indicate_for_RNAseq, if_square = input$Tsne_if_square_for_RNAseq, perplexity = input$Tsne_perplexity_for_RNAseq, theseed = input$Tsne_theseed_for_RNAseq)
        })

        Pearson_correlation_input_for_RNAseq <- reactive({
          DEP2::plot_cor(deg(), pal = input$Pearson_pal_for_RNAseq, pal_rev = input$Pearson_pal_rev_for_RNAseq,
                         lower = input$Pearson_lower_for_RNAseq, upper = input$Pearson_upper_for_RNAseq,
                         add_values = input$add_values_for_RNAseq_person, value_size = input$value_size_for_RNAseq_person,
                         digits = input$value_digits_for_RNAseq_person)
        })

        Gowers_distance_input_for_RNAseq <- reactive({
        DEP2::plot_dist(x = deg(), pal = input$Gower_pal_for_RNAseq, pal_rev = input$Gower_pal_rev_for_RNAseq,
                        add_values = input$add_values_for_RNAseq_gower,
                        value_size = input$value_size_for_RNAseq_gower, digits = input$value_digits_for_RNAseq_gower)
        })

        diagno_dispests_input_for_RNAseq <- reactive({
          DESeq2::plotDispEsts(dds())
        })

        MA_plot_input_for_RNAseq <- reactive({
          # Create a Progress object
          progress <- shiny::Progress$new()
          progress$set(message = "Plotting", value = 0.66)
          # Close the progress when this reactive exits (even if there's an error)
          on.exit(progress$close())

          DEP2::plot_ma_RNA(deg(), contrast = input$contrast_for_RNAseq_MA, hlines = input$lfc,
                            intgenes = if(input$labeled_for_RNAseq_MAplot & input$selected_genes_for_RNAseq_MAplot) {input$selected_proteins_for_RNAseq_MAplot} else {NULL},
                            add_rug = input$add_rug_for_RNAseq
          )
        })

        selected_genes_input_for_RNAseq <- reactive({
          # Create a Progress object
          # progress <- shiny::Progress$new()
          # progress$set(message = "Plotting", value = 0.66)
          # # Close the progress when this reactive exits (even if there's an error)
          # on.exit(progress$close())
          if("symbol" %in% colnames(as.data.frame(diff()@test_result))) {
            Annotation_obj <- as.data.frame(diff()@test_result)
            Gene <- rownames(as.data.frame(diff()@test_result))[match(input$selected_genes_for_RNAseq, as.data.frame(diff()@test_result)$symbol)]
          } else {
            Annotation_obj <- NULL
            Gene <- input$selected_genes_for_RNAseq
          }
          print(ggplot_Counts(dds = dds(), gene = Gene, annotation_obj = Annotation_obj, intgroup = input$selected_group_for_RNAseq, ylimZero = input$ylimZero_for_RNAseq, add_labels = input$add_labels_for_RNAseq, labels_repel = input$labels_repel_for_RNAseq))
        })


        output$heatmap_for_DERNAseq <- renderPlot({
          heatmap_input_for_RNAseq()
        })

        output$custom_volcano_for_RNAseq <- renderPlot({
          print(custom_volcano_input_for_RNAseq())
        })

        output$pca_for_RNAseq <- renderPlot({
          pca_input_for_RNAseq()
        })


        output$umap_for_RNAseq <- renderPlot({
          umap_input_for_RNAseq()
        })

        output$Tsne_for_RNAseq <- renderPlot({
          Tsne_input_for_RNAseq()
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

        output$downloadUMAP_for_RNAseq <- downloadHandler(
          filename = 'UMAP.pdf',
          content = function(file) {
            pdf(file, width = input$umap_Width_for_RNAseq, height = input$umap_Height_for_RNAseq)
            print(umap_input_for_RNAseq())
            dev.off()
          }
        )

        output$downloadTSNE_for_RNAseq <- downloadHandler(
          filename = 'TSNE.pdf',
          content = function(file) {
            pdf(file, width = input$Tsne_Width_for_RNAseq, height = input$Tsne_Height_for_RNAseq)
            print(Tsne_input_for_RNAseq())
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
            DESeq2::plotDispEsts(dds())
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

        # return(res_for_RNAseq_1)
      })

      ###"question for RNA"
      #### help ----
      observeEvent(input$help_format_RNA, {
        showModal(modalDialog(
          title = "Format specifications for DEG-RNAseq",
          includeMarkdown(system.file("extdata", "DEG_RNAseq.md", package = "DEP2")),
          h4("Example:"),
          tags$img(
            # src = base64enc::dataURI(file = system.file("extdata", "DEG_RNAseq.png", package = "DEP2"), mime = "image/png"),
            src = base64enc::dataURI(file = "www/DEG_RNAseq.png", mime = "image/png"),
            width = 750
          ),
          easyClose = TRUE,
          footer = NULL,
          size = "l"
        ))
      })

      return(returnval)

    }
  )
}
