library(QFeatures)

DEP_pep_UI <- function(id, labelname){
  tagList(
    sidebarLayout(
      # sidebarPanel = sidebarPanel(),
      DEP_pep_sidebar_mod(id = id,labelname=labelname), ## DEP sidebar module
      DEP_pep_body_mod(id = id)  ## DEP pagebody module
      # mainPanel = mainPanel()
    ))
}


DEP_pep_sidebar_mod <-  function(id, labelname){
  ns = NS(id)
  cat(paste0("DEP_pep_sidebar ns is",ns(""),"\n"))

  tagList(
    # sidebarLayout(
    sidebarPanel(
      width = 3,
      h3(labelname),

      h4("Aggregation option"),
      bsCollapse( id = "pep_sidebar",
                  open = "Peptides Files",
                  multiple = F,
                  bsCollapsePanel("Peptides Files", selected = TRUE,
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

                                  # file input options, from Quantity table
                                  conditionalPanel(
                                    condition = paste0("input['",ns("uploadmode"),"']","== 'fromTable'"), ## condition in mod
                                    fluidRow(
                                      column(width = 9,
                                             fileInput(ns('file1'),width = "300px",
                                                       'Peptides.txt',
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
                                                                  '.csv','.gz','.tsv'))
                                             )
                                             # fileInput(ns('resultRData'),'Load saved aggregated result RData',width = "300px")
                                      ),
                                      column(
                                        width = 1, h4(),
                                        br(),
                                        actionButton(ns("help_format_pep"),
                                                     label = "", icon = icon("question-circle"),
                                                     style = "color: #f6f6f6; background-color: #2c3b41; border-color: #2c3b41"
                                        ),
                                        shinyBS::bsTooltip(
                                          ns("help_format_pep"),
                                          "How to provide your input data",
                                          "top",
                                          options = list(container = "body")
                                        )
                                      )
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

                                  radioButtons(ns("anno"),
                                               "Sample annotation",
                                               choices = list("Parse from columns" = "columns",
                                                              "Use Experimental Design" = "expdesign"),
                                               selected = "columns"),
                                  shinyBS::bsTooltip(ns("anno"), "Sample annotation type", "top", options = list(container = "body"))
                  ),
                  bsCollapsePanel("Data proccess",
                                  style = "primary",
                                  uiOutput(ns("name")),
                                  uiOutput(ns("id")),
                                  selectizeInput(ns("delim"), "Delimiter", choices = c(";", "|", ","), selected = ";"),
                                  uiOutput(ns("filt")),
                                  numericInput(ns("thr"), "Allowed max.miss.num at least one condition", min = 0, max = 20, value = 0),

                                  uiOutput(ns("intensitycols")),

                                  prettySwitch(
                                    inputId = ns("imp_norm_order"),
                                    label = "Imputation before normalization",
                                    value = T,
                                    status = "success",
                                    fill = T
                                  ),

                                  selectizeInput(ns("norm_method"), "Peptide normalization method",
                                                 choices = c("quantiles","quantiles.robust", "diff.median",
                                                             "vsn"
                                                             # , "center.mean","center.median",
                                                             # "div.mean","div.median", "diff.median"
                                                 ),
                                                 selected = "quantiles.robust", multiple = FALSE),
                                  radioButtons(ns("imputation"),
                                               "Imputation type",
                                               choices = c(c("man", "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "min", "zero"), "mixed on features"),
                                               selected = "MinDet"),
                                  shinyBS::bsTooltip(ns("id"), "Name of the column containing feature IDs", "right", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("delim"), "Set the delimiter separating the name and id column", "right", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("intensitycols"), "Choose the expression columns of your data", "right", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("remove_prefix"), "remove the prefix character of expression columns names, such as 'LFQ intensity' or 'intensity' in Maxquant result",
                                                     "right", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("remove_suffix"), "remove the suffix character of expression columns names",
                                                     "right", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("imputation"), "Choose an imputation method for peptide missing values", "right", options = list(container = "body"))
                  ),
                  bsCollapsePanel("Aggregation",
                                  style = "primary",
                                  # uiOutput(ns("fname")),
                                  selectizeInput(ns("aggregate_Peptide_Type"), "Peptides for quantification",
                                                 choices = c("Unique","Unique + Razor"),
                                                 selected = "Unique + Razor", multiple = FALSE),
                                  selectizeInput(ns("aggregation_method"), "Proteingroups aggregation algorithm",
                                                 choices = c("totalSum", "totalMean", "medianPolish","RobustSummary"), selected = "RobustSummary"),

                                  shinyBS::bsTooltip(ns("aggregate_Peptide_Type"), "Summarize proteingroups only using unique peptides or including Razor peptide","right",
                                                     options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("aggregation_method"), "The aggregation algorithm to summarize peptide quantity values","right",
                                                     options = list(container = "body"))
                  )
      ),
      actionButton(ns("aggregate_button"), "Aggregate"),
      # disabled(actionButton(ns("aggregate_button"), "Aggregate")),
      br(),br(),
      uiOutput(outputId = ns("save_aggregate_button")),
      br(),
      uiOutput(ns("aggregation_check")),
      br(),

      h4("Differential test options"),
      bsCollapse( id = "pro_sidebar",
                  open = "Columns",
                  multiple = F,
        bsCollapsePanel("Columns",
                        style = "primary",
                        # uiOutput(ns("intensitycols")),
                        # uiOutput(ns("id")),
                        shinyBS::bsTooltip(ns("name"), "Name of the column containing feature names", "right", options = list(container = "body")),

                        uiOutput(ns("peptidescol")),
                        # numericInput(ns("thr"), "allowed max.miss.num at least one condition", min = 0, max = 20, value = 0),
                        uiOutput(ns("control")),
                        uiOutput(ns("order")),
                        uiOutput(ns("test_manual")),
                        # selectizeInput(ns("norm_method"), "normalized method",
                        #                choices = c("quantiles","quantiles.robust",
                        #                            "vsn","center.mean",
                        #                            "center.median","div.mean",
                        #                            "div.median", "diff.median"),
                        #                selected = "quantiles", multiple = FALSE),

                        shinyBS::bsTooltip(ns("peptidescol"), "Choose the peptides column of your data. And this column should be numeric. It is only used to plot heatmap and custom volcano to color the peptides number 1 and 2 to help user check the accuracy of the quantification. Of course, if your data does not have this column, You can leave it blank", "right", options = list(container = "body")),
                        shinyBS::bsTooltip(ns("filt"), "The filtered columns based on", "right", options = list(container = "body")),
                        shinyBS::bsTooltip(ns("control"), "Choose your control condition", "right", options = list(container = "body")),
                        shinyBS::bsTooltip(ns("test_manual"), "Choose the contrasts that you want to test", "right", options = list(container = "body"))
        ),
        bsCollapsePanel("FDR correction",  ## DEP FDR correct options
                        style = "primary",
                        radioButtons(ns("FDR_type"),
                                     "type of FDR value",
                                     choices = c("Strimmer's qvalue(t)","Strimmer's qvalue(p)","BH","Storey's qvalue"),
                                     selected = "Strimmer's qvalue(t)"),
                        shinyBS::bsTooltip(ns("FDR_type"), "Choose the method of pvalue adjustment", "right", options = list(container = "body"))
        ),
        bsCollapsePanel("Threshold method",
                        style = "primary",
                        radioButtons(ns("threshold_method"),
                                     "Threshold method",
                                     choices = c("intersect", "curve"),
                                     selected = "intersect"),
                        shinyBS::bsTooltip(ns("threshold_method"), "Choose the method that the cutoff of significant proteins based on. [intersect] means by adjusted pvalue and log2 fold change cutoff and [curve] means by curvature and x0 fold cutoff", "right", options = list(container = "body"))
        )
      ),
      disabled(actionButton(ns("analyze"), "Analyze")),
      shinyBS::bsTooltip(ns("analyze"), "Click on it to analyze your data", "right", options = list(container = "body")),
      tags$hr(),
      tags$style(type="text/css", "#downloadData {background-color:white;color: black;font-family: Source Sans Pro}"),
      uiOutput(ns("downloadTable")),
      uiOutput(ns("downloadButton")),
      h6(),
      tags$style(type="text/css", "#Save_RData {background-color:white;color: black;font-family: Source Sans Pro}"),
      # uiOutput(ns("downloadButton_for_save_RData")),
      hr(),
      uiOutput(ns("downloadButton_for_log")),
      shinyBS::bsTooltip(ns("downloadTable"), "Choose a dataset to save, and here we offer four forms of datasets for downloading including results.txt, significant_proteins.txt, displayed_subset.txt, full_dataset.txt", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("downloadButton"), "Click on it to download the selected dataset", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("downloadButton_for_save_RData"), "Click on it to Save the RData, and when you upload the results.RData in the menuItem [Files] next time, you can get the same result as this time", "right", options = list(container = "body"))

    )
  )
}


DEP_pep_body_mod <- function(id){
  ns = NS(id)
  cat(paste0("DEP_pep_body_mod ns is",ns(""),"\n"))

  tagList(
    mainPanel(
      #** DEP tabItem----
      # tabName = "DEP_tab",
      fluidRow(tags$head(tags$style(type="text/css", "
                        #loadmessage {
                        top: 0px; left: 0px;
                        width: 100%; padding: 5px 0px 5px 0px;
                        text-align: center; font-weight: bold;
                        font-size: 100%; color: #000000;
                        background-color: #FFC1C1; z-index: 105;}")), ## 提示条的样式
               conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                tags$div("Calculating...Please wait...",id="loadmessage"))
      ),
      helpText("Please cite: "),
      fluidRow(

        conditionalPanel(condition = paste0("input['",ns("threshold_method"),"']","=='intersect'"),
                         box(numericInput(ns("p"),
                                          "adj. P value",
                                          min = 0.0001, max = 0.1, value = 0.05),
                             width = 2),
                         box(numericInput(ns("lfc"),
                                          "Log2 fold change",
                                          min = 0, max = 10, value = 1),
                             width = 2)
        ),
        conditionalPanel(condition = paste0("input['",ns("threshold_method"),"']","=='curve'"),
                         box(numericInput(ns("curvature"),
                                          "curvature",
                                          min = 0.1, max = 6, step = 0.1, value = 0.6),
                             width = 2),
                         box(numericInput(ns("x0_fold"),
                                          "x0 fold",
                                          min = 1, max = 10, value = 1),
                             width = 2)
        ),

        infoBoxOutput(ns("significantBox")),
        box(radioButtons(ns("pres"),
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
        shinyBS::bsTooltip(ns("curvature"), "Set the curvature threshold", "top", options = list(container = "body")),
        shinyBS::bsTooltip(ns("x0_fold"), "Set the x0 fold threshold, it decides the x0 ratio to the standard deviations of L2FC", "top", options = list(container = "body")),
        shinyBS::bsTooltip(ns("pres"), "The type of data scaling used for heatmap plotting. Either the fold change (contrast) or the centered log2-intensity (centered)", "top", options = list(container = "body")),
        shinyBS::bsTooltip(ns("contrasts"), "The type of contrasts that will be tested. This can be all possible pairwise comparisons (all), limited to the comparisons versus the control (control), or manually defined contrasts (manual)", "top", options = list(container = "body"))
      ),
      fluidRow(
        column(width = 7,
               box(title = "Top Table",
                   box(uiOutput(ns("select")), width = 6),
                   box(uiOutput(ns("exclude")), width = 6),
                   DT::dataTableOutput(ns("table")), width = 12),
               shinyBS::bsTooltip(ns("select"), "Select the rows that are significant of the selected contrasts", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("exclude"), "Exclude the rows that are significant of the selected contrasts", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("table"), "Table of the significant results", "top", options = list(container = "body"))
        ),
        column(width = 5,

               ## Result tabs ----
               tabBox(title = "Result Plots", width = 12,
                      # tabsetPanel(id = ns("DEP_results_tabs"),
                      tabPanel(title = "Selected Protein",
                               fluidRow(
                                 box(numericInput(ns("selected_plot_Width"),
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                                 box(numericInput(ns("selected_plot_Height"),
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                               ),
                               fluidRow(
                                 plotOutput(ns("selected_plot")),
                                 downloadButton(ns('downloadPlot'), 'Save plot')
                               ),
                               shinyBS::bsTooltip(ns("selected_plot_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("selected_plot_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Heatmap",
                               fluidRow(
                                 column(width = 6,
                                        # box(checkboxInput(ns("manual_heatmap"),
                                        #                   "Manual heatmap",
                                        #                   value = FALSE), width = 12)
                                 ),
                                 column(width = 6,
                                        box(selectizeInput(ns("colorbar"),
                                                           "colorbar",
                                                           choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                                                           selected = c("RdBu"), multiple = FALSE), width = 12)
                                 )
                               ),
                               fluidRow(
                                 box(uiOutput(ns("heatmap_cntrst")), width = 12)
                               ),

                               fluidRow(

                                 box(numericInput(ns("k"),
                                                  "Kmeans",
                                                  min = 0, max = 15, value = 6),
                                     width = 4),
                                 box(numericInput(ns("limit"),
                                                  "Color",
                                                  min = 0, max = 16, value = 4),
                                     width = 4),
                                 box(numericInput(ns("size"),
                                                  "Heat size(4-30)",
                                                  min = 4, max = 30, value = 10),
                                     width = 4)

                               ),
                               fluidRow(
                                 dropdownButton(

                                   conditionalPanel(condition = paste0("input['",ns("peptidescol"),"']"," !== ''"), # "input.peptidescol !== ''"
                                                    checkboxInput(ns("if_rowname_color"),
                                                                  "row color",
                                                                  value = FALSE)
                                   ),

                                   uiOutput(ns("chooseToshow")),

                                   numericInput(ns("row_font_size"),
                                                "row font size",
                                                min = 0, max = 15, value = 6),
                                   numericInput(ns("col_font_size"),
                                                "col font size",
                                                min = 0, max = 16, value = 10),
                                   numericInput(ns("Width"),
                                                "heatmap width",
                                                min = 1, max = 30, value = 7),
                                   numericInput(ns("Height"),
                                                "heatmap height",
                                                min = 1, max = 30, value = 10),
                                   checkboxInput(ns("cluster_columns"),
                                                 "Cluster columns",
                                                 value = TRUE),
                                   uiOutput(ns("Custom_columns_order")),
                                   uiOutput(ns("mysplit")),

                                   tooltip = tooltipOptions(title = "More parameters"),
                                   icon = icon("gear"),
                                   # status = "info",
                                   width = "300px"),
                                 br("")
                               ),
                               fluidRow(
                                 uiOutput(ns("plot")),
                                 downloadButton(ns('downloadHeatmap'), 'Save heatmap')),
                               shinyBS::bsTooltip(ns("colorbar"), "Choose the color palette of the colorbar", "top", options = list(container = "body")),
                               # shinyBS::bsTooltip(ns("manual_heatmap"), "Whether plot the manual heatmap", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("heatmap_cntrst"), "Choose the contrasts that you want to show, act when [Manual heatmap] is TRUE, can be one or more contrasts", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("k"), "Set the number of k-means clusters", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("limit"), "Set the outer limits of the color scale", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("size"), "Set the height of the plot region", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("if_rowname_color"), "Whether set the color of the row names,(purple: Peptides == 1, blue: Peptides == 2)", "top", options = list(container = "body")),
                               # shinyBS::bsTooltip(ns("if_chooseToshow"), "Whether only show the rows that you want to show", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("chooseToshow"), "Choose the rows that you want to show", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("row_font_size"), "Set the size of row labels", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("col_font_size"), "Set the size of column labels", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Width"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Height"), "Height of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("cluster_columns"), "Whether make cluster on columns", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Custom_columns_order"), "Order of column, act when [Cluster columns] is FALSE", "top", options = list(container = "body")),
                               # shinyBS::bsTooltip(ns("if_mysplit"), "Whether order the splices manually", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("mysplit"), "The ordered splices of the heatmap (from top to bottom), only act when both [If mysplit] and [Kmeans] are TRUE", "top", options = list(container = "body"))
                      ),
                      tabPanel(title = "Custom Volcano",
                               fluidRow(
                                 column(width = 4,
                                        colourpicker::colourInput(inputId = ns("up_color"),
                                                                  label = "up",
                                                                  showColour = "both",
                                                                  palette = "square",
                                                                  value = "#CD8500")
                                 ),
                                 column(width = 4,
                                        colourpicker::colourInput(inputId = ns("down_color"),
                                                                  label = "down",
                                                                  showColour = "both",
                                                                  palette = "square",
                                                                  value = "#0099CC")
                                 ),
                                 column(width = 4,
                                        colourpicker::colourInput(inputId = ns("stable_color"),
                                                                  label = "not significant",
                                                                  showColour = "both",
                                                                  palette = "square",
                                                                  value = "#BEBEBE")
                                 )
                               ),
                               fluidRow(
                                 box(uiOutput(ns("Volcano_cntrst")), width = 8),
                                 box(selectizeInput(ns("labelWay"),
                                                    "label way",
                                                    choices = c("all", "up", "down"),
                                                    selected = c("all significant"), multiple = FALSE),width = 4)
                               ),
                               fluidRow(
                                 dropdownButton(
                                   # checkboxInput(ns("if_label_rectangle"),
                                   #               "label with rectangle",
                                   #               value = FALSE),
                                   conditionalPanel(condition = paste0("input['",ns("peptidescol"),"']"," !== ''"), #"input.peptidescol !== ''",#&& input.peptidescol !== '' typeof input.peptidescol !== null
                                                    checkboxInput(ns("if_peptide_color"),
                                                                  "peptide color",
                                                                  value = FALSE),
                                                    # column(width = 6,
                                                    uiOutput(ns("Peptides1")),
                                                    uiOutput(ns("Peptides2"))
                                   ),

                                   uiOutput(ns("selected_proteins")),
                                   numericInput(ns("showNum"),
                                                "show number",
                                                min = 0, max = 100000, value = 20),
                                   numericInput(ns("fontSize"),
                                                "font size",
                                                min = 0, max = 20, value = 5),
                                   numericInput(ns("dotsize"),
                                                "dot size",
                                                min = 0, max = 20, value = 2),
                                   numericInput(ns("custom_volcano_Width"),
                                                "width",
                                                min = 1, max = 30, value = 7),
                                   numericInput(ns("custom_volcano_Height"),
                                                "height",
                                                min = 1, max = 30, value = 7),
                                   tooltip = tooltipOptions(title = "More parameters"),
                                   icon = icon("gear"),
                                   width = "300px"),
                                 conditionalPanel(condition = paste0("input['",ns("threshold_method"),"']"," == 'intersect'"), #"input.threshold_method == 'intersect'",
                                                  box(checkboxInput(ns("P_adj"),
                                                                    "Adjusted p values",
                                                                    value = FALSE),
                                                      checkboxInput(ns("Same_width"),
                                                                    "Same width",
                                                                    value = FALSE),
                                                      width = 4)
                                 ),
                                 conditionalPanel(condition = paste0("input['",ns("threshold_method"),"']"," == 'curve'"), #"input.threshold_method == 'curve'",
                                                  box(checkboxInput(ns("Same_width"),
                                                                    "Same width",
                                                                    value = FALSE),
                                                      width = 4)
                                 )

                               ),
                               fluidRow(
                                 plotOutput(ns("custom_volcano"), height = 600),
                                 downloadButton(ns('download_custom_volcano'), 'Save custom volcano')
                               ),
                               shinyBS::bsTooltip(ns("up_color"), "Set the color of the points those are up-regulated", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("down_color"), "Set the color of the points those are down-regulated", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("stable_color"), "Set the color of the points those are not significant", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Volcano_cntrst"), "Choose the contrast to plot", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("labelWay"), "Choose the way of adding labels, one of all significant, up, down, and selected proteins", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("if_label_rectangle"), "Whether add rectangle underneath the text, making it easier to read", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("if_peptide_color"), "Whether set the color of the points,(default, purple: Peptides == 1, blue: Peptides == 2)", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Peptides1"), "Set the color of the points which the number of Peptides == 1", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Peptides2"), "Set the color of the points which the number of Peptides == 2", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("stroke"), "Set the thickness of black line around the point", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("selected_proteins"), "Choose the point labels to show, act when [label Way] is selected proteins", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("showNum"), "Set the number of the labels to add, act when [label way] is significant, up or down. Note that, when it is larger than the number of significant result, it is seted to the number of significant result", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("fontSize"), "Set the font size of the label", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("dotsize"), "Set the dot size of the label", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("custom_volcano_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("custom_volcano_Height"), "Height of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("P_adj"), "Whether or not to use adjusted p values", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(ns("Same_width"), "Whether the x axis to have the same width from 0", "top", options = list(container = "body"))
                      )
                      # )
               ),

               ## QC tabs ----
               # tabBox(title = "QC Plots", width = 12,
               tabsetPanel(id = ns("DEP_QC_tabs"),

                           tabPanel(title = "Pca plot",
                                    fluidRow(
                                      box(selectizeInput(ns("Indicate"),
                                                         "Color and shape",
                                                         choices = c("condition", "replicate", "Condition", "Replicate"),
                                                         selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                      box(checkboxInput(ns("if_square"),
                                                        "if square",
                                                        value = FALSE),
                                          width = 3)),
                                    fluidRow(
                                      box(uiOutput(ns("pca_top_n")),width = 9)
                                    ),
                                    fluidRow(
                                      box(numericInput(ns("pca_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7), width = 6),
                                      box(numericInput(ns("pca_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 6)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("pca"), height = 600),
                                      downloadButton(ns('downloadPca'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("Indicate"), "Set the color, shape and facet_wrap of the plot", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("if_square"), "Whether x limit is equal to y limit", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("pca_top_n"), "Set the number of top variable proteins to consider", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("pca_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("pca_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "UMAP plot",
                                    fluidRow(
                                      box(selectizeInput(ns("umap_Indicate"),
                                                         "Color and shape",
                                                         choices = c("condition", "replicate", "Condition", "Replicate"),
                                                         selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                      box(checkboxInput(ns("umap_if_square"),
                                                        "if square",
                                                        value = FALSE),
                                          width = 6)),
                                    fluidRow(
                                      box(uiOutput(ns("umap_n_neighbors")),width = 10)
                                    ),
                                    fluidRow(
                                      box(numericInput(ns("umap_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7), width = 6),
                                      box(numericInput(ns("umap_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 6)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("umap"), height = 600),
                                      downloadButton(ns('downloadUMAP'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("umap_Indicate"), "Set the color, shape and facet_wrap of the plot", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("umap_if_square"), "Whether x limit is equal to y limit", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("umap_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("umap_Height"), "Height of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("umap_n_neighbors"), "number of nearest neighbors for umap", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "t-SNE plot",
                                    fluidRow(
                                      box(selectizeInput(ns("Tsne_Indicate"),
                                                         "Color and shape",
                                                         choices = c("condition", "replicate", "Condition", "Replicate"),
                                                         selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                      box(checkboxInput(ns("Tsne_if_square"),
                                                        "if square",
                                                        value = FALSE),
                                          width = 6)
                                    ),
                                    fluidRow(
                                      box(uiOutput(ns("Tsne_perplexity")),width = 8),
                                      box(numericInput(ns("Tsne_theseed"),
                                                       "set seed for t-SNE",
                                                       value = 42),
                                          width = 4)
                                    ),
                                    fluidRow(
                                      box(numericInput(ns("Tsne_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7), width = 6),
                                      box(numericInput(ns("Tsne_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 6)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("Tsne"), height = 600),
                                      downloadButton(ns('downloadTSNE'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("Tsne_Indicate"), "Set the color, shape and facet_wrap of the plot", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Tsne_if_square"), "Whether x limit is equal to y limit", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Tsne_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Tsne_Height"), "Height of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Tsne_perplexity"), "perplexity for t-SNE", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Statistical plot",

                                    fluidRow(
                                      box(selectInput(ns("stastical_plot_x"),
                                                      "plot statistics",
                                                      choices = c("diff","t.stastic","p.val","p.adj"),
                                                      multiple = T,
                                                      selected = c("diff","t.stastic","p.val","p.adj")),
                                          width = 6),

                                      box(uiOutput(ns("plot_stat_contrasts")),

                                          width = 6)
                                    ),
                                    fluidRow(
                                      box(numericInput(ns("stastical_plot_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7), width = 6),
                                      box(numericInput(ns("stastical_plot_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 6)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("stastical_plot"),height=600),
                                      downloadButton(ns('downloadstastical_plot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("stastical_plot_x"), "Select the statistic to plot", "top", options = list(container = "body")),
                                    # shinyBS::bsTooltip(ns("stastical_plot_y"), "Set the vector that y axis represents", "top", options = list(container = "body")),
                                    # shinyBS::bsTooltip(ns("stastical_plot_type"), "Set the plot type: x-y represents that X-axis represents your selected x, and y-axis represents your selected y; histogram represents the histogram plot of your selected x", "top", options = list(container = "body")),
                                    # shinyBS::bsTooltip(ns("Volcano_cntrst_stastical_plot"), "Choose the contrast to plot", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("stastical_plot_Width"), "Width of the figure to export when save", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("stastical_plot_Height"), "Height of the figure to export when save", "top", options = list(container = "body"))
                           ),
                           #tabPanel(title = "Pca plot",
                           #         plotOutput("pca", height = 600),
                           #         downloadButton('downloadPca', 'Save')
                           #)
                           tabPanel(title = "meanSdPlot",
                                    fluidRow(
                                      box(numericInput(ns("meanSdPlot_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7), width = 6),
                                      box(numericInput(ns("meanSdPlot_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 6)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("meanSdPlot"), height = 600),
                                      downloadButton(ns('downloadmeanSdPlot'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("meanSdPlot_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("meanSdPlot_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Pearson correlation",
                                    fluidRow(
                                      box(selectizeInput(ns("Pearson_pal"),
                                                         "color panel",
                                                         choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Blues",  "BuGn", "BuPu", "GnBu", "Greens" , "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
                                                         selected = c("PRGn"), multiple = FALSE),
                                          checkboxInput(ns("Pearson_pal_rev"),
                                                        "pal rev",
                                                        value = FALSE), width = 4),
                                      box(numericInput(ns("Pearson_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7),
                                          numericInput(ns("Pearson_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 4),
                                      box(numericInput(ns("Pearson_lower"),
                                                       "lower",
                                                       min = -1, max = 1, value = -1),
                                          numericInput(ns("Pearson_upper"),
                                                       "upper",
                                                       min = -1, max = 1, value = 1), width = 4)
                                    ),
                                    fluidRow(
                                      box(checkboxInput(ns("add_values_for_DEP_person"),
                                                        "Add values",
                                                        value = FALSE),
                                          numericInput(ns("value_size_for_DEP_person"),
                                                       "Value size",
                                                       min = 1, max = 30, value = 10),
                                          numericInput(ns("value_digits_for_DEP_person"),
                                                       "Value digits",
                                                       min = 1, max = 30, value = 2), width = 12, collapsible = TRUE, collapsed = TRUE)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("Pearson_correlation"), height = 600),
                                      downloadButton(ns('download_Pearson_correlation'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("Pearson_pal"), "Set the color panel (from RColorBrewer)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Pearson_pal_rev"), "Whether or not to invert the color palette", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Pearson_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Pearson_Height"), "Height of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Pearson_lower"), "Set the lower limit of the color scale", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Pearson_upper"), "Set the upper limit of the color scale", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Gower's distance",
                                    fluidRow(
                                      box(selectizeInput(ns("Gower_pal"),
                                                         "color panel",
                                                         choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Blues",  "BuGn", "BuPu", "GnBu", "Greens" , "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
                                                         selected = c("YlOrRd"), multiple = FALSE),
                                          checkboxInput(ns("Gower_pal_rev"),
                                                        "pal rev",
                                                        value = TRUE), width = 6),
                                      box(numericInput(ns("Gower_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7),
                                          numericInput(ns("Gower_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 6)
                                    ),
                                    fluidRow(
                                      box(checkboxInput(ns("add_values_for_DEP_gower"),
                                                        "Add values",
                                                        value = FALSE),
                                          numericInput(ns("value_size_for_DEP_gower"),
                                                       "Value size",
                                                       min = 1, max = 30, value = 10),
                                          numericInput(ns("value_digits_for_DEP_gower"),
                                                       "Value digits",
                                                       min = 1, max = 30, value = 2), width = 12, collapsible = TRUE, collapsed = TRUE)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("Gowers_distance"), height = 600),
                                      downloadButton(ns('download_Gowers_distance'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("Gower_pal"), "Set the color panel (from RColorBrewer)", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Gower_pal_rev"), "Whether or not to invert the color palette", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Gower_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Gower_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Sample CVs",
                                    fluidRow(
                                      box(numericInput(ns("Sample_CVs_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7),
                                          width = 6),
                                      box(
                                        numericInput(ns("Sample_CVs_Height"),
                                                     "height",
                                                     min = 1, max = 30, value = 7), width = 6
                                      )
                                    ),
                                    fluidRow(
                                      plotOutput(ns("Sample_CVs"), height = 600),
                                      downloadButton(ns('download_Sample_CVs'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("Sample_CVs_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("Sample_CVs_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Protein Numbers",
                                    fluidRow(
                                      box(numericInput(ns("numbers_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7), width = 6),
                                      box(numericInput(ns("numbers_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 6)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("numbers"), height = 600),
                                      downloadButton(ns('downloadNumbers'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("numbers_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("numbers_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Sample coverage",
                                    fluidRow(
                                      box(numericInput(ns("coverage_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7), width = 6),
                                      box(numericInput(ns("coverage_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 6)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("coverage"), height = 600),
                                      downloadButton(ns('downloadCoverage'), 'Save')

                                    ),
                                    shinyBS::bsTooltip(ns("coverage_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("coverage_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Normalization",
                                    fluidRow(
                                      box(numericInput(ns("norm_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7), width = 6),
                                      box(numericInput(ns("norm_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 6)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("norm"), height = 600),
                                      downloadButton(ns('downloadNorm'), 'Save')

                                    ),
                                    shinyBS::bsTooltip(ns("norm_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("norm_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           tabPanel(title = "Missing values - Quant",
                                    fluidRow(
                                      box(numericInput(ns("detect_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7), width = 6),
                                      box(numericInput(ns("detect_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 6)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("detect"), height = 600),
                                      downloadButton(ns('downloadDetect'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("detect_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("detect_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                           ),
                           # tabPanel(title = "Missing values - Heatmap",  ## It might is too large to draw a peptide missing values heatmap
                           #          fluidRow(
                           #            box(
                           #              checkboxInput(ns("cluster_columns_for_missval_heatmap"),
                           #                            "Cluster columns",
                           #                            value = TRUE), width = 4),
                           #            box(uiOutput(ns("Custom_columns_order_for_missval_heatmap")),
                           #                width = 8)
                           #          ),
                           #          fluidRow(
                           #            box(numericInput(ns("missval_heatmap_Width"),
                           #                             "width",
                           #                             min = 1, max = 30, value = 7), width = 6),
                           #            box(numericInput(ns("missval_heatmap_Height"),
                           #                             "height",
                           #                             min = 1, max = 30, value = 10), width = 6)),
                           #          fluidRow(
                           #            plotOutput(ns("missval"), height = 600),
                           #            downloadButton(ns('downloadMissval'), 'Save')
                           #          ),
                           #          shinyBS::bsTooltip(ns("cluster_columns_for_missval_heatmap"), "Whether make cluster on columns", "top", options = list(container = "body")),
                           #          shinyBS::bsTooltip(ns("Custom_columns_order_for_missval_heatmap"), "Order of column, act when [Cluster columns] is FALSE", "top", options = list(container = "body")),
                           #          shinyBS::bsTooltip(ns("missval_heatmap_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                           #          shinyBS::bsTooltip(ns("missval_heatmap_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                           # ),
                           tabPanel(title = "Imputation",
                                    fluidRow(
                                      box(numericInput(ns("imputation_Width"),
                                                       "width",
                                                       min = 1, max = 30, value = 7), width = 6),
                                      box(numericInput(ns("imputation_Height"),
                                                       "height",
                                                       min = 1, max = 30, value = 7), width = 6)
                                    ),
                                    fluidRow(
                                      plotOutput(ns("imputation"), height = 600),
                                      downloadButton(ns('downloadImputation'), 'Save')
                                    ),
                                    shinyBS::bsTooltip(ns("imputation_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                    shinyBS::bsTooltip(ns("imputation_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                           )
                           # )
               )
        )
      )
    )
  )
}


DEP_pep_server_module <- function(id){
  moduleServer(
    id,
    function(input,output,session){
      ns = session$ns

      observe({
        cat("Check aggregate input... ")
        req(input$id, input$name, input$intensitycols, input$imputation, input$aggregate_Peptide_Type)
        cat("pass\n")
        # shinyjs::enable("aggregate_button") ## check all input for aggregation is ready. is exist

        cat("Check anlyze input... ")
        req(input$peptidescol)
        req(data_QF())
        cat("pass\n")
        shinyjs::enable("analyze") ## check all input for Analyze is ready. and aggregation is finish
      })

      ### UI functions of aggregation or analyze options --------------------------------------------------------
      output$aggregation_check <- renderUI({
        if(!is.null(data_QF())){
          HTML("Aggregation finished! please perform differential analysis in the below pannels")
        }else{
          HTML("Please perform aggregation first in the above panels!")
        }
      })


      output$name <- renderUI({
        if(input$uploadmode == "fromTable"){
          shiny::validate(need(!is.null(data()), "Please upload quantity file in 'Files' firstly!"))
        }else{
          validate(need(!is.null(data()), "Please upload log in 'Files' firstly!"))
        }

        # shiny::validate(need(!is.null(data_QF()), " "))

        selectizeInput(inputId = ns("name"),
                       label = "Protein name column",
                       choices = colnames(peptide()),#
                       selected = "Gene.names")
      })

      output$id <- renderUI({
        shiny::validate(need(!is.null(peptide()), " "))
        # shiny::validate(need(!is.null(data_QF()), " "))
        selectizeInput(inputId = ns("id"),
                       label = "Protein ID column",
                       choices = colnames(peptide()),
                       selected = "Proteins")
      })

      output$peptidescol <- renderUI({
        shiny::validate(need(!is.null(peptide()), " "))
        shiny::validate(need(!is.null(data_QF()), " "))

        if(".n" %in% colnames(rowData(data_QF()[["protein"]]))){
          return(selectizeInput(inputId = ns("peptidescol"),
                                label = "peptide column",
                                choices=c(colnames(rowData(data_QF()[["protein"]])), ""),
                                selected = ".n"))
        }else{
          return(selectizeInput(inputId = ns("peptidescol"),
                                label = "peptide column",
                                choices=c("", colnames(rowData(data_QF()[["protein"]]))),
                                selected = NULL))
        }
      })

      output$intensitycols <- renderUI({
        LFQcols = grep("^Intensity\\.|^intensity\\.",colnames(peptide()),value = T)
        if(length(LFQcols)==0) grep("^Intensity",colnames(peptide()),value = T)
        if(length(LFQcols)==0) grep("^intensity",colnames(peptide()),value = T)
        if(length(LFQcols)==0) grep("^intensity",colnames(peptide()),value = T)

        bsCollapse(id = session$ns("IntensityCols"), open = "Expression columns",
                   bsCollapsePanel("Expression columns",
                                   selectizeInput(ns("intensitycols"),
                                                  label = "Choose the expression columns",
                                                  choices=colnames(peptide()),
                                                  multiple = TRUE,
                                                  selected = grep("^Intensity\\.|^intensity\\.",LFQcols,value = T), width = '100%'),
                                   checkboxInput(inputId = session$ns("remove_prefix"),"remove prefix of label", value = T),
                                   checkboxInput(inputId = session$ns("remove_suffix"),"remove suffix of label", value = T)
                                   , style = "primary"))


      })

      output$filt <- renderUI({
        shiny::validate(need(!is.null(peptide()), " "))
        # shiny::validate(need(!is.null(data_QF()), " "))
        selectizeInput(inputId = ns("filt"),
                       label = "Filter on columns" ,
                       choices = c("", colnames(peptide())),
                       multiple = TRUE,
                       selected = c("Reverse","Potential.contaminant"))
      })


      output$control <- renderUI({
        if(input$uploadmode == "fromTable"){
          validate(need(!is.null(input$file1), "Please upload quantity file in Files"))
        }else{
          validate(need(!is.null(input$file_log), "Please upload log in Files"))
        }

        validate(need(!is.null(input$intensitycols), "Please select the Expression columns"))
        validate(need(length(input$intensitycols) > 1 , "More expression columns is required"))
        validate(need(!(is.null(data_QF()) && is.null(input$resultRData)), "Please finish aggregate first."))

        if (input$anno == "columns" & !is.null(data_QF()) & input$contrasts == "control") {
          my_data <- data_QF()
          selectizeInput(ns("control"), "Control",
                         choices = colData(my_data[["protein"]])$condition %>% unique,
                         selected = NULL)

        } else if (input$anno == "expdesign" & !is.null(expdesign()) & input$contrasts == "control"){
            selectizeInput(ns("control"), "Control",
                           choices = make.names(expdesign()$condition)%>% unique() )

        }else if (input$anno == "logexpdesign" & !is.null(upload_log()) & !is.null(expdesign()) ){
          # selected111 <<- upload_log()$inputVals()[["control"]]
          # condition11 <<- make.names(expdesign()$condition) %>% unique()
          selectizeInput(ns("control"), "Control",
                         choices = make.names(expdesign()$condition)%>% unique(),
                         selected = upload_log()$inputVals()[["control"]]
                         )
        }
      })

      output$order <- renderUI({
        if(input$uploadmode == "fromTable"){
          validate(need(!is.null(input$file1), "Please upload quantity file in Files"))
        }else{
          validate(need(!is.null(input$file_log), "Please upload log in Files"))
        }

        validate(need(!is.null(input$intensitycols), "Please select the Expression columns"))
        validate(need(length(input$intensitycols) > 1 , "More expression columns is required"))
        shiny::validate(need(!is.null(peptide()), " "))
        shiny::validate(need(!(is.null(data_QF()) && is.null(input$resultRData)), "Please finish aggregate first."))

        # rdataFile <- input$resultRData
        # if(is.null(rdataFile)){
        #   data <- data_QF()[["protein"]]
        # }else{
        #   load(rdataFile)
        #   data <- data_QF(data)
        # }
        # condition = colData(data)$condition
        condition = filt0()@colData$condition %>% unique()
        selectizeInput(session$ns("order"),
                       "Set the groups order",
                       choices = condition, multiple = TRUE,
                       selected = NULL, width = '100%')
      })

      output$test_manual <- renderUI({
        validate(need(!is.null(input$file1), ""))
        if(!is.null(data_QF()) & input$contrasts == "manual"){

          test_manual_name <- colData(data_QF()[["protein"]])$condition %>% unique
          test_manual_name <- cbind(combn(test_manual_name,2),combn(test_manual_name,2, FUN = rev))
          test_manual_name <- apply(test_manual_name, 2, function(i){paste(i[1], i[2], sep = "_vs_")})
          selectizeInput(ns("test_manual"), "Manual test",
                         choices=test_manual_name, selected = NULL, multiple = TRUE)
        }
      })

      output$Peptides1 <- renderUI({
        if(input$if_peptide_color) {
          colourpicker::colourInput(inputId = ns("Peptides_1"),
                                    label = "Peptides1",
                                    showColour = "both",
                                    palette = "square",
                                    value = "#A020F0")
        }
      })

      output$Peptides2 <- renderUI({
        if(input$if_peptide_color) {
          colourpicker::colourInput(inputId = ns("Peptides_2"),
                                    label = "Peptides2",
                                    showColour = "both",
                                    palette = "square",
                                    value = "#0000FF")#
        }
      })


      # observeEvents upon uploadmode and log_file. fromLog -----------------------------------------------------------------

      ## updata options according uploadmode and upload_log
      Columns_to_updata = reactiveVal(value = 0)
      Columns_rendered = reactiveVal(value = 0)

      observeEvent(upload_log(),{


        cat("Update options 1 \n")
        Option_to_updata <- c("imputation", "FDR_type", "threshold_method","contrasts")
        for(i in Option_to_updata){
          # message("update ",i," to '", upload_log()$inputVals()[[i]],"'\n", sep = "")
          updateRadioButtons(session = session, inputId = i,
                             selected = upload_log()$inputVals()[[i]] )
        }

        updatePrettySwitch(session = session, inputId = "imp_norm_order",
                           value = upload_log()$inputVals()[["imp_norm_order"]] )

        Option_to_updata2 <- c("remove_prefix", "remove_suffix")
        for(i in Option_to_updata2){
          updateCheckboxInput(session = session, inputId = i,
                              value = upload_log()$inputVals()[[i]] )
        }

        Option_to_updata3 <- c("norm_method", "aggregate_Peptide_Type", "aggregation_method")
        for(i in Option_to_updata3){
          updateSelectizeInput(session = session, inputId = i,
                               selected = upload_log()$inputVals()[[i]] )
        }

        message("Update options 1 finished \n")

        data_QF <- data_QF({
          data_QF = upload_log()$resultVals()$data
          # data_QF11 <<- data_QF
        })

        Columns_to_updata = Columns_to_updata(1)

      },ignoreNULL = T,ignoreInit = T)

      observeEvent(
        input$peptidescol,{
          Columns_rendered = Columns_rendered(1)
        },ignoreNULL = T,ignoreInit = T
      )

      observeEvent(
        c(Columns_to_updata(),Columns_rendered()),
        {
          if(Columns_to_updata() == 1 & Columns_rendered() == 1){
            cat("Update options in 'Colmuns' \n")

            column_opt = c("name", "id", "delim", "intensitycols",
                           "filt", "order", "control",
                           "norm_method", "aggregate_Peptide_Type", "aggregation_method"
                           )
            for(i in column_opt){
              # message("update ",i," to '", upload_log()$inputVals()[[i]],"'\n", sep = "")
              updateSelectizeInput(session = session, inputId = i,
                                   selected = upload_log()$inputVals()[[i]] )
            }

            Threshold_opt = c("lfc","p","curvature","x0_fold")
            for(i in column_opt){
              # message("update ",i," to '", upload_log()$inputVals()[[i]],"'\n", sep = "")
              updateNumericInput(session = session, inputId = i,
                                 value = upload_log()$inputVals()[[i]] )
            }

            Columns_rendered = Columns_rendered(0) # back to 0 once restore the options
            Columns_to_updata = Columns_to_updata(0) # back to 0 once restore the options
            message("Update 'Colmuns' finnished \n")
          }

        },
        ignoreNULL = T,ignoreInit = T
      )

      ## update anno option according uploadmode
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
                             choices = list(
                               # "Parse from columns" = "columns",
                               "Use ExpDesign in log" = "logexpdesign"),
                             selected = "logexpdesign"
          )
        }

      })

      # if intensitycols is changed from log, update the anno to columns
      observeEvent(input$intensitycols,{
        if(input$uploadmode == "fromLog" && !identical(input$intensitycols, upload_log()$inputVals()[["intensitycols"]]) )
          updateRadioButtons(session = session,
                             inputId = "anno",
                             label = "Sample annotation",
                             choices = list("Parse from columns" = "columns",
                                            "Use ExpDesign in log" = "logexpdesign"),
                             selected = "columns"
          )
      })

      # observeEvents upon threshold_method -----------------------------------------------------------------
      observe(
        if(input$threshold_method == "curve") {
          insertTab(inputId = ns("DEP_QC_tabs"),
                    tabPanel(title = "Normal distribution - Fit",
                             fluidRow(
                               box(uiOutput(ns("cntrst_norm_distribution")), width = 8)
                             ),
                             fluidRow(
                               box(numericInput(ns("norm_distribution_Width"),
                                                "width",
                                                min = 1, max = 30, value = 7), width = 6),
                               box(numericInput(ns("norm_distribution_Height"),
                                                "height",
                                                min = 1, max = 30, value = 7), width = 6)
                             ),
                             fluidRow(
                               plotOutput(ns("norm_distribution_plot"),height=600),
                               downloadButton(ns('download_norm_distribution_plot'), 'Save')
                             ),
                             shinyBS::bsTooltip(ns("cntrst_norm_distribution"), "Choose the contrast to plot", "top", options = list(container = "body")),
                             shinyBS::bsTooltip(ns("norm_distribution_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                             shinyBS::bsTooltip(ns("norm_distribution_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                    ),

                    target = "Pca plot", position = "before"
          )
        }
      )

      observe(
        if(input$threshold_method == "intersect") {
          removeTab(inputId = ns("DEP_QC_tabs"), target = "DEP curve 2")
        }
      )

      # observeEvents upon imp_norm_order -----------------
      observe(
        if(!input$imp_norm_order){
          updateSelectizeInput(session = session,
                               inputId = "norm_method",
                               choices = c("quantiles","quantiles.robust", "diff.median","vsn"),
                               selected = "quantiles"
                               )
        }else if(input$imp_norm_order){
          updateSelectizeInput(session = session,
                               inputId = "norm_method",
                               choices = c("quantiles","quantiles.robust", "diff.median","vsn"),
                               selected = "quantiles.robust"
          )
        }
      )

      ### check cutoff input ### --------------------------------------------------------
      iv <- InputValidator$new()
      iv$add_rule("p", sv_between(0.0001, 0.1))
      iv$add_rule("lfc", sv_between(0, 10))
      iv$add_rule("curvature", sv_between(0.1, 6))
      iv$add_rule("x0_fold", sv_between(1, 10))
      iv$enable()

      ### Reactive functions of peptide ### --------------------------------------------------
      upload_log <- reactive({
        logFile <- input$file_log

        if( !is.null(logFile) & input$uploadmode == "fromLog" ){

          cat("Read infile_log\n")
          upload_log <- readRDS(logFile$datapath)

          if( upload_log$appVersion == DEP2:::app_version && inherits(upload_log,"pgpep_log") ){
            sendSweetAlert(
              session = shiny::getDefaultReactiveDomain(),
              title = "Load log successfully!",
              text = "Pass version check and load.",
              type = "success"
            )
            message("Read successfully\n")
            return(upload_log)

          }else if(!inherits(upload_log,"pgpep_log")){
            sendSweetAlert(
              session = shiny::getDefaultReactiveDomain(),
              title = "Error !",
              text = "The uploaded file is not a log for peptide-group aggregating analysis, unload the log.\n
              Please ensure that: The log file is a log from protein-group analysis. And the app version is correct,
              DEP2 app is only compatible the log from same version. Old log may not be used on a updated app.
              ",
              type = "warning"
            )
            return(NULL)

          }else if(upload_log$appVersion != DEP2:::app_version){
            sendSweetAlert(
              session = shiny::getDefaultReactiveDomain(),
              title = "Error !",
              text = "The version of uploaded file don't match current app version.\n
              Please ensure that: The log file is a log from peptide-group aggregating analysis. And the app version is correct,
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

      expdesign <- reactive({

        # from log
        if(input$uploadmode == "fromLog" && !is.null(peptide()) ){
          upload_log = upload_log()

          if(input$anno == "logexpdesign"){
            message("Load expdesign from log")
            expdesign <- upload_log$resultVals()$expdesign

          }else if(input$anno == "columns"){
            cols <- which(colnames(peptide()) %in% input$intensitycols)

            if(length(cols) == 0){
              return(NULL)
            }else{
              message("Creat expdesign from colnames parse")
              label <- colnames(peptide())[cols]
              expdesign <- get_exdesign_parse(label, mode = "delim", sep = "_",
                                              remove_prefix = input$remove_prefix, remove_suffix = input$remove_suffix)
            }
          }
        }else if(input$uploadmode == "fromLog" && is.null(peptide())){
          return(NULL)
        }

        # from table
        if(input$uploadmode == "fromTable"){
          inFile <- input$file2
          if ( input$anno == "columns" & (!is.null(peptide())) ){
            cols <- which(colnames(peptide()) %in% input$intensitycols)

            if(length(cols) == 0){
              return(NULL)
            }else{
              message("Creat expdesign from colnames parse")
              label <- colnames(peptide())[cols]
              expdesign <- get_exdesign_parse(label, mode = "delim", sep = "_",
                                              remove_prefix = input$remove_prefix, remove_suffix = input$remove_suffix)
            }
          }else if(input$anno == "expdesign" & (!is.null(inFile)) & (!is.null(peptide())) ){
            message("Load expdesign from table")
            expdesign <- read.csv(inFile$datapath, header = TRUE,
                                  sep = "\t", stringsAsFactors = FALSE) %>%
              mutate(id = row_number())
          }else{
            return(NULL)
          }
        }

        # my_expdesign <<- expdesign
        return(expdesign)


        # inFile <- input$file2
        # if (is.null(inFile) & (!is.null(data_QF())) ){
        #   cols <- which(colnames(data_QF()) %in% input$intensitycols)
        #   if(length(cols) == 0){
        #     return(NULL)
        #   }else{
        #     label <- colnames(data_QF())[cols]
        #     expdesign <- get_exdesign_parse(label, mode = "delim", sep = "_",
        #                                     remove_prefix = input$remove_prefix, remove_suffix = input$remove_suffix)
        #
        #   }
        # }else{
        #   print(inFile)
        #   print(data_QF)
        #   read.csv(inFile$datapath, header = TRUE,
        #            sep = "\t", stringsAsFactors = FALSE) %>%
        #     mutate(id = row_number())
        # }
      })

      peptide <- reactive({
        cat("peptide1 \n")
        # from log
        upload_log = upload_log()
        if(!is.null(upload_log) & input$uploadmode == "fromLog"){
          my_data <- upload_log$resultVals()$peptide
          my_data111 <- my_data
          return(my_data)
        }else if( is.null(upload_log) & input$uploadmode == "fromLog"){
          # cat("Initiate: log file is null \n")
          return(NULL)
        }

        # from table
        inFile <- input$file1 ## the peptides file input
        if (is.null(inFile))
          return(NULL)
        my_data <- data.table::fread(inFile$datapath, header = TRUE,
                                     sep = "\t", stringsAsFactors = FALSE, integer64 = "numeric") %>%
          mutate(Pepid = row_number())
        coln.int64 <- names(which(sapply(my_data, bit64::is.integer64))) ## if read in as int64, transform to numeric
        if (length(coln.int64) > 0L)
          my_data[, c(coln.int64) := lapply(.SD, as.numeric), .SDcols = coln.int64]
        my_data = as.data.frame(my_data)
        ## replace the [xxx] characters in colnames
        colnames(my_data) = gsub("^\\[(.*)\\] ","",colnames(my_data))
        colnames(my_data) = make.names(colnames(my_data))


        cat("Read in peptide \n")
        return(my_data)

      })

      ## Creat a Qfeatures objects by make_pe
      pe <- reactive({
        # rdataFile <- input$resultRData
        # if(is.null(rdataFile)){
        Peptide <- peptide()
        need(all(input$intensitycols %in% colnames(Peptide)) & length(input$intensitycols) > 1 , "Please select the Expression columns in the sidebar")
        Ecols <- which(colnames(Peptide) %in% input$intensitycols)

        ## make pe
        # if (input$anno == "columns") {
        #   pe <- DEP2::make_pe_parse(Peptide, columns = Ecols,
        #                             # protein_col = "Pepid",
        #                             mode = "delim", sep = "_",
        #                             remove_prefix = input$remove_prefix,
        #                             log2transform = T)
        # }
        # if (input$anno == "expdesign") {
        #   pe <- DEP2::make_pe(Peptide = Peptide, columns = Ecols,
        #                       # protein_col = "Pepid",
        #                       expdesign = expdesign(),
        #                       log2transform = T)
        #   colData(pe)$replicate = as.character(colData(se)$replicate)
        # }

        # expdesign111 <<- expdesign()
        pe <- DEP2::make_pe(Peptide = Peptide, columns = Ecols,
                            # protein_col = "Pepid",
                            expdesign = expdesign(),
                            log2transform = T)
        cat("Construct object finish \n")
        return(pe)

        # }else{
        #   load(file = rdataFile$datapath)
        #   return(pe)
        # }
      })

      ## peptide-level filter,normalization, imputation.
      ### filter
      pe_filt <- reactive({
        pe <- pe()
        cat("11")
        if(length(input$filt)==0 || (length(input$filt)==1 && input$filt == "")){
          filtcols = NULL
          filter_formula = NULL
        }else{
          cat("22")
          filtcols = input$filt
          filtcols = filtcols[filtcols!=""]
          for(i in 1:length(filtcols)){
            filtcol = filtcols[i]
            if(i == 1){
              filter_formula = paste( "~(is.na(",filtcol,")|",filtcol,"=='')")
            }else{
              filter_formula = paste(filter_formula, paste( "&(is.na(",filtcol,")|",filtcol,"=='')"))
            }
          }
          filter_formula = as.formula(filter_formula)
        }

        cat("33")
        pe_filt <- DEP2::filter_pe(pe, thr = input$thr,
                                   filter_formula = filter_formula)

        cat("Filter finished.\n")
        return(pe_filt)
      })

      ### imputation and normalization.The order could change
      pe_norm <- reactive({
        # rdataFile <- input$resultRData
        # if(is.null(rdataFile)){
        upload_log = upload_log()
        if( !is.null(upload_log) & input$uploadmode == "fromLog" ){
          options_before_imp <- c("name", "id", "delim",
                                  "intensitycols", "filt",
                                  "imputation", "imp_norm_order", "norm_method")
          # input11 <<- input
          check_opt_change <- options_before_imp %>% sapply(.,function(x){
            identical(input[[x]], upload_log()$inputVals()[[x]])
          })

          # check_opt_change2 <<- check_opt_change

          if(all(check_opt_change)){ # if the option before imputation unchange, directly use the old imputation result
            imp_norm = upload_log()$resultVals()$pe_norm
            message("Use the imputation in log file!")
            return(imp_norm)
          }

        }

        pe_filt = pe_filt()
        # cat("begin normalization ...")
        if(input$imp_norm_order){
          withProgress(message = "Imputate & normalize", value = 0.7,{
            cat("Impute before normalize ...\n")
            pe_norm <- addAssay(pe_filt,
                                DEP2::impute(pe_filt[["peptideRaw"]], fun = input$imputation),
                                name = "peptideImp")
            # pe_norm <- QFeatures::impute(pe_filt, method = input$imputation, i = "peptideRaw", "peptideImp")
            message("\tImpute finished")
            pe_norm <- normalize_pe(pe_norm, method = input$norm_method, i = "peptideImp", name = "peptideNorm")
            message("\tNormalize finished")
          })
        }else{
          withProgress(message = "Normalize & imputate", value = 0.7,{
            cat("Impute after normalize ...\n")
            pe_norm <- normalize_pe(pe_filt, method = input$norm_method, i = "peptideRaw", name = "peptideNorm")
            message("\tNormalize finished")

            pe_norm <- addAssay(pe_norm,
                                DEP2::impute(pe_norm[["peptideNorm"]], fun = input$imputation),
                                name = "peptideImp")
            # pe_norm <-  QFeatures::impute(pe_norm, method = input$imputation, i = "peptideNorm")
            message("\tImpute finished")
          })
        }

        return(pe_norm)
        # }else{
        #   load(rdataFile$datapath)
        #   return(pe_norm)
        # }
      })


      data_QF <- reactiveVal(NULL) ## data store the aggregated QF object.


      filt0 <- reactive({
        # rdataFile <- input$resultRData
        # if(is.null(rdataFile)){
        data <- data_QF()[["protein"]]
        # }else{
        #   data_QF <- data_QF({
        #     rdataFile <- input$resultRData
        #     load(rdataFile$datapath)
        #     cat("load rdataFile in filt0 \n")
        #     data
        #   })
        #   data <- data_QF()[["protein"]]
        # }

        # req(input$id, input$name)
        validate(need(!(input$name == "" & input$id == ""), "Please ensure that: at least one of your name column and id column is non-empty!"))
        if(input$name == "" & input$id == "") {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = "Please ensure that: at least one of your name column and id column is non-empty!",
            type = "warning"
          )
        }

        # unique_names <- make_unique(proteins = filtered, names = ifelse(input$name == "", input$id, input$name), ids = ifelse(input$id == "", input$name, input$id), delim = input$delim)
        rowData(data) <- DEP2::make_unique(rowData(data) %>% as.data.frame(),
                                           names =  ifelse(input$name == "",
                                                           "smallestProteingroups",
                                                           # ifelse(input$aggregate_Peptide_Type == "Unique + Razor", "smallestProteingroups", input$id),
                                                           input$name),
                                           ids = "smallestProteingroups",
                                           # ifelse(input$aggregate_Peptide_Type == "Unique + Razor", "smallestProteingroups", input$id),
                                           delim = input$delim)
        rownames(data) <- rowData(data)$name
        cat("Make_unique finished")
        data
      })

      iv1 <- InputValidator$new()
      iv1$add_rule("order",
                   sv_required(message = "Complete order is required, sort all groups", test = function(x){
                     is.null(x) || length(x) == (length(filt0()@colData$condition %>% unique))
                   })
      )
      iv1$enable()

      the_order <- reactive({
        req(iv1$is_valid())

        if(length(input$order) == length(filt0()@colData$condition %>% unique))
          return(input$order)
        return(NULL)
      })

      filt <- reactive({
        filt = filt0()

        if(is.null(the_order())){
          return(filt)
        }else{

          filt@colData$condition = factor(filt@colData$condition, levels = the_order())
          filt = filt[,order(filt@colData$condition)]
          return(filt)
        }
      })

      norm <- reactive({
        norm <- normalize_vsn(filt())
        cat("VSN finished\n")
        return(norm)
      })


      imp <- reactive({
        norm()
      })



      FDR_type <- reactive({
        FDR_type <- input$FDR_type
        # FDR_type <- ifelse(length(grep(FDR_type,"qval"))>1,"qval","lfdr")
      })

      df <- reactive({
        cat("begin test_diff\n")
        if(input$contrasts == "control") {
          validate(
            need(input$control != "", "Please select a control condition under menuItem Columns in the DEP-LFQ options of the sidebar")
            # need(all(input$intensitycols %in% colnames(rowData(data_QF()[["protein"]]))), "Please select the Expression columns in the sidebar")
          )
        }

        if(input$contrasts == "manual") {
          validate(
            need(input$test_manual != "", "Please select manual contrasts to test under menuItem Columns")
            # need(all(input$intensitycols %in% colnames(rowData(data_QF()[["protein"]]))), "Please select the Expression columns in the sidebar")
          )
        }

        inFile1 <- NULL
        if(is.null(inFile1)){
          if(input$contrasts == "control"){
            df <- DEP2::test_diff(se = imp(), type = input$contrasts, control = input$control,
                                  fdr.type = FDR_type())
          }

          if(input$contrasts == "all") {
            df <- DEP2::test_diff(se = imp(), type = "all",
                                  fdr.type = FDR_type())
          }

          if(input$contrasts == "manual") {
            df <- DEP2::test_diff(se = imp(), type = input$contrasts, test = input$test_manual,
                                  fdr.type = FDR_type())
          }
        } else {
          #load your saved RData in order to get the same result (imp is the key of if result from two analysis being the same)
          load(file = inFile1$datapath)
          if(input$contrasts == "control"){
            df <- DEP2::test_diff(se = my_imp, type = input$contrasts, control = input$control,
                                  fdr.type = FDR_type())
          }

          if(input$contrasts == "all") {
            df <- DEP2::test_diff(se = my_imp, type = "all",
                                  fdr.type = FDR_type())
          }

          if(input$contrasts == "manual") {
            df <- DEP2::test_diff(se = my_imp, type = input$contrasts, test = input$test_manual,
                                  fdr.type = FDR_type())
          }
        }
        df
      })

      dep <- reactive({
        if(input$threshold_method == "intersect") {
          my_dep <- add_rejections(diff = df(), alpha = input$p, lfc = input$lfc, thresholdmethod = input$threshold_method)
        }
        if(input$threshold_method == "curve") {
          my_dep <- add_rejections(diff = df(), thresholdmethod = input$threshold_method, curvature = input$curvature, x0_fold = input$x0_fold)
        }

        if(is.null(input$peptidescol)|input$peptidescol==""){
          rowData(my_dep)$Peptides = 3
        }else{
          rowData(my_dep)$Peptides = rowData(my_dep)[,input$peptidescol]
        }
        return(my_dep)
      })


      ## All object and functions upon 'Aggregation' button input ----------------------
      observeEvent(input$aggregate_button,{

        print("aggregation begin")
        # req(input$id, input$name, input$intensitycols, input$imputation, input$aggregate_Peptide_Type)
        if(is.null(input$intensitycols)) {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = "Please ensure the expression columns in 'Data proccess' opt",
            type = "warning"
          )
        }
        req(input$intensitycols)


        if(is.null(input$imputation)) {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = "Please check Imputation type in the 'Data proccess' opt",
            type = "warning"
          )
        }
        req(input$imputation)

        if(is.null(input$aggregate_Peptide_Type)) {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = "Please check options in the 'Aggregation' opt",
            type = "warning"
          )
        }
        req(input$aggregate_Peptide_Type)

        # if(is.null(input$imputation)) {
        #   sendSweetAlert(
        #     session = shiny::getDefaultReactiveDomain(),
        #     title = "warning !",
        #     text = "Please ensure that: select expression columns in the Columns opt",
        #     type = "warning"
        #   )
        # }
        # req(input$imputation)

        if(is.null(input$filt)) {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "Note :",
            text = "Your filter column is empty! If needed, please choose the filter columns. Otherwise, please ignore this message",
            type = "info"
          )
        }

        data_QF <- data_QF({
          input$aggregate_button
          isolate({
            pe_norm = pe_norm()
            fcol = input$id

            cat("Begin aggregate... \n")
            aggrefun = switch(input$aggregation_method,
                              totalSum = base::colSums,
                              totalMean = base::colMeans,
                              medianPolish = MsCoreUtils::medianPolish,
                              RobustSummary = MsCoreUtils::robustSummary)
            i = ifelse(input$imp_norm_order, "peptideNorm","peptideImp")

            withProgress(message = 'Aggregeting, please wait for a while', {
              protein <- DEP2::aggregate_pe(pe_norm, aggrefun = aggrefun, aggregate_Peptide_Type = input$aggregate_Peptide_Type,
                                            fcol = fcol, i = i, reserve = input$name)
            })
            cat("Finished \n")
            protein
          })
        })

        output$save_aggregate_button <- renderUI({
          downloadButton(outputId = session$ns("save_aggregate_button2"), label = "save aggregation result")
        })

      })

      output$save_aggregate_button2 <- downloadHandler(
        filename = function() {
          paste0("aggregation_result.RData")
        },
        content = function(file) {
          peptide = peptide()
          pe = pe()
          pe_norm = pe_norm()
          data = data_QF()
          save(peptide, pe, pe_norm, data,file = file)
        }
      )

      output$Save_RData <- downloadHandler(
        filename = function() { paste("results", ".RData", sep = "") },
        content = function(file) {
          withProgress(message = 'Please wait ...', value = 0.66, {
            peptide = peptide()
            pe = pe()
            pe_norm = pe_norm()
            data = data_QF()
            filt = filt(); norm = norm();imp = imp(); dep = dep()
            save(peptide, pe, pe_norm, data, filt = filt, norm, imp, dep, file=file)})}
      )

      output$Save_log <- downloadHandler(
        filename = function() { paste("Export_log_pep", format(Sys.time(), "%m%d_%H%M"),
                                      ".logrds", sep = "") },
        content = function(file) {
          withProgress(message = 'Compressing result. Please wait ...', value = 0.66, {
            message("Export log...")
            thelog = pgpep_log$new(input = input,
                                 reactive_vals = list(
                                   peptide = peptide(), # required
                                   expdesign = expdesign(),  # required
                                   pe = pe(),
                                   pe_filt = pe_filt(),
                                   pe_norm = pe_norm(),
                                   data_QF = data_QF(),
                                   filt0 = filt0(),
                                   filt = filt(),
                                   norm = norm(),
                                   imp = imp(),
                                   dep = dep()
                                 ),
                                 app_version = DEP2:::app_version
            )
            thelog111 <<- thelog
            saveRDS(thelog, file = file)
          })
        }
      )

      ## All object and functions upon 'Analyze' button input ----------------------
      observeEvent(input$analyze, {

        ### Interactive UI functions ### ------------------------------------------
        output$downloadTable <- renderUI({
          selectizeInput(ns("dataset"),
                         "Choose a dataset to save" ,
                         c("results","significant_proteins",
                           "displayed_subset","full_dataset"))
        })

        output$downloadButton <- renderUI({
          downloadButton(ns('downloadData'), 'Save', class = "downloadData")
        })

        output$downloadButton_for_save_RData <- renderUI({
          downloadButton(ns("Save_RData"), "save result RData", class = "Save_RData")
        })

        output$downloadButton_for_log <- renderUI({ ## download log file
          downloadButton(session$ns("Save_log"), "Export log file")
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
          selectizeInput(ns("select"),
                         "Select direct comparisons",
                         choices=names,
                         multiple = TRUE)
        })

        output$exclude <- renderUI({
          row_data <- rowData(dep())
          cols <- grep("_significant", colnames(row_data))
          names <- colnames(row_data)[cols]
          names <- gsub("_significant","",names)
          selectizeInput(ns("exclude"),
                         "Exclude direct comparisons",
                         choices = names,
                         multiple = TRUE)
        })


        output$heatmap_cntrst <- renderUI({
          if (!is.null(selected())) {
            df <- rowData(selected())
            cols <- grep("_significant$",colnames(df))
            selectizeInput(ns("heatmap_cntrst"),
                           "Contrast",
                           choices = gsub("_significant", "", colnames(df)[cols]), multiple = TRUE)
          }
        })


        output$mysplit <- renderUI({
          selectizeInput(ns("mysplit"),
                         "my split",
                         choices = c(1 : input$k), multiple = TRUE)

        })

        output$Custom_columns_order <- renderUI({
          selectizeInput(ns("Custom_columns_order"),
                         "Custom columns order",
                         choices = colnames(assay(dep())), multiple = TRUE)

        })

        output$Custom_columns_order_for_missval_heatmap <- renderUI({
          selectizeInput(ns("Custom_columns_order_for_missval_heatmap"),
                         "Custom columns order",
                         choices = colnames(assay(dep())), multiple = TRUE)
        })


        output$volcano_cntrst <- renderUI({
          if (!is.null(selected())) {
            df <- rowData(selected())
            cols <- grep("_significant$",colnames(df))
            selectizeInput(ns("volcano_cntrst"),
                           "Contrast",
                           choices = gsub("_significant", "", colnames(df)[cols]))
          }
        })

        output$Volcano_cntrst <- renderUI({
          if (!is.null(selected())) {
            df <- rowData(selected())
            cols <- grep("_significant$",colnames(df))
            selectizeInput(ns("Volcano_cntrst"),
                           "Contrast",
                           choices = gsub("_significant", "", colnames(df)[cols]))
          }
        })

        output$Volcano_cntrst_stastical_plot <- renderUI({
          if (!is.null(selected())) {
            df <- rowData(selected())
            cols <- grep("_significant$",colnames(df))
            selectizeInput(ns("Volcano_cntrst_stastical_plot"),
                           "Contrast",
                           choices = gsub("_significant", "", colnames(df)[cols]))
          }
        })

        output$cntrst_norm_distribution <- renderUI({
          if (!is.null(selected())) {
            df <- rowData(selected())
            cols <- grep("_significant$",colnames(df))
            selectizeInput(ns("cntrst_norm_distribution"),
                           "Contrast",
                           choices = gsub("_significant", "", colnames(df)[cols]))
          }
        })

        output$selected_proteins <- renderUI({
          row_Data <- rowData(dep())
          selectizeInput(ns("selected_proteins"),
                         "selected proteins",
                         choices = row_Data$name, selected = NULL, multiple = TRUE)
        })

        output$chooseToshow <- renderUI({
          selectizeInput(ns("chooseToshow"),
                         "selected proteins",
                         choices = name_for_choose(), selected = NULL, multiple = TRUE)
        })

        output$plot_stat_contrasts <- renderUI({
          selectInput(session$ns("plot_stat_contrasts"),
                      "plot contrasts",
                      choices = get_contrast(selected()),
                      multiple = T,
                      selected = get_contrast(selected()))
        })


        output$umap_n_neighbors <- renderUI(
          sliderInput(inputId = ns("umap_n_neighbors"),"number of nearest neighbors",min=2,max=min(15,ncol(dep())),value=4)
        )

        output$pca_top_n  <- renderUI(
          sliderInput(inputId = ns("pca_top_n"),"number of top variable proteins to consider",min = 2, max = nrow(dep()), value=500)
        )

        output$Tsne_perplexity <- renderUI(
          sliderInput(inputId = ns("Tsne_perplexity"),"perplexity",min=1,max= (ncol(dep())-1)/3 - 2^(-31) ,value=2)
        )

        ### Reactive functions for plots ------------------------------------------------
        name_for_choose <- reactive({
          # heatmap_name(dep = selected(), manual = input$manual_heatmap, manual_name = input$heatmap_cntrst)
          DEP2::get_signicant(selected(), contrasts = input$heatmap_cntrst, return_type = "names")
        })

        excluded <- reactive({
          DEP2:::exclude_contrasts(dep(), input$exclude)
        })

        selected <- reactive({
          DEP2:::select_contrasts(excluded(), input$select)
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
            if(input$cluster_columns || is.null(input$Custom_columns_order) || length(input$Custom_columns_order) < 2) {
              selected = selected()
            } else {
              cols <- input$Custom_columns_order
              selected <- selected()[,cols]
            }
            DEP2::plot_heatmap(object = selected,
                               type = input$pres,
                               manual_contrast = input$heatmap_cntrst,
                               kmeans = TRUE,
                               k = input$k,
                               color = input$colorbar,
                               col_limit = input$limit,
                               row_font_size = input$row_font_size,
                               col_font_size = input$col_font_size,
                               cluster_columns = input$cluster_columns,
                               split_order = input$mysplit,
                               label_few_peptide_rows = input$if_rowname_color,
                               chooseToshow = input$chooseToshow
            )
            # plot_heatmap(dep = selected(),
            #              type = input$pres,
            #              manual = input$manual_heatmap,
            #              manual_name = input$heatmap_cntrst,
            #              kmeans = TRUE,
            #              k = input$k,
            #              color = input$colorbar,
            #              col_limit = input$limit,
            #              #same_trend = input$same_trend,
            #              row_font_size = input$row_font_size,
            #              col_font_size = input$col_font_size,
            #              cluster_columns = input$cluster_columns,
            #              if_mysplit = input$if_mysplit,
            #              mysplit = input$mysplit,
            #              if_rowname_color = input$if_rowname_color,
            #              if_chooseToshow = input$if_chooseToshow,
            #              chooseToshow = input$chooseToshow,
            #              column_order = if(input$cluster_columns) {NULL} else {input$Custom_columns_order})
          })
        })

        # volcano_input <- reactive({
        #   if(!is.null(input$volcano_cntrst)) {
        #     plot_volcano(selected(),
        #                  input$volcano_cntrst,
        #                  input$fontsize,
        #                  input$check_names,
        #                  input$p_adj,
        #                  same_width = input$same_width,
        #                  my_breaks = input$my_breaks,
        #                  mybreaks = as.numeric(input$mybreaks))
        #   }
        # })

        custom_volcano_input <- reactive({
          if(!is.null(input$Volcano_cntrst)) {
            withProgress(message = 'Plotting', value = 0.66, {
              DEP2::plot_volcano(object = selected(),
                                 contrast = input$Volcano_cntrst,
                                 add_threshold_line = input$threshold_method,
                                 adjusted = ifelse(input$threshold_method == "intersect", input$P_adj, TRUE),
                                 label_size = input$fontSize,
                                 label_number = input$showNum,
                                 label_trend = input$labelWay,
                                 chooseTolabel = input$selected_proteins,
                                 x_symmetry = input$Same_width,
                                 fcCutoff = input$lfc,
                                 pCutoff = input$p,
                                 curvature = input$curvature,
                                 x0_fold = input$x0_fold,
                                 dot_size = input$dotsize,
                                 highlight_PGs_with_few_peptides = input$if_peptide_color,
                                 # label.rectangle = input$if_label_rectangle,
                                 # stroke = input$stroke,
                                 down_color = input$down_color,
                                 stable_color = input$stable_color,
                                 up_color = input$up_color,
                                 peptide_1_color = input$Peptides_1,
                                 peptide_2_color = input$Peptides_2,
                                 # threshold_method = input$threshold_method,
                                 # PTM = FALSE
              )
            })
          }
        })

        stastical_plot_input <- reactive({
          if(!is.null(input$plot_stat_contrasts)){
            # statistics_plot(dep = selected(),
            #                 x = input$stastical_plot_x, y = input$stastical_plot_y,
            #                 plottype = input$stastical_plot_type,
            #                 contrast=input$Volcano_cntrst_stastical_plot)
            DEP2::plot_statistics(selected(), statistic = input$stastical_plot_x ,contrasts = input$plot_stat_contrasts)
          }
        })

        norm_distribution_input <- reactive({
          if(!is.null(input$cntrst_norm_distribution)){
            DEP2::plot_norm_distribution(dep = dep(),
                                         contrast=input$cntrst_norm_distribution)
          }
        })

        pca_input <- reactive({
          DEP2::plot_pca(object = dep(), indicate = input$Indicate, if_square = input$if_square, n = input$pca_top_n, features = "proteins")
        })

        umap_input <- reactive({
          DEP2::plot_umap(dep(), indicate = input$umap_Indicate, if_square = input$umap_if_square, n_neighbors = input$umap_n_neighbors, features = "proteins")
        })

        Tsne_input <- reactive({
          DEP2::plot_Tsne(dep(), indicate = input$Tsne_Indicate, if_square = input$Tsne_if_square, perplexity = input$Tsne_perplexity, theseed = input$Tsne_theseed, features = "proteins")
        })

        Pearson_correlation_input <- reactive({
          DEP2::plot_cor(dep(), pal = input$Pearson_pal, pal_rev = input$Pearson_pal_rev, lower = input$Pearson_lower, upper = input$Pearson_upper, add_values = input$add_values_for_DEP_person, value_size = input$value_size_for_DEP_person, digits = input$value_digits_for_DEP_person)
        })

        meanSdPlot_input <- reactive({
          DEP2::meanSdPlot(dep())
        })

        Gowers_distance_input <- reactive({
          DEP2::plot_dist(dep(), pal = input$Gower_pal, pal_rev = input$Gower_pal_rev, add_values = input$add_values_for_DEP_gower, value_size = input$value_size_for_DEP_gower, digits = input$value_digits_for_DEP_gower)
        })

        Sample_CVs_input <- reactive({
          DEP2::plot_cvs(dep())
        })


        norm_input <- reactive({
          protein <- data_QF()
          thenorm = which(names(protein) == "peptideNorm")
          Peptide_before = protein[[thenorm-1]]
          Peptide_normalized = protein[[thenorm-1]]
          DEP2:: plot_normalization(Peptide_before,
                                    Peptide_normalized)
        })

        # missval_input <- reactive({
        #   withProgress(message = 'Plotting missing heatmap', value = 0.3, {
        #     plot_missval(pe_filt()[["peptideRaw"]],
        #                  # norm(),
        #                  cluster_columns = input$cluster_columns_for_missval_heatmap,
        #                  column_order = if(input$cluster_columns_for_missval_heatmap) {NULL} else {input$Custom_columns_order_for_missval_heatmap}
        #     )
        #   })
        # })

        detect_input <- reactive({
          DEP2::plot_detect(data_QF()[["peptideRaw"]])
        })

        imputation_input <- reactive({
          peptide_raw <- data_QF()[["peptideRaw"]]
          peptide_norm <- data_QF()[["peptideNorm"]]
          peptide_imp <- data_QF()[["peptideImp"]]

          if(input$imp_norm_order){
            DEP2::plot_imputation(peptide_raw,
                                  peptide_imp,
                                  peptide_norm)
          }else{
            DEP2::plot_imputation(peptide_raw,
                                  peptide_norm,
                                  peptide_imp)
          }
        })

        numbers_input <- reactive({
          DEP2::plot_numbers(se = norm(), plot = TRUE, features_type = "proteins")
        })

        coverage_input <- reactive({
          DEP2::plot_coverage(se = norm(), plot = TRUE, features_type = "proteins")
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

        # output$volcano <- renderPlot({
        #   volcano_input()
        # })

        output$custom_volcano <- renderPlot({
          custom_volcano_input()
        })

        output$stastical_plot <- renderPlot({
          stastical_plot_input()
        })

        output$norm_distribution_plot <- renderPlot({
          norm_distribution_input()
        })

        output$pca <- renderPlot({
          pca_input()
        })

        output$umap <- renderPlot({
          umap_input()
        })

        output$Tsne <- renderPlot({
          Tsne_input()
        })

        output$Pearson_correlation <- renderPlot({
          Pearson_correlation_input()
        })

        output$meanSdPlot <- renderPlot({
          meanSdPlot_input()
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

        # output$missval <- renderPlot({
        #   missval_input()
        # })

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
            plotOutput(ns("heatmap"), height = (100 * as.numeric(input$size)))
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

        output$downloadstastical_plot <- downloadHandler(
          filename = function() {
            paste0("statistical_plot_", input$Volcano_cntrst_stastical_plot, ".pdf")
          },
          content = function(file) {
            pdf(file, width = input$stastical_plot_Width, height = input$stastical_plot_Height)
            print(stastical_plot_input())
            dev.off()
          }
        )

        output$download_norm_distribution_plot <- downloadHandler(
          filename = 'norm_distribution.pdf',
          content = function(file) {
            pdf(file, width = input$norm_distribution_Width, height = input$norm_distribution_Height)
            print(norm_distribution_input())
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

        output$downloadUMAP <- downloadHandler(
          filename = 'UMAP.pdf',
          content = function(file) {
            pdf(file, width = input$umap_Width, height = input$umap_Height)
            print(umap_input())
            dev.off()
          }
        )

        output$downloadTSNE <- downloadHandler(
          filename = 'TSNE.pdf',
          content = function(file) {
            pdf(file, width = input$Tsne_Width, height = input$Tsne_Height)
            print(Tsne_input())
            dev.off()
          }
        )

        output$downloadmeanSdPlot <- downloadHandler(
          filename = 'meanSdPlot.pdf',
          content = function(file) {
            pdf(file, width = input$meanSdPlot_Width, height = input$meanSdPlot_Height)
            print(meanSdPlot_input())
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

        # output$downloadMissval <- downloadHandler(
        #   filename = "missing_values_heatmap.pdf",
        #   content = function(file) {
        #     pdf(file, width = input$missval_heatmap_Width, height = input$missval_heatmap_Height)
        #     print(missval_input())
        #     dev.off()
        #   }
        # )

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
      }) ## end of analyze event

      ###"question for ptm"
      observeEvent(input$help_format_pep, {
        showModal(modalDialog(
          title = "Format specifications",
          # includeMarkdown(system.file("extdata", "DEP_LFQ.md", package = "DEP2")),
          includeMarkdown("www/DEP_pep.md"),
          br(),
          easyClose = TRUE,
          footer = NULL,
          size = "l"
        ))
      })

      return(df)

    })
}
