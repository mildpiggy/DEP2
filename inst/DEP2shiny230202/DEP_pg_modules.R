## the sidebar_mod body_mod and server mod of DEP analysis from proteingroups quantification

get_DEPglobal_env <- function () {
  if (!exists(".DEPglobal_data_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".DEPglobal_data_Env", new.env(), envir=envir)
  }
  get(".DEPglobal_data_Env", envir = .GlobalEnv)
}

DEP_pg_UI <- function(id, labelname){
  tagList(
    sidebarLayout(
      DEP_pg_sidebar_mod(id = id,labelname=labelname), ## DEP sidebar module
      DEP_pg_body_mod(id = id)  ## DEP pagebody module
    ))
}

DEP_pg_sidebar_mod <-  function(id,labelname){
  ns = NS(id)
  cat(paste0("DEP_pg_sidebar ns is",ns(""),"\n"))

  tagList(
    # sidebarLayout(
      sidebarPanel(
        width = 3,
        h3(labelname),
        bsCollapse(
          id = ns("DEPpg_sidebar"),
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

                          # file input options, from Quantity table
                          conditionalPanel(
                            condition = paste0("input['",ns("uploadmode"),"']","== 'fromTable'"), ## condition in mod
                            fluidRow(
                              column(width = 9,
                                     fileInput(ns('file1'),width = "300px",
                                               'ProteinGroups.txt',
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
                                                          '.csv','.gz',
                                                          '.tsv'))
                                     )
                                     # ,
                                     # fileInput(ns('resultRData'),'Load saved result RData',width = "300px")
                              ),
                              column(
                                width = 1, h4(),
                                br(),
                                actionButton(ns("help_format_DEP"),
                                             label = "", icon = icon("question-circle"),
                                             style = "color: #f6f6f6; background-color: #2c3b41; border-color: #2c3b41"
                                ),
                                shinyBS::bsTooltip(
                                  ns("help_format_DEP"),
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

                          ## how to generate experiment design
                          radioButtons(ns("anno"),
                                       "Sample annotation",
                                       choices = list("Parse from columns" = "columns",
                                                      "Use Experimental Design" = "expdesign"),
                                       selected = "columns")

          ),
          # menuItemOutput("columns"),
          bsCollapsePanel("Columns",
                          style = "primary",
                          uiOutput(outputId = ns("name")),
                          uiOutput(outputId = ns("id")),
                          selectizeInput(inputId = ns("delim"), "Delimiter", choices = c(";", "|"), selected = ";"),
                          uiOutput(ns("intensitycols")),
                          uiOutput(ns("peptidescol")),
                          # menuItem("Expression columns", icon = icon("th"),
                          #   uiOutput("intensitycols")),
                          uiOutput(ns("filt")),
                          numericInput(ns("thr"), "Allowed max.miss.num at least one condition", min = 0, max = 20, value = 0),
                          uiOutput(ns("order")),
                          uiOutput(ns("control")),
                          uiOutput(ns("test_manual")),
                          shinyBS::bsTooltip(ns("order"), "Set order of groups", "top", options = list(container = "body"))
          ),
          bsCollapsePanel("Imputation options",
                          style = "primary",
                          radioButtons(inputId = ns("imputation"),
                                       label = "Imputation type",
                                       choices = c(c("man", "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "min", "zero"), "mixed on proteins",
                                                   # "mixed on samples",  # 'mixed on samples' impute control and other condition with different method
                                                   "RF"),  ## increase "RF" opt 20211023
                                       selected = "MinProb"),
                          shinyBS::bsTooltip(ns("imputation"), "Choose an imputation method", "right", options = list(container = "body")),
                          actionButton(ns("help_imputation"),
                                       label = "Detailed information", #icon = icon("question-circle"),
                                       style = "color: #f6f6f6; background-color: #2c3b41; border-color: #2c3b41; text-align: left; margin: 0 0 0px 0px",
                                       #rgba(255, 255, 255, 0.1) rgba(0, 0, 0, 0.1)
                          ),
                          shinyBS::bsTooltip(
                            ns("help_imputation"),
                            "The detailed information of imputation methods",
                            "right",
                            options = list(container = "body")
                          )
                          # p(a("Detailed information link",
                          #     href = "https://www.rdocumentation.org/packages/MSnbase/versions/1.20.7/topics/impute-methods",
                          #     target="_blank"))
          ),
          bsCollapsePanel("FDR correction",  ## DEP FDR correct options
                          style = "primary",
                          radioButtons(ns("FDR_type"),
                                       "type of FDR value",
                                       choices = c("Strimmer's qvalue(t)","Strimmer's qvalue(p)","BH","Storey's qvalue")),
                          shinyBS::bsTooltip("fdr_correction", "Choose the method of pvalue adjustment", "right", options = list(container = "body"))
          ),
          bsCollapsePanel("Threshold method",
                          style = "primary",
                          radioButtons(ns("threshold_method"),
                                       "Threshold method",
                                       choices = c("intersect", "curve"),
                                       selected = "intersect"),
                          shinyBS::bsTooltip(ns("threshold_method"), "Choose the method that the cutoff of significant proteins based on. [intersect] means by adjusted pvalue and log2 fold change cutoff and [curve] means by curvature and x0 fold cutoff", "right", options = list(container = "body"))
          ),
          open = "Files"
        ),

        shinyBS::bsTooltip(ns("anno"), "Sample annotation type", "top", options = list(container = "body")),
        shinyBS::bsTooltip(ns("name"), "Name of the column containing feature names", "right", options = list(container = "body")),
        shinyBS::bsTooltip(ns("id"), "Name of the column containing feature IDs", "right", options = list(container = "body")),
        shinyBS::bsTooltip(ns("delim"), "Set the delimiter separating the name and id column", "right", options = list(container = "body")),
        shinyBS::bsTooltip(ns("intensitycols"), "Choose the expression columns of your data", "right", options = list(container = "body")),
        shinyBS::bsTooltip(ns("remove_prefix"), "remove the prefix character of expression columns names, such as 'LFQ intensity' or 'intensity' in Maxquant result",
                           "right", options = list(container = "body")),
        shinyBS::bsTooltip(ns("remove_suffix"), "remove the suffix character of expression columns names",
                           "right", options = list(container = "body")),
        shinyBS::bsTooltip(ns("peptidescol"), "Choose the peptides column of your data. And this column should be numeric. It is only used to plot heatmap and custom volcano to color the peptides number 1 and 2 to help user check the accuracy of the quantification. Of course, if your data does not have this column, You can leave it blank", "right", options = list(container = "body")),
        shinyBS::bsTooltip(ns("filt"), "The filtered columns based on", "right", options = list(container = "body")),
        shinyBS::bsTooltip(ns("thr"), "Set the threshold for the allowed number of missing values in at least one condition", "right", options = list(container = "body")),
        shinyBS::bsTooltip(ns("control"), "Choose your control condition", "right", options = list(container = "body")),
        shinyBS::bsTooltip(ns("test_manual"), "Choose the contrasts that you want to test", "right", options = list(container = "body")),

        actionButton(ns("analyze"), "Analyze"),

        shinyBS::bsTooltip("analyze", "Click on it to analyze your data", "right", options = list(container = "body")),
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
        shinyBS::bsTooltip(ns("downloadButton_for_save_RData"), "Click on it to Save the RData, and when you upload the results.RData in the menuItem [Files] next time, you can get the same result as this time", "right", options = list(container = "body")),

        # ),
        # mainPanel("bbb")
      )
  )
}

DEP_pg_body_mod <- function(id){
  ns = NS(id)
  cat(paste0("DEP_pg_body_mod ns is",ns(""),"\n"))
  # cat(ls(envir = parent.frame()))
  # # cat(parent.frame()$modid)
  # modid = get("modid",envir = parent.frame())
  # assign(modid(id),"modid",envir = parent.frame())
  # parent.frame()$modid = modid(id)
  tagList(
    mainPanel(#** DEP tabItem----
              tabName = id,
              width = 9,
              helpText("Please cite: "),
              fluidRow(

                conditionalPanel(condition = paste0("input['",ns("threshold_method"),"']","== 'intersect'"), ## condition in mod
                                 box(numericInput(ns("p"),
                                                  "adj. P value",
                                                  min = 0.0001, max = 0.1, value = 0.05),
                                     width = 2),
                                 box(numericInput(ns("lfc"),
                                                  "Log2 fold change",
                                                  min = 0, max = 10, value = 1),
                                     width = 2)
                ),
                conditionalPanel(condition = paste0("input['",ns("threshold_method"),"']","== 'curve'"),
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
                                 selected = "centered", inline = TRUE),
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
                       tabBox(
                         title = "Result Plots", width = 12, id = ns("DEP_results_tabs"),

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
                                    # conditionalPanel(condition = paste0("!input['",ns("if_chooseToshow"),"']"),#&& "!input.if_chooseToshow" input.peptidescol !== '' typeof input.peptidescol !== null
                                    #                  box(numericInput(ns("k"),
                                    #                                   "Kmeans",
                                    #                                   min = 0, max = 15, value = 6),
                                    #                      width = 4)
                                    # ),
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

                                      width = "300px"),
                                    br("")
                                  ),
                                  fluidRow(
                                    uiOutput(ns("plot")),
                                    downloadButton(ns('downloadHeatmap'), 'Save heatmap')
                                  ),

                                  shinyBS::bsTooltip(ns("colorbar"), "Choose the color palette of the colorbar", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("heatmap_cntrst"), "Choose the contrasts that you want to show, act when [Manual heatmap] is TRUE, can be one or more contrasts", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("k"), "Set the number of k-means clusters", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("limit"), "Set the outer limits of the color scale", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("size"), "Set the height of the plot region", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("if_rowname_color"), "Whether set the color of the row names,(purple: Peptides == 1, blue: Peptides == 2)", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("chooseToshow"), "Choose the rows that you want to show", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("row_font_size"), "Set the size of row labels", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("col_font_size"), "Set the size of column labels", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("Height"), "Height of the figure to export", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("cluster_columns"), "Whether make cluster on columns", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("Custom_columns_order"), "Order of column, act when [Cluster columns] is FALSE", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("if_mysplit"), "Whether order the splices manually", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("mysplit"), "The ordered splices of the heatmap (from top to bottom), only act when both [If mysplit] and [Kmeans] are TRUE", "top", options = list(container = "body"))
                         ),

                         tabPanel(title = "Custom Volcano",
                                  fluidRow(
                                    column(width = 4,
                                           colourpicker::colourInput(inputId = ns("up_color"),
                                                                     label = "up",
                                                                     showColour = "both",
                                                                     palette = "square",
                                                                     value = "#B2182B")
                                    ),
                                    column(width = 4,
                                           colourpicker::colourInput(inputId = ns("down_color"),
                                                                     label = "down",
                                                                     showColour = "both",
                                                                     palette = "square",
                                                                     value = "#2166AC")
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
                                                       choices = c("all", "up","down"),
                                                       selected = c("all"), multiple = FALSE),width = 4)
                                  ),
                                  fluidRow(
                                    dropdownButton(
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
                                  # shinyBS::bsTooltip(ns("if_label_rectangle"), "Whether add rectangle underneath the text, making it easier to read", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("if_peptide_color"), "Whether set the color of the points,(default, purple: Peptides == 1, blue: Peptides == 2)", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("Peptides1"), "Set the color of the points which the number of Peptides == 1", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("Peptides2"), "Set the color of the points which the number of Peptides == 2", "top", options = list(container = "body")),
                                  # shinyBS::bsTooltip(ns("stroke"), "Set the thickness of black line around the point", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("selected_proteins"), "Choose the point labels to show, act when [label Way] is selected proteins", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("showNum"), "Set the number of the labels to add, act when [label way] is significant, up or down. Note that, when it is larger than the number of significant result, it is seted to the number of significant result", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("fontSize"), "Set the font size of the label", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("dotsize"), "Set the dot size of the label", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("custom_volcano_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("custom_volcano_Height"), "Height of the figure to export", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("P_adj"), "Whether or not to use adjusted p values", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("Same_width"), "Whether the x axis to have the same width from 0", "top", options = list(container = "body"))
                         )
                       ),
                       hr(),

                       ## QC tabs ----
                       tabsetPanel(id = ns("DEP_QC_tabs"),
                       # tabBox(title = "QC Plots", width = 12, id = ns("DEP_QC_tabs"),
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
                         #tabPanel(title = "Pca plot",
                         #         plotOutput("pca", height = 600),
                         #         downloadButton('downloadPca', 'Save')
                         #)
                         tabPanel(title = "Statistical plot",
                                  # fluidRow(
                                  #   box(uiOutput(ns("Volcano_cntrst_stastical_plot")), width = 8),
                                  #   box(radioButtons(ns("stastical_plot_type"),
                                  #                    "plot type",
                                  #                    choices = c("x-y","histogram")), width = 4)
                                  # ),
                                  fluidRow(
                                    box(selectInput(ns("stastical_plot_x"),
                                                    "plot statistics",
                                                    choices = c("diff","t.stastic","p.val","p.adj"),
                                                    multiple = T,
                                                    selected = c("diff","t.stastic","p.val","p.adj")),
                                        width = 6),
                                    # box(selectInput(ns("stastical_plot_y"),
                                    #                 "y",
                                    #                 choices = c("-log10.padj", "-log10.pval","padj","pval","t")), width = 6)
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
                                    downloadButton(ns('downloadmeanSdPlot'), 'Save dataset')
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
                         tabPanel(title = "Missing values - Heatmap",
                                  fluidRow(
                                    box(
                                      checkboxInput(ns("cluster_columns_for_missval_heatmap"),
                                                    "Cluster columns",
                                                    value = TRUE), width = 4),
                                    box(uiOutput(ns("Custom_columns_order_for_missval_heatmap")),
                                        width = 8)
                                  ),
                                  fluidRow(
                                    box(numericInput(ns("missval_heatmap_Width"),
                                                     "width",
                                                     min = 1, max = 30, value = 7), width = 6),
                                    box(numericInput(ns("missval_heatmap_Height"),
                                                     "height",
                                                     min = 1, max = 30, value = 10), width = 6)),
                                  fluidRow(
                                    plotOutput(ns("missval"), height = 600),
                                    downloadButton(ns('downloadMissval'), 'Save')
                                  ),
                                  shinyBS::bsTooltip(ns("cluster_columns_for_missval_heatmap"), "Whether make cluster on columns", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("Custom_columns_order_for_missval_heatmap"), "Order of column, act when [Cluster columns] is FALSE", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("missval_heatmap_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                                  shinyBS::bsTooltip(ns("missval_heatmap_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                         ),
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

                       )
                )
              ))
    # ,change_reactive_ID(id)
  )

}

DEP_pg_server_module <- function(id){
  moduleServer(
    id,
    function(input,output,session){

      cat("DEP_server_module load")

      ### UI functions for analyze options ### --------------------------------------------------------

      output$name <- renderUI({
        if(input$uploadmode == "fromTable"){
          shiny::validate(need(!is.null(data()), "Please upload quantity file in 'Files' firstly!"))
        }else{
          validate(need(!is.null(data()), "Please upload log in 'Files' firstly!"))
        }

        selectizeInput(session$ns("name"),
                       "Name column",
                       choices=colnames(data()),
                       selected = "Gene.names")
      })

      output$id <- renderUI({
        selectizeInput(session$ns("id"),
                       "ID column",
                       choices=colnames(data()),
                       selected = "Protein.IDs")
      })

      output$peptidescol <- renderUI({
        # cat("peptidescol")
        if("Peptides" %in% colnames(data())){
          return(selectizeInput(session$ns("peptidescol"),
                                "Peptide column",
                                choices=c(colnames(data()), ""),
                                selected = "Peptides"))
        }else{
          return(selectizeInput(session$ns("peptidescol"),
                                "Peptide column",
                                choices=c("", colnames(data())),
                                selected = NULL))
        }
      })

      output$intensitycols <- renderUI({
        LFQcols = grep("^LFQ",colnames(data()),value = T)
        if(length(LFQcols)==0) LFQcols = grep("^Intensity",colnames(data()),value = T)
        if(length(LFQcols)==0) LFQcols = grep("^intensity",colnames(data()),value = T)
        if(length(LFQcols)==0) LFQcols = grep(".Quantity",colnames(data()),value = T)
        # if(length(LFQcols)==0) LFQcols = grep("^quantity|Quantity",colnames(data()),value = T)
        # selectizeInput("intensitycols",
        #                " ",#Expression columns
        #                choices=colnames(data()),
        #                multiple = TRUE,
        #                selected = grep("^LFQ",LFQcols,value = T), width = '100%')
        # shinydashboardPlus::box(
        #   selectizeInput(session$ns("intensitycols"),
        #                                        label = "Choose the expression columns",
        #                                        choices=colnames(data()),
        #                                        multiple = TRUE,
        #                                        selected = grep("^LFQ",LFQcols,value = T), width = '100%')
        #   , title = "Expression columns", width = 12, icon = icon("th"), collapsible = TRUE, collapsed = TRUE, background = "light-blue")
        bsCollapse(id = session$ns("IntensityCols"), open = "Expression columns",
                   bsCollapsePanel("Expression columns",
                                   selectizeInput(session$ns("intensitycols"),
                                                  label = "Choose the expression columns",
                                                  choices=colnames(data()),
                                                  multiple = TRUE,
                                                  selected = LFQcols, width = '100%'),
                                   checkboxInput(inputId = session$ns("remove_prefix"),"remove prefix of label", value = T),
                                   checkboxInput(inputId = session$ns("remove_suffix"),"remove suffix of label", value = T)
                                   , style = "primary"))
      })

      output$filt <- renderUI({
        selectizeInput(session$ns("filt"),
                       "Filter on columns" ,
                       colnames(data()),
                       multiple = TRUE,
                       selected = c("Reverse","Potential.contaminant"))
      })

      output$order <- renderUI({
        if(input$uploadmode == "fromTable"){
          validate(need(!is.null(input$file1), "Please upload quantity file in Files"))
        }else{
          validate(need(!is.null(input$file_log), "Please upload log in Files"))
        }

        validate(need(!is.null(input$intensitycols), "Please select the Expression columns"))
        validate(need(length(input$intensitycols) > 1 , "More expression columns is required"))
        # groups = exp_design()[,input$groupby]
        # if (input$anno == "columns" & !is.null(data()) ) {
        #   condition = expdesign()$condition
        # } else {
        #   if (input$anno == "expdesign" & !is.null(expdesign()) ) {
        #     condition = make.names(expdesign()$condition)
        #   }
        # }

        condition = make.names(expdesign()$condition) %>% unique()

        selectizeInput(session$ns("order"),
                       "Set the groups order",
                       choices = condition, multiple = TRUE,
                       selected = NULL, width = '100%')
      })

      output$control <- renderUI({
        if(input$uploadmode == "fromTable"){
          validate(need(!is.null(input$file1), "Please upload quantity file in Files"))
        }else{
          validate(need(!is.null(input$file_log), "Please upload log in Files"))
        }

        validate(need(!is.null(input$intensitycols), "Please select the Expression columns"))
        validate(need(length(input$intensitycols) > 1 , "More expression columns is required"))

        if(input$contrasts == "control"){
          # if (input$anno == "columns" & !is.null(data()) ) {
          #   condition = make.names(expdesign()$condition)
          # } else {
          #   if (input$anno == "expdesign" & !is.null(expdesign()) ) {
          #     condition = make.names(expdesign()$condition)
          #   }
          # }

          condition = make.names(expdesign()$condition) %>% unique()

          selectizeInput(session$ns("control"),
                         "Control",
                         choices = condition,
                         selected = NULL)
        }
      })

      output$test_manual <- renderUI({
        validate(need(!is.null(input$file1), ""))
        if(!is.null(data()) & input$contrasts == "manual"){
          cols <- which(colnames(data()) %in% input$intensitycols)#according to intensitycols
          if(input$remove_prefix){
            prefix <- get_prefix(data()[,cols] %>% colnames())
            label = colnames(data())[cols] %>% gsub(prefix,"",.)
          }else label = colnames(data())[cols]

          test_manual_name <- unique(make.names(unlist(lapply(label %>% strsplit(., split = "_"), function(x){x[1]}))))
          test_manual_name <- cbind(combn(test_manual_name,2),combn(test_manual_name,2, FUN = rev))
          test_manual_name <- apply(test_manual_name, 2, function(i){paste(i[1], i[2], sep = "_vs_")})
          selectizeInput(session$ns("test_manual"), "Manual test",
                         choices=test_manual_name, selected = NULL, multiple = TRUE)
        }
      })

      output$Peptides1 <- renderUI({
        if(input$if_peptide_color) {
          colourpicker::colourInput(inputId = session$ns("Peptides_1"),
                                    label = "Peptides1",
                                    showColour = "both",
                                    palette = "square",
                                    value = "#A020F0")
        }
      })

      output$Peptides2 <- renderUI({
        if(input$if_peptide_color) {
          colourpicker::colourInput(inputId = session$ns("Peptides_2"),
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

        Columns_to_updata = Columns_to_updata(1)
        message("Update options 1 finished \n")
      },ignoreNULL = T,ignoreInit = T)

      observeEvent(
        # input$contrl,
        # input$order,
        input$peptidescol,
        {
          Columns_rendered = Columns_rendered(1)
        },ignoreNULL = T,ignoreInit = T
      )

      observeEvent(
        c(Columns_to_updata(),Columns_rendered()),
        {
          if(Columns_to_updata() == 1 & Columns_rendered() == 1){
            cat("Update options in 'Colmuns' \n")

            column_opt = c("name", "id", "delim", "intensitycols", "peptidescol",
                           "filt", "order", "control")
            for(i in column_opt){
              # message("update ",i," to '", upload_log()$inputVals()[[i]],"'\n", sep = "")
              upload_log22 <<- upload_log
              if(i == "intensitycols") temp11 <<- upload_log()$inputVals()[[i]]
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
                             choices = list("Parse from columns" = "columns",
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
          insertTab(inputId = "DEP_QC_tabs",
                    tabPanel(title = "Normal distribution - Fit",
                             fluidRow(
                               box(uiOutput(session$ns("cntrst_norm_distribution")), width = 8)
                             ),
                             fluidRow(
                               box(numericInput(session$ns("norm_distribution_Width"),
                                                "width",
                                                min = 1, max = 30, value = 7), width = 6),
                               box(numericInput(session$ns("norm_distribution_Height"),
                                                "height",
                                                min = 1, max = 30, value = 7), width = 6)
                             ),
                             fluidRow(
                               plotOutput(session$ns("norm_distribution_plot"),height=600),
                               downloadButton(session$ns('download_norm_distribution_plot'), 'Save')
                             ),
                             shinyBS::bsTooltip(session$ns("cntrst_norm_distribution"), "Choose the contrast to plot", "top", options = list(container = "body")),
                             shinyBS::bsTooltip(session$ns("norm_distribution_Width"), "Width of the figure to export", "top", options = list(container = "body")),
                             shinyBS::bsTooltip(session$ns("norm_distribution_Height"), "Height of the figure to export", "top", options = list(container = "body"))
                    ),

                    target = "Pca plot", position = "before"
          )
        }
      )

      observe(
        if(input$threshold_method == "intersect") {
          print("observe2")
          removeTab(inputId = "DEP_QC_tabs", target = "Normal distribution - Fit")
        }
      )

      ### check input ### --------------------------------------------------------
      iv <- InputValidator$new()
      iv$add_rule("p", sv_between(0.0001, 0.1))
      iv$add_rule("lfc", sv_between(0, 10))
      iv$add_rule("curvature", sv_between(0.1, 6))
      iv$add_rule("x0_fold", sv_between(1, 10))
      iv$enable()


      ### Reactive functions ### --------------------------------------------------
      upload_log <- reactive({
        logFile <- input$file_log

        if( !is.null(logFile) ){

          cat("Read infile_log\n")
          upload_log <- readRDS(logFile$datapath)
          upload_log11 <<- upload_log

          if( upload_log$appVersion == DEP2:::app_version && inherits(upload_log,"pg_log") ){
            sendSweetAlert(
              session = shiny::getDefaultReactiveDomain(),
              title = "Load log successfully!",
              text = "Pass version check and load.",
              type = "success"
            )
            message("Read successfully\n")
            return(upload_log)

          }else if(!inherits(upload_log,"pg_log")){
            sendSweetAlert(
              session = shiny::getDefaultReactiveDomain(),
              title = "Error !",
              text = "The uploaded file is not a log for protein-group analysis, unload the log.\n
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
              text = "The version of uploaded file don't match current app version now\n
              Please ensure that: The log file is a log from protein-group analysis. And the app version is correct,
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

      ## generate a experiement design table
      expdesign <- reactive({
        # from log
        if(input$uploadmode == "fromLog" && !is.null(data()) ){
          upload_log = upload_log()

          if(input$anno == "logexpdesign"){
            cat("Use expdesign from log\n")
            expdesign <- upload_log$resultVals()$expdesign
            # expdesign11 <<- expdesign
          }else if(input$anno == "columns"){
            cols <- which(colnames(data()) %in% input$intensitycols)

            if(length(cols) == 0){
              return(NULL)
            }else{
              cat("Parse from columns\n")
              label <- colnames(data())[cols]
              expdesign <- get_exdesign_parse(label, mode = "delim", sep = "_",
                                              remove_prefix = input$remove_prefix, remove_suffix = input$remove_suffix)
            }
          }
        }else if(input$uploadmode == "fromLog" && is.null(data())){
          return(NULL)
        }

        # from table
        if(input$uploadmode == "fromTable"){
          inFile <- input$file2
          if ( input$anno == "columns" & (!is.null(data())) ){
            cat("Parse from columns.\n")
            cols <- which(colnames(data()) %in% input$intensitycols)

            if(length(cols) == 0){
              return(NULL)
            }else{
              label <- colnames(data())[cols]
              expdesign <- get_exdesign_parse(label, mode = "delim", sep = "_",
                                              remove_prefix = input$remove_prefix, remove_suffix = input$remove_suffix)
            }
          }else if(input$anno == "expdesign" & (!is.null(inFile)) & (!is.null(data())) ){
            cat("Use expdesign from upload table\n")
            expdesign <- read.csv(inFile$datapath, header = TRUE,
                                  sep = "\t", stringsAsFactors = FALSE) %>%
              mutate(id = row_number())
          }else{
            return(NULL)
          }
        }

        my_expdesign <<- expdesign
        return(expdesign)
      })


      data <- reactive({

        # from log
        upload_log = upload_log()
        if(!is.null(upload_log) & input$uploadmode == "fromLog"){
          my_data <<- upload_log$resultVals()$data
          return(my_data)
        }else if( is.null(upload_log) & input$uploadmode == "fromLog"){
          # cat("Initiate: log file is null \n")
          return(NULL)
        }

        # from table
        inFile <- input$file1
        if (is.null(inFile)){
          # cat("Initiate: infile1 is null \n")
          return(NULL)
        }
        # cat("infile1 is not null now \n")
        my_data <- data.table::fread(inFile$datapath, header = TRUE,
                                     sep = "\t", stringsAsFactors = FALSE, integer64 = "numeric") %>%
          mutate(id = row_number())
        coln.int64 <- names(which(sapply(my_data, bit64::is.integer64)))
        if (length(coln.int64) > 0L)
          my_data[, c(coln.int64) := lapply(.SD, as.numeric), .SDcols = coln.int64]
        my_data = as.data.frame(my_data)
        ## replace the [xxx] characters in colnames
        colnames(my_data) = gsub("^\\[(.*)\\] ","",colnames(my_data))
        colnames(my_data) = make.names(colnames(my_data))

        my_data
      })



      # filt0
      filt <- reactive({
        data <- data()

        cols <- which(colnames(data) %in% input$intensitycols)
        data[,cols] = DEP2::clean_character(data[,cols])

        validate(need(!(input$name == "" && input$id == ""), "Please ensure that: at least one of your name column and id column is non-empty!"))

        if(input$name == "" & input$id == "") {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = "Please ensure that: at least one of your name column and id column is non-empty!",
            type = "warning"
          )
        }

        unique_names <- DEP2::make_unique(proteins = data, names = ifelse(input$name == "", input$id, input$name), ids = ifelse(input$id == "", input$name, input$id), delim = input$delim)

        ind_empty = c(grep("^\\.\\d*$", unique_names$name), which(unique_names$name == ""))
        if(length(ind_empty) > 0) {
          unique_names = unique_names[-ind_empty, ]
        }

        # if (input$anno == "columns") {
        #   se <- DEP2::make_se_parse(unique_names, cols, mode = "delim", sep = "_",
        #                             remove_prefix = input$remove_prefix, remove_suffix = input$remove_suffix)
        # }
        # if (input$anno == "expdesign") {
        #   se <- DEP2::make_se(unique_names, cols, expdesign())
        #   colData(se)$replicate = as.character(colData(se)$replicate)
        # }

        # make se
        se <- DEP2::make_se(unique_names, cols, expdesign()) ## when input$anno == "columns", use expdesign() and make_se instead of make_se_parse
        colData(se)$replicate = as.character(colData(se)$replicate)

        # filter by identify info
        filtered <- se

        if(is.null(input$filt)){
          message("the filter column is empty! Do not filter with column")
        }else{
          cols_filt <- intersect(input$filt, colnames(data))
          if(length(cols_filt) > 0){
            for(i in cols_filt){
              filter_formula = paste( "~(is.na(",i,")|",i,"=='')") %>% as.formula()
              message(filter_formula)
              filtered <- DEP2:::filter_se(filtered, filter_formula = filter_formula)
            }
          }
        }

        # filter by missing number
        thr <- ifelse(is.na(input$thr), 0, input$thr)
        my_filt <- filter_se(filtered, thr = thr)

        return(my_filt)
      })

      iv1 <- InputValidator$new()
      iv1$add_rule("order",
                   sv_required(message = "Complete order is required, please sort all groups", test = function(x){
                     is.null(x) || length(x) == (length(filt()@colData$condition %>% unique))
                   })
      )
      iv1$enable()

      the_order <- reactive({
        req(iv1$is_valid())

        if(length(input$order) == length(filt()@colData$condition %>% unique))
          return(input$order)
        return(NULL)
      })

      # filt <- reactive({
      #   filt = filt0()
      #
      #   if(is.null(the_order())){
      #     return(filt)
      #   }else{
      #     filt = DEP2::Order_cols(filt,the_order())
      #
      #     return(filt)
      #   }
      # })

      norm <- reactive({
        my_norm <- try({normalize_vsn(filt())})
        if(class(my_norm) == "try-error"){
          ## if vsn failed for a small matrix, skip vsn normalization
          my_norm = filt()
        }
        my_norm
      })

      imp <- reactive({
        # inFile1 <- input$resultRData

        upload_log = upload_log()
        if( !is.null(upload_log) & input$uploadmode == "fromLog" ){
          options_before_imp <- c("name", "id", "delim", "intensitycols",
                                  "peptidescol", "filt", "imputation")

          check_opt_change <- options_before_imp %>% sapply(.,function(x){
            identical(input[[x]], upload_log()$inputVals()[[x]])
          })

          if(all(check_opt_change)){ # if the option before imputation unchange, directly use the old imputation result
            my_imp = upload_log()$resultVals()$imp
            message("Use the imputation in log file!")
            return(my_imp)
          }

        }


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
          my_imp <- DEP2::impute(
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
            my_imp <- sample_specific_imputation
            validate(
              need(length(which(is.na((assay(my_imp))))) == 0, "Please select another imputation method")
            )
            my_imp
          } else {
            my_imp <- impute(norm(), input$imputation)
          }
        }

        message("Use the imputation:", input$imputation)
        return(my_imp)
      })



      FDR_type <- reactive({
        FDR_type <- input$FDR_type
        # FDR_type <- ifelse(length(grep(FDR_type,"qval"))>1,"qval","lfdr")
      })

      df <- reactive({
        if(input$contrasts == "control") {
          validate(
            need(input$control != "", "Please select a Control condition under menuItem Columns in the DEP-LFQ options of the sidebar"),
            need( !is.null(input$intensitycols) && all(input$intensitycols %in% colnames(data())), "Please select the Expression columns in the sidebar")
          )
        }

        if(input$contrasts == "manual") {
          validate(
            need(input$test_manual != "", "Please select manual contrasts to test under menuItem Columns"),
            need(!is.null(input$intensitycols) && all(input$intensitycols %in% colnames(data())), "Please select the Expression columns in the sidebar")
          )
        }

        inFile1 <- input$resultRData
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
            load(file = inFile1$datapath)
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


      dep0 <- reactive({
        req(input$threshold_method %in% c("intersect", "curve"))
        if(input$threshold_method == "intersect") {
          my_dep <- DEP2::add_rejections(diff = df(), alpha = input$p, lfc = input$lfc, thresholdmethod = input$threshold_method)
        }
        if(input$threshold_method == "curve") {
          my_dep <- DEP2::add_rejections(diff = df(), thresholdmethod = input$threshold_method, curvature = input$curvature, x0_fold = input$x0_fold)
        }

        if(is.null(input$peptidescol)|input$peptidescol==""){
          rowData(my_dep)$Peptides = 3
        }else{
          rowData(my_dep)$Peptides = rowData(my_dep)[,input$peptidescol]
        }

        return(my_dep)
      })

      dep <- reactive({
        dep = dep0()

        if(is.null(the_order())){
          return(dep)
        }else{
          dep = DEP2::Order_cols(dep,the_order())

          return(dep)
        }
      })

      output$Save_RData <- downloadHandler(
        filename = function() { paste("results", ".RData", sep = "") },
        content = function(file) {
          withProgress(message = 'Please wait ...', value = 0.66, {
            data = data(); filt = filt(); norm = norm();imp = imp(); dep = dep()
            save(data = data, filt = filt, norm, imp, dep, file=file)})}
      )

      output$Save_log <- downloadHandler(
        filename = function() { paste("Export_log_pg", format(Sys.time(), "%m%d_%H%M"),
                                      ".logrds", sep = "") },
        content = function(file) {
          withProgress(message = 'Compressing result. Please wait ...', value = 0.66, {
            thelog = pg_log$new(input = input,
                                reactive_vals = list(
                                  data = data(),
                                  expdesign = expdesign(),
                                  filt = filt(),
                                  norm = norm(),
                                  imp = imp(),
                                  dep = dep()
                                ),
                                app_version = DEP2:::app_version
            )
            thelog11 <<- thelog
            saveRDS(thelog, file = file)
          })
        }
      )


      ## All object and functions upon 'Analyze': observeevent  ### ----------------------
      observeEvent(input$analyze, {

        if(is.null(input$intensitycols)) {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = "Please ensure that: select expression columns in the Columns opt",
            type = "warning"
          )
        }
        req(!is.null(input$intensitycols))

        if(is.null(input$filt)) {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "Note :",
            text = "Your filter column is empty! If needed, please choose the filter columns. Otherwise, please ignore this message",
            type = "info"
          )
        }
        print(input$name)


        ### Interactive UI functions in analyze, download, output & plots### ------------------------------------------
        output$downloadTable <- renderUI({
          selectizeInput(session$ns("dataset"),
                         "Choose a dataset to save" ,
                         c("results","significant_proteins",
                           "displayed_subset","full_dataset"))
        })

        output$downloadButton <- renderUI({
          downloadButton(session$ns('downloadData'), 'Save', class = "downloadData")
        })

        output$downloadButton_for_save_RData <- renderUI({
          downloadButton(session$ns("Save_RData"), "save result RData", class = "Save_RData")
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
            info_box <- infoBox("Significant proteins",
                                paste0(num_signif,
                                       " out of ",
                                       num_total),
                                paste0(signif(frac * 100, digits = 3),
                                       "% of proteins differentially expressed"),
                                icon = icon("thumbs-up", lib = "glyphicon"),
                                # color = "green",
                                width = 4)
          }
          info_box
        })

        output$select <- renderUI({
          row_data <- rowData(dep())
          cols <- grep("_significant", colnames(row_data))
          names <- colnames(row_data)[cols]
          names <- gsub("_significant", "", names)
          selectizeInput(session$ns("select"),
                         "Select direct comparisons",
                         choices=names,
                         multiple = TRUE)
        })

        output$exclude <- renderUI({
          row_data <- rowData(dep())
          cols <- grep("_significant", colnames(row_data))
          names <- colnames(row_data)[cols]
          names <- gsub("_significant","",names)
          selectizeInput(session$ns("exclude"),
                         "Exclude direct comparisons",
                         choices = names,
                         multiple = TRUE)
        })


        output$heatmap_cntrst <- renderUI({
          if (!is.null(selected())) {
            df <- rowData(selected())
            cols <- grep("_significant$",colnames(df))
            selectizeInput(session$ns("heatmap_cntrst"),
                           "Contrast",
                           choices = gsub("_significant", "", colnames(df)[cols]), multiple = TRUE)
          }
        })


        output$mysplit <- renderUI({
          selectizeInput(session$ns("mysplit"),
                         "my split",
                         choices = c(1 : input$k), multiple = TRUE)

        })

        output$Custom_columns_order <- renderUI({
          selectizeInput(session$ns("Custom_columns_order"),
                         "Custom columns order",
                         choices = colData(filt())$condition %>% unique(), multiple = TRUE)

        })

        output$Custom_columns_order_for_missval_heatmap <- renderUI({
          selectizeInput(session$ns("Custom_columns_order_for_missval_heatmap"),
                         "Custom columns order",
                         choices = colnames(assay(dep())), multiple = TRUE)
        })


        output$Volcano_cntrst <- renderUI({
          if (!is.null(selected())) {
            df <- rowData(selected())
            cols <- grep("_significant$",colnames(df))
            selectizeInput(session$ns("Volcano_cntrst"),
                           "Contrast",
                           choices = gsub("_significant", "", colnames(df)[cols]))
          }
        })

        output$cntrst_norm_distribution <- renderUI({
          if (!is.null(selected())) {
            df <- rowData(selected())
            cols <- grep("_significant$",colnames(df))
            selectizeInput(session$ns("cntrst_norm_distribution"),
                           "Contrast",
                           choices = gsub("_significant", "", colnames(df)[cols]))
          }
        })

        output$selected_proteins <- renderUI({
          row_Data <- rowData(dep())
          selectizeInput(session$ns("selected_proteins"),
                         "selected proteins",
                         choices = row_Data$name, selected = NULL, multiple = TRUE)
        })

        output$chooseToshow <- renderUI({
          selectizeInput(session$ns("chooseToshow"),
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

        output$umap_n_neighbors <- renderUI({
          cat("umap_n_neighbors")
          sliderInput(inputId = session$ns("umap_n_neighbors"),"number of nearest neighbors",min=2,max=min(15,ncol(dep())),value=4)}
        )

        output$pca_top_n  <- renderUI(
          sliderInput(inputId = session$ns("pca_top_n"),"number of top variable proteins to consider",min = 2, max = nrow(dep()), value=500)
        )

        output$Tsne_perplexity <- renderUI(
          sliderInput(inputId = session$ns("Tsne_perplexity"),"perplexity",min=1,max= (ncol(dep())-1)/3 - 2^(-31) ,value=2)
        )

        ### Reactive functions ### ------------------------------------------------
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
          DEP2::get_results(selected())
        })

        table <- reactive({
          DEP2::get_results(selected())
        })

        selected_plot_input <- reactive ({
          if(!is.null(input$table_rows_selected)) {
            selected_id <- table()[input$table_rows_selected,1]
            DEP2::plot_single(selected(), selected_id, input$pres) ## function from DEP package
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
                               # manual = input$manual_heatmap,
                               manual_contrast = input$heatmap_cntrst,
                               kmeans = TRUE,
                               k = input$k,
                               color = input$colorbar,
                               col_limit = input$limit,
                               #same_trend = input$same_trend,
                               row_font_size = input$row_font_size,
                               col_font_size = input$col_font_size,
                               cluster_columns = input$cluster_columns,
                               # if_mysplit = input$if_mysplit,
                               split_order = input$mysplit,
                               label_few_peptide_rows = input$if_rowname_color,
                               # if_chooseToshow = input$if_chooseToshow,
                               chooseToshow = input$chooseToshow
                               # column_order = if(input$cluster_columns) {NULL} else {input$Custom_columns_order}
            )
          })
        })

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
                                 down_color = input$down_color,
                                 stable_color = input$stable_color,
                                 up_color = input$up_color,
                                 peptide_1_color = input$Peptides_1,
                                 peptide_2_color = input$Peptides_2
              )
            })
          }
        })

        stastical_plot_input <- reactive({
          if(!is.null(input$plot_stat_contrasts)){
            DEP2::plot_statistics(selected(), statistic = input$stastical_plot_x ,contrasts = input$plot_stat_contrasts)
          }
        })

        norm_distribution_input <- reactive({
          if(!is.null(input$cntrst_norm_distribution)){
            DEP2::plot_norm_distribution(object = dep(),
                                         contrast = input$cntrst_norm_distribution)
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
          DEP2::plot_normalization(filt(),
                                   norm())
        })

        missval_input <- reactive({
          DEP2::plot_missval(norm(),
                             cluster_columns = input$cluster_columns_for_missval_heatmap,
                             column_order = if(input$cluster_columns_for_missval_heatmap) {NULL} else {input$Custom_columns_order_for_missval_heatmap}
          )
        })

        detect_input <- reactive({
          DEP2::plot_detect(norm())
        })

        imputation_input <- reactive({
          DEP2::plot_imputation(norm(),
                                df())
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
            plotOutput(session$ns("heatmap"), height = (100 * as.numeric(input$size)))
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

        output$download_custom_volcano <- downloadHandler(
          filename = function() {
            paste0("custom_Volcano_", input$Volcano_cntrst, ".pdf")
          },
          content = function(file) {
            pdf(file, width = input$custom_volcano_Width, height = input$custom_volcano_Height)
            print(custom_volcano_input())
            dev.off()
          }
        )

        output$downloadstastical_plot <- downloadHandler(
          filename = function() {
            paste0("statistical_plot_", paste0(input$plot_stat_contrasts, collapse = "_"), ".pdf")
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

      ###"question for DEP"
      observeEvent(input$help_format_DEP, {
        showModal(modalDialog(
          title = "Format specifications",
          # includeMarkdown(system.file("extdata", "DEP_LFQ.md", package = "DEP2")),
          includeMarkdown("www/DEP_LFQ.md"),
          br(),
          # h4("Example:"),
          # tags$img(
          #   # src = base64enc::dataURI(file = system.file("extdata", "DEP_LFQ.png", package = "DEP2"), mime = "image/png"),
          #   src = base64enc::dataURI(file = "www/DEP_LFQ.png", mime = "image/png"),
          #   width = 750
          # ),
          easyClose = TRUE,
          footer = NULL,
          size = "l"
        ))
      })

      ### help imputation
      observeEvent(input$help_imputation, {
        showModal(modalDialog(
          title = "The detailed information of imputation methods",
          # includeMarkdown(system.file("extdata", "impute.md", package = "DEP2")),
          includeMarkdown("www/impute.md"),
          easyClose = TRUE,
          footer = NULL,
          size = "l"
        ))
      })

      return(df)


    }
  )
}


