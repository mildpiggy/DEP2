## the sidebar_mod body_mod and server mod of DEP-ptm analysis from Phospho (STY)Sites quantification

DEP_ptm_UI <- function(id){
  tagList(
    sidebarLayout(
      DEP_ptm_sidebar_mod(id = id), ## DEP-ptm sidebar module
      DEP_ptm_body_mod(id = id)  ## DEP-ptm pagebody module
    ))
}

DEP_ptm_sidebar_mod <-  function(id){
  ns = NS(id)
  cat(paste0("DEP_ptm_sidebar ns is",ns(""),"\n"))

  tagList(
    # sidebarLayout(
    sidebarPanel(
      width = 3,
      h3(id),
      bsCollapse(
        bsCollapsePanel("Files",
                        style = "primary",
                        fluidRow(
                          #column(width = 2,
                          #actionButton("help_format_DEP",width = 12,
                          #              label = "", icon = icon("bell")
                          #      )),
                          column(width = 9,
                                 fileInput(ns('file1_ptm'),width = "300px",
                                           'PTM Sites.txt',
                                           accept=c('text/csv',
                                                    'text/comma-separated-values,text/plain',
                                                    '.csv')),
                                 fileInput(ns('file2_ptm'),width = "300px",
                                           'ExperimentalDesign.txt',
                                           accept=c('text/csv',
                                                    'text/comma-separated-values,text/plain',
                                                    '.csv')),
                                 fileInput(ns('resultRData_ptm'),'load saved result RData',width = "300px")
                          ),
                          column(
                            width = 1, h4(),
                            br(),
                            actionButton(ns("help_format_DEP_ptm"),
                                         label = "", icon = icon("question-circle"),
                                         style = "color: #f6f6f6; background-color: #2c3b41; border-color: #2c3b41"
                            ),
                            shinyBS::bsTooltip(
                              ns("help_format_DEP_ptm"),
                              "How to provide your input data",
                              "top",
                              options = list(container = "body")
                            )
                          )),
                        radioButtons(ns("anno_ptm"),
                                     "Sample annotation",
                                     choices = list("Parse from columns" = "columns",
                                                    "Use Experimental Design" = "expdesign"),
                                     selected = "columns")

        ),
        # menuItemOutput("columns"),
        bsCollapsePanel("Columns",
                        style = "primary",
                        uiOutput(ns("name_ptm")),
                        uiOutput(ns("id_ptm")),
                        uiOutput(ns("amino_acid_ptm")),
                        uiOutput(ns("aa_position_ptm")),
                        selectizeInput(ns("delim_ptm"), "Delimiter", choices = c(";", "|"), selected = ";"),
                        uiOutput(ns("intensitycols_ptm")),
                        uiOutput(ns("filt_ptm")),#input.peptidescol !== ''
                        uiOutput(ns("filt_num_cutoff_ptm")),
                        conditionalPanel(condition = paste0("input['",ns("filt_num_cutoff_ptm"),"']","!== ''"),#&& "input.filt_num_cutoff_ptm !== ''"
                                         numericInput(ns("filt_num_ptm"),
                                                      "Cutoff",
                                                      min = 0, max = 1, value = 0.7),

                        ),
                        numericInput(ns("thr_ptm"), "Allowed max.miss.num at least one condition", min = 0, max = 20, value = 0),
                        uiOutput(ns("order")),
                        uiOutput(ns("control_ptm")),
                        uiOutput(ns("test_manual_ptm"))
        ),
        bsCollapsePanel("Imputation options",
                        style = "primary",
                        radioButtons(inputId = ns("imputation_ptm"),
                                     label = "Imputation type",
                                     choices = c(c("man", "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "min", "zero"), "mixed on proteins", "mixed on samples","RF"),  ## increase "RF" opt 20211023
                                     selected = "MinProb"),
                        shinyBS::bsTooltip(ns("imputation_ptm"), "Choose an imputation method", "right", options = list(container = "body")),
                        actionButton(ns("help_imputation_ptm"),
                                     label = "Detailed information", #icon = icon("question-circle"),
                                     style = "color: #f6f6f6; background-color: #2c3b41; border-color: #2c3b41; text-align: left; margin: 0 0 0px 0px",
                                     #rgba(255, 255, 255, 0.1) rgba(0, 0, 0, 0.1)
                        ),
                        shinyBS::bsTooltip(
                          ns("help_imputation_ptm"),
                          "The detailed information of imputation methods",
                          "right",
                          options = list(container = "body")
                        )
                        # p(a("Detailed information link",
                        #     href = "https://www.rdocumentation.org/packages/MSnbase/versions/1.20.7/topics/impute-methods",
                        #     target="_blank"))
        ),
        bsCollapsePanel("Protein quantification correction",
                        style = "primary",
                        prettySwitch(
                          inputId = ns("Protein_level_correction"),
                          label = "Correct quantity by protein-level quantity",
                          status = "success",
                          fill = F
                        ),
                        uiOutput(ns("protein_group_choose")),
                        uiOutput(ns("Protein_level_correction_options")),
                        shinyBS::bsTooltip(ns("Protein_level_correction"),
                                           "If on, enriched(modified)-peptide intensity will be Corrected by protein-level qutification before differential test")

        ),
        bsCollapsePanel("FDR correction",  ## DEP FDR correct options
                        style = "primary",
                        radioButtons(ns("FDR_type_ptm"),
                                     "type of FDR value",
                                     choices = c("Strimmer's qvalue(t)","Strimmer's qvalue(p)","BH","Storey's qvalue")),
                        shinyBS::bsTooltip("FDR_type_ptm", "Choose the method of pvalue adjustment", "right", options = list(container = "body"))
        ),
        bsCollapsePanel("Threshold method",
                        style = "primary",
                        radioButtons(ns("threshold_method_ptm"),
                                     "Threshold method",
                                     choices = c("intersect", "curve"),
                                     selected = "intersect"),
                        shinyBS::bsTooltip(ns("threshold_method_ptm"), "Choose the method that the cutoff of significant peptides based on. [intersect] means by adjusted pvalue and log2 fold change cutoff and [curve] means by curvature and x0 fold cutoff", "right", options = list(container = "body"))
        ),

        open = "Files"
      ),

      shinyBS::bsTooltip(ns("anno_ptm"), "Sample annotation type", "top", options = list(container = "body")),
      shinyBS::bsTooltip(ns("id_ptm"), "Name of the column containing feature IDs", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("delim_ptm"), "Set the delimiter separating the name and id column", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("intensitycols_ptm"), "Choose the expression columns of your data", "right", options = list(container = "body")),
      # shinyBS::bsTooltip("peptidescol", "Choose the peptides column of your data. And this column should be numeric. It is only used to plot heatmap and custom volcano to color the peptides number 1 and 2 to help user check the accuracy of the quantification. Of course, if your data does not have this column, You can leave it blank", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("filt_ptm"), "The filtered columns based on", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("filt_num_cutoff_ptm"), "Set the column that cutoff based on (such as the location prosibility/score), and it should be numeric.", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("filt_num_ptm"), "Set the cutoff [0 ~ 1], the rows that greater than or equal to it will be kept", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("thr_ptm"), "Set the threshold for the allowed number of missing values in at least one condition", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("control_ptm"), "Choose your control condition", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("test_manual_ptm"), "Choose the contrasts that you want to test", "right", options = list(container = "body")),
      actionButton(ns("analyze_ptm"), "Analyze"),

      shinyBS::bsTooltip("analyze_ptm", "Click on it to analyze your data", "right", options = list(container = "body")),
      tags$hr(),
      tags$style(type="text/css", "#downloadData_ptm {background-color:white;color: black;font-family: Source Sans Pro}"),
      uiOutput(ns("downloadTable_ptm")),
      uiOutput(ns("downloadButton_ptm")),
      h6(),
      tags$style(type="text/css", "#Save_RData_ptm {background-color:white;color: black;font-family: Source Sans Pro}"),
      uiOutput(ns("downloadButton_for_save_RData_ptm")),
      shinyBS::bsTooltip(ns("downloadTable_ptm"), "Choose a dataset to save, and here we offer four forms of datasets for downloading including results.txt, significant_proteins.txt, displayed_subset.txt, full_dataset.txt", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("downloadButton_ptm"), "Click on it to download the selected dataset", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("downloadButton_for_save_RData_ptm"), "Click on it to Save the RData, and when you upload the results.RData in the menuItem [Files] next time, you can get the same result as this time", "right", options = list(container = "body"))



      # ),
      # mainPanel("bbb")
    )
  )
}

DEP_ptm_body_mod <- function(id){
  ns = NS(id)
  cat(paste0("DEP_ptm_body_mod ns is",ns(""),"\n"))
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
                # box(numericInput("p",
                #                  "adj. P value",
                #                  min = 0.0001, max = 0.1, value = 0.05),
                #     width = 2),
                # box(numericInput("lfc",
                #                  "Log2 fold change",
                #                  min = 0, max = 10, value = 1),
                #     width = 2),
                conditionalPanel(condition = paste0("input['",ns("threshold_method_ptm"),"']","== 'intersect'"), ## condition in mod
                                 box(numericInput(ns("p_ptm"),
                                                  "adj. P value",
                                                  min = 0.0001, max = 0.1, value = 0.05),
                                     width = 2),
                                 box(numericInput(ns("lfc_ptm"),
                                                  "Log2 fold change",
                                                  min = 0, max = 10, value = 1),
                                     width = 2)
                ),
                conditionalPanel(condition = paste0("input['",ns("threshold_method_ptm"),"']","== 'curve'"),
                                 box(numericInput(ns("curvature_ptm"),
                                                  ns("curvature"),
                                                  min = 0.1, max = 6, step = 0.1, value = 0.6),
                                     width = 2),
                                 box(numericInput(ns("x0_fold_ptm"),
                                                  "x0 fold",
                                                  min = 1, max = 10, value = 1),
                                     width = 2)
                ),
                # box(uiOutput("p"),
                #     uiOutput("curvature"), width = 2),
                # box(uiOutput("lfc"),
                #   uiOutput("x0_fold"), width = 2),
                # uiOutput("curvature"),
                # uiOutput("x0_fold"),
                infoBoxOutput(ns("significantBox_ptm")),
                box(radioButtons(ns("pres_ptm"),
                                 "Data presentation",
                                 c("contrast", "centered"),
                                 selected = "centered"),
                    width = 2),
                box(radioButtons(inputId = ns("contrasts_ptm"),
                                 label = "Contrasts",
                                 choices = c("control", "all", "manual"),
                                 selected = "control", inline = TRUE),
                    width = 2),
                shinyBS::bsTooltip(ns("p_ptm"), "Set the false discovery rate threshold", "top", options = list(container = "body")),
                shinyBS::bsTooltip(ns("lfc_ptm"), "Set the log2 fold change threshold", "top", options = list(container = "body")),
                shinyBS::bsTooltip(ns("curvature_ptm"), "Set the curvature threshold", "top", options = list(container = "body")),
                shinyBS::bsTooltip(ns("x0_fold_ptm"), "Set the x0 fold threshold, it decides the x0 ratio to the standard deviations of L2FC", "top", options = list(container = "body")),
                shinyBS::bsTooltip(ns("pres_ptm"), "The type of data scaling used for heatmap plotting. Either the fold change (contrast) or the centered log2-intensity (centered)", "top", options = list(container = "body")),
                shinyBS::bsTooltip(ns("contrasts_ptm"), "The type of contrasts that will be tested. This can be all possible pairwise comparisons (all), limited to the comparisons versus the control (control), or manually defined contrasts (manual)", "top", options = list(container = "body"))
              ),
              fluidRow(
                column(width = 7,
                       box(title = "Top Table",
                           box(uiOutput(ns("select_ptm")), width = 6),
                           box(uiOutput(ns("exclude_ptm")), width = 6),
                           DT::dataTableOutput(ns("table_ptm")), width = 12),
                       shinyBS::bsTooltip(ns("select_ptm"), "Select the rows that are significant of the selected contrasts", "top", options = list(container = "body")),
                       shinyBS::bsTooltip(ns("exclude_ptm"), "Exclude the rows that are significant of the selected contrasts", "top", options = list(container = "body")),
                       shinyBS::bsTooltip(ns("table_ptm"), "Table of the significant results", "top", options = list(container = "body"))
                ),
                column(width = 5,
                       tabBox(title = "Result Plots", width = 12,
                              tabsetPanel(id = ns("DEP_ptm_results_tabs"),
                                          tabPanel(title = "Selected Peptide",
                                                   fluidRow(
                                                     box(numericInput(ns("selected_plot_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7), width = 6),
                                                     box(numericInput(ns("selected_plot_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 6)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("selected_plot_ptm")),
                                                     downloadButton(ns('downloadPlot_ptm'), 'Save plot')
                                                   ),
                                                   shinyBS::bsTooltip(ns("selected_plot_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("selected_plot_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "Heatmap",
                                                   fluidRow(
                                                     column(width = 6,
                                                            box(checkboxInput(ns("manual_heatmap_ptm"),
                                                                              "Manual heatmap",
                                                                              value = FALSE), width = 12)
                                                     ),
                                                     column(width = 6,
                                                            box(selectizeInput(ns("colorbar_ptm"),
                                                                               "colorbar",
                                                                               choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                                                                               selected = c("RdBu"), multiple = FALSE), width = 12)
                                                     )
                                                   ),
                                                   fluidRow(
                                                     box(uiOutput(ns("heatmap_cntrst_ptm")), width = 12)
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
                                                     conditionalPanel(condition = paste0("!input['",ns("if_chooseToshow_ptm"),"']"),#&& "!input.if_chooseToshow" input.peptidescol !== '' typeof input.peptidescol !== null
                                                                      box(numericInput(ns("k_ptm"),
                                                                                       "Kmeans",
                                                                                       min = 0, max = 15, value = 6),
                                                                          width = 4)
                                                     ),
                                                     # box(numericInput("k",
                                                     #                  "Kmeans",
                                                     #                  min = 0, max = 15, value = 6),
                                                     #     width = 4),
                                                     box(numericInput(ns("limit_ptm"),
                                                                      "Color",
                                                                      min = 0, max = 16, value = 4),
                                                         width = 4),
                                                     box(numericInput(ns("size_ptm"),
                                                                      "Heat size(4-30)",
                                                                      min = 4, max = 30, value = 10),
                                                         width = 4)

                                                   ),
                                                   fluidRow(
                                                     dropdownButton(
                                                       # selectizeInput("colorbar",
                                                       #      "colorbar",
                                                       #      choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                                                       #      selected = c("RdBu"), multiple = FALSE),
                                                       checkboxInput(ns("if_chooseToshow_ptm"),
                                                                     "only display selected rows",
                                                                     value = FALSE),
                                                       uiOutput(ns("chooseToshow_ptm")),

                                                       numericInput(ns("row_font_size_ptm"),
                                                                    "row font size",
                                                                    min = 0, max = 15, value = 6),
                                                       numericInput(ns("col_font_size_ptm"),
                                                                    "col font size",
                                                                    min = 0, max = 16, value = 10),
                                                       numericInput(ns("Width_ptm"),
                                                                    "heatmap width",
                                                                    min = 1, max = 30, value = 7),
                                                       numericInput(ns("Height_ptm"),
                                                                    "heatmap height",
                                                                    min = 1, max = 30, value = 10),
                                                       checkboxInput(ns("cluster_columns_ptm"),
                                                                     "Cluster columns",
                                                                     value = TRUE),
                                                       uiOutput(ns("Custom_columns_order_ptm")),
                                                       conditionalPanel(condition = paste0("!input","['",ns("if_chooseToshow_ptm"),"']"), #"!input.if_chooseToshow"
                                                                        # checkboxInput(ns("if_mysplit_ptm"),
                                                                        #               "If mysplit",
                                                                        #               value = FALSE),
                                                                        uiOutput(ns("mysplit_ptm"))
                                                       ),
                                                       # checkboxInput("if_mysplit",
                                                       #               "If mysplit",
                                                       #               value = FALSE),
                                                       # uiOutput("mysplit"),
                                                       tooltip = tooltipOptions(title = "More parameters"),
                                                       icon = icon("gear"),
                                                       # status = "info",
                                                       width = "300px"),
                                                     br("")
                                                   ),
                                                   fluidRow(
                                                     uiOutput(ns("plot_ptm")),
                                                     downloadButton(ns('downloadHeatmap_ptm'), 'Save heatmap')),
                                                   shinyBS::bsTooltip(ns("colorbar_ptm"), "Choose the color palette of the colorbar", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("manual_heatmap_ptm"), "Whether plot the manual heatmap", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("heatmap_cntrst_ptm"), "Choose the contrasts that you want to show, act when [Manual heatmap] is TRUE, can be one or more contrasts", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("k_ptm"), "Set the number of k-means clusters", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("limit_ptm"), "Set the outer limits of the color scale", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("size_ptm"), "Set the height of the plot region", "top", options = list(container = "body")),
                                                   # shinyBS::bsTooltip("if_rowname_color", "Whether set the color of the row names,(purple: Peptides == 1, blue: Peptides == 2)", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("if_chooseToshow_ptm"), "Whether only show the rows that you want to show", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("chooseToshow_ptm"), "Choose the rows that you want to show", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("row_font_size_ptm"), "Set the size of row labels", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("col_font_size_ptm"), "Set the size of column labels", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Height_ptm"), "Height of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("cluster_columns_ptm"), "Whether make cluster on columns", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Custom_columns_order_ptm"), "Order of column, act when [Cluster columns] is FALSE", "top", options = list(container = "body")),
                                                   # shinyBS::bsTooltip(ns("if_mysplit_ptm"), "Whether order the splices manually", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("mysplit_ptm"), "The ordered splices of the heatmap (from top to bottom), only act when both [If mysplit] and [Kmeans] are TRUE", "top", options = list(container = "body"))

                                          ),
                                          tabPanel(title = "Custom Volcano",
                                                   fluidRow(
                                                     column(width = 4,
                                                            colourpicker::colourInput(inputId = ns("up_color_ptm"),
                                                                                      label = "up",
                                                                                      showColour = "both",
                                                                                      palette = "square",
                                                                                      value = "#B2182B")
                                                     ),
                                                     column(width = 4,
                                                            colourpicker::colourInput(inputId = ns("down_color_ptm"),
                                                                                      label = "down",
                                                                                      showColour = "both",
                                                                                      palette = "square",
                                                                                      value = "#2166AC")
                                                     ),
                                                     column(width = 4,
                                                            colourpicker::colourInput(inputId = ns("stable_color_ptm"),
                                                                                      label = "not significant",
                                                                                      showColour = "both",
                                                                                      palette = "square",
                                                                                      value = "#BEBEBE")
                                                     )
                                                   ),
                                                   fluidRow(
                                                     box(uiOutput(ns("volcano_cntrst_ptm")), width = 8),
                                                     box(selectizeInput(ns("labelWay_ptm"),
                                                                        "label way",
                                                                        choices = c("all", "up","down"),
                                                                        selected = c("all"), multiple = FALSE),width = 4)
                                                   ),
                                                   fluidRow(
                                                     box(
                                                       # checkboxInput(ns("if_label_rectangle_ptm"),
                                                       #               "label with rectangle",
                                                       #               value = FALSE),

                                                       # conditionalPanel(condition = paste0("input['",ns("peptidescol"),"']"," !== ''"), #"input.peptidescol !== ''",#&& input.peptidescol !== '' typeof input.peptidescol !== null
                                                       #                  checkboxInput(ns("if_peptide_color"),
                                                       #                                "peptide color",
                                                       #                                value = FALSE),
                                                       #                  # column(width = 6,
                                                       #                  uiOutput(ns("Peptides1")),
                                                       #                  uiOutput(ns("Peptides2"))
                                                       # ),
                                                       # checkboxInput("if_peptide_color",
                                                       #               "peptide color",
                                                       #               value = FALSE),
                                                       # # column(width = 6,
                                                       # uiOutput("Peptides1"),
                                                       # uiOutput("Peptides2"),
                                                       # ),
                                                       # column(width = 6, uiOutput("Peptides2")),
                                                       # numericInput(ns("stroke_ptm"),
                                                       #              "point outside width",
                                                       #              min = 0, max = 5, value = 0.6),
                                                       uiOutput(ns("selected_proteins_ptm")),
                                                       numericInput(ns("showNum_ptm"),
                                                                    "show number",
                                                                    min = 0, max = 100000, value = 20),
                                                       numericInput(ns("fontSize_ptm"),
                                                                    "font size",
                                                                    min = 0, max = 20, value = 5),
                                                       numericInput(ns("dotsize_ptm"),
                                                                    "dot size",
                                                                    min = 0, max = 20, value = 2),
                                                       numericInput(ns("custom_volcano_Width_ptm"),
                                                                    "width",
                                                                    min = 1, max = 30, value = 7),
                                                       numericInput(ns("custom_volcano_Height_ptm"),
                                                                    "height",
                                                                    min = 1, max = 30, value = 7), width = 8, collapsible = TRUE, collapsed = TRUE),
                                                     conditionalPanel(condition = paste0("input['",ns("threshold_method_ptm"),"']"," == 'intersect'"), #"input.threshold_method == 'intersect'",
                                                                      box(checkboxInput(ns("P_adj_ptm"),
                                                                                        "Adjusted p values",
                                                                                        value = FALSE),
                                                                          checkboxInput(ns("Same_width_ptm"),
                                                                                        "Same width",
                                                                                        value = FALSE),
                                                                          width = 4)
                                                     ),
                                                     conditionalPanel(condition = paste0("input['",ns("threshold_method_ptm"),"']"," == 'curve'"), #"input.threshold_method == 'curve'",
                                                                      box(checkboxInput(ns("Same_width_ptm"),
                                                                                        "Same width",
                                                                                        value = FALSE),
                                                                          width = 4)
                                                     )
                                                     # box(checkboxInput("P_adj",
                                                     #                   "Adjusted p values",
                                                     #                   value = FALSE),
                                                     #     checkboxInput("Same_width",
                                                     #                   "Same width",
                                                     #                   value = FALSE),
                                                     #     width = 4)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("custom_volcano_ptm"), height = 600),
                                                     downloadButton(ns('download_custom_volcano_ptm'), 'Save custom volcano')
                                                   ),
                                                   shinyBS::bsTooltip(ns("up_color_ptm"), "Set the color of the points those are up-regulated", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("down_color_ptm"), "Set the color of the points those are down-regulated", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("stable_color_ptm"), "Set the color of the points those are not significant", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("volcano_cntrst_ptm"), "Choose the contrast to plot", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("labelWay_ptm"), "Choose the way of adding labels, one of all significant, up, down, and selected peptides", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("if_label_rectangle_ptm"), "Whether add rectangle underneath the text, making it easier to read", "top", options = list(container = "body")),
                                                   # shinyBS::bsTooltip("if_peptide_color_ptm", "Whether set the color of the points,(default, purple: Peptides == 1, blue: Peptides == 2)", "top", options = list(container = "body")),
                                                   # shinyBS::bsTooltip("Peptides1_ptm", "Set the color of the points which the number of Peptides == 1", "top", options = list(container = "body")),
                                                   # shinyBS::bsTooltip("Peptides2_ptm", "Set the color of the points which the number of Peptides == 2", "top", options = list(container = "body")),
                                                   # shinyBS::bsTooltip(ns("stroke_ptm"), "Set the thickness of black line around the point", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("selected_proteins_ptm"), "Choose the point labels to show, act when [label Way] is selected peptides", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("showNum_ptm"), "Set the number of the labels to add, act when [label way] is significant, up or down. Note that, when it is larger than the number of significant result, it is seted to the number of significant result", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("fontSize_ptm"), "Set the font size of the label", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("dotsize_ptm"), "Set the dot size of the label", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("custom_volcano_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("custom_volcano_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("P_adj_ptm"), "Whether or not to use adjusted p values", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Same_width_ptm"), "Whether the x axis to have the same width from 0", "top", options = list(container = "body"))
                                          )
                              )),
                       tabBox(title = "QC Plots", width = 12,
                              tabsetPanel(id = ns("DEP_QC_tabs_ptm"),
                                          tabPanel(title = "Pca plot",
                                                   fluidRow(
                                                     box(selectizeInput(ns("Indicate_ptm"),
                                                                        "Color and shape",
                                                                        choices = c("condition", "replicate", "Condition", "Replicate"),
                                                                        selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                                     box(checkboxInput(ns("if_square_ptm"),
                                                                       "if square",
                                                                       value = FALSE),
                                                         width = 3)),
                                                   fluidRow(
                                                     box(uiOutput(ns("pca_top_n_ptm")),width = 9)
                                                   ),
                                                   fluidRow(
                                                     box(numericInput(ns("pca_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7), width = 6),
                                                     box(numericInput(ns("pca_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 6)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("pca_ptm"), height = 600),
                                                     downloadButton(ns('downloadPca_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("Indicate_ptm"), "Set the color, shape and facet_wrap of the plot", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("if_square_ptm"), "Whether x limit is equal to y limit", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("pca_top_n_ptm"), "Set the number of top variable peptides to consider", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("pca_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("pca_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "UMAP plot",
                                                   fluidRow(
                                                     box(selectizeInput(ns("umap_Indicate_ptm"),
                                                                        "Color and shape",
                                                                        choices = c("condition", "replicate", "Condition", "Replicate"),
                                                                        selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                                     box(checkboxInput(ns("umap_if_square_ptm"),
                                                                       "if square",
                                                                       value = FALSE),
                                                         width = 6)),
                                                   fluidRow(
                                                     box(uiOutput(ns("umap_n_neighbors_ptm")),width = 10)
                                                   ),
                                                   fluidRow(
                                                     box(numericInput(ns("umap_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7), width = 6),
                                                     box(numericInput(ns("umap_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 6)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("umap_ptm"), height = 600),
                                                     downloadButton(ns('downloadUMAP_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("umap_Indicate_ptm"), "Set the color, shape and facet_wrap of the plot", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("umap_if_square_ptm"), "Whether x limit is equal to y limit", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("umap_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("umap_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("umap_n_neighbors_ptm"), "number of nearest neighbors for umap", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "t-SNE plot",
                                                   fluidRow(
                                                     box(selectizeInput(ns("Tsne_Indicate_ptm"),
                                                                        "Color and shape",
                                                                        choices = c("condition", "replicate", "Condition", "Replicate"),
                                                                        selected = c("condition", "replicate"), multiple = TRUE), width = 6),
                                                     box(checkboxInput(ns("Tsne_if_square_ptm"),
                                                                       "if square",
                                                                       value = FALSE),
                                                         width = 6)
                                                   ),
                                                   fluidRow(
                                                     box(uiOutput(ns("Tsne_perplexity_ptm")),width = 8),
                                                     box(numericInput(ns("Tsne_theseed_ptm"),
                                                                      "set seed for t-SNE",
                                                                      value = 42),
                                                         width = 4)
                                                   ),
                                                   fluidRow(
                                                     box(numericInput(ns("Tsne_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7), width = 6),
                                                     box(numericInput(ns("Tsne_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 6)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("Tsne_ptm"), height = 600),
                                                     downloadButton(ns('downloadTSNE_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("Tsne_Indicate_ptm"), "Set the color, shape and facet_wrap of the plot", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Tsne_if_square_ptm"), "Whether x limit is equal to y limit", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Tsne_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Tsne_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Tsne_perplexity_ptm"), "perplexity for t-SNE", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "Statistical plot",
                                                   # fluidRow(
                                                   #   box(uiOutput(ns("Volcano_cntrst_stastical_plot_ptm")), width = 8),
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
                                                     plotOutput(ns("stastical_plot_ptm"),height=600),
                                                     downloadButton(ns('downloadstastical_plot_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("stastical_plot_x"), "Select the statistic to plot", "top", options = list(container = "body")),
                                                   # shinyBS::bsTooltip(ns("stastical_plot_y"), "Set the vector that y axis represents", "top", options = list(container = "body")),
                                                   # shinyBS::bsTooltip(ns("stastical_plot_type"), "Set the plot type: x-y represents that X-axis represents your selected x, and y-axis represents your selected y; histogram represents the histogram plot of your selected x", "top", options = list(container = "body")),
                                                   # shinyBS::bsTooltip(ns("Volcano_cntrst_stastical_plot_ptm"), "Choose the contrast to plot", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("stastical_plot_Width"), "Width of the figure to export when save", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("stastical_plot_Height"), "Height of the figure to export when save", "top", options = list(container = "body"))
                                          ),
                                          #tabPanel(title = "Pca plot",
                                          #         plotOutput("pca", height = 600),
                                          #         downloadButton('downloadPca', 'Save')
                                          #)
                                          tabPanel(title = "meanSdPlot",
                                                   fluidRow(
                                                     box(numericInput(ns("meanSdPlot_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7), width = 6),
                                                     box(numericInput(ns("meanSdPlot_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 6)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("meanSdPlot_ptm"), height = 600),
                                                     downloadButton(ns('downloadmeanSdPlot_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("meanSdPlot_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("meanSdPlot_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "Pearson correlation",
                                                   fluidRow(
                                                     box(selectizeInput(ns("Pearson_pal_ptm"),
                                                                        "color panel",
                                                                        choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Blues",  "BuGn", "BuPu", "GnBu", "Greens" , "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
                                                                        selected = c("PRGn"), multiple = FALSE),
                                                         checkboxInput(ns("Pearson_pal_rev_ptm"),
                                                                       "pal rev",
                                                                       value = FALSE), width = 4),
                                                     box(numericInput(ns("Pearson_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7),
                                                         numericInput(ns("Pearson_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 4),
                                                     box(numericInput(ns("Pearson_lower_ptm"),
                                                                      "lower",
                                                                      min = -1, max = 1, value = -1),
                                                         numericInput(ns("Pearson_upper_ptm"),
                                                                      "upper",
                                                                      min = -1, max = 1, value = 1), width = 4)
                                                   ),
                                                   fluidRow(
                                                     box(checkboxInput(ns("add_values_for_DEP_person_ptm"),
                                                                       "Add values",
                                                                       value = FALSE),
                                                         numericInput(ns("value_size_for_DEP_person_ptm"),
                                                                      "Value size",
                                                                      min = 1, max = 30, value = 10),
                                                         numericInput(ns("value_digits_for_DEP_person_ptm"),
                                                                      "Value digits",
                                                                      min = 1, max = 30, value = 2), width = 12, collapsible = TRUE, collapsed = TRUE)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("Pearson_correlation_ptm"), height = 600),
                                                     downloadButton(ns('download_Pearson_correlation_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("Pearson_pal_ptm"), "Set the color panel (from RColorBrewer)", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Pearson_pal_rev_ptm"), "Whether or not to invert the color palette", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Pearson_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Pearson_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Pearson_lower_ptm"), "Set the lower limit of the color scale", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Pearson_upper_ptm"), "Set the upper limit of the color scale", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "Gower's distance",
                                                   fluidRow(
                                                     box(selectizeInput(ns("Gower_pal_ptm"),
                                                                        "color panel",
                                                                        choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Blues",  "BuGn", "BuPu", "GnBu", "Greens" , "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"),
                                                                        selected = c("YlOrRd"), multiple = FALSE),
                                                         checkboxInput(ns("Gower_pal_rev_ptm"),
                                                                       "pal rev",
                                                                       value = TRUE), width = 6),
                                                     box(numericInput(ns("Gower_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7),
                                                         numericInput(ns("Gower_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 6)
                                                   ),
                                                   fluidRow(
                                                     box(checkboxInput(ns("add_values_for_DEP_gower_ptm"),
                                                                       "Add values",
                                                                       value = FALSE),
                                                         numericInput(ns("value_size_for_DEP_gower_ptm"),
                                                                      "Value size",
                                                                      min = 1, max = 30, value = 10),
                                                         numericInput(ns("value_digits_for_DEP_gower_ptm"),
                                                                      "Value digits",
                                                                      min = 1, max = 30, value = 2), width = 12, collapsible = TRUE, collapsed = TRUE)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("Gowers_distance_ptm"), height = 600),
                                                     downloadButton(ns('download_Gowers_distance_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("Gower_pal_ptm"), "Set the color panel (from RColorBrewer)", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Gower_pal_rev_ptm"), "Whether or not to invert the color palette", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Gower_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Gower_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "Sample CVs",
                                                   fluidRow(
                                                     box(numericInput(ns("Sample_CVs_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7),
                                                         width = 6),
                                                     box(
                                                       numericInput(ns("Sample_CVs_Height_ptm"),
                                                                    "height",
                                                                    min = 1, max = 30, value = 7), width = 6
                                                     )
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("Sample_CVs_ptm"), height = 600),
                                                     downloadButton(ns('download_Sample_CVs_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("Sample_CVs_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Sample_CVs_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "Peptide Numbers",
                                                   fluidRow(
                                                     box(numericInput(ns("numbers_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7), width = 6),
                                                     box(numericInput(ns("numbers_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 6)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("numbers_ptm"), height = 600),
                                                     downloadButton(ns('downloadNumbers_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("numbers_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("numbers_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "Sample coverage",
                                                   fluidRow(
                                                     box(numericInput(ns("coverage_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7), width = 6),
                                                     box(numericInput(ns("coverage_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 6)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("coverage_ptm"), height = 600),
                                                     downloadButton(ns('downloadCoverage_ptm'), 'Save')

                                                   ),
                                                   shinyBS::bsTooltip(ns("coverage_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("coverage_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "Normalization",
                                                   fluidRow(
                                                     box(numericInput(ns("norm_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7), width = 6),
                                                     box(numericInput(ns("norm_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 6)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("norm_ptm"), height = 600),
                                                     downloadButton(ns('downloadNorm_ptm'), 'Save')

                                                   ),
                                                   shinyBS::bsTooltip(ns("norm_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("norm_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "Missing values - Quant",
                                                   fluidRow(
                                                     box(numericInput(ns("detect_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7), width = 6),
                                                     box(numericInput(ns("detect_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 6)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("detect_ptm"), height = 600),
                                                     downloadButton(ns('downloadDetect_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("detect_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("detect_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "Missing values - Heatmap",
                                                   fluidRow(
                                                     box(
                                                       checkboxInput(ns("cluster_columns_for_missval_heatmap_ptm"),
                                                                     "Cluster columns",
                                                                     value = TRUE), width = 4),
                                                     box(uiOutput(ns("Custom_columns_order_for_missval_heatmap_ptm")),
                                                         width = 8)
                                                   ),
                                                   fluidRow(
                                                     box(numericInput(ns("missval_heatmap_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7), width = 6),
                                                     box(numericInput(ns("missval_heatmap_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 10), width = 6)),
                                                   fluidRow(
                                                     plotOutput(ns("missval_ptm"), height = 600),
                                                     downloadButton(ns('downloadMissval_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("cluster_columns_for_missval_heatmap_ptm"), "Whether make cluster on columns", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("Custom_columns_order_for_missval_heatmap_ptm"), "Order of column, act when [Cluster columns] is FALSE", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("missval_heatmap_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("missval_heatmap_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                                          ),
                                          tabPanel(title = "Imputation",
                                                   fluidRow(
                                                     box(numericInput(ns("imputation_Width_ptm"),
                                                                      "width",
                                                                      min = 1, max = 30, value = 7), width = 6),
                                                     box(numericInput(ns("imputation_Height_ptm"),
                                                                      "height",
                                                                      min = 1, max = 30, value = 7), width = 6)
                                                   ),
                                                   fluidRow(
                                                     plotOutput(ns("imputation_ptm"), height = 600),
                                                     downloadButton(ns('downloadImputation_ptm'), 'Save')
                                                   ),
                                                   shinyBS::bsTooltip(ns("imputation_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                                                   shinyBS::bsTooltip(ns("imputation_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                                          )

                              ))
                )
              ))
    # ,change_reactive_ID(id)
  )

}

DEP_ptm_server_module <- function(id, Omics_res){
  moduleServer(
    id,
    function(input,output,session){

      Omics_res_list <- reactive({Omics_res_list_save <<- reactiveValuesToList(Omics_res)})

      # observeEvent(input$threshold_method,{

        # req(!is.null(input$threshold_method))
        cat("DEP_ptm_server_module load")
        thedf <- reactive(NULL)

        ### UI functions ### --------------------------------------------------------
        output$name_ptm <- renderUI({
          shiny::validate(need(!is.null(data_ptm()), " "))
          selectizeInput(session$ns("name_ptm"),
                         "Gene name column",
                         choices=colnames(data_ptm()),
                         selected = "Gene.names")
        })

        output$id_ptm <- renderUI({
          shiny::validate(need(!is.null(data_ptm()), " "))
          selectizeInput(session$ns("id_ptm"),
                         "Protin ID column",
                         choices=colnames(data_ptm()),
                         selected = "Protein")
        })

        output$amino_acid_ptm <- renderUI({
          selectizeInput(session$ns("amino_acid_ptm"),
                         "Amino acid column",
                         choices=colnames(data_ptm()),
                         selected = "Amino.acid")
        })

        output$aa_position_ptm <- renderUI({
          selectizeInput(session$ns("aa_position_ptm"),
                         "Amino acid position column",
                         choices=colnames(data_ptm()),
                         selected = "Position")
        })

       # output$peptidescol_ptm <- renderUI({ input$peptidescol_input_ptm})
        # peptidescol_input_ptm
        # output$peptidescol_ptm <- renderUI({
        #   if("Peptides" %in% colnames(data_ptm())){
        #     return(selectizeInput("peptidescol_ptm",
        #                    "peptide column",
        #                    choices=c(colnames(data_ptm()), ""),
        #                    selected = "Peptides"))
        #   }else{
        #     return(selectizeInput("peptidescol_ptm",
        #                    "peptide column",
        #                    choices=c("", colnames(data_ptm())),
        #                    selected = NULL))
        #   }
        # })
        output$intensitycols_ptm <- renderUI({
          # LFQcols = grep("^LFQ",colnames(data_ptm()),value = T)
          # if(length(LFQcols)==0) grep("^Intensity",colnames(data_ptm()),value = T)
          # if(length(LFQcols)==0) grep("^intensity",colnames(data_ptm()),value = T)
          # if(length(LFQcols)==0) grep("^intensity",colnames(data_ptm()),value = T)
          # selectizeInput("intensitycols_ptm",
          #                " ",#Expression columns
          #                choices=colnames(data_ptm()),
          #                multiple = TRUE,
          #                selected = grep("^LFQ",LFQcols,value = T), width = '100%')
          shinydashboardPlus::box(selectizeInput(session$ns("intensitycols_ptm"),
                                                 label = "Choose the expression columns",
                                                 choices=colnames(data_ptm()),
                                                 multiple = TRUE,
                                                 selected = colnames(data_ptm())[grep("^Intensity. *",colnames(data_ptm()),value = F) %>% .[!.%in%grep("__",colnames(data_ptm()),value = F)]], width = '100%'), title = "Expression columns", width = 12, icon = icon("th"), collapsible = TRUE, collapsed = TRUE, background = "light-blue")
        })



        output$filt_ptm <- renderUI({
          selectizeInput(session$ns("filt_ptm"),
                         "Filter on columns" ,
                         colnames(data_ptm()),
                         multiple = TRUE,
                         selected = c("Reverse","Potential.contaminant"))
        })

        output$filt_num_cutoff_ptm <- renderUI({
          selectizeInput(session$ns("filt_num_cutoff_ptm"),
                         "Cutoff based on" ,
                         colnames(data_ptm()),
                         multiple = FALSE,
                         selected = c("Localization.prob"))
        })

        output$order <- renderUI({
          validate(need(!is.null(input$file1_ptm), "Please upload expression ProteinGroups in Files"))
          validate(need(!is.null(input$intensitycols_ptm), "Please select the Expression columns"))
          validate(need(length(input$intensitycols_ptm) > 1 , "More expression columns is required"))
          # groups = exp_design()[,input$groupby]
          if (input$anno_ptm == "columns" & !is.null(data_ptm()) & input$contrasts_ptm == "control") {
            print("aaaa")
            my_data <- data_ptm()
            cols <- which(colnames(data_ptm()) %in% input$intensitycols_ptm)#according to intensitycols
            # if(input$remove_prefix){
              prefix <- DEP2::get_prefix(data_ptm()[,cols] %>% colnames())
              label = colnames(data_ptm())[cols] %>% gsub(prefix,"",.)
            # }else label = colnames(data_ptm())[cols]
            # groups
            condition = make.names(unlist(lapply(label %>% strsplit(., split = "_"), function(x){x[1]})))
          } else {
            print("bbb")
            if (input$anno_ptm == "expdesign" & !is.null(expdesign_ptm()) & input$contrasts_ptm == "control") {
              condition = make.names(expdesign_ptm()$condition)
            }
          }

          selectizeInput(session$ns("order"),
                         "Set the groups order",
                         choices = condition, multiple = TRUE,
                         selected = NULL, width = '100%')
        })

        output$control_ptm <- renderUI({
          validate(need(!is.null(input$file1_ptm), ""))
          if (input$anno_ptm == "columns" & !is.null(data_ptm()) & input$contrasts_ptm == "control") {
            my_data <- data_ptm()
            # cols <- grep("^LFQ", colnames(data_ptm()))
            cols <- which(colnames(data_ptm()) %in% input$intensitycols_ptm)#according to intensitycols
            temp <- data_ptm()[,cols] %>% colnames()
            prefix <- DEP2::get_prefix(temp)
            selectizeInput(session$ns("control_ptm"), "Control",
                           choices=make.names(unlist(lapply(colnames(data_ptm())[cols] %>% gsub(prefix,"",.) %>% strsplit(., split = "_"), function(x){x[1]}))), selected = NULL)
          } else {
            if (input$anno_ptm == "expdesign" & !is.null(expdesign_ptm()) & input$contrasts_ptm == "control") {
              selectizeInput(session$ns("control_ptm"),
                             "Control",
                             choices = make.names(expdesign_ptm()$condition))
            }
          }
        })


        output$test_manual_ptm <- renderUI({
          validate(need(!is.null(input$file1_ptm), ""))
          if(!is.null(data_ptm()) & input$contrasts_ptm == "manual"){
            # cols <<- grep("^LFQ", colnames(data_ptm()))
            cols <- which(colnames(data_ptm()) %in% input$intensitycols_ptm)#according to intensitycols
            prefix <<- get_prefix(data_ptm()[,cols] %>% colnames())
            test_manual_name <- unique(make.names(unlist(lapply(colnames(data_ptm())[cols] %>% gsub(prefix,"",.) %>% strsplit(., split = "_"), function(x){x[1]}))))
            test_manual_name <- cbind(combn(test_manual_name,2),combn(test_manual_name,2, FUN = rev))
            test_manual_name <- apply(test_manual_name, 2, function(i){paste(i[1], i[2], sep = "_vs_")})
            selectizeInput(session$ns("test_manual_ptm"), "Manual test",
                           choices=test_manual_name, selected = NULL, multiple = TRUE)
          }
        })

        output$protein_group_choose <- renderUI({
          Omics_res_list = Omics_res_list()
          Omics_res_list2 = Omics_res_list[grep("^Proteome",names(Omics_res_list))]

          if(length(Omics_res_list2) == 0) {
            Omics_res_list2 = NULL
          }else{
            Omics_res_list2 = Omics_res_list2[
              which((Omics_res_list2 %>%
                       lapply(., function(x){
                         x1 <<- x
                         temp = try(x(), silent = T)
                         if(class(temp) == "try-error")
                           return(FALSE)
                         return(TRUE)
                         # is.na(x1 <<- x)
                       }) %>%
                       unlist()))
            ]}

          if(is.null(Omics_res_list2) ||length(Omics_res_list2) < 1)
            return(p("Required a exit proteomics analysis. Please perform a proteomics analysis first!"))
          conditionalPanel(
            condition = paste0("input['", session$ns("Protein_level_correction"),"']"),
            selectInput(inputId = session$ns("protein_group_choose"),
                        label = "Corresponding proteomics", choices = names(Omics_res_list2))
          )
        })

        output$Protein_level_correction_options <- renderUI({
          options <- tagList(
            radioButtons(
              inputId = session$ns("correct_key"),
              label = "The column correction based on",
              choices = list("Gene name" = "name", "Protein ID" = "ID"),
              inline = TRUE
              # status = "success",
              # fill = TRUE
            ),
            radioButtons(
              inputId = session$ns("unidentified_protein_treatment"),
              label = "How to treat peptide whose proteins unidentified in corresponding proteomics",
              choices = c("retain", "remove"),
              inline = TRUE,
              # status = "success",
              # fill = TRUE
            ),
            radioButtons(
              inputId = session$ns("correct_level"),
              label = "Corresponding relationship between",
              choices = list("condition to condition" = "condition", "repetition to repetition" = "label"),
              inline = TRUE,
              # status = "success",
              # fill = TRUE
            ),
            shinyBS::bsTooltip(session$ns("correct_key"),
                               "Match enriched-peptides quantity and corresponding proteomics based on Gene name or protein ID (both were set in the Columns options)",
                               "top", options = list(container = "body")),
            shinyBS::bsTooltip(session$ns("unidentified_protein_treatment"),
                               "Some enriched peptides' proteins may lack effective quantity in corresponding proteomics. /n
                               Choose remove to eliminate these peptides or choose retain to transmit the original peptide quantity(without correction)",
                               "top", options = list(container = "body"))

          )
          conditionalPanel(
            condition = paste0("input['", session$ns("Protein_level_correction"),"']"),
            options)

        })


      # output$Peptides1_ptm <- renderUI({
      #   if(input$if_peptide_color_ptm) {
      #     colourInput(inputId = "Peptides_1_ptm",
      #                 label = "Peptides1",
      #                 showColour = "both",
      #                 palette = "square",
      #                 value = "#A020F0")
      #   }
      # })

      # output$Peptides2_ptm <- renderUI({
      #   if(input$if_peptide_color_ptm) {
      #     colourpicker::colourInput(inputId = "Peptides_2_ptm",
      #                               label = "Peptides2",
      #                               showColour = "both",
      #                               palette = "square",
      #                               value = "#0000FF")#
      #   }
      # })




      # observe(
      #   if(input$threshold_method_ptm == "intersect") {
      #     # removeTab(inputId = "DEP_ptm_results_tabs", target = "DEP curve 1")
      #     insertTab(inputId = "DEP_ptm_results_tabs",
      #               tabPanel(title = "Volcano plot",
      #                        fluidRow(
      #                          box(uiOutput(session$ns("volcano_cntrst_ptm")), width = 9),
      #                          box(numericInput(session$ns("fontsize_ptm"),
      #                                           "Font size",
      #                                           min = 0, max = 8, value = 4),
      #
      #                              width = 3)
      #                        ),
      #                        fluidRow(
      #                          box(checkboxInput(session$ns("my_breaks_ptm"),
      #                                            "My breaks",
      #                                            value = FALSE),
      #                              multiInput(session$ns("mybreaks_ptm"),"Mybreaks",choices = seq(-40, 40, 0.5)),
      #                              numericInput(session$ns("Volcano_Width_ptm"),
      #                                           "width",
      #                                           min = 1, max = 30, value = 7),
      #                              numericInput(session$ns("Volcano_Height_ptm"),
      #                                           "height",
      #                                           min = 1, max = 30, value = 7),
      #                              #selectizeInput("mybreaks_ptm","Mybreaks",choices = seq(-40, 40, 1), multiple = TRUE, size = 10),
      #                              width = 8, collapsible = TRUE, collapsed = TRUE),
      #                          box(checkboxInput(session$ns("check_names_ptm"),
      #                                            "Display names",
      #                                            value = FALSE),
      #                              checkboxInput(session$ns("p_adj_ptm"),
      #                                            "Adjusted p values",
      #                                            value = FALSE),
      #                              checkboxInput(session$ns("same_width_ptm"),
      #                                            "Same width",
      #                                            value = FALSE),
      #                              width = 4)),
      #                        # fluidRow(
      #                        #   plotOutput(session$ns("volcano_ptm"), height = 600),
      #                        #   downloadButton(session$ns('downloadVolcano_ptm'), 'Save volcano')
      #                        # ),
      #                        shinyBS::bsTooltip(session$ns("volcano_cntrst_ptm"), "Choose the contrast to plot", "top", options = list(container = "body")),
      #                        shinyBS::bsTooltip(session$ns("fontsize_ptm"), "Set the size of name labels", "top", options = list(container = "body")),
      #                        shinyBS::bsTooltip(session$ns("my_breaks_ptm"), "Whether set the breaks of x axis by yourself", "top", options = list(container = "body")),
      #                        shinyBS::bsTooltip(session$ns("mybreaks_ptm"), "Choose the breaks of x axis, act when [My breaks] is TRUE", "top", options = list(container = "body")),
      #                        shinyBS::bsTooltip(session$ns("Volcano_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
      #                        shinyBS::bsTooltip(session$ns("Volcano_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body")),
      #                        shinyBS::bsTooltip(session$ns("check_names_ptm"), "Whether or not to plot names", "top", options = list(container = "body")),
      #                        shinyBS::bsTooltip(session$ns("p_adj_ptm"), "Whether or not to use adjusted p values", "top", options = list(container = "body")),
      #                        shinyBS::bsTooltip(session$ns("same_width_ptm"), "Whether the x axis to have the same width from 0", "top", options = list(container = "body"))
      #               ),
      #               target = "Heatmap", position = "after"
      #     )
      #
      #   }
      # )



        observe(
          if(input$threshold_method_ptm == "curve") {
            removeTab(inputId = "DEP_ptm_results_tabs", target = "Volcano plot")
            # removeTab(inputId = "DEP_ptm_results_tabs", target = "Custom Volcano")
            # insertTab(inputId = "DEP_ptm_results_tabs",
            #           tabPanel("DEP curve 1", "This a DEP curve method added tab"),
            #           target = "Heatmap", position = "after"
            # )
          }
        )


        observe(
          if(input$threshold_method_ptm == "curve") {
            insertTab(inputId = "DEP_QC_tabs_ptm",
                      tabPanel(title = "Normal distribution - Fit",
                               fluidRow(
                                 box(uiOutput(session$ns("cntrst_norm_distribution_ptm")), width = 8)
                               ),
                               fluidRow(
                                 box(numericInput(session$ns("norm_distribution_Width_ptm"),
                                                  "width",
                                                  min = 1, max = 30, value = 7), width = 6),
                                 box(numericInput(session$ns("norm_distribution_Height_ptm"),
                                                  "height",
                                                  min = 1, max = 30, value = 7), width = 6)
                               ),
                               fluidRow(
                                 plotOutput(session$ns("norm_distribution_plot_ptm"),height=600),
                                 downloadButton(session$ns('download_norm_distribution_plot_ptm'), 'Save')
                               ),
                               shinyBS::bsTooltip(session$ns("cntrst_norm_distribution_ptm"), "Choose the contrast to plot", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(session$ns("norm_distribution_Width_ptm"), "Width of the figure to export", "top", options = list(container = "body")),
                               shinyBS::bsTooltip(session$ns("norm_distribution_Height_ptm"), "Height of the figure to export", "top", options = list(container = "body"))
                      ),

                      target = "Pca plot", position = "before"
            )
          }
        )

        observe(
         if(input$threshold_method_ptm == "intersect") {
           removeTab(inputId = "DEP_QC_tabs_ptm", target = "DEP curve 2")
         }
       )

        ### check input ### --------------------------------------------------------
        iv <- InputValidator$new()
        iv$add_rule("p_ptm", sv_between(0.0001, 0.1))
        iv$add_rule("lfc_ptm", sv_between(0, 10))
        iv$add_rule("curvature_ptm", sv_between(0.1, 6))
        iv$add_rule("x0_fold_ptm", sv_between(1, 10))
        iv$add_rule("filt_num_ptm", sv_between(0, 1))
        iv$enable()

        ### Reactive functions ### --------------------------------------------------
        expdesign_ptm <- reactive({
          inFile <- input$file2_ptm
          if (is.null(inFile))
            return(NULL)
          read.csv(inFile$datapath, header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE) %>%
            mutate(id = row_number())
        })

        data_ptm <- reactive({
          inFile <- input$file1_ptm
          if (is.null(inFile))
            return(NULL)
          my_data <<- read.csv(inFile$datapath, header = TRUE,
                               sep = "\t", stringsAsFactors = FALSE) %>%
            mutate(id = row_number())
        })

      #   data_ptm <- reactive({
      #     inFile <- input$file1_ptm
      #     if (is.null(inFile)){
      #       cat("infile1 is null /n")
      #       return(NULL)
      #     }
      #     cat("infile1 is not null /n")
      #     my_data <- data.table::fread(inFile$datapath, header = TRUE,
      #                                  sep = "\t", stringsAsFactors = FALSE, integer64 = "numeric") %>%
      #       mutate(id = row_number())
      # coln.int64 <- names(which(sapply(my_data, bit64::is.integer64)))
      # if (length(coln.int64) > 0L)
      #   my_data[, c(coln.int64) := lapply(.SD, as.numeric), .SDcols = coln.int64]
      #     my_data = as.data.frame(my_data)
      #     colnames(my_data) = gsub("^\\[(.*)\\] ","",colnames(my_data))
      #     colnames(my_data) = make.names(colnames(my_data))
      #     my_data <- my_data
      #   })


        filt0_ptm <- reactive({
          data <- data_ptm()
          # aaaa <<- data
          # order_save <<- the_order()
          # if(is.null(input$peptidescol_ptm)|input$peptidescol_ptm==""){
          #   data$Peptides = 3
          # }else{
          #   data$Peptides = data[,input$peptidescol_ptm]
          # }

          cols <- which(colnames(data) %in% input$intensitycols_ptm)
          data[,cols] = apply(data[,cols], 2, function(x){
              x[!(!grepl("[A-z]",x) & grepl("\\d",x))] = 0
              return(as.numeric(x))
            })


          if(is.null(input$filt_ptm)){
            cat("the filter column is empty! Do not filter with column")
          }

          # req(input$id_ptm, input$name_ptm)
          validate(need(!(input$name_ptm == "" & input$id_ptm == ""), "Please ensure that: at least one of your name column and id column is non-empty!"))
        if(input$name_ptm == "" & input$id_ptm == "") {
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = "Please ensure that: at least one of your name column and id column is non-empty!",
            type = "warning"
          )
        }
          validate(need(input$amino_acid_ptm != "", "Please ensure that: Amino acid column is non-empty"))
          validate(need(input$aa_position_ptm != "", "Please ensure that: Amino acid position column is non-empty"))

          filter_column_names = input$filt_ptm
          cols_filt <- intersect(filter_column_names, colnames(data))
          data_save1 <<- data
          if(length(cols_filt) > 0){
            message("Filtering based on '", paste(cols_filt,
                                                  collapse = "', '"), "' column(s)")

            filtered = data[apply((data[, cols_filt,drop = F] == "" | is.na(data[, cols_filt, drop = F])), 1, all), ]
          }
          filtered_save1 <<- filtered
          filtered = make_unique_ptm(filtered, gene_name = ifelse(input$name_ptm == "", input$id_ptm, input$name_ptm),
                                     protein_ID = ifelse(input$id_ptm == "", input$name_ptm, input$id_ptm),
                                     modified_name = NULL,
                                     aa = input$amino_acid_ptm, pos = input$aa_position_ptm,
                                     delim = input$delim_ptm
                                     )
          filtered_save2 <<- filtered

          # if (!is.null(cols_filt)) {
          #   NAs <- is.na(data[, cols_filt])
          #   data[, cols_filt][NAs] <- ""
          #   if (length(cols_filt) == 1) {
          #     data <- dplyr::filter(data, data[, cols_filt] !=
          #                             "+")
          #   }
          #   else if (length(cols_filt) > 1) {
          #     data <- dplyr::filter(data, !apply(data[, cols_filt] ==
          #                                          "+", 1, any))
          #   }
          # }

          # filtered <- try(filter_peptide(PTMdata = data,
          #                                gene.name = ifelse(input$name_ptm == "", input$id_ptm, input$name_ptm),
          #                                protein.ID = ifelse(input$id_ptm == "", input$name_ptm, input$id_ptm),
          #                                aa = input$amino_acid_ptm,
          #                                pos = input$aa_position_ptm,
          #                                filter_column_names = input$filt_ptm,
          #                                delim = input$delim_ptm,
          #                                cutoff_based = input$filt_num_cutoff_ptm,
          #                                cutoff = input$filt_num_ptm),
          #                 silent = TRUE)
          # validate(need(filtered, "The function filter_peptide get an error"))
          # a133 <<- input$filt
          a122 <<- input$name_ptm
          a133 <<- input$id_ptm

          ind_empty = c(grep("^_", filtered$name), which(filtered$name == ""))
          if(length(ind_empty) > 0) {
            filtered = filtered[-ind_empty, ]
          }

          if (input$anno_ptm == "columns") {
            se <- DEP2::make_se_parse(filtered, cols, mode = "delim", sep = "_")
          }
          if (input$anno_ptm == "expdesign") {
            se <- DEP2::make_se(filtered, cols, expdesign_ptm())
            colData(se)$replicate = as.character(colData(se)$replicate)
          }
          se_save <<- se
          cutoff_col = input$filt_num_cutoff_ptm; cutoff_val = input$filt_num_ptm
          filter_formula = paste( "~(!is.na(",cutoff_col,")&",cutoff_col," >= ",cutoff_val, ")") %>% as.formula()
          filt <- DEP2::filter_se(se_save,thr = input$thr_ptm, filter_formula = filter_formula) ## filter the location prosibility or other scores
          my_filt <<- filt
        })

        iv1 <- InputValidator$new()
        iv1$add_rule("order",
                     sv_required(message = "Complete order is required, sort all groups", test = function(x){
                       is.null(x) || length(x) == (length(filt0_ptm()@colData$condition %>% unique))
                     })
        )
        iv1$enable()

        the_order <- reactive({
          req(iv1$is_valid())
          order_save1 <<- input$order
          if(length(input$order) == length(filt0_ptm()@colData$condition %>% unique))
            return(input$order)
          return(NULL)
        })

        filt_ptm <- reactive({
          filt_ptm = filt0_ptm()
          order_save <<- the_order()
          my_filt_save <<- filt_ptm
          if(is.null(the_order())){
            return(filt_ptm)
          }else{
            order_save2 <<- the_order()
            filt_ptm@colData$condition = factor(filt_ptm@colData$condition, levels = the_order())
            filt_ptm = filt_ptm[,order(filt_ptm@colData$condition)]
            # filt@colData = filt@colData %>% as.data.frame() %>% arrange(., condition) %>% DataFrame()
            return(filt_ptm)
          }

        })

        norm_ptm <- reactive({
           my_norm <<- normalize_vsn(filt_ptm())
        })


        imp_ptm <-  reactive({
          inFile1 <- input$resultRData_ptm

          if(is.null(inFile1)){
            if(input$imputation_ptm == "mixed on proteins") { ## mixed imputation
              # Extract protein names with missing values
              # in all replicates of at least one condition
              proteins_MNAR <- get_df_long(norm_ptm()) %>%
                dplyr::group_by(name, condition) %>%
                dplyr::summarize(NAs = all(is.na(intensity))) %>%
                dplyr::filter(NAs) %>%
                dplyr::pull(name) %>%
                unique()

              # Get a logical vector
              MNAR <- names(norm_ptm()) %in% proteins_MNAR

              # Perform a mixed imputation
              my_imp <- DEP2::impute(
                norm_ptm(),
                fun = "mixed",
                randna = !MNAR, # we have to define MAR which is the opposite of MNAR
                mar = "knn", # imputation function for MAR
                mnar = "MinDet") # imputation function for MNAR
              validate(
                need(length(which(is.na((assay(my_imp)[which(MNAR),])))) == 0, "Please select another imputation method")
              )
              # my_imp

            } else {
              if(input$imputation_ptm == "mixed on samples") { ## perfrom MNAR impute for control and MAR for treatment. May suit for AP-MS experiment
                sample_specific_imputation <- norm_ptm()
                MSnSet <- as(sample_specific_imputation, "MSnSet")

                # Impute differently for two sets of samples
                MSnSet_imputed1 <- MSnbase::impute(MSnSet[, which(gsub("_\\d$", "", colnames(assay(sample_specific_imputation))) %in% input$control_ptm)], method = "MinProb") # control samples that need not random (MNAR) left-censored imputation
                MSnSet_imputed2 <- MSnbase::impute(MSnSet[, -which(gsub("_\\d$", "", colnames(assay(sample_specific_imputation))) %in% input$control_ptm)], method = "knn")# treatment samples that need random(MAR) imputation

                # Combine into the SummarizedExperiment object
                assay(sample_specific_imputation) <- cbind(
                  MSnbase::exprs(MSnSet_imputed1),
                  MSnbase::exprs(MSnSet_imputed2))
                my_imp <- sample_specific_imputation
                validate(
                  need(length(which(is.na((assay(my_imp))))) == 0, "Please select another imputation method")
                )
                # my_imp
              } else {
                #set.seed(12345)
                my_imp <- DEP2::impute(norm_ptm(), input$imputation_ptm)
                # my_imp <<- DEP2::impute(norm_ptm(), input$imputation_ptm)
              }
            }
            my_imp_ptm <<- my_imp
            return(my_imp)
          } else {
            load(file = inFile1$datapath)
            my_imp <- my_imp
            return(my_imp)
          }
        })

        Correct_ptm <- reactive({
          imp_ptm = imp_ptm()
          PG_choosed = input$protein_group_choose
          omictype <- strsplit(PG_choosed,"_")[[1]][1]

          the_res <- Omics_res[[PG_choosed]]()
          # PG_data <- assay(the_res)
          correct_key_save <<- input$correct_key

          Correct_ptm = DEP2::correct_PTM_by_Protein(enriched_peptide = imp_ptm,
                                                     relative_protein = the_res,
                                                     correct_key = input$correct_key,
                                                     correct_level = input$correct_level,
                                                     unidentified_treatment = input$unidentified_protein_treatment)
          Correct_ptm_save <<- Correct_ptm
        })

        FDR_type_ptm <- reactive({
          FDR_type_ptm <- input$FDR_type_ptm
          # FDR_type_ptm <- ifelse(length(grep(FDR_type_ptm,"qval"))>1,"qval","lfdr")
        })

        # df_all_ptm = reactive({
        #   inFile1 <- input$resultRData_ptm
        #   # if(is.null(inFile1)){
        #   #   df_all <- DEP2::test_diff(se = imp_ptm(), type = "all",FDRmethod=input$fdr_correction_ptm)
        #   # }else{
        #   #   load(file = inFile1$datapath)
        #   #   df_all <- DEP2::test_diff(se = my_imp, type = "all",FDRmethod=input$fdr_correction_ptm)
        #   # }
        #
        #   if(is.null(inFile1)){
        #     df_all <- test_diff(se = imp_ptm(), type = "all", fdr.type = FDR_type_ptm())
        #   }else{
        #     load(file = inFile1$datapath)
        #     df_all <- test_diff(se = my_imp, type = "all", fdr.type = FDR_type_ptm())
        #   }
        #   # df_all2 <<- df_all
        #   df_all
        # })

        df_ptm <- reactive({
          if(input$contrasts_ptm == "control") {
            validate(
              need(input$control_ptm != "", "Please select a control condition under menuItem Columns in the DEP-PTM options of the sidebar"),
              need(all(input$intensitycols_ptm %in% colnames(data_ptm())), "Please select the Expression columns in the sidebar")
            )
          }

        if(input$contrasts_ptm == "manual") {
          validate(
            need(input$test_manual_ptm != "", "Please select manual contrasts to test under menuItem Columns"),
            need(all(input$intensitycols_ptm %in% colnames(data_ptm())), "Please select the Expression columns in the sidebar")
          )
        }

        inFile1 <- input$resultRData_ptm
        if(input$Protein_level_correction){
          my_imp = Correct_ptm()
        }else{
          my_imp = imp_ptm()
        }
        # if(is.null(inFile1)){
          if(input$contrasts_ptm == "control"){
            df <- test_diff(se = my_imp, type = input$contrasts_ptm, control = input$control_ptm,
                             fdr.type = FDR_type_ptm())
          }

          if(input$contrasts_ptm == "all") {
            df <- test_diff(se = my_imp, type = "all", fdr.type = FDR_type_ptm())
          }

          if(input$contrasts_ptm == "manual") {
            df <- test_diff(se = my_imp, type = input$contrasts_ptm, test = input$test_manual_ptm,
                            fdr.type = FDR_type_ptm())
          }
        # } else {
          # load your saved RData in order to get the same result (imp is the key of if result from two analysis being the same)
        #   load(file = inFile1$datapath)
        #   my_imp <- my_imp
        #   if(input$contrasts_ptm == "control"){
        #     df <- test_diff(se = my_imp, type = input$contrasts_ptm, control = input$control_ptm,
        #                      fdr.type = FDR_type_ptm())
        #   }
        #
        #   if(input$contrasts_ptm == "all") {
        #     df <- test_diff(se = my_imp, type = "all", fdr.type = FDR_type_ptm())
        #   }
        #
        #   if(input$contrasts_ptm == "manual") {
        #     df <- test_diff(se = my_imp, type = input$contrasts_ptm, test = input$test_manual_ptm,
        #                      fdr.type = FDR_type_ptm())
        #   }
        # }
        # my_df <<- df
        return(df)
      })

        thedf <- df_ptm

        # dep_all_ptm <- reactive({
        #   if(all( input$intensitycols_ptm %in% colnames(data_ptm())) ){
        #     req(input$threshold_method_ptm)
        #     if(input$threshold_method_ptm == "intersect") {
        #       req(input$p_ptm,input$lfc_ptm)
        #       dep_all <- add_rejections(diff = df_all_ptm(), alpha = input$p_ptm, lfc = input$lfc_ptm, thresholdmethod = input$threshold_method_ptm)
        #     }
        #     if(input$threshold_method_ptm == "curve") {
        #       req(input$curvature_ptm,input$x0_fold_ptm)
        #       dep_all <- add_rejections(diff = df_all_ptm(), thresholdmethod = input$threshold_method_ptm, curvature = input$curvature_ptm, x0_fold = input$x0_fold_ptm)
        #     }
        #     # dep_all3 <<- dep_all
        #     req(input$p_ptm,input$lfc_ptm)
        #     # get_dep_ptm_siglist(dep = dep_all,alpha = input$p_ptm , lfc = input$lfc_ptm, threshold_method = input$threshold_method_ptm, curvature = input$curvature_ptm, x0_fold = input$x0_fold_ptm)
        #     # get_dep_listforgsea(dep_all)
        #     # Genelist_server_module("bbb",dep_all = dep_all,alpha = input$p_ptm , lfc = input$lfc_ptm)
        #     return(dep_all)
        #   }
        # })

          dep_ptm <- reactive({
            aaa <<- input$curvature_ptm
            if(input$threshold_method_ptm == "intersect") {
              my_dep <- add_rejections(diff = df_ptm(), alpha = input$p_ptm, lfc = input$lfc_ptm, thresholdmethod = input$threshold_method_ptm)
            }
            if(input$threshold_method_ptm == "curve") {
              my_dep <- add_rejections(diff = df_ptm(), thresholdmethod = input$threshold_method_ptm, curvature = input$curvature_ptm, x0_fold = input$x0_fold_ptm)
            }

            my_dep_ptm_save <<- my_dep
            # if(is.null(input$peptidescol_ptm)|input$peptidescol_ptm==""){
            #     rowData(my_dep)$Peptides = 3
            # }else{
            #    rowData(my_dep)$Peptides = rowData(my_dep)[,input$peptidescol_ptm]
            # }
            return(my_dep)
          })

        output$Save_RData_ptm <- downloadHandler(
          filename = function() { paste("results", ".RData", sep = "") },
          content = function(file) {
            withProgress(message = 'Please wait ...', value = 0.66, {
              save(my_data, my_filt, my_norm, my_imp, my_dep, file=file)})}
        )

        ## All object and functions upon 'Analyze' input  ### ----------------------
         observeEvent(input$analyze_ptm, {
  if(is.null(input$filt_ptm)) {
    sendSweetAlert(
      session = shiny::getDefaultReactiveDomain(),
      title = "Note :",
      text = "Your filter column is empty! If needed, please choose the filter columns. Otherwise, please ignore this message",
      type = "info"
    )
  }

    # validate(
    #   need(input$filt_ptm,"filter column is empty\n")
    # )



    # temp <- observe(dep_all_ptm())
    ### Interactive UI functions ### ------------------------------------------
    output$downloadTable_ptm <- renderUI({
      selectizeInput(session$ns("dataset_ptm"),
                     "Choose a dataset to save" ,
                     c("results","significant_peptides",
                       "displayed_subset","full_dataset"))
    })

    output$downloadButton_ptm <- renderUI({
      downloadButton(session$ns('downloadData_ptm'), 'Save', class = "downloadData_ptm")
    })

    output$downloadButton_for_save_RData_ptm <- renderUI({
      downloadButton(session$ns("Save_RData_ptm"), "save result RData", class = "Save_RData_ptm")
    })

    output$significantBox_ptm <- renderInfoBox({
      num_total <- dep_ptm() %>%
        nrow()
      num_signif <- dep_ptm() %>%
        .[rowData(.)$significant, ] %>%
        nrow()
      frac <- num_signif / num_total

      if(frac > 0.2) {
        info_box <- infoBox("Significant peptides",
                            paste0(num_signif,
                                   " out of ",
                                   num_total),
                            paste0("Too large fraction (",
                                   signif(frac * 100, digits = 3),
                                   "%) of peptides differentially expressed"),
                            icon = icon("minus", lib = "glyphicon"),
                            color = "orange",
                            width = 4)
      }
      if(frac == 0) {
        info_box <- infoBox("Significant peptides",
                            paste0(num_signif,
                                   " out of ",
                                   num_total),
                            "No peptides differentially expressed",
                            icon = icon("thumbs-down", lib = "glyphicon"),
                            color = "red",
                            width = 4)
      }
      if(frac > 0 & frac <= 0.2) {
        info_box <-     infoBox("Significant peptides",
                              paste0(num_signif,
                                     " out of ",
                                     num_total),
                              paste0(signif(frac * 100, digits = 3),
                                     "% of peptides differentially expressed"),
                              icon = icon("thumbs-up", lib = "glyphicon"),
                              color = "green",
                              width = 4)
      }
      info_box
    })

    output$select_ptm <- renderUI({
      row_data <- rowData(dep_ptm())
      cols <- grep("_significant", colnames(row_data))
      names <- colnames(row_data)[cols]
      names <- gsub("_significant", "", names)
      selectizeInput(session$ns("select_ptm"),
                     "Select direct comparisons",
                     choices=names,
                     multiple = TRUE)
    })


    output$exclude_ptm <- renderUI({
      row_data <- rowData(dep_ptm())
      cols <- grep("_significant", colnames(row_data))
      names <- colnames(row_data)[cols]
      names <- gsub("_significant","",names)
      selectizeInput(session$ns("exclude_ptm"),
                     "Exclude direct comparisons",
                     choices = names,
                     multiple = TRUE)
    })


    output$heatmap_cntrst_ptm <- renderUI({
      if (!is.null(selected_ptm())) {
        df <- rowData(selected_ptm())
        cols <- grep("_significant$",colnames(df))
        selectizeInput(session$ns("heatmap_cntrst_ptm"),
                       "Contrast",
                       choices = gsub("_significant", "", colnames(df)[cols]), multiple = TRUE)
      }
    })


    output$mysplit_ptm <- renderUI({
      selectizeInput(session$ns("mysplit_ptm"),
                     "my split",
                     choices = c(1 : input$k_ptm), multiple = TRUE)

    })

    output$Custom_columns_order_ptm <- renderUI({
      selectizeInput(session$ns("Custom_columns_order_ptm"),
                     "Custom columns order",
                     choices = colnames(assay(dep_ptm())), multiple = TRUE)

    })

    output$Custom_columns_order_for_missval_heatmap_ptm <- renderUI({
      selectizeInput(session$ns("Custom_columns_order_for_missval_heatmap_ptm"),
                     "Custom columns order",
                     choices = colnames(assay(dep_ptm())), multiple = TRUE)
    })


    # output$volcano_cntrst_ptm <- renderUI({
    #   if (!is.null(selected_ptm())) {
    #     df <- rowData(selected_ptm())
    #     cols <- grep("_significant$",colnames(df))
    #     selectizeInput(session$ns("volcano_cntrst_ptm"),
    #                    "Contrast",
    #                    choices = gsub("_significant", "", colnames(df)[cols]))
    #   }
    # })

    output$volcano_cntrst_ptm <- renderUI({
      if (!is.null(selected_ptm())) {
        df <- rowData(selected_ptm())
        cols <- grep("_significant$",colnames(df))
        selectizeInput(session$ns("volcano_cntrst_ptm"),
                       label = "Contrast",
                       choices = gsub("_significant", "", colnames(df)[cols]))
      }
    })

    # output$Volcano_cntrst_stastical_plot_ptm <- renderUI({
    #   if (!is.null(selected_ptm())) {
    #     df <- rowData(selected_ptm())
    #     cols <- grep("_significant$",colnames(df))
    #     selectizeInput(session$ns("Volcano_cntrst_stastical_plot_ptm"),
    #                    "Contrast",
    #                    choices = gsub("_significant", "", colnames(df)[cols]))
    #   }
    # })

    output$cntrst_norm_distribution_ptm <- renderUI({
      if (!is.null(selected_ptm())) {
        df <- rowData(selected_ptm())
        cols <- grep("_significant$",colnames(df))
        selectizeInput(session$ns("cntrst_norm_distribution_ptm"),
                       "Contrast",
                       choices = gsub("_significant", "", colnames(df)[cols]))
      }
    })


    output$selected_proteins_ptm <- renderUI({
        row_Data <- rowData(dep_ptm())
        selectizeInput(session$ns("selected_proteins_ptm"),
                       "selected peptides",
                       choices = row_Data$name, selected = NULL, multiple = TRUE)
    })

    output$chooseToshow_ptm <- renderUI({
      selectizeInput(session$ns("chooseToshow_ptm"),
                     "selected peptides",
                     choices = name_for_choose_ptm(), selected = NULL, multiple = TRUE)
    })

    output$plot_stat_contrasts <- renderUI({
      selectInput(session$ns("plot_stat_contrasts"),
                  "plot contrasts",
                  choices = get_contrast(selected_ptm()),
                  multiple = T,
                  selected = get_contrast(selected_ptm()))
    })

    output$umap_n_neighbors_ptm <- renderUI({
      cat("umap_n_neighbors")
      sliderInput(inputId = session$ns("umap_n_neighbors_ptm"),"number of nearest neighbors",min=2,max=min(15,ncol(dep_ptm())),value=4)}
    )

    output$pca_top_n_ptm  <- renderUI(
      sliderInput(inputId = session$ns("pca_top_n_ptm"),"number of top variable proteins to consider",min = 2, max = nrow(dep_ptm()), value=500)
    )


    output$Tsne_perplexity_ptm <- renderUI(
      sliderInput(inputId = session$ns("Tsne_perplexity_ptm"),"perplexity",min=1,max= (ncol(dep_ptm())-1)/3 - 2^(-31) ,value=2)
    )



    ### Reactive functions ### ------------------------------------------------
    name_for_choose_ptm <- reactive({
      heatmap_name(dep = selected_ptm(), manual = input$manual_heatmap_ptm, manual_name = input$heatmap_cntrst_ptm)
      })

    excluded_ptm <- reactive({
      DEP2:::exclude_contrasts(dep_ptm(), input$exclude_ptm)
    })

    selected_ptm <- reactive({
      DEP2:::select_contrasts(excluded_ptm(), input$select_ptm)
    })

    res_ptm <- reactive({
      bbbb <<- selected_ptm()
      PTM_get_results(selected_ptm())
    })

    table_ptm <- reactive({
      aaaaa <<- res_ptm()
      PTM_get_table(res_ptm(), input$pres_ptm)
    })

    #### plots ----
    selected_plot_input_ptm <- reactive ({
      if(!is.null(input$table_ptm_rows_selected)) {
        selected_id <- table_ptm()[input$table_ptm_rows_selected,1]
        DEP2::plot_single(selected_ptm(), selected_id, input$pres_ptm)
      }
    })

    heatmap_input_ptm <- reactive({

      withProgress(message = 'Plotting', value = 0.66, {
        if(input$cluster_columns_ptm || is.null(input$Custom_columns_order_ptm) || length(input$Custom_columns_order_ptm) < 2) {
          selected = selected_ptm()
        } else {
          cols <- input$Custom_columns_order_ptm
          selected <- selected_ptm()[,cols]
        }

        DEP2::plot_heatmap(object = selected,
                           type = input$pres_ptm,
                           manual_contrast = input$heatmap_cntrst_ptm,
                           kmeans = TRUE,
                           k = input$k_ptm,
                           color = input$colorbar_ptm,
                           col_limit = input$limit_ptm,
                           row_font_size = input$row_font_size_ptm,
                           col_font_size = input$col_font_size_ptm,
                           cluster_columns = input$cluster_columns_ptm,
                           split_order = input$mysplit_ptm,
                           label_few_peptide_rows = F,
                           chooseToshow = input$chooseToshow_ptm,
                           )
        # plot_heatmap(dep = selected_ptm(),
        #              type = input$pres_ptm,
        #              manual = input$manual_heatmap_ptm,
        #              manual_name = input$heatmap_cntrst_ptm,
        #              kmeans = TRUE,
        #              k = input$k_ptm,
        #              color = input$colorbar_ptm,
        #              col_limit = input$limit_ptm,
        #              #same_trend = input$same_trend_ptm,
        #              row_font_size = input$row_font_size_ptm,
        #              col_font_size = input$col_font_size_ptm,
        #              cluster_columns = input$cluster_columns_ptm,
        #              if_mysplit = input$if_mysplit_ptm,
        #              mysplit = input$mysplit_ptm,
        #              if_rowname_color = FALSE,
        #              if_chooseToshow = input$if_chooseToshow_ptm,
        #              chooseToshow = input$chooseToshow_ptm,
        #              column_order = if(input$cluster_columns_ptm) {NULL} else {input$Custom_columns_order_ptm})
      })
    })


    # volcano_input_ptm <- reactive({
    #   if(!is.null(input$volcano_cntrst_ptm)) {
    #     plot_volcano(dep = selected_ptm(),
    #                  contrast = input$volcano_cntrst_ptm,
    #                  label_size = input$fontsize_ptm,
    #                  add_names = input$check_names_ptm,
    #                  adjusted = input$p_adj_ptm,
    #                  same_width = input$same_width_ptm,
    #                  my_breaks = input$my_breaks_ptm,
    #                  mybreaks = as.numeric(input$mybreaks_ptm))
    #   }
    # })


    custom_volcano_input_ptm <- reactive({
      if(!is.null(input$volcano_cntrst_ptm)) {
        withProgress(message = 'Plotting', value = 0.66, {
          DEP2::plot_volcano(object = selected_ptm(),
                             contrast = input$volcano_cntrst_ptm,
                             add_threshold_line = input$threshold_method_ptm,
                             adjusted = ifelse(input$threshold_method_ptm == "intersect", input$P_adj_ptm, TRUE),
                             label_size = input$fontSize_ptm,
                             label_number = input$showNum_ptm,
                             label_trend = input$labelWay_ptm,
                             chooseTolabel = input$selected_proteins_ptm,
                             x_symmetry = input$Same_width_ptm,
                             fcCutoff = input$lfc_ptm,
                             pCutoff = input$p_ptm,
                             curvature = input$curvature_ptm,
                             x0_fold = input$x0_fold_ptm,
                             dot_size = input$dotsize_ptm,
                             highlight_PGs_with_few_peptides = F,
                             down_color = input$down_color_ptm,
                             stable_color = input$stable_color_ptm,
                             up_color = input$up_color_ptm,
          )

          # my_plot_volcano(dep = selected_ptm(),
          #                 contrast = input$Volcano_cntrst_ptm,
          #                 adjusted = ifelse(input$threshold_method_ptm == "intersect", input$P_adj_ptm, TRUE),
          #                 labelWay = switch (input$labelWay_ptm,
          #                         "all significant" = "all significant",
          #                         "up" = "up",
          #                         "down" = "down",
          #                         "selected peptides" = "selected proteins"),
          #                 showNum = input$showNum_ptm,
          #                 chooseTolabel = input$selected_proteins_ptm,
          #                 fontSize = input$fontSize_ptm,
          #                 dotsize = input$dotsize_ptm,
          #                 same_width = input$Same_width_ptm,
          #                 fcCutoff = input$lfc_ptm,
          #                 adjpCutoff = input$p_ptm,
          #                 if_peptide_color = FALSE,
          #                 label.rectangle = input$if_label_rectangle_ptm,
          #                 stroke = input$stroke_ptm,
          #                 down_color = input$down_color_ptm,
          #                 stable_color = input$stable_color_ptm,
          #                 up_color = input$up_color_ptm,
          #                 threshold_method = input$threshold_method_ptm,
          #                 curvature = input$curvature_ptm ,
          #                 x0_fold = input$x0_fold_ptm,
          #                 PTM = TRUE
          # )
        })
      }
    })


    # stastical_plot_input_ptm <- reactive({
    #   if(!is.null(input$Volcano_cntrst_stastical_plot_ptm)){
    #     statistics_plot(dep = selected_ptm(),
    #                     x = input$stastical_plot_x_ptm, y = input$stastical_plot_y_ptm,
    #                     plottype = input$stastical_plot_type_ptm,
    #                     contrast=input$Volcano_cntrst_stastical_plot_ptm)
    #   }
    # })
    stastical_plot_input <- reactive({
      if(!is.null(input$plot_stat_contrasts)){
        # statistics_plot(dep = selected(),
        #                 x = input$stastical_plot_x, y = input$stastical_plot_y,
        #                 plottype = input$stastical_plot_type,
        #                 contrast=input$Volcano_cntrst_stastical_plot)
        DEP2::plot_statistics(selected_ptm(), statistic = input$stastical_plot_x ,contrasts = input$plot_stat_contrasts)
      }
    })

    norm_distribution_input_ptm <- reactive({
      if(!is.null(input$cntrst_norm_distribution_ptm)){
        DEP2::plot_norm_distribution(dep = dep_ptm(),
                                     contrast=input$cntrst_norm_distribution_ptm)
      }
    })

    pca_input_ptm <- reactive({
      DEP2::plot_pca(object = dep_ptm(), indicate = input$Indicate_ptm, if_square = input$if_square_ptm, n = input$pca_top_n_ptm, features = "peptides")
    })

    umap_input_ptm <- reactive({
      DEP2::plot_umap(dep_ptm(), indicate = input$umap_Indicate_ptm, if_square = input$umap_if_square_ptm, n_neighbors = input$umap_n_neighbors_ptm, features = "peptides")
    })

    Tsne_input_ptm <- reactive({
      DEP2::plot_Tsne(dep_ptm(), indicate = input$Tsne_Indicate_ptm, if_square = input$Tsne_if_square_ptm, perplexity = input$Tsne_perplexity_ptm, theseed = input$Tsne_theseed_ptm, features = "peptides")
    })

    Pearson_correlation_input_ptm <- reactive({
      DEP2::plot_cor(dep_ptm(), pal = input$Pearson_pal_ptm, pal_rev = input$Pearson_pal_rev_ptm, lower = input$Pearson_lower_ptm, upper = input$Pearson_upper_ptm, add_values = input$add_values_for_DEP_person_ptm, value_size = input$value_size_for_DEP_person_ptm, digits = input$value_digits_for_DEP_person_ptm)
    })

    meanSdPlot_input_ptm <- reactive({
      DEP2::meanSdPlot(dep_ptm())
    })

    Gowers_distance_input_ptm <- reactive({
      DEP2::plot_dist(dep_ptm(), pal = input$Gower_pal_ptm, pal_rev = input$Gower_pal_rev_ptm, add_values = input$add_values_for_DEP_gower_ptm, value_size = input$value_size_for_DEP_gower_ptm, digits = input$value_digits_for_DEP_gower_ptm)
    })

    Sample_CVs_input_ptm <- reactive({
      DEP2::plot_cvs(dep_ptm())
    })


    norm_input_ptm <- reactive({
      DEP2::plot_normalization(filt_ptm(),
                               norm_ptm())
    })

    missval_input_ptm <- reactive({
      DEP2::plot_missval(norm_ptm(),
                         cluster_columns = input$cluster_columns_for_missval_heatmap_ptm,
                         column_order = if(input$cluster_columns_for_missval_heatmap_ptm) {NULL} else {input$Custom_columns_order_for_missval_heatmap_ptm}
      )
    })

    detect_input_ptm <- reactive({
      DEP2::plot_detect(norm_ptm())
    })

    imputation_input_ptm <- reactive({
      DEP2::plot_imputation(norm_ptm(),
                            df_ptm())
    })

    numbers_input_ptm <- reactive({
      DEP2::plot_numbers(se = norm_ptm(), plot = TRUE, features_type = "peptides")
    })

    coverage_input_ptm <- reactive({
      DEP2::plot_coverage(se = norm_ptm(), plot = TRUE, features_type = "peptides")
    })

    ### Output functions ### --------------------------------------------------
    output$table_ptm <- DT::renderDataTable({
      table_ptm()
    }, options = list(pageLength = 25, scrollX = T),
    selection = list(selected = c(1)))

    output$selected_plot_ptm <- renderPlot({
      selected_plot_input_ptm()
    })

    output$heatmap_ptm <- renderPlot({
      print(heatmap_input_ptm())
    })

    # output$volcano_ptm <- renderPlot({
    #   volcano_input_ptm()
    # })

    output$custom_volcano_ptm <- renderPlot({
      custom_volcano_input_ptm()
    })

    output$stastical_plot_ptm <- renderPlot({
      stastical_plot_input()
    })

    output$norm_distribution_plot_ptm <- renderPlot({
      norm_distribution_input_ptm()
      })

    output$pca_ptm <- renderPlot({
      pca_input_ptm()
    })

    output$umap_ptm <- renderPlot({
      umap_input_ptm()
    })

    output$Tsne_ptm <- renderPlot({
      Tsne_input_ptm()
    })

    output$Pearson_correlation_ptm <- renderPlot({
      Pearson_correlation_input_ptm()
    })

    output$meanSdPlot_ptm <- renderPlot({
      meanSdPlot_input_ptm()
    })

    output$Gowers_distance_ptm <- renderPlot({
      Gowers_distance_input_ptm()
    })

    output$Sample_CVs_ptm <- renderPlot({
      Sample_CVs_input_ptm()
    })

    output$norm_ptm <- renderPlot({
      norm_input_ptm()
    })

    output$missval_ptm <- renderPlot({
      missval_input_ptm()
    })

    output$detect_ptm <- renderPlot({
      detect_input_ptm()
    })

    output$imputation_ptm <- renderPlot({
      imputation_input_ptm()
    })

    output$numbers_ptm <- renderPlot({
      numbers_input_ptm()
    })

    output$coverage_ptm <- renderPlot({
      coverage_input_ptm()
    })

    observe({
      output$plot_ptm <- renderUI({
        plotOutput(session$ns("heatmap_ptm"), height = (100 * as.numeric(input$size_ptm)))
      })
    })

    ### Download objects and functions ### ------------------------------------
    datasetInput_ptm <- reactive({
      switch(input$dataset_ptm,
             "results" = PTM_get_results(dep_ptm()),
             "significant_peptides" = PTM_get_results(
              dep_ptm()) %>%
               dplyr::filter(significant) %>%
               dplyr::select(-significant),
             "displayed_subset" = res_ptm() %>%
               dplyr::filter(significant) %>%
               dplyr::select(-significant),
             "full_dataset" = get_df_wide(dep_ptm()))
    })

    output$downloadData_ptm <- downloadHandler(
      filename = function() { paste(input$dataset_ptm, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput_ptm(),
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") }
    )

    output$downloadPlot_ptm <- downloadHandler(
      filename = function() {
        paste0("Barplot_", table_ptm()[input$table_ptm_rows_selected,1], ".pdf")
      },
      content = function(file) {
        pdf(file, width = input$selected_plot_Width_ptm, height = input$selected_plot_Height_ptm)
        print(selected_plot_input_ptm())
        dev.off()
      }
    )

    output$downloadHeatmap_ptm <- downloadHandler(
      filename = 'Heatmap.pdf',
      content = function(file) {
        pdf(file, width = input$Width_ptm, height = input$Height_ptm)
        print(heatmap_input_ptm())
        dev.off()
      }
    )

    output$downloadVolcano_ptm <- downloadHandler(
      filename = function() {
        paste0("Volcano_", input$volcano_cntrst_ptm, ".pdf")
      },
      content = function(file) {
        pdf(file, width = input$Volcano_Width_ptm, height = input$Volcano_Height_ptm)
        print(volcano_input_ptm())
        dev.off()
      }
    )

    output$download_custom_volcano_ptm <- downloadHandler(
      filename = function() {
        paste0("custom_Volcano_", input$volcano_cntrst_ptm, ".pdf")
      },
      content = function(file) {
        pdf(file, width = input$custom_volcano_Width_ptm, height = input$custom_volcano_Height_ptm)
        print(custom_volcano_input_ptm())
        dev.off()
      }
    )

    output$downloadstastical_plot_ptm <- downloadHandler(
      filename = function() {
        paste0("statistical_plot_", paste0(input$plot_stat_contrasts,collapse = "_"), ".pdf")
      },
      content = function(file) {
        pdf(file, width = input$stastical_plot_Width_ptm, height = input$stastical_plot_Height_ptm)
        print(stastical_plot_input())
        dev.off()
      }
    )

    output$download_norm_distribution_plot_ptm <- downloadHandler(
      filename = 'norm_distribution.pdf',
      content = function(file) {
        pdf(file, width = input$norm_distribution_Width_ptm, height = input$norm_distribution_Height_ptm)
        print(norm_distribution_input_ptm())
        dev.off()
      }
    )

    output$downloadPca_ptm <- downloadHandler(
      filename = 'Pca.pdf',
      content = function(file) {
        pdf(file, width = input$pca_Width_ptm, height = input$pca_Height_ptm)
        print(pca_input_ptm())
        dev.off()
      }
    )

    output$downloadUMAP_ptm <- downloadHandler(
      filename = 'UMAP.pdf',
      content = function(file) {
        pdf(file, width = input$umap_Width_ptm, height = input$umap_Height_ptm)
        print(umap_input_ptm())
        dev.off()
      }
    )

    output$downloadTSNE_ptm <- downloadHandler(
      filename = 'TSNE.pdf',
      content = function(file) {
        pdf(file, width = input$Tsne_Width_ptm, height = input$Tsne_Height_ptm)
        print(Tsne_input_ptm())
        dev.off()
      }
    )

    output$downloadmeanSdPlot_ptm <- downloadHandler(
      filename = 'meanSdPlot.pdf',
      content = function(file) {
        pdf(file, width = input$meanSdPlot_Width_ptm, height = input$meanSdPlot_Height_ptm)
        print(meanSdPlot_input_ptm())
        dev.off()
      }
    )

    output$download_Pearson_correlation_ptm <- downloadHandler(
      filename = 'Pearson_correlation.pdf',
      content = function(file) {
        pdf(file, width = input$Pearson_Width_ptm, height = input$Pearson_Height_ptm)
        print(Pearson_correlation_input_ptm())
        dev.off()
      }
    )

    output$download_Gowers_distance_ptm <- downloadHandler(
      filename = 'Gowers_distance.pdf',
      content = function(file) {
        pdf(file, width = input$Gower_Width_ptm, height = input$Gower_Height_ptm)
        print(Gowers_distance_input_ptm())
        dev.off()
      }
    )

    output$download_Sample_CVs_ptm <- downloadHandler(
      filename = 'Sample_CVs.pdf',
      content = function(file) {
        pdf(file, width = input$Sample_CVs_Width_ptm, height = input$Sample_CVs_Height_ptm)
        print(Sample_CVs_input_ptm())
        dev.off()
      }
    )

    output$downloadNorm_ptm <- downloadHandler(
      filename = "normalization.pdf",
      content = function(file) {
        pdf(file, width = input$norm_Width_ptm, height = input$norm_Height_ptm)
        print(norm_input_ptm())
        dev.off()
      }
    )

    output$downloadMissval_ptm <- downloadHandler(
      filename = "missing_values_heatmap.pdf",
      content = function(file) {
        pdf(file, width = input$missval_heatmap_Width_ptm, height = input$missval_heatmap_Height_ptm)
        print(missval_input_ptm())
        dev.off()
      }
    )

    output$downloadDetect_ptm <- downloadHandler(
      filename = "missing_values_quant.pdf",
      content = function(file) {
        pdf(file, width = input$detect_Width_ptm, height = input$detect_Height_ptm)
        gridExtra::grid.arrange(detect_input_ptm())
        dev.off()
      }
    )

    output$downloadImputation_ptm <- downloadHandler(
      filename = "imputation.pdf",
      content = function(file) {
        pdf(file, width = input$imputation_Width_ptm, height = input$imputation_Height_ptm)
        print(imputation_input_ptm())
        dev.off()
      }
    )

    output$downloadNumbers_ptm <- downloadHandler(
      filename = "numbers.pdf",
      content = function(file) {
        pdf(file, width = input$numbers_Width_ptm, height = input$numbers_Height_ptm)
        print(numbers_input_ptm())
        dev.off()
      }
    )

    output$downloadCoverage_ptm <- downloadHandler(
      filename = "coverage.pdf",
      content = function(file) {
        pdf(file, width = input$coverage_Width_ptm, height = input$coverage_Height_ptm)
        print(coverage_input_ptm())
        dev.off()
      }

    )
  })


      ###"question for DEP-PTM"
        observeEvent(input$help_format_DEP_ptm, {
          showModal(modalDialog(
            title = "Format specifications for DEP-PTM",
            includeMarkdown(system.file("extdata", "DEP_PTM.md", package = "DEP2")),
            h4("Example:"),
            tags$img(
              src = base64enc::dataURI(file = system.file("extdata", "DEP_PTM.png", package = "DEP2"), mime = "image/png"),
              width = 750
            ),
            easyClose = TRUE,
            footer = NULL,
            size = "l"
          ))
        })
        ### help imputation for ptm
        observeEvent(input$help_imputation_ptm, {
          showModal(modalDialog(
            title = "The detailed information of imputation methods",
            includeMarkdown(system.file("extdata", "impute.md", package = "DEP2")),
            easyClose = TRUE,
            footer = NULL,
            size = "l"
          ))
        })

        return(df_ptm)

      # },
      # ignoreNULL = TRUE,
      # once = T
      # )


    })
  }


