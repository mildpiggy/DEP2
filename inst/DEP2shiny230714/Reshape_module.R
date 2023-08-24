
Reshape_UI <- function(id){
  tagList(
    sidebarLayout(
      Reshape_sidebar_mod(id = id),
      Reshape_body_mod(id = id)
    ))
}

#* Reshape sidebar ----
Reshape_sidebar_mod <-  function(id,label="Reshape_sidabar"){
  ns = NS(id)
  cat(paste0("Reshape_sidebar ns is",ns(""),"\n"))

  tagList(
    sidebarPanel(
      width = 4,
      h3(id),
      bsCollapse(
        multiple = T,
        open = c("Files","Reshape options"),
        bsCollapsePanel("Files",
                        style = "primary",
                        fluidRow(
                          column(width = 12,
                                 fileInput(ns('table_for_reshape'),width = "300px",
                                           'Long table',
                                           accept=c('text/csv',
                                                    'text/comma-separated-values,text/plain',
                                                    '.csv','.tsv',".gz",".GZ")),
                          )
                        ),
                        fluidRow(
                          column(width = 12,
                                 selectInput(ns('sep'),
                                           "The field separator/delimiter character of table",
                                           choices = c("\\t" = "\t", ",", ";","|"),
                                           selected = "\\t"
                                 )
                          )
                        )
        ),
        bsCollapsePanel("filter", ## filter options pannel
                        style = "primary",
                        tagList(
                          fluidRow(
                            column(width = 12,
                                   # shiny::div(HTML(
                                   #   paste0("<h5>",
                                   #          "If need to filter the input 'Long table' before, ",
                                   #          "please fill in the filter conditions in the table below.",
                                   #          "Else please skip these step.",
                                   #          "</h5>")
                                   # )),
                                   h5("If need to filter the input 'Long table' before filter,
                                   please fill in the filter conditions in the table below and click filter button.
                                   Else please skip these step."),
                                   h6("Each complete row corresponds to a filter rule")
                            )
                          ),
                          fluidRow(column( width = 12,rhandsontable::rHandsontableOutput(ns("filtOpt_table")) )),
                          fluidRow(column( width = 12,actionButton(ns("filt_butt"),"Perform filter") )),
                          # fluidRow(shiny::verbatimTextOutput(ns("test_opt"))),
                          # fluidRow(shiny::verbatimTextOutput(ns("test_opt2"))),
                          fluidRow(column( width = 12,p("Filter based on rule(s):") )),
                          fluidRow(column( width = 12,shiny::verbatimTextOutput(ns("test_opt3")) )),
                          fluidRow(column( width = 12,shiny::verbatimTextOutput(ns("filter_message")) )),
                        )
        # )
                        # )

        ),
        bsCollapsePanel("Reshape options",
                        style = "primary",
                        uiOutput(ns("Reshape_opts")),
                        shiny::actionButton(ns("reshape_butt"),"Reshape table"),
                        br(),br(),
                        uiOutput(ns("down_butt"))
        )
      ),
      shinyBS::bsPopover(ns("sep"), "Delimiter",
                         content = paste0("<p>",
                                          "The delimiter of input table, ",
                                          "character used to separate individual columns within a table, such as:<br>",
                                          "CSV(Comma-Separated): <b>,</b> <br>",
                                          "TSV(Tab-Separated) or TXT:  <b>\\\\t</b> <br>",
                                          "Semicolon-Delimited: <b>;</b> <br>",
                                          "You check check the table <b>delimiter</b> through text editors like Notepad.",
                                          "</p>"), trigger = 'hover')

    ),
    shinyBS::bsTooltip(ns("rm_prefix"), "Remove the prefix of variables in the Sample column, like the file path",
                       "top", options = list(container = "body")),
    shinyBS::bsTooltip(ns("rm_suffix"), "Remove the suffix of variables in the Sample column",
                       "top", options = list(container = "body"))
  )
}



#* Reshape body ----
Reshape_body_mod <- function(id, label = "Reshape_body") {
  ns = NS(id)
  tagList(
    mainPanel(
      tabName = id,
      width = 8,
      ## ** The style of busy tooltip ----
      tags$head(tags$style(type="text/css", "
                        #loadmessage {
                        top: 0px; left: 0px;
                        width: 100%; padding: 5px 0px 5px 0px;
                        text-align: center; font-weight: bold;
                        font-size: 100%; color: #000000;
                        background-color: #FFC1C1; z-index: 105;}")), ##  提示条的样式
      conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                       tags$div("calculating...please wait...",id="loadmessage")),
      tabsetPanel(
        tabPanel(
          "Input table",
          dataTableOutput(ns("input_table")),
          # verbatimTextOutput(ns("input_table_str")),
          # verbatimTextOutput(ns("input_table_classes"))
        ),
        tabPanel(
          "Filtered table",
          dataTableOutput(ns("filtered_table")),
          # verbatimTextOutput(ns("input_table_str")),
          # verbatimTextOutput(ns("input_table_classes"))
        ),
        tabPanel(
          "Unique table",
          dataTableOutput(ns("unique_table"))
        ),
        tabPanel(
          "Output table",
          dataTableOutput(ns("output_table"))
        )
      )

    )

  )
}

#* Reshape server ----
Reshape_Server <- function(id, Omics_res) {
  ID = id
  print("Reshape server loading")
  moduleServer(
    ID,
    function(input, output, session) {
      ns = session$ns

      shinyBS::addPopover(session, "distPlot", "Data", content = paste0("<p>Waiting time between ",
                                                               "eruptions and the duration of the eruption for the Old Faithful geyser ",
                                                               "in Yellowstone National Park, Wyoming, USA.</p><p>Azzalini, A. and ",
                                                               "Bowman, A. W. (1990). A look at some data on the Old Faithful geyser. ",
                                                               "Applied Statistics 39, 357-365.</p>"), trigger = 'click')

      long_table <- reactive({
        inFile <- input$table_for_reshape
        if (is.null(inFile))
          return(NULL)
        fread(inFile$datapath, sep = input$sep, stringsAsFactors = F, header = T) %>% as.data.frame()
      })

      filt_table <- reactive({
        filt_condition = isolate(filt_condition())

        long_table = long_table()
        filt_butt = input$filt_butt

        if(is.null(filt_butt)||filt_butt == 0) return(NULL)
        if(is.null(long_table)) return(NULL)

        if(is.null(filt_condition) || length(filt_condition) < 1){
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "No filter rules.",
            text = paste0("Input the filter rules first"),
            type = "warning"
          )
          return(NULL) # no filter rules

        }else{
          filt_table <- long_table
          for(i in filt_condition){
            filt_table = try({
              filt_table %>% filter(!!rlang::parse_expr(i))
            })
            if(any(class(filt_table) == "try-error")){ # if filter error
              sendSweetAlert(
                session = shiny::getDefaultReactiveDomain(),
                title = "Error in filter",
                text = paste0(
                  "Something go wrong. Please check your input filter rules:\n",
                  i),
                type = "error"
              )
              return(NULL)
            }

            if(nrow(filt_table) < 1){ # if all rows is filtered out
              sendSweetAlert(
                session = shiny::getDefaultReactiveDomain(),
                title = "No row is remained after filter",
                text = paste0(
                  "All rows are filtered out. Please check your input filter rules:",
                  i),
                type = "warning"
              )
              return(NULL)
            }

          }

          return(filt_table)
        }

      })

      classes <- reactive({
        if(!is.null(long_table())){
          thedata <- long_table()
          classes <- 1:ncol(thedata) %>% sapply(function(x){
            class(thedata[,x])
          })
          classes
        }else{
          NULL
        }
      })

      unique_table <- reactive({
        a = input$reshape_butt
        unique_table <- NULL
        isolate({

          if(!(is.null(long_table())||
               is.null(input$feature_col)||all(input$feature_col == "")||
               is.null(input$sample_col)||all(input$sample_col == "")||
               is.null(input$expression_col)||all(input$expression_col == "")
          )){
            if(!is.null(filt_table()) && # filter is finnished
               all.equal(colnames(filt_table()), colnames(long_table())) # filtered table is from long table
            ){
              thedata = filt_table()
            }else thedata <- long_table()
            # thesampcol <- input$sample_col
            unique_table <- thedata[,c(input$feature_col,input$sample_col,input$expression_col)] %>% unique()
          }else{
            unique_table <- NULL
          }
        })
        return(unique_table)
      })

      out_table <- reactive({
        a = input$reshape_butt
        out_tb <- NULL
        isolate({
          if(!(is.null(long_table())||
               is.null(unique_table())
          )){
            if(is.null(input$shrink_col)|| all(input$shrink_col == "")){
              shrink_ident_cols = NULL
            }else{shrink_ident_cols = input$shrink_col}

            if(is.null(input$extend_col)|| all(input$extend_col == "")){
              extend_ident_cols = NULL
            }else{extend_ident_cols = input$extend_col}

            out_tb <- reshape_long2wide(long_table(),feature_col = input$feature_col,
                                        expression_col = input$expression_col,sample_col = input$sample_col,
                                        remove_sample_prefix = input$rm_prefix,remove_sample_suffix = input$rm_suffix,
                                        shrink_ident_cols = shrink_ident_cols,
                                        extend_ident_cols = extend_ident_cols)
          }else{out_tb <- NULL}
        })

        return(out_tb)
      })

      # The filter condition input
      filt_opt_table <- reactive({

        ## accept the value from rhandsontable
        opt_table = hot_to_r(input$filtOpt_table)
        if(is.null(opt_table)) return(NULL)

        ## only retain the complete rows from input filter opt table
        keep_rows = which(!MatrixGenerics::rowAnys(opt_table == ""))
        if(length(keep_rows) < 1) return(NULL)
        opt_table = opt_table[keep_rows,]
        return(opt_table)
      })

      filt_condition <- reactive({
        ## validate rules
        filt_opt_table <- filt_opt_table()

        if(!is.null(filt_opt_table)){
          Boolean_operators = "&"
          filt_condition = vector()

          # when long_table is changed
          if(any(!filt_opt_table$filteronColumn %in% colnames(long_table()))){
            cat
            return(NULL)
          }

          for (i in 1:nrow(filt_opt_table)) {
            # if(filt_opt_table$filteronColumn[i] %in% colnames(long_table()))
            if(filt_opt_table$judgmentRule[i] == "containCharecter" ){
              filt_condition[i] = paste0("grepl('",
                                         filt_opt_table$criteriaValue[i],
                                         "',",
                                         filt_opt_table$filteronColumn[i],
                                         ")"
              )
            }else if(filt_opt_table$judgmentRule[i] == "excludeCharecter" ){
              filt_condition[i] = paste0("!grepl('",
                                         filt_opt_table$criteriaValue[i],
                                         "',",
                                         filt_opt_table$filteronColumn[i],
                                         ")"
              )
            }else if(filt_opt_table$judgmentRule[i] %in% c(">","<","=","!=") ){
              filt_condition[i] = paste0(filt_opt_table$filteronColumn[i],
                                         filt_opt_table$judgmentRule[i],
                                         filt_opt_table$criteriaValue[i])
            }
          }
          return(filt_condition)
          # theformula = paste(filt_condition, collapse = Boolean_operators)
        }else{
          return(NULL)
        }
      })



      # The output rhandsontable for filter condition input, up to 7 rules
      output$filtOpt_table <- rhandsontable::renderRHandsontable({
        temp_table = data.frame(a= 1:4,b = 4:1,c = letters[1:4])
        if(is.null(long_table())) return(NULL)
        # cols = colnames(temp_table)
        cols = colnames(long_table())
        opt_table = data.frame(filteronColumn = factor(replicate(7,""),levels = c("",cols)),
                               judgmentRule = factor(replicate(7,""),level=c("",">","<","=","!=","containCharecter","excludeCharecter")),
                               criteriaValue = replicate(7,""),
                               stringsAsFactors = F)
        rhandsontable::rhandsontable(opt_table,
                                     colHeaders = c("Filter on Which Column","Judgment Rule", "Criteria Value"),
                                     rowHeaders = paste("Rule",1:7))
      })

      output$test_opt <- renderText({
        input$filtOpt_table -> temp
        hot_to_r(temp)[,1] %>% as.character()
      })

      output$test_opt2 <- renderText({
        filt_opt_table()
      })

      output$test_opt3 <- renderText({
        filt_condition()
      })

      output$filter_message <- renderText({
        if(!is.null(long_table()) && !is.null(filt_table())){
          paste0("Before filter: ",nrow(long_table())," rows.\n",
                 "After filter: ",nrow(filt_table())," rows."
          )
        }else NULL
      })

      output$Reshape_opts <- renderUI({
        tagList(
          pickerInput(
            ns("feature_col"),
            "Feature column",
            choices = colnames(long_table())[which(classes() == "character")],
            options = list(
              style = "warning",
              title = "Select the column of feature identifer")
          ),
          pickerInput(
            ns("sample_col"),
            "Sample column",
            choices = colnames(long_table())[which(classes() == "character")],
            options = list(
              style = "warning",
              title = "Select the column of samples")
          ),
          shiny::checkboxInput(ns("rm_prefix"),"Remove prefix of samples",value = T),
          shiny::checkboxInput(ns("rm_suffix"),"Remove suffix of samples",value = T),
          pickerInput(
            ns("expression_col"),
            "Expression column",
            choices = colnames(long_table()),
            options = list(
              style = "warning",
              title = "Select the column of expression/abundance values")
          ),
          pickerInput(
            ns("shrink_col"),
            "The column to shrink",
            choices = colnames(long_table()),
            options = list(
              style = "warning",
              title = "Select the column that will be shrinked in wide table"),
            multiple = T
          ),
          pickerInput(
            ns("extend_col"),
            "The column to extend",
            choices = colnames(long_table()),
            options = list(
              style = "warning",
              title = "Select the column that will be extended in wide table"),
            multiple = T
          ),

        )
      })

      output$down_butt <- renderUI({
        if(!is.null(out_table())){
          downloadButton(ns("down_button"),"down output table")
        }else{NULL}
      })

      output$down_button <- downloadHandler(
        filename = function() {
          paste0("out_table_", id, ".txt")
        },
        content = function(file) {
          write.table(out_table(),
                      file,
                      col.names = TRUE,
                      row.names = FALSE,
                      sep ="\t")
        }
      )

      output$input_table <- renderDataTable({
        if(!is.null(long_table())){
          thedata <- long_table()
          render_tb <- trim_table_character(thedata, length_lim = 14)
          suppressWarnings(render_tb)
        }else{
          NULL
        }
      },options = list(autoWidth = F, scrollX = T)
      )

      output$filtered_table <- renderDataTable({
        if(!is.null(filt_table())){
          thedata <- filt_table()
          suppressWarnings(trim_table_character(thedata, length_lim = 14))
        }else{
          NULL
        }
      },options = list(autoWidth = F, scrollX = T)
      )

      output$unique_table <- renderDataTable({
        if(!is.null(unique_table())){
          thedata <- unique_table()
          suppressWarnings(trim_table_character(thedata, length_lim = 14))
        }else{
          NULL
        }
      },options = list(autoWidth = F, scrollX = T)
      )

      output$output_table <- renderDataTable({
        if(!(is.null(long_table())||is.null(out_table()))){
          trim_table_character(out_table(), length_lim = 14)
        }else{
          NULL
        }
      },options = list(autoWidth = F, scrollX = T)
      )

    }
  )
}


trim_table_character <- function(thedata,length_lim = 14){
  classes <- 1:ncol(thedata) %>% sapply(function(x){
    class(thedata[,x])
  })
  character_cols = which(classes == "character")


  if(length(character_cols) > 1){
    thedata[,character_cols] = thedata[,character_cols] %>% apply(.,2,function(x,len_limit){
      # x = thedata[,1]
      t_df <- data.frame(ori = unique(x))
      t_df$trans_code = stringr::str_conv(t_df$ori,"UTF-8")
      t_df$trans_char = ifelse(stringr::str_length(t_df$trans_code) > len_limit,
                               paste0(stringr::str_sub(t_df$trans_code,1,len_limit),"..."),
                               t_df$trans_code)
      x <- t_df$trans_char[match(x,t_df$ori)]

      return(x)
    },len_limit = length_lim)
  }else if(length(character_cols) == 1){
    x = thedata[,character_cols]

    t_df <- data.frame(ori = unique(x))
    t_df$trans_code = stringr::str_conv(t_df$ori,"UTF-8")
    t_df$trans_char = ifelse(stringr::str_length(t_df$trans_code) > len_limit,
                             paste0(stringr::str_sub(t_df$trans_code,1,len_limit),"..."),
                             t_df$trans_code)
    x <- t_df$trans_char[match(x,t_df$ori)]
    thedata[,character_cols] = x

    # x = stringr::str_conv(x,"UTF-8")
    # thedata[,character_cols] =  ifelse(stringr::str_length(x) > len_limit, paste0(stringr::str_sub(x,1,len_limit),"..."),x)
  }

  return(thedata)
}
