
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
                                           "The field separator character",
                                           choices = c("\\t" = "\t", ",", ";"),
                                           selected = "\\t"
                                 )
                          )
                        )
        ),
        bsCollapsePanel("Reshape options",
                        style = "primary",
                        uiOutput(ns("Reshape_opts")),
                        shiny::actionButton(ns("reshape_butt"),"Reshape table"),
                        br(),br(),
                        uiOutput(ns("down_butt"))
        )
      )
    )
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

      long_table <- reactive({
        inFile <- input$table_for_reshape
        if (is.null(inFile))
          return(NULL)
        read.csv(inFile$datapath, sep = "\t", stringsAsFactors = F, header = T) %>% as.data.frame()
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
               is.null(input$feature_col)||input$feature_col == ""||
               is.null(input$sample_col)||input$sample_col == ""||
               is.null(input$expression_col)||input$expression_col == ""
          )){
            # thesampcol <- input$sample_col
            thedata <- long_table()
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
            if(is.null(input$shrink_col)|| input$shrink_col == ""){
              shrink_ident_cols = NULL
            }else{shrink_ident_cols = input$shrink_col}

            if(is.null(input$extend_col)|| input$extend_col == ""){
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
      # x111 <<- x
      # x = x111
      x = thedata[,1]
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
