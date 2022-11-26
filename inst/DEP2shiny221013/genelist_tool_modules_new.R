# library(shiny)
library(dndselectr)
library(ggVennDiagram)
library(plotly)

# utility function
`%||%` <- function(a, b) if (is.null(a)) b else a

# Module: UI for exctract genelist item from choosed omics. drag and down
exctract_genelist_UI_suffixFun <- function(suffix=""){
  exctract_genelist_UI <- function(id, server=FALSE, session = getDefaultReactiveDomain()){
    ID = id
    id <- paste0(id, "-", suffix)
    ns <- NS(session$ns(id))
    ## The session here is the session of genelist_tool_Server. Session$ns is the ns in genelist_tool_UI testaaa here
    if(!server){
      print(paste0("the unservered exctract_genelist_UI import id is ",id,". the nsid is ",session$ns(id)))
    }

    genelistUI <- paste("genelistUI",id)

    if (server) {
      print(paste0("the servered exctract_genelist_UI import id is ",id,". the nsid is ",session$ns(id)))
      genelistUI <- uiOutput(ns("genelistUI"))
    }

    tagList(
      HTML("<style>
         .ds-dropzone .ds-dropoption { cursor: default; }
         .ds-handle { cursor: move; }
         </style>"),
      div(class = "blue",
          icon("bars", class="ds-handle"),
          br(),
          strsplit(ID,"-ds-")[[1]][1],
          genelistUI
          # switch(as.character(server),
          #        "TRUE" = verbatimTextOutput(ns('pcutoffval')),
          #        "FALSE" = NULL),
          # dataTableOutput(ns("restest"))
      )
    )
  }
  return(exctract_genelist_UI)
}

# Module: Server exctract genelist item from choosed omics
# ID transmit from the dragInput choose. ns is (dragInput id)-(dragInput choose). Omics_serv_res to render
exctract_genelist_Server <- function(id, ID, Omics_Serv_res, suffix = ""){
  ID = id
  id <- paste0(id, "-", suffix)
  moduleServer(id,
               function(input,output,session){
                 ns = session$ns
                 # omics_info = get(".OmicsModule",envir = get_DEPglobal_env())
                 cat(paste0("the exctract_genelist_Server ns is ", ns(''), "\n"))
                 cat(paste0("the exctract_genelist_Server ID is ", ID, "\n"))

                 # OmicsModules = get(".OmicsModule",envir = get_DEPglobal_env())
                 # sta
                 omictype <- strsplit(ID,"_")[[1]][1]
                 if( omictype == "Proteomepg" | omictype == "Proteomepep"| omictype == "PTM"){
                   ## genelistUI in proteome
                   output$genelistUI <- renderUI({
                     tagList(
                       selectizeInput(ns("contrast"),
                                      "Contrast",
                                      choices = get_contrast(test_res()),
                                        # grep("_p.adj$",colnames(rowData(test_res())),value = T) %>% gsub("_p.adj", "", .),
                                      multiple = T),
                       fluidRow(
                         column(width = 6, radioButtons(ns("trend"),
                                                        "Threshold method",
                                                        choices = c("both", "upregulate", "downregulate"),
                                                        selected = "both")),
                         column(width = 6, radioButtons(ns("method"),
                                                        "Threshold method",
                                                        choices = c("intersect", "curve"),
                                                        selected = "intersect"))
                       ),
                       conditionalPanel(condition = paste0("input['",ns("method"),"']","== 'intersect'"),
                                        fluidRow(
                                          column(width = 6,
                                                 numericInput(ns('pcutoff'),
                                                              "pcutoff",
                                                              min = 0, max = 0.2,
                                                              value =  0.05,step = 0.01)),
                                          column(width = 6, numericInput(ns('l2fc'),
                                                                         "l2fc",
                                                                         min = 0, max = 8,
                                                                         value = session$input[[ns('slider')]] %||% 1,step = 0.1))
                                        )
                       ) ,
                       conditionalPanel(condition = paste0("input['",ns("method"),"']","== 'curve'"),
                                        fluidRow(
                                          column(width = 6, numericInput(ns('curvature'),
                                                                         "curvature",
                                                                         min = 0, max = 1,
                                                                         value = session$input[[ns('slider')]] %||% 1,step = 0.01),),
                                          column(width = 6, numericInput(ns('x0_fold'),
                                                                         "x0 fold",
                                                                         min = 0, max = 8,
                                                                         value = session$input[[ns('slider')]] %||% 1,step = 0.1))
                                        )
                       ),
                       uiOutput(ns("table_bttn"))
                     )
                   })
                 }else if( omictype == "RNAseq"){
                   ## genelistUI in RNAseq
                   output$genelistUI <- renderUI({
                     tagList(
                       selectizeInput(ns("contrast"),
                                      "Contrast",
                                      choices = get_contrast(test_res()), multiple = T),
                       fluidRow(
                         column(width = 6, radioButtons(ns("trend"),
                                                        "Threshold method",
                                                        choices = c("both", "upregulate", "downregulate"),
                                                        selected = "both"))
                       ),
                       fluidRow(
                         column(width = 6,
                                numericInput(ns('pcutoff'),
                                             "pcutoff",
                                             min = 0, max = 0.2,
                                             value = 0.05,step = 0.01)),
                         column(width = 6, numericInput(ns('l2fc'),
                                                        "l2fc",
                                                        min = 0, max = 8,
                                                        value = session$input[[ns('slider')]] %||% 1,step = 0.1))
                       ),
                       uiOutput(ns("table_bttn"))
                     )
                   })
                 }else if(omictype == "Timecourse"){
                   output$genelistUI <- renderUI({
                     test_res_save <<- test_res()
                     tagList(
                       fluidRow(
                         column(width = 6,
                                selectizeInput(ns("cluster"),
                                               "Cluster",
                                               choices = test_res()$res$Timecourse_cluster, multiple = T))
                       ),
                       uiOutput(ns("table_bttn"))
                     )
                   })
                 }

                 output$pcutoffval <- renderPrint({
                   pcutoff()
                 })

                 output$exctractValues <- renderPrint({
                   reactiveValuesToList(input)
                 })

                 output$table_bttn <- renderUI({
                   sig_res_save3 <<- sig_res()
                   if(omictype == "RNAseq" && nrow(test_res()@geneinfo) == 0 ){ ## check id transform of RNAseq data
                     return(
                       tagList(
                         bsButton(ns("table_bttn"),label = "ID Transformation in RNAseq analysis module is necessary for post analysis",
                                  style = "warning", block = TRUE, type = "action", disabled = T)
                       )
                     )
                   }

                   if(!is.null(sig_res()) && !is.na(sig_res()) &&nrow(sig_res()) > 0){
                     return(
                       tagList(
                         bsButton(ns("table_bttn"),label = paste0("candidate list: ",nrow(sig_res())),
                                  style = "primary", block = TRUE, type = "action", disabled = F)
                       )
                     )
                   }else{
                     return(
                       tagList(
                         bsButton(ns("table_bttn"),label = paste0("no candidate match certain condition!"),
                                  style = "warning", block = TRUE, type = "action", disabled = T)
                       )
                     )
                   }
                 })


                 omics_id <- strsplit(ID, split = "-ds-")[[1]][1]
                 test_res = reactive({
                   # Omics_Serv_res_save <<- Omics_Serv_res
                   theomics <- Omics_Serv_res()[[omics_id]]()
                   test_res_save <<- theomics
                 })

                 iv1 <- InputValidator$new()
                 iv1$add_rule("pcutoff", sv_between(0.00001, 0.2))
                 iv1$add_rule("l2fc", sv_between(0, 10))
                 iv1$enable()

                 pcutoff <- reactive({
                   req(test_res)
                   req(iv1$is_valid())
                   input$pcutoff
                 })
                 l2fc <- reactive({
                   req(test_res)
                   req(iv1$is_valid())
                   input$l2fc
                 })


                 ## three type of omics data
                 if(omictype == "Proteomepg" | omictype == "Proteomepep" | omictype == "PTM"){
                   iv2 <- InputValidator$new()
                   iv2$add_rule("curvature", sv_between(0.1, 6))
                   iv2$add_rule("x0_fold", sv_between(1, 10))
                   iv2$enable()
                   curvature <- reactive({
                     req(iv2$is_valid())
                     input$curvature
                   })
                   x0_fold <- reactive({
                     req(iv2$is_valid())
                     input$x0_fold
                   })

                   sig_res <- reactive({
                     req(!is.null(input$method))  ## check input UI initialization was finished
                     dep_res <- add_rejections(diff = test_res(),
                                               alpha = pcutoff(), lfc = l2fc(),
                                               thresholdmethod = input$method,
                                               curvature = curvature(), x0_fold = x0_fold())
                     # dep_res_save <<- dep_res
                     # if(nrow(dep_res) > 0){
                     sigcol <- paste0(input$contrast,"_significant")
                     req( all(sigcol %in% colnames( rowData(dep_res))) )
                     sigdf <- rowData(dep_res)[,sigcol] %>% as.data.frame()
                     # sigcol_save <<- sigcol
                     sig <- dep_res[rowAnys( sigdf == TRUE, na.rm = T ), ]
                     sig_save <<- sig
                     sig_res <- get_results(sig)
                     sig_res_save <<- sig_res
                     if(omictype == "PTM"){
                       identifier_cols = 1:4 ## "name", "ID", "gene_name", "protein_ID"
                     }else{
                       identifier_cols = 1:2 ## "name", "ID"
                     }

                     if(nrow(sig) > 0){
                       if(length(sigcol) == 1){
                         sig_res <- sig_res %>% .[,c(identifier_cols,grep(paste0("(^",input$contrast,"_)", collapse = "|"),colnames(.)))]
                         if(input$trend == "both"){
                           return(sig_res_save2 <<- sig_res)
                         }else if(input$trend == "upregulate"){
                           return(sig_res_save2 <<- sig_res[(sig_res %>% dplyr::select(paste0(input$contrast,"_ratio"))) > 0,])
                         }else if(input$trend == "downregulate"){
                           return(sig_res_save2 <<- sig_res[(sig_res %>% dplyr::select(paste0(input$contrast,"_ratio"))) < 0,])
                         }
                       }else if(length(sigcol) > 1){
                         return(sig_res_save2 <<- sig_res %>% .[,c(identifier_cols,grep(paste0("(^",input$contrast,"_)", collapse = "|"),colnames(.)))])
                       }

                     }else{
                       return(NULL)
                     }
                   })

                 }else if(omictype == "RNAseq"){
                   sig_res <- reactive({
                     the_res = test_res()@test_result %>% as.data.frame()
                     gene_info = test_res()@geneinfo
                     if(nrow(gene_info) == 0)
                       return(NULL)
                     the_res$name = gene_info$SYMBOL

                     the_res_save <- the_res
                     thecontrast = input$contrast
                     pcol <- paste0(thecontrast,"_p.adj")
                     L2Fcol <- paste0(thecontrast,"_diff")
                     req( all(pcol %in% colnames( the_res)) && all(L2Fcol %in% colnames( the_res)) )

                     if(length(thecontrast) > 1){
                       sigdf = sapply(thecontrast, function(x){
                         pcol <- paste0(x,"_p.adj")
                         L2Fcol <- paste0(x,"_diff")
                         return(the_res[,pcol] < pcutoff() & abs(the_res[,L2Fcol] )> l2fc())
                       })
                       sig_res = the_res[rowAnys(sigdf, na.rm = T),] %>% dplyr::select(name,starts_with(thecontrast)) %>% rownames_to_column("ID")
                       if(nrow(sig_res) == 0) return(NULL)
                       return(sig_res)

                     }else if(length(thecontrast) == 1){
                       sig_rows = which(the_res[,pcol] < pcutoff() & abs(the_res[,L2Fcol] )> l2fc())
                       sig_res = the_res[sig_rows,] %>% dplyr::select(name,starts_with(thecontrast)) %>% rownames_to_column("ID")
                       sig_res_save <<- sig_res
                       if(nrow(sig_res) == 0) return(NULL)
                       if(input$trend == "both"){
                         return(sig_res)
                       }else if(input$trend == "upregulate"){
                         return(sig_res[(sig_res %>% dplyr::select(paste0(thecontrast,"_diff"))) > 0,])
                       }else if(input$trend == "downregulate"){
                         return(sig_res[(sig_res %>% dplyr::select(paste0(thecontrast,"_diff"))) < 0,])

                       }
                     }
                   })

                 }else if(omictype == "Timecourse"){
                   sig_res <- reactive({
                     req(!is.null(input$cluster))
                     sig_res = test_res()$res
                     sig_res = sig_res[which(sig_res$Timecourse_cluster %in% input$cluster),]
                     sig_res_save <<- sig_res
                   })
                 }

                 observeEvent(input$table_bttn, {
                   showModal(
                     modalDialog(footer = modalButton("close"),
                                 size = "m",
                                 easyClose = T,
                                 DT::DTOutput(ns("sig_table"))
                     )
                   )
                 })

                 output$sig_table <- DT::renderDT({
                   # sig_res_save <- sig_res()
                   sig_res()
                 }, options = list(pageLength = 10, scrollX = T))
                 # sig_res_saved <<- sig_res
                 return(sig_res)
               })
}


# Module: UI of whole genelist tool
genelist_tool_UI <- function(id,labelname = "GenelisttoolUI"){
  ns = NS(id)
  tagList(
    fluidRow(
      column(width = 2,
             uiOutput(ns("dragzone")),
             # actionBttn(
             #   inputId = ns("import_genelist"),
             #   label = "Import gene list",
             #   style = "bordered",
             #   color = "primary",
             #   size = "sm",
             #   icon = icon("retweet")
             # ),
             shinyBS::bsTooltip(ns("refresh"), "Click to refresh dragzone options and clear dropzone inputs", "top", options = list(container = "body")),
             # h5("Dropzone Values"),
             # verbatimTextOutput(ns("showme")),
             # textInput("item", label = h3("Item to append or select"), value = ""),
             # actionButton("append", label = "Append to dropzone"),
             # actionButton("select", label = "Select"),
             # actionButton("remove", label = "Remove selected"),
             # checkboxGroupInput("new_presets", label="Presets",
             #                    choices = c("one", "two"),
             #                    selected = c("one"),
             #                    inline = TRUE),
             # actionButton("update_presets", "Update Dropzone Presets")
             uiOutput(ns("inputlist_view")),
             # div(style="overflow-x: auto; height:400px;overflow-y: auto",
             #     textOutput(ns("Genelist1"))
             #     # verbatimTextOutput(ns("modules2"))
             # ),

      ),
      # column(width = 6,textOutput(ns("Genelist1"))),
      splitLayout(


        # uiOutput(ns("inputlist_view")),
        # div(
        #   uiOutput(ns("dragzone")),
        #   h5("Dropzone Values"),
        #   verbatimTextOutput(ns("showme")),
        #   # textInput("item", label = h3("Item to append or select"), value = ""),
        #   # actionButton("append", label = "Append to dropzone"),
        #   # actionButton("select", label = "Select"),
        #   # actionButton("remove", label = "Remove selected"),
        #   # checkboxGroupInput("new_presets", label="Presets",
        #   #                    choices = c("one", "two"),
        #   #                    selected = c("one"),
        #   #                    inline = TRUE),
        #   # actionButton("update_presets", "Update Dropzone Presets")
        #   div(style="overflow-x: auto; height:400px;overflow-y: auto",
        #
        #       verbatimTextOutput(ns("modules2"))
        #   )
        # ),
        div(
          column(width = 12, h4("Genelist 1"),
                 prettySwitch(
                   inputId = ns("if_paste_1"),
                   label = "Import by paste",
                   status = "primary",
                   fill = TRUE
                 )
          ),
          div(uiOutput(ns("ui_dropzone1")))
        ),
        div(
          column(width = 12, h4("Genelist 2"),
                 prettySwitch(
                   inputId = ns("if_paste_2"),
                   label = "Import by paste",
                   status = "primary",
                   fill = TRUE
                 )),
          uiOutput(ns("ui_dropzone2")),
        ),
        div(
          column(width = 12, h4("Genelist 3"),
                 prettySwitch(
                   inputId = ns("if_paste_3"),
                   label = "Import by paste",
                   status = "primary",
                   fill = TRUE
                 )),
          uiOutput(ns("ui_dropzone3")),
        ),
        div(
          column(width = 12, h4("Genelist 4"),
                 prettySwitch(
                   inputId = ns("if_paste_4"),
                   label = "Import by paste",
                   status = "primary",
                   fill = TRUE
                 )),
          uiOutput(ns("ui_dropzone4")),
        ),
        cellWidths = "20%"
      ),

    ),
    # fluidRow(verbatimTextOutput(ns("modules1"))),
    fluidRow(
      column(width = 4,
             h4("Venn plotly"),
             plotlyOutput(ns("venn_plotly"),height = "600px"),
             downloadButton(ns('downloadlist'), 'Save list')
      ),
      column(width = 8,
             # verbatimTextOutput(ns("venn_selected")), ## to show selected venn_click()
             # verbatimTextOutput(ns("venn_selected2")), ## to show selected venn_click2()
             h4("Subset Heatmap"),
             fluidRow(
               column(width = 1,
                      numericInput(ns("k"),
                                   "kmanes",
                                   min = 1, max = 20, value = 5)),
               column(width = 2,
                      selectizeInput(ns("heatmap_color"),
                                     "color panel",
                                     choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                                     selected = c("RdBu"), multiple = FALSE),
                      ),
               column(width = 2,
                      numericInput(ns("heatmap_color_limit"),
                                   "color limit",
                                   min = 0, max = 16, value = 4),
                      ),
               column(width = 2,
                      numericInput(ns("row_font_size"),
                                   "row font size",
                                   min = 0, max = 15, value = 6)),
               column(width = 2,
                      numericInput(ns("col_font_size"),
                                   "col font size",
                                   min = 0, max = 16, value = 10)),
               column(width = 2,
                      numericInput(ns("heatmap_width"),
                                   "heatmap width",
                                   min = 1, max = 30, value = 5)),
               column(width = 2,
                      numericInput(ns("heatmap_height"),
                                   "heatmap height",
                                   min = 1, max = 30, value = 8)),
             ),
             uiOutput(ns("Venn_heatmap")),
             downloadButton(ns('downloadHeatmap'), 'Save heatmap')
      ),
      # bsCollapse(id = ns("collapseExample"),
      #            multiple = T,
      #            bsCollapsePanel("venn plotly",
      #                            "f",
      #                            # uiOutput(ns("venn_plot")),
      #                            plotlyOutput(ns("venn_plotly")),
      #                            style = "success"),
      #            bsCollapsePanel("venn selected",
      #                            "f",
      #                            # uiOutput(ns("venn_plot")),
      # uiOutput(ns("overlayheatmap_plot"))
    ), ## render heatmap according dropzone inputs assay and event_data of venn plot plotly_click.
    # fluidRow(
    #   column(verbatimTextOutput(ns("inputvals1")),width = 3),
    #   actionButton(ns("save_modsession"),"save mod session"))
  )
}

# Module: Server of whole genelist tool. Omics_res is the test result of omics analysis
genelist_tool_Server <- function(id, Omics_res){

  # Omics_res <- reactiveValuesToList(omics_Severs_saved)
  # Omics_res <- reactiveValues(name1 = "awf", name2 = "awefwa")
  # cat()
  moduleServer(
    id,
    function(input,output,session){
      print(class(Omics_res))

      check_omics <- function(omics_res){
        Omics_res = reactiveValuesToList(Omics_res)
        check = Omics_res %>% lapply(., function(x){
          return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
        }) %>% unlist
        Omics_res = Omics_res[check]
      }

      # Omics_res_list <- reactiveVal(check_omics(isolate(omics_res)))
      Omics_res_list <- reactive( check_omics(omics_res) )

      observeEvent(input$refresh,{
        # Omics_res_list2 <- Omics_res_list(check_omics(isolate(omics_res)))
        updateDropZoneInput(session = session, inputId = "dropzone1", presets = list(value = character(0))
                            ,thechoices <- isolate(drop_choose1()) )
        updateDropZoneInput(session = session, inputId = "dropzone2", presets = list(value = character(0))
                            ,thechoices <- isolate(drop_choose2()) )
        updateDropZoneInput(session = session, inputId = "dropzone3", presets = list(value = character(0))
                            ,thechoices <- isolate(drop_choose3()) )
        updateDropZoneInput(session = session, inputId = "dropzone4", presets = list(value = character(0))
                            ,thechoices <- isolate(drop_choose4()) )
      })

      dropzone_input_list <- reactive({
        # drop_record <- drop_record()
        # print(paste0("drop_record: ",drop_record()))
        dropzone1_list <- layer_modules1 %>%reactiveValuesToList()
        dropzone2_list <- layer_modules2 %>% reactiveValuesToList()
        dropzone3_list <- layer_modules3 %>% reactiveValuesToList()
        dropzone4_list <- layer_modules4 %>% reactiveValuesToList()

        dropzone_input_list <- list(
          Dropzone1 = dropzone1_list,
          Dropzone2 = dropzone2_list,
          Dropzone3 = dropzone3_list,
          Dropzone4 = dropzone4_list
        )
        dropzone_input_list_save <<- dropzone_input_list

        dropzone_input_list2 = dropzone_input_list %>% lapply(., function(x){
          temp = c()
          if(length(x)>0){
            for(i in 1:length(x)){
              if(!is.null(x[[i]])){
                temp = c(temp,i)
              }
            }
            if(length(temp)>0)
              return(x[temp])
          }
          return(NULL)
        })
        dropzone_input_list2_save <<- dropzone_input_list2

        to_upper = input$to_upper
        if( length(dropzone_input_list2) > 0 ){
          dropzone_input_list3 = dropzone_input_list2 %>% lapply(., function(x){
            the_results = c()
            if(length(x) > 0 ){
              for (i in 1:length(x)) {
                the_res = x[[i]]()
                if(!is.na(the_res) && !is.null(the_res) && nrow(the_res) > 0 ){
                  the_results = c(the_results,the_res$name)
                }
              }
            }
            the_results = the_results %>% unique()
            if(length(the_results) > 0){
              if(to_upper)
                the_results = toupper(the_results) %>% make.unique()
              return(the_results)
            }else{
              return(NULL)
            }
          })
        }else{return(NULL)}
        dropzone_input_list3_save <<- dropzone_input_list3
        return(dropzone_input_list3)
      })

      venn_data <- reactive({
        dropzone_input_list3 = dropzone_input_list()

        if( length(dropzone_input_list3) > 1 ){
          venn <- ggVennDiagram::Venn(dropzone_input_list3)
          data <- ggVennDiagram::process_data(venn)
          return(data)
        }else{
          return(NULL)
        }

      })

      venn_output_df <- reactive({
        data <- venn_data()
        venn_list = data@region$item
        region_names = data@region$name
        df <- data.frame(sapply(venn_list, "[", i = 1:max(sapply(venn_list, length))))
        df[is.na(df)] <- ""
        colnames(df) <- region_names
        df_save <<- df
        return(df)
      })

      venn_items <- reactive({
        data <- venn_data()
        data_save <<- data
        req(!is.null(data))
        items <- data@region %>% dplyr::rowwise() %>%
          dplyr::mutate(text = stringr::str_wrap(paste0(.data$item, collapse = " "), width = 40)) %>%
          # sf::st_as_sf() %>%
          dplyr::mutate(ratio = round(count/sum(data@region$count),3)) %>%
          dplyr::mutate(count2 = paste(count,"(",ratio*100,"%)",sep = ""))
        label_coord = sf::st_centroid(items$geometry) %>% sf::st_coordinates()
        items$textx = label_coord[,1]
        items$texty = label_coord[,2]
        return(items)
      })

      venn_click <- reactive({
        # req(!is.null(venn_data()))
        click_list <- event_data(event = "plotly_selected",source = session$ns("venn_plotly"))
        return(click_list)
      })

      venn_click2 = reactiveVal()
      observeEvent(venn_click()$key,{
        # venn_click2 = venn_click2(NULL)
        # if(!is.null(venn_click()$key))  ## venn_click2 changed whenever the venn_click$key changed and wasn't null
          venn_click2 = venn_click2(venn_click()$key)

      }, ignoreNULL = F)


      # render dragzone according exited omics results
      output$dragzone <- renderUI({
        input$refresh
        tagList(
          h4("Dragzone"),
          actionBttn(
            inputId = session$ns("refresh"),
            label = "refresh",
            style = "bordered",
            color = "primary",
            size = "sm",
            icon = icon("retweet")
          ),
          br(),br(),
          dragZone(session$ns("dragzone"),
                   choices = names(isolate(Omics_res_list())),
                   # list(one = "One",
                   #              two = "Two"),
                   style = "text-align: center;
                      font-size: 12px;
                      background-color: 	#F8F8FF;
                      border-radius: 6px;
                      min-width: 80px;
                      color: 	#6495ED;"
          ),
          br(),
          prettySwitch(
            inputId = session$ns("to_upper"),
            label = "Transform to upper case",
            status = "primary",
            fill = TRUE
          )
        )
      })
      print("render dragzone finished")

      ## dropzone1
      exctract_genelist_UI1 <- exctract_genelist_UI_suffixFun("dropzoneA")
      drop_choose1 <- reactive({
        Omics_res_list = isolate(Omics_res_list())
        x = lapply(names(Omics_res_list), exctract_genelist_UI1)
        names(x) = names(Omics_res_list)
        return(x <- x)
      })
      output$ui_dropzone1 <- renderUI({
        temp <- tagList(conditionalPanel(condition = paste0("!input['",session$ns("if_paste_1"),"']"),
                                         dropZoneInput(session$ns("dropzone1"),
                                                       choices = isolate(drop_choose1()),
                                                       selectable = F,
                                                       multivalued = TRUE,
                                                       server = exctract_genelist_UI1,
                                                       style = "height: 600px;
                                    overflow-y: auto;
                                    color: #6495ED;
                                    font-color: black;
                                    border-radius: 6px;"
                                         )),
                        conditionalPanel(condition = paste0("input['",session$ns("if_paste_1"),"']"),
                                         column(width = 12,
                                                textAreaInput(inputId = session$ns("text_input_1"),
                                                              label = "Paste your gene list",
                                                              placeholder = "TP53\nPTEN\n\nor\n\nTP53  6.21\nPTEN  -1.53", rows = 27, resize = "none")
                                         ))
        )
        dropZoneServer(session, "dropzone1", exctract_genelist_UI1)
        return(temp)
      })
      print("render dropzone1 finished")

      # (2) This stores returned reactives from layer modules
      layer_modules1 <- reactiveValues()

      # (3) Call server module for new dropzone items
      observeEvent(
        list(a = input$dropzone1,
             b = input$if_paste_1),
        {
          if(input$if_paste_1){
            purrr::map(
              setdiff(names(layer_modules1), "paste"),
              ~ { layer_modules1[[.]] <- NULL })
            layer_modules1[["paste"]] = reactive({
              if(input$if_paste_1){
                temp <- input$text_input_1
                if(is.null(temp) || temp == "") return(NULL)
                genelist <- read.table(text = temp)
                if(nrow(genelist) == 0) return(NULL)
                # gene_name <- genelist[,1] %>% sapply( ., function(x){strsplit(x,";")[[1]][1] })
                gene_name <- sapply(genelist[,1],the1stname)
                gene_name = data.frame(name = gene_name %>% sapply(.,function(x){strsplit(x,"\\.")[[1]][1]})  %>% unlist)
                return(gene_name)
              }
              return(NULL)
            })

          }else{
            layer_modules1[["paste"]] = NULL
            layer_modules1_list = reactiveValuesToList(layer_modules1)
            check = layer_modules1_list %>% lapply(., function(x){
              return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
            }) %>% unlist
            purrr::map(
              setdiff(input$dropzone1, names(layer_modules1)[check]),
              ~ { layer_modules1[[.]] <- temp2 <- exctract_genelist_Server(id = ., Omics_Serv_res = Omics_res_list, suffix = "dropzoneA") }
            )
            cat(paste0("outside ns is ",session$ns(""),"\n"))
            purrr::map(
              setdiff(names(layer_modules1)[check], input$dropzone1),
              ~ { layer_modules1[[.]] <- NULL })
            print("exctract_genelist_Server load")
            # layer_modules1_saved <<- layer_modules1
          }
        }, priority = 1, ignoreNULL = F, ignoreInit = T)

      ## dropzone2
      exctract_genelist_UI2 <- exctract_genelist_UI_suffixFun("dropzoneB")
      drop_choose2 <- reactive({
        x = lapply(names(Omics_res_list()), exctract_genelist_UI2)
        names(x) = names(Omics_res_list())
        return(x <- x)
      })
      output$ui_dropzone2 <- renderUI({
        temp <- tagList(conditionalPanel(condition = paste0("!input['",session$ns("if_paste_2"),"']"),
                                         dropZoneInput(session$ns("dropzone2"),
                                                       choices = isolate(drop_choose2()),
                                                       selectable = F,
                                                       multivalued = TRUE,
                                                       server = exctract_genelist_UI2,
                                                       style = "height: 600px;
                                    overflow-y: auto;
                                    color: #6495ED;
                                    font-color: black;
                                    border-radius: 6px;"
                                         )),
                        conditionalPanel(condition = paste0("input['",session$ns("if_paste_2"),"']"),
                                         column(width = 12,
                                                textAreaInput(inputId = session$ns("text_input_2"),
                                                              label = "Paste your gene list",
                                                              placeholder = "TP53\nPTEN\n\nor\n\nTP53  6.21\nPTEN  -1.53", rows = 27, resize = "none")
                                         ))
        )
        dropZoneServer(session, "dropzone2", exctract_genelist_UI2)
        return(temp)
      })
      print("render dropzone2 finished")

      layer_modules2 <- reactiveValues()
      # gene_name <- reactive({
      #   if(input$if_paste_2){
      #     temp <<- input$text_input_2
      #     if(is.null(temp)) return(NULL)
      #     genelist <- read.table(text = temp)
      #     if(is.null(genelist) || length(genelist) == 0) return(NULL)
      #     gene_name <- genelist[,1] %>% sapply( ., function(x){strsplit(x,";")[[1]][1] })
      #
      #     gene_name = data.frame(name = gene_name %>% sapply(.,function(x){strsplit(x,"\\.")[[1]][1]})  %>% unlist)
      #     gene_name_save <<- gene_name
      #     return(gene_name)
      #   }
      #   return(NULL)
      # })
      observeEvent(
        list(a = input$dropzone2,
             b = input$if_paste_2),
        {
          if(input$if_paste_2){
            layer_modules2[["paste"]] = reactive({
              if(input$if_paste_2){
                temp <- input$text_input_2
                if(is.null(temp) || temp == "") return(NULL)
                genelist <- read.table(text = temp)
                if(nrow(genelist) == 0) return(NULL)
                gene_name <- sapply(genelist[,1],the1stname)
                gene_name = data.frame(name = gene_name %>% sapply(.,function(x){strsplit(x,"\\.")[[1]][1]})  %>% unlist)
                gene_name_save <<- gene_name
                return(gene_name)
              }
              return(NULL)
            })

          }else{
            layer_modules2[["paste"]] = NULL
            layer_modules2_list = reactiveValuesToList(layer_modules2)
            check = layer_modules2_list %>% lapply(., function(x){
              return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
            }) %>% unlist
            dropzone2_input = input$dropzone2
            # cat(dropzone2_input)
            purrr::map(
              setdiff(dropzone2_input, names(layer_modules2)[check]),
              ~ { layer_modules2[[.]] <- exctract_genelist_Server(id = ., Omics_Serv_res = Omics_res_list, suffix = "dropzoneB")  }
            )
            cat(paste0("outside ns is ",session$ns(""),"\n"))
            purrr::map(
              setdiff(names(layer_modules2)[check], dropzone2_input),
              ~ { layer_modules2[[.]] <- NULL })
            # layer_modules2_saved <<- layer_modules2
          }
        },
        priority = 1, ignoreNULL = F, ignoreInit = T)

      ## dropzone3
      exctract_genelist_UI3 <- exctract_genelist_UI_suffixFun("dropzoneC")
      drop_choose3 <- reactive({
        x = lapply(names(Omics_res_list()), exctract_genelist_UI3)
        names(x) = names(Omics_res_list())
        return(x <- x)
      })
      output$ui_dropzone3 <- renderUI({
        temp <- tagList(conditionalPanel(condition = paste0("!input['",session$ns("if_paste_3"),"']"),
                                         dropZoneInput(session$ns("dropzone3"),
                                                       choices = isolate(drop_choose3()),
                                                       selectable = F,
                                                       multivalued = TRUE,
                                                       server = exctract_genelist_UI3,
                                                       style = "height: 600px;
                                    overflow-y: auto;
                                    color: #6495ED;
                                    font-color: black;
                                    border-radius: 6px;"
                                         )),
                        conditionalPanel(condition = paste0("input['",session$ns("if_paste_3"),"']"),
                                         column(width = 12,
                                                textAreaInput(inputId = session$ns("text_input_3"),
                                                              label = "Paste your gene list",
                                                              placeholder = "TP53\nPTEN\n\nor\n\nTP53  6.21\nPTEN  -1.53", rows = 27, resize = "none")
                                         ))
        )
        dropZoneServer(session, "dropzone3", exctract_genelist_UI3)
        return(temp)
      })
      print("render dropzone3 finished")

      layer_modules3 <- reactiveValues()
      observeEvent(
        list(a = input$dropzone3,
             b = input$if_paste_3),
        {
          if(input$if_paste_3){
            print("pate input")
            layer_modules3[["paste"]] = reactive({
              if(input$if_paste_3){
                temp <- input$text_input_3
                if(is.null(temp) || temp == "") return(NULL)
                genelist <- read.table(text = temp)
                if(nrow(genelist) == 0) return(NULL)
                # gene_name <- genelist[,1] %>% sapply( ., function(x){strsplit(x,";")[[1]][1] })
                gene_name <- sapply(genelist[,1],the1stname)
                gene_name = data.frame(name = gene_name %>% sapply(.,function(x){strsplit(x,"\\.")[[1]][1]})  %>% unlist)
                return(gene_name)
              }
              return(NULL)
            })

          }else{
            layer_modules3[["paste"]] = NULL
            layer_modules3_list = reactiveValuesToList(layer_modules3)
            check = layer_modules3_list %>% lapply(., function(x){
              return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
            }) %>% unlist
            purrr::map(
              setdiff(input$dropzone3, names(layer_modules3)[check]),
              ~ { layer_modules3[[.]] <- exctract_genelist_Server(id = ., Omics_Serv_res = Omics_res_list, suffix = "dropzoneC")  }
            )
            cat(paste0("outside ns is ",session$ns(""),"\n"))
            purrr::map(
              setdiff(names(layer_modules3)[check], input$dropzone3),
              ~ { layer_modules3[[.]] <- NULL })
          }
        }, priority = 1, ignoreNULL = F, ignoreInit = T)

      ## dropzone4
      exctract_genelist_UI4 <- exctract_genelist_UI_suffixFun("dropzoneD")
      drop_choose4 <- reactive({
        x = lapply(names(Omics_res_list()), exctract_genelist_UI4)
        names(x) = names(Omics_res_list())
        return(x <- x)
      })
      output$ui_dropzone4 <- renderUI({
        temp <- tagList(conditionalPanel(condition = paste0("!input['",session$ns("if_paste_4"),"']"),
                                         dropZoneInput(session$ns("dropzone4"),
                                                       choices = isolate(drop_choose4()),
                                                       selectable = F,
                                                       multivalued = TRUE,
                                                       server = exctract_genelist_UI4,
                                                       style = "height: 600px;
                                    overflow-y: auto;
                                    color: #6495ED;
                                    font-color: black;
                                    border-radius: 6px;"
                                         )),
                        conditionalPanel(condition = paste0("input['",session$ns("if_paste_4"),"']"),
                                         column(width = 12,
                                                textAreaInput(inputId = session$ns("text_input_4"),
                                                              label = "Paste your gene list",
                                                              placeholder = "TP53\nPTEN\n\nor\n\nTP53  6.21\nPTEN  -1.53", rows = 27, resize = "none")
                                         ))
        )
        dropZoneServer(session, "dropzone4", exctract_genelist_UI4)
        return(temp)
        # if(!input$if_paste_4){
        #   print("render dropZoneInput")
        #   drop_input <- dropZoneInput(session$ns("dropzone4"),
        #                               choices = isolate(drop_choose4()),
        #                               selectable = F,
        #                               multivalued = TRUE,
        #                               server = exctract_genelist_UI4,
        #                               style = "height: 600px;
        #                             overflow-y: auto;
        #                             color: #6495ED;
        #                             font-color: black;
        #                             border-radius: 6px;")
        #   dropZoneServer(session, "dropzone4", exctract_genelist_UI4)
        #   return(drop_input)
        # }else {
          # column(width = 12,
          #        textAreaInput(inputId = session$ns("text_input_4"),
          #                      label = "Paste your gene list",
          #                      placeholder = "TP53\nPTEN\n\nor\n\nTP53  6.21\nPTEN  -1.53", rows = 27, resize = "none")
          # )
        #
        # }
        # })
      })
      print("render dropzone4 finished")

      layer_modules4 <- reactiveValues()

      observeEvent(
        list(a = input$dropzone4,
             b = input$if_paste_4),
        {
          if(input$if_paste_4){
            print("pate input")
            layer_modules4[["paste"]] = reactive({
              if(input$if_paste_4){
                temp <- input$text_input_4
                if(is.null(temp) || temp == "") return(NULL)
                genelist <- read.table(text = temp)
                if(nrow(genelist) == 0) return(NULL)
                # gene_name <- genelist[,1] %>% sapply( ., function(x){strsplit(x,";")[[1]][1] })
                gene_name <- sapply(genelist[,1],the1stname)
                gene_name = data.frame(name = gene_name %>% sapply(.,function(x){strsplit(x,"\\.")[[1]][1]})  %>% unlist)
                return(gene_name)
              }
              return(NULL)
            })

          }else{
            layer_modules4[["paste"]] = NULL
            layer_modules4_list = reactiveValuesToList(layer_modules4)
            check = layer_modules4_list %>% lapply(., function(x){
              return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
            }) %>% unlist
            purrr::map(
              setdiff(input$dropzone4, names(layer_modules4)[check]),
              ~ { layer_modules4[[.]] <- exctract_genelist_Server(id = ., Omics_Serv_res = Omics_res_list, suffix = "dropzoneD")  }
            )
            purrr::map(
              setdiff(names(layer_modules4)[check], input$dropzone4),
              ~ { layer_modules4[[.]] <- NULL })
            # layer_modules4_saved <<- layer_modules4
          }
        }, priority = 1, ignoreNULL = F, ignoreInit = T)

      # observeEvent(input$save_modsession, {
      #   session_saved2 <<- session
      #   session_saved3 <<- session$parent
      #   session_saved4 <<- .subset2(session, "parent")
      #   re_domain <<- getDefaultReactiveDomain() ## the session is the getDefaultReactiveDomain()
      #   # cat(aaa())
      #   # aaa_saved <<- aaa
      # })


      output$showme <- renderText({
        paste("Dropzone:", paste0(input$dropzone, collapse = " "))
      })

      # output$inputvals1 <- renderPrint({
      #   reactiveValuesToList(input)
      # })

      # output$venn_selected <- renderPrint({
      #   venn_click()
      # })

      # output$venn_selected2 <- renderPrint({
      #   venn_click2()
      # })

      heatmap_plot = reactive({
        library(ComplexHeatmap)

        venn_click2 <- venn_click2() %>% unique()
        # venn_click2_save <<- venn_click2
        if(is.null(venn_click2))
          return(NULL)

        items = venn_items()
        data = venn_data()
        validate(need(!is.null(venn_data()), "At least need two list"))
        selected_proteins <- data@region$item[which(data@region$id %in% venn_click2)] %>% unlist() %>% unique()

        print(input$dropzone1)
        dropzone_selected_omics <- c(input$dropzone1, input$dropzone2, input$dropzone3, input$dropzone4)
        dropzone_selected_omics2 <- dropzone_selected_omics %>% sapply(., function(x){strsplit(x,"-ds-")[[1]][1]}) %>% unique

        # omics_list <- Omics_res %>% reactiveValuesToList() %>% .[dropzone_selected_omics2] %>% lapply(function(x) x())
        omics_list <- Omics_res_list() %>% .[dropzone_selected_omics2] %>% lapply(function(x) x())
        multi_ht <- plot_multi_heatmap(omics_list = omics_list, choose_name = selected_proteins, to_upper = input$to_upper,
                                       color = input$heatmap_color,col_limit = input$heatmap_color_limit,
                                       width = input$heatmap_width,
                                       height = input$heatmap_height,
                                       km = input$k)

        multi_ht
      })

      output$Venn_heatmap <- renderUI({
        plotOutput(session$ns("heatmap_selected"), height = (155 * as.numeric(input$heatmap_height)))
      })

      output$heatmap_selected <- renderPlot({
        req(!is.null(heatmap_plot()))
        heatmap_plot()
      })

      output$downloadHeatmap <- downloadHandler(
        filename = 'Heatmap.pdf',
        content = function(file) {
          req(!is.null(heatmap_plot()))
          pdf(file, width = input$heatmap_width *2.5, height = input$heatmap_height*2.5)
          print(heatmap_plot())
          dev.off()
        }
      )

      output$downloadlist <- downloadHandler(
        filename = 'Venn_lists.txt',
        content = function(file) {
          write.table(venn_output_df(),file,sep = "\t",
                      col.names = TRUE,
                      row.names = FALSE)
        }
      )

      output$modules1 <- renderPrint({
        input$dropzone2
      })

      output$inputlist_view <- renderUI({
        dropzone_input_list <- dropzone_input_list()
        # dropzone_input_list <- dropzone_input_list_save
        if (!is.null(dropzone_input_list[["Dropzone1"]]) && length(dropzone_input_list[["Dropzone1"]]) >0) {
          tag1 <- tagList(
            p("proteins(genes) in Dropzone1",style = "font-weight:bold"),
            div( style=paste0("overflow-x: auto; height:100px;overflow-y: auto; ",
                              "font-size:xx-small; ") ,
                 textOutput(session$ns("Genelist1"))
            ),
            br()
          )
        }else tag1 <- NULL

        if (!is.null(dropzone_input_list[["Dropzone2"]]) && length(dropzone_input_list[["Dropzone2"]]) >0) {
          tag2 <- tagList(
            p("proteins(genes) in Dropzone2",style = "font-weight:bold"),
            div( style=paste0("overflow-x: auto; height:100px;overflow-y: auto; ",
                              "font-size:xx-small; ") ,
                 textOutput(session$ns("Genelist2"))
            ),
            br()
          )
        }else tag2 <- NULL

        if (!is.null(dropzone_input_list[["Dropzone3"]]) && length(dropzone_input_list[["Dropzone3"]]) >0) {
          tag3 <- tagList(
            p("proteins(genes) in Dropzone3",style = "font-weight:bold"),
            div( style=paste0("overflow-x: auto; height:100px;overflow-y: auto; ",
                              "font-size:xx-small; ") ,
                 textOutput(session$ns("Genelist3"))
            ),
            br()
          )
        }else tag3 <- NULL

        if (!is.null(dropzone_input_list[["Dropzone4"]]) && length(dropzone_input_list[["Dropzone4"]]) >0) {
          tag4 <- tagList(
            p("proteins(genes) in Dropzone4",style = "font-weight:bold"),
            div( style=paste0("overflow-x: auto; height:100px;overflow-y: auto; ",
                              "font-size:xx-small; ") ,
                 textOutput(session$ns("Genelist4"))
            ),
            br()
          )
        }else tag4 <- NULL

        tagList(
          # br(),
          h5("Drop input viewer",style = "font-weight:bold"),
          tag1, tag2,tag3,tag4
        )

      })

      output$Genelist1 <- renderPrint({
        # dropzone_input_list <- dropzone_input_list_save
        dropzone_input_list <- dropzone_input_list()
        if(!is.null(dropzone_input_list[["Dropzone1"]]) && length(dropzone_input_list[["Dropzone1"]] >0)){
          cat(paste0(dropzone_input_list[["Dropzone1"]],collapse = " "))
        }else{
          ""
        }
      })

      output$Genelist2 <- renderPrint({
        # dropzone_input_list <- dropzone_input_list_save
        dropzone_input_list <- dropzone_input_list()
        if(!is.null(dropzone_input_list[["Dropzone2"]]) && length(dropzone_input_list[["Dropzone2"]] >0)){
          cat(paste0(dropzone_input_list[["Dropzone2"]],collapse = " "))
        }else{
          NULL
        }
      })

      output$Genelist3 <- renderPrint({
        # dropzone_input_list <- dropzone_input_list_save
        dropzone_input_list <- dropzone_input_list()
        if(!is.null(dropzone_input_list[["Dropzone3"]]) && length(dropzone_input_list[["Dropzone3"]] >0)){
          cat(paste0(dropzone_input_list[["Dropzone3"]],collapse = " "))
        }else{
          ""
        }
      })

      output$Genelist4 <- renderPrint({
        # dropzone_input_list <- dropzone_input_list_save
        dropzone_input_list <- dropzone_input_list()
        if(!is.null(dropzone_input_list[["Dropzone4"]]) && length(dropzone_input_list[["Dropzone4"]] >0)){
          cat(paste0(dropzone_input_list[["Dropzone4"]],collapse = " "))
        }else{
          ""
        }
      })


      output$venn_plotly <- renderPlotly({
        data <- venn_data()
        items <- venn_items()
        shiny::validate(need(!is.null(venn_data()), "Must provide two input"))
        # items <<- data@region %>% dplyr::rowwise() %>%
        #   dplyr::mutate(text = stringr::str_wrap(paste0(.data$item, collapse = " "), width = 40)) %>%
        #   # sf::st_as_sf() %>%
        #   dplyr::mutate(ratio = round(count/sum(data@region$count),3)) %>%
        #   dplyr::mutate(count2 = paste(count,"(",ratio*100,"%)",sep = ""))
        # label_coord = sf::st_centroid(items$geometry) %>% sf::st_coordinates()
        # items$textx = label_coord[,1]
        # items$texty = label_coord[,2]

        p1 <- items %>% highlight_key(~id) %>% {
          ggplot(.) +
            geom_sf(aes_string(fill = "count"),lwd=0.5,color ="grey80") +
            geom_sf_text(aes_string(label = "name"), data = data@setLabel,inherit.aes = F) +
            geom_text(aes_string(label = "count2", text = "text", x = "textx", y = "texty"), show.legend = FALSE,size=4) +
            scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
            theme_void()} %>%
          ggplotly(tooltip = "text",
                   # layout.dragmode = 'select',
                   source  = session$ns("venn_plotly")) %>%
          highlight(on = c("plotly_selected"),off = "plotly_deselect",layout.dragmode = 'select')
        # %>%
        #   layout(dragmode = "select")

        return(p1)
      })
    }
  )
}
