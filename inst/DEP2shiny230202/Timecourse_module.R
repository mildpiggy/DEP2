library(TCseq)

Timecourse_UI <- function(id){
  tagList(
    sidebarLayout(
      Timecourse_sidebar_mod(id = id),
      Timecourse_body_mod(id = id)
    ))
}

#* Timecourse sidebar ----
Timecourse_sidebar_mod <-  function(id,label="Timecourse_sidabar"){
  ns = NS(id)
  cat(paste0("Timecourse_sidebar ns is",ns(""),"\n"))

  tagList(
    sidebarPanel(
      width = 3,
      h3(id),

      uiOutput(ns("Timecourse_drag_drop")),
      br(),
      actionButton(ns("analyze_for_Timecourse"), "Analyze")
    )
  )
}

#* Timecourse body ----
Timecourse_body_mod <- function(id, label = "Timecourse_body") {
  ns = NS(id)
  tagList(
    mainPanel(
      tabName = id,
      width = 9,
      fluidRow(

        column(width = 3,
               uiOutput(ns("groupby"))),
        column(width = 6,
               uiOutput(ns("Custom_columns_order"))),
        # column(width = 2,
        #        selectizeInput(ns("algo"),
        #                       "Clustering method",
        #                       choices = c("cm", "km", "hc"),
        #                       selected = c("cm"), multiple = FALSE)),
        column(width = 3,
               numericInput(ns("k"),
                            "k",
                            min = 1, max = 30, value = 10)
        )
      ),
      fluidRow(
        column(width = 3,
               numericInput(ns("row_font_size"),
                            "row font size",
                            min = 0, max = 15, value = 6)),
        column(width = 3,
               numericInput(ns("col_font_size"),
                            "col font size",
                            min = 0, max = 16, value = 10)),
        column(width = 3,
               numericInput(ns("Width"),
                            "heatmap width",
                            min = 1, max = 30, value = 7)),
        column(width = 3,
               numericInput(ns("Height"),
                            "heatmap height",
                            min = 1, max = 30, value = 10)),
      ),
      fluidRow(
        # column(numericInput(ns("size"),
        #               "Heat size(4-30)",
        #               min = 4, max = 30, value = 10),
        #  width = 3),
        column(width = 3,
               selectizeInput(ns("colorbar"),
                              "colorbar",
                              choices = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                              selected = c("RdBu"), multiple = FALSE)
        ),
        column(numericInput(ns("limit"),
                            "Color limit",
                            min = 0, max = 11, value = 3),
               width = 3),
        column(numericInput(ns("seed"),
                            "random seed for clustering",
                            min = 0, max = 100000, value = 10086),
               width = 3)
      ),
      fluidRow(
        column(uiOutput(ns("Timecourse_heatmap")),
               downloadButton(ns('downloadHeatmap'), 'Save heatmap'),
               width = 8),
        column(DT::dataTableOutput(ns("Timecourse_table")),
               downloadButton(ns('downloadTable'), 'Save table'),
               width = 4)
      ),

      shinyBS::bsTooltip(ns("row_font_size"), "Set the size of row labels", "top", options = list(container = "body")),
      shinyBS::bsTooltip(ns("col_font_size"), "Set the size of column labels", "top", options = list(container = "body")),
      shinyBS::bsTooltip(ns("Width"), "Width of the figure to export", "top", options = list(container = "body")),
      shinyBS::bsTooltip(ns("Height"), "Height of the figure to export", "top", options = list(container = "body")),
      shinyBS::bsTooltip(ns("Custom_columns_order"), "Set order of condition", "top", options = list(container = "body"))

    )


  )
}

#* Timecourse server ----
Timecourse_server_module <- function(id, Omics_res) {
  ID = id
  print("Timecourse server loading")
  moduleServer(
    ID,
    function(input, output, session) {
      ns = session$ns

      dropzone_input_list <- reactive({
        # drop_genelist_modules <- drop_genelist_modules_save
        # temp <- drop_genelist_modules()
        dropzone_input_list <- drop_genelist_modules %>% reactiveValuesToList()
        if(length(dropzone_input_list) == 0 ) return(NULL)
        dropzone_input_list_save <- dropzone_input_list[!unlist(map(dropzone_input_list,is.null))]
        if(length(dropzone_input_list_save) == 0 ) return(NULL)
        return(dropzone_input_list_save)
      })

      selected_omics <- reactive({
        if(is.null(dropzone_input_list())) return(NULL)
        selected_omics <- names(dropzone_input_list()) %>% sapply(., function(x){strsplit(x,"-ds-")[[1]][1]}) %>% unique
      })

      gene_name <- reactive({
        shiny::validate(need(length(selected_omics()) == 1, "Only perform on one omics result"))
        genelist <- dropzone_input_list()
        if(length(genelist)>0){
          gene_name <- genelist %>% lapply(., function(x){
            x()$name
          }) %>% unlist() %>% unique()
        }else{
          gene_name<- NULL
        }
        gene_name = data.frame(name = gene_name)

        return(gene_name)
      })

      selected_omics_result <- reactive({

        Omics_res[[selected_omics()]]()
      })

      ht_mat <- reactive({
        shiny::validate(need(length(selected_omics()) == 1, "Only perform on one omics result"))

        omictype = strsplit(selected_omics(),"_")[[1]][1]
        ht_mat <- assay(selected_omics_result())
        if(omictype == "RNAseq"){
          geneinfo = selected_omics_result()@geneinfo
          ht_mat = ht_mat[which(geneinfo$SYMBOL %in% gene_name()$name),]
          rownames(ht_mat) = make.names(geneinfo$SYMBOL[which(geneinfo$SYMBOL %in% gene_name()$name)], unique = T)
        }else if(omictype == "Proteomepg" | omictype == "Proteomepep"){
          ht_mat = ht_mat[ which(rownames(ht_mat) %in% gene_name()$name),]
        }
          return(ht_mat)
      })

      exp_design <- reactive({
        shiny::validate(need(length(selected_omics()) == 1, "Only perform on one omics result"))
        exp_design <- selected_omics_result()@colData
      })

      Omics_res_list <- reactive({
        # Omics_res = isolate(reactiveValuesToList(Omics_res)) ## break the interactive between omics results and drag input
        # fresh <- input$fresh_drag
        Omics_res <- reactiveValuesToList(Omics_res)
        Omics_res <- Omics_res[names(Omics_res) != ID]
        # Omics_res = reactiveValuesToList(omics_Servers_saved)
        check <- Omics_res %>% lapply(., function(x){
          return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
        }) %>% unlist
        Omics_res <- Omics_res[check] ## remove
      })

      # Omics_res_list()

      output$Timecourse_drag_drop <- renderUI({
        tagList(
          uiOutput(ns("ui_dragzone")),
          uiOutput(ns("ui_dropzone"))
        )
      })

      #** dragzone UI from genelist_tool ----
      output$ui_dragzone <- renderUI({
        # if(input$import_from_text_paste) {
        tagList(
          h5("Dragzone"),
          dragZone(ns("dragzone"),
                   choices = names(Omics_res_list()),
                   # list(one = "One",
                   #              two = "Two"),
                   style = "text-align: center;
                      font-size: 12px;
                      background-color:   #F8F8FF;
                      border-radius: 6px;
                      min-width: 80px;
                      color:  #6495ED;"
          )
          # ,actionButton(ns("fresh_drag"),"Fresh"),
          # shinyBS::bsTooltip(ns("fresh_drag"), "Fresh Drag-Drop input selections", "down", options = list(container = "body"))
        )
      })

      #** drop UI from genelist_tool ----
      exctract_genelist_UI_Timecourse <- exctract_genelist_UI_suffixFun("Timecourse")

      drop_choose <- reactive({
        x = lapply(names(Omics_res_list()), exctract_genelist_UI_Timecourse)
        names(x) = names(Omics_res_list())
        return(x <- x)
      })
      output$ui_dropzone <- renderUI({
        print("render dropZoneInput")
        # thedrop_choose <- isolate(drop_choose())
        thedrop_choose <- drop_choose()
        drop_input <- dropZoneInput(ns("Timecourse_dropzone"),
                                    choices = thedrop_choose,
                                    # list(one = exctract_genelist_UI("one"),
                                    #              two = exctract_genelist_UI("two")),
                                    # presets = list(values = names(drop_genelist_modules)[check]),
                                    selectable = F,
                                    multivalued = F,
                                    server = exctract_genelist_UI_Timecourse,
                                    style = "height: 450px;
                                    overflow-y: auto;
                                    color: #6495ED;
                                    font-color: black;
                                    border-radius: 6px;"
        )
        dropZoneServer(session, "Timecourse_dropzone", exctract_genelist_UI_Timecourse)
        return(tagList(h5("Dropzone"), drop_input))
      })

      # print("render Timecourse_dropzone finished")

      drop_genelist_modules <- reactiveValues()
      observeEvent(input$Timecourse_dropzone, {
        drop_genelist_modules_list = reactiveValuesToList(drop_genelist_modules)
        check = drop_genelist_modules_list %>% lapply(., function(x){
          return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
        }) %>% unlist
        purrr::map(
          setdiff(input$Timecourse_dropzone, names(drop_genelist_modules)[check]),
          ~ { drop_genelist_modules[[.]] <- exctract_genelist_Server(id = ., Omics_Serv_res = Omics_res_list, suffix = "Timecourse")
            # exctract_genelist_Server2(id = ., Omics_Serv_res = Omics_res_list, suffix = "Timecourse")
            }
        )
        cat(paste0("outside ns is ",session$ns(""),"\n"))
        purrr::map(
          setdiff(names(drop_genelist_modules)[check], input$Timecourse_dropzone),
          ~ { drop_genelist_modules[[.]] <- NULL })
        # print(input$Timecourse_dropzone)
        print("exctract_genelist_Server2 load")
        drop_genelist_modules
      }, priority = 1, ignoreNULL = F)


      output$groupby <- renderUI({
        shiny::validate(need(length(selected_omics()) == 1, "Only perform on one omics result"))
        selectizeInput(session$ns("groupby"),
                       "Set the groups by",
                       choices = colnames(exp_design()), multiple = F,
                       selected = "condition")
      })

      output$Custom_columns_order <- renderUI({
        shiny::validate(need(length(input$groupby) == 1, "Only perform on one omics result"))
        groups = exp_design()[,input$groupby]
        selectizeInput(session$ns("Custom_columns_order"),
                       "Set the groups order",
                       choices = groups, multiple = TRUE,
                       selected = NULL, width = '100%')
      })


      Timecourse_heatmap_res <- reactiveVal()
      observeEvent(input$analyze_for_Timecourse,{
        print("heatmap")
        ht_mat_1 <- ht_mat()
        exp_design <- exp_design()

        # Timecourse_heatmap_res <- Timecourse_heatmap_res(
        #   withProgress(message = 'Plotting', value = 0.66, {
        #     Timecourse_heatmap(ht_mat = ht_mat_save,
        #                        exp_design = exp_design_save,
        #                        groupby = input$groupby,
        #                        columns_order = input$Custom_columns_order,
        #                        algo = input$algo,
        #                        k = input$k,
        #                        color = input$colorbar,
        #                        col_limit = input$limit,
        #                        row_font_size = input$row_font_size,
        #                        col_font_size = input$col_font_size,
        #                        heatmap_height = input$Height,
        #                        heatmap_width = input$Width,
        #                        seed = input$seed
        #     )
        #   })
        # )
        Timecourse_heatmap_input <- reactive({
          shiny::validate( need( length(unique(exp_design[,input$groupby])) >2 , "Cluster is only meaningful for experiement with more the two groups!") )
          withProgress(message = 'Plotting', value = 0.66, {
            req(nrow(ht_mat_1) > 0)

            res = get_tc_cluster(ht_mat = ht_mat(),
                                 exp_design = exp_design,
                                 groupby = input$groupby,
                                 group_order = input$Custom_columns_order,
                                 # algo = input$algo,
                                 k = input$k,
                                 color = input$colorbar,
                                 col_limit = input$limit,
                                 row_font_size = input$row_font_size,
                                 col_font_size = input$col_font_size,
                                 heatmap_height = input$Height,
                                 heatmap_width = input$Width,
                                 seed = input$seed
            )

            Timecourse_heatmap_res = Timecourse_heatmap_res(res)
            return(res)
          })
        })

        output$heatmap <- renderPlot({
          req(!is.null(Timecourse_heatmap_input()))
          Timecourse_heatmap_input()$ht
        })
        output$Timecourse_heatmap <- renderUI({
          plotOutput(session$ns("heatmap"), height = (150 * as.numeric(input$Height)))
        })

        output$downloadHeatmap <- downloadHandler(
          filename = 'Heatmap.pdf',
          content = function(file) {
            req(!is.null(Timecourse_heatmap_input()))
            pdf(file, width = input$Width *1.9, height = input$Height*1.9)
            print(Timecourse_heatmap_input()$ht)
            dev.off()
          }
        )

        output$Timecourse_table <- DT::renderDataTable({
          req(!is.null(Timecourse_heatmap_input()))
          Timecourse_heatmap_input()$res
          # datatable(Timecourse_heatmap_res()$res)
        }, options = list(pageLength = 25, scrollX = T),style = "default")

        output$downloadTable <- downloadHandler(
          filename = 'Cluster_table.txt',
          content = function(file) {
            write.table(Timecourse_heatmap_input()$res, file,
                        col.names = TRUE, row.names = FALSE, sep ="\t")
          }
        )

      })

      return(Timecourse_heatmap_res)
    }
  )
}


# Module: Server exctract genelist item from choosed omics
# ID transmit from the dragInput choose. ns is (dragInput id)-(dragInput choose).
exctract_genelist_Server2 <- function(id, ID, Omics_Serv_res, suffix = ""){
  ID = id
  id <- paste0(id, "-", suffix)
  moduleServer(id,
               function(input,output,session){
                 ns = session$ns
                 # omics_info = get(".OmicsModule",envir = get_DEPglobal_env())
                 cat(paste0("the exctract_genelist_Server2 ns is ", ns(''), "\n"))
                 cat(paste0("the exctract_genelist_Server2 ID is ", ID, "\n"))

                 # OmicsModules = get(".OmicsModule",envir = get_DEPglobal_env())
                 # sta
                 omictype <- strsplit(ID,"_")[[1]][1]
                 if( omictype == "Proteomepg" | omictype == "Proteomepep"){
                   ## genelistUI in proteome
                   output$genelistUI <- renderUI({
                     # h3(id)
                     thechoices <- grep("_p.adj$",colnames(rowData(test_res())),value = T) %>% gsub("_p.adj", "", .)
                     tagList(
                       selectizeInput(ns("contrast"),
                                      "Contrast",
                                      choices = thechoices,
                                      multiple = T,
                                      selected = thechoices[1]),
                       fluidRow(
                         uiOutput(ns("trend")),
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
                     thechoices = grep("_padj$",colnames(test_res()@test_result),value = T) %>% gsub("_padj", "", .)
                     tagList(
                       selectizeInput(ns("contrast"),
                                      "Contrast",
                                      choices = thechoices,
                                      multiple = T,
                                      selected = thechoices[1]),
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
                 }

                 output$trend <- renderUI({
                   if(length(input$contrast) < 2){
                     column(width = 6, radioButtons(ns("trend"),
                                                    "Threshold method",
                                                    choices = c("both", "upregulate", "downregulate"),
                                                    selected = "both"))
                   }else return(NULL)

                 })

                 output$pcutoffval <- renderPrint({
                   pcutoff()
                 })

                 output$exctractValues <- renderPrint({
                   reactiveValuesToList(input)
                 })

                 output$table_bttn <- renderUI({
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
                   theomics <- Omics_Serv_res()[[omics_id]]()
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

                 if(omictype == "Proteomepg" | omictype == "Proteomepep"){
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

                     # if(nrow(dep_res) > 0){
                     sigcol <- paste0(input$contrast,"_significant")
                     req( all(sigcol %in% colnames( rowData(dep_res))) )
                     sigdf <- rowData(dep_res)[,sigcol] %>% as.data.frame()

                     sig <- dep_res[rowAnys( sigdf == TRUE, na.rm = T ), ]
                     if(nrow(sig) > 0){
                       if(length(sigcol) == 1){
                         sig_res <- get_results(sig) %>% .[,c(1:2,grep(paste0("(^",input$contrast,"_)", collapse = "|"),colnames(.)))]
                         if(input$trend == "both"){
                           return(sig_res_save2 <- sig_res)
                         }else if(input$trend == "upregulate"){
                           return(sig_res_save2 <- sig_res[(sig_res %>% dplyr::select(paste0(input$contrast,"_ratio"))) > 0,])
                         }else if(input$trend == "downregulate"){
                           return(sig_res_save2 <- sig_res[(sig_res %>% dplyr::select(paste0(input$contrast,"_ratio"))) < 0,])
                         }
                       }else if(length(sigcol) > 1){
                         return(sig_res_save2 <- get_results(sig) %>% .[,c(1:2,grep(paste0("(^",input$contrast,"_)", collapse = "|"),colnames(.)))])
                       }

                     }else{
                       return(NULL)
                     }
                   })
                 }else if(omictype == "RNAseq"){
                   sig_res = reactive({
                     the_res = test_res()@test_result %>% as.data.frame()
                     gene_info = test_res()@geneinfo
                     if(nrow(gene_info) == 0)
                       return(NULL)
                     the_res$name = gene_info$SYMBOL

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
                       sig_res_save <- sig_res
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

                 return(sig_res)
               })
}
