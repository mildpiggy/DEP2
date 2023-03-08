source("GSEA_submodule.R")
source("source_enrich_ana.R")


#* GSEA UI ----
GSEA_UI <- function(id){
  tagList(
    sidebarLayout(
      GSEA_sidebar_mod(id = id),
      GSEA_body_mod(id = id)
    ))
}

#* GSEA sidebar ----
GSEA_sidebar_mod <-  function(id,label="GSEA_sidabar"){
  ns = NS(id)
  # cat(paste0("GSEA_sidebar ns is",ns(""),"\n"))

  tagList(
    sidebarPanel(
      width = 3,
      h3(id),
      checkboxInput(ns("import_from_text_paste"), "import from text paste", value = FALSE),

      conditionalPanel(condition = paste0("!input['",ns("import_from_text_paste"),"']"),#&& "input.filt_num_cutoff_ptm !== ''"
                       tagList(
                         uiOutput(ns("GSEA_drag_drop"))
                         # uiOutput(ns("GSEA_dragzone")),
                         # uiOutput(ns("ui_dropzone"))
                       )
      ),
      conditionalPanel(condition = paste0("input['",ns("import_from_text_paste"),"']"),#&& "input.filt_num_cutoff_ptm !== ''"
                       textAreaInput(inputId = ns("text_input_for_GSEA"), label = "Please paste your gene list", placeholder = "TP53  6.21\nPTEN  -1.53", rows = 8, width = "100%")
      ),

      # uiOutput(ns("text_input_for_GSEA")),
      numericInput(ns("topn"),"Top n genes", value = NULL),
      selectizeInput(ns("type_for_GSEA"), "Select database", choices = c("GO", "KEGG", "Reactome", "Msigdb"), selected = "GO", multiple = TRUE),
      conditionalPanel(condition = paste0("input['",ns("type_for_GSEA"),"'] &&","input['",ns("type_for_GSEA"),"'].indexOf('Msigdb') > -1"), #&& Msigdb in the selectize of type_for_GSEA
                       tagList(checkboxGroupInput(ns("Msigdb_opt"), label = "Msigdb subsets",
                                                  choices = c("c2.cp (Canonical pathways)" = "C2.CP",
                                                              "c2.cgp (chemical and genetic perturbations)" = "C2.CGP",
                                                              "c2.cp.biocarta" = "C2.BIOCARTA",
                                                              "c3.all (motif gene sets)" = "C3", "c3.tft (transcription factor targets)" = "C3.TFT",
                                                              "c6.all (oncogenic signatures)" = "C6", "c7.all (immunologic signatures)" = "C7",
                                                              "h.all.v6.1 (hallmark gene sets)" = "H"
                                                              # , "other collections" = "other"
                                                              ) ),
                               shinyBS::bsTooltip(
                                 ns("Msigdb_opt"),
                                 "Selected geneset collection in Msigdb",
                                 "right",
                                 options = list(container = "body")
                               )),
                       # conditionalPanel(condition = paste0("input['",ns("Msigdb_opt"),"'] &&","input['",ns("Msigdb_opt"),"'].indexOf('other') > -1"), #&& Msigdb in the selectize of type_for_GSEA
                       #                  div(style = "word-wrap:break-word;overflow-wrap: break-word;width:100%;overflow-x: auto;",
                       #                      shinyTree(ns("Msigdb_opt2"), checkbox =T,theme="proton", themeIcons = FALSE, themeDots = T),
                       #                      # verbatimTextOutput(ns("shinyTreeselect")),
                       #                      # verbatimTextOutput(ns("shinyTreeselect2")),
                       #                      # verbatimTextOutput(ns("shinyTreeselect3")),
                       #                  )
                       # ),
                       p(a("  MsigDB link",
                           href = "http://www.gsea-msigdb.org/gsea/msigdb/index.jsp",
                           target="_blank"))
      ),

      uiOutput(ns("organism_for_GSEA")),
      uiOutput(ns("df_with_lg2fc")),
      actionButton(ns("analyze_for_GSEA"), "Analyze"),
      shinyBS::bsTooltip(ns("topn"), "Reserve the topn genes in the list(sorted by absolute values of L2FC)  to GSEA", "top", options = list(container = "body")),
      shinyBS::bsTooltip(ns("import_from_text_paste"), "Import genes by paste gene list", "top", options = list(container = "body")),
      shinyBS::bsTooltip(ns("genelist_tool_for_GSEA"), "Choose the gene list from gene list tool options", "top", options = list(container = "body")),
      shinyBS::bsTooltip(ns("text_input_for_GSEA"), "Paste your gene list here", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("type_for_GSEA"), "Select database, supporting GO, KEGG, reactome pathway and MsigDB database", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("organism_for_GSEA"), "Select the organism. Note that: GO and KEGG support 21 organisms, while reactome supporting 7 organisms. When your selected database includes reactome, it will only appears 7 organisms those are supported by reactome. The others support 3 organisms", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("df_with_lg2fc"), "Whether your gene list with log2 fold change", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("analyze_for_GSEA"), "Click on it to analyze", "right", options = list(container = "body"))
    )

  )
}

#* GSEA body ----
GSEA_body_mod <- function(id, label = "GSEA_body") {
  ns = NS(id)
  tagList(
    mainPanel(
      tabName = id,
      width = 9,
      fluidRow(tags$head(tags$style(type="text/css", "
                        #loadmessage {
                        top: 0px; left: 0px;
                        width: 100%; padding: 5px 0px 5px 0px;
                        text-align: center; font-weight: bold;
                        font-size: 100%; color: #000000;
                        background-color: #FFC1C1; z-index: 105;}")), ## 提示条的样式
               conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                tags$div("calculating...please wait...",id="loadmessage"))
               ),
      uiOutput(ns("GSEA_tabsetPanel"))
    )
  )
}



#* GSEA server ----
GSEA_server_module <- function(id, Omics_res) {
  ID = id
  print("GSEA server loading")
  moduleServer(
    ID,
    function(input, output, session) {
      library(enrichplot)
      library(ggplot2)

      library(DOSE)
      library(GO.db)
      library(dplyr)
      library(clusterProfiler)
      library(clusterProfiler.dplyr)
      library(ReactomePA)

      ns = session$ns

      annoSpecies_df <- DEP2:::annoSpecies_df()
      # annoSpecies_df_for_reactome <- readRDS("annoSpecies_df_for_reactome.rds")

      #** organism selection ----
      output$organism_for_GSEA <- renderUI({

        if(is.element("Msigdb", input$type_for_GSEA)){
          selectizeInput(ns("organism_for_GSEA"), "Select organism", choices=c("Human", "Mouse", "Rat"), selected="Human")
        } else {
          if(is.element("Reactome", input$type_for_GSEA)){
            selectizeInput(ns("organism_for_GSEA"), "Select organism", choices=annoSpecies_df$species, selected="Human")
          } else {
            selectizeInput(ns("organism_for_GSEA"), "Select organism", choices=annoSpecies_df$species, selected="Human")
          }
        }
      })

      #** drag-drop input ----
      output$GSEA_drag_drop <- renderUI({
        tagList(
          uiOutput(ns("ui_dragzone")),
          uiOutput(ns("ui_dropzone"))
        )
      })

      Omics_res_list <- reactive({
        Omics_res = reactiveValuesToList(Omics_res)
        omictype <- names(Omics_res) %>% sapply(function(x){strsplit(x,"_")[[1]][1]})
        Omics_res <- Omics_res[which(!omictype %in% c("Timecourse","PTM"))]
        check = Omics_res %>% lapply(., function(x){
          return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
        }) %>% unlist


        Omics_res = Omics_res[check] ## remove
      })

      #*** dragzone UI from genelist_tool ----
      output$ui_dragzone <- renderUI({
        # if(input$import_from_text_paste) {
        if(is.null(Omics_res_list())||length(Omics_res_list()) < 1){
          mes <- HTML("No  omics data, please perform differential analyze in omics modules.")
        }else{
          mes <- NULL
        }
        tagList(
          h5("Dragzone"),
          mes,
          dragZone(ns("dragzone"),
                   choices = names(Omics_res_list()),
                   style = "text-align: center;
                      font-size: 12px;
                      background-color:   #F8F8FF;
                      border-radius: 6px;
                      min-width: 80px;
                      color:  #6495ED;"
          )
        )

        # }
      })

      #*** drop UI from genelist_tool ----
      exctract_genelist_UI_GSEA <- exctract_genelist_UI_suffixFun("GSEAgenelist")

      drop_choose <- reactive({
        x = lapply(names(Omics_res_list()), exctract_genelist_UI_GSEA)
        names(x) = names(Omics_res_list())
        return(x <- x)
      })
      output$ui_dropzone <- renderUI({
        print("render dropZoneInput")
        drop_choose <- isolate(drop_choose())
        drop_input <- dropZoneInput(ns("GSEA_dropzone"),
                                    choices = drop_choose,
                                    # list(one = exctract_genelist_UI("one"),
                                    #              two = exctract_genelist_UI("two")),
                                    # presets = list(values = "one-ds-1"),
                                    selectable = F,
                                    multivalued = F,
                                    server = exctract_genelist_UI_GSEA,
                                    style = "height: 450px;
                                    overflow-y: auto;
                                    color: #6495ED;
                                    font-color: black;
                                    border-radius: 6px;"
        )
        dropZoneServer(session, "GSEA_dropzone", exctract_genelist_UI_GSEA)
        return(tagList(h5("Dropzone"), drop_input))
      })
      # print("render GSEA_dropzone finished")

      drop_genelist_modules <- reactiveValues() # drop_genelist_modules store the reactive genelist value dou
      observeEvent(input$GSEA_dropzone, {
        check = reactiveValuesToList(drop_genelist_modules) %>% lapply(., function(x){
          return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
        }) %>% unlist
        purrr::map(
          setdiff(input$GSEA_dropzone, names(drop_genelist_modules)[check]),
          ~ { drop_genelist_modules[[.]] <- exctract_genelist_Server_GSEA(id = ., Omics_Serv_res = Omics_res_list, suffix = "GSEAgenelist") }
        )
        cat(paste0("outside ns is ",session$ns(""),"\n"))
        purrr::map(
          setdiff(names(drop_genelist_modules)[check], input$GSEA_dropzone),
          ~ { drop_genelist_modules[[.]] <- NULL })
        print("exctract_genelist_Server_GSEA load")

      }, priority = 1, ignoreNULL = F)


      #** Msigdb_options shinyTree input ----
      # output$Msigdb_opt2 <- renderTree({
      #   list("C1 (Gene sets corresponding to human chromosome and cytogenetic band)" = "C1",
      #        "C2 (Gene sets curated from online pathway databases and the biomedical literature)" =
      #          list("CGP (chemical and genetic perturbations)" = "CGP",
      #               "CP (Canonical pathways)" =
      #                 list("CP:BIOCARTA (Gene sets from BioCarta pathway database)" = "CP:BIOCARTA",
      #                      "CP:PID (Gene sets from PID pathway database)" = "CP:PID",
      #                      "CP:WIKIPATHWAYS (Gene sets from WikiPathways database)" = "CP:WIKIPATHWAYS")
      #          ),
      #        "C3 (Gene sets representing potential targets of regulation by transcription factors or microRNAs)" =
      #          list("MIR (microRNA targets)" =
      #                 list("MIR:MIRDB (Human miRNA targets as catalogued by miRDB)" = "MIR:MIRDB",
      #                      "MIR:MIR_Legacy (Older gene sets that sharing putative target sites of human miRNA in 3'-UTRs)" = "MIR:MIR_Legacy")
      #          ),
      #        "C4 (Computational gene sets defined by mining large collections of cancer-oriented microarray data)" =
      #          list("CGN (cancer gene neighborhoods by Subramanian, Tamayo et al. 2005.)" = "CGN",
      #               "CM (cancer modules Segal et al. 2004.)" = "CM"),
      #        "C5 (ontology gene sets)" =
      #          list("HPO (Human Phenotype Ontology)" = "HPO"),
      #        "C8 (Gene sets curated cluster markers for cell types identified in single-cell sequencing studies of human tissue)" = "C8")
      #
      # })

      # Msigdb_opt_selected2 <- reactive({
      #   tree_sel <- get_selected(input$Msigdb_opt2 , format = "slices")
      #   # shinytreeInput_treatment(tree_sel %>% unlist() %>% names)
      #   tree_sel_save2 <- shinytreeInput_treatment(tree_sel %>% unlist() %>% names)
      # })

      #** input handle ----
      # genelist <- reactive({ strsplit(input$text_input_for_GSEA,'\n')[[1]] })
      gene_df <- reactive({
        if(input$import_from_text_paste){

          temp <- input$text_input_for_GSEA
          # check the nrow of text input
          check <- strsplit(temp,'\n')[[1]]
          if(is.null(check) || length(check) == 0)
            return(NULL)

          genelist <- fread(text = temp, header = F) %>% as.data.frame()
          if(ncol(genelist) < 2){
            sendSweetAlert(
              session = shiny::getDefaultReactiveDomain(),
              title = "warning !",
              text = paste0("Please check! GSEA need text input with two columns. First column is name, the second is lf2c."),
              type = "warning"
            )
            return(NULL)
          }
          colnames(genelist)[1:2] = c("name","fc")
          genelist$name <-  sapply(genelist$name,the1stname)
          genelist$fc = as.numeric(genelist$fc)
        }else{
          print("aa1")
          genelist <- drop_genelist_modules %>% reactiveValuesToList() # change genelist reactivevalues into List
          if(is.null(genelist) || length(genelist) < 1)
            return(NULL)

          genelist = genelist[[1]]()
          genelist <- genelist[,c("name","L2FC")]
          genelist$name <-  sapply(genelist$name,the1stname)
          colnames(genelist) = c("name","fc")
        }
        genelist <- as.data.frame(genelist)

        return(genelist)
      })

      GSEA_organism = reactive(input$organism_for_GSEA)

      gene_list <- reactive({
        gene_df = gene_df()
        if(is.null(gene_df) || nrow(gene_df) == 0) return(NULL)
        # print("a1")
        # shiny::validate(
        #   need(nrow(gene_df) != 0, message = "No genes meet your requirements, and can not do the GO analysis")
        # )

        print("a2")
        pkg = annoSpecies_df$pkg[annoSpecies_df$species == GSEA_organism()]
        # req( try(require(pkg, character.only = TRUE)) )
        if(!require(pkg, character.only = TRUE)){
          sendSweetAlert(
            session = shiny::getDefaultReactiveDomain(),
            title = "warning !",
            text = paste0("It needs the annotation data package ",pkg," for the species '", ORA_organism(),"'. But it is not installed!",
                          "Please install it first."),
            type = "warning"
          )
          return(NULL)
        }

        print("a3")
        orgDB <- get(pkg)
        gene_id_table <- DEP2:::map_to_entrezid(as.character(gene_df$name), orgDB = orgDB)
        ids <- gene_id_table %>% tibble::rownames_to_column() %>%
          dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID)

        ids <- ids %>% inner_join(., as.data.frame(gene_df), by = "name")

        ids <- ids[!is.na(ids$ENTREZID) & !is.na(ids$fc), ]

        de = ids$fc
        names(de) = unlist(ids$ENTREZID)
        de = sort(de, decreasing = T)

        if((!is.na(input$topn)) && (!is.null(input$topn)) && (input$topn != "")){
          de = de[order(abs(de),decreasing = T)[1:min(input$topn,length(de))]]
          de = sort(de, decreasing = T)
        }
        de = de[which(!duplicated(names(de)))] ## remove duplicated genes, only retain the first one(with a highest l2fc)

        print("a4")
        return(de)
      })



      #* GSEA_tabsetPanel display ----
      output$GSEA_tabsetPanel <- renderUI({
        tabsetPanel(id = session$ns("GSEA"),
                    GSEA_GO_UI(id = session$ns("GO"), tab_name = "GO"),
                    GSEA_other_UI(id = session$ns("KEGG"), tab_name = "KEGG"),
                    GSEA_other_UI(id = session$ns("Reactome") , tab_name = "Reactome"),
                    GSEA_Msigdb_UI(id = session$ns("Msigdb"))
        )
      })

      observeEvent(input$GSEA,{
        hideTab(inputId = "GSEA", target = "KEGG"); hideTab(inputId = "GSEA", target = "Reactome"); hideTab(inputId = "GSEA", target = "Msigdb")
      },once = T, ignoreNULL = T, ignoreInit = F)

      observeEvent(input$type_for_GSEA,{
        if(length(input$type_for_GSEA) > 0) {
          for(i in input$type_for_GSEA) showTab(inputId = "GSEA", target = i)
        }
        GSEA_database_unselected <- setdiff(c("GO","KEGG","Reactome","Msigdb"),input$type_for_GSEA)
        if(length(GSEA_database_unselected) > 0) {
          for(i in GSEA_database_unselected) hideTab(inputId = "GSEA", target = i)
        }
      }, ignoreNULL = T, ignoreInit = T)


      #* Event of analyze ----
      GSEA_database_selected <- reactive(input$type_for_GSEA)

      observeEvent(input$analyze_for_GSEA,{
        # warning if select no database
        type_for_GSEA_Pass1 <- ifelse((is.null(input$type_for_GSEA) || input$type_for_GSEA == ""), F, T)
        shinyFeedback::feedbackWarning("type_for_GSEA", !type_for_GSEA_Pass1, "Please select at least one database")
        req(type_for_GSEA_Pass1)

        # print(gene_name())
        print(GSEA_database_selected())
        print(GSEA_organism())



        GSEA_organism = GSEA_organism()

        # ** import enrichment analyze module,  ----
        if("GO" %in% GSEA_database_selected()){
          print("start gseaGO")
          GSEA_GO_server_module(id = "GO",
                               # gene_df = gene_df,
                               gene_list = gene_list(), organism_for_GSEA = GSEA_organism,
                               annoSpecies_df = annoSpecies_df)
        }

        if("KEGG" %in% GSEA_database_selected()){
          print("loading KEGG server")
          GSEA_other_server_module(id = "KEGG",
                                  # gene_df = gene_df,
                                  gene_list = gene_list(),  organism_for_GSEA = GSEA_organism,
                                  annoSpecies_df = annoSpecies_df, type = "kegg")
          print("loaded KEGG server")
        }

        if("Reactome" %in% GSEA_database_selected()){
          print("loading Reactome server")
          GSEA_other_server_module(id = "Reactome",
                                  # gene_df = gene_df,
                                  gene_list = gene_list(),  organism_for_GSEA = GSEA_organism,
                                  annoSpecies_df = annoSpecies_df, type = "reactome")
          print("loaded Reactome server")
        }

        if("Msigdb" %in% GSEA_database_selected()){
          GSEA_Msigdb_server_module(id = "Msigdb",
                                   # gene_df = gene_df,
                                   gene_list = gene_list(),  organism_for_GSEA = GSEA_organism,
                                   annoSpecies_df = annoSpecies_df,
                                   Msigdb_selection = input$Msigdb_opt, Msigdb_selection2 = NULL
                                     # Msigdb_opt_selected2()
                                   ) # id here should not add the namespace prefix. Don't use ns("Msigdb")
        }

      })

    }
  )
}




# *functions for drag-drop input in GSEA ------------------------------------------
# Module: Server exctract genelist item from choosed omics
# ID transmit from the dragInput choose. ns is (dragInput id)-(dragInput choose). Omics_serv_res to render
exctract_genelist_Server_GSEA <- function(id, ID, Omics_Serv_res, suffix = ""){
  ID = id
  id <- paste0(id, "-", suffix)
  moduleServer(
    id,
    function(input,output,session){
      ns = session$ns
      # omics_info = get(".OmicsModule",envir = get_DEPglobal_env())
      cat(paste0("the exctract_genelist_Server_GSEA ns is ", ns(''), "\n"))
      cat(paste0("the exctract_genelist_Server_GSEA ID is ", ID, "\n"))

      # OmicsModules = get(".OmicsModule",envir = get_DEPglobal_env())
      # sta
      omictype <- strsplit(ID,"_")[[1]][1]
      if( omictype == "Proteomepg" | omictype == "Proteomepep"){
        ## genelistUI in proteome
        output$genelistUI <- renderUI({
          tagList(
            selectizeInput(ns("contrast"),
                           "Contrast",
                           choices = get_contrast(test_res()), multiple = F),
            uiOutput(ns("table_bttn"))
          )
        })
      }else if( omictype == "RNAseq"){
        ## genelistUI in RNAseq
        output$genelistUI <- renderUI({
          tagList(
            selectizeInput(ns("contrast"),
                           "Contrast",
                           choices = grep("_p.adj$",colnames(test_res()@test_result),value = T) %>% gsub("_p.adj", "", .), multiple = F),
            uiOutput(ns("table_bttn"))
          )
        })
      }

      # output$restest <- renderDataTable({
      #   (sig_res())
      # }, options = list(pageLength = 10, scrollX = T))

      output$table_bttn <- renderUI({
        # omictype <- strsplit(ID,"_")[[1]][1]
        theres <- test_res()
        thegenelist <- genelist()

        if((omictype == "RNAseq") && (nrow(theres@geneinfo) == 0) ){ ## check id transform of RNAseq data
          return(
            tagList(
              bsButton(ns("table_bttn"),label = "ID Transformation is necessary for enrichment analysis",
                       style = "warning", block = TRUE, type = "action", disabled = T)
            )
          )
        }else if(!is.null(thegenelist) &&nrow(thegenelist) > 0){
          return(
            tagList(
              bsButton(ns("table_bttn"),label = paste0("Gene list: ",nrow(thegenelist)),
                       style = "primary", block = TRUE, type = "action", disabled = F)
            )
          )
        }else{
          return(
            tagList(
              bsButton(ns("table_bttn"),label = paste0("Can't extract gene list"),
                       style = "warning", block = TRUE, type = "action", disabled = T)
            )
          )
        }
      })


      omics_id <- strsplit(ID, split = "-ds-")[[1]][1]
      test_res <-  reactive({
        test_res <- Omics_Serv_res()[[omics_id]]()
        test_res
      })

      iv1 <- InputValidator$new()
      iv1$add_rule("pcutoff", sv_between(0.00001, 0.2))
      iv1$add_rule("l2fc", sv_between(0, 10))
      iv1$enable()


      if(omictype == "Proteomepg" | omictype == "Proteomepep"){

        genelist <- reactive({
          req(test_res)
          req(!is.null(input$contrast))  ## check input UI initialization was finished
          dep_res <- rowData(test_res()) %>% as.data.frame()
          genelist = dep_res[,c("name", "ID", paste0(input$contrast,c("_diff","_p.adj")))]
          colnames(genelist) = c("name", "ID", "L2FC", "padj")

          return(genelist)
        })
      }else if(omictype == "RNAseq"){
        genelist <- reactive({
          the_res = test_res()@test_result %>% as.data.frame()
          gene_info = test_res()@geneinfo

          if(nrow(gene_info) == 0)
            return(NULL)

          thecontrast = input$contrast
          genelist <- data.frame(name = gene_info$SYMBOL,L2FC =  the_res[, paste0(thecontrast,"_diff")], padj = the_res[paste0(thecontrast,"_p.adj")])

          return(genelist)
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
        genelist()
      }, options = list(pageLength = 10, scrollX = T))

      return(genelist)
    })
}

