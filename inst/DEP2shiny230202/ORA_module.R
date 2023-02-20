source("ORA_submodule.R")
source("source_enrich_ana.R")
# library(enrichplot)
# library(ggplot2)

# library(DOSE)
# library(GO.db)
# library(dplyr)
# library(clusterProfiler.dplyr)
# library(ReactomePA)

#* ORA UI ----
ORA_UI <- function(id){
  tagList(
    sidebarLayout(
      ORA_sidebar_mod(id = id),
      ORA_body_mod(id = id)
    ))
}

#* ORA sidebar ----
ORA_sidebar_mod <-  function(id,label="ORA_sidabar"){
  ns = NS(id)
  cat(paste0("ORA_sidebar ns is",ns(""),"\n"))

  tagList(
    sidebarPanel(
      width = 3,
      h3(id),
      checkboxInput(ns("import_from_text_paste"), "import from text paste", value = FALSE),

      conditionalPanel(condition = paste0("!input['",ns("import_from_text_paste"),"']"),#&& "input.filt_num_cutoff_ptm !== ''"
                       tagList(
                         uiOutput(ns("ORA_drag_drop"))
                         # uiOutput(ns("ORA_dragzone")),
                         # uiOutput(ns("ui_dropzone"))
                       )
      ),
      conditionalPanel(condition = paste0("input['",ns("import_from_text_paste"),"']"),#&& "input.filt_num_cutoff_ptm !== ''"
                       textAreaInput(inputId = ns("text_input_for_ORA"), label = "Please paste your gene list", resize = "none",
                                     placeholder = "TP53\nPTEN\n\nor\n\nTP53  6.21\nPTEN  -1.53", rows = 8, width = "100%")
      ),

      # uiOutput(ns("text_input_for_ORA")),
      selectizeInput(ns("type_for_ORA"), "Select database", choices = c("GO", "KEGG", "Reactome", "Msigdb"), selected = "GO", multiple = TRUE),
      conditionalPanel(condition = paste0("input['",ns("type_for_ORA"),"'] &&","input['",ns("type_for_ORA"),"'].indexOf('Msigdb') > -1"), #&& Msigdb in the selectize of type_for_ORA
                       tagList(checkboxGroupInput(ns("Msigdb_opt"), label = "Msigdb subsets",
                                    choices = c("c2.cp (Canonical pathways )" = "C2.CP",
                                                "c2.cgp (chemical and genetic perturbations)" = "C2.CGP",
                                                "c2.cp.biocarta" = "C2.BIOCARTA",
                                                "c3.all (motif gene sets)" = "C3", "c3.tft (transcription factor targets)" = "C3.TFT",
                                                "c6.all (oncogenic signatures)" = "C6", "c7.all (immunologic signatures)" = "C7",
                                                "h.all.v6.1 (hallmark gene sets)" = "H", "other collections" = "other") ),
                               shinyBS::bsTooltip(
                                 ns("Msigdb_opt"),
                                 "Selected geneset collection in Msigdb",
                                 "right",
                                 options = list(container = "body")
                               )),
                       # conditionalPanel(condition = paste0("input['",ns("Msigdb_opt"),"'] &&","input['",ns("Msigdb_opt"),"'].indexOf('other') > -1"), #&& Msigdb in the selectize of type_for_ORA
                       #                  div(style = "word-wrap:break-word;overflow-wrap: break-word;width:100%;overflow-x: auto;",
                       #                      shinyTree(ns("Msigdb_opt2"), checkbox =T,theme="proton", themeIcons = FALSE, themeDots = T),
                       #                      verbatimTextOutput(ns("shinyTreeselect")),
                       #                      verbatimTextOutput(ns("shinyTreeselect2")),
                       #                      verbatimTextOutput(ns("shinyTreeselect3")),
                       #                  )
                       # ),
                       p(a("  MsigDB link",
                           href = "http://www.gsea-msigdb.org/gsea/msigdb/index.jsp",
                           target="_blank"))
      ),

      uiOutput(ns("organism_for_ORA")),
      uiOutput(ns("df_with_lg2fc")),
      actionButton(ns("analyze_for_ORA"), "Analyze"),
      shinyBS::bsTooltip(ns("import_from_text_paste"), "Import genes by paste gene list", "top", options = list(container = "body")),
      shinyBS::bsTooltip(ns("genelist_tool_for_ORA"), "Choose the gene list from gene list tool options", "top", options = list(container = "body")),
      shinyBS::bsTooltip(ns("text_input_for_ORA"), "Paste your gene list here", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("type_for_ORA"), "Select database, supporting GO, KEGG, reactome pathway and MsigDB database", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("organism_for_ORA"), "Select the organism. Note that: GO and KEGG support 21 organisms, while reactome supporting 7 organisms. When your selected database includes reactome, it will only appears 7 organisms those are supported by reactome. The others support 3 organisms", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("df_with_lg2fc"), "Whether your gene list with log2 fold change", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("analyze_for_ORA"), "Click on it to analyze", "right", options = list(container = "body"))
    )

  )}



#* ORA body ----
ORA_body_mod <- function(id, label = "ORA_body") {
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
      uiOutput(ns("ORA_tabsetPanel"))
    )

  )
}

#* ORA server ----
ORA_server_module2 <- function(id, Omics_res) {
  ID = id
  print("ORA server loading")
  moduleServer(
    ID,
    function(input, output, session) {
      library(enrichplot)
      library(ggplot2)

      library(DOSE)
      library(GO.db)
      library(dplyr)
      library(clusterProfiler.dplyr)
      library(ReactomePA)

      ns = session$ns

      annoSpecies_df <- DEP2:::annoSpecies_df()
      # annoSpecies_df_for_reactome <- DEP2:::annoSpecies_df()

      #* organism selection ----
      output$organism_for_ORA <- renderUI({

        if(is.element("Msigdb", input$type_for_ORA)){
          selectizeInput(ns("organism_for_ORA"), "Select organism", choices=c("Human", "Mouse", "Rat"), selected="Human")
        } else {
          if(is.element("Reactome", input$type_for_ORA)){
            selectizeInput(ns("organism_for_ORA"), "Select organism", choices=annoSpecies_df$species, selected="Human") ## Now reactome and GO have the same selection
          } else {
            selectizeInput(ns("organism_for_ORA"), "Select organism", choices=annoSpecies_df$species, selected="Human")
          }
        }
      })

      #* drag-drop input ----
      output$ORA_drag_drop <- renderUI({
        tagList(
          uiOutput(ns("ui_dragzone")),
          uiOutput(ns("ui_dropzone"))
        )
      })
      # Omics_res_list <- reactive(reactiveValuesToList(Omics_res))
      Omics_res_list <- reactive({
        Omics_res = reactiveValuesToList(Omics_res)
        # Omics_res = reactiveValuesToList(omics_Servers_saved)
        check = Omics_res %>% lapply(., function(x){
          return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
        }) %>% unlist
        Omics_res = Omics_res[check] ## remove
      })

      #** dragzone UI from genelist_tool ----
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

      #** dropzone UI from genelist_tool ----
      exctract_genelist_UI_ORA <- exctract_genelist_UI_suffixFun("ORAgenelist")

      drop_choose <- reactive({
        x = lapply(names(Omics_res_list()), exctract_genelist_UI_ORA)
        names(x) = names(Omics_res_list())
        return(x <- x)
      })
      output$ui_dropzone <- renderUI({
        print("render dropZoneInput")
        drop_input <- dropZoneInput(ns("ORA_dropzone"),
                                    choices = drop_choose(),
                                    # list(one = exctract_genelist_UI("one"),
                                    #              two = exctract_genelist_UI("two")),
                                    # presets = list(values = "one-ds-1"),
                                    selectable = F,
                                    multivalued = F,
                                    server = exctract_genelist_UI_ORA,
                                    style = "height: 450px;
                                    overflow-y: auto;
                                    color: #6495ED;
                                    font-color: black;
                                    border-radius: 6px;"
        )
        dropZoneServer(session, "ORA_dropzone", exctract_genelist_UI_ORA)
        return(tagList(h5("Dropzone"), drop_input))
      })
      # print("render ORA_dropzone finished")

      #** dropzone server, return the significant data.frame ----
      drop_genelist_modules <- reactiveValues()
      observeEvent(input$ORA_dropzone, {
        # drop_genelist_modules_list = reactiveValuesToList(drop_genelist_modules)
        check = reactiveValuesToList(drop_genelist_modules) %>% lapply(., function(x){
          return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
        }) %>% unlist
        purrr::map(
          setdiff(input$ORA_dropzone, names(drop_genelist_modules)[check]),
          ~ { drop_genelist_modules[[.]] <- exctract_genelist_Server(id = ., Omics_Serv_res = Omics_res_list, suffix = "ORAgenelist") }
        )
        cat(paste0("outside ns is ",session$ns(""),"\n"))
        purrr::map(
          setdiff(names(drop_genelist_modules)[check], input$ORA_dropzone),
          ~ { drop_genelist_modules[[.]] <- NULL })
        print("exctract_genelist_Server load")
      }, priority = 1, ignoreNULL = F)

      #* Msigdb_options shinyTree input ----
      # output$Msigdb_opt2 <- renderTree({
      #   list("C1 (Gene sets corresponding to human chromosome and cytogenetic band)" = "C1",
      #        "C2 (Gene sets curated from online pathway databases and the biomedical literature)" =
      #          list("CGP (chemical and genetic perturbations)" = "CGP",
      #               "CP (Canonical pathways)" =
      #                 list("CP:BIOCARTA (Gene sets from BioCarta pathway database)" = "CP:BIOCARTA",
      #                      "CP:PID (Gene sets from PID pathway database)" = "CP:PID",
      #                      "CP:WIKIPATHWAYS (Gene sets from WikiPathways database)" = "CP:WIKIPATHWAYS")
      #               ),
      #        "C3 (Gene sets representing potential targets of regulation by transcription factors or microRNAs)" =
      #          list("MIR (microRNA targets)" =
      #                 list("MIR:MIRDB (Human miRNA targets as catalogued by miRDB)" = "MIR:MIRDB",
      #                      "MIR:MIR_Legacy (Older gene sets that sharing putative target sites of human miRNA in 3'-UTRs)" = "MIR:MIR_Legacy")
      #               ),
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
      # })

      #* input handle ----
      # genelist <- reactive({ strsplit(input$text_input_for_ORA,'\n')[[1]] })
      gene_name <- reactive({
        if(input$import_from_text_paste){
          temp <- input$text_input_for_ORA
          # check the nrow of text input
          check <- strsplit(temp,'\n')[[1]]
          if(is.null(check) || length(check) == 0)
            return(NULL)
          # if(is.null(temp)) return(NULL)
          # genelist <- strsplit(input$text_input_for_ORA,'\n')[[1]]
          # if(nrow(genelist) == 0) return(NULL)
          # gene_name <-unlist(strsplit(genelist[,1],";")[[1]])

          genelist <- fread(text = temp, header = F) %>% as.data.frame()

          if(nrow(genelist) == 0) return(NULL)
          gene_name <- genelist[,1] %>% sapply(., the1stname) %>% unique()

        }else{
          genelist <- drop_genelist_modules %>% reactiveValuesToList()
          if(length(genelist)>0){
            gene_name = vector()
            for(i in 1:length(genelist)){
              theres <- genelist[[i]]()
              if(is.null(theres)||nrow(theres) < 1){
                gene_name = gene_name; next
              }
              if("gene_name" %in% colnames(theres)){
                gene_name <- c(gene_name, theres$gene_name)
              }else {gene_name <- c(gene_name, theres$name)}
            }
          }else{
            gene_name<- NULL
          }
        }

        if(is.null(gene_name)||length(gene_name) < 1) return(NULL) ## when no sig res (no gene as input)
        gene_name = data.frame(name = gene_name %>% DEP2:::rm_digit_end()  %>% unique()) %>%
          dplyr::filter(!is.na(name))

        return(gene_name)
      })

      ORA_organism = reactive(input$organism_for_ORA)

      gene_id <- reactive({
        gene_df = gene_name()
        if(is.null(gene_df) || nrow(gene_df) == 0) return(NULL)
        # shiny::validate(
        #   need(nrow(gene_df) != 0, message = "No genes meet your requirements, and can not do the GO analysis")
        # )
        pkg = annoSpecies_df$pkg[annoSpecies_df$species == ORA_organism()]
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
        print("loading species annotation data")
        orgDB <- get(pkg)
        gene_id_table <- DEP2:::map_to_entrezid(as.character(gene_df$name), orgDB = orgDB)
        ids <- gene_id_table %>% tibble::rownames_to_column() %>%
          dplyr::rename(., name = rowname, ENTREZID = id) %>% dplyr::select(name, ENTREZID) ## 2 cols table: name(from) & ENTREZID(to)
        # %>% filter(!is.na(ENTREZID)) %>% filter(!duplicated(ENTREZID))

        # ids = my_to_entrezid(orgDB = orgDB, gene = as.character(gene_df$name))
        ids <- ids %>% inner_join(., gene_df, by = "name") %>% filter(!is.na(ENTREZID))
        ids = unlist(ids$ENTREZID)
        names(ids) = ids
        print("Mapping ENTREZID finished!")
        return(ids)
      })



      #* ORA_tabsetPanel display ----
      output$ORA_tabsetPanel <- renderUI({
        tabsetPanel(id = session$ns("ORA"),
                    ORA_GO_UI(session$ns("GO")),
                    ORA_other_UI(id = session$ns("KEGG"), tab_name = "KEGG"),
                    ORA_other_UI(id = session$ns("Reactome") , tab_name = "Reactome"),
                    ORA_Msigdb_UI(id = session$ns("Msigdb"))
        )
      })

      observeEvent(input$ORA,{
        hideTab(inputId = "ORA", target = "KEGG"); hideTab(inputId = "ORA", target = "Reactome"); hideTab(inputId = "ORA", target = "Msigdb")

      },once = T, ignoreNULL = T, ignoreInit = F)

      observeEvent(input$type_for_ORA,{
        if(length(input$type_for_ORA) > 0) {
          for(i in input$type_for_ORA) showTab(inputId = "ORA", target = i)
        }
        ORA_database_unselected <- setdiff(c("GO","KEGG","Reactome","Msigdb"),input$type_for_ORA)
        if(length(ORA_database_unselected) > 0) {
          for(i in ORA_database_unselected) hideTab(inputId = "ORA", target = i)
        }
      }, ignoreNULL = T, ignoreInit = T)

      #* Event of analyze ----
      ORA_database_selected <- reactive(input$type_for_ORA)

      observeEvent(input$analyze_for_ORA,{
        # warning if select no database
        type_for_ORA_Pass1 <- ifelse((is.null(input$type_for_ORA) || input$type_for_ORA == ""), F, T)
        shinyFeedback::feedbackWarning("type_for_ORA", !type_for_ORA_Pass1, "Please select at least one database")
        req(type_for_ORA_Pass1)

        # print(gene_name())
        print(ORA_database_selected())
        print(ORA_organism())



        ORA_organism = ORA_organism()

        # ** import enrichment analyze module,  ----
        if("GO" %in% ORA_database_selected()){
          ORA_GO_server_module(id = "GO",
                               # gene_df = gene_df,
                               gene_id = gene_id(), organism_for_ORA = ORA_organism,
                               annoSpecies_df = annoSpecies_df)
        }

        if("KEGG" %in% ORA_database_selected()){
          print("loading KEGG server")
          ORA_other_server_module(id = "KEGG",
                                  # gene_df = gene_df,
                                  gene_id = gene_id(),  organism_for_ORA = ORA_organism,
                                  annoSpecies_df = annoSpecies_df, type = "kegg")
          print("loaded KEGG server")
        }

        if("Reactome" %in% ORA_database_selected()){
          print("loading Reactome server")
          ORA_other_server_module(id = "Reactome",
                                  # gene_df = gene_df,
                                  gene_id = gene_id(),  organism_for_ORA = ORA_organism,
                                  annoSpecies_df = annoSpecies_df_for_reactome, type = "reactome")
          print("loaded Reactome server")
        }

        if("Msigdb" %in% ORA_database_selected()){
          ORA_Msigdb_server_module(id = "Msigdb",
                                   # gene_df = gene_df,
                                   gene_id = gene_id(),  organism_for_ORA = ORA_organism,
                                   annoSpecies_df = annoSpecies_df,
                                   Msigdb_selection = input$Msigdb_opt, Msigdb_selection2 = NULL
                                     # Msigdb_opt_selected2()
                                   ) # id here should not add the namespace prefix. Don't use ns("Msigdb")
        }

      })
    })
  }











