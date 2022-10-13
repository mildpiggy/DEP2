
#* PPI UI ----
PPI_UI <- function(id){
  tagList(
    sidebarLayout(
      PPI_sidebar_mod(id = id), ## DEP-ptm sidebar module
      PPI_body_mod(id = id)  ## DEP-ptm pagebody module
    )
  )
}

#* PPI sidebar ----
PPI_sidebar_mod <-  function(id,label="PPI_sidabar"){
  ns = NS(id)
  cat(paste0("ORA_sidebar ns is",ns(""),"\n"))

  tagList(
    sidebarPanel(
      width = 3,
      h3(id),
      # checkboxInput(ns("import_from_genelist_tool_for_ppi"), "Import from gene list tool", value = FALSE),
      # uiOutput(ns("genelist_tool_for_ppi")),
      # uiOutput(ns("input_text")),
      checkboxInput(ns("import_from_text_paste"), "import from text paste", value = FALSE),
      conditionalPanel(condition = paste0("!input['",ns("import_from_text_paste"),"']"),#&& "input.filt_num_cutoff_ptm !== ''"
                       tagList(
                         uiOutput(ns("PPI_drag_drop"))
                       )
      ),
      conditionalPanel(condition = paste0("input['",ns("import_from_text_paste"),"']"),#&& "input.filt_num_cutoff_ptm !== ''"
                       textAreaInput(inputId = ns("text_input_for_PPI"), label = "Please paste your gene list", placeholder = "TP53\nPTEN\n\nor\n\nTP53  6.21\nPTEN  -1.53", rows = 8, width = "100%")
      ),
      selectInput(ns("organism"),
                  "Select organism",
                  choices=
                    c("Anopheles", "Arabidopsis", "Bovine", "Worm", "Canine", "Fly", "Zebrafish", "E coli strain K12", "Chicken", "Human", "Mouse", "Rhesus", "Malaria", "Chimp", "Rat", "Yeast", "Streptomyces coelicolor", "Pig", "Toxoplasma gondii", "Xenopus")
                  ,
                  selected="Human"),
      # actionButton("load_PPIdata","Load Datas First!"),
      helpText("choose organism string data to load", br(),"If it doesn't exist, app will give messages on the body", br(),"please download files from string following them first"),
      selectizeInput(ns("chooseScore"),
                     "Select concerned scores" ,
                     choices=c("combined_score" ,"neighborhood","fusion","cooccurence","coexpression","experimental","database","textmining"),
                     multiple = TRUE,
                     selected = c("combined_score" ,"neighborhood","fusion","cooccurence","coexpression","experimental","database","textmining")),
      fluidRow(
        column(width = 10,
               numericInput(ns("scorecutoff"),
                            "scores cutoff",
                            min = 400, max = 1000, value = 400))),
      actionButton(ns("String_annalysis"),"String It!"),
      # radioButtons("reactome_color",
      #                "colorBy",
      #                c("pvalue", "p.adjust"),
      #                selected = "p.adjust"),
      # h4("\n"),
      tags$hr(),
      tags$style(type="text/css", "#downloadPPITable {background-color:white;color: black;font-family: Source Sans Pro}"),
      tags$style(type="text/css", "#downloadPPInetwork {background-color:white;color: black;font-family: Source Sans Pro}"),
      # h4("\n"),
      uiOutput(ns("download_PPITable")),
      h6(),
      uiOutput(ns("download_PPInetwork")),
      # downloadButton("downloadPPITable", "Save table"),
      # downloadButton("downloadPPInetwork", "Save network")
      # shinyBS::bsTooltip(ns("import_from_genelist_tool_for_ppi"), "Import genes from the gene list tool options", "top", options = list(container = "body")),
      # shinyBS::bsTooltip(ns("genelist_tool_for_ppi"), "Choose the gene list from gene list tool options", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("text_input_for_PPI"), "Paste gene list here", "top", options = list(container = "body")),
      shinyBS::bsTooltip(ns("organism"), "Select the organism", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("chooseScore"), "Select which type of evidence will contribute to the prediction of the score under the active interaction sources", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("scorecutoff"), "Minimum required interaction score. The minimum required interaction score puts a threshold on the confidence score, such that only interaction above this score are included in the predicted network. Confidence limits are as follows: low confidence - 0.15 (or better); medium confidence - 0.4; high confidence - 0.7; highest confidence - 0.9", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("String_annalysis"), "Click on it to analyze", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("download_PPITable"), "Click on it to download the Protein-Protein Interaction table", "right", options = list(container = "body")),
      shinyBS::bsTooltip(ns("download_PPInetwork"), "Click on it to download the PPI(Protein-Protein Interaction Network)", "right", options = list(container = "body"))
    ))
}

#* PPI body ----
PPI_body_mod <- function(id, label = "PPI_body") {
  ns = NS(id)
  tagList(
    mainPanel(
      tabName = id,
      width = 9,
      fluidRow(
        column(width = 12,
               # tabBox(width = 12,
               # tabPanel(

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
               ## ** Network options ----
               fluidRow(height=30,
                        column(width = 3,
                               colourpicker::colourInput(inputId = ns("nodecolor"),
                                                         label = "input node color(like:#2EA9DF)",
                                                         showColour = "both",
                                                         palette = "square",
                                                         value = "#2EA9DF")
                               # textInput("nodecolor", "input node color(like:#2EA9DF)", "#2EA9DF")
                        ),
                        column(width = 3,
                               checkboxInput(ns("changesize"),
                                             "whether change nodesizes by connection",
                                             value = FALSE)),
                        column(width = 3,
                               numericInput(ns("nodesize"),
                                            "basic node size",
                                            min = 10, max = 150, value = 20)),
                        column(width = 3,
                               numericInput(ns("PPIfontsize"),
                                            "fonts size",
                                            min = 10, max = 150, value = 25))
               ),
               fluidRow(height=30,
                        column(width = 3,
                               colourpicker::colourInput(inputId = ns("linecolor"),
                                                         label = "input line color(like:  #ADD8E6)",
                                                         showColour = "both",
                                                         palette = "square",
                                                         value = "#ADD8E6")
                               # textInput("linecolor", "input line color(like:  #ADD8E6)", "#ADD8E6")
                        ),

                        column(width = 3,
                               checkboxInput(ns("changewidth"),
                                             "whether change link width by connection",
                                             value = FALSE)),
                        column(width = 3,
                               numericInput(ns("linewidth"),
                                            "basic line width",
                                            min = 3, max = 80, value = 5)),
                        column(width = 3,
                               selectizeInput(ns("layoutway"),
                                              "choose layout style of network",
                                              choices=c("layout_as_star","layout_components","layout_in_circle", "layout_nicely" ,"layout_on_grid","layout_on_sphere",
                                                        "layout_randomly","layout_with_dh","layout_with_drl" ,"layout_with_fr","layout_with_gem","layout_with_graphopt",
                                                        "layout_with_kk","layout_with_lgl"  ),
                                              selected="layout_components"))
               ),
               fluidRow(height=30,
                        column(width = 2,
                               selectizeInput(ns("nodeshape"),
                                              "choose node shape",
                                              choices=c("square", "triangle", "box", "circle", "dot", "star",
                                                        "ellipse", "database", "text", "diamond"),
                                              selected="dot")),
                        column(width = 2,
                               checkboxInput(ns("highlightkey"),
                                             "highlight key nodes?",
                                             value = TRUE)),
                        column(width = 2,
                               checkboxInput(ns("smoothline"),
                                             "smooth lines",
                                             value = FALSE)),
                        column(width = 2,
                               selectizeInput(ns("smoothtype"),
                                              "type to smooth the line",
                                              choices=c('dynamic', 'continuous', 'discrete', 'diagonalCross', 'straightCross', 'horizontal',
                                                        'vertical', 'curvedCW', 'curvedCCW', 'cubicBezier'),
                                              selected="continuous")),
                        column(width = 3,
                               numericInput(ns("plothigh"),
                                            "plot high",
                                            min = 10, max = 200, value = 75))
               ),
               fluidRow(
                 # visNetworkOutput("network",height=input$plothigh*10)
                 uiOutput(ns("PPI_plot"))
               ),
               fluidRow(
                 DT::dataTableOutput(ns("String_Table"))
               ),
               shinyBS::bsTooltip(ns("nodecolor"), "choose the color of input node", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("changesize"), "whether change nodesizes by connection", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("nodesize"), "basic node size", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("PPIfontsize"), "fonts size", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("linecolor"), "choose the color of input line", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("changewidth"), "whether change link width by connection", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("linewidth"), "basic line width", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("layoutway"), "choose layout style of network", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("nodeshape"), "choose node shape", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("highlightkey"), "highlight key nodes?", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("smoothline"), "smooth lines", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("smoothtype"), "type to smooth the line", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("plothigh"), "the height of the plot region", "top", options = list(container = "body")),
               shinyBS::bsTooltip(ns("String_Table"), "the table of PPI (Protein-Protein Interaction)", "top", options = list(container = "body"))
               # )
               # )

        )
      )

    )

  )
}


#* PPI server ----
PPI_server_module <- function(id, Omics_res) {
  ID = id
  print("PPI server loading")
  moduleServer(
    ID,
    function(input, output, session) {

      ns = session$ns


      #* drag-drop input ----
      # Omics_res_list <- reactive(reactiveValuesToList(Omics_res))
      Omics_res_list <- reactive({
        Omics_res = reactiveValuesToList(Omics_res)
        # Omics_res = reactiveValuesToList(omics_Servers_saved)
        check = Omics_res %>% lapply(., function(x){
          return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
        }) %>% unlist
        Omics_res = Omics_res[check] ## remove
      })

      # Omics_res_list_save <<- Omics_res_list()
      output$PPI_drag_drop <- renderUI({
        validate(need(length(Omics_res_list()) > 0 && !all(sapply(Omics_res_list(), is.null)), message = "There are no omics analysis result. Please do omics analysis in . \n\n"))
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
                   style =
                     "text-align: center;
                      font-size: 12px;
                      background-color:   #F8F8FF;
                      border-radius: 6px;
                      min-width: 80px;
                      color:  #6495ED;"
          )
        )
      })

      #** drop UI from genelist_tool ----
      exctract_genelist_UI_PPI <- exctract_genelist_UI_suffixFun("PPIgenelist")

      drop_choose <- reactive({
        x = lapply(names(Omics_res_list()), exctract_genelist_UI_PPI)
        names(x) = names(Omics_res_list())
        return(x <- x)
      })
      output$ui_dropzone <- renderUI({
        print("render dropZoneInput")
        drop_input <- dropZoneInput(ns("PPI_dropzone"),
                                    choices = drop_choose(),
                                    # list(one = exctract_genelist_UI("one"),
                                    #              two = exctract_genelist_UI("two")),
                                    # presets = list(values = "one-ds-1"),
                                    selectable = F,
                                    multivalued = F,
                                    server = exctract_genelist_UI_PPI,
                                    style = "height: 450px;
                                    overflow-y: auto;
                                    color: #6495ED;
                                    font-color: black;
                                    border-radius: 6px;"
        )
        dropZoneServer(session, "PPI_dropzone", exctract_genelist_UI_PPI)
        return(tagList(h5("Dropzone"), drop_input))
      })
      # print("render PPI_dropzone finished")

      drop_genelist_modules <- reactiveValues()
      observeEvent(input$PPI_dropzone, {
        check = reactiveValuesToList(drop_genelist_modules) %>% lapply(., function(x){
          return( class(try(x(), silent = T)) != "try-error" && !is.null(x()) )
        }) %>% unlist
        purrr::map(
          setdiff(input$PPI_dropzone, names(drop_genelist_modules)[check]),
          ~ { drop_genelist_modules[[.]] <- exctract_genelist_Server(id = ., Omics_Serv_res = Omics_res_list, suffix = "PPIgenelist") } ## defined in genelist_tool_modules_new
        )
        # cat(paste0("outside ns is ",session$ns(),"\n"))
        purrr::map(
          setdiff(names(drop_genelist_modules)[check], input$PPI_dropzone),
          ~ { drop_genelist_modules[[.]] <- NULL })
        print("exctract_genelist_Server load")
        # drop_genelist_modules_save <<- drop_genelist_modules
      }, priority = 1, ignoreNULL = F)


      #* handle gene_name input ----
      gene_name <- reactive({
        if(input$import_from_text_paste){ ## text input
          temp <- input$text_input_for_PPI
          # check the nrow of text input
          check <- strsplit(temp,'\n')[[1]]
          if(is.null(check) || length(check) == 0)
            return(NULL)

          genelist <- fread(text = temp, header = F) %>% as.data.frame()
          # genelist_save <<- genelist

          # if(ncol(genelist) < 1) return(NULL)
          gene_name = genelist[,1] %>% sapply(., the1stname) %>% unique()
          # gene_name_save1 <<- gene_name
          # colnames(genelist)[1] = c("name")
          # genelist$name <-  sapply(genelist$name,the1stname)
          # genelist$fc = as.numeric(genelist$fc)

        }else{ ## drag-drop input
          genelist <- drop_genelist_modules %>% reactiveValuesToList()
          # genelist_save <<- genelist
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
        # gene_name_save1 <<- gene_name
        if(is.null(gene_name)||length(gene_name) < 1) return(NULL) ## when no sig res (no gene as input)
        gene_name = data.frame(name = gene_name %>% DEP2:::rm_digit_end()  %>% unique()) %>%
          dplyr::filter(!is.na(name))
        # gene_name = name = gene_name %>% DEP2:::rm_digit_end()  %>% unique()
        # gene_name = gene_name[!is.na(gene_name)]
        # gene_name_save2 <<- gene_name
        return(gene_name)
      })

      speciesname <- input$organism

      taxonomy_df <- data.frame(
        species = c(
          "Anopheles", "Arabidopsis", "Bovine", "Worm",
          "Canine", "Fly", "Zebrafish", "E coli strain K12",
          "Chicken", "Human", "Mouse",
          "Rhesus", "Malaria", "Chimp", "Rat",
          "Yeast", "Streptomyces coelicolor", "Pig", "Toxoplasma gondii",
          "Xenopus"
        ),
        taxonomy = c(
          "7165", "3702", "9913", "6239", "9612", "7227", "7955", "511145", "9031", "9606", "10090", "9544", "5833", "9598", "10116", "4932", "100226", "9823", "5811", "8364"
        ),
        pkg = c(
          "org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db",
          "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db",
          "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db",
          "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db", "org.Rn.eg.db",
          "org.Sc.sgd.db", "org.Sco.eg.db", "org.Ss.eg.db", "org.Tgondii.eg.db",
          "org.Xl.eg.db"
        ),
        stringsAsFactors = FALSE
      )

      script_dir <- system.file("PPIdata",package = "DEP2")
      version_download <- "11.0"
      species_ID <- taxonomy_df$taxonomy[taxonomy_df$species == speciesname]

      proteinsDf_url <- paste("https://stringdb-static.org/download/protein.info.v", version_download, "/", species_ID, ".protein.info.v", version_download, ".txt.gz", sep = "")
      aliases_url <- paste("https://stringdb-static.org/download/protein.aliases.v", version_download, "/", species_ID, ".protein.aliases.v", version_download, ".txt.gz", sep = "")
      links_detail_url <- paste("https://stringdb-static.org/download/protein.links.detailed.v", version_download, "/", species_ID, ".protein.links.detailed.v", version_download, ".txt.gz", sep = "")
      proteinsDf_file <- paste(script_dir,"/",speciesname, "/", species_ID, ".protein.info.v", version_download, ".txt.gz", sep="")
      aliases_file = paste(script_dir,"/",speciesname, "/", species_ID, ".protein.aliases.v", version_download, ".txt.gz", sep="")
      links_detail_file = paste(script_dir,"/",speciesname, "/", species_ID, ".protein.links.detailed.v", version_download, ".txt.gz", sep="")

      strEnv = get_string_Env()


      #* observeRvent for PPITools ----
      observeEvent(input$String_annalysis,{
        withProgress(message = 'Please wait', value = 0, {
          memory.limit(2000000)#2GB


          #** download ----
          output$downloadPPITable <- downloadHandler(
            filename = "stringResult.txt",
            content = function(file) {
              write.table(links2(),
                          file,
                          col.names = TRUE,
                          row.names = FALSE,
                          sep ="\t") }
          )

          output$downloadPPInetwork <- downloadHandler(
            filename = "stringNetwork.html",
            content = function(file) {
              aplot <- aplot()
              htmlwidgets::saveWidget(aplot, file)}
          )

          # validate(need(nrow(gene_name()) != 0, "No genes meet your requirements, and can not do the PPI(Protein-Protein Interaction Network) analysis"))

          speciesname = input$organism
          # print(speciesname)
          orgDbname = reactive({
            the_annoSpecies_df = DEP2:::annoSpecies_df()
            orgDB = the_annoSpecies_df$pkg[the_annoSpecies_df$species == species]
            # require(orgDB, character.only = TRUE)
            # switch(speciesname,
            #        "Anopheles" = "org.Ag.eg.db",
            #        "Arabidopsis" = "org.At.tair.db",
            #        "Bovine" = "org.Bt.eg.db",
            #        "Worm" = "org.Ce.eg.db",
            #        "Canine" = "org.Cf.eg.db",
            #        "Fly" = "org.Dm.eg.db",
            #        "Zebrafish" = "org.Dr.eg.db",
            #        "E coli strain K12" = "org.EcK12.eg.db",
            #        "Chicken" = "org.Gg.eg.db",
            #        "Human" = "org.Hs.eg.db",
            #        "Mouse" = "org.Mm.eg.db",
            #        "Rhesus" = "org.Mmu.eg.db",
            #        "Malaria" = "org.Pf.plasmo.db",
            #        "Chimp" = "org.Pt.eg.db",
            #        "Rat" = "org.Rn.eg.db",
            #        "Yeast" = "org.Sc.sgd.db",
            #        "Streptomyces coelicolor" = "org.Sco.eg.db",
            #        "Pig" = "org.Ss.eg.db",
            #        "Toxoplasma gondii" = "org.Tgondii.eg.db",
            #        "Xenopus" = "org.Xl.eg.db")
          })


          # shiny::validate(
          #   need(
          #     require(orgDbname, character.only = TRUE),
          #     paste0("The package ", orgDbname, " is not installed/available. Try installing it with BiocManager::install('", orgDbname, "')")
          #   )
          # )


          check_ppi <- reactive({
            print("check ppi data")
            # spe_save <<- input$organism
            try(DEP2::load_PPIdata(input$organism), silent = TRUE)
          })

          links <- reactive({
            print("Searching links.")
            links <<- DEP2::test_PPI(gene_name()$name, species = input$organism,
                                          choose_scores = input$chooseScore, score_cutoff = 400)
            links_save <<- links
          })
          links2 <- reactive({
            print("Filter links.")
            links2 <<- dplyr::filter(links(), combined_score >= input$scorecutoff)
            links2_save <<- links2
          })



          # geneID <- reactive({
          #   # isol
          #   require(orgDbname, character.only = TRUE)
          #   geneID <- my_to_entrezid(orgDB = get(orgDbname), gene = as.character(gene_name()$name)) %>%
          #     tibble::rownames_to_column() %>% dplyr::rename(SYMBOL = rowname, ENTREZID = id) %>%
          #     dplyr::select(SYMBOL, ENTREZID) %>% dplyr::filter(!is.na(ENTREZID))
          #
          #   # geneID_save <<- geneID
          # })
          #
          # incProgress(0.3)
          #
          # data_mapped <- reactive({
          #   check_ppi()
          #   aliasDf <- get("aliasDf",envir = strEnv)
          #   proteinsDf <- get("proteinsDf",envir = strEnv)
          #   protein_links_detail <- get("protein_links_detail",envir = strEnv)
          #   data_mapped <- geneID() %>% mymap(my_data_frame_id_col_names = "ENTREZID",
          #                                     aliasDf=get("aliasDf",envir = strEnv),
          #                                     proteinsDf=get("proteinsDf",envir = strEnv), ## data_mapped 3cols ENTREZID SYMBOL STRING_id
          #                                     removeUnmappedRows = TRUE,reload = F)
          # })
          # link_table <- reactive({
          #   data_mapped <- data_mapped()
          #   myget_interactions(protein_detail = get("protein_links_detail",envir = strEnv),mapped_data = data_mapped,choose_scores = input$chooseScore, score_cutoff = 400)
          # })
          # links <- reactive({
          #   data_mapped <- data_mapped()
          #   link_table() %>%                                                ## links 3cols from to weight from&to is SYMBOL
          #     mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>%
          #     mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"])
          # })
          # links2 <- reactive({filterlink(links(), scorecutoff=input$scorecutoff) })


          output$String_Table <- DT::renderDataTable({
            # print("b1")
            # if(nrow(geneID()) == 0){
            #   # return(gene_name())
            #   require(orgDbname, character.only = TRUE)
            #   return(my_to_entrezid(orgDB = get(orgDbname), gene = as.character(gene_name()$name)) %>%
            #            tibble::rownames_to_column() %>% dplyr::rename(SYMBOL = rowname, ENTREZID = id) %>%
            #            dplyr::select(SYMBOL, ENTREZID))
            # }
            validate(need(nrow(gene_name()) > 1, "No genes (or less than 2) meet your requirements, and can not do the PPI(Protein-Protein Interaction Network) analysis!"))

            # require(orgDbname, character.only = TRUE)
            # orgDB = get(orgDbname)
            # geneID = map_to_entrezid(gene_name(), orgDB = orgDB)
            # validate(need(nrow(geneID) > 1, "No genes (or less than 2) meet your requirements, and can not do the PPI(Protein-Protein Interaction Network) analysis!"))

            validate(
              need(class(check_ppi()) != "try-error",
                   paste("The local STRING data of select organism is absent! You need to download file of organism by runing DEP2::load_PPIdata")
              )
            )
            print("b4")
            validate(need(nrow(links2()) != 0, "There is no Protein-Protein Interaction Network"))
            print("b5")
            validate(
              need(
                require(orgDbname, character.only = TRUE),
                paste0("The package ", orgDbname, " is not installed/available. Try installing it with BiocManager::install('", orgDbname, "')")
              )
            )

            links2()
          }, filter = 'top', options = list( autoWidth = F,scrollX = TRUE) )

          #** output ----
          output$download_PPITable <- renderUI({
            downloadButton(ns("downloadPPITable"), "Save table", class = "downloadPPITable")
          })

          output$download_PPInetwork <- renderUI({
            downloadButton(ns("downloadPPInetwork"), "Save network", class = "downloadPPInetwork")
          })

          output$PPI_plot <- renderUI({
            print("a1")
            validate(need(nrow(gene_name()) != 0, "No genes meet your requirements, and can not do the PPI(Protein-Protein Interaction Network) analysis!"))
            print("a2")
            # validate(need(nrow(geneID()) != 0, "No genes meet your requirements, and can not do the PPI(Protein-Protein Interaction Network) analysis"))
            print("a3")
            validate(
              need(class(check_ppi()) != "try-error",
                   paste("You need to download file of organism '" , speciesname, "' from string first.\nPlease try to download them in R with:\n",
                         "download.file(", "'", proteinsDf_url, "' , ", "'", proteinsDf_file, "' , mode = 'wb')\n",
                         "download.file(", "'", aliases_url, "' , ", "'", aliases_file, "' , mode = 'wb')\n",
                         "download.file(", "'", links_detail_url, "' , ", "'", links_detail_file, "' , mode = 'wb')\n", sep = "")
              )
            )
            print("a4")
            validate(need(nrow(links2()) != 0, "There is no Protein-Protein Interaction Network"))
            print("a5")
            # validate(
            #   need(
            #     require(orgDbname, character.only = TRUE),
            #     paste0("The package ", orgDbname, " is not installed/available. Try installing it with BiocManager::install('", orgDbname, "')")
            #   )
            # )
            print("a6")
            visNetworkOutput(ns("network"),height=input$plothigh*10)
          })
          output$network <- renderVisNetwork({
            # validate(need(nrow(links2()) != 0, "There is no Protein-Protein Interaction Network"))
            aplot()
          })

          aplot <- reactive({
            # isolate(linksTable = links2())
            shiny::validate(need(nrow(gene_name()) > 1, "PPI should perform on at least 2 proteins/genes"))
            aplot <- stringNetwork(linksTable = links2() , layoutway=input$layoutway , nodecolor= input$nodecolor , nodeshape=input$nodeshape, linecolor=input$linecolor , nodesize=input$nodesize , changesize=input$changesize , fontsize=input$PPIfontsize , changewidth=input$changewidth , linewidth=input$linewidth,smoothline = input$smoothline, smoothtype = input$smoothtype,highlightkey=input$highlightkey)
            # aplot_save <<- aplot_save
          })

          incProgress(0.5)

          incProgress(0.2)
        })

      })
    }
  )
}




