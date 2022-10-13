#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


# source("app_can_use_2.1.20_source.R")
source("source.R")
# source("genelist_tool_modules.R")
source("ORA_module.R")
source("GSEA_module.R")
source("DEP_pg_modules.R")
source("DEG_modules.R")
source("genelist_tool_modules_new.R")
source("DEP_pep_modules.R")
source("DEP_ptm_modules.R")
source("PPI_module.R")
source("Timecourse_module.R")

library(BiocStyle)
library(dplyr)
library(tibble)
library(purrr)
library(SummarizedExperiment)
library(DEP2)
require(shiny)
# require(shinydashboard)
require(shinyBS)
# library(shinydashboardPlus)
library(shinyvalidate)
library(shinyWidgets)
library(shinyjs)
# library(bslib)


rm(list = ls(envir = get_DEPglobal_env(), all.names = T), envir = get_DEPglobal_env())
# initialization
# initialized_omicsmod = data.frame(name = c("Proteomic","PTM","RNA"),
#                                   type = c("Pg","PTM","RNAseq"),
#                                   ID = c("Pg_1","PTM_1","RNAseq_1"))
initialized_omicsmod <- data.frame(
  name = c(""),
  type = c(""),
  ID = c("")
)

# initialized_postmod = data.frame(name = c("Annotation","ORA","GSEA","PPI tool","Instegrated analysis"),
#                                  type = c("Annotation","ORA","GSEA","PPI","Instegrated"),
#                                  ID = c("Annotation_1","ORA_1","GSEA_1","PPI_1","Instegrated_1"))
initialized_postmod <- data.frame(
  name = c(""),
  type = c(""),
  ID = c("")
)
assign(".OmicsModule", initialized_omicsmod, envir = get_DEPglobal_env())
assign(".PostModule", initialized_postmod, envir = get_DEPglobal_env())
# rm(list = ls(envir = get_DEPglobal_env(),all.names = T), envir = get_DEPglobal_env())

# Define UI
ui <- navbarPage(
  title = "DEP2",
  theme = shinythemes::shinytheme("spacelab"),
  id = "DEPnavbar",
  useShinyjs() ,
  tabPanel(
    title = "Welcome",
    header = h1("asdf"),
    fluidRow(
      column(
        style = "background-color:white;",
        width = 12,
        # includeMarkdown("Tutorial.md")
        # uiOutput("htmlpage")
        # includeHTML("../rmd_test/test2.html")
        # includeCSS("../rmd_test/test2.html") ## use includeCSS to insert html vignette with css. includeHTML cause
        # includeMarkdown("../rmd_test/test2.md")
      )
    )
  ),
  ## add mod button
  tags$script(
    HTML("var header = $('.navbar > .container-fluid');
                              header.append('<div style=\"float:right; padding-top: 8px\"><button id=\"addmod\" type=\"button\" class=\"btn btn-primary action-button\" onclick=\"AddModule()\">Add module</button></div>')")
  ),
  tags$head(tags$style(HTML("
                               body {
                                  width: 100% !important;
                                  max-width: 100% !important;
                               }
                               ")))
)





# Define global server logic for navbarPage
server <- function(input, output,session = session) {
  options(shiny.maxRequestSize = 5000 * 1024^2)
  ## store the omics analysis servers result, for the crosstalk of different modules
  omics_Servers <- reactiveValues()
  # bs_themer()
  # DEP_env = reactive({get_DEPglobal_env()})
  ## input of addmod
  output$htmlpage <- renderUI({

    # HTML(markdown::markdownToHTML(
    #   knitr::knit("../rmd_test/test2.Rmd", quiet = TRUE),
    #   # "../rmd_test/test2.html",
    #   fragment.only=T))
    tags$iframe(
      seamless="seamless",
      src="../rmd_test/test2.html",
      width="100%", height="400", frameBorder="0")

  })

  observeEvent(input$addmod, {
    # id <- paste0("Dropdown", input$append, "a")
    showModal(modalDialog(
      footer = NULL,
      size = "m",
      h2("Add analysis module"),
      tabsetPanel(
        tabPanel(
          "Omics analysis pipeline",
          tags$h4("Add omics analysis module"),
          fluidRow(
            column(
              width = 12,
              awesomeRadio(
                inputId = "OmicsType",
                label = "Select omics type",
                width = "100%",
                choices = list(
                  "Quantificated proteome (based on proteingroup quantification)" = "Proteomepg",
                  "Quantificated proteome (based on peptide quantification)" = "Proteomepep",
                  "Modification-specific proteome (modified peptide quantification)" = "PTM",
                  "Transcriptome (RNA counts)" = "RNAseq"
                )
              )
            ),
          ),
          shinyFeedback::useShinyFeedback(),
          textInput("OmicsName", "Input the name of omics"),
          tagList(
            actionButton("add_omic", "Add"),
            modalButton("Cancel")
          )
        ),
        tabPanel(
          "Post analysis",
          tags$h4("Add post-analysis module"),
          fluidRow(
            column(
              width = 12,
              awesomeRadio(
                inputId = "PostanalysisType",
                label = "Select post-analysis type",
                width = "100%",
                choices = c(
                  "Over-representation analysis" = "ORA",
                  "GSEA" = "GSEA",
                  "Protein-protein interaction" = "PPI",
                  "Timecourse cluster" = "Timecourse",
                  "Integrated analysis" = "Instegrated"
                )
              )
            )
          ),
          tagList(
            actionButton("add_postanalysis", "Add"),
            modalButton("cancel")
          )
        )
      ),
    ))
  })

  ## add omic modules
  observeEvent(input$add_omic,
    {
      ID <- NULL
      theOmicsName_Pass1 <- ifelse((is.null(input$OmicsName) | input$OmicsName == ""), F, T)
      shinyFeedback::feedbackWarning("OmicsName", !theOmicsName_Pass1, "Please input an name")
      req(theOmicsName_Pass1)

      if (exists(".OmicsModule", envir = get_DEPglobal_env())) {
        OmicsModules <- get(".OmicsModule", envir = get_DEPglobal_env())
        theOmicsName_Pass2 <- ifelse(input$OmicsName %in% OmicsModules$name, F, T)
        shinyFeedback::feedbackWarning("OmicsName", !theOmicsName_Pass2, "This name was used, please change")
        req(theOmicsName_Pass2)

        if (length(which(OmicsModules$type == input$OmicsType)) > 0) {
          ID <- OmicsModules[which(OmicsModules$type == input$OmicsType), "ID"] %>%
            .[length(.)] %>%
            sapply(., function(x) {
              strsplit(x, "_")[[1]] %>% .[length(.)]
            }) %>%
            as.numeric()
        } else {
          ID <- 0
        }
        ID <- paste(input$OmicsType, ID + 1, sep = "_")
        OmicsModules <- rbind(OmicsModules, c(name = input$OmicsName, type = input$OmicsType, ID = ID))
        assign(".OmicsModule", OmicsModules, envir = get_DEPglobal_env())
      } else {
        ID <- paste(input$OmicsType, 1, sep = "_")
        assign(".OmicsModule", data.frame(name = input$OmicsName, type = input$OmicsType, ID = ID), envir = get_DEPglobal_env())
      }

      ns <- reactive({
        NS(ID)
      })

      #*import omics UI modules and Server modules ----
      if (input$OmicsType == "Proteomepg") {
        ## import DEPpg analyze module
        appendTab(
          inputId = "DEPnavbar",
          tabPanel(
            ID,
            DEP_pg_UI(id = ID)
          ),
          # target = "Welcome",
          select = T
        )
        cat(ID)
        cat(" load UI \n")
        print(paste0("proteomics ns is ", ns()("")))
        cat(ID)
        cat(" load server \n")
        # omics_Servers[[ID]] <- DEP_server_module(ID)

        observeEvent(input[[ns()("threshold_method")]],
          {
            if (!is.null(ID)) {
              cat(ID)
              cat(" load server \n")
              # DEP_server_module(ID)
              omics_Servers[[ID]] <- DEP_pg_server_module(ID) ## DEP Server. The test values are stored in omics_Servers for postanalysis modules to call
              omics_Servers_saved <<- omics_Servers
              # omics_Servers_list <- reactive({reactiveValuesToList(omics_Servers)})
            }
          },
          priority = 0,
          ignoreNULL = T,
          ignoreInit = F,
          once = T
        )
      } else if (input$OmicsType == "RNAseq") {
        ## import DEG test module
        appendTab(
          inputId = "DEPnavbar",
          tabPanel(
            ID,
            sidebarLayout(
              DEG_sidebar_mod(id = ID, labelname = input$OmicsName), ## DEG sidebar module
              # mainPanel()
              DEG_body_mod(id = ID, labelname = input$OmicsName) ## DEG pagebody module
            ),
          ),
          # target = "Welcome",
          select = T
        )
        cat(ID)
        cat(" load UI \n")
        observeEvent(input[[ns()("filter_rowsum")]],
          {
            if (!is.null(ID)) {
              cat(ID)
              cat(" load server \n")
              omics_Servers[[ID]] <- DEG_server_module(ID) ## DEG Server. The test values are stored in omics_Servers for postanalysis modules to call
              omics_Servers_saved <<- omics_Servers
            }
          },
          priority = 0,
          ignoreNULL = T,
          ignoreInit = F,
          once = T
        )
      } else if (input$OmicsType == "Proteomepep") {
        ## import DEPpep analyze module
        appendTab(
          inputId = "DEPnavbar",
          tabPanel(
            ID,
            DEP_pep_UI(ID),
          ),
          # target = "Welcome",
          select = T
        )
        cat(ID)
        cat(" load UI \n")
        observeEvent(input[[ns()("threshold_method")]],
          {
            if (!is.null(ID)) {
              cat(ID)
              cat(" load server \n")
              omics_Servers[[ID]] <- DEP_pep_server_module(ID) ## DEG Server. The test values are stored in omics_Servers for postanalysis modules to call
              omics_Servers_saved <<- omics_Servers
            }
          },
          priority = 0,
          ignoreNULL = T,
          ignoreInit = F,
          once = T
        )
      }else if (input$OmicsType == "PTM") { ## import DEP_ptm test module
        appendTab(inputId = "DEPnavbar",
                  tabPanel(ID,
                           DEP_ptm_UI(id = ID)
                  ),
                  # target = "Welcome",
                  select = T
        )
        cat(ID)
        cat(" load UI \n")
        print(paste0("PTM ns is ",ns()("")))
        cat(ID)
        cat(" load server \n")
        # omics_Servers[[ID]] <- DEP_server_module(ID)

        observeEvent(input[[ns()("threshold_method_ptm")]],{
          if(!is.null(ID)){
            cat(ID)
            cat(" load server \n")
            # DEP_server_module(ID)
            omics_Servers[[ID]] <- DEP_ptm_server_module(ID, Omics_res = omics_Servers)   ## DEP_ptm Server. The test values are stored in omics_Servers for postanalysis modules to call
            omics_Servers_saved <<- omics_Servers
            # omics_Servers_list <- reactive({reactiveValuesToList(omics_Servers)})
          }
        },priority = 0, ignoreNULL =T, ignoreInit = F,once = T)
      }


      removeModal() ## close modalDialog
    },
    priority = 10,
    ignoreInit = F
  )


  #* import post analysis modules -----
  observeEvent(input$add_postanalysis, {
    if (exists(".OmicsModule", envir = get_DEPglobal_env())) {
      Postmodules <- get(".PostModule", envir = get_DEPglobal_env())
      if (length(which(Postmodules$type == input$PostanalysisType)) > 0) {
        ID <- Postmodules[which(Postmodules$type == input$PostanalysisType), "ID"] %>%
          .[length(.)] %>%
          sapply(., function(x) {
            strsplit(x, "_")[[1]][2]
          }) %>%
          as.numeric()
      } else {
        ID <- 0
      }
      ID <- paste(input$PostanalysisType, ID + 1, sep = "_")
      PostModule <- rbind(Postmodules, c(name = ID, type = input$PostanalysisType, ID = ID))
      assign(".PostModule", PostModule, envir = get_DEPglobal_env())
    } else {
      ID <- paste(input$PostanalysisType, 1, sep = "_")
      assign(".PostModule", data.frame(name = ID, type = input$PostanalysisType, ID = ID), envir = get_DEPglobal_env())
    }

    # add instegrated analysis. Or genelist tool
    if (input$PostanalysisType == "Instegrated") {
      appendTab(
        inputId = "DEPnavbar",
        tabPanel(ID, genelist_tool_UI(ID)), ## inserttab of genelis tool UI
        # target = "Home",
        select = T
      )
      #  call genelist tool modules tool, transmit omics module infor & test result
      # omics_Servers_list <- reactiveValuesToList(omics_Servers)
      # omics_Servers_list_saved <<- omics_Servers_list
      # cat(names(omics_Servers_list))
      # ID_saved <<- ID
      genelist_tool_Server(ID, Omics_res = omics_Servers)
    }else if(input$PostanalysisType == "ORA") {
      appendTab(inputId = "DEPnavbar",
                tabPanel(ID,
                         ORA_UI(id = ID)
                ),
                # target = "Welcome",
                select = T
      )
      # session_save3 <<- session
      ORA_server_module2(ID, Omics_res = omics_Servers)
    } else if(input$PostanalysisType == "GSEA") {
      appendTab(inputId = "DEPnavbar",
                tabPanel(ID,
                         GSEA_UI(id = ID)
                ),
                # target = "Welcome",
                select = T
      )
      GSEA_server_module(ID, Omics_res = omics_Servers)
    }else if(input$PostanalysisType == "PPI"){
      appendTab(inputId = "DEPnavbar",
                tabPanel(ID,
                         PPI_UI(id = ID)
                ),
                # target = "Welcome",
                select = T
      )
      PPI_server_module(ID, Omics_res = omics_Servers)
    }else if(input$PostanalysisType == "Timecourse"){
      appendTab(inputId = "DEPnavbar",
                tabPanel(ID,
                         Timecourse_UI(id = ID)
                ),
                # target = "Welcome",
                select = T
      )
      ## the time-course cluster result is also saved in omics_Servers for following extract.
      omics_Servers[[ID]] <- Timecourse_server_module(ID, Omics_res = omics_Servers)
    }

    removeModal() ## close modalDialog
  })

}

# Run the application
shinyApp(ui = ui, server = server)
