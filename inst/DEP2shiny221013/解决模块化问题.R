library(shiny)
### create a tabset Panel in module environment
tabModuleUI <- function(id, title){
  ns = NS(id)
  return(
    tabPanel(title = title, value = id,
             shinydashboard::box(width = NULL,
                                 tabsetPanel(id = ns("Product_module"),
                                             tabPanel(title = "Product_module 0", value = "Product_module 0") #,
                                             #tabPanel(title = "Product 1", value = "Product 1")
                                 ),
                                 radioButtons(ns("threshold_method"),
                                              "Threshold method",
                                              choices = c("intersect", "curve"),
                                              selected = "intersect"),
                                 ## button to add new tab
                                 actionButton(inputId = ns("button"), label = "Click Module")))
  )
}

# tabModuleUI1 <- function(id){
#   ns = NS(id)
#   return(
#     sidebarPanel(radioButtons(ns("threshold_method"),
#                               "Threshold method",
#                               choices = c("intersect", "curve"),
#                               selected = "intersect"))
#   )
# }
# tabModuleUI2 <- function(id){
#   ns = NS(id)
#   return(
#     mainPanel(actionButton(inputId = ns("button"), label = "Click Module"))
#   )
# }

tabModuleServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      id_num = reactiveVal(1)
      parentSession <- get("session", envir = parent.frame(1))
      observe({
        # cat("dsf")
        if((threshold_method1 <<- input$threshold_method) == "intersect") {
          # cat("observe1")
          aasefwds <- reactive({input$threshold_method})
        }

      })
      observeEvent(input$button, {
        # cat("a")
        # cat("load ui")
        id_chr = paste0("Product_module", id_num())
        # server know which tab is shown
        # print(input[["Product_module"]])
        # print(input$button)
        ### but cannot add a newtab in it
        appendTab(inputId = ns("Product_module"), tab = tabPanel(title = id_chr, value = id_chr,
                                                                 shinydashboard::box()), session = parentSession)
        # ### even when i add a ns to inputId
        # appendTab(inputId = ns("Product_module"), tab = tabPanel(title = id_chr, value = id_chr,
        #                                                          shinydashboard::box()), session = session)
        ### but it can add a newtab in global environment
        # appendTab(inputId = "Product", tab = tabPanel(title = id_chr, value = id_chr,
        #                                               shinydashboard::box()), session = parentSession)
        id_num = id_num(id_num()+1)
      })
      
      
    }
  )}


## create a tabsetPanel in global environment

ui <- fluidPage(title="",
                tagList(tabsetPanel(
                  id = "Product",
                  tabPanel(title = "Product 0", value = "Product 0"))),
                ## button to add new tab
                actionButton(inputId = "button", label = "Click Main")
)

server <- function(input, output, session){
  id_numm = reactiveVal(1)
  id_numm2 = reactiveVal()
  # ns = reactiveVal(NS("test"))
  observeEvent(input[["button"]],{
    id_chr = paste0("Product", id_numm())
    # server know which tab is choosen
    # print(input[["Product"]])
    # add add a newtab to it
    appendTab(inputId = "Product", tab = tabModuleUI(id = id_chr, title = id_chr), select = T)
    # appendTab(inputId = "Product", tab = tabPanel(id_chr, 
    #                                               sidebarLayout(tabModuleUI1(id = id_chr),
    #                                                             # mainPanel()
    #                                                             tabModuleUI2(id = id_chr)
    #                                               )), select = T)
    
    # callModule(tabModuleServer, id = id_chr, idModule = id_chr)
    # tabModuleServer(id_chr)
    id_numm = id_numm(id_numm()+1)
    # id_numm2 = reactive({id_chr})
    
    
    
    id_numm2 = id_numm2({id_chr})
    # id_numm2_save <<- id_numm2()
    
    ns = reactive({NS(id_chr)})
    ns_save <<- ns()
    tm = reactive({input[[ns()("threshold_method")]]})
    
    # input2 <<- input
    # input3 <<- tm()
    cat(id_numm2())
    cat(" load ui \n")
    
    # observeEvent(input[[ns()("threshold_method")]],{
    #   cat(id_numm2())
    #   cat(" load mod \n")
    #   if(!is.null(id_numm2()))
    #     tabModuleServer(id_numm2())
    # },priority = 0, ignoreNULL =T, ignoreInit = T)
    
    invalidateLater(1000)
  },priority = 1, ignoreInit = F)
  
  
  
  thedata <- reactivePoll(1000,
                          session,
                          checkFunc = function(){
                            cat(id_numm())
                            cat("\n")
                            id_numm()
                            },
                          valueFunc = function(){
                            print("value \n")
                            if(!is.null(id_numm2())){
                            cat(id_numm2())
                            cat(" load mod \n")
                            # tabModuleServer(id_numm2())
                            thedata_save <<- "changed"
                            return("changed ")
                            }
                            "change"
                          }
                          )
  
  observe({
    invalidateLater(3000)
    cat("change2 \n")
    a <- reactive(thedata())
    # print(thedata())
  })
  
  
  # observeEvent(input[["button"]],{
  #   cat(id_numm2())
  #   cat(" load mod \n")
  #   if(!is.null(id_numm2()))
  #     tabModuleServer(id_numm2())
  # },priority = 0, ignoreNULL =F, ignoreInit = T)
  # observeEvent(id_numm2(),{
  #   cat(id_numm2())
  #   if(!is.null(ns())) temp <<- input[ns()("threshold_method")]
  #   
  #   if( (!is.null(id_numm2())) && (!is.null(ns())) && (!is.null(input[ns()("threshold_method")]))){
  #     cat(" called ")
  #     tabModuleServer(id_numm2())
  #   }
  #     
  # })
  # observe({
  #   cat(id_numm2())
  #   if(!is.null(id_numm2()))
  #     tabModuleServer(id_numm2())
  # })
  
}

shinyApp(ui,server)