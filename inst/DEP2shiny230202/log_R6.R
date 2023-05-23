module_log = R6::R6Class(
  classname = "module_log",

  private = list(
    shiny_input = NULL,
    reactive_values = NULL
  ),

  public = list(
    initialize = function(input, reactive_vals, app_version = NULL){
      private$shiny_input = input
      private$reactive_values = reactive_vals
      self$appVersion = app_version
    },

    appVersion = NULL,

    inputVals = function(){
      private$shiny_input
    },

    resultVals = function(){
      private$reactive_values
    }
  )

)

pg_log = R6::R6Class(
  classname = "pg_log",
  inherit = module_log
)

ptm_log = R6::R6Class(
  classname = "ptm_log",
  inherit = module_log
)

pgpep_log = R6::R6Class(
  classname = "pgpep_log",
  inherit = module_log
)

transcriptome_log = R6::R6Class(
  classname = "transcriptome_log",
  inherit = module_log
)
