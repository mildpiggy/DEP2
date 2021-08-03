.onLoad <- function(libname = find.package("DEP2"), pkgname = "DEP2"){

  # CRAN Note avoidance
  if(getRversion() >= "2.15.1")
    utils::globalVariables(
      c( # functions.R globalVariables
        "ID", ".", "condition", "label",
        "value", "rowname", "miss_val",
        "samples", "logFC", "qval", "comparison",
        "val", "name", "significant", "contrasts",
        "CI.L", "CI.R", "P.Value", "variable", "temp",

        # iBAQ.R globalVariables
        "Unique..Groups.", "Protein.group.IDs",
        "mean_ctrl", "stoichiometry",

        # plot_functions.R globalVariables
        "bin", "Freq", "Var1", "missval", "num",
        "cs", "cs_frac", "cluster", "index",
        "iBAQ_value", "Peptides", "lfc", "ibaq",
        "LFC", "ymin", "ymax", "x", "y", "z",
        "qnorm", "error", "var",

        # plot_functions2.R globalVariables
        "n_con", "conditions", "category",

        # enrichR_functions.R globalVariables
        "Adjusted.P.value", "Term", "log_odds",
        "Overlap", "bg_IN", "bg_OUT", "IN", "OUT",
        "installed.packages", "contrast", "valid", "frac"
        )
    )
  invisible()
}
