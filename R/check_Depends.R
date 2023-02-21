## functions for package install
checek_req_depends <- function(need_packages,need_version =NULL, install = T, mess = "", pak_install = T, bioc_pkg = NULL){
  need_packages2 = gsub(".*/","",need_packages)

  ## check whether need_packages installed
  lacks = if(is.null(need_version)){
    !sapply(need_packages2, rlang::is_installed,version = need_version)
  }else{
    suppressWarnings(!mapply(rlang::is_installed, pkg = need_packages2, version = need_version, SIMPLIFY = T))
  }

  lack_packages = need_packages[lacks]
  need_version = need_version[lacks]

  lack_packages_info = paste0(lack_packages,ifelse(is.na(need_version),"",paste0("(>",need_version,")")))
  if(length(lack_packages)>0){
    cat("Lack such required packages ", paste0(lack_packages_info,collapse = ", "), " ", mess,".\n",
        sep = "")
    if(install){ # install lacked packages
      # rlang::check_installed(lack_packages,version = need_version)
      check_installed(lack_packages,version = need_version, pak_install = pak_install,bioc_pkg = bioc_pkg[bioc_pkg %in%lack_packages])
    }else{
      return(lack_packages)
    }

    lack_packages = ifelse(is.null(need_version),
                           need_packages[!sapply(need_packages2, rlang::is_installed,version = need_version)],
                           need_packages[!mapply(rlang::is_installed, pkg = need_packages2, version = need_version, SIMPLIFY = T)]
    )
    if(length(lack_packages) > 0 && !is.na(lack_packages)){
      return(lack_packages)
    }else{
      return(T)
    }
  }else{
    return(T)
  }
}

## functions to check required packages for DESeq2
#' Check required packages for DESeq2.
#'
#' Check the depend packages for RNA-Seq analysis workflow based on DESeq2
#'
#' @param install logic(1). Whether install the missing packages.
#'
#' @export
#'
check_RNAseq_depends <- function(install = T){
  need_packages = c(
    "DESeq2", "IHW", "ashr", "apeglm"
  )
  bioc_pkg = c( "DESeq2", "IHW", "apeglm")
  checek_req_depends(need_packages,install=install, mess = "for RNAseq workflow",bioc_pkg = bioc_pkg)
}

## functions to check required packages for shiny application
check_shiny_depends = function(install = T){
 need_packages = c("shinyvalidate",
                   "shinyWidgets",
                   "shinythemes",
                   "shinyBS",
                   "shinyFeedback",
                   "shinydashboard",
                   # shinydashboardPlus,
                   "markdown",
                   "serenity-r/dndselectr",
                   "rmarkdown",
                   "shiny",
                   "colourpicker"
                   )
 bioc_pkg = c("serenity-r/dndselectr")
 checek_req_depends(need_packages,install = install, mess = "for the shiny application", bioc_pkg = bioc_pkg)

}

## functions for ORA and GSEA anaysis. ClusterProfiler
#' Check required packages for enrichment analysis.
#'
#' Check the depend packages for functional enrichment. Including clusterProfiler, ReactomePA and msigdbr
#'
#' @param install logic(1). Whether install the missing packages.
#'
#' @export
#'
check_enrichment_depends = function(install = T){
  need_packages = c(
    "clusterProfiler", "YuLab-SMU/clusterProfiler.dplyr", "ReactomePA", "msigdbr",
    "org.Hs.eg.db", "org.Mm.eg.db","dbplyr"
  )
  bioc_pkg = need_packages
  need_version = c("4.3.1.900",rep(NA, times = 6))
  checek_req_depends(need_packages,need_version = need_version, install = install, mess = "for enrichment analysis", bioc_pkg = bioc_pkg)
}

#' Check organism annotation package
#'
#' Check and install organism annotation package(OrgDb from Bioconductor).
#'
#' @param organism Character(1), the organism name, could be one of c(
#' "Chicken", "Human", "Mouse", "Fly", "Rat", "Yeast"). Other species with
#' \href{http://bioconductor.org/packages/release/BiocViews.html#___OrgDb}{OrgDb} are also optional.
#' @param install logic(1). Whether install if organism annotation db package is missed.
#'
#' @export
#'
check_organismDB_depends <- function(organism,install = T){
  organism_df = annoSpecies_df()
  need_packages = organism_df$pkg[which(organism_df$species == organism)]
  bioc_pkg = need_packages
  checek_req_depends(need_packages, install = install, bioc_pkg = bioc_pkg)
}

#' Check required packages for PPI analysis.
#'
#' Check the depend packages for PPI analysis.
#'
#' @param install logic(1). Whether install the missing packages.
#'
#' @export
#'
check_PPI_depends <- function(install = T){
  need_packages = c(
    # "dbplyr",
    "igraph", "visNetwork","R.utils"
  )
  checek_req_depends(need_packages, install = install, mess = "for PPI enrichment")
}

# pkg = "serenity-r/dndselectr"
# check_installed("serenity-r/dndselectr")
# reason = NULL; version = NULL; compare = NULL;
# action = NULL; call = rlang::caller_env(); pak_install = T
check_installed <- function (pkg, reason = NULL, ..., version = NULL, compare = NULL,
          action = NULL, call = rlang::caller_env(), pak_install = T,
          bioc_pkg = NULL)
{
  rlang::check_dots_empty0(...)
  rlang:::check_action(action)
  abort <- rlang::env_get(rlang::ns_env("rlang"), "abort")
  info <- rlang:::pkg_version_info(pkg, version = version, compare = compare)
  needs_install <- !rlang:::detect_installed(info)
  pkg <- info$pkg
  version <- info$ver
  compare <- info$cmp
  missing_pkgs <- pkg[needs_install]
  missing_vers <- version[needs_install]
  missing_cmps <- compare[needs_install]
  if (!length(missing_pkgs)) {
    return(invisible(NULL))
  }
  cnd <- rlang:::new_error_package_not_found(missing_pkgs, missing_vers,
                                     missing_cmps, reason = reason, call = call)
  restart <- rlang::peek_option("rlib_restart_package_not_found") %||%
    TRUE
  if (!rlang::is_bool(restart)) {
    abort("`rlib_restart_package_not_found` must be a logical value.")
  }
  if (!rlang::is_interactive() || !restart || any(missing_cmps %in%
                                           c("<", "<="))) {
    stop(cnd)
  }
  if (rlang:::signal_package_not_found(cnd)) {
    return(with_options(rlib_restart_package_not_found = FALSE,
                        check_installed(pkg, reason, version = version, call = call)))
  }
  header <- rlang::cnd_header(cnd)
  n <- length(missing_pkgs)
  question <- rlang:::pluralise(n, "Would you like to install it?",
                        "Would you like to install them?")
  question <- rlang:::paste_line(paste0(rlang:::ansi_info(), " ", header),
                         paste0(rlang:::ansi_cross(), " ", question), .trailing = TRUE)
  if (rlang::is_true(rlang::peek_option("rlang:::check_installed_test_hook"))) {
    return(question)
  }
  cat(question)
  if (utils::menu(c("Yes", "No")) != 1) {
    invokeRestart("abort", cnd)
  }
  if (!rlang::is_null(action)) {
    action(missing_pkgs)
  }else if (rlang::is_installed("pak") && pak_install) {
    missing_pkgs2 = missing_pkgs[!missing_pkgs %in% bioc_pkg]
    missing_bioc_pkgs = missing_pkgs[missing_pkgs %in% bioc_pkg]
    if(length(missing_pkgs2) > 0){
      pkg_install <- rlang::env_get(rlang::ns_env("pak"), "pkg_install")
      pkg_install(missing_pkgs2, ask = FALSE)
    }
    if(length(missing_bioc_pkgs) > 0){
      BiocManager::install(missing_bioc_pkgs,update = FALSE)
    }
  }else {
    utils::install.packages(missing_pkgs)
  }
}

