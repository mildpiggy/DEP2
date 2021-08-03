context("9a - Wrapper functions")

test_that("import_MaxQuant throws error without valid input", {
  expect_error(import_MaxQuant("test_data", UbiLength_ExpDesign, filter = c("Reverse", "Potential.contaminant"),
                               intensities = "LFQ", names = "Gene.names", ids = "Protein.IDs", delim = ";"))
  expect_error(import_MaxQuant(test_data, "UbiLength_ExpDesign", filter = c("Reverse", "Potential.contaminant"),
                               intensities = "LFQ", names = "Gene.names", ids = "Protein.IDs", delim = ";"))
  expect_error(import_MaxQuant(test_data, UbiLength_ExpDesign, filter = test_data,
                               intensities = "LFQ", names = "Gene.names", ids = "Protein.IDs", delim = ";"))
  expect_error(import_MaxQuant(test_data, UbiLength_ExpDesign, filter = c("Reverse", "Potential.contaminant"),
                               intensities = LFQ, names = "Gene.names", ids = "Protein.IDs", delim = ";"))
  expect_error(import_MaxQuant(test_data, UbiLength_ExpDesign, filter = c("Reverse", "Potential.contaminant"),
                               intensities = "LFQ", names = Gene.names, ids = "Protein.IDs", delim = ";"))
  expect_error(import_MaxQuant(test_data, UbiLength_ExpDesign, filter = c("Reverse", "Potential.contaminant"),
                               intensities = "LFQ", names = "Gene.names", ids = Protein.IDs, delim = ";"))
  expect_error(import_MaxQuant(test_data, UbiLength_ExpDesign, filter = c("Reverse", "Potential.contaminant"),
                               intensities = "LFQ", names = "Gene.names", ids = "Protein.IDs", delim = test_data))
  expect_error(import_MaxQuant(test_data, UbiLength_ExpDesign, filter = "Bla",
                               intensities = "LFQ", names = "Gene.names", ids = "Protein.IDs", delim = ";"))
  expect_error(import_MaxQuant(test_data, UbiLength_ExpDesign, filter = c("Reverse", "Potential.contaminant"),
                               intensities = "Bla", names = "Gene.names", ids = "Protein.IDs", delim = ";"))
  expect_error(import_MaxQuant(test_data, UbiLength_ExpDesign, filter = c("Reverse", "Potential.contaminant"),
                               intensities = "Protein", names = "Gene.names", ids = "Protein.IDs", delim = ";"))
})

test_that("import_MaxQuant returns a SummarizedExperiment object", {
  expect_is(import_MaxQuant(test_data, UbiLength_ExpDesign), "SummarizedExperiment")
})

test_that("process throws error without valid input", {
  expect_error(process("test_se", thr = 0, fun = "MinProb"))
  expect_error(process(test_se, thr = "0", fun = "MinProb"))
  expect_error(process(test_se, thr = 0, fun = MinProb))
  expect_error(process(test_se, thr = 0, fun = "Bla"))
})

test_that("process returns a SummarizedExperiment object", {
  expect_is(process(test_se), "SummarizedExperiment")
})

test_that("analyze_dep throws error without valid input", {
  expect_error(analyze_dep("test_impute", "control", "Ctrl", alpha = 0.05,
                           lfc = 1, test = NULL, incl_repl = FALSE))
  expect_error(analyze_dep(test_impute, control, "Ctrl", alpha = 0.05,
                           lfc = 1, test = NULL, incl_repl = FALSE))
  expect_error(analyze_dep(test_impute, "control", Ctrl, alpha = 0.05,
                           lfc = 1, test = NULL, incl_repl = FALSE))
  expect_error(analyze_dep(test_impute, "control", "Ctrl", alpha = "0.05",
                           lfc = 1, test = NULL, incl_repl = FALSE))
  expect_error(analyze_dep(test_impute, "control", "Ctrl", alpha = 0.05,
                           lfc = "1", test = NULL, incl_repl = FALSE))
  expect_error(analyze_dep(test_impute, type = "manual", alpha = 0.05,
                           lfc = 1, test = bla, incl_repl = FALSE))
  expect_error(analyze_dep(test_impute, "control", "Ctrl", alpha = 0.05,
                           lfc = 1, test = NULL, incl_repl = "FALSE"))
  expect_error(analyze_dep(test_impute, type = "Bla", alpha = 0.05,
                           lfc = 1, test = NULL, incl_repl = FALSE))
})

test_that("analyze_dep returns a SummarizedExperiment object", {
  expect_is(analyze_dep(test_impute, "control", "Ctrl"), "SummarizedExperiment")
})

test_that("plot_all throws error without valid input", {
  expect_error(plot_all("test_sign", "volcano"))
  expect_error(plot_all(test_sign, volcano))
  expect_error(plot_all(test_sign, "Bla"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(1)]
  expect_error(plot_all(test_sign_error, "volcano"))

  test_sign_error2 <- test_sign
  SummarizedExperiment::rowData(test_sign_error2) <- SummarizedExperiment::rowData(test_sign_error2)[,-c(41:43)]
  expect_error(plot_all(test_sign_error2, "volcano"))

  test_sign_error3 <- test_sign
  SummarizedExperiment::rowData(test_sign_error3) <- SummarizedExperiment::rowData(test_sign_error3)[,-c(28:30,33:35,38:40)]
  expect_error(plot_all(test_sign_error3, "volcano"))
})

test_that("plot_all output ", {
  expect_message(plot_all(test_sign, "volcano"), c("Plot volcano:\n"))
  expect_message(plot_all(test_sign, "heatmap"), c("Plot heatmap:\n"))
  expect_message(plot_all(test_sign, "comparison"), c("Plot comparison graphs"))
  expect_message(plot_all(test_sign, "freq"), c("Plot frequency graphs"))
})

test_that("filter_MaxQuant throws error without valid input", {
  expect_error(filter_MaxQuant("test_data", c("Reverse", "Potential.contaminant")))
  expect_error(filter_MaxQuant(test_data, Potential.contaminant))
})

test_that("filter_MaxQuant returns a data.frame", {
  expect_is(filter_MaxQuant(test_data, c("Reverse", "Potential.contaminant")), "data.frame")
  expect_is(filter_MaxQuant(test_data, "Reverse"), "data.frame")
})

test_that("exclude_deps throws error without valid input", {
  expect_error(exclude_deps("test_sign", "Ubi6_vs_Ctrl"))
  expect_error(exclude_deps(test_sign, "Bla"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-c(41:43)]
  expect_error(exclude_deps(test_sign_error, "Ubi6_vs_Ctrl"))
})

test_that("exclude_deps returns a SummarizedExperiment object", {
  expect_is(exclude_deps(test_sign, "Ubi6_vs_Ctrl"), "SummarizedExperiment")
  expect_is(exclude_deps(test_sign, c("Ubi6_vs_Ctrl", "Ubi1_vs_Ctrl")), "SummarizedExperiment")
  expect_is(exclude_deps(test_sign, NULL), "SummarizedExperiment")
})

test_that("select_deps throws error without valid input", {
  expect_error(select_deps("test_sign", "Ubi6_vs_Ctrl"))
  expect_is(select_deps(test_sign, c("Ubi6_vs_Ctrl", "Ubi4_vs_Ctrl")), "SummarizedExperiment")
  expect_error(select_deps(test_sign, "Bla"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-c(41:43)]
  expect_error(select_deps(test_sign_error, "Ubi6_vs_Ctrl"))
})

test_that("select_deps returns a SummarizedExperiment object", {
  expect_is(select_deps(test_sign, "Ubi6_vs_Ctrl"), "SummarizedExperiment")
  expect_is(select_deps(test_sign, NULL), "SummarizedExperiment")
})

test_that("get_table throws error without valid input", {
  expect_error(get_table("test_results", "contrast"))
  expect_error(get_table(test_results, contrast))
  expect_error(get_table(test_results, "bla"))
  expect_error(get_table(test_results[,-c(10:12)], "contrast"))
  expect_error(get_table(test_results[,-c(13:16)], "contrast"))
})

test_that("get_table returns a data.frame", {
  expect_is(get_table(test_results, "contrast"), "data.frame")
  expect_is(get_table(test_results, "centered"), "data.frame")
})
