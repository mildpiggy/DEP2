context("5 - Plot functions results")

test_that("plot_single throws error without valid input", {
  expect_error(plot_single("test_sign", "USP19", "centered"))
  expect_error(plot_single(test_sign, USP19, "centered"))
  expect_error(plot_single(test_sign, "USP19", centered))
  expect_error(plot_single(test_sign, "USP19", "test"))
  expect_error(plot_single(test_sign, "USP15", "centered"))
  expect_error(plot_single(test_sign, "test", "centered"))

  test_sign_error1 <- test_sign
  SummarizedExperiment::colData(test_sign_error1) <- SummarizedExperiment::colData(test_sign_error1)[,-(3)]
  expect_error(plot_single(test_sign_error1, "USP19", "centered"))

  test_sign_error2 <- test_sign
  SummarizedExperiment::rowData(test_sign_error2) <- SummarizedExperiment::rowData(test_sign_error2)[,-(1)]
  expect_error(plot_single(test_sign_error2, "USP19", "centered"))

  test_sign_error3 <- test_sign
  SummarizedExperiment::rowData(test_sign_error3) <- SummarizedExperiment::rowData(test_sign_error3)[,-c(28:30,33:35,38:40)]
  expect_error(plot_single(test_sign_error3, "USP19", "centered"))
})

test_that("plot_single returns a ggplot object", {
  expect_is(plot_single(test_sign, "USP19", "centered"), "ggplot")
  expect_is(plot_single(test_sign, "USP19", "contrast"), "ggplot")
})

test_that("plot_heatmap throws error without valid input", {
  expect_error(plot_heatmap("test_sign", "centered", k = 6, col_limit = 4, labelsize = 8))
  expect_error(plot_heatmap(test_sign, centered, k = 6, col_limit = 4, labelsize = 8))
  expect_error(plot_heatmap(test_sign, "centered", k = "6", col_limit = 4, labelsize = 8))
  expect_error(plot_heatmap(test_sign, "centered", k = 6, col_limit = "4", labelsize = 8))
  expect_error(plot_heatmap(test_sign, "centered", k = 6, col_limit = 4, labelsize = "8"))
  expect_error(plot_heatmap(test_sign, "test", k = 6, col_limit = 4, labelsize = 8))
  expect_error(plot_heatmap(test_sign, "centered", k = 6, col_limit = 4, labelsize = 8, indicate = "bla"))
  expect_error(plot_heatmap(test_sign, "centered", k = 6, col_limit = 4, labelsize = 8, indicate = bla))

  test_sign_error1 <- test_sign
  SummarizedExperiment::colData(test_sign_error1) <- SummarizedExperiment::colData(test_sign_error1)[,-(3)]
  expect_error(plot_heatmap(test_sign_error1, "centered", k = 6, col_limit = 4))

  test_sign_error2 <- test_sign
  SummarizedExperiment::rowData(test_sign_error2) <- SummarizedExperiment::rowData(test_sign_error2)[,-(44)]
  expect_error(plot_heatmap(test_sign_error2, "centered", k = 6, col_limit = 4))

  test_sign_error3 <- test_sign
  SummarizedExperiment::rowData(test_sign_error3) <- SummarizedExperiment::rowData(test_sign_error3)[,-c(28:30,33:35,38:40)]
  expect_error(plot_heatmap(test_sign_error3, "centered", k = 6, col_limit = 4))
})

test_that("plot_heatmap returns a HeatmapList object", {
  expect_is(plot_heatmap(test_sign, "centered", col_limit = 4, row_font_size = 3), "HeatmapList")
  expect_is(plot_heatmap(test_sign, "contrast", col_limit = 8, row_font_size = 3), "HeatmapList")
  expect_is(plot_heatmap(test_sign, "centered", kmeans = TRUE, col_limit = 4, row_font_size = 3), "HeatmapList")
  expect_is(plot_heatmap(test_sign, "contrast", kmeans = TRUE, col_limit = 8, row_font_size = 3), "HeatmapList")
  expect_is(plot_heatmap(test_sign, "centered", kmeans = TRUE, col_limit = 8, row_font_size = 3, indicate = "replicate"), "HeatmapList")
})

test_that("plot_heatmap returns a data frame object", {
  expect_is(plot_heatmap(test_sign, "centered", col_limit = 4, row_font_size = 3, plot = FALSE), "data.frame")
  expect_is(plot_heatmap(test_sign, "contrast", col_limit = 8, row_font_size = 3, plot = FALSE), "data.frame")
  expect_is(plot_heatmap(test_sign, "centered", kmeans = TRUE, col_limit = 4, row_font_size = 3, plot = FALSE), "data.frame")
  expect_is(plot_heatmap(test_sign, "contrast", kmeans = TRUE, col_limit = 8, row_font_size = 3, plot = FALSE), "data.frame")
})

test_that("plot_volcano throws error without valid input", {
  expect_error(plot_volcano("test_sign", "Ubi6_vs_Ctrl", label_size = 5, add_names = TRUE))
  expect_error(plot_volcano(test_sign, Ubi6_vs_Ctrl, label_size = 5, add_names = TRUE))
  expect_error(plot_volcano(test_sign, "Ubi6_vs_Ctrl", label_size = 5, add_names = "TRUE"))
  expect_error(plot_volcano(test_sign, "Ubi6_vs_Ctrl", label_size = "5", add_names = TRUE))
  expect_error(plot_volcano(test_sign, "Ubi5_vs_Ctrl", label_size = 5, add_names = TRUE))

  test_sign_error1 <- test_sign
  SummarizedExperiment::rowData(test_sign_error1) <- SummarizedExperiment::rowData(test_sign_error1)[,-(1)]
  expect_error(plot_volcano(test_sign_error1, "Ubi6_vs_Ctrl", label_size = 5, add_names = TRUE))

  test_sign_error2 <- test_sign
  SummarizedExperiment::rowData(test_sign_error2) <- SummarizedExperiment::rowData(test_sign_error2)[,-c(28:30,33:35,38:40)]
  expect_error(plot_volcano(test_sign_error2, "Ubi6_vs_Ctrl", label_size = 5, add_names = TRUE))

  test_sign_error3 <- test_sign
  SummarizedExperiment::rowData(test_sign_error3) <- SummarizedExperiment::rowData(test_sign_error3)[,-(41:43)]
  expect_error(plot_volcano(test_sign_error3, "Ubi6_vs_Ctrl", label_size = 5, add_names = TRUE))
})

test_that("plot_volcano returns a ggplot object", {
  expect_is(plot_volcano(test_sign, "Ubi6_vs_Ctrl", label_size = 5, add_names = TRUE), "ggplot")
  expect_is(plot_volcano(test_sign, "Ubi6_vs_Ctrl", add_names = FALSE), "ggplot")
  expect_is(plot_volcano(test_sign, "Ubi6_vs_Ctrl", add_names = FALSE, adjusted = TRUE), "ggplot")
})
