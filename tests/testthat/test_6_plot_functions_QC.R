context("6 - Plot functions QC")

test_that("meanSdPlot returns a list object", {
  expect_is(meanSdPlot(test_vsn), "list")
})

test_that("plot_normalization throws error without valid input", {
  expect_error(plot_normalization(test_data))

  test_filter_error <- test_filter
  SummarizedExperiment::colData(test_filter_error) <- SummarizedExperiment::colData(test_filter_error)[,-(3)]
  expect_error(plot_normalization(test_filter_error, test_vsn))

  test_vsn_error <- test_vsn
  SummarizedExperiment::colData(test_vsn_error) <- SummarizedExperiment::colData(test_vsn_error)[,-(3)]
  expect_error(plot_normalization(test_filter, test_vsn_error))
})

test_that("plot_normalization returns a ggplot object", {
  expect_is(plot_normalization(test_filter, test_vsn), "ggplot")
})

test_that("plot_imputation throws error without valid input", {
  expect_error(plot_imputation(test_data))

  test_vsn_error <- test_vsn
  SummarizedExperiment::colData(test_vsn_error) <- SummarizedExperiment::colData(test_vsn_error)[,-(3)]
  expect_error(plot_imputation(test_vsn_error, test_impute))

  test_impute_error <- test_impute
  SummarizedExperiment::colData(test_impute_error) <- SummarizedExperiment::colData(test_impute_error)[,-(3)]
  expect_error(plot_imputation(test_vsn, test_impute_error))
})

test_that("plot_imputation returns a ggplot object", {
  expect_is(plot_imputation(test_vsn, test_impute), "ggplot")
})

test_that("plot_detect throws error without valid input", {
  expect_error(plot_detect("test_filter"))

  NAs <- apply(SummarizedExperiment::assay(test_filter), 1, function(x) any(is.na(x)))
  no_NAs <- test_filter[!NAs,]
  expect_error(plot_detect(no_NAs))
})

test_that("plot_detect returns a grob object", {
  expect_is(plot_detect(test_filter), "grob")
})

test_that("plot_missval throws error without valid input", {
  expect_error(plot_missval("test_filter"))

  NAs <- apply(SummarizedExperiment::assay(test_filter), 1, function(x) any(is.na(x)))
  no_NAs <- test_filter[!NAs,]
  expect_error(plot_missval(no_NAs))
})

test_that("plot_missval returns a HeatmapList object", {
  expect_is(plot_missval(test_filter), "HeatmapList")
})

test_that("plot_p_hist throws error without valid input", {
  expect_error(plot_p_hist("test_sign", FALSE, FALSE))
  expect_error(plot_p_hist(test_sign, "FALSE", FALSE))
  expect_error(plot_p_hist(test_sign, FALSE, "FALSE"))

  expect_error(plot_p_hist(test_filter))

  test_filter_error <- test_sign
  SummarizedExperiment::rowData(test_filter_error) <- SummarizedExperiment::rowData(test_filter_error)[,-(1)]
  expect_error(plot_p_hist(test_filter_error))

  test_filter_error2 <- test_sign
  SummarizedExperiment::rowData(test_filter_error2) <- SummarizedExperiment::rowData(test_filter_error2)[,-c(30,35,40)]
  expect_error(plot_p_hist(test_filter_error2))

})

test_that("plot_p_hist returns a ggplot object", {
  expect_is(plot_p_hist(test_sign), "ggplot")
  expect_is(plot_p_hist(test_sign, TRUE, FALSE), "ggplot")
  expect_is(plot_p_hist(test_sign, FALSE, TRUE), "ggplot")
})
