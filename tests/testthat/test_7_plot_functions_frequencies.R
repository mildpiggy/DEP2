context("7 - Plot functions frequencies")

test_that("plot_numbers throws error without valid input", {
  expect_error(plot_numbers("test_filter"))
})

test_that("plot_numbers returns a ggplot object", {
  expect_is(plot_numbers(test_filter), "ggplot")
})

test_that("plot_coverage throws error without valid input", {
  expect_error(plot_coverage("test_filter"))
})

test_that("plot_coverage returns a ggplot object", {
  expect_is(plot_coverage(test_filter), "ggplot")
})

test_that("plot_frequency throws error without valid input", {
  expect_error(plot_frequency("test_filter"))
})

test_that("plot_frequency returns a ggplot object", {
  expect_is(plot_frequency(test_filter), "ggplot")
})


test_that("plot_cond_freq throws error without valid input", {
  expect_error(plot_cond_freq("test_sign"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(44)]
  expect_error(plot_cond_freq(test_sign_error))

  test_sign_error2 <- test_sign
  SummarizedExperiment::rowData(test_sign_error2) <- SummarizedExperiment::rowData(test_sign_error2)[,-(41:43)]
  expect_error(plot_cond_freq(test_sign_error2))
})

test_that("plot_cond_freq returns a ggplot object", {
  expect_is(plot_cond_freq(test_sign), "ggplot")
})

test_that("plot_cond_overlap throws error without valid input", {
  expect_error(plot_cond_overlap("test_sign"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(44)]
  expect_error(plot_cond_overlap(test_sign_error))

  test_sign_error2 <- test_sign
  SummarizedExperiment::rowData(test_sign_error2) <- SummarizedExperiment::rowData(test_sign_error2)[,-(41:43)]
  expect_error(plot_cond_overlap(test_sign_error2))
})

test_that("plot_cond_overlap returns a grob object", {
  expect_is(plot_cond_overlap(test_sign), "grob")
})

test_that("plot_cond throws error without valid input", {
  expect_error(plot_cond("test_sign"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(44)]
  expect_error(plot_cond(test_sign_error))

  test_sign_error2 <- test_sign
  SummarizedExperiment::rowData(test_sign_error2) <- SummarizedExperiment::rowData(test_sign_error2)[,-(41:43)]
  expect_error(plot_cond(test_sign_error2))
})

test_that("plot_cond returns a grob object", {
  expect_is(plot_cond(test_sign), "grob")
})
