context("2 - Missing values filtering and data normalization")

test_that("filter_missval throws error without valid input", {
  expect_error(filter_missval("test_se", 0))
  expect_error(filter_missval(test_se, "0"))
  expect_error(filter_missval(test_se, 4))
  expect_error(filter_missval(test_se, -1))

  test_se_error <- test_se
  SummarizedExperiment::rowData(test_se_error) <- SummarizedExperiment::rowData(test_se_error)[,-(24:25)]
  expect_error(filter_missval(test_se_error, 0))

  test_se_error2 <- test_se
  SummarizedExperiment::colData(test_se_error2) <- SummarizedExperiment::colData(test_se_error2)[,-(3)]
  expect_error(filter_missval(test_se_error2, 0))
})

test_that("filter_missval returns a SummarizedExperiment", {
  expect_is(filter_missval(test_se, 0), "SummarizedExperiment")
})

test_that("filter_missval returns correct number of rows", {
  expect_equal(filter_missval(test_se, 0) %>% nrow(), 203)
  expect_equal(filter_missval(test_se, 1) %>% nrow(), 243)
  expect_equal(filter_missval(test_se, 2) %>% nrow(), 393)
})

test_that("filter_proteins throws error without valid input", {
  expect_error(filter_proteins("test_se", "complete"))
  expect_error(filter_proteins(test_se, complete))
  expect_error(filter_proteins(test_se, "abc"))
  expect_error(filter_proteins(test_se, "condition", thr = NULL))
  expect_error(filter_proteins(test_se, "condition", thr = -1))
  expect_error(filter_proteins(test_se, "condition", thr = 4))
  expect_error(filter_proteins(test_se, "percentage", min = NULL))
  expect_error(filter_proteins(test_se, "percentage", min = -0.5))
  expect_error(filter_proteins(test_se, "percentage", min = 1.5))

  test_se_error <- test_se
  SummarizedExperiment::rowData(test_se_error) <- SummarizedExperiment::rowData(test_se_error)[,-(24:25)]
  expect_error(filter_proteins(test_se_error, "complete"))

  test_se_error2 <- test_se
  SummarizedExperiment::colData(test_se_error2) <- SummarizedExperiment::colData(test_se_error2)[,-(3)]
  expect_error(filter_proteins(test_se_error2, "complete"))
})

test_that("filter_proteins returns a SummarizedExperiment", {
  expect_is(filter_proteins(test_se, "complete"), "SummarizedExperiment")
  expect_is(filter_proteins(test_se, "condition", thr = 0), "SummarizedExperiment")
  expect_is(filter_proteins(test_se, "fraction", min = 0.5), "SummarizedExperiment")
})

test_that("filter_proteins returns correct number of rows", {
  expect_equal(filter_proteins(test_se, "condition", thr = 0) %>% nrow(), 203)
  expect_equal(filter_proteins(test_se, "condition", thr = 1) %>% nrow(), 243)
  expect_equal(filter_proteins(test_se, "condition", thr = 2) %>% nrow(), 393)
  expect_equal(filter_proteins(test_se, "complete") %>% nrow(), 75)
  expect_equal(filter_proteins(test_se, "fraction", min = 0.5) %>% nrow(), 152)
  expect_equal(filter_proteins(test_se, "fraction", min = 0.25) %>% nrow(), 214)
})

test_that("normalize_vsn throws error without valid input", {
  expect_error(normalize_vsn("test_filter"))
})

test_that("normalize_vsn returns a SummarizedExperiment", {
  expect_is(normalize_vsn(test_filter), "SummarizedExperiment")
})
