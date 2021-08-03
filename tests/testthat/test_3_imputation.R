context("3 - Missing data imputation")

test_that("se2msn throws error without valid input", {
  expect_error(se2msn("test_se"))
})

test_that("se2msn returns a MSnSet object", {
  expect_is(se2msn(test_se), "MSnSet")
})

test_that("manual_impute throws error without valid input", {
  expect_error(manual_impute("test_vsn", 0.3, 1.8))
  expect_error(manual_impute(test_vsn, "0.3", 1.8))
  expect_error(manual_impute(test_vsn, 0.3, "1.8"))

  NAs <- apply(SummarizedExperiment::assay(test_vsn), 1, function(x) any(is.na(x)))
  no_NAs <- test_vsn[!NAs,]
  expect_error(manual_impute(no_NAs, 0.3, 1.8))
})

test_that("manual_impute returns a MSnSet object", {
  expect_is(manual_impute(test_vsn), "SummarizedExperiment")
})

test_that("manual_impute returns an object without missing values", {
  result <- SummarizedExperiment::assay(manual_impute(test_vsn))
  expect_true(all(!is.na(result)))
})

test_that("impute throws error without valid input", {
  expect_error(impute("test_vsn", "QRILC"))
  expect_error(impute(test_vsn, QRILC))
  expect_error(impute(test_vsn, "FOO"))

  test_vsn_error <- test_vsn
  SummarizedExperiment::rowData(test_vsn_error) <- SummarizedExperiment::rowData(test_vsn_error)[,-(24:25)]
  expect_error(impute(test_vsn_error, "QRILC"))

  NAs <- apply(SummarizedExperiment::assay(test_vsn), 1, function(x) any(is.na(x)))
  no_NAs <- test_vsn[!NAs,]
  expect_warning(impute(no_NAs, "QRILC"))
})

test_that("impute returns a MSnSet object", {
  expect_is(impute(test_vsn, "QRILC"), "SummarizedExperiment")
  expect_is(impute(test_vsn, "man"), "SummarizedExperiment")
})

test_that("impute returns an object without missing values", {
  result <- SummarizedExperiment::assay(impute(test_vsn, "QRILC"))
  expect_true(all(!is.na(result)))
})
