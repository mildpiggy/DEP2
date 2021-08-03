context("1 - Make unique names and generate a SummarizedExperiment")

test_that("make_unique throws error without valid input", {
  expect_error(make_unique("test_data", "Gene.names", "Protein.IDs"))
  expect_error(make_unique(test_data, "Gene.name", "Protein.IDs"))
  expect_error(make_unique(test_data, "Gene.names", "Protein.ID"))

  error_test_data <- test_data
  error_test_data$Gene.names[1:20] <- NA
  error_test_data$Protein.IDs[11:20] <- NA
  expect_error(make_unique(error_test_data, "Gene.names", "Protein.IDs"))
})

test_that("make_unique returns a data.frame", {
  expect_is(make_unique(test_data, "Gene.names", "Protein.IDs"), "data.frame")
  expect_is(make_unique(tibble::as_tibble(test_data), "Gene.names", "Protein.IDs"), "data.frame")
})

test_that("make_unique returns unique names", {
  expect_false(any(duplicated(make_unique(test_data, "Gene.names", "Protein.IDs")$name)))
})

test_that("make_se throws error without valid input", {
  expect_error(make_se("test_unique", 21:32, UbiLength_ExpDesign))
  expect_error(make_se(test_unique, "21:32", UbiLength_ExpDesign))
  expect_error(make_se("test_unique", 21:32, "UbiLength_ExpDesign"))

  expect_error(make_se(test_unique, 1:10, UbiLength_ExpDesign))
  expect_error(make_se(test_unique[,-(36:37)], 21:32, UbiLength_ExpDesign))
  expect_error(make_se(test_unique, 21:32, UbiLength_ExpDesign[,-(2)]))

  error_ExpDesign <- UbiLength_ExpDesign
  error_ExpDesign$label[1:3] <- "Bla"
  expect_error(make_se(test_unique, 21:32, error_ExpDesign))
})

test_that("make_se returns a SummarizedExperiment", {
  expect_is(make_se(test_unique, 21:32, UbiLength_ExpDesign[c(7:12,1:6),]), "SummarizedExperiment")
  expect_is(make_se(test_unique, 21:32, UbiLength_ExpDesign), "SummarizedExperiment")
  expect_is(make_se(tibble::as_tibble(test_unique), 21:32, UbiLength_ExpDesign), "SummarizedExperiment")
  expect_is(make_se(test_unique, 21:32, tibble::as_tibble(UbiLength_ExpDesign)), "SummarizedExperiment")
})

test_that("get_prefix throws error without valid input", {
  words <- c("AB", "ABC", "ABCD")
  expect_error(get_prefix(10))
  expect_error(get_prefix("words"))
  expect_error(get_prefix(c(words, NA)))
})

test_that("get_prefix returns a character object", {
  expect_is(get_prefix(c("AB", "ABC", "ABCD")), "character")
})

test_that("get_prefix return the right prefix", {
  expect_equal(get_prefix(c("AB", "ABC", "ABCD")), "AB")
})

test_that("make_se_parse throws error without valid input", {
  expect_error(make_se_parse("test_unique", 21:32))
  expect_error(make_se_parse(test_unique, "21:32"))
  expect_error(make_se_parse(test_unique, 1:10))
  expect_error(make_se_parse(test_unique[,-(36:37)], 21:32))
  expect_error(make_se_parse(test_unique, 21:32, mode = "bla"))
  expect_error(make_se_parse(test_unique, 21:32, mode = bla))
})

test_that("make_se_parse returns a SummarizedExperiment", {
  expect_is(make_se_parse(test_unique, 21:32, mode = "char"), "SummarizedExperiment")
  expect_is(make_se_parse(test_unique, 21:32, mode = "delim"), "SummarizedExperiment")
  expect_is(make_se_parse(tibble::as_tibble(test_unique), 21:32), "SummarizedExperiment")
})
