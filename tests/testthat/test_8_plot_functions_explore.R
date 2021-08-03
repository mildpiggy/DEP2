context("8 - Plot functions explore")

test_that("plot_pca throws error without valid input", {
  expect_error(plot_pca("test_sign", x = 1, y = 2, n = 100))
  expect_error(plot_pca(test_sign, x = "1", y = 2, n = 100))
  expect_error(plot_pca(test_sign, x = 1, y = "2", n = 100))
  expect_error(plot_pca(test_sign, x = 1, y = 2, n = "100"))
  expect_error(plot_pca(test_sign, x = 20, y = 2, n = 100))
  expect_error(plot_pca(test_sign, x = 1, y = 2, n = 500))
  expect_error(plot_pca(test_sign, x = 1, y = 2, n = 100, indicate = c("a", "b", "c", "d")))
  expect_error(plot_pca(test_sign, x = 1, y = 2, n = 100, indicate = c("a", "b", "c")))
})

test_that("plot_pca returns a ggplot object", {
  expect_is(plot_pca(test_sign, x = 1, y = 2, n = 100), "ggplot")
  expect_is(plot_pca(test_sign, x = 1, y = 2, n = 100, label = TRUE), "ggplot")
  expect_is(plot_pca(test_sign, x = 1, y = 2, n = 100, point_size = 2), "ggplot")
  expect_is(plot_pca(test_sign, x = 1, y = 2, n = 100, indicate = "condition"), "ggplot")
  expect_is(plot_pca(test_sign, x = 1, y = 2, n = 100, indicate = c("label", "replicate", "condition")), "ggplot")
})

test_that("plot_cor throws error without valid input", {
  expect_error(plot_cor("test_sign", TRUE, -1, 1, "PRGn", FALSE))
  expect_error(plot_cor(test_sign, "TRUE", -1, 1, "PRGn", FALSE))
  expect_error(plot_cor(test_sign, TRUE, "-1", 1, "PRGn", FALSE))
  expect_error(plot_cor(test_sign, TRUE, -1, "1", "PRGn", FALSE))
  expect_error(plot_cor(test_sign, TRUE, -1, 1, PRGn, FALSE))
  expect_error(plot_cor(test_sign, TRUE, -1, 1, "PRGn", "FALSE"))
  expect_error(plot_cor(test_sign, TRUE, -1, 1, "Accent", FALSE))
  expect_error(plot_cor(test_sign, TRUE, -1, 2, "PRGn", FALSE))
  expect_error(plot_cor(test_sign, TRUE, -1, 1, "PRGn", FALSE, indicate = replicate))
  expect_error(plot_cor(test_sign, TRUE, -1, 1, "PRGn", FALSE, indicate = "bla"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(44)]
  expect_error(plot_cor(test_sign_error, TRUE, -1, 1, "PRGn", FALSE))
})

test_that("plot_cor_matrix returns a ComplexHeatmap object", {
  expect_is(plot_cor(test_sign, TRUE, -1, 1, "PRGn", FALSE), "HeatmapList")
  expect_is(plot_cor(test_sign, TRUE, -1, 1, "RdBu", TRUE), "HeatmapList")
  expect_is(plot_cor(test_sign, TRUE, -1, 1, "PRGn", FALSE, indicate = "replicate"), "HeatmapList")
})

test_that("plot_dist throws error without valid input", {
  expect_error(plot_dist("test_sign", TRUE, "PRGn", FALSE))
  expect_error(plot_dist(test_sign, "TRUE", "PRGn", FALSE))
  expect_error(plot_dist(test_sign, TRUE, PRGn, FALSE))
  expect_error(plot_dist(test_sign, TRUE, "PRGn", "FALSE"))
  expect_error(plot_dist(test_sign, TRUE, "Accent", FALSE))
  expect_error(plot_dist(test_sign, TRUE, "PRGn", FALSE, indicate = replicate))
  expect_error(plot_dist(test_sign, TRUE, "PRGn", FALSE, indicate = "bla"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(44)]
  expect_error(plot_dist(test_sign_error, TRUE, "PRGn", FALSE))
})

test_that("plot_dist_matrix returns a ComplexHeatmap object", {
  expect_is(plot_dist(test_sign, TRUE, "PRGn", FALSE), "HeatmapList")
  expect_is(plot_dist(test_sign, TRUE, "RdBu", TRUE), "HeatmapList")
  expect_is(plot_dist(test_sign, TRUE, "PRGn", FALSE, indicate = "replicate"), "HeatmapList")
})
