context("4 - Differential testing")

test_that("test_diff throws error without valid input", {
  expect_error(test_diff("test_impute", "control", "Ctrl"))
  expect_error(test_diff(test_impute, "control", Ctrl))
  expect_error(test_diff(test_impute, control, "Ctrl"))
  expect_error(test_diff(test_impute, "CONTROL", "Ctrl"))
  expect_error(test_diff(test_impute, "control", "control"))
  expect_error(test_diff(test_impute, "manual"))
  expect_error(test_diff(test_impute, "manual", test = Ubi4_vs_Ctrl))
  expect_error(test_diff(test_impute, "manual", test = "test"))
  expect_error(test_diff(test_impute, "manual", test = "Ubi4_vs_Ctrl", incl_repl = "TRUE"))

  test_impute_error <- test_impute
  SummarizedExperiment::colData(test_impute_error) <- SummarizedExperiment::colData(test_impute_error)[,-(3)]
  expect_error(test_diff(test_impute_error, "control", "Ctrl"))

  test_impute_error2 <- test_impute
  SummarizedExperiment::rowData(test_impute_error2) <- SummarizedExperiment::rowData(test_impute_error2)[,-(24:25)]
  expect_error(test_diff(test_impute_error2, "control", "Ctrl"))
})

test_that("test_diff returns a SummarizedExperiment object", {
  expect_is(test_diff(test_impute, "control", "Ctrl"), "SummarizedExperiment")
  expect_is(test_diff(test_impute, "manual", test = "Ubi4_vs_Ctrl"), "SummarizedExperiment")
  expect_is(test_diff(test_impute, "all", design_formula = formula(~ 0 + condition + replicate)), "SummarizedExperiment")
})

test_that("test_diff returns an object with diff, p.val and p.adj columns", {
  result <- SummarizedExperiment::rowData(test_diff(test_impute, "control", "Ctrl"))
  expect_equal(grep("_diff$", colnames(result)), c(28,33,38))
  expect_equal(grep("_p.adj$", colnames(result)), c(29,34,39))
  expect_equal(grep("_p.val$", colnames(result)), c(30,35,40))
})

test_that("add_rejections throws error without valid input", {
  expect_error(add_rejections("test_lm", 0.05, 1))
  expect_error(add_rejections(test_lm, "0.05", 1))
  expect_error(add_rejections(test_lm, 0.05, "1"))
  expect_error(add_rejections(test_impute, 0.05, 1))

  test_lm_error <- test_lm
  SummarizedExperiment::rowData(test_lm_error) <- SummarizedExperiment::rowData(test_lm_error)[,-(24:25)]
  expect_error(add_rejections(test_lm_error, 0.05, 1))
})

test_that("add_rejections returns a SummarizedExperiment object", {
  expect_is(add_rejections(test_lm, 0.05, 1), "SummarizedExperiment")

  test_lm2 <- test_lm
  SummarizedExperiment::rowData(test_lm2) <- SummarizedExperiment::rowData(test_lm2)[,-c(33:34,38:39)]
  expect_is(add_rejections(test_lm2, 0.05, 1), "SummarizedExperiment")
})

test_that("add_rejections returns an object with significance columns", {
  result <- SummarizedExperiment::rowData(add_rejections(test_lm, 0.05, 1))
  expect_equal(grep("^significant$", colnames(result)), 44)
  expect_equal(grep("_significant$", colnames(result)), 41:43)
  expect_equal(nrow(result[result$significant,]), 45)
})

test_that("get_results throws error without valid input", {
  expect_error(get_results("test_sign"))
  expect_error(get_results(test_impute))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(24:25)]
  expect_error(get_results(test_sign_error))
})

test_that("get_results returns a data.frame", {
  expect_is(get_results(test_sign), "data.frame")
})

test_that("get_results returns a data.frame with the expected column", {
  result <- get_results(test_sign)
  expect_equal(nrow(result), 203)
  expect_equal(ncol(result), 19)
  expect_equal(grep("_centered$", colnames(result)), 16:19)
  expect_equal(grep("_ratio$", colnames(result)), 13:15)
  expect_equal(nrow(result[result$significant,]), 45)
})

test_that("get_df_wide throws error without valid input", {
  expect_error(get_df_wide("test_sign"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(1)]
  expect_error(get_df_wide(test_sign_error))
})

test_that("get_df_wide returns a data.frame", {
  expect_is(get_df_wide(test_sign), "data.frame")
})

test_that("get_df_long throws error without valid input", {
  expect_error(get_df_long("test_sign"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(1)]
  expect_error(get_df_long(test_sign_error))
})

test_that("get_df_long returns a data.frame", {
  expect_is(get_df_long(test_sign), "data.frame")
})


test_that("p-values are ordered correctly", {
  
  shuffled_order <- sample(seq_len(nrow(test_impute)))
  test_shuffled <- test_impute[shuffled_order, ]
  pval1 <- rowData(test_diff(test_impute, "manual", test = "Ubi4_vs_Ctrl"))$Ubi4_vs_Ctrl_p.val
  pval2 <- rowData(test_diff(test_shuffled, "manual", test = "Ubi4_vs_Ctrl"))$Ubi4_vs_Ctrl_p.val
  
  expect_false(all(pval1 == pval2))
  expect_equal(pval1[shuffled_order], pval2)
  
  # Compare ordered against limma output
  design <- model.matrix(~ 0 + test_impute$condition)
  colnames(design) <- levels(as.factor(test_impute$condition))
  contrasts <- makeContrasts(Ubi4 - Ctrl, levels=design)
  fit1 <- lmFit(assay(test_impute), design) 
  fit2 <- contrasts.fit(fit1, contrasts)
  fit3 <- eBayes(fit2)
  expect_equal(pval1, topTable(fit3, sort.by = "none", number=Inf)$P.Value)  
  
  # Compare unordered against limma output
  design <- model.matrix(~ 0 + test_shuffled$condition)
  colnames(design) <- levels(as.factor(test_impute$condition))
  contrasts <- makeContrasts(Ubi4 - Ctrl, levels=design)
  fit1 <- lmFit(assay(test_shuffled), design) 
  fit2 <- contrasts.fit(fit1, contrasts)
  fit3 <- eBayes(fit2)
  expect_equal(pval2, topTable(fit3, sort.by = "none", number=Inf)$P.Value)  
  
})