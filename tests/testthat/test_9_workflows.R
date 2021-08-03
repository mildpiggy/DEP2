context("9 - Workflows")

test_that("LFQ throws error without valid input", {
  expect_error(LFQ("test_data", UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, "UbiLength_ExpDesign", "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, QRILC, "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", control, "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", Ctrl, filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = test,
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = Gene.names, ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = Protein.IDs, alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = "0.05", lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = "1"))

  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = "test",
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "Test", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.Names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.ids", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign[,-(2)], "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "test", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
})

test_that("LFQ returns a list object", {
  expect_is(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1), "list")
  expect_is(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1L), "list")
  expect_is(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = "Reverse",
                name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1), "list")
})

test_that("LFQ returns a list object with all expected objects", {
  expect_equal(names(LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)), c("data","se","filt","norm","imputed","diff","dep","results","param"))
  expect_equal(names(LFQ(tibble::as_tibble(test_data), UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                         name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)), c("data","se","filt","norm","imputed","diff","dep","results","param"))
  expect_equal(names(LFQ(test_data, tibble::as_tibble(UbiLength_ExpDesign), "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
                         name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)), c("data","se","filt","norm","imputed","diff","dep","results","param"))
})

test_that("report throws error without valid input", {
  res <- LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
             name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)
  expect_error(report("res"))
  expect_error(report(res[-1]))
})

test_that("report output ", {
  res <- LFQ(test_data, UbiLength_ExpDesign, "QRILC", "control", "Ctrl", filter = c("Reverse", "Potential.contaminant"),
             name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)
  expect_message(report(res), c("Render reports\n"))
})

