context("shiny app")

test_that("run_app throws error without valid input", {
  expect_error(run_app("test"))
})

