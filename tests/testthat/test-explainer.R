test_that("ubiq_explainer works for regression and classification", {
  skip_if_not_installed("ranger")
  
  # 1. Regression test (mtcars)
  rf_reg <- ranger::ranger(mpg ~ ., data = mtcars, num.trees = 20, seed = 42)
  exp_reg <- ubiq_explainer(rf_reg, mtcars[, -1], mtcars$mpg, seed = 42)
  res_reg <- explain_uncertainty(exp_reg, mtcars[1, -1], n_samples = 50)
  
  expect_s3_class(res_reg, "ubiq_explanation")
  expect_equal(nrow(res_reg$summary), ncol(mtcars) - 1)
  expect_true(all(res_reg$summary$sign_stability >= 0 & res_reg$summary$sign_stability <= 1))
  
  # 2. Multi-class classification test (iris)
  rf_cls <- ranger::ranger(Species ~ ., data = iris, num.trees = 20, probability = TRUE, seed = 42)
  exp_cls <- ubiq_explainer(rf_cls, iris[, -5], iris$Species, seed = 42)
  res_cls <- explain_uncertainty(exp_cls, iris[1, -5], n_samples = 50, class_idx = 1L)
  
  expect_s3_class(res_cls, "ubiq_explanation")
  expect_equal(res_cls$class_name, "setosa")
  
  # 3. Verify ggplot2 visualizations build without error
  p1 <- plot_uncertainty_bars(res_cls)
  p2 <- plot_uncertainty_distribution(res_cls, feature = "Petal.Length")
  p3 <- plot_uncertainty_comparison(res_cls)
  
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})