library(testthat)

test_that(".num_samples computes correctly", {
  # Vector input
  expect_equal(.num_samples(c(1, 2, 3)), 3)

  # Matrix input
  mat <- matrix(1:9, nrow = 3)
  expect_equal(.num_samples(mat), 3)

  # Empty input
  expect_error(.num_samples(numeric(0)), "Input must not be empty")

  # Data frame input
  df <- data.frame(a = 1:3, b = 4:6)
  expect_equal(.num_samples(df), 3)

  # Scalar input
  expect_equal(.num_samples(1), 1)
})

test_that(".fourier_frequencies computes correctly", {
  # Basic input
  expect_equal(.fourier_frequencies(4), c(0, 0.25, 0.5, -0.25))
})


