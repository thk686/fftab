library(testthat)

# tests/testthat/test-fft.R

test_that(".fft computes correct FFT without normalization", {
  x <- c(1, 2, 3, 4)
  expected <- stats::fft(x)
  result <- .fft(x)
  expect_equal(result, expected)
})

test_that(".fft computes correct FFT with normalization", {
  x <- c(1, 2, 3, 4)
  expected <- stats::fft(x) / length(x)
  result <- .fft(x, norm = TRUE)
  expect_equal(result, expected)
})

test_that(".fft handles edge cases with empty input", {
  x <- numeric(0)
  expected <- stats::fft(x)
  result <- .fft(x)
  expect_equal(result, expected)
})

test_that(".fft handles single-element input", {
  x <- c(42)
  expected <- stats::fft(x)
  result <- .fft(x)
  expect_equal(result, expected)
})

test_that(".fft handles complex input", {
  x <- c(1+1i, 2+2i, 3+3i)
  expected <- stats::fft(x)
  result <- .fft(x)
  expect_equal(result, expected)
})

test_that(".fft normalization works with complex input", {
  x <- c(1+1i, 2+2i, 3+3i)
  expected <- stats::fft(x) / length(x)
  result <- .fft(x, norm = TRUE)
  expect_equal(result, expected)
})

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
