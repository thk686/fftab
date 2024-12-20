library(testthat)

# Helper Functions for Testing
generate_random_data <- function(n) rnorm(n)
generate_time_series <- function(n, frequency) ts(rnorm(n), frequency = frequency)

# Test Suite
test_that("cross_spec.default computes the cross FFT correctly", {
  a <- generate_random_data(8)
  b <- generate_random_data(8)

  result <- cross_spec(a, b)
  expect_s3_class(result, "tidy_fft")
  expect_true("fx" %in% colnames(result))
})

test_that("cross_spec.ts preserves time series attributes", {
  a <- generate_time_series(8, frequency = 4)
  b <- generate_time_series(8, frequency = 4)

  result <- cross_spec(a, b)
  expect_s3_class(result, "tidy_fft")
  expect_equal(attr(result, ".tsp"), attr(a, "tsp"))
})

test_that("cross_spec.ts fails for mismatched frequencies", {
  a <- generate_time_series(8, frequency = 4)
  b <- generate_time_series(8, frequency = 5)

  expect_error(cross_spec(a, b), "frequency")
})

test_that("cross_spec.array preserves array dimensions", {
  a <- array(runif(8), dim = c(2, 4))
  b <- array(runif(8), dim = c(2, 4))

  result <- cross_spec(a, b)
  expect_s3_class(result, "tidy_fft")
  expect_equal(attr(result, ".dim"), dim(a))
})

test_that("cross_spec.array fails for mismatched dimensions", {
  a <- array(runif(8), dim = c(2, 4))
  b <- array(runif(6), dim = c(3, 2))

  expect_error(cross_spec(a, b), "dim")
})

test_that("cross_spec.tidy_fft respects conjugate and normalization options", {
  a <- tidy_fft(generate_random_data(8))
  b <- tidy_fft(generate_random_data(8))

  result_conj <- cross_spec(a, b, conj = TRUE)
  result_no_conj <- cross_spec(a, b, conj = FALSE)
  expect_false(identical(result_conj$fx, result_no_conj$fx))

  result_norm <- cross_spec(a, b, norm = TRUE)
  expect_false(identical(result_norm$fx, result_conj$fx))
})

test_that("cross_spec handles edge cases correctly", {
  a <- rep(0, 8)
  b <- rep(0, 8)
  result <- cross_spec(a, b)
  expect_true(all(result$fx == 0))

  a <- generate_random_data(8)
  b <- a
  result <- cross_spec(a, b)
  expect_true(all(Re(result$fx) >= 0))
})

test_that("cross_spec produces expected results for known inputs", {
  n <- 8

  # Test 1: Delta function
  a <- c(1, rep(0, 7))
  b <- c(1, rep(0, 7))
  result <- cross_spec(a, b)
  expected <- rep(1, n)
  expect_equal(get_re(result), expected, tolerance = 1e-6)

  # Test 2: Sine waves (same frequency)
  freq <- 1
  a <- sin(2 * pi * freq * (0:(n - 1)) / n)
  b <- sin(2 * pi * freq * (0:(n - 1)) / n)
  result <- cross_spec(a, b)
  expected <- fft(a) * Conj(fft(b))
  expect_equal(get_fx(result), expected, tolerance = 1e-6)

  # Test 3: Orthogonal sine and cosine waves
  a <- sin(2 * pi * freq * (0:(n - 1)) / n)
  b <- cos(2 * pi * freq * (0:(n - 1)) / n)
  result <- cross_spec(a, b)

  # Exclude the positive and negative frequencies of synchrony
  synchrony_index <- which(Mod(get_fx(result)) == max(Mod(get_fx(result))))
  remaining_indices <- setdiff(seq_along(get_fx(result)), synchrony_index)

  # Assert modulus is near zero for all other frequencies
  expect_true(all(Mod(get_fx(result)[remaining_indices]) < 1e-6))
})
