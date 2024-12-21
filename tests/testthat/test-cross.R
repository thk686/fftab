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

test_that("Can find lag and correlation for sinusoid", {
  t <- seq(-pi, pi, len = 256)
  a <- sin(8 * t)
  b <- sin(8 * t + pi / 4)

  # Perform FFT
  fft_a <- stats::fft(a)
  fft_b <- stats::fft(b)

  # Identify the dominant frequency (skip DC component)
  freq_idx <- which.max(Mod(fft_a)[2:(length(t)/2)]) + 1

  # Calculate phase difference at dominant frequency
  phase_diff <- Arg(fft_b[freq_idx]) - Arg(fft_a[freq_idx])
  phase_diff <- (phase_diff + pi) %% (2 * pi) - pi  # Normalize to [-π, π]

  # Calculate lag in radians
  lag_ab <- phase_diff

  # Calculate normalized cross-power spectrum for correlation
  cross_power <- fft_a[freq_idx] * Conj(fft_b[freq_idx])
  cor_ab <- Mod(cross_power) / (Mod(fft_a[freq_idx]) * Mod(fft_b[freq_idx]))

  expect_equal(lag_ab, pi / 4, tolerance = 0.01)
  expect_equal(cor_ab, 1, tolerance = 0.01)
})

test_that("Can calculate unshifted correlation function in frequency domain", {
  a <- rnorm(256)
  b <- rnorm(256)

  fft_a <- stats::fft(a)
  fft_b <- stats::fft(b)

  # Compute mean of signals from DC component (first FFT coefficient)
  mean_a <- Re(fft_a[1]) / length(a)
  mean_b <- Re(fft_b[1]) / length(b)

  # Remove mean (DC component set to zero)
  fft_a[1] <- 0
  fft_b[1] <- 0

  # Compute variance using Parseval's theorem
  var_a <- sum(Mod(fft_a)^2) / length(a)^2
  var_b <- sum(Mod(fft_b)^2) / length(b)^2

  # Compute cross-correlation (inner product in frequency domain)
  cross_corr <- sum(fft_a * Conj(fft_b)) / length(a)^2

  # Compute correlation coefficient
  cor_ab <- Re(cross_corr) / sqrt(var_a * var_b)

  expect_equal(cor_ab, cor(a, b), tolerance = 0.01)
})

test_that("Can calculate shifted correlation function in frequency domain", {
  set.seed(42)  # For reproducibility

  t <- seq(-pi, pi, len = 256)
  a <- sin(8 * t) + rnorm(length(t), sd = 0.25)
  b <- sin(8 * t + pi / 4) + rnorm(length(t), sd = 0.25)

  # Perform FFT
  fft_a <- stats::fft(a)
  fft_b <- stats::fft(b)

  # Cross-Power Spectrum
  cross_power <- fft_a * Conj(fft_b)

  # Compute weighted average phase difference across all frequencies
  weights <- Mod(cross_power)  # Use magnitude as weights
  phase_diffs <- Arg(cross_power)  # Phase differences

  # Avoid DC component (index 1)
  weights <- weights[-1]
  phase_diffs <- phase_diffs[-1]

  # Calculate weighted phase shift
  weighted_phase_diff <- sum(weights * phase_diffs) / sum(weights)

  # Apply global phase shift to fft_b
  shifted_fft_b <- fft_b * exp(1i * weighted_phase_diff)

  # Compute variance using Parseval's theorem
  var_a <- sum((Mod(fft_a)[-1])^2) / length(a)^2
  var_b <- sum((Mod(shifted_fft_b)[-1])^2) / length(b)^2

  # Compute cross-correlation (inner product in frequency domain)
  cross_corr <- sum(fft_a * Conj(shifted_fft_b)) / length(a)^2

  # Compute correlation coefficient
  cor_ab <- Re(cross_corr) / sqrt(var_a * var_b)

  # Inverse FFT to reconstruct shifted signal
  shifted_b <- Re(stats::fft(shifted_fft_b, inverse = TRUE) / length(a))

  expect_equal(cor_ab, cor(a, shifted_b), tolerance = 0.01)
})

