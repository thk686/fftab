library(testthat)

library(testthat)

# Helper Functions for Testing
generate_random_data <- function(n) {
  rnorm(n)
}

generate_time_series <- function(n, frequency) {
  ts(rnorm(n), frequency = frequency)
}

# Begin Test Suite
test_that("cross_fft.default correctly computes the cross FFT", {
  a <- generate_random_data(8)
  b <- generate_random_data(8)

  result <- cross_fft(a, b)
  expect_s3_class(result, "tidy_fft")
  expect_true("fx" %in% colnames(result))
})

test_that("cross_fft.ts checks and preserves time series attributes", {
  a <- generate_time_series(8, frequency = 4)
  b <- generate_time_series(8, frequency = 4)

  result <- cross_fft(a, b)
  expect_s3_class(result, "tidy_fft")
  expect_equal(attr(result, ".tsp"), attr(a, "tsp"))
})

test_that("cross_fft.ts fails for mismatched frequencies", {
  a <- generate_time_series(8, frequency = 4)
  b <- generate_time_series(8, frequency = 5)

  expect_error(cross_fft(a, b), "frequency")
})

test_that("cross_fft.array checks and preserves array attributes", {
  a <- array(runif(8), dim = c(2, 4))
  b <- array(runif(8), dim = c(2, 4))

  result <- cross_fft(a, b)
  expect_s3_class(result, "tidy_fft")
  expect_equal(attr(result, ".dim"), dim(a))
})

test_that("cross_fft.array fails for mismatched dimensions", {
  a <- array(runif(8), dim = c(2, 4))
  b <- array(runif(6), dim = c(3, 2))

  expect_error(cross_fft(a, b), "dim")
})

test_that("cross_fft.tidy_fft computes the cross FFT correctly", {
  a <- tidy_fft(generate_random_data(8))
  b <- tidy_fft(generate_random_data(8))

  result <- cross_fft(a, b)
  expect_s3_class(result, "tidy_fft")
  expect_true("fx" %in% colnames(result))
})

test_that("cross_fft.tidy_fft respects the conjugate option", {
  a <- tidy_fft(generate_random_data(8))
  b <- tidy_fft(generate_random_data(8))

  result_conj <- cross_fft(a, b, conj = TRUE)
  result_no_conj <- cross_fft(a, b, conj = FALSE)

  expect_false(identical(result_conj$fx, result_no_conj$fx))
})

test_that("cross_fft.default handles edge cases gracefully", {
  # Edge case: All zeros
  a <- rep(0, 8)
  b <- rep(0, 8)
  result <- cross_fft(a, b)
  expect_true(all(result$fx == 0))

  # Edge case: Identical inputs
  a <- generate_random_data(8)
  b <- a
  result <- cross_fft(a, b)
  expect_true(all(Re(result$fx) >= 0))  # Power spectrum is non-negative
})

test_that("%cfft% operator works as expected", {
  a <- generate_random_data(8)
  b <- generate_random_data(8)

  result <- a %cfft% b
  expect_s3_class(result, "tidy_fft")
})

test_that("cross_fft.tidy_fft validates normalization and complexity", {
  a <- tidy_fft(generate_random_data(8), norm = TRUE)
  b <- tidy_fft(generate_random_data(8), norm = FALSE)
  expect_error(cross_fft(a, b), ".is_normalized")
})

test_that("cross_fft produces expected results for known inputs", {
  # Test 1: Delta function
  a <- c(1, 0, 0, 0)
  b <- c(1, 0, 0, 0)
  result <- cross_fft(a, b)
  expected_fx <- rep(1, length(a))
  expect_equal(get_re(result), expected_fx)

  # Test 2: Sine waves with same frequency
  n <- 8
  freq <- 1
  a <- sin(2 * pi * freq * (0:(n - 1)) / n)
  b <- sin(2 * pi * freq * (0:(n - 1)) / n)
  result <- cross_fft(a, b)
  expected_fx <- fft(a) * Conj(fft(b))
  expect_equal(get_fx(result), expected_fx)

  # Test 3: Sine waves with different frequencies
  freq_a <- 1
  freq_b <- 2
  a <- sin(2 * pi * freq_a * (0:(n - 1)) / n)
  b <- sin(2 * pi * freq_b * (0:(n - 1)) / n)
  result <- cross_fft(a, b)
  expect_true(all(Mod(get_fx(result)) < 1e-6))  # Near-zero for all frequencies

  # Test 4: Real and constant inputs
  a <- rep(1, 8)
  b <- rep(1, 8)
  result <- cross_fft(a, b)
  expected_fx <- fft(a) * Conj(fft(b))
  expect_equal(get_fx(result), expected_fx)

  # Test 5: Orthogonal sine and cosine waves
  freq <- 1
  a <- sin(2 * pi * freq * (0:(n - 1)) / n)
  b <- cos(2 * pi * freq * (0:(n - 1)) / n)
  result <- cross_fft(a, b)
  expect_true(all(Mod(get_fx(result)[-c(2, n)]) < 1e-6))
})

# fourier_stdev <- function(x) {
#   sqrt((sum(Mod(x) ^ 2) - Re(x[1]) ^ 2) / length(x) ^ 2)
# }
#
# test_that("Can compute forier stats using base fft",{
#   x <- rnorm(100, mean = 5, sd = 2)
#   fx <- stats::fft(x)
#   expect_equal(fourier_stdev(fx), sd(x), tolerance = 0.1)
#   expect_equal(Re(fx[1])/length(fx), mean(x))
# })
#
# test_that("Can compute correlation using base fft without inverse transform", {
#   n <- 256
#   freq <- 8
#   t <- seq(-pi, pi, len = n)
#   a <- sin(freq * t)
#   b <- sin(freq * t + pi / 4)
#
#   # Standardize signals
#   a <- scale(a)
#   b <- scale(b)
#
#   # Compute FFT of signals
#   fft_a <- fft(a)
#   fft_b <- fft(b)
#
#   # Compute the cross-spectrum (element-wise product)
#   cross_spectrum <- fft_a * Conj(fft_b)
#
#   # Compute the zero-lag correlation directly in the frequency domain
#   computed_correlation <- Re(sum(cross_spectrum)) / n^2
#
#   # Compare with cor() result
#   expect_equal(computed_correlation, drop(cor(a, b)), tolerance = 0.01)
# })
#
# test_that("Can compute maximum correlation and phase shift using base fft", {
#   n <- 256
#   freq <- 8
#   t <- seq(-pi, pi, len = n)
#   a <- sin(freq * t)
#   b <- sin(freq * t + pi / 4)
#
#   # Standardize signals
#   a <- scale(a)
#   b <- scale(b)
#
#   # Compute FFT of signals
#   fft_a <- fft(a)
#   fft_b <- fft(b)
#
#   # Compute cross-correlation via FFT
#   cross_correlation <- fft(fft_a * Conj(fft_b), inverse = TRUE) / n^2
#   cross_correlation <- Re(cross_correlation)
#
#   # Find the maximum correlation and its index
#   max_correlation <- max(cross_correlation)
#   max_index <- which.max(cross_correlation)
#
#   # Compute the lag and phase shift
#   lag <- ifelse(max_index > n / 2, max_index - n, max_index - 1) # Adjust for cyclic shifts
#   phase_shift <- 2 * pi * lag * freq / n
#
#   # Compare results
#   expect_equal(max_correlation, 1, tolerance = 0.01)
#   expect_equal(phase_shift, pi / 4, tolerance = 0.01)
# })
#
# test_that("Can compute maximum correlation and phase shift using base fft without inverse transform", {
#   n <- 256
#   freq <- 8
#   t <- seq(-pi, pi, len = n)
#   a <- sin(freq * t)
#   b <- sin(freq * t + pi / 4)
#
#   # Standardize signals
#   a <- scale(a)
#   b <- scale(b)
#
#   # Compute FFT of signals
#   fft_a <- fft(a)
#   fft_b <- fft(b)
#
#   # Compute the cross-spectrum
#   cross_spectrum <- fft_a * Conj(fft_b)
#
#   # Find the maximum correlation and corresponding phase
#   magnitudes <- Mod(cross_spectrum)  # Magnitudes of the cross-spectrum
#   max_index <- which.max(magnitudes)  # Index of the maximum magnitude
#   max_correlation <- 2 * magnitudes[max_index] / n^2  # Normalize by n^2 and scale by 2
#   phase_shift <- Arg(cross_spectrum[max_index])  # Phase of the maximum value
#
#   # Adjust the phase shift for the frequency and lag
#   lag <- ifelse(max_index > n / 2, max_index - n, max_index - 1)  # Adjust for cyclic shifts
#   adjusted_phase_shift <- 8 * pi * lag / n  # Corrected phase shift formula
#
#   # Compare results
#   expect_equal(max_correlation, 1, tolerance = 0.01)
#   expect_equal(adjusted_phase_shift, pi / 4, tolerance = 0.01)
# })

# test_that("filter_max_mod finds the correct frequency", {
#   # Input: simple tidy_fft object
#   a <- tidy_fft(c(1, 0, -1, 0))
#
#   # Compute cross-spectrum (self-correlation for simplicity)
#   cross_spec_result <- cross_spec(a, a)
#
#   # Find the frequency with the maximum magnitude
#   max_mod_result <- filter_max_mod(cross_spec_result)
#
#   # Verify the number of rows (should only return the maximum frequency)
#   expect_equal(nrow(max_mod_result), 1)
#
#   # Verify that the returned `mod` matches the maximum `mod` in the original data
#   expect_equal(get_mod(max_mod_result), max(get_mod(cross_spec_result)))
#
#   # Verify that the frequency dimension values are correct
#   expected_dim <- cross_spec_result |>
#     change_repr("polr") |>
#     dplyr::filter(dim_1 > 0) |>
#     dplyr::filter(mod == max(mod)) |>
#     dplyr::select(starts_with("dim_"))
#
#   expect_equal(max_mod_result |> dplyr::select(starts_with("dim_")), expected_dim)
# })
#
# test_that("max_correlation_phase_units computes correct correlation and phase shift", {
#   # Parameters
#   n <- 256
#   freq <- 10
#   phase_shift <- pi / 4
#   t <- seq(0, 1, length.out = n)
#
#   # Generate signals
#   x1 <- sin(2 * pi * freq * t)
#   x2 <- sin(2 * pi * freq * t + phase_shift)
#
#   # FFT
#   f1 <- tidy_fft(x1)
#   f2 <- tidy_fft(x2)
#
#   # Expected results
#   expected_correlation <- 1.0
#   expected_shift_units <- phase_shift / (2 * pi * freq)
#
#   # Compute correlation and phase shift
#   result <- max_correlation_phase_units(f1, f2, freq_scale = n)
#   expect_equal(result["max_correlation"], expected_correlation, tolerance = 1e-6)
#   expect_equal(result["shift_units"], expected_shift_units, tolerance = 1e-6)
# })
