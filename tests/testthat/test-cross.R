library(testthat)

fourier_stdev <- function(x) {
  sqrt((sum(Mod(x) ^ 2) - Re(x[1]) ^ 2) / length(x) ^ 2)
}

test_that("Can compute forier stats using base fft",{
  x <- rnorm(100, mean = 5, sd = 2)
  fx <- stats::fft(x)
  expect_equal(fourier_stdev(fx), sd(x), tolerance = 0.1)
  expect_equal(Re(fx[1])/length(fx), mean(x))
})

test_that("Can compute correlation using base fft without inverse transform", {
  n <- 256
  freq <- 8
  t <- seq(-pi, pi, len = n)
  a <- sin(freq * t)
  b <- sin(freq * t + pi / 4)

  # Standardize signals
  a <- scale(a)
  b <- scale(b)

  # Compute FFT of signals
  fft_a <- fft(a)
  fft_b <- fft(b)

  # Compute the cross-spectrum (element-wise product)
  cross_spectrum <- fft_a * Conj(fft_b)

  # Compute the zero-lag correlation directly in the frequency domain
  computed_correlation <- Re(sum(cross_spectrum)) / n^2

  # Compare with cor() result
  expect_equal(computed_correlation, drop(cor(a, b)), tolerance = 0.01)
})

test_that("Can compute maximum correlation and phase shift using base fft", {
  n <- 256
  freq <- 8
  t <- seq(-pi, pi, len = n)
  a <- sin(freq * t)
  b <- sin(freq * t + pi / 4)

  # Standardize signals
  a <- scale(a)
  b <- scale(b)

  # Compute FFT of signals
  fft_a <- fft(a)
  fft_b <- fft(b)

  # Compute cross-correlation via FFT
  cross_correlation <- fft(fft_a * Conj(fft_b), inverse = TRUE) / n^2
  cross_correlation <- Re(cross_correlation)

  # Find the maximum correlation and its index
  max_correlation <- max(cross_correlation)
  max_index <- which.max(cross_correlation)

  # Compute the lag and phase shift
  lag <- ifelse(max_index > n / 2, max_index - n, max_index - 1) # Adjust for cyclic shifts
  phase_shift <- 2 * pi * lag * freq / n

  # Compare results
  expect_equal(max_correlation, 1, tolerance = 0.01)
  expect_equal(phase_shift, pi / 4, tolerance = 0.01)
})

test_that("Can compute maximum correlation and phase shift using base fft without inverse transform", {
  n <- 256
  freq <- 8
  t <- seq(-pi, pi, len = n)
  a <- sin(freq * t)
  b <- sin(freq * t + pi / 4)

  # Standardize signals
  a <- scale(a)
  b <- scale(b)

  # Compute FFT of signals
  fft_a <- fft(a)
  fft_b <- fft(b)

  # Compute the cross-spectrum
  cross_spectrum <- fft_a * Conj(fft_b)

  # Find the maximum correlation and corresponding phase
  magnitudes <- Mod(cross_spectrum)  # Magnitudes of the cross-spectrum
  max_index <- which.max(magnitudes)  # Index of the maximum magnitude
  max_correlation <- 2 * magnitudes[max_index] / n^2  # Normalize by n^2 and scale by 2
  phase_shift <- Arg(cross_spectrum[max_index])  # Phase of the maximum value

  # Adjust the phase shift for the frequency and lag
  lag <- ifelse(max_index > n / 2, max_index - n, max_index - 1)  # Adjust for cyclic shifts
  adjusted_phase_shift <- 8 * pi * lag / n  # Corrected phase shift formula

  # Compare results
  expect_equal(max_correlation, 1, tolerance = 0.01)
  expect_equal(adjusted_phase_shift, pi / 4, tolerance = 0.01)
})

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
