# Helper Functions for Testing
generate_random_data <- function(n) rnorm(n)
generate_time_series <- function(n, frequency) ts(rnorm(n), frequency = frequency)

# Test Suite
test_that("cross_spec.default computes the cross FFT correctly", {
  a <- generate_random_data(8)
  b <- generate_random_data(8)

  result <- cross_spec(a, b)
  expect_s3_class(result, "fftab")
  expect_true("fx" %in% colnames(result))
})

test_that("cross_spec.ts preserves time series attributes", {
  a <- generate_time_series(8, frequency = 4)
  b <- generate_time_series(8, frequency = 4)

  result <- cross_spec(a, b)
  expect_s3_class(result, "fftab")
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
  expect_s3_class(result, "fftab")
  expect_equal(attr(result, ".dim"), dim(a))
})

test_that("cross_spec.array fails for mismatched dimensions", {
  a <- array(runif(8), dim = c(2, 4))
  b <- array(runif(6), dim = c(3, 2))

  expect_error(cross_spec(a, b), "dim")
})

test_that("cross_spec.fftab respects conjugate and normalization options", {
  a <- fftab(generate_random_data(8))
  b <- fftab(generate_random_data(8))

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
  freq_idx <- which.max(Mod(fft_a)[2:(length(t) / 2)]) + 1

  # Calculate phase difference at dominant frequency
  phase_diff <- Arg(fft_b[freq_idx]) - Arg(fft_a[freq_idx])
  phase_diff <- (phase_diff + pi) %% (2 * pi) - pi # Normalize to [-π, π]

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
  set.seed(42) # For reproducibility

  n <- 256
  nn <- ceiling(n / 2)
  pdiff <- pi / 4
  t <- seq(-pi, pi, len = n)
  a <- sin(8 * t) + rnorm(length(t), sd = 0.1)
  b <- sin(8 * t + pdiff) + rnorm(length(t), sd = 0.1)

  # Perform FFT
  fft_a <- stats::fft(a)
  fft_b <- stats::fft(b)

  # Cross-Power Spectrum
  cross_power <- fft_b * Conj(fft_a)

  # Compute weighted average phase difference across positive frequencies
  weights <- Mod(cross_power[2:ceiling(n / 2)]) # Use magnitude as weights
  phase_diffs <- Arg(cross_power[2:ceiling(n / 2)]) # Phase differences

  # Calculate weighted phase shift
  weighted_phase_diff <- sum(weights * phase_diffs) / sum(weights)

  expect_equal(weighted_phase_diff, pdiff, tolerance = 0.1)

  # Apply global phase shift to fft_b
  shifted_fft_b <- complex(
    modulus = Mod(fft_b),
    argument = Arg(fft_b) + weighted_phase_diff
  )
  shifted_fft_b[2:nn] <- Conj(shifted_fft_b[2 + n - 2:nn])

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

test_that("phase_diff correlation matches time-domain correlation after alignment", {
  # Parameters for signal generation
  freq <- 2 # Frequency in Hz
  sample_rate <- 100 # Samples per second
  duration <- 1 # Duration in seconds
  phase_shift <- pi / 4 # Known phase shift in radians (45°)

  # Time vector
  t <- seq(0, duration, length.out = sample_rate * duration)

  # Generate two sine waves with known phase difference
  signal_a <- sin(2 * pi * freq * t)
  signal_b <- sin(2 * pi * freq * t + phase_shift)

  # Add small noise for realism
  set.seed(42) # Ensure reproducibility
  noise_a <- rnorm(length(t), sd = 0.001)
  noise_b <- rnorm(length(t), sd = 0.001)
  signal_a <- signal_a + noise_a
  signal_b <- signal_b + noise_b

  # Remove DC components
  signal_a <- signal_a - mean(signal_a)
  signal_b <- signal_b - mean(signal_b)

  # Compute phase difference and correlation using `phase_diff`
  result <- phase_diff(signal_a, signal_b)
  computed_phase_diff <- result[["phase_diff"]]
  frequency_domain_correlation <- result[["correlation"]]

  # Apply the phase shift in the time domain
  aligned_signal_b <- sin(2 * pi * freq * t + phase_shift - computed_phase_diff)
  aligned_signal_b <- aligned_signal_b - mean(aligned_signal_b)

  # Compute time-domain correlation
  time_domain_correlation <- cor(signal_a, aligned_signal_b)

  # Validate frequency and time-domain correlation match
  expect_equal(
    frequency_domain_correlation,
    time_domain_correlation,
    tolerance = 0.01
  )
})

test_that(".correlation computes correct normalized correlation between signals", {
  # Parameters for signal generation
  freq <- 2 # Frequency in Hz
  sample_rate <- 100 # Samples per second
  duration <- 1 # Duration in seconds
  phase_shift <- pi / 4 # Phase shift in radians (45°)

  # Time vector
  t <- seq(0, duration, length.out = sample_rate * duration)

  # Generate two sine waves
  signal_a <- sin(2 * pi * freq * t)
  signal_b <- sin(2 * pi * freq * t + phase_shift)

  # Add some noise to make the test more realistic
  noise_a <- rnorm(length(t), sd = 0.1)
  noise_b <- rnorm(length(t), sd = 0.1)
  signal_a <- signal_a + noise_a
  signal_b <- signal_b + noise_b

  # Compute FFT of the signals
  fft_a <- fftab(signal_a)
  fft_b <- fftab(signal_b)

  # Compute correlation
  computed_corr <- .correlation(fft_a, fft_b)

  # Estimate expected correlation from time domain
  expected_corr <- cor(signal_a, signal_b)

  # Compare with a tolerance for numerical stability
  expect_equal(computed_corr, expected_corr, tolerance = 0.05)
})

test_that(".shift_phase with arbitrary phase shift works in the time domain", {
  # Generate a sine wave signal
  n <- 1024
  t <- seq(0, 2 * pi, length.out = n)
  x <- sin(t) # Original sine wave signal

  # Perform FFT
  fft_x <- fftab(x)

  # Shift phase by pi/4
  shifted_fft <- .shift_phase(fft_x, pi / 4)

  # Inverse FFT to get back to time domain
  shifted_x <- ifftab(shifted_fft)

  # Expected signal (shift sine by pi/4 manually)
  expected_x <- sin(t + pi / 4)

  # Compare time-domain signals
  expect_equal(shifted_x, expected_x, tolerance = 0.01)
})
