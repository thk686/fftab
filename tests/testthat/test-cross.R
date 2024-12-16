library(testthat)

test_that("filter_max_mod finds the correct frequency", {
  # Input: simple tidy_fft object
  a <- tidy_fft(c(1, 0, -1, 0))

  # Compute cross-spectrum (self-correlation for simplicity)
  cross_spec_result <- cross_spec(a, a)

  # Find the frequency with the maximum magnitude
  max_mod_result <- filter_max_mod(cross_spec_result)

  # Verify the number of rows (should only return the maximum frequency)
  expect_equal(nrow(max_mod_result), 1)

  # Verify that the returned `mod` matches the maximum `mod` in the original data
  expect_equal(get_mod(max_mod_result), max(get_mod(cross_spec_result)))

  # Verify that the frequency dimension values are correct
  expected_dim <- cross_spec_result |>
    change_repr("polr") |>
    dplyr::filter(dim_1 > 0) |>
    dplyr::filter(mod == max(mod)) |>
    dplyr::select(starts_with("dim_"))

  expect_equal(max_mod_result |> dplyr::select(starts_with("dim_")), expected_dim)
})

test_that("max_correlation_phase_units computes correct correlation and phase shift", {
  # Parameters
  n <- 256
  freq <- 10
  phase_shift <- pi / 4
  t <- seq(0, 1, length.out = n)

  # Generate signals
  x1 <- sin(2 * pi * freq * t)
  x2 <- sin(2 * pi * freq * t + phase_shift)

  # FFT
  f1 <- tidy_fft(x1)
  f2 <- tidy_fft(x2)

  # Expected results
  expected_correlation <- 1.0
  expected_shift_units <- phase_shift / (2 * pi * freq)

  # Compute correlation and phase shift
  result <- max_correlation_phase_units(f1, f2, freq_scale = n)
  expect_equal(result["max_correlation"], expected_correlation, tolerance = 1e-6)
  expect_equal(result["shift_units"], expected_shift_units, tolerance = 1e-6)
})
