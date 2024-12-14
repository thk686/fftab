test_that("fourier_frequencies.default computes correctly", {
  # Basic case
  expect_equal(fourier_frequencies.default(4), c(0, 0.25, 0.5, -0.25))

  # Larger input
  expect_equal(
    fourier_frequencies.default(6),
    c(0, 1/6, 2/6, 3/6, -2/6, -1/6)
  )

  # Invalid input
  expect_error(fourier_frequencies.default(1), "Minimum length is 2")
})

test_that("fourier_frequencies.ts scales correctly with time series frequency", {
  ts_obj <- ts(c(1, 2, 3, 4, 5, 6), start = c(2000, 1), frequency = 12)

  # Check correct scaling
  freqs <- fourier_frequencies.ts(ts_obj)
  expect_equal(freqs, fourier_frequencies.default(6) * 12)

  # Edge case: Single-point time series
  ts_invalid <- ts(1, start = c(2000, 1), frequency = 12)
  expect_error(fourier_frequencies.ts(ts_invalid), "Minimum length is 2")
})

test_that("tidy_fft.default produces correct structure", {
  x <- c(1, 0, -1, 0)

  # Test "cplx" representation
  res <- tidy_fft.default(x, repr = "cplx")
  expect_true(inherits(res, "tidy_fft"))
  expect_equal(res$frequency, fourier_frequencies.default(x))
  expect_equal(res$fx, stats::fft(x))

  # Test "rect" representation
  res_rect <- tidy_fft.default(x, repr = "rect")
  expect_equal(res_rect$re, Re(stats::fft(x)))
  expect_equal(res_rect$im, Im(stats::fft(x)))

  # Test "polr" representation
  res_polr <- tidy_fft.default(x, repr = "polr")
  expect_equal(res_polr$mod, Mod(stats::fft(x)))
  expect_equal(res_polr$arg, Arg(stats::fft(x)))
})

test_that("tidy_ifft correctly reconstructs the original signal", {
  x <- c(1, 0, -1, 0)
  fft_res <- tidy_fft.default(x, repr = "cplx")

  # Reconstruct signal
  recon <- tidy_ifft(fft_res)
  expect_equal(recon, x)

  # Check real-only case
  expect_true(is.numeric(recon))

  # Test with complex input
  complex_input <- tidy_fft.default(complex(real = x, imaginary = x), repr = "cplx")
  recon_complex <- tidy_ifft(complex_input)
  expect_true(is.complex(recon_complex))
})

test_that("change_repr correctly switches representations", {
  x <- c(1, 0, -1, 0)
  fft_res <- tidy_fft.default(x, repr = "cplx")

  # Switch to rectangular
  rect <- change_repr(fft_res, "rect")
  expect_equal(rect$re, Re(fft_res$fx))
  expect_equal(rect$im, Im(fft_res$fx))

  # Switch to polar
  polr <- change_repr(rect, "polr")
  expect_equal(polr$mod, Mod(fft_res$fx))
  expect_equal(polr$arg, Arg(fft_res$fx))
})

