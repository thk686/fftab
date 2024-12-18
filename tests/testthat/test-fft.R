library(testthat)

base_fft_norm <- function(x) stats::fft(x) / length(x)

test_that("tidy_fft.default computes FFT results correctly", {
  x <- c(1, 0, -1, 0)

  # Complex representation
  res <- tidy_fft(x)
  expect_true(inherits(res, "tidy_fft"))
  expect_equal(res$dim_1, .fourier_frequencies(x))
  expect_equal(get_fx(res), base_fft_norm(x))

  # Rectangular representation
  res_rect <- tidy_fft(x) |> to_rect()
  expect_equal(res_rect$re, Re(base_fft_norm(x)))
  expect_equal(res_rect$im, Im(base_fft_norm(x)))

  # Polar representation
  res_polr <- tidy_fft(x) |> to_polr()
  expect_equal(res_polr$mod, Mod(base_fft_norm(x)))
  expect_equal(res_polr$arg, Arg(base_fft_norm(x)))
})

test_that("tidy_fft.ts computes FFT results with proper scaling", {
  ts_obj <- ts(c(1, 0, -1, 0), frequency = 4)

  # Complex representation
  res <- tidy_fft(ts_obj)
  expect_true(inherits(res, "tidy_fft"))
  expect_equal(res$dim_1, .fourier_frequencies(ts_obj) * 4)  # Scaled by frequency
  expect_equal(res$fx, base_fft_norm(as.vector(ts_obj)))

  # Check original attributes are retained
  expect_equal(attr(res, ".tsp"), attr(ts_obj, "tsp"))
})

test_that("tidy_fft.array computes FFT results correctly", {
  arr <- array(1:8, dim = c(2, 2, 2))

  # Complex representation
  res <- tidy_fft(arr)
  expect_true(inherits(res, "tidy_fft"))
  expect_equal(res$fx, as.vector(base_fft_norm(arr)))

  # Rectangular representation
  res_rect <- tidy_fft(arr) |> to_rect()
  expect_equal(res_rect$re, Re(as.vector(base_fft_norm(arr))))
  expect_equal(res_rect$im, Im(as.vector(base_fft_norm(arr))))

  # Original dimensions should be retained
  expect_equal(attr(res, ".dim"), dim(arr))
})

test_that("tidy_ifft reconstructs original signal accurately", {
  x <- c(1, 0, -1, 0)
  fft_res <- tidy_fft(x)

  # Reconstruct signal
  recon <- tidy_ifft(fft_res)
  expect_equal(Re(recon), x)  # Ensures rounding errors are accounted for

  # Test with array
  arr <- array(1:8, dim = c(2, 2, 2))
  fft_arr_res <- tidy_fft(arr)
  recon_arr <- tidy_ifft(fft_arr_res)
  expect_equal(Re(recon_arr), arr)
})

test_that("error handling works as expected", {
  expect_error(tidy_fft(list(1, 2, 3)), "is.numeric(x) || is.complex(x) is not TRUE")
  # expect_error(tidy_fft(c()), "is.vector(x) is not TRUE")
})
