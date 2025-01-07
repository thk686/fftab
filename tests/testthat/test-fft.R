test_that("fftab.default computes FFT results correctly", {
  x <- c(1, 0, -1, 0)

  # Complex representation
  res <- fftab(x)
  expect_true(inherits(res, "fftab"))
  expect_equal(res$.dim_1, .fourier_frequencies(x))
  expect_equal(get_fx(res), stats::fft(x))

  # Rectangular representation
  res_rect <- fftab(x) |> to_rect()
  expect_equal(res_rect$re, Re(stats::fft(x)))
  expect_equal(res_rect$im, Im(stats::fft(x)))

  # Polar representation
  res_polr <- fftab(x) |> to_polr()
  expect_equal(res_polr$mod, Mod(stats::fft(x)))
  expect_equal(res_polr$arg, Arg(stats::fft(x)))
})

test_that("fftab.ts computes FFT results with proper scaling", {
  ts_obj <- ts(c(1, 0, -1, 0), frequency = 4)

  # Complex representation
  res <- fftab(ts_obj)
  expect_true(inherits(res, "fftab"))
  expect_equal(res$.dim_1, .fourier_frequencies(ts_obj) * 4) # Scaled by frequency
  expect_equal(res$fx, stats::fft(as.vector(ts_obj)))

  # Check original attributes are retained
  expect_equal(attr(res, ".tsp"), attr(ts_obj, "tsp"))
})

test_that("fftab.array computes FFT results correctly", {
  arr <- array(1:8, dim = c(2, 2, 2))

  # Complex representation
  res <- fftab(arr)
  expect_true(inherits(res, "fftab"))
  expect_equal(res$fx, as.vector(stats::fft(arr)))

  # Rectangular representation
  res_rect <- fftab(arr) |> to_rect()
  expect_equal(res_rect$re, Re(as.vector(stats::fft(arr))))
  expect_equal(res_rect$im, Im(as.vector(stats::fft(arr))))

  # Original dimensions should be retained
  expect_equal(attr(res, ".dim"), dim(arr))
})

test_that("ifftab reconstructs original signal accurately", {
  x <- c(1, 0, -1, 0)
  fft_res <- fftab(x)

  # Reconstruct signal
  recon <- ifftab(fft_res)
  expect_equal(Re(recon), x) # Ensures rounding errors are accounted for

  # Test with array
  arr <- array(1:8, dim = c(2, 2, 2))
  fft_arr_res <- fftab(arr)
  recon_arr <- ifftab(fft_arr_res)
  expect_equal(Re(recon_arr), arr)
})

test_that("error handling works as expected", {
  expect_error(fftab(list(1, 2, 3)))
  expect_error(fftab(integer(0)))
  expect_error(fftab(NULL))
})
