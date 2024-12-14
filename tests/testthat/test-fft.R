library(testthat)

test_that(".num_samples computes correctly", {
  # Vector input
  expect_equal(.num_samples(c(1, 2, 3)), 3)

  # Matrix input
  mat <- matrix(1:9, nrow = 3)
  expect_equal(.num_samples(mat), 3)

  # Empty input
  expect_error(.num_samples(numeric(0)), "Input must not be empty")
})

test_that(".fourier_frequencies computes normalized frequencies", {
  # Vector input
  expect_equal(.fourier_frequencies(4), c(0, 0.25, 0.5, -0.25))

  # Larger vector input
  expect_equal(.fourier_frequencies(6), c(0, 1 / 6, 2 / 6, 3 / 6, -2 / 6, -1 / 6))

  # Invalid input
  expect_error(.fourier_frequencies(1), "Minimum length is 2")
})

test_that("fourier_frequencies.default computes correctly", {
  # Vector input
  res <- fourier_frequencies.default(4)
  expect_equal(res$dim_1, c(0, 0.25, 0.5, -0.25))

  # Larger vector input
  res <- fourier_frequencies.default(6)
  expect_equal(res$dim_1, c(0, 1 / 6, 2 / 6, 3 / 6, -2 / 6, -1 / 6))
})

test_that("fourier_frequencies.ts scales frequencies correctly", {
  ts_obj <- ts(c(1, 2, 3, 4, 5, 6), start = c(2000, 1), frequency = 12)

  # Check scaled frequencies
  freqs <- fourier_frequencies.ts(ts_obj)
  expect_equal(freqs, fourier_frequencies(6) |> dplyr::mutate(dim_1 = 12 * dim_1))

  # Edge case: Single-point time series
  ts_invalid <- ts(1, start = c(2000, 1), frequency = 12)
  expect_error(fourier_frequencies.ts(ts_invalid), "Minimum length is 2")
})

test_that("fourier_frequencies.array computes frequencies for arrays", {
  # 2D matrix input
  mat <- matrix(1:9, nrow = 3, ncol = 3)
  res <- fourier_frequencies(mat)
  ex <- tidyr::expand_grid(dim_2 = c(0, 1 / 3, -1 / 3),
                           dim_1 = c(0, 1 / 3, -1 / 3))
  expect_equal(res, ex)

  # 3D array input
  arr <- array(1:27, dim = c(3, 3, 3))
  res <- fourier_frequencies(arr)
  ex <- tidyr::expand_grid(dim_3 = c(0, 1 / 3, -1 / 3),
                           dim_2 = c(0, 1 / 3, -1 / 3),
                           dim_1 = c(0, 1 / 3, -1 / 3))
  expect_equal(res, ex)
})

test_that("tidy_fft.default computes FFT results", {
  x <- c(1, 0, -1, 0)

  # Complex representation
  res <- tidy_fft(x, repr = "cplx")
  expect_true(inherits(res, "tidy_fft"))
  expect_equal(res$dim_1, fourier_frequencies(x)$dim_1)
  expect_equal(res$fx, stats::fft(x))

  # Rectangular representation
  res_rect <- tidy_fft.default(x, repr = "rect")
  expect_equal(res_rect$re, Re(stats::fft(x)))
  expect_equal(res_rect$im, Im(stats::fft(x)))

  # Polar representation
  res_polr <- tidy_fft.default(x, repr = "polr")
  expect_equal(res_polr$mod, Mod(stats::fft(x)))
  expect_equal(res_polr$arg, Arg(stats::fft(x)))
})

test_that("tidy_fft.matrix computes 2D FFT results", {
  mat <- matrix(c(1, 0, -1, 0, 1, 2, 3, 4), nrow = 2)

  # Rectangular representation
  res <- tidy_fft(mat, repr = "rect")
  expect_equal(ncol(res), 4)  # dim_1, dim_2, re, im
  expect_true("dim_1" %in% colnames(res))
  expect_true("dim_2" %in% colnames(res))

  # Check values
  fft_res <- stats::fft(mat)
  expect_equal(res$re, Re(as.vector(fft_res)))
  expect_equal(res$im, Im(as.vector(fft_res)))
})

test_that("tidy_ifft reconstructs original signal", {
  x <- c(1, 0, -1, 0)
  fft_res <- tidy_fft.default(x, repr = "cplx")

  # Reconstruct signal
  recon <- tidy_ifft(fft_res)
  expect_equal(recon, x)

  # Test with matrix
  mat <- matrix(c(1, 0, -1, 0, 1, 2, 3, 4), nrow = 2)
  fft_mat_res <- tidy_fft.array(mat, repr = "cplx")
  recon_mat <- tidy_ifft(fft_mat_res)
  expect_equal(recon_mat, mat)
})

test_that("fft(x) matches fourier frequencies for arrays of varying sizes", {
  # Test with various dimensions
  dims_list <- list(
    c(2, 3),       # 2D array
    c(3, 4, 5),    # 3D array
    c(1, 5, 2, 3)  # 4D array
  )

  for (dims in dims_list) {
    # Create an array with the specified dimensions
    arr <- array(1:prod(dims), dim = dims)

    # Compute FFT and frequencies
    fft_res <- as.vector(fft(arr))
    dim_grid <- fourier_frequencies(arr)

    # Check that the order matches
    expect_equal(length(fft_res), nrow(dim_grid),
                 info = paste("Length mismatch for dims", toString(dims)))

    # Verify that the frequencies match expected combinations
    expected_freqs <- do.call(tidyr::expand_grid, lapply(rev(dims[dims > 1]), .fourier_frequencies))
    names(expected_freqs) <- names(dim_grid)
    expect_equal(dim_grid, expected_freqs,
                 info = paste("Frequency mismatch for dims", toString(dims)))
  }
})

test_that("Impulse position matches expected frequency for row sine wave", {
  dims <- c(32, 40)  # Row-major input dimensions
  dim_row <- 8       # Frequency of the sine wave along rows
  t <- seq(0, 2 * pi, length.out = dims[1])  # Time vector for rows
  x <- matrix(sin(dim_row * t), nrow = dims[1], ncol = dims[2])  # Input sine wave

  # Perform FFT
  y <- tidy_fft(x, repr = "polr")

  # Find the frequencies grid
  dim_grid <- fourier_frequencies(x)

  # Find max impulse in the FFT result
  impulse <- y[which.max(y$mod), ]

  # Assertions
  expect_equal(impulse$dim_2, 0, tolerance = 1e-6)             # No variation along columns
  expect_equal(impulse$dim_1, dim_row / dims[1], tolerance = 1e-6)  # Expected frequency along rows
})

test_that("Impulse position matches expected frequency for column sine wave", {
  dims <- c(32, 40)  # Matrix dimensions
  dim_col <- 10      # Frequency of sine wave along columns
  t <- seq(0, 2 * pi, length.out = dims[2])  # Time vector for columns
  x <- matrix(sin(dim_col * t), nrow = dims[1], ncol = dims[2], byrow = TRUE)  # Input sine wave

  # Perform FFT
  y <- tidy_fft(x, repr = "polr")

  # Find max impulse
  impulse <- y[which.max(y$mod), ]

  # Assertions
  expect_equal(impulse$dim_1, 0, tolerance = 1e-6)             # No variation along rows
  expect_equal(impulse$dim_2, dim_col / dims[2], tolerance = 1e-6)  # Expected column frequency
})

test_that("Impulse position matches expected frequency for array", {
  dims <- c(16, 32, 64)
  clen <- c(2, 8, 2)
  x <- array(0, dim = dims)
  for (i in 1:3) {
    x <- x + slice.index(x, i) / clen[i]
  }
  y <- tidy_fft(sin(2 * pi * x), repr = "polr")

  # Find max impulse
  y |>
    dplyr::mutate(i = 1:dplyr::n()) |>
    dplyr::filter(dim_1 > 0, dim_2 > 0, dim_3 > 0) |>
    dplyr::filter(mod == max(mod)) -> impl

  # Assertions
  expect_equal(impl$dim_1, 1 / clen[1], tolerance = 1e-6)
  expect_equal(impl$dim_2, 1 / clen[2], tolerance = 1e-6)
  expect_equal(impl$dim_3, 1 / clen[3], tolerance = 1e-6)
})

