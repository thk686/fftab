# tests/testthat/test-fourier_frequencies.R

library(testthat)
library(dplyr)
library(tidyr)

# ---------------------------
# .fourier_frequencies Tests
# ---------------------------
test_that(".fourier_frequencies computes normalized frequencies", {
  # Vector input
  expect_equal(.fourier_frequencies(4), c(0, 0.25, 0.5, -0.25))

  # Larger vector input
  expect_equal(.fourier_frequencies(6), c(0, 1 / 6, 2 / 6, 3 / 6, -2 / 6, -1 / 6))
})


# ---------------------------
# fourier_frequencies.default Tests
# ---------------------------
test_that("fourier_frequencies.default computes correctly", {
  # Vector input
  res <- fourier_frequencies.default(4)
  expect_equal(res$.dim_1, c(0, 0.25, 0.5, -0.25))

  # Larger vector input
  res <- fourier_frequencies.default(6)
  expect_equal(res$.dim_1, c(0, 1 / 6, 2 / 6, 3 / 6, -2 / 6, -1 / 6))
})


# ---------------------------
# fourier_frequencies.ts Tests
# ---------------------------
test_that("fourier_frequencies.ts scales frequencies correctly", {
  ts_obj <- ts(c(1, 2, 3, 4, 5, 6),
               start = c(2000, 1),
               frequency = 12
  )

  # Check scaled frequencies
  freqs <- fourier_frequencies.ts(ts_obj)
  expected <- fourier_frequencies(6) |> dplyr::mutate(.dim_1 = 12 * .dim_1)
  expect_equal(freqs, expected)
})


# ---------------------------
# fourier_frequencies.array Tests
# ---------------------------
test_that("fourier_frequencies.array computes frequencies for arrays", {
  # 2D matrix input
  mat <- matrix(1:9, nrow = 3, ncol = 3)
  res <- fourier_frequencies(mat)
  ex <- rev(tidyr::expand_grid(
    .dim_2 = c(0, 1 / 3, -1 / 3),
    .dim_1 = c(0, 1 / 3, -1 / 3)
  ))
  expect_equal(res, ex)

  # 3D array input
  arr <- array(1:27, dim = c(3, 3, 3))
  res <- fourier_frequencies(arr)
  ex <- rev(tidyr::expand_grid(
    .dim_3 = c(0, 1 / 3, -1 / 3),
    .dim_2 = c(0, 1 / 3, -1 / 3),
    .dim_1 = c(0, 1 / 3, -1 / 3)
  ))
  expect_equal(res, ex)
})


# ---------------------------
# FFT vs fourier_frequencies Tests
# ---------------------------
test_that("fft(x) matches fourier frequencies for arrays of varying sizes", {
  dims_list <- list(
    c(2, 3),   # 2D array
    c(3, 4, 5), # 3D array
    c(1, 5, 2, 3) # 4D array
  )

  for (dims in dims_list) {
    # Create an array with the specified dimensions
    arr <- array(1:prod(dims), dim = dims)

    # Compute FFT and frequencies
    fft_res <- as.vector(fft(arr))
    dim_grid <- fourier_frequencies(arr)

    # Check that the order matches
    expect_equal(length(fft_res),
                 nrow(dim_grid),
                 info = paste("Length mismatch for dims", toString(dims))
    )

    # Verify that the frequencies match expected combinations
    expected_freqs <- rev(do.call(tidyr::expand_grid, lapply(rev(dims), .fourier_frequencies)))
    names(expected_freqs) <- names(dim_grid)
    expect_equal(dim_grid,
                 expected_freqs,
                 info = paste("Frequency mismatch for dims", toString(dims))
    )
  }
})


# ---------------------------
# Row Sine Wave Impulse Position
# ---------------------------
test_that("Impulse position matches expected frequency for row sine wave", {
  dims <- c(32, 40) # Row-major input dimensions
  dim_row <- 8      # Frequency of the sine wave along rows
  t <- seq(0, 2 * pi, length.out = dims[1]) # Time vector for rows
  x <- matrix(sin(dim_row * t), nrow = dims[1], ncol = dims[2])

  # Perform FFT
  y <- tidy_fft(x) |> to_polr()

  # Find max impulse in the FFT result
  impulse <- y[which.max(y$mod), ]

  # Assertions
  expect_equal(impulse$.dim_2, 0, tolerance = 1e-6)
  expect_equal(impulse$.dim_1, dim_row / dims[1], tolerance = 1e-6)
})


# ---------------------------
# Column Sine Wave Impulse Position
# ---------------------------
test_that("Impulse position matches expected frequency for column sine wave", {
  dims <- c(32, 40)
  dim_col <- 10
  t <- seq(0, 2 * pi, length.out = dims[2])
  x <- matrix(sin(dim_col * t), nrow = dims[1], ncol = dims[2], byrow = TRUE)

  # Perform FFT
  y <- tidy_fft(x) |> to_polr()

  # Find max impulse
  impulse <- y[which.max(y$mod), ]

  # Assertions
  expect_equal(impulse$.dim_1, 0, tolerance = 1e-6)
  expect_equal(impulse$.dim_2, dim_col / dims[2], tolerance = 1e-6)
})


# ---------------------------
# Array Impulse Position
# ---------------------------
test_that("Impulse position matches expected frequency for array", {
  dims <- c(16, 32, 64)
  clen <- c(2, 8, 2)
  x <- array(0, dim = dims)
  for (i in 1:3) {
    x <- x + slice.index(x, i) / clen[i]
  }
  y <- tidy_fft(sin(2 * pi * x)) |> to_polr()

  # Find max impulse
  impl <- y |>
    dplyr::filter(.dim_1 > 0, .dim_2 > 0, .dim_3 > 0) |>
    dplyr::filter(mod == max(mod))

  # Assertions
  expect_equal(impl$.dim_1, 1 / clen[1], tolerance = 1e-6)
  expect_equal(impl$.dim_2, 1 / clen[2], tolerance = 1e-6)
  expect_equal(impl$.dim_3, 1 / clen[3], tolerance = 1e-6)
})
