library(testthat)

# tests/testthat/test-fft.R

test_that(".fft computes correct FFT without normalization", {
  x <- c(1, 2, 3, 4)
  expected <- stats::fft(x)
  result <- .fft(x)
  expect_equal(result, expected)
})

test_that(".fft computes correct FFT with normalization", {
  x <- c(1, 2, 3, 4)
  expected <- stats::fft(x) / length(x)
  result <- .fft(x, norm = TRUE)
  expect_equal(result, expected)
})

test_that(".fft handles edge cases with empty input", {
  x <- numeric(0)
  expected <- stats::fft(x)
  result <- .fft(x)
  expect_equal(result, expected)
})

test_that(".fft handles single-element input", {
  x <- c(42)
  expected <- stats::fft(x)
  result <- .fft(x)
  expect_equal(result, expected)
})

test_that(".fft handles complex input", {
  x <- c(1+1i, 2+2i, 3+3i)
  expected <- stats::fft(x)
  result <- .fft(x)
  expect_equal(result, expected)
})

test_that(".fft normalization works with complex input", {
  x <- c(1+1i, 2+2i, 3+3i)
  expected <- stats::fft(x) / length(x)
  result <- .fft(x, norm = TRUE)
  expect_equal(result, expected)
})

test_that(".get_dim_cols selects columns starting with 'dim_'", {
  df <- data.frame(
    .dim_1 = 1:3,
    .dim_2 = 4:6,
    value = 7:9
  )
  result <- .get_dim_cols(df)
  expect_equal(names(result), c(".dim_1", ".dim_2"))
})

test_that(".get_dim_cols returns an empty data frame if no 'dim_' columns exist", {
  df <- data.frame(
    value = 1:3,
    category = letters[1:3]
  )
  result <- .get_dim_cols(df)
  expect_equal(ncol(result), 0)
  expect_s3_class(result, "data.frame")
})

test_that(".get_dim_cols handles tibbles correctly", {
  library(tibble)
  df <- tibble::tibble(
    .dim_1 = 1:3,
    .dim_2 = 4:6,
    other = 7:9
  )
  result <- .get_dim_cols(df)
  expect_equal(names(result), c(".dim_1", ".dim_2"))
  expect_s3_class(result, "tbl_df")
})

test_that(".get_dim_cols handles empty data frames gracefully", {
  df <- data.frame()
  result <- .get_dim_cols(df)
  expect_equal(ncol(result), 0)
  expect_s3_class(result, "data.frame")
})

# tests/testthat/test-is_normalized.R

test_that(".is_normalized detects normalization attribute correctly", {
  # Object with .is_normalized attribute set to TRUE
  obj1 <- 1:10
  attr(obj1, ".is_normalized") <- TRUE
  expect_true(.is_normalized(obj1))

  # Object with .is_normalized attribute set to FALSE
  obj2 <- 1:10
  attr(obj2, ".is_normalized") <- FALSE
  expect_false(.is_normalized(obj2))

  # Object without .is_normalized attribute
  obj3 <- 1:10
  expect_null(.is_normalized(obj3))

  # Object with non-logical .is_normalized attribute
  obj4 <- 1:10
  attr(obj4, ".is_normalized") <- "not_logical"
  expect_equal(.is_normalized(obj4), "not_logical")
})

test_that(".num_samples computes correctly", {
  # Vector input
  expect_equal(.num_samples(c(1, 2, 3)), 3)

  # Matrix input
  mat <- matrix(1:9, nrow = 3)
  expect_equal(.num_samples(mat), 3)

  # Empty input
  expect_error(.num_samples(numeric(0)), "Input must not be empty")

  # Data frame input
  df <- data.frame(a = 1:3, b = 4:6)
  expect_equal(.num_samples(df), 3)

  # Scalar input
  expect_equal(.num_samples(1), 1)
})

test_that(".fourier_frequencies computes correctly", {
  # Basic input
  expect_equal(.fourier_frequencies(4), c(0, 0.25, 0.5, -0.25))
})
