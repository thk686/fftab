library(testthat)

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
  x <- c(1 + 1i, 2 + 2i, 3 + 3i)
  expected <- stats::fft(x)
  result <- .fft(x)
  expect_equal(result, expected)
})

test_that(".fft normalization works with complex input", {
  x <- c(1 + 1i, 2 + 2i, 3 + 3i)
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

test_that(".fourier_frequencies handles edge cases", {
  # Single-element input
  expect_equal(.fourier_frequencies(1), 0)

  # Invalid input
  expect_error(.fourier_frequencies("invalid"), "non-numeric")
})

test_that(".is_complex detects complex attribute correctly", {
  obj1 <- 1:10
  attr(obj1, ".is_complex") <- TRUE
  expect_true(.is_complex(obj1))

  obj2 <- 1:10
  attr(obj2, ".is_complex") <- FALSE
  expect_false(.is_complex(obj2))

  obj3 <- 1:10
  expect_null(.is_complex(obj3))
})

test_that(".dim retrieves dimension attribute", {
  obj1 <- 1:10
  attr(obj1, ".dim") <- c(2, 5)
  expect_equal(.dim(obj1), c(2, 5))

  obj2 <- 1:10
  expect_null(.dim(obj2))
})

test_that(".is_array detects array attribute", {
  obj1 <- 1:10
  attr(obj1, ".dim") <- c(2, 5)
  expect_true(.is_array(obj1))

  obj2 <- 1:10
  expect_false(.is_array(obj2))
})

test_that(".size retrieves size attribute", {
  obj1 <- 1:10
  attr(obj1, ".size") <- 10
  expect_equal(.size(obj1), 10)

  obj2 <- 1:10
  expect_null(.size(obj2))
})

test_that(".variance computes correct variance for real input", {
  x <- rnorm(1024)
  fft_x <- fftab(x)
  expect_equal(.variance(fft_x), var(x), tolerance = 0.01)

  x <- rnorm(1024, mean = -2, sd = 2)
  fft_x <- fftab(x)
  expect_equal(.variance(fft_x), var(x), tolerance = 0.01)

  x <- rnorm(1024, mean = 3, sd = 5)
  fft_x <- fftab(x, norm = TRUE)
  expect_equal(.variance(fft_x), var(x), tolerance = 0.01)
})

test_that(".shift_phase correctly shifts phase for real input", {
  # Generate a simple sine wave
  x <- sin(seq(0, 2 * pi, length.out = 128))
  fft_x <- fftab(x)

  # Shift phase by pi/2
  shifted <- .shift_phase(fft_x, pi / 2)

  # Expected phase shift: arg should increase by pi/2 for positive frequencies
  expected <- to_polr(fft_x) |>
    dplyr::mutate(
      arg = dplyr::case_when(
        .dim_1 == 0.0 ~ arg,
        .dim_1 == .nyquist(fft_x) ~ arg,
        .dim_1 > 0.0 ~ arg + pi / 2,
        .dim_1 < 0.0 ~ arg - pi / 2,
        TRUE ~ arg
      )
    )

  expect_equal(shifted$arg, expected$arg, tolerance = 1e-6)
})

test_that(".shift_phase correctly shifts phase for complex input", {
  # Generate a complex sine wave
  x <- complex(
    real = cos(seq(0, 2 * pi, length.out = 128)),
    imaginary = sin(seq(0, 2 * pi, length.out = 128))
  )
  fft_x <- fftab(x)

  # Shift phase by pi/4
  shifted <- .shift_phase(fft_x, pi / 4)

  # Expected phase shift: Add pi/4 to all arguments
  expected <- to_polr(fft_x) |> dplyr::mutate(arg = arg + pi / 4)

  expect_equal(shifted$arg, expected$arg, tolerance = 1e-6)
})

test_that(".shift_phase correctly shifts phase in the time domain", {
  # Generate a sine wave signal
  n <- 1024
  t <- seq(0, 2 * pi, length.out = n)
  x <- sin(t) # Original sine wave signal

  # Perform FFT
  fft_x <- fftab(x)

  # Shift phase by pi/2
  shifted_fft <- .shift_phase(fft_x, pi / 2)

  # Inverse FFT to get back to time domain
  shifted_x <- ifftab(shifted_fft)

  # Expected signal after phase shift (pi/2 phase shift -> cos wave)
  expected_x <- cos(t)

  # Compare time-domain signals
  expect_equal(shifted_x, expected_x, tolerance = 0.01)
})

test_that(".shift_phase shifts a cosine wave correctly in the time domain", {
  # Generate a cosine wave signal
  n <- 1024
  t <- seq(0, 2 * pi, length.out = n)
  x <- cos(t) # Original cosine wave signal

  # Perform FFT
  fft_x <- fftab(x)

  # Shift phase by pi/2
  shifted_fft <- .shift_phase(fft_x, pi / 2)

  # Inverse FFT to get back to time domain
  shifted_x <- ifftab(shifted_fft)

  # Expected signal after phase shift (pi/2 phase shift -> -sin wave)
  expected_x <- -sin(t)

  # Compare time-domain signals
  expect_equal(shifted_x, expected_x, tolerance = 0.01)
})
