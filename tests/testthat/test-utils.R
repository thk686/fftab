library(testthat)

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

  # Minimum size input
  expect_error(.fourier_frequencies(numeric(0)), "Input must not be empty")
  expect_error(.fourier_frequencies(1), "Minimum length is 2")
})

test_that("can_repr and get_repr work correctly", {
  tft <- tibble::tibble(fx = complex(real = c(1, 2), imaginary = c(0, 1)))

  # Check representations
  expect_true(can_repr(tft, "cplx"))
  expect_false(can_repr(tft, "polr"))
  expect_equal(get_repr(tft), "cplx")
})

test_that("get_repr works correctly for multiple representations", {
  # Example with only 'cplx' representation
  tft_cplx <- tibble::tibble(
    dim_1 = c(0, 0.25, -0.5, -0.25),
    fx = complex(real = c(1, 0, -1, 0), imaginary = c(0, 1, 0, -1))
  )
  expect_equal(get_repr(tft_cplx), "cplx")

  # Example with 'polr' representation only
  tft_polr <- tibble::tibble(
    dim_1 = c(0, 0.25, -0.5, -0.25),
    mod = c(1, 1, 1, 1),
    arg = c(0, pi/2, -pi, -pi/2)
  )
  expect_equal(get_repr(tft_polr), "polr")

  # Example with 'rect' representation only
  tft_rect <- tibble::tibble(
    dim_1 = c(0, 0.25, -0.5, -0.25),
    re = c(1, 0, -1, 0),
    im = c(0, 1, 0, -1)
  )
  expect_equal(get_repr(tft_rect), "rect")

  # Example with multiple representations: 'cplx' and 'rect'
  tft_multi <- tibble::tibble(
    dim_1 = c(0, 0.25, -0.5, -0.25),
    fx = complex(real = c(1, 0, -1, 0), imaginary = c(0, 1, 0, -1)),
    re = c(1, 0, -1, 0),
    im = c(0, 1, 0, -1)
  )
  expect_true(setequal(get_repr(tft_multi), c("rect", "cplx")))

  # Example with all three representations
  tft_all <- tibble::tibble(
    dim_1 = c(0, 0.25, -0.5, -0.25),
    fx = complex(real = c(1, 0, -1, 0), imaginary = c(0, 1, 0, -1)),
    re = c(1, 0, -1, 0),
    im = c(0, 1, 0, -1),
    mod = c(1, 1, 1, 1),
    arg = c(0, pi/2, -pi, -pi/2)
  )
  expect_setequal(get_repr(tft_all), c("cplx", "rect", "polr"))
})

test_that("can_repr works correctly for single representation checks", {
  # Example with only 'cplx' representation
  tft_cplx <- tibble::tibble(
    dim_1 = c(0, 0.25, -0.5, -0.25),
    fx = complex(real = c(1, 0, -1, 0), imaginary = c(0, 1, 0, -1))
  )
  expect_true(can_repr(tft_cplx, "cplx"))
  expect_false(can_repr(tft_cplx, "polr"))
  expect_false(can_repr(tft_cplx, "rect"))

  # Example with multiple representations
  tft_multi <- tibble::tibble(
    dim_1 = c(0, 0.25, -0.5, -0.25),
    fx = complex(real = c(1, 0, -1, 0), imaginary = c(0, 1, 0, -1)),
    re = c(1, 0, -1, 0),
    im = c(0, 1, 0, -1)
  )
  expect_true(can_repr(tft_multi, "cplx"))
  expect_true(can_repr(tft_multi, "rect"))
  expect_false(can_repr(tft_multi, "polr"))

  # Example with all three representations
  tft_all <- tibble::tibble(
    dim_1 = c(0, 0.25, -0.5, -0.25),
    fx = complex(real = c(1, 0, -1, 0), imaginary = c(0, 1, 0, -1)),
    re = c(1, 0, -1, 0),
    im = c(0, 1, 0, -1),
    mod = c(1, 1, 1, 1),
    arg = c(0, pi/2, -pi, -pi/2)
  )
  expect_true(can_repr(tft_all, "cplx"))
  expect_true(can_repr(tft_all, "rect"))
  expect_true(can_repr(tft_all, "polr"))
})

test_that("can_repr and get_repr handle edge cases", {
  # Empty tibble
  tft_empty <- tibble::tibble()
  expect_equal(get_repr(tft_empty), character(0))
  expect_false(can_repr(tft_empty, "cplx"))
  expect_false(can_repr(tft_empty, "rect"))
  expect_false(can_repr(tft_empty, "polr"))

  # Tibble with unrelated columns
  tft_unrelated <- tibble::tibble(dim_1 = c(0, 0.25, -0.5, -0.25), extra = 1:4)
  expect_equal(get_repr(tft_unrelated), character(0))
  expect_false(can_repr(tft_unrelated, "cplx"))
  expect_false(can_repr(tft_unrelated, "rect"))
  expect_false(can_repr(tft_unrelated, "polr"))
})

test_that("get_fx extracts complex Fourier coefficients correctly from a polar representation", {
  # Input signal
  x <- c(1, 0, -1, 0)

  # Generate tidy_fft object
  fft_result <- tidy_fft(x) |> to_polr()  # Convert to polar representation

  # Extract complex Fourier coefficients using get_fx
  fx_values <- get_fx(fft_result)

  # Expected FFT result
  expected_fx <- fft(x)

  # Test if extracted coefficients match the expected FFT result
  expect_equal(fx_values, expected_fx)

  # Test that the result is complex
  expect_true(is.complex(fx_values))
})

