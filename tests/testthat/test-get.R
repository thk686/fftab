fft_obj <- tidy_fft(rnorm(32)) |> set_repr(c("cplx", "rect", "polr"))

# Tests start here
test_that("get_fx returns raw Fourier coefficients", {
  result <- get_fx(fft_obj)
  expect_type(result, "complex")
  expect_equal(result, fft_obj$fx)
})

test_that("get_fx_norm normalizes Fourier coefficients correctly", {
  result_norm <- get_fx_norm(fft_obj, norm = TRUE)
  result_non_norm <- get_fx_norm(fft_obj, norm = FALSE)

  expect_type(result_norm, "complex")
  expect_type(result_non_norm, "complex")

  # Check normalization logic
  expect_equal(result_norm, fft_obj$fx / nrow(fft_obj))
  expect_equal(result_non_norm, fft_obj$fx)
})

test_that("get_re extracts real parts", {
  result <- get_re(fft_obj)
  expect_type(result, "double")
  expect_equal(result, fft_obj$re)
})

test_that("get_im extracts imaginary parts", {
  result <- get_im(fft_obj)
  expect_type(result, "double")
  expect_equal(result, fft_obj$im)
})

test_that("get_mod extracts magnitudes", {
  result <- get_mod(fft_obj)
  expect_type(result, "double")
  expect_equal(result, fft_obj$mod)
})

test_that("get_arg extracts phase angles", {
  result <- get_arg(fft_obj)
  expect_type(result, "double")
  expect_equal(result, fft_obj$arg)
})

test_that("get_rect returns matrix with real and imaginary parts", {
  result <- get_rect(fft_obj)
  expect_true(is.matrix(result))
  expect_equal(colnames(result), c("re", "im"))
})

test_that("get_polr returns matrix with magnitude and phase", {
  result <- get_polr(fft_obj)
  expect_true(is.matrix(result))
  expect_equal(colnames(result), c("mod", "arg"))
})

# Error handling
test_that("functions handle invalid inputs", {
  expect_error(get_fx(list(1, 2, 3)))
  expect_error(get_fx_norm(NULL))
  expect_error(get_re(NA))
  expect_error(get_im(NULL))
  expect_error(get_mod(list()))
  expect_error(get_arg(NULL))
  expect_error(get_rect(NULL))
  expect_error(get_polr(NULL))
})
