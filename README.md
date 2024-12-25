
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidyfft <a href="https://github.com/thk686/tidyfft"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/thk686/tidyfft/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/thk686/tidyfft/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of tidyfft is to make working with fft’s in R easier and more
consistent. It follows the tidy philosophy by storing output in a
tibble.

Features:

1.  Provides frequencies in each dimension in cycles per
    sample-interval.
2.  Caches attributes sufficient to allow reconstruction of original
    object when computing the inverse transform.
3.  Easy conversion between complex, rectangular, and polar
    representations.
4.  Tracks complex versus real valued input.

## Installation

You can install the development version of tidyfft from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("thk686/tidyfft")
```

## Example

Using tidy_fft with ggplot.

``` r
tidy_fft(sunspot.month, norm = TRUE) |>
  to_rect(.keep = "all") |>
  to_polr(.keep = "all") |>
  print() ->
  ssm.fft
#> # A tibble: 3,177 × 6
#>     .dim_1 fx              re     im   mod    arg
#>      <dbl> <cpl>        <dbl>  <dbl> <dbl>  <dbl>
#>  1 0       51.96+0.00i 52.0    0     52.0   0    
#>  2 0.00378  4.37+4.99i  4.37   4.99   6.63  0.852
#>  3 0.00755 -0.86+5.08i -0.860  5.08   5.15  1.74 
#>  4 0.0113  -2.65-5.70i -2.65  -5.70   6.29 -2.01 
#>  5 0.0151  -4.64-0.59i -4.64  -0.586  4.68 -3.02 
#>  6 0.0189  -2.51+4.07i -2.51   4.07   4.79  2.12 
#>  7 0.0227  -0.68+2.26i -0.684  2.26   2.36  1.87 
#>  8 0.0264  -1.07-2.22i -1.07  -2.22   2.46 -2.02 
#>  9 0.0302   0.18-1.36i  0.180 -1.36   1.37 -1.44 
#> 10 0.0340   2.46-1.10i  2.46  -1.10   2.70 -0.422
#> # ℹ 3,167 more rows
```

``` r
ggplot(fortify(sunspot.month)) +
  geom_line(aes(x = Index, y = Data)) +
  ylab("Sunspot count") +
  xlab("Year") +
  theme_bw() ->
  p1

xlocs <- c(1, 0.1, 0.01)
xlabs <- c("1", "10", "100")

ssm.fft |>
  dplyr::filter(.dim_1 > 0) |>
  ggplot() +
  geom_point(aes(x = .dim_1, y = mod)) +
  geom_smooth(aes(x = .dim_1, y = mod)) +
  scale_y_continuous(trans = "log", labels = function(y) signif(y, 1)) +
  scale_x_continuous(trans = "log", breaks = xlocs, labels = xlabs) +
  xlab("Cycle duration (years)") +
  ylab("Mean amplitude") +
  theme_bw() ->
  p2

print(p1 / p2)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />
