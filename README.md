
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidyfft

<!-- badges: start -->

[![R-CMD-check](https://github.com/thk686/tidyfft/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/thk686/tidyfft/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of tidyfft is to make working with fft’s in R easier and more
consistent. It follows the tidy philosophy by storing output in a
tibble.

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
library(tidyfft)
(
    tidy_fft(sunspot.month) |>
    change_repr("rect", .keep = "all") |>
    change_repr("cplx", .keep = "all") -> sunspot.month.fft
)
#> # A tibble: 3,177 × 6
#>      dim_1     mod    arg      re      im fx                    
#>      <dbl>   <dbl>  <dbl>   <dbl>   <dbl> <cpl>                 
#>  1 0       165092.  0     165092.      0  165092.2000+    0.000i
#>  2 0.00378  21066.  0.852  13876.  15850.  13876.3375+15850.464i
#>  3 0.00755  16362.  1.74   -2734.  16132.  -2733.7209+16131.802i
#>  4 0.0113   19976. -2.01   -8428. -18111.  -8428.2640-18110.753i
#>  5 0.0151   14872. -3.02  -14755.  -1861. -14755.1410- 1860.878i
#>  6 0.0189   15207.  2.12   -7987.  12941.  -7986.8153+12940.843i
#>  7 0.0227    7490.  1.87   -2173.   7168.  -2173.4245+ 7168.128i
#>  8 0.0264    7830. -2.02   -3396.  -7055.  -3395.5271- 7055.326i
#>  9 0.0302    4367. -1.44     573.  -4329.    572.7342- 4329.374i
#> 10 0.0340    8564. -0.422   7812.  -3510.   7811.7219- 3510.364i
#> # ℹ 3,167 more rows
```

``` r
library(ggfortify)
library(patchwork)

ggplot(fortify(sunspot.month)) + 
  aes(x = Index, y = Data) + 
  geom_line() +
  xlab("Year") +
  ylab("Sunspot count") +
  theme_bw() -> p1

xlocs <- c(1, 0.1, 0.01)
xlabs <- c("1", "10", "100")

sunspot.month.fft |>
  dplyr::filter(dim_1 > 0) |>
  ggplot() +
    aes(x = dim_1, y = mod) +
    geom_point() +
    scale_y_continuous(trans = "log", labels = function(y) round(y)) +
    scale_x_continuous(trans = "log", breaks = xlocs, labels = xlabs) +
    xlab("Cycle duration (years)") +
    ylab("Cycle amplitude") +
    geom_smooth() +
    theme_bw() -> p2

print(p1 / p2)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
