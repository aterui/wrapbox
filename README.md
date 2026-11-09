wrapbox
================

# wrapbox

`wrapbox` provides R wrappers and utilities for working with
hydrological rasters and stream networks using **terra**, **sf**, and
**WhiteboxTools**. The package includes tools for converting
flow-direction data, extracting stream networks in raster and vector
formats, and delineating nested or unnested watersheds from outlet
locations.

## Installation

Install the development/stable version from GitHub with `remotes`:

``` r
install.packages("remotes")

## development version
remotes::install_github("aterui/wrapbox")

## version-specific installation; the latest sable version is v0.2.1
remotes::install_github("aterui/wrapbox@v0.2.1")
```

Load the package:

``` r
library(wrapbox)
library(terra)
library(sf)
```

The package uses WhiteboxTools through the R package `whitebox`. If
WhiteboxTools is not yet installed on your computer, install it using:

``` r
install.packages("whitebox")
whitebox::install_whitebox()
```

## Example

See **Articles \> Getting Started** for example.

## Other utilities

See **Reference** or use `?function_name` in R to see details and
examples for each function.

## Development

Bug reports and feature requests are welcome through the [GitHub issue
tracker](https://github.com/aterui/wrapbox/issues).
