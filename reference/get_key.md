# Return file name keys for MERIT Hydro rasters

Generates a regex-ready key string identifying all MERIT Hydro raster
tiles that overlap with a given polygon. Tile codes follow the MERIT
Hydro naming convention, e.g. `"n30e135"` for the tile covering 30–35°N,
135–140°E.

## Usage

``` r
get_key(shape)
```

## Arguments

- shape:

  Polygon object of class `sf` defining the extent.

## Value

A single character string of pipe-separated tile codes, e.g.
`"n30e135|n30e140"`. Suitable for use with
[`grep`](https://rdrr.io/r/base/grep.html) or
[`str_detect`](https://stringr.tidyverse.org/reference/str_detect.html).

## References

MERIT Hydro: <https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_Hydro/>

## Author

Akira Terui, <hanabi0111@gmail.com>
