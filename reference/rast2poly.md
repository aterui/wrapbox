# Convert raster to polygon(s)

Convert raster to polygon(s)

## Usage

``` r
rast2poly(x, simplify = TRUE, keep = 0.05)
```

## Arguments

- x:

  Raster of class `SpatRaster`

- simplify:

  Logical. Whether output polygons are simplified or not.

- keep:

  Numeric. Proportion of vertices kept after polygon simplification.
  Must be between 0 and 1 (exclusive). Ignored if `simplify = FALSE`.

## Value

An `sf` object of class `MULTIPOLYGON`.

## Author

Akira Terui, <hanabi0111@gmail.com>
