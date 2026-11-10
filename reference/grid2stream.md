# Convert stream grid to vector stream

Converts a raster stream grid to a vector line layer using
[`whitebox::wbt_raster_streams_to_vector()`](https://whiteboxr.gishub.org/reference/wbt_raster_streams_to_vector.html).

## Usage

``` r
grid2stream(f_dir, str_grid, output = NULL, set_crs = TRUE)
```

## Arguments

- f_dir:

  Flow direction raster of class `SpatRaster`.

- str_grid:

  Stream grid raster of class `SpatRaster`. Required if `snap = TRUE`.

- output:

  Character. Optional file path to save the output stream vector. If
  `NULL` (default), the result is returned in memory only.

- set_crs:

  Logical. Whether the output should inherit the CRS from `f_dir`.
  Defaults to `TRUE`.

## Value

A vector line layer of class `sf`.

## Author

Akira Terui, <hanabi0111@gmail.com>
