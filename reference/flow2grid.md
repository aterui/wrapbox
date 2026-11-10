# Convert flow accumulation raster to stream grid

Extracts a stream grid from a flow accumulation raster by applying a
minimum drainage area threshold via
[`whitebox::wbt_extract_streams()`](https://whiteboxr.gishub.org/reference/wbt_extract_streams.html).

## Usage

``` r
flow2grid(f_acc, threshold, output)
```

## Arguments

- f_acc:

  Flow accumulation raster of class `SpatRaster`.

- threshold:

  Numeric. Minimum drainage area threshold for stream initiation. The
  unit inherits from the flow accumulation layer.

- output:

  Character. File path for the output stream grid raster.

## Value

A stream grid raster of class `SpatRaster`.

## Author

Akira Terui, <hanabi0111@gmail.com>
