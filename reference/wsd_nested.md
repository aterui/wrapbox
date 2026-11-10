# Delineate nested watersheds

Delineate nested watersheds for multiple outlet points using a D8 flow
direction raster.

## Usage

``` r
wsd_nested(
  outlet,
  id_col = NULL,
  f_dir,
  str_grid = NULL,
  snap = TRUE,
  snap_dist = 5
)
```

## Arguments

- outlet:

  Outlet point layer of class `sf`.

- id_col:

  Column name specifying outlet id. This column information will be
  appended to the output polygon layer.

- f_dir:

  Flow direction raster of class `SpatRaster`.

- str_grid:

  Stream grid raster of class `SpatRaster`. Required if `snap = TRUE`.

- snap:

  Logical. Whether snapping of outlet points should be performed.
  Defaults to `TRUE`.

- snap_dist:

  Numeric. Distance threshold for snapping points to stream grid.
  Measured in the unit of input raster files. Defaults to `5`.

## Value

A polygon layer of class `sf` with one row per outlet point. If `id_col`
is supplied, the output includes the identifier column and coordinates
of the original (`x0`, `y0`) and snapped (`x`, `y`) outlet points.

## Author

Akira Terui, <hanabi0111@gmail.com>
