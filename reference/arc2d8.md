# Convert ArcGIS flow direction to D8 flow direction

Reclassifies a flow direction raster from ArcGIS D8 encoding to
WhiteboxTools D8 (flow pointer) encoding. The two formats use different
powers-of-2 assignments for the 8 cardinal and intercardinal directions:

## Usage

``` r
arc2d8(x)
```

## Arguments

- x:

  Flow direction raster of class `SpatRaster` in ArcGIS D8 encoding.
  Values 247 and 255 (ArcGIS no-data codes) are converted to `NA`.

## Value

A `SpatRaster` of integer type with WhiteboxTools D8 encoding.

## Details

**ArcGIS D8:**

     32  64 128
     16   0   1
      8   4   2

**WhiteboxTools D8:**

     64 128   1
     32   0   2
     16   8   4

## Author

Akira Terui, <hanabi0111@gmail.com>
