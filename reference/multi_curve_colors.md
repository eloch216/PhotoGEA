# Set of colors for plotting multiple curves

`multi_curve_colors` returns a vector of color specifications that work
reasonably well for plotting multiple curves on the same axes.

`multi_curve_line_colors` returns the same vector, but with the first
color set to be transparent. `multi_curve_point_colors` also returns the
same vector, but with all colors except the first set to transparent.
These color specifications can be helpful when plotting measured data
along with fits, allowing the data to be displayed as points and the
fits as lines.

## Usage

``` r
multi_curve_colors()

  multi_curve_line_colors()

  multi_curve_point_colors()
```

## Details

The color set was originally formed by calling the following:

`multi_curve_colors <- c( "#000000", RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Paired")[c(1:10,12)], RColorBrewer::brewer.pal(8, "Dark2") )`

## Value

A character vector with 28 elements, each of which is a hexadecimal
color specification.

## Examples

``` r
multi_curve_colors()
#>  [1] "#000000" "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F"
#>  [8] "#E5C494" "#B3B3B3" "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99"
#> [15] "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#B15928" "#1B9E77"
#> [22] "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"

multi_curve_line_colors()
#>  [1] "#00000000" "#66C2A5"   "#FC8D62"   "#8DA0CB"   "#E78AC3"   "#A6D854"  
#>  [7] "#FFD92F"   "#E5C494"   "#B3B3B3"   "#A6CEE3"   "#1F78B4"   "#B2DF8A"  
#> [13] "#33A02C"   "#FB9A99"   "#E31A1C"   "#FDBF6F"   "#FF7F00"   "#CAB2D6"  
#> [19] "#6A3D9A"   "#B15928"   "#1B9E77"   "#D95F02"   "#7570B3"   "#E7298A"  
#> [25] "#66A61E"   "#E6AB02"   "#A6761D"   "#666666"  

multi_curve_point_colors()
#>  [1] "#000000"   "#66C2A500" "#FC8D6200" "#8DA0CB00" "#E78AC300" "#A6D85400"
#>  [7] "#FFD92F00" "#E5C49400" "#B3B3B300" "#A6CEE300" "#1F78B400" "#B2DF8A00"
#> [13] "#33A02C00" "#FB9A9900" "#E31A1C00" "#FDBF6F00" "#FF7F0000" "#CAB2D600"
#> [19] "#6A3D9A00" "#B1592800" "#1B9E7700" "#D95F0200" "#7570B300" "#E7298A00"
#> [25] "#66A61E00" "#E6AB0200" "#A6761D00" "#66666600"
```
