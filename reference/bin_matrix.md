# Build the bin (composite-link) matrix C

Builds the composite-link matrix \\C\\ such that the wide-bin
probabilities are \\\gamma = C \pi\\, where \\\pi\\ is the vector of
fine-grid probabilities. Handles both aligned and misaligned wide bins
via the rectangle method.

## Usage

``` r
bin_matrix(wide_breaks, fine_breaks)
```

## Arguments

- wide_breaks:

  Numeric matrix or data frame with two columns giving the lower and
  upper limits of each wide bin, or a sorted numeric vector of length
  \\J + 1\\ when the wide bins partition \\(a, b)\\ contiguously.

- fine_breaks:

  Numeric vector of length \\I + 1\\ giving the fine-grid breakpoints.

## Value

A \\J \times I\\ matrix `C`, where `C[j, i]` is the proportion of fine
interval \\i\\ that lies in wide bin \\j\\.

## Examples

``` r
fine  <- seq(0, 10, by = 0.5)
wide  <- seq(0, 10, by = 2)
C     <- bin_matrix(wide, fine)
dim(C)
#> [1]  5 20
```
