# Build the bin (composite-link) matrix C

Given the wide bin boundaries \\(L_j, U_j)\\ (\\j = 1,\dots,J\\) and the
fine-grid breakpoints, this function constructs the \\J \times I\\
composite-link matrix \\C\\ such that the wide-bin probabilities are
\\\gamma = C \pi\\, where \\\pi\\ is the vector of fine-grid
probabilities.

## Usage

``` r
bin_matrix(wide_breaks, fine_breaks)
```

## Arguments

- wide_breaks:

  Numeric matrix or data frame with two columns giving the lower and
  upper limit of each wide bin, or a sorted numeric vector of length
  \\J + 1\\ when the wide bins partition \\(a, b)\\ contiguously.

- fine_breaks:

  Numeric vector of length \\I + 1\\ giving the fine-grid breakpoints.

## Value

A \\J \times I\\ numeric matrix `C`.

## Details

Element \\c\_{ji}\\ is the proportion of the \\i\\th fine interval that
falls inside the \\j\\th wide bin. When wide-bin boundaries align with
the fine grid, this gives 0/1 entries (as in Section 2 of Lambert and
Eilers, 2009). Misaligned or overlapping wide bins are handled by the
rectangle method (Section 4 of the same paper).

## Examples

``` r
fine_breaks <- seq(0, 10, by = 0.5)
# five wide bins of width 2
wide_breaks <- seq(0, 10, by = 2)
C <- bin_matrix(wide_breaks, fine_breaks)
dim(C)            # 5 x 20
#> [1]  5 20
rowSums(C)        # all 1: every wide bin fully covered
#> [1] 4 4 4 4 4
```
