
<!-- README.md is generated from README.Rmd. Please edit that file -->
tuckerR.mmgg
============

This package performs Three-Mode Principal Components using Tuckers Models and plot interactive Biplot.Some experiment design generated three-way or three-mode data, repeated observations of a set of attributes for a set of individuals in different conditions. The information was displayed in a three-dimensional array, and the structure of the data was explored using Three-Mode Principal Component Analysis, the Tucker-2 Model.

Installation
------------

You can install tuckerR.mmgg from github with:

``` r
# install.packages("devtools")
devtools::install_github("gusart/tuckerR_mmgg")
```

Also you can install from CRAN as a packages

``` pack_cran
install.packages("tuckerR.mmgg"")
```

Important contribution of this package
--------------------------------------

The most important contribution of this package are the interactive biplot graphics and the application of the `diffit()` function to find the best combination of components to retain.

Example
-------

The data are part of collection of local populations corresponding to different races that are conserved in the Active Germplasm Bank of INTA Pergamino Experimental Station, Argentina.A data frame with 10 characteristics of 31 maize populations in two different conditions corresponding to production areas of Buenos Aires. Since the variables are repeated in both places the data frame has a total of 20 variables, 10 for an environment and evaluated them in the other conditions.

### Compute the number of components to retain using diffit methods

The diffit method is used to apply when we need to know the axis number to be gathered in the P mode, and Q mode. The third mode, K it is related to the environment numbers. The diffit method consist on fitting each value with the Tuckle algorithm.

``` r
library(tuckerR.mmgg)
#> 
#> Attaching package: 'tuckerR.mmgg'
#> The following object is masked from 'package:graphics':
#> 
#>     plot
data(maize_pop,package = "tuckerR.mmgg")
dif_sal <- diffit(maize_pop,amb=2)
print(dif_sal)
#> $models
#>      P  Q K  S SCE_diffit iter_diffit diffit coc_diffit
#> 1    2  1 2  5      35.15         836  35.15  14.831224
#> 33   2  2 2  6      37.52           8   2.37         NA
#> 65   3  3 2  8      51.91          14  14.39   4.496875
#> 96   3  4 2  9      55.11           8   3.20         NA
#> 97   4  4 2 10      63.32          90   8.21   1.765591
#> 98   5  4 2 11      67.97          17   4.65         NA
#> 129  5  5 2 12      74.80           5   6.83   1.352475
#> 130  6  5 2 13      79.85           5   5.05   1.608280
#> 131  7  5 2 14      82.99           7   3.14         NA
#> 132  8  5 2 15      86.79           4   3.80   2.317073
#> 163  8  6 2 16      88.43           6   1.64         NA
#> 164  9  6 2 17      91.03           4   2.60   2.765957
#> 195  9  7 2 18      91.97           4   0.94         NA
#> 196 10  7 2 19      92.93           9   0.96         NA
#> 227 10  8 2 20      94.84           5   1.91   1.540323
#> 228 11  8 2 21      96.08           4   1.24   1.097345
#> 229 12  8 2 22      97.21           3   1.13   1.066038
#> 260 12  9 2 23      98.27           2   1.06   1.796610
#> 261 13  9 2 24      98.86           2   0.59   1.685714
#> 262 14  9 2 25      99.21           2   0.35   1.093750
#> 263 15  9 2 26      99.53           2   0.32   1.600000
#> 264 16  9 2 27      99.73           2   0.20   1.666667
#> 265 17  9 2 28      99.85           2   0.12   1.200000
#> 266 18  9 2 29      99.95           2   0.10   2.500000
#> 297 18 10 2 30      99.99           1   0.04   4.000000
#> 298 19 10 2 31     100.00           1   0.01        Inf
#> 299 20 10 2 32     100.00           1   0.00         NA
#> 
#> $critic_value
#> [1] 3.703704
#> 
#> attr(,"class")
#> [1] "diff"
```

In determining the number of components to be retained, the first row must be ignored, since it is the first value in the algorithm. In this case, the best combination is 3 − 3 − 2.

This is a basic example which shows you how to solve a common problem:

``` r
data(maize_pop)
output <- tucker2R(maize_pop,amb=2,stand=TRUE,nc1=3,nc2=3)
```

### Extract the core matrix.

to obtain the core matrix, type:

``` r
output$matrizG  
#>           [,1]     [,2]      [,3]      [,4]     [,5]       [,6]
#> [1,] 10.260719 1.847900  3.553432  8.380775 3.021522 -0.5999851
#> [2,] -2.014825 3.989558  3.306571 -1.322206 3.332721 -4.2685767
#> [3,] -1.290695 3.355101 -3.429868  1.325232 3.341179  3.2866310
#>           [,1]     [,2]      [,3]      [,4]     [,5]       [,6]
#> [1,] 10.260719 1.847900  3.553432  8.380775 3.021522 -0.5999851
#> [2,] -2.014825 3.989558  3.306571 -1.322206 3.332721 -4.2685767
#> [3,] -1.290695 3.355101 -3.429868  1.325232 3.341179  3.2866310
```

### The plot from output of function

``` r
plot(output) 
```

![](./tools/graphics-1.png)
