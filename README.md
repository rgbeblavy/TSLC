
# TSLC

<!-- badges: start -->
<!-- badges: end -->

The goal of TSLC is to take into account information on both the taxa
and genes/gene pathways in order to find differentially abundant taxa
that are associated with a potential difference in function across two
populations of interest. This is a two-stage approach that looks at
differences in function across the population in the first stage, and
then looks at differences in the taxa associated with each
differentially abundant function in the second stage.

## Installation

You can install the development version of TSLC like so:

``` r
install_github(repo = "rgbeblavy/TSLC")
```

## Example

Here, we will illustrate how to use this package to implement the
two-stagewise approach to find differentially abundant taxa given
differentially abundant functions. For the ease of illustration, we use
a very low-dimensional data set consisting of 10 subjects, 10 taxa, and
5 gene pathways.

First, we will load in the file ‘wmg_demo_data.csv’, which gives
synthetic data on the functions, as well as the taxa contributing to
each function. In this data frame, the rows that do not have a vertical
bar “\|” give the aggregated abundances of each of the synthetic
pathways, while the rows with a vertical bar tell us the contribution of
the taxon to each function. The following chunk of code loads this file.

``` r
library(tidyverse)
library(TSLC)

wmg.dat <- read.csv(file = "wmg_demo_data.csv", header = TRUE)
wmg.dat[,1] <- NULL # Remove the index column
head(wmg.dat, 7)
```

    ##   X..Pathway Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
    ## 1       pwy1       8       9       1      10       8       0       8       8
    ## 2  pwy1|s__A       8       0       0       5       0       0       6       8
    ## 3  pwy1|s__B       0       9       1       5       8       0       2       0
    ## 4       pwy2      10       0      10       6      15       5       8      18
    ## 5  pwy2|s__A       5       0       5       0       9       0       1       0
    ## 6  pwy2|s__C       5       0       5       6       6       0       7      10
    ## 7  pwy2|s__D       0       0       0       0       0       5       0       8
    ##   Sample9 Sample10
    ## 1       2        5
    ## 2       1        4
    ## 3       1        1
    ## 4       0        9
    ## 5       0        2
    ## 6       0        0
    ## 7       0        7

Now, we just need a data frame that contains the outcome and any
covariates of interest for each of the subjects. We will store this
information in the variable ‘meta.data’. Note that the names in the
SampleID column correspond to the column names in the previous data
frame. For illustration purposes, we randomly generate this dataset.
However, this would be obtained from real data in practice.

``` r
set.seed(2025)
y <- c(rep(1,5), rep(0,5)) # response variable
x <- rnorm(n = 10, mean = 0, sd = 1) # covariate

meta.data <- data.frame("SampleID" = colnames(wmg.dat)[-1],
                        "y" = y,
                        "X" = x)
head(meta.data)
```

    ##   SampleID y          X
    ## 1  Sample1 1  0.6207567
    ## 2  Sample2 1  0.0356414
    ## 3  Sample3 1  0.7731545
    ## 4  Sample4 1  1.2724891
    ## 5  Sample5 1  0.3709754
    ## 6  Sample6 0 -0.1628543

The function get.brp.taxa() implementing our two stage approach takes
these two data frames as arguments.

``` r
res <- get.brp.taxa(wmgs.dat = wmg.dat, meta.dat = meta.data)
```

The outputs of the function get.brp.taxa() contains a list of the
following information:

1.  results for DA analysis on functions

``` r
fun.res <- res$fun.res
head(fun.res)
```

    ##                   Est        SD        pval       adjp
    ## Intercept  0.44221195 0.7864032 0.573896523 0.80345513
    ## pwy1       0.28562041 0.2337765 0.221795677 0.77628487
    ## pwy2       0.02052461 0.1939548 0.915723883 0.91572388
    ## pwy3       0.31918108 0.3977769 0.422314473 0.80345513
    ## pwy4       0.07020848 0.1928187 0.715770901 0.83506605
    ## pwy5      -0.69746487 0.2688715 0.009485406 0.06639784

2.  results for identifying BRp taxa given DA functions

``` r
tax.res <- res$tax.res
head(tax.res)
```

    ## $pwy5
    ##                  Est        SD        pval       adjp
    ## pwy5|s__G -0.2056811 0.2334870 0.378366354 0.73834437
    ## pwy5|s__H -0.1219142 0.1852466 0.510462064 0.73834437
    ## pwy5|s__I  0.1155514 0.2148366 0.590675497 0.73834437
    ## pwy5|s__J  0.4573189 0.1686515 0.006695531 0.03347766
    ## X         -0.1354806 0.5655684 0.810680951 0.81068095

3.  names of all taxa identified as BRp

``` r
brp.tax <- res$brp.tax
head(brp.tax)
```

    ## [1] "pwy5|s__J"
