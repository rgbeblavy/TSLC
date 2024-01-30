
<!-- README.md is generated from README.Rmd. Please edit that file -->

# daTaxa

<!-- badges: start -->
<!-- badges: end -->

The goal of daTaxa is to take into account information on both the taxa
and genes/gene pathways in order to find differentially abundant taxa
that are associated with a difference in phenotype across two
populations of interest. This method is a two-stagewise approach that
uses an L1 penalized probit log-constrast model in each stage. The first
stage of this method is to find differentially abundant functions. Then
for each function identified as differentially abundant in the first
stage, the second stage involves subsetting the taxonomic profile to
look only at taxa associated with the given function, and applies the
same model to find differentially abundant taxa.

## Installation

You can install the development version of daTaxa like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

Here, we will illustrate how to use this package to implement the
two-stagewise approach to find differentially abundant taxa given
differentially abundant functions.

First, we will load in the file ‘tax_prof_demo.csv’, which gives
synthetic data designed to mimic the output from Metaphlan. The
following chunk of code loads this in. Note that the rows correspond to
sample IDs and the columns correspond to taxa.

``` r
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.4.4     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
#> ✔ purrr     1.0.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(daTaxa)

metaphlan.out <- read.csv(file = "tax_prof_demo.csv", header = TRUE)
rownames(metaphlan.out) <- metaphlan.out[,1]
metaphlan.out[,1] <- NULL
head(metaphlan.out)
#>         s__A s__B s__C s__D s__E s__F s__G s__H s__I s__J
#> Sample1   13    6   15    0   10    6    2    0    8    1
#> Sample2    0    9    9    8   11    0    0    0    0    2
#> Sample3    5    2    9    4    9   19    0    0    0    3
#> Sample4    5    5    9    6    0    0    0    0    4    4
#> Sample5    9    8   11    0   10    0    0    5    0    5
#> Sample6    0    0    8    5   10   13    0    6    0   20
```

Now, we will load in the file ‘functional_prof_demo.csv’, which contains
synthetic data designed to mimic the output from Humann. This will give
us information about the abundances of the pathways and the taxa
associated with each pathway. The following chunk of code does this.

``` r
humann.out <- read.csv(file = "functional_prof_demo.csv", header = TRUE)
humann.out[,1] <- NULL # Remove the index column
head(humann.out)
#>   X..Pathway Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
#> 1       pwy1       8       9       1      10       8       0       8       8
#> 2  pwy1|s__A       8       0       0       5       0       0       6       8
#> 3  pwy1|s__B       0       9       1       5       8       0       2       0
#> 4       pwy2      10       0      10       6      15       5       8      18
#> 5  pwy2|s__A       5       0       5       0       9       0       1       0
#> 6  pwy2|s__C       5       0       5       6       6       0       7      10
#>   Sample9 Sample10
#> 1       2        5
#> 2       1        4
#> 3       1        1
#> 4       0        9
#> 5       0        2
#> 6       0        0
```

In the data frame above, the rows that do not have a vertical bar “\|”
give the aggregated abundances of each of the synthetic pathways, while
the rows with a vertical bar give the abundance of the pathway
associated with a given species. With the functional profile in this
form, we can apply the function ‘get.taxa’ to find the list of taxa
associated with each pathway. The following chunk of code does this.

``` r
Fjs <- get.taxa(humann.out)
print(Fjs)
#> $pwy1
#> [1] "s__A" "s__B"
#> 
#> $pwy2
#> [1] "s__A" "s__C" "s__D"
#> 
#> $pwy3
#> [1] "s__B" "s__C" "s__D" "s__E" "s__F"
#> 
#> $pwy4
#> [1] "s__C" "s__E" "s__F"
#> 
#> $pwy5
#> [1] "s__G" "s__H" "s__I" "s__J"
```

We will need the above list for running the function to implement the
two-stagewise approach to differential abundance analysis. Before doing
this, we will need to get the functional profile in the proper form
where the samples are rows and the columns are the pathways. The
following chunk of code will help us convert the Humann output into the
correct form.

``` r
humann.out.t <- humann.out %>%
  dplyr::filter(!grepl("|", X..Pathway, fixed = TRUE)) # Select only pathways giving total counts
pwy_names <- humann.out.t$X..Pathway
rownames(humann.out.t) <- pwy_names

humann.out.t <- humann.out.t %>%
  dplyr::select(-X..Pathway) %>%
  t() %>%
  as.matrix() # Must be a matrix to run the function for the two-stagewise method

print(humann.out.t)
#>          pwy1 pwy2 pwy3 pwy4 pwy5
#> Sample1     8   10   15   17   11
#> Sample2     9    0   16   12    2
#> Sample3     1   10   24   13    3
#> Sample4    10    6    9    0    8
#> Sample5     8   15    9    6   10
#> Sample6     0    5   10   21   26
#> Sample7     8    8   15    0   42
#> Sample8     8   18   10   14   55
#> Sample9     2    0   13   13   50
#> Sample10    5    9   10    3   62
```

Now that the pathway profile is in the correct form, we will add a
pseudocount and convert this to proportions.

``` r
#humann.out.t <- transform(humann.out.t, as.numeric()) # Ensure that we have numerics
humann.out.t <- humann.out.t + 0.5 # Add pseudocount
humann.out.t <- proportions(humann.out.t, margin = 1) # Convert to proportions
print(humann.out.t)
#>                 pwy1       pwy2       pwy3        pwy4       pwy5
#> Sample1  0.133858268 0.16535433 0.24409449 0.275590551 0.18110236
#> Sample2  0.228915663 0.01204819 0.39759036 0.301204819 0.06024096
#> Sample3  0.028037383 0.19626168 0.45794393 0.252336449 0.06542056
#> Sample4  0.295774648 0.18309859 0.26760563 0.014084507 0.23943662
#> Sample5  0.168316832 0.30693069 0.18811881 0.128712871 0.20792079
#> Sample6  0.007751938 0.08527132 0.16279070 0.333333333 0.41085271
#> Sample7  0.112582781 0.11258278 0.20529801 0.006622517 0.56291391
#> Sample8  0.079069767 0.17209302 0.09767442 0.134883721 0.51627907
#> Sample9  0.031055901 0.00621118 0.16770186 0.167701863 0.62732919
#> Sample10 0.060109290 0.10382514 0.11475410 0.038251366 0.68306011
```

Now that the functional profile is in the correct form and converted to
proportions, we will add a pseudocount and convert the taxonomic profile
to proportions.

``` r
metaphlan.out <- as.matrix(metaphlan.out) # Must have a matrix to run the function for the two-stagewise method
metaphlan.out <- metaphlan.out + .5 # Add pseudocount
metaphlan.out <- proportions(metaphlan.out, margin = 1) # Convert to proportions
print(metaphlan.out)
#>                 s__A        s__B       s__C        s__D        s__E        s__F
#> Sample1  0.204545455 0.098484848 0.23484848 0.007575758 0.159090909 0.098484848
#> Sample2  0.011363636 0.215909091 0.21590909 0.193181818 0.261363636 0.011363636
#> Sample3  0.098214286 0.044642857 0.16964286 0.080357143 0.169642857 0.348214286
#> Sample4  0.144736842 0.144736842 0.25000000 0.171052632 0.013157895 0.013157895
#> Sample5  0.179245283 0.160377358 0.21698113 0.009433962 0.198113208 0.009433962
#> Sample6  0.007462687 0.007462687 0.12686567 0.082089552 0.156716418 0.201492537
#> Sample7  0.096153846 0.096153846 0.09615385 0.032051282 0.006410256 0.108974359
#> Sample8  0.077272727 0.004545455 0.16818182 0.077272727 0.095454545 0.059090909
#> Sample9  0.018072289 0.018072289 0.04216867 0.066265060 0.114457831 0.114457831
#> Sample10 0.069148936 0.015957447 0.03723404 0.079787234 0.111702128 0.005319149
#>                 s__G        s__H        s__I       s__J
#> Sample1  0.037878788 0.007575758 0.128787879 0.02272727
#> Sample2  0.011363636 0.011363636 0.011363636 0.05681818
#> Sample3  0.008928571 0.008928571 0.008928571 0.06250000
#> Sample4  0.013157895 0.013157895 0.118421053 0.11842105
#> Sample5  0.009433962 0.103773585 0.009433962 0.10377358
#> Sample6  0.007462687 0.097014925 0.007462687 0.30597015
#> Sample7  0.108974359 0.006410256 0.057692308 0.39102564
#> Sample8  0.004545455 0.068181818 0.077272727 0.36818182
#> Sample9  0.006024096 0.006024096 0.006024096 0.60843373
#> Sample10 0.005319149 0.026595745 0.005319149 0.64361702
```

Now, we just need a vector indicating group membership and a matrix of
the covariate(s) of interest to perform the differential abundance
analysis. The file ‘meta_dat_demo.csv’ gives the simulated response
variable and a simulated covariate. This can be loaded in using the
following chunk of code.

``` r
meta.dat <- read.csv(file = "meta_dat_demo.csv", header = TRUE)
meta.dat[,1] <- NULL # Remove the index column since it's redundant
head(meta.dat)
#>   Group  Covariate
#> 1     0 -1.2432097
#> 2     0  0.9918590
#> 3     0  0.8639651
#> 4     0 -1.7550385
#> 5     0  0.1488082
#> 6     1  1.1100746
```

We can then use this to extract the vector indicating group membership
and matrix containing the covariate.

``` r
y <- meta.dat[,"Group"] %>% as.vector()
X <- meta.dat[,"Covariate"] %>% as.matrix()
```

Finally, we can use these quantities to determine differentially
abundant taxa given differentially abundant functions.

``` r
res <- two.stage.da.analysis(y = y, f_prof = humann.out.t, tax_prof = metaphlan.out, Fjs = Fjs, x = X)
print(res)
#> $da.fun
#> [1] "pwy5"
#> 
#> $da.taxa
#> character(0)
#> 
#> $result
#>      Taxon   P.value Adj.p.value
#> s__A  s__A 1.0000000   1.0000000
#> s__B  s__B 1.0000000   1.0000000
#> s__C  s__C 1.0000000   1.0000000
#> s__D  s__D 1.0000000   1.0000000
#> s__E  s__E 1.0000000   1.0000000
#> s__F  s__F 1.0000000   1.0000000
#> s__G  s__G 0.7015549   0.8257163
#> s__H  s__H 0.8257163   0.8257163
#> s__I  s__I 0.4116526   0.6860877
#> s__J  s__J 0.0327565   0.1637825
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.
