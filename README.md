
<!-- README.md is generated from README.Rmd. Please edit that file -->

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
5 gene pathways. The first 5 subjects belong to one population, while
the last 5 subjects belong to another population. The first 4 pathways
are not differentially abundant across the two populations, while the
last pathway is differentially abundant across the two populations. The
taxa are labeled A through J, and only species J is truly differentially
abundant. In practice, this method is designed to be applied to a
high-dimensional taxonomic and functional profile.

First, we will load in the file ‘tax_prof_demo.csv’, which gives
synthetic data designed to mimic the output from Metaphlan. The
following chunk of code loads this in. Note that the rows correspond to
sample IDs and the columns correspond to taxa.

``` r
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
#> ✔ purrr     1.0.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(TSLC)

### Load demo taxonomic profile
metaphlan.out <- read.csv(file = "tax_prof_demo.csv", header = TRUE) %>%
  column_to_rownames(var = "X") %>%
  as.matrix()

### Add psuedocount and convert to proportions
metaphlan.out <- ifelse(metaphlan.out == 0, 0.5, metaphlan.out) %>%
  proportions(margin = 1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample_Name")
head(metaphlan.out)
#>   Sample_Name       s__A       s__B      s__C        s__D       s__E       s__F
#> 1     Sample1 0.20967742 0.09677419 0.2419355 0.008064516 0.16129032 0.09677419
#> 2     Sample2 0.01204819 0.21686747 0.2168675 0.192771084 0.26506024 0.01204819
#> 3     Sample3 0.09523810 0.03809524 0.1714286 0.076190476 0.17142857 0.36190476
#> 4     Sample4 0.14285714 0.14285714 0.2571429 0.171428571 0.01428571 0.01428571
#> 5     Sample5 0.18000000 0.16000000 0.2200000 0.010000000 0.20000000 0.01000000
#> 6     Sample6 0.00781250 0.00781250 0.1250000 0.078125000 0.15625000 0.20312500
#>         s__G        s__H       s__I       s__J
#> 1 0.03225806 0.008064516 0.12903226 0.01612903
#> 2 0.01204819 0.012048193 0.01204819 0.04819277
#> 3 0.00952381 0.009523810 0.00952381 0.05714286
#> 4 0.01428571 0.014285714 0.11428571 0.11428571
#> 5 0.01000000 0.100000000 0.01000000 0.10000000
#> 6 0.00781250 0.093750000 0.00781250 0.31250000
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
  dplyr::filter(!grepl("|", X..Pathway, fixed = TRUE)) # select only pathways giving total counts
pwy_names <- humann.out.t$X..Pathway
rownames(humann.out.t) <- pwy_names

humann.out.t <- humann.out.t %>%
  dplyr::select(-X..Pathway) %>% t()

### Add psuedocount and convert to proportions
tmp <- humann.out.t %>% as.matrix()
tmp <- ifelse(tmp == 0, 0.5, tmp) %>% proportions(margin = 1)

humann.out.t <- tmp %>%
  as.data.frame() %>%
  dplyr::mutate(Sample_Name = rownames(.)) %>%
  dplyr::select(Sample_Name, everything())

print(humann.out.t)
#>          Sample_Name       pwy1        pwy2      pwy3        pwy4       pwy5
#> Sample1      Sample1 0.13114754 0.163934426 0.2459016 0.278688525 0.18032787
#> Sample2      Sample2 0.22784810 0.012658228 0.4050633 0.303797468 0.05063291
#> Sample3      Sample3 0.01960784 0.196078431 0.4705882 0.254901961 0.05882353
#> Sample4      Sample4 0.29850746 0.179104478 0.2686567 0.014925373 0.23880597
#> Sample5      Sample5 0.16666667 0.312500000 0.1875000 0.125000000 0.20833333
#> Sample6      Sample6 0.00800000 0.080000000 0.1600000 0.336000000 0.41600000
#> Sample7      Sample7 0.10884354 0.108843537 0.2040816 0.006802721 0.57142857
#> Sample8      Sample8 0.07619048 0.171428571 0.0952381 0.133333333 0.52380952
#> Sample9      Sample9 0.02547771 0.006369427 0.1656051 0.165605096 0.63694268
#> Sample10    Sample10 0.05617978 0.101123596 0.1123596 0.033707865 0.69662921
```

Now, we just need a vector indicating group membership and a matrix of
the covariate(s) of interest to perform the differential abundance
analysis. The file ‘meta_dat_demo.csv’ gives the simulated response
variable and a simulated covariate. This can be loaded in using the
following chunk of code.

``` r
meta.dat <- read.csv(file = "meta_dat_demo.csv", header = TRUE)
head(meta.dat)
#>   X Group  Covariate
#> 1 1     0 -1.2432097
#> 2 2     0  0.9918590
#> 3 3     0  0.8639651
#> 4 4     0 -1.7550385
#> 5 5     0  0.1488082
#> 6 6     1  1.1100746
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
res <- two.stage.da.analysis(y = y, fun_prof = humann.out.t, tax_prof = metaphlan.out, Fjs = Fjs, x = X)
print(res)
#> [1] "Error: Functional Profile Must Be a Matrix."
```

The output above tells us that pathway 5 is the only differentially
abundant pathway, and only taxon 4 is differentially abundant.
