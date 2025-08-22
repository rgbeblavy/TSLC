
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TSLC

<!-- badges: start -->
<!-- badges: end -->

The goal of TSLC is to take into account information on both the taxa and genes/gene pathways in order to find differentially abundant taxa that are associated with a potential difference in function across two populations of interest. This is a two-stage approach that looks at differences in function across the population in the first stage, and then looks at differences in the taxa associated with each differentially abundant function in the second stage.

## Installation

You can install the development version of TSLC like so:

``` r
install_github(repo = "rgbeblavy/TSLC")
```

## Example

Here, we will illustrate how to use this package to implement the two-stagewise approach to find differentially abundant taxa given differentially abundant functions. For the ease of illustration, we use a very low-dimensional data set consisting of 10 subjects, 10 taxa, and 5 gene pathways. The first 5 subjects belong to one population, while the last 5 subjects belong to another population. The first 4 pathways are not differentially abundant across the two populations, while the last pathway is differentially abundant across the two populations. The taxa are labeled A through J, and only species J is truly differentially abundant. In practice, this method is designed to be applied to a high-dimensional taxonomic and functional profile.

First, we will load in the file 'tax_prof_demo.csv', which gives synthetic data designed to mimic the output from Metaphlan. The following chunk of code loads this in. Note that the rows correspond to sample IDs and the columns correspond to taxa.

```{r example}
library(tidyverse)
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
```


Now, we will load in the file 'functional_prof_demo.csv', which contains synthetic data designed to mimic the output from Humann. This will give us information about the abundances of the pathways and the taxa associated with each pathway. The following chunk of code does this.

```{r}
humann.out <- read.csv(file = "functional_prof_demo.csv", header = TRUE)
humann.out[,1] <- NULL # Remove the index column
head(humann.out)
```

In the data frame above, the rows that do not have a vertical bar "|" give the aggregated abundances of each of the synthetic pathways, while the rows with a vertical bar give the abundance of the pathway associated with a given species. With the functional profile in this form, we can apply the function 'get.taxa' to find the list of taxa associated with each pathway. The following chunk of code does this.

```{r}
Fjs <- get.taxa(humann.out)
print(Fjs)
```

We will need the above list for running the function to implement the two-stagewise approach to differential abundance analysis. Before doing this, we will need to get the functional profile in the proper form where the samples are rows and the columns are the pathways. The following chunk of code will help us convert the Humann output into the correct form.

```{r}
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
```


Now, we just need a vector indicating group membership and a matrix of the covariate(s) of interest to perform the differential abundance analysis. The file 'meta_dat_demo.csv' gives the simulated response variable and a simulated covariate. This can be loaded in using the following chunk of code.

```{r}
meta.dat <- read.csv(file = "meta_dat_demo.csv", header = TRUE)
head(meta.dat)
```

We can then use this to extract the vector indicating group membership and matrix containing the covariate.

```{r}
y <- meta.dat[,"Group"] %>% as.vector()
X <- meta.dat[,"Covariate"] %>% as.matrix()
```

Finally, we can use these quantities to determine differentially abundant taxa given differentially abundant functions.

```{r}
res <- two.stage.da.analysis(y = y, fun_prof = humann.out.t, tax_prof = metaphlan.out, Fjs = Fjs, x = X)
print(res)
```

The output above tells us that pathway 5 is the only differentially abundant pathway, and only taxon 4 is differentially abundant.






