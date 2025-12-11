#' Run The Two-Stagewise Method for Differential Abundance Analysis
#'
#' This is a function that takes a data frame containing the abundances of the functions along with their associated taxa and
#' a data frame containing the outcome variable and any covariates of interest. It runs the two-stagewise approach to finding
#' BRp taxa and returns the results from both stages along with a list of all identified BRp taxa.
#'
#' @param wmgs.dat A data frame containing information on the functions and taxa associated with each function for all samples.
#' The first column should be labeled X..Pathway and contain the names of the functions and taxa linked with each function. The
#' remaining columns should contain the sample IDs and list the counts for each function and linked taxon.
#' @param meta.dat A data frame containing information on the outcome variable and any covariates for each subject. This must have one
#' column labeled "SampleID" to denote the participants and another column labeled "y" to denote the outcome variable. The remaining
#' columns should include information on any covariates of interest.
#' @param f.cutoff Adjusted p-value threshold to use for the first stage to identify DA functions.
#' @param tax.cutoff Adjusted p-value threshold to use for the second stage to identify BRp taxa.
#'
#' @return This function returns a list containing the following:
#'
#' fun.res: A data frame containing the results from running the first stage of the method where we identify DA functions.
#' tax.res: A list containing the results from the second stage where we try to identify BRp taxa given each function. The names
#' of each element of the list correspond to the functions, and the data frame contains the results from the second stage.
#' brp.tax: A vector containing the names of all identified BRp taxa.
#'
get.brp.taxa <- function (wmgs.dat, meta.dat, f.cutoff = 0.1, tax.cutoff = 0.1) {
  ### wmgs.dat: Data frame containing information on the functions and the taxa linked to them
  ### meta.dat: Data frame containing data on response and covariate(s)

  ############################## First Stage: DA Functions ##############################
  ### Data preprocessing on functional profile
  F.prof <- wmgs.dat %>%
    dplyr::filter(!grepl(pattern = "|", x = wmgs.dat[,"X..Pathway"], fixed = TRUE)) %>%
    column_to_rownames(var = "X..Pathway") %>%
    as.matrix() %>%
    t()
  pc.f <- min(F.prof[which(F.prof > 0)]) / 10 # pseudo count for functional profile
  F.prof <- ifelse(F.prof == 0, pc.f, F.prof) %>% proportions(margin = 1) # convert functional profile to proportions

  ### Complete data needed for first stage
  full.dat <- F.prof %>%
    as.data.frame() %>%
    rownames_to_column(var = "SampleID") %>%
    dplyr::inner_join(meta.data, by = "SampleID")
  full.dat %>% head()

  ### Extract outcome variable, functional profile, and covariate(s)
  y <- full.dat[,"y"]
  fun.dat <- full.dat[,colnames(F.prof)] %>% as.matrix()
  x <- full.dat %>% dplyr::select(-SampleID, -y, -colnames(fun.dat)) %>% as.matrix()

  ### Run first stage of our method
  rslt.s1 <- comp.debias.probit(y = y, M = fun.dat, x = x)
  rownames(rslt.s1) <- c("Intercept", colnames(fun.dat), colnames(x))

  ### Store names of DA pathways in a vector
  da.pwys <- rslt.s1[which(rslt.s1[,"adjp"] < f.cutoff),] %>% rownames() # DA pathway names



  ############################## Second Stage: Taxa Associated with Significant Functions ##############################
  brp.res <- vector(mode = "list", length = length(da.pwys)) # initialize list to store results from second stage
  names(brp.res) <- da.pwys
  brp.tax.names <- NULL # initialize vector to store names of BRp taxa

  for (da.pwy in da.pwys) {
    ### Data preprocessing on taxonomic profile
    tax.prof <- wmgs.dat %>%
      dplyr::filter(grepl(pattern = da.pwy, x = wmgs.dat[,"X..Pathway"], fixed = TRUE)) %>%
      column_to_rownames(var = "X..Pathway") %>%
      as.matrix() %>%
      t()
    tax.prof[,-1] <- ifelse(tax.prof[,-1] == 0, 0.5, tax.prof[,-1]) %>% proportions(margin = 1) # add pseudo count and convert to proportions
    tax.prof[,da.pwy] <- log(F.prof[,da.pwy]) # use log of function abundance as outcome
    head(tax.prof)

    ### Complete data needed for second stage
    full.dat2 <- tax.prof %>%
      as.data.frame() %>%
      rownames_to_column(var = "SampleID") %>%
      dplyr::inner_join(meta.data, by = "SampleID")
    full.dat2 %>% head()

    ### Extract outcome variable, taxonomic profile, and covariate(s)
    y <- full.dat2[,da.pwy]
    tax.dat <- full.dat2 %>% dplyr::select(contains("|")) %>% as.matrix()
    x <- full.dat2 %>% dplyr::select(-SampleID, -y, -colnames(tax.prof)) %>% as.matrix()

    ### Run second stage of the proposed method
    rslt.s2 <- comp.debias(y = y, M = tax.dat, x = x)
    rownames(rslt.s2) <- c(colnames(tax.dat), colnames(x))
    brp.res[[da.pwy]] <- rslt.s2

    ### Store names of identified BRp taxa in a vector
    brp.tax.names <- c(brp.tax.names, rownames(rslt.s2[which(rslt.s2$adjp < tax.cutoff),]))
    brp.tax.names <- brp.tax.names[grepl(pattern = "|", x = brp.tax.names, fixed = TRUE)]
  }

  return(list("fun.res" = rslt.s1,
              "tax.res" = brp.res,
              "brp.tax" = brp.tax.names))
}




