#' Obtain Taxa Associated With Each Pathway
#'
#' This is a function that takes the output from Humann as an input and returns a list of all
#' the species associated with each of the pathways.
#'
#' @param func A data frame where the first column contains the pathways and associated taxa and the remaining
#' columns contain the sample IDs. The first column must have the name X..Pathway.
#'
#' @return This function returns a list where the names are the gene pathways and each element of the list
#' is a vector containing all species associated with the corresponding pathway.
get.taxa <- function (func) {
  ### Here func is the Humann output

  ### Step 1: Get vector with just function names func
  func.t <- func %>% dplyr::filter(!grepl("|", X..Pathway, fixed = TRUE))
  func_names <- func.t$X..Pathway

  ### Step 2: Get vector with func|g__.s__
  all_func_names <- func$X..Pathway
  taxa_names <- all_func_names[grepl("|", all_func_names, fixed = TRUE)]

  ### Step 3:
  species_df <- data.frame("Species_Name" = 0,
                           "X..Pathway" = taxa_names)

  for (i in 1:length(taxa_names)) {

    if (grepl("s__", taxa_names[i], fixed = TRUE) == TRUE) {
      #species_df[i, 1] <- strsplit(taxa_names[i], split = ".s__", fixed = TRUE)[[1]][2]
      species_df[i, 1] <- strsplit(taxa_names[i], split = "s__", fixed = TRUE)[[1]][2]
    }

    if (grepl("|unclassified", taxa_names[i], fixed = TRUE) == TRUE) {
      species_df[i, 1] <- "unclassified"
    }

  }

  species_df[,"Species_Name"] <- paste0("s__", species_df[,"Species_Name"])

  ### Step 4: Construct the list Fjs
  Fjs <- vector(mode = "list", length = length(func_names))
  names(Fjs) <- func_names

  for (i in 1:length(func_names)) {

    if (grepl(pattern = "[", x = func_names[i], fixed = TRUE)) {
      func_names[i] <- gsub(pattern = "[", replacement = "", x = func_names[i], fixed = TRUE)
    }

    if (grepl(pattern = "]", x = func_names[i], fixed = TRUE)) {
      func_names[i] <- gsub(pattern = "]", replacement = "", x = func_names[i], fixed = TRUE)
    }

    Fjs[[i]] <- species_df[which(grepl(pattern = func_names[i], x = species_df[,2]) == TRUE), 1]
  }

  return(Fjs)
}

#' Run The Two-Stagewise Method for Differential Abundance Analysis
#'
#' This is a function that takes a vector indicating the group membership, the functional profile,
#' the taxonomic profile, a list containing all taxa associated with each function, and a matrix of covariates
#' and runs the two-stagewise approach to finding differentially abundant taxa given differentially
#' abundant functions.
#'
#' @param y A vector with elements from {0, 1} to indicate which population each sample belongs to. The order
#' of the elements of y must match the order of the samples in the rows of the functional profile, taxonomic
#' profile, and matrix of covariates.
#' @param fun_prof A matrix with the samples as rows and pathways as columns. Note that the abundances of
#' each sample must be converted to proportions in order for this to work. Also, the order of the rows
#' must match the order of the samples in the response vector.
#' @param tax_prof A matrix with the samples as rows and taxa as columns. Note that the abundances of each sample
#' must be converted to proportions in order for this to work. The rownames of the taxonomic profile must
#' match up to the rownames of the functional profile.
#' @param Fjs A list where the names of the list are the same as the column names of the functional profile,
#' and each element of the list is a vector indicating which microbes in the taxonomic profile are associated
#' with the corresponding function. The elements of the vector can be given as either column indices or names.
#' @param x A matrix where the samples are rows and the columns are covariates to be adjusted for in the model.
#' The rownames must match up to those of the functional and taxonomic profile, and any NA values must be
#' removed from this matrix before running the method.
#' @param cutoff A number to indicate which threshold the adjusted p-value should be cutoff at. All taxa with
#' an adjusted p-value below this threshold will be detected as differentially abundant. To correct for
#' multiple comparisons, we use the Benjamini-Hochberg method.
#'
#' @return This function returns a list with the following elements.
#'
#' da.taxa: A list where the names are each function identified as differentially abundant and the elements
#' of the list are vectors containing all associated differentially abundant taxa.
#'
#' result: A data frame where the columns are each taxa name, their nominal p-value, and adjusted p-value.
two.stage.da.analysis <- function (y, fun_prof, tax_prof, Fjs, x, cutoff = 0.20) {
  num.fun <- ncol(fun_prof); num.taxa <- ncol(tax_prof)
  pvals <- rep(1, num.taxa); adjp <- rep(1, num.taxa)
  names(pvals) <- colnames(tax_prof); names(adjp) <- colnames(tax_prof)

  if (!is.matrix(fun_prof)) {
    return("Error: Functional Profile Must Be a Matrix.")
  }

  if (sum(rowSums(fun_prof)) != nrow(fun_prof)) {
    return("Error: Each row of functional profile must be a composition.")
  }

  if (!is.matrix(tax_prof)) {
    return("Error: Taxonomic Profile Must Be a Matrix.")
  }

  if (sum(rowSums(tax_prof)) != nrow(tax_prof)) {
    return("Error: Each row of taxonomic profile must be a composition.")
  }

  if (!is.null(x)) {
    if (!is.matrix(x)) {
      return("Error: Covariates must be given as a matrix.")
    }
  }

  if (sum(sort(unique(y)) == c(0,1)) != 2) {
    return("Error: Response Y must be a binary variable with elements {0,1}.")
  }

  if (length(y) != nrow(fun_prof)) {
    return("Error: Number of samples in functional profile do not match length of response variable.")
  }

  if (length(y) != nrow(tax_prof)) {
    return("Error: Number of samples in taxonomic profile do not match length of response variable.")
  }

  if (nrow(fun_prof) != nrow(tax_prof)) {
    return("Error: Number of samples in functional profile do not match number of variables in taxonomic profile.")
  }

  ### First Stage: Identify DA Functions
  rslt.F <- comp.debias.probit(y = y, M = fun_prof, x = x)
  da.fun.nums <- which(rslt.F[,"adjp"] < cutoff) - 1 # The -1 is because of the intercept
  da.fun.nums <- da.fun.nums[da.fun.nums > 0]
  da.fun.nums <- da.fun.nums[which(da.fun.nums <= num.fun)]

  if (length(da.fun.nums) == 0) {
    return("No taxa are identified as differentially abundant.")
  }

  da.fun.names <- colnames(fun_prof)[da.fun.nums]

  ### Second Stage: Identify DA Taxa
  ### For each function identified as DA, we perform DA analysis on the taxa associated with it
  my.list <- vector(mode = "list", length = length(da.fun.names))
  names(my.list) <- da.fun.names

  for (fun.name in da.fun.names) {
    taxa.associated <- Fjs[[fun.name]]

    if (length(taxa.associated) < 2) {
      next
    }

    MM <- tax_prof[,taxa.associated]
    MM <- sweep(MM, 1, rowSums(MM), "/")

    rslt.M <- comp.debias.probit(y = y, MM, x = NULL) # Perform second stage

    rownames(rslt.M) <- c("Intercept", taxa.associated)
    my.list[[fun.name]] <- taxa.associated[which(rslt.M[2:nrow(rslt.M),"adjp"] < cutoff)]

    pvals[taxa.associated] <- rslt.M[1+(1:length(taxa.associated)),"pval"] # The 1+ is for the intercept
    adjp[taxa.associated] <- rslt.M[1+(1:length(taxa.associated)),"adjp"] # The 1+ is for the intercept
  }

  da.taxa.nums <- which(adjp < cutoff)
  da.taxa.names <- colnames(tax_prof)[da.taxa.nums]

  final.result.df <- data.frame("Taxon" = colnames(tax_prof),
                                "P.value" = pvals,
                                "Adj.p.value" = adjp)
  final.result.df <- final.result.df[-which(final.result.df[,"Adj.p.value"] == 1),] # Remove taxa not associated with a DA function

  return(list(da.taxa = my.list,
              result = final.result.df))
}


