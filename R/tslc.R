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

  ### Step 3: Construct data frame with species from each pathway name
  species_df <- data.frame("Species_Name" = 0,
                           "X..Pathway" = taxa_names)

  for (i in 1:length(taxa_names)) {

    if (grepl("s__", taxa_names[i], fixed = TRUE) == TRUE) {
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
#' the taxonomic profile, a list containing all taxa associated with each function, and a matrix of covariates.
#' It runs the two-stagewise approach to finding DA taxa given DA functions and returns a list of DA taxa
#' assocaiated with each DA function.
#'
#' @param y A vector with elements from {0, 1} to indicate which population each sample belongs to. The order
#' of the elements of y must match the order of the samples in the rows of the functional profile, taxonomic
#' profile, and matrix of covariates.
#' @param fun_prof A matrix with the data for the functional profile. The rows must be samples with a column
#' Sample_Name to indicate the sample IDs. The remaining columns must be for the pathways. Note that the
#' abundances of each sample must be converted to proportions in order for this to work. Also, the order of
#' the rows must match the order of the samples in the response vector.
#' @param tax_prof A matrix with the data for the taxonomic profile. The rows must be samples with a column
#' Sample_Name to indicate the sample IDs. The remaining columns must be for the pathways. Note that the
#' abundances of each sample must be converted to proportions in order for this to work. Also, the order of
#' the rows must match the order of the samples in the response vector.
#' @param Fjs A list where the names of the list are the same as the column names of the functional profile,
#' and each element of the list is a vector indicating which microbes in the taxonomic profile are associated
#' with the corresponding function. The elements of the vector can be given as either column indices or names.
#' @param x A matrix where the samples are rows and the columns are covariates to be adjusted for in the model.
#' The rownames must match up to those of the functional and taxonomic profile, and any NA values must be
#' removed from this matrix before running the method.
#' @param pwy_cutoff A number to indicate which threshold the p-value for the pathways should be cutoff at. All
#' pathways with a p-value below this threshold will be detected as DA.
#' @param tax_cutoff A number to indicate which threshold the p-value for the taxa should be cutoff at. All
#' taxa with a p-value below this threshold will be detected as DA. To correct for
#' multiple comparisons, we use the Benjamini-Hochberg method.
#'
#' @return This function returns the following:
#'
#' associated.da.taxa: A list where the names are each function identified as differentially abundant and the elements
#' of the list are vectors containing all associated differentially abundant taxa.
#'
two.stage.da.analysis <- function (y, fun_prof, tax_prof, Fjs, x, pwy_cutoff = 0.05, tax_cutoff = 0.10) {

  if (!is.null(x)) {
    if (!is.matrix(x)) {
      x <- as.matrix(x)
    }

    if (is.null(colnames(x))) {
      colnames(x) <- paste0("X", 1:ncol(x))
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


  ####################### First Stage: Identify DA Functions #######################
  F.prof <- fun_prof %>%
    dplyr::select(-Sample_Name) %>%
    as.matrix()

  if (sum(rowSums(F.prof)) != nrow(F.prof)) {
    return("Error: Each row of functional profile must be a composition.")
  }

  rslt <- comp.debias.probit(y, F.prof, x)
  rownames(rslt) <- c("Intercept", colnames(F.prof), colnames(x))

  ### Initialize a list to store DA species conditional on each DA pathway
  da.pwy <- rslt[which(rslt$pval < pwy_cutoff),] %>% rownames()

  if (length(da.pwy) == 0) {
    return("No functions are identified as differentially abundant.")
  }


  ####################### Second Stage: DA Taxa Given DA Functions #######################
  associated.da.species <- vector(mode = "list", length = length(da.pwy))
  names(associated.da.species) <- da.pwy

  for (pwy.name in da.pwy) {
    ### Get species names associated with pathway
    species.names <- Fjs[[pwy.name]]
    species.names <- species.names[which(species.names %in% colnames(tax_prof))] # keep only those species that met filtering criterion

    if (length(species.names) < 2) {
      associated.da.species[[pwy.name]] <- 0
      next
    }

    ### Create a data frame with function abundance and the abundance of associated taxa
    tmp2 <- tax_prof %>%
      dplyr::select(Sample_Name, species.names) %>%
      inner_join(fun_prof[,c("Sample_Name", pwy.name)], by = "Sample_Name")
    #tmp2[,pwy.name] <- log(tmp2[,pwy.name]) # take log of function abundance

    ### Run the second stage
    y.f <- tmp2[,pwy.name] %>% log()
    M <- tmp2[,species.names] %>%
      as.matrix() %>%
      proportions(margin = 1)

    rslt2 <- comp.debias(y = y.f, M = M)

    ### Add species associated with pathway to the list
    da.tax <- rslt2[which(rslt2$adjp < tax_cutoff),] %>% rownames()
    associated.da.species[[pwy.name]] <- da.tax

  }

  return(associated.da.species)
}

