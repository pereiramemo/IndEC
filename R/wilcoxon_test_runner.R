###############################################################################
### 1. Def wilcoxon subfunction
###############################################################################

# X <- X 
# G <- WQI_GROUPS_TBL 
# NSLOTS <- 6 
# P <- 5e-2

wilcoxon_test_runner_subfunction <- function(X_g, X_r, P, j) {
  
  X_g_j <- dplyr::select(X_g, all_of(j))[,1]
  X_r_j <- dplyr::select(X_r, all_of(j))[,1]

  X_g_j_sum <- sum(X_g_j)
  X_r_j_sum <- sum(X_r_j)

  if (X_g_j_sum != 0 || X_r_j_sum != 0) {
    
    x <- stats::wilcox.test(X_g_j, X_r_j, alternative = "two.sided")
    
    if (x$p.value <= P) {
      return(data.frame(asv = j, pval = x$p.value))
    } 
  } 
}

###############################################################################
### 3. Def fun wilcoxon for loop
###############################################################################

#' Wilcoxon rank sum test runner
#'
#' This function tests the differential abundance of features (e.g., ASVs, OTUs, or OPUs) between groups using the Wilcoxon rank sum test.
#' @param X Abundance table (i.e., data.frame) formatted as samples x features. Sample names should be set as rownames.
#' @param G Groups table (i.e., data.frame). This table should have two columns: the first column with the sample names and the second column with groups numbers. 
#' @param P P-value used to select singificant features. Default 1e-3.
#' @param NSLOTS Number of cores used. Default 4. 
#' @return This function returns a list where each element corresponds to a group (tested against all the other groups) and consists of a vector of significant features.
#' @importFrom foreach %dopar%
#' @importFrom dplyr %>%
#' @export

wilcoxon_test_runner <- function(X, G, P = 1e-3, NSLOTS = 4){
  
  doParallel::registerDoParallel(cores=NSLOTS)
  print(paste("Using", NSLOTS, "cores"))
  
  # format tables
  X <- tibble::rownames_to_column(X, "sample")
  X$sample <- as.character(X$sample)
  colnames(G) <- c("sample","group")
  G$sample <- as.character(G$sample)
  
  # join tables
  X_ext <- dplyr::left_join(X, G, by = "sample")
  
  # Get dimensions
  ASVS <- colnames(dplyr::select(X_ext, starts_with("ASV"))) 
          
  # Define groups
  i <- order(unique(G$group))
  GROUPS <- unique(G$group)[i]
  
  # run wilcoxon test
  SELECTED_ASVS <- list()
  
  for (i in GROUPS) {
    
    X_g <- dplyr::filter(X_ext, group == i)
    
    X_r <- dplyr::filter(X_ext, group %in% GROUPS[(! GROUPS %in% i )])
    
    x <- foreach::foreach(j = ASVS) %dopar% 
           wilcoxon_test_runner_subfunction(X_g, X_r, P, j)
  
    output_df <- do.call("rbind",x) %>%
                 as.data.frame() %>%
                 arrange(pval)
  
    SELECTED_ASVS[[i]] <- output_df
  
  }
  
  return(SELECTED_ASVS)
}
