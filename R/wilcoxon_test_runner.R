###############################################################################
### 1. Set env
###############################################################################

library(tidyverse)
library(vegan)
library(doParallel)

###############################################################################
### 2. Def wilcox subfunction
###############################################################################

wilcoxon_test_runner_subfunction <- function(X_g, X_r, P, j) {
  
  X_g_j <- X_g %>% dplyr::select(j) %>% .[,1]
  X_r_j <- X_r %>% dplyr::select(j) %>% .[,1]

  X_g_j_sum <- sum(X_g_j)
  X_r_j_sum <- sum(X_r_j)

  if (X_g_j_sum != 0 || X_r_j_sum != 0) {
    x <- wilcox.test(X_g_j, X_r_j, alternative = "two.sided")
    
    if (x$p.value <= P) {
      return(j)
    } 
  } 
}

###############################################################################
### 3. Def fun wilcox for loop
###############################################################################

#' Wilcoxon rank sum test runner
#'
#' This function tests the differential abundance of features (e.g., ASVs, OTUs, or OPUs) between groups using the Wilcoxon rank sum test. The function includes the rarification of the abundance matrix and the filtering of singletons.
#' @param X Abundance table (i.e., data.frame or matrix) formatted as samples x features. Sample names should be set as rownames.
#' @param G Groups table (i.e., data.frame). This table should have two columns: the first column with the sample names and the second column with groups numbers. 
#' @param P P-value used to select singificant features. Default 1e-3.
#' @param NSLOTS Number of cores used. Default 4. 
#' @return This function returns a list where each element corresponds to a group (tested against all the other groups) and consists of a vector of significant features. 
#' @export

wilcoxon_test_runner <- function(X, G, P = 1e-3, NSLOTS = 4){
  
  registerDoParallel(cores=NSLOTS)
  print(paste("Using", NSLOTS, "cores"))
  
  # rarefy
  rare_min <- rowSums(X) %>% min()
  X_rare <- rrarefy(x = X, sample = rare_min)
  print(c("min count to rarefy:", rare_min))
  
  # remove singletones
  i <- colSums(X_rare) > 1
  X_rare_filt <- X_rare[,i] %>%
                 as.data.frame() %>%
                 rownames_to_column("sample")
  print(c("num of singletons:", sum(!i)))
  print(c("dim (rarefied and no singletons):", dim(X_rare_filt)))
  
  # join tables
  colnames(G) <- c("sample","group")
  X_rare_file_ext <- left_join(X_rare_filt, G, by = "sample")
  
  # Get dimensions
  ASVS <- X_rare_file_ext %>% 
          select(starts_with("ASV")) %>% 
          colnames()
  
  # Define groups
  i <- unique(G$group) %>% order()
  GROUPS <- unique(G$group)[i]
  
  # run wilcox test
  SELECTED_ASVS <- list()
  
  for (i in GROUPS) {
    
    X_g <- X_rare_file_ext %>%
           filter(group == i)
    
    X_r <- X_rare_file_ext %>%
            filter(group %in% GROUPS[(! GROUPS %in% i )])
    
  x <- foreach(j = ASVS) %dopar% 
       wilcoxon_test_runner_subfunction(X_g, X_r, P, j)
  
  SELECTED_ASVS[[i]]<- unlist(x)
  
  }
  return(SELECTED_ASVS)
}
