#' Paired Swap Permutation Test
#'
#' Evaluates whether two dependent predictors (X1, X2) provide equal 
#' explanatory power regarding a target variable (Y).
#'
#' @param Y A vector containing the target variable.
#' @param X1 A vector containing the first predictor.
#' @param X2 A vector containing the second predictor.
#' @param statistic Character string specifying the test statistic. 
#'        Options: "mutual_information" (discrete), "pearson" (continuous), 
#'        "kendall" (continuous), "dcov" (continuous).
#' @param B Integer specifying the number of permutation iterations. Default is 10000.
#' @param standardize Logical. If TRUE, continuous variables are automatically 
#'        standardized (scaled to zero mean and unit variance) before test evaluation. 
#'        Highly recommended for Distance Covariance (dcov). Ignored for discrete metrics.
#' @return A list with class "paired_swap" containing the test results.
#' @export
paired_swap_test <- function(Y, X1, X2, 
                             statistic = c("mutual_information", "pearson", "kendall", "dcov"), 
                             B = 10000,
                             standardize = TRUE) {
  
  # 1. Match argument and validate inputs
  statistic <- match.arg(statistic)
  n <- length(Y)
  
  if (length(X1) != n || length(X2) != n) {
    stop("Vectors Y, X1, and X2 must all have the same length.")
  }
  
  # Require energy package if dcov is selected
  if (statistic == "dcov" && !requireNamespace("energy", quietly = TRUE)) {
    stop("Package 'energy' is required for distance covariance. Please install it.")
  }
  
  # 2. Standardize continuous variables if requested
  if (statistic != "mutual_information" && standardize) {
    Y  <- as.vector(scale(Y))
    X1 <- as.vector(scale(X1))
    X2 <- as.vector(scale(X2))
  }
  
  # 3. Define the internal function for calculating the chosen statistic
  calc_stat <- function(X, Y_t, stat) {
    if (stat == "pearson") return(cor(X, Y_t, method = "pearson"))
    if (stat == "kendall") return(cor(X, Y_t, method = "kendall"))
    if (stat == "dcov")    return(energy::dcov(X, Y_t))
    if (stat == "mutual_information") {
      # Calculate empirical Mutual Information I(X;Y)
      tbl <- table(Y_t, X)
      p_xy <- tbl / n
      p_x <- colSums(tbl) / n
      p_y <- rowSums(tbl) / n
      p_x_p_y <- outer(p_y, p_x, "*")
      
      # Sum only over non-zero probabilities to avoid log2(0)
      valid <- p_xy > 0
      mi <- sum(p_xy[valid] * log2(p_xy[valid] / p_x_p_y[valid]))
      return(mi)
    }
  }
  
  # 4. Calculate observed statistics
  stat1_obs <- calc_stat(X1, Y, statistic)
  stat2_obs <- calc_stat(X2, Y, statistic)
  diff_obs <- stat1_obs - stat2_obs
  
  diff_sim <- numeric(B)
  
  # 5. Perform the Permutation Loop based on the Domain
  if (statistic == "mutual_information") {
    # DISCRETE DOMAIN: Symmetric Direct Swap
    for (i in seq_len(B)) {
      swap_mask <- sample(c(TRUE, FALSE), n, replace = TRUE)
      
      X1_perm <- X1
      X2_perm <- X2
      
      # Swap values directly
      X1_perm[swap_mask] <- X2[swap_mask]
      X2_perm[swap_mask] <- X1[swap_mask]
      
      diff_sim[i] <- calc_stat(X1_perm, Y, statistic) - calc_stat(X2_perm, Y, statistic)
    }
    method_name <- "Paired Swap Permutation Test (Categorical/Symmetric)"
    
  } else {
    # CONTINUOUS DOMAIN: ECDF-Mapped Swap
    X1_sorted <- sort(X1)
    X2_sorted <- sort(X2)
    rank_X1 <- rank(X1, ties.method = "random")
    rank_X2 <- rank(X2, ties.method = "random")
    
    for (i in seq_len(B)) {
      swap_mask <- sample(c(TRUE, FALSE), n, replace = TRUE)
      
      X1_perm <- X1
      X2_perm <- X2
      
      # Swap relative ranks via ECDF mapping
      X1_perm[swap_mask] <- X1_sorted[rank_X2[swap_mask]]
      X2_perm[swap_mask] <- X2_sorted[rank_X1[swap_mask]]
      
      diff_sim[i] <- calc_stat(X1_perm, Y, statistic) - calc_stat(X2_perm, Y, statistic)
    }
    method_name <- sprintf("ECDF-Mapped Paired Swap Test (%s)", tools::toTitleCase(statistic))
  }
  
  # 6. Compute Exact P-value (accounting for the observed data as one permutation)
  p_value <- (sum(abs(diff_sim) >= abs(diff_obs)) + 1) / (B + 1)
  
  # 7. Construct and return result object
  res <- list(
    method = method_name,
    statistic_name = statistic,
    B = B,
    standardized = (statistic != "mutual_information" && standardize),
    stat1_obs = stat1_obs,
    stat2_obs = stat2_obs,
    diff_obs = diff_obs,
    p_value = p_value,
    null_distribution = diff_sim
  )
  class(res) <- "paired_swap"
  return(res)
}

#' Print method for paired_swap objects
#' @export
print.paired_swap <- function(x, ...) {
  cat("\n\t", x$method, "\n\n")
  cat("Data: Evaluated", x$B, "permutations")
  if (x$standardized) cat(" (Continuous variables were automatically standardized)")
  cat("\n")
  
  cat(sprintf("Observed %s (X1): %.4f\n", x$statistic_name, x$stat1_obs))
  cat(sprintf("Observed %s (X2): %.4f\n", x$statistic_name, x$stat2_obs))
  cat(sprintf("Observed Difference (X1 - X2): %.4f\n", x$diff_obs))
  cat(sprintf("Exact p-value = %s\n\n", format.pval(x$p_value, digits = 4, eps = 1/(x$B + 1))))
  
  if(x$p_value < 0.05) {
    if(x$diff_obs > 0) cat("Alternative hypothesis: X1 has significantly stronger explanatory power than X2.\n")
    else cat("Alternative hypothesis: X2 has significantly stronger explanatory power than X1.\n")
  } else {
    cat("Null hypothesis: X1 and X2 have equal explanatory power.\n")
  }
}