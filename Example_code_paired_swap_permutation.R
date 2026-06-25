# ==============================================================================
# Paired Swap Permutation Test - Usage Examples
# ==============================================================================
# This script provides minimal reproducible examples for applying the 
# Paired Swap Permutation Test in both continuous and discrete domains.
# Ensure that the 'paired_swap_test.R' function is loaded before running.

source("paired_swap_test.R")

# ==============================================================================
# EXAMPLE 1: CONTINUOUS DOMAIN (Distance Covariance)
# ==============================================================================
cat("--------------------------------------------------\n")
cat("EXAMPLE 1: Continuous Domain (Distance Covariance)\n")
cat("--------------------------------------------------\n")

# Require MASS for multivariate normal generation
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
if (!requireNamespace("energy", quietly = TRUE)) install.packages("energy")

set.seed(23)

# 1. Simulate pathological continuous data (Heavy-tailed Y, Asymmetric X1/X2)
# Here, X1 is designed to be a significantly stronger predictor of Y than X2.
n_cont <- 50
rho <- 0.5
r1 <- 0.7  # Strong relation between X1 and latent Y
r2 <- 0.3  # Weak relation between X2 and latent Y

Sigma <- matrix(c(
  1,  r1,  r2,
  r1,  1, rho,
  r2, rho,  1
), nrow = 3)

Z <- MASS::mvrnorm(n_cont, mu = c(0, 0, 0), Sigma = Sigma)
U <- pnorm(Z)

# Map to asymmetric marginals
Y_cont  <- qlnorm(U[,1], meanlog = 0, sdlog = 1.5)  # Heavy right tail
X1_cont <- qnorm(U[,2], mean = 10, sd = 2)          # Normally distributed
X2_cont <- qlnorm(U[,3], meanlog = 2, sdlog = 0.5)  # Log-normally distributed

# 2. Run the ECDF-Mapped Paired Swap Test using Distance Covariance
cat("Running ECDF-Mapped Swap Test using 'dcov'...\n")
result_cont <- paired_swap_test(Y = Y_cont, X1 = X1_cont, X2 = X2_cont, 
                                statistic = "dcov", B = 5000)

# 3. View the results
print(result_cont)


# ==============================================================================
# EXAMPLE 2: CATEGORICAL DOMAIN (Mutual Information)
# ==============================================================================
cat("\n--------------------------------------------------\n")
cat("EXAMPLE 2: Categorical Domain (Mutual Information)\n")
cat("--------------------------------------------------\n")

set.seed(123)

# 1. Simulate categorical data (e.g., Linguistic compounds)
# Y is the semantic type. X1 and X2 share the same vocabulary (X_1 to X_10).
# X1 perfectly predicts the type, while X2 is heavily distorted by noise.
n_cat <- 100
y_levels <- c("Type_A", "Type_B", "Type_C")
x_levels <- paste0("X_", 1:10)

Y_cat <- character(n_cat)
X1_cat <- character(n_cat)
X2_cat <- character(n_cat)

for (i in 1:n_cat) {
  # Base latent root
  z <- sample(1:10, 1)
  
  # Y depends directly on the region of the root
  if (z <= 3) Y_cat[i] <- "Type_A"
  else if (z <= 7) Y_cat[i] <- "Type_B"
  else Y_cat[i] <- "Type_C"
  
  # X1 is a clean predictor (no noise)
  X1_cat[i] <- x_levels[z]
  
  # X2 is a weak predictor (heavy local noise applied)
  shift <- sample(c(-3, -2, 2, 3), 1)
  x2_idx <- max(1, min(10, z + shift))
  X2_cat[i] <- x_levels[x2_idx]
}

Y_cat <- factor(Y_cat, levels = y_levels)
X1_cat <- factor(X1_cat, levels = x_levels)
X2_cat <- factor(X2_cat, levels = x_levels)

# 2. Run the Symmetric Paired Swap Test using Mutual Information
cat("Running Symmetric Swap Test using 'mutual_information'...\n")
result_cat <- paired_swap_test(Y = Y_cat, X1 = X1_cat, X2 = X2_cat, 
                               statistic = "mutual_information", B = 5000)

# 3. View the results
print(result_cat)

#Optional: Plot the null distribution for the categorical test
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  df_plot <- data.frame(Diff = result_cat$null_distribution)
  ggplot(df_plot, aes(x = Diff)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = result_cat$diff_obs, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Null Distribution: Mutual Information Difference",
          x = "Simulated Difference I(X1;Y) - I(X2;Y)", y = "Frequency") +
     theme_minimal()
}
