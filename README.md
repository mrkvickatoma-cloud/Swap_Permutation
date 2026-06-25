# Paired Swap Permutation Test

An exact, non-parametric methodological framework for comparing the relative explanatory power of two dependent predictors ($X_1$, $X_2$) regarding a common target variable ($Y$).

## Overview

Comparing the predictive strength of two dependent variables is a fundamental statistical challenge. Classical asymptotic procedures (e.g., Vuong's closeness test, Hotelling-Williams test) frequently collapse under pathological data conditions such as heavy-tailed distributions or extreme categorical sparsity. Traditional resampling methods also fail:
* **Naive permutations** independently shuffle predictors, destroying their covariance structure and inflating Type I error.
* **Paired bootstrap** evaluates variance around the alternative hypothesis ($H_1$) rather than the null ($H_0$).

The **Paired Swap Permutation Test** overcomes these limitations by evaluating the test statistic strictly under the exact null hypothesis of functional exchangeability.

## Key Features

* **Categorical Domains (The Symmetric Swap):** For discrete data sharing the same state space, the algorithm uses a symmetric within-subject swapping mechanism evaluated via **Mutual Information**. It is immune to sparsity-induced entropy compression.
* **Continuous Domains (The ECDF-Mapped Swap):** For continuous data with differing marginal distributions, the algorithm introduces an Empirical Cumulative Distribution Function (ECDF) mapping step. It exchanges relative ranks rather than raw values, preventing marginal density corruption and avoiding metric space compression.
* **Supported Statistics:** Mutual Information (discrete), Pearson's $r$, Kendall's $\tau$, and Distance Covariance (`dcov`).
* **Auto-Standardization:** Safely standardizes continuous variables to prevent scale dominance (critical for Distance Covariance).

## Installation & Requirements

The function relies on base R for most operations. If you wish to use Distance Covariance, the `energy` package is required. No complex installation is needed—simply source the script in your R environment.

```R
# Install 'energy' if using distance covariance
# install.packages("energy")

# Source the main function
source("paired_swap_test.R")