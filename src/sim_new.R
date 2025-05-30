# Rashomon Effect Simulation
# Based on Breiman (2001): "Statistical Modeling: The Two Cultures"

library(leaps)  # for regsubsets (best subset selection)
library(data.table)
library(MASS)

set.seed(123)  # for reproducibility

# Generate population data
n_pop <- 10000
n_vars <- 30

# Generate correlated covariates (this creates more similar models)
# Create correlation matrix with some structure
Sigma <- matrix(0.3, n_vars, n_vars)  # baseline correlation
diag(Sigma) <- 1

# Add some stronger correlations between groups of variables
# Group 1: variables 1-5
Sigma[1:5, 1:5] <- 0.7
diag(Sigma) <- 1

# Group 2: variables 6-10
Sigma[6:10, 6:10] <- 0.6
diag(Sigma) <- 1

# Group 3: variables 11-15
Sigma[11:15, 11:15] <- 0.65
diag(Sigma) <- 1

# Generate correlated predictors
X_pop <- mvrnorm(n_pop, mu = rep(0, n_vars), Sigma = Sigma)
colnames(X_pop) <- paste0("X", 1:n_vars)

# Generate coefficients - many small effects rather than few large ones
# This creates a situation where many subsets can perform similarly
true_coefs <- rnorm(n_vars, mean = 0, sd = 1)
# Make some coefficients larger
important_vars <- sample(1:n_vars, 10)
true_coefs[important_vars] <- true_coefs[important_vars] * 2

# Generate outcome with moderate noise
# Higher noise relative to signal creates more model uncertainty
signal <- X_pop %*% true_coefs
noise_sd <- sd(signal) * 0.8  # noise is 80% of signal SD
y_pop <- signal + rnorm(n_pop, sd = noise_sd)

# Combine into data frame
pop_data <- data.frame(y = y_pop, X_pop)

# Sample sizes to test
sample_sizes <- seq(from=100, to=2000, by=100)
sample_sizes = c(100)

# Function to find models within threshold of best model
find_rashomon_set <- function(data, n_vars_select = 5, threshold = 0.01) {
  # Get all possible 5-variable combinations
  all_vars <- names(data)[-1]  # exclude 'y'
  all_combos <- combn(all_vars, n_vars_select)
  
  # Calculate RSS for each combination
  n_models <- ncol(all_combos)
  rss_values <- numeric(n_models)
  
  cat(sprintf("  Evaluating %d possible 5-variable models...\n", n_models))
  
  # Progress indicator for large computations
  pb <- txtProgressBar(min = 0, max = n_models, style = 3)
  
  for (i in 1:n_models) {
    vars <- all_combos[, i]
    formula <- as.formula(paste("y ~", paste(vars, collapse = " + ")))
    lm_fit <- lm(formula, data = data)
    rss_values[i] <- sum(residuals(lm_fit)^2)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Find minimum RSS
  min_rss <- min(rss_values)
  
  # Count models within threshold
  threshold_rss <- min_rss * (1 + threshold)
  n_models_in_set <- sum(rss_values <= threshold_rss)
  
  # Also calculate for different thresholds to see distribution
  thresholds <- c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50)
  counts <- sapply(thresholds, function(t) {
    sum(rss_values <= min_rss * (1 + t))
  })
  
  cat(sprintf("  Models within various thresholds of best:\n"))
  for (i in seq_along(thresholds)) {
    cat(sprintf("    Within %d%%: %d models\n", thresholds[i]*100, counts[i]))
  }
  
  return(list(
    n_models_1pct = n_models_in_set,
    all_counts = data.frame(
      threshold = thresholds,
      n_models = counts
    ),
    min_rss = min_rss,
    rss_values = rss_values
  ))
}

# Run simulation for each sample size
results <- data.table(sample_size = integer(), 
                      n_models_rashomon = integer())

# Store detailed results for analysis
detailed_results <- list()

cat("Running Rashomon Effect simulation...\n")
cat("Data generated with correlated predictors and distributed effects\n\n")

for (n in sample_sizes) {
  cat(sprintf("Sample size: %d\n", n))
  
  # Take sample from population
  if (n <= n_pop) {
    sample_idx <- sample(1:n_pop, n, replace = FALSE)
  } else {
    # For samples larger than population, sample with replacement
    sample_idx <- sample(1:n_pop, n, replace = TRUE)
  }
  
  sample_data <- pop_data[sample_idx, ]
  
  # Find Rashomon set size
  rashomon_results <- find_rashomon_set(sample_data)
  
  # Store results
  results <- rbind(results, 
                   data.table(sample_size = n, 
                              n_models_rashomon = rashomon_results$n_models_1pct))
  
  # Store detailed results
  detailed_results[[as.character(n)]] <- rashomon_results
  
  cat("\n")
}

# Save main results
write.csv(results, file = '../data/result.csv', row.names = FALSE)

# Save detailed results for different thresholds
all_threshold_results <- data.table()
for (n in names(detailed_results)) {
  temp <- detailed_results[[n]]$all_counts
  temp$sample_size <- as.integer(n)
  all_threshold_results <- rbind(all_threshold_results, temp)
}
write.csv(all_threshold_results, file = '../data/detailed_results.csv', row.names = FALSE)

cat("\nResults saved to '../data/result.csv' and '../data/detailed_results.csv'\n")
print(results)

# Show correlation structure impact
cat("\n\nCorrelation structure of predictors:\n")
cat("Average correlation: ", mean(Sigma[upper.tri(Sigma)]), "\n")
cat("This creates many competing models with similar performance\n")

