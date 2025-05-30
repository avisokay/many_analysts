library(tidyverse)
library(combinat)

### This simulation demonstrates the Rashomon effect
### More sample size means fewer models in the rashomon set for the same epsilon
### We create one large dataset and sample from it at different sizes

# Function to create the master dataset
create_master_dataset <- function(n_obs = 30000, predictors = 30, seed = 42) {
  set.seed(seed)
  
  # Generate true coefficients - make them realistic but not too strong
  true_coeffs <- rep(0, predictors)
  # Make first 12 variables truly predictive with varying strengths
  true_coeffs[1:12] <- c(1.2, -0.8, 1.5, -0.6, 0.9, -1.1, 0.7, -0.4, 
                         1.0, -0.9, 0.5, -0.7)
  
  # Generate predictors with correlation structure
  X <- matrix(rnorm(n_obs * predictors), nrow = n_obs, ncol = predictors)
  
  # Add correlation structure to make variable selection more ambiguous
  correlation_strength <- 0.4
  for(i in 1:8) {
    for(j in 1:2) {
      if(i + j*8 <= predictors) {
        X[, i + j*8] <- correlation_strength * X[, i] + 
          sqrt(1 - correlation_strength^2) * X[, i + j*8]
      }
    }
  }
  
  # Generate response with fixed noise level
  noise_sd <- 2.5
  y <- X %*% true_coeffs + rnorm(n_obs, sd = noise_sd)
  
  cat("Created master dataset with", n_obs, "observations and", predictors, "predictors\n")
  cat("True signal variables: 1-12, with noise SD =", noise_sd, "\n\n")
  
  return(list(X = X, y = y, true_coeffs = true_coeffs))
}

# Function to evaluate Rashomon set for a given sample from master dataset
evaluate_rashomon_set <- function(master_data, sample_size, threshold, 
                                  n_var = 5, max_models = 15000, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Sample from master dataset
  n_total <- nrow(master_data$X)
  sample_indices <- sample(n_total, sample_size, replace = FALSE)
  
  X_sample <- master_data$X[sample_indices, , drop = FALSE]
  y_sample <- master_data$y[sample_indices]
  
  # Create all possible 5-variable combinations
  all_combinations <- combn(ncol(X_sample), n_var, simplify = FALSE)
  
  # Sample combinations if too many (but keep it large enough for good coverage)
  if (length(all_combinations) > max_models) {
    sampled_indices <- sample(length(all_combinations), max_models)
    all_combinations <- all_combinations[sampled_indices]
  }
  
  cat("Sample size:", sample_size, "- Evaluating", length(all_combinations), "models\n")
  
  # Evaluate each model
  results <- map_dfr(all_combinations, function(vars) {
    # Fit linear regression
    X_subset <- X_sample[, vars, drop = FALSE]
    
    # Use QR decomposition for numerical stability
    qr_fit <- qr(X_subset)
    coeffs <- qr.coef(qr_fit, y_sample)
    
    # Handle potential NA coefficients from rank deficiency
    if (any(is.na(coeffs))) {
      return(tibble(
        variables = list(vars),
        rss = Inf,
        coefficients = list(rep(NA, length(vars))),
        intercept = NA
      ))
    }
    
    # Calculate intercept and RSS
    intercept <- mean(y_sample) - mean(X_subset %*% coeffs)
    y_pred <- intercept + X_subset %*% coeffs
    rss <- sum((y_sample - y_pred)^2)
    
    tibble(
      variables = list(vars),
      rss = rss,
      coefficients = list(as.numeric(coeffs)),
      intercept = intercept
    )
  })
  
  # Remove any infinite RSS models
  results <- results[is.finite(results$rss), ]
  
  if (nrow(results) == 0) {
    return(list(
      sample_size = sample_size,
      threshold = threshold,
      total_models = 0,
      rashomon_set_size = 0,
      rashomon_proportion = 0,
      best_rss = NA,
      threshold_rss = NA
    ))
  }
  
  # Find best RSS and threshold
  best_rss <- min(results$rss)
  threshold_rss <- best_rss * (1 + threshold)
  rashomon_set <- results[results$rss <= threshold_rss, ]
  
  return(list(
    sample_size = sample_size,
    threshold = threshold,
    total_models = nrow(results),
    rashomon_set_size = nrow(rashomon_set),
    rashomon_proportion = nrow(rashomon_set) / nrow(results),
    best_rss = best_rss,
    threshold_rss = threshold_rss,
    rashomon_models = rashomon_set
  ))
}

# Main simulation function
run_rashomon_simulation <- function() {
  # Create master dataset once
  master_data <- create_master_dataset(n_obs = 1e5, predictors = 30, seed = 42)
  
  # Sample sizes to test - geometric progression
  sample_sizes <- c(1000, 5000, 10000, 25000, 50000)
  
  # Different thresholds to test
  thresholds <- c(0.01)
  
  # Run simulation across all combinations
  simulation_results <- expand_grid(
    sample_size = sample_sizes,
    threshold = thresholds
  ) %>%
    mutate(
      # Use different seeds for each sample to avoid bias
      seed = row_number() + 100
    ) %>%
    mutate(
      results = pmap(list(sample_size, threshold, seed), function(ss, th, sd) {
        evaluate_rashomon_set(master_data, ss, th, seed = sd, max_models = 15000)
      })
    ) %>%
    mutate(
      rashomon_size = map_dbl(results, ~.x$rashomon_set_size),
      rashomon_proportion = map_dbl(results, ~.x$rashomon_proportion),
      best_rss = map_dbl(results, ~.x$best_rss),
      total_models = map_dbl(results, ~.x$total_models)
    )
  
  return(simulation_results)
}

# Function to extract different "stories" from Rashomon set
extract_different_stories <- function(rashomon_models, n_stories = 3) {
  if (nrow(rashomon_models) < n_stories) {
    selected_indices <- 1:nrow(rashomon_models)
  } else {
    selected_indices <- sample(nrow(rashomon_models), n_stories)
  }
  
  stories <- map(selected_indices, function(idx) {
    vars <- rashomon_models$variables[[idx]]
    coeffs <- rashomon_models$coefficients[[idx]]
    intercept <- rashomon_models$intercept[[idx]]
    rss <- rashomon_models$rss[[idx]]
    
    # Create equation string
    terms <- paste0(round(coeffs, 2), "*x", vars)
    equation <- paste("y =", round(intercept, 2), "+", paste(terms, collapse = " + "))
    equation <- str_replace_all(equation, "\\+ -", "- ")
    
    list(
      model_index = idx,
      variables = vars,
      equation = equation,
      rss = rss
    )
  })
  
  return(stories)
}

# Function to analyze how variable selection changes with sample size
analyze_variable_stability <- function(results_data) {
  # Look at variable selection patterns across sample sizes for one threshold
  threshold_data <- results_data %>%
    filter(threshold == 0.02) %>%  # Use 2% threshold
    arrange(sample_size)
  
  # Extract variable importance for each sample size
  var_analysis <- map_dfr(1:nrow(threshold_data), function(i) {
    sample_info <- threshold_data[i, ]
    rashomon_models <- sample_info$results[[1]]$rashomon_models
    
    if (nrow(rashomon_models) == 0) {
      return(tibble(sample_size = sample_info$sample_size))
    }
    
    # Count variable frequency
    var_counts <- rep(0, 30)
    for(j in 1:nrow(rashomon_models)) {
      vars <- rashomon_models$variables[[j]]
      var_counts[vars] <- var_counts[vars] + 1
    }
    
    # Top 10 variables
    top_vars <- order(var_counts, decreasing = TRUE)[1:10]
    
    tibble(
      sample_size = sample_info$sample_size,
      rashomon_size = nrow(rashomon_models),
      top_var_1 = top_vars[1],
      top_var_2 = top_vars[2],
      top_var_3 = top_vars[3],
      top_var_1_count = var_counts[top_vars[1]],
      top_var_2_count = var_counts[top_vars[2]],
      top_var_3_count = var_counts[top_vars[3]]
    )
  })
  
  return(var_analysis)
}

# ===== MAIN EXECUTION =====

cat("=== RASHOMON EFFECT SIMULATION ===\n")
cat("Creating fixed dataset and sampling at different sizes\n")
cat("This demonstrates: Larger samples → Fewer models in Rashomon set\n\n")

# Run the simulation
cat("Running simulation (this may take several minutes)...\n")
results <- run_rashomon_simulation()

# Display key results
cat("\n=== RESULTS SUMMARY ===\n")
results_summary <- results %>%
  select(sample_size, threshold, total_models, rashomon_size, rashomon_proportion) %>%
  arrange(threshold, sample_size) %>%
  mutate(
    rashomon_proportion = paste0(round(rashomon_proportion * 100, 1), "%")
  )

print(results_summary, n = Inf)

# Show the key insight: larger samples -> smaller Rashomon sets
cat("\n=== KEY INSIGHT: Sample Size vs Rashomon Set Size ===\n")
key_insight <- results %>%
  select(sample_size, threshold, rashomon_size) %>%
  arrange(threshold, sample_size)

print(key_insight, n = Inf)

# Analyze variable selection stability
cat("\n=== VARIABLE SELECTION STABILITY ===\n")
var_stability <- analyze_variable_stability(results)
print(var_stability)

# ===== SAVE RESULTS TO CSV =====

# Create data directory if it doesn't exist
if (!dir.exists("../data")) {
  dir.create("../data", recursive = TRUE)
}

# Prepare main results for CSV export
results_for_csv <- results %>%
  select(sample_size, threshold, rashomon_size, rashomon_proportion, 
         best_rss, total_models) %>%
  mutate(
    rashomon_proportion = round(rashomon_proportion, 4),
    best_rss = round(best_rss, 2)
  )

write.csv(results_for_csv, file = '../data/rashomon_results.csv', row.names = FALSE)
cat("\nMain results saved to '../data/rashomon_results.csv'\n")

# Save variable stability analysis
write.csv(var_stability, file = '../data/variable_stability.csv', row.names = FALSE)
cat("Variable stability analysis saved to '../data/variable_stability.csv'\n")

# Save detailed models for a few key sample sizes (2% threshold)
key_sample_sizes <- c(1000, 4000, 12000, 24000)
detailed_models_list <- list()

for(sample_sz in key_sample_sizes) {
  result_row <- results %>%
    filter(sample_size == sample_sz, threshold == 0.02)
  
  if(nrow(result_row) > 0) {
    rashomon_models <- result_row$results[[1]]$rashomon_models
    
    if(nrow(rashomon_models) > 0) {
      models_flat <- rashomon_models %>%
        mutate(
          sample_size = sample_sz,
          threshold = 0.02,
          model_id = row_number(),
          var1 = map_dbl(variables, ~.x[1]),
          var2 = map_dbl(variables, ~.x[2]), 
          var3 = map_dbl(variables, ~.x[3]),
          var4 = map_dbl(variables, ~.x[4]),
          var5 = map_dbl(variables, ~.x[5]),
          coef1 = map_dbl(coefficients, ~.x[1]),
          coef2 = map_dbl(coefficients, ~.x[2]),
          coef3 = map_dbl(coefficients, ~.x[3]),
          coef4 = map_dbl(coefficients, ~.x[4]),
          coef5 = map_dbl(coefficients, ~.x[5])
        ) %>%
        select(sample_size, threshold, model_id, var1:var5, coef1:coef5, intercept, rss)
      
      detailed_models_list[[length(detailed_models_list) + 1]] <- models_flat
    }
  }
}

if(length(detailed_models_list) > 0) {
  all_detailed_models <- bind_rows(detailed_models_list)
  write.csv(all_detailed_models, file = '../data/rashomon_detailed_models.csv', row.names = FALSE)
  cat("Detailed models saved to '../data/rashomon_detailed_models.csv'\n")
}

# Show examples of different "stories" for small vs large samples
cat("\n=== EXAMPLE: DIFFERENT STORIES FROM SAME DATA ===\n")

# Small sample stories
small_sample_result <- results %>%
  filter(sample_size == 1000, threshold == 0.02) %>%
  pull(results) %>%
  .[[1]]

if(nrow(small_sample_result$rashomon_models) >= 3) {
  cat("Small sample (n=1000) - Multiple good stories:\n")
  stories_small <- extract_different_stories(small_sample_result$rashomon_models, 3)
  
  for (i in seq_along(stories_small)) {
    cat("Story", i, ":", stories_small[[i]]$equation, "\n")
    cat("RSS:", round(stories_small[[i]]$rss, 2), "\n")
    cat("Variables: x", paste(stories_small[[i]]$variables, collapse = ", x"), "\n\n")
  }
}

# Large sample stories
large_sample_result <- results %>%
  filter(sample_size == 24000, threshold == 0.02) %>%
  pull(results) %>%
  .[[1]]

if(nrow(large_sample_result$rashomon_models) >= 1) {
  cat("Large sample (n=24000) - Fewer good stories:\n")
  n_stories_large <- min(3, nrow(large_sample_result$rashomon_models))
  stories_large <- extract_different_stories(large_sample_result$rashomon_models, n_stories_large)
  
  for (i in seq_along(stories_large)) {
    cat("Story", i, ":", stories_large[[i]]$equation, "\n")
    cat("RSS:", round(stories_large[[i]]$rss, 2), "\n")
    cat("Variables: x", paste(stories_large[[i]]$variables, collapse = ", x"), "\n\n")
  }
}

cat("=== SIMULATION COMPLETE ===\n")
cat("Key finding: As sample size increases from the same population,\n")
cat("the Rashomon set shrinks, reducing model ambiguity.\n")
cat("\nFiles saved:\n")
cat("- ../data/rashomon_results.csv (main results)\n")
cat("- ../data/variable_stability.csv (variable selection patterns)\n")
cat("- ../data/rashomon_detailed_models.csv (detailed model info)\n")