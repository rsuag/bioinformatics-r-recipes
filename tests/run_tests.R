source("recipes/explore_counts.R")
source("recipes/adjust_pvalues.R")
source("recipes/beta_binomial.R")
source("recipes/logistic_classifier.R")
source("recipes/fasta_gc.R")
source("recipes/metabolomics_missingness.R")

counts <- matrix(c(1, 2, 3, 2, 4, 1, 4, 1, 2), nrow = 3)
explored <- explore_counts(counts)
stopifnot(length(explored$library_sizes) == 3L,
          identical(dim(explored$pca), c(3L, 2L)),
          all(abs(colSums(2^explored$log_cpm - 1) - 1e6) < 1e-6))

adjusted <- adjust_pvalues(c(a = 0.01, b = 0.02, c = 0.5))
stopifnot(identical(adjusted$test, c("a", "b", "c")),
          all(diff(adjusted$q_value) >= 0), adjusted$q_value[1] == 0.03)

posterior <- beta_binomial(0, 1)
stopifnot(abs(posterior["posterior_mean"] - 1 / 3) < 1e-12,
          posterior["lower"] < posterior["posterior_mean"],
          posterior["upper"] > posterior["posterior_mean"])

gc <- fasta_gc("data/toy_sequences.fasta")
stopifnot(nrow(gc) == 3L, gc$gc_fraction[1] == 0.5,
          gc$called_bases[2] == 6L)

intensity <- matrix(c(2, NA, 4, 8), nrow = 2)
missing <- metabolomics_missingness(intensity)
stopifnot(identical(as.numeric(missing$feature_missing), c(0, 0.5)),
          is.na(missing$log2_intensity[2, 1]))

set.seed(8)
classifier_data <- data.frame(x = rnorm(80), z = rnorm(80))
classifier_data$case <- rbinom(80, 1, plogis(classifier_data$x - 0.3 * classifier_data$z))
model_result <- evaluate_logistic(classifier_data, "case", c("x", "z"), seed = 8)
stopifnot(!length(intersect(model_result$train_rows, model_result$test_rows)),
          length(model_result$test_rows) == nrow(model_result$predictions),
          sum(model_result$confusion) == length(model_result$test_rows),
          all(model_result$metrics >= 0 & model_result$metrics <= 1))

cat("All recipe tests passed.\n")
