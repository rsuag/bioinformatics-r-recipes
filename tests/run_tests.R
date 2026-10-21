source("recipes/explore_counts.R")
source("recipes/adjust_pvalues.R")
source("recipes/beta_binomial.R")
source("recipes/logistic_classifier.R")
source("recipes/fasta_gc.R")
source("recipes/metabolomics_missingness.R")
source("recipes/kmer_frequencies.R")
source("recipes/variant_allele_frequency.R")
source("recipes/differential_expression.R")
source("recipes/gene_set_enrichment.R")
source("recipes/cluster_samples.R")
source("recipes/markov_transitions.R")
source("recipes/metabolomics_fold_change.R")
source("recipes/bootstrap_mean.R")

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

kmers <- kmer_frequencies("ACACNAC", 2)
stopifnot(unname(kmers["AC"]) == 3L, sum(kmers) == 4L)

genotypes <- matrix(c(0, 1, 2, NA, NA, NA), nrow = 2, byrow = TRUE)
frequencies <- variant_allele_frequency(genotypes)
stopifnot(frequencies$alternate_frequency[1] == 0.5,
          is.na(frequencies$alternate_frequency[2]),
          frequencies$missing_fraction[2] == 1)

expression <- matrix(c(1, 3, 2, 4, 4, 2, 5, 3), nrow = 2)
de <- differential_expression(expression, c("control", "control", "treated", "treated"))
stopifnot(de$log2_difference[1] == 3, de$log2_difference[2] == -1,
          all(de$q_value >= de$p_value))

enrichment <- gene_set_enrichment(c("a", "b"),
                                  list(both = c("a", "b"), none = "d"), letters[1:6])
stopifnot(enrichment$overlap[1] == 2L, enrichment$p_value[2] == 1,
          enrichment$p_value[1] < enrichment$p_value[2])

profiles <- matrix(c(1, 2, 3, 4, 5, 2, 3, 4, 5, 6), nrow = 2)
clusters <- cluster_samples(profiles)
stopifnot(nrow(clusters$distance) == 5L,
          all(diag(clusters$distance) == 0),
          isTRUE(all.equal(clusters$distance, t(clusters$distance))))

transitions <- markov_transitions(c("A", "B", "A", "A"), c("A", "B", "C"))
stopifnot(transitions$counts["A", "A"] == 1L,
          transitions$counts["A", "B"] == 1L,
          all(rowSums(transitions$probabilities)[1:2] == 1),
          sum(transitions$probabilities["C", ]) == 0)

metabolites <- matrix(c(10, 5, 10, NA, 20, 10, 20, 20), nrow = 2)
fold_change <- metabolomics_fold_change(metabolites,
                                        c("control", "control", "treated", "treated"))
stopifnot(abs(fold_change$fold_change[1] - 2) < 1e-12,
          fold_change$n_first[2] == 1L)

boot <- bootstrap_mean(rep(5, 4), replicates = 100, seed = 4)
stopifnot(all(boot == 5))

cat("All recipe tests passed.\n")
