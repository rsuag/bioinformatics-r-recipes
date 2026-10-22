# Render one result figure for every standalone recipe.
# Run from the repository root:
# Rscript scripts/render_recipe_figures.R [output_directory]

recipe_files <- c(
  "adjust_pvalues.R", "beta_binomial.R", "bootstrap_mean.R",
  "cluster_samples.R", "differential_expression.R", "explore_counts.R",
  "fasta_gc.R", "gene_set_enrichment.R", "kmer_frequencies.R",
  "logistic_classifier.R", "markov_transitions.R",
  "metabolomics_fold_change.R", "metabolomics_missingness.R",
  "variant_allele_frequency.R"
)
invisible(lapply(file.path("recipes", recipe_files), source))

blue <- "#2369A8"
orange <- "#D65F35"
teal <- "#2A9D8F"
gray <- "#747474"

with_png <- function(path, code) {
  grDevices::png(path, width = 960, height = 600, res = 120)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(4.5, 4.5, 3, 1), las = 1, bty = "l")
  force(code)
}

render_recipe_figures <- function(output_dir = "docs/figures/recipes") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  adjusted <- adjust_pvalues(c(gene_a = 0.001, gene_b = 0.008, gene_c = 0.02,
                               gene_d = 0.08, gene_e = 0.3, gene_f = 0.7))
  with_png(file.path(output_dir, "adjust_pvalues.png"), {
    values <- t(as.matrix(adjusted[c("p_value", "q_value")]))
    barplot(values, beside = TRUE, names.arg = adjusted$test,
            col = c(blue, orange), ylab = "Probability", main = "Multiple-testing adjustment")
    abline(h = 0.05, lty = 2, col = gray)
    legend("topleft", c("p-value", "BH q-value"), fill = c(blue, orange), bty = "n")
  })

  posterior <- beta_binomial(7, 20)
  with_png(file.path(output_dir, "beta_binomial.png"), {
    x <- seq(0, 1, length.out = 400)
    y <- dbeta(x, 8, 14)
    plot(x, y, type = "l", lwd = 3, col = blue, xlab = "Proportion",
         ylab = "Posterior density", main = "Beta-binomial posterior")
    polygon(c(x[x >= posterior["lower"] & x <= posterior["upper"]],
              rev(x[x >= posterior["lower"] & x <= posterior["upper"]])),
            c(y[x >= posterior["lower"] & x <= posterior["upper"]],
              rep(0, sum(x >= posterior["lower"] & x <= posterior["upper"]))),
            col = grDevices::adjustcolor(blue, 0.25), border = NA)
    abline(v = posterior["posterior_mean"], col = orange, lwd = 2)
  })

  measurements <- c(4.2, 5.1, 3.9, 4.8, 5.3, 4.5, 5.0, 4.1)
  boot <- bootstrap_mean(measurements, replicates = 2000, seed = 11)
  with_png(file.path(output_dir, "bootstrap_mean.png"), {
    stripchart(measurements, method = "jitter", vertical = FALSE, pch = 19,
               col = blue, xlab = "Measurement", main = "Mean and 95% bootstrap interval")
    segments(boot["lower"], 1.18, boot["upper"], 1.18, lwd = 5, col = orange)
    points(boot["mean"], 1.18, pch = 18, cex = 1.6, col = "black")
  })

  profiles <- matrix(c(2, 3, 8, 9, 2, 4, 7, 9, 8, 7, 3, 2,
                       7, 8, 2, 1, 3, 4, 8, 7), nrow = 4,
                     dimnames = list(paste0("feature_", 1:4), paste0("sample_", 1:5)))
  clustered <- cluster_samples(profiles)
  with_png(file.path(output_dir, "cluster_samples.png"), {
    plot(clustered$tree, hang = -1, col = blue, main = "Hierarchical sample clustering",
         xlab = "", sub = "Complete linkage on standardized features")
  })

  set.seed(31)
  expression <- matrix(rnorm(80 * 8, 6, 0.45), nrow = 80,
                       dimnames = list(paste0("gene_", 1:80), NULL))
  expression[1:10, 5:8] <- expression[1:10, 5:8] + 1.8
  de <- differential_expression(expression, rep(c("control", "treated"), each = 4))
  with_png(file.path(output_dir, "differential_expression.png"), {
    highlighted <- !is.na(de$q_value) & de$q_value < 0.05 & abs(de$log2_difference) >= 1
    plot(de$log2_difference, -log10(de$p_value), pch = 19,
         col = ifelse(highlighted, orange, gray), xlab = "Mean log2 difference",
         ylab = "-log10(p-value)", main = "Normalized-expression screen")
    abline(v = c(-1, 1), lty = 2, col = "#BBBBBB")
  })

  set.seed(4)
  counts <- matrix(rpois(60 * 6, 50), nrow = 60,
                   dimnames = list(paste0("gene_", 1:60), paste0("sample_", 1:6)))
  counts[1:15, 4:6] <- counts[1:15, 4:6] + 65
  explored <- explore_counts(counts)
  with_png(file.path(output_dir, "explore_counts.png"), {
    colors <- rep(c(blue, orange), each = 3)
    plot(explored$pca, pch = 19, cex = 1.5, col = colors, xlab = "PC1", ylab = "PC2",
         main = "Count-matrix sample structure")
    text(explored$pca, labels = rownames(explored$pca),
         pos = rep(c(1, 3, 4, 2), length.out = nrow(explored$pca)), cex = 0.8)
  })

  gc <- fasta_gc("data/toy_sequences.fasta")
  with_png(file.path(output_dir, "fasta_gc.png"), {
    barplot(gc$gc_fraction, names.arg = gc$id, col = blue, ylim = c(0, 1),
            ylab = "GC fraction", main = "GC content by FASTA sequence")
  })

  universe <- paste0("g", 1:100)
  selected <- paste0("g", 1:15)
  enrichment <- gene_set_enrichment(selected,
    list(nitrogen_cycle = paste0("g", c(1:10, 41:50)),
         carbon_cycle = paste0("g", c(8:20, 51:60)),
         stress_response = paste0("g", 61:85)), universe)
  with_png(file.path(output_dir, "gene_set_enrichment.png"), {
    par(mar = c(4.5, 9, 3, 1))
    score <- -log10(pmax(enrichment$q_value, .Machine$double.xmin))
    barplot(score, names.arg = enrichment$set, horiz = TRUE, las = 1, col = teal,
            xlab = "-log10(BH q-value)", main = "Gene-set overrepresentation")
  })

  kmers <- kmer_frequencies("ACGTACGTACGTTTACGCGCGTATATACG", 2)
  with_png(file.path(output_dir, "kmer_frequencies.png"), {
    barplot(kmers, col = blue, ylab = "Count", xlab = "2-mer",
            main = "Observed DNA k-mer frequencies")
  })

  set.seed(24)
  classifier_data <- data.frame(signal_a = rnorm(100), signal_b = rnorm(100))
  classifier_data$case <- rbinom(100, 1,
    plogis(-0.2 + 1.4 * classifier_data$signal_a - 0.7 * classifier_data$signal_b))
  classifier <- evaluate_logistic(classifier_data, "case", c("signal_a", "signal_b"), seed = 24)
  with_png(file.path(output_dir, "logistic_classifier.png"), {
    jittered <- jitter(classifier$predictions$actual, amount = 0.08)
    plot(jittered, classifier$predictions$probability, pch = 19,
         col = ifelse(classifier$predictions$actual == 1, orange, blue),
         xaxt = "n", xlab = "Actual class", ylab = "Predicted probability",
         main = "Held-out classifier predictions", xlim = c(-0.35, 1.35))
    axis(1, at = 0:1, labels = c("control", "case"))
    abline(h = 0.5, lty = 2, col = gray)
  })

  transitions <- markov_transitions(c("A", "C", "A", "A", "G", "C", "A", "G", "G", "C"),
                                    c("A", "C", "G"))
  with_png(file.path(output_dir, "markov_transitions.png"), {
    matrix_for_plot <- transitions$probabilities[nrow(transitions$probabilities):1, , drop = FALSE]
    image(seq_len(ncol(matrix_for_plot)), seq_len(nrow(matrix_for_plot)), t(matrix_for_plot),
          col = grDevices::colorRampPalette(c("white", blue))(20), zlim = c(0, 1),
          xaxt = "n", yaxt = "n", xlab = "To state", ylab = "From state",
          main = "First-order transition probabilities")
    axis(1, at = 1:3, labels = colnames(matrix_for_plot))
    axis(2, at = 1:3, labels = rev(rownames(matrix_for_plot)))
  })

  metabolite_values <- matrix(c(10, 20, 5, 12, 11, 19, 6, 11,
                                22, 10, 10, 24, 20, 11, 9, 25), nrow = 4,
                              dimnames = list(paste0("metabolite_", 1:4), NULL))
  metabolite_fc <- metabolomics_fold_change(metabolite_values,
                                             c("control", "control", "treated", "treated"))
  with_png(file.path(output_dir, "metabolomics_fold_change.png"), {
    colors <- ifelse(metabolite_fc$log2_difference >= 0, orange, blue)
    barplot(metabolite_fc$log2_difference, names.arg = metabolite_fc$feature,
            col = colors, ylab = "Median log2 difference", main = "Metabolite fold changes")
    abline(h = 0, col = "black")
  })

  intensity <- matrix(c(100, 200, NA, 80, NA, 320, 90, 210,
                        400, 120, 190, NA, 130, NA, 220, 300), nrow = 4,
                      dimnames = list(paste0("feature_", 1:4), paste0("sample_", 1:4)))
  missingness <- metabolomics_missingness(intensity)
  with_png(file.path(output_dir, "metabolomics_missingness.png"), {
    map <- is.na(intensity)[nrow(intensity):1, , drop = FALSE]
    image(seq_len(ncol(map)), seq_len(nrow(map)), t(map), col = c("#E8EEF4", orange),
          xaxt = "n", yaxt = "n", xlab = "Sample", ylab = "Feature",
          main = "Missing metabolomics measurements")
    axis(1, at = seq_len(ncol(map)), labels = colnames(map))
    axis(2, at = seq_len(nrow(map)), labels = rev(rownames(map)))
    legend("topright", c("Observed", "Missing"), fill = c("#E8EEF4", orange), bty = "n")
  })

  genotypes <- matrix(c(0, 1, 2, 0, 1, 1, 2, 2, NA, 0,
                        0, 0, 1, 0, 0, 2, 2, 2, 1, NA), nrow = 4, byrow = TRUE,
                      dimnames = list(paste0("variant_", 1:4), paste0("sample_", 1:5)))
  variants <- variant_allele_frequency(genotypes)
  with_png(file.path(output_dir, "variant_allele_frequency.png"), {
    mids <- barplot(variants$alternate_frequency, names.arg = variants$variant,
                    col = blue, ylim = c(0, 1), ylab = "Alternate allele frequency",
                    main = "Diploid variant summary")
    points(mids, variants$missing_fraction, pch = 18, col = orange, cex = 1.3)
    legend("topleft", c("Allele frequency", "Missing fraction"),
           fill = c(blue, NA), border = c(NA, NA), pch = c(NA, 18),
           col = c(NA, orange), bty = "n")
  })

  invisible(normalizePath(output_dir))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- if (length(args)) args[1] else "docs/figures/recipes"
  render_recipe_figures(output_dir)
  message("Rendered 14 recipe figures in ", output_dir)
}
