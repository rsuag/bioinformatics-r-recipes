# Task: check a genes-by-samples count matrix and explore sample structure.
# Input: numeric matrix with genes in rows, samples in columns, and nonnegative counts.
# Output: library sizes, log2(CPM + 1), and sample PCA coordinates.
# Exploratory only: this is not a differential-expression analysis.

explore_counts <- function(counts) {
  if (!is.matrix(counts) || !is.numeric(counts) || nrow(counts) < 2L || ncol(counts) < 3L ||
      anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop("counts must be a finite, nonnegative numeric matrix with >=2 genes and >=3 samples")
  }
  library_sizes <- colSums(counts)
  if (any(library_sizes == 0)) stop("every sample must have at least one count")
  cpm <- sweep(counts, 2L, library_sizes, "/") * 1e6
  log_cpm <- log2(cpm + 1)
  pca <- prcomp(t(log_cpm), center = TRUE, scale. = FALSE)
  list(library_sizes = library_sizes, log_cpm = log_cpm, pca = pca$x[, 1:2, drop = FALSE])
}

if (sys.nframe() == 0L) {
  counts <- matrix(c(20, 8, 0, 4, 13, 2, 9, 7,
                     23, 6, 1, 5, 15, 2, 8, 6,
                     19, 9, 0, 3, 12, 3, 10, 8,
                     4, 18, 8, 16, 2, 12, 3, 13,
                     5, 16, 9, 17, 3, 11, 2, 14,
                     3, 19, 7, 15, 1, 13, 4, 12), nrow = 8)
  rownames(counts) <- paste0("gene_", 1:8)
  colnames(counts) <- paste0("sample_", 1:6)
  result <- explore_counts(counts)
  print(result$library_sizes)
  print(round(result$pca, 3))
}
