# Task: group samples by their standardized molecular profiles.
# Input: features x samples finite numeric matrix; at least 2 each.
# Output: sample distance matrix and hierarchical clustering object.
# Constant features are omitted; clustering is descriptive, not validation.

cluster_samples <- function(values) {
  if (!is.matrix(values) || !is.numeric(values) || nrow(values) < 2L ||
      ncol(values) < 2L || anyNA(values) || any(!is.finite(values))) {
    stop("values must be a finite features-by-samples matrix with >=2 rows and columns")
  }
  variable <- apply(values, 1L, sd) > 0
  if (!any(variable)) stop("at least one feature must vary across samples")
  standardized <- t(scale(t(values[variable, , drop = FALSE])))
  distances <- dist(t(standardized))
  list(distance = as.matrix(distances), tree = hclust(distances, method = "complete"),
       retained_features = which(variable))
}

if (sys.nframe() == 0L) {
  x <- matrix(c(1, 2, 3, 4, 1, 2, 4, 5, 8, 8, 7, 6), nrow = 3)
  print(cluster_samples(x)$tree)
}
