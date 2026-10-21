# Task: compare observed metabolite intensities between two groups.
# Input: positive features x samples matrix (NA allowed) and two group labels.
# Output: per-group observed counts, median log2 difference, and fold change.
# Descriptive only; no imputation or significance test is performed.

metabolomics_fold_change <- function(intensity, group) {
  if (!is.matrix(intensity) || !is.numeric(intensity) || !nrow(intensity) ||
      length(group) != ncol(intensity) || anyNA(group) ||
      !is.character(group) && !is.factor(group) ||
      any(!is.finite(intensity[!is.na(intensity)])) ||
      any(intensity[!is.na(intensity)] <= 0)) {
    stop("supply positive features-by-samples intensities (NA allowed) and group labels")
  }
  group <- factor(group)
  if (nlevels(group) != 2L || any(tabulate(group) < 2L)) {
    stop("exactly two groups with >=2 samples each are required")
  }
  log_values <- log2(intensity)
  by_group <- lapply(levels(group), function(label) {
    x <- log_values[, group == label, drop = FALSE]
    list(n = rowSums(!is.na(x)), median = apply(x, 1L, function(v) {
      if (all(is.na(v))) NA_real_ else median(v, na.rm = TRUE)
    }))
  })
  delta <- by_group[[2]]$median - by_group[[1]]$median
  ids <- rownames(intensity)
  if (is.null(ids)) ids <- paste0("feature_", seq_len(nrow(intensity)))
  data.frame(feature = ids, n_first = by_group[[1]]$n, n_second = by_group[[2]]$n,
             log2_difference = delta, fold_change = 2^delta, row.names = NULL)
}

if (sys.nframe() == 0L) {
  x <- matrix(c(10, 5, 12, NA, 20, 10, 22, 12), nrow = 2)
  print(metabolomics_fold_change(x, c("control", "control", "treated", "treated")))
}
