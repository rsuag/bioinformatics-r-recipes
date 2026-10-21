# Task: summarize missing measurements in a metabolomics intensity matrix.
# Input: numeric matrix with features in rows and samples in columns.
# Output: missingness by feature and sample, plus log2 of observed positive values.
# Missing values remain missing; this recipe does not impute or infer mechanisms.

metabolomics_missingness <- function(intensity) {
  if (!is.matrix(intensity) || !is.numeric(intensity) || !nrow(intensity) || !ncol(intensity) ||
      any(!is.finite(intensity[!is.na(intensity)])) || any(intensity[!is.na(intensity)] <= 0)) {
    stop("intensity must be a numeric matrix of positive observed values and NA")
  }
  feature_missing <- rowMeans(is.na(intensity))
  sample_missing <- colMeans(is.na(intensity))
  log2_intensity <- log2(intensity)
  list(feature_missing = feature_missing, sample_missing = sample_missing,
       log2_intensity = log2_intensity)
}

if (sys.nframe() == 0L) {
  intensity <- matrix(c(100, 200, NA, 80, NA, 320, 90, 210, 400, 120, 190, NA),
                      nrow = 4, dimnames = list(paste0("feature_", 1:4), paste0("sample_", 1:3)))
  result <- metabolomics_missingness(intensity)
  print(result$feature_missing)
  print(result$sample_missing)
}
