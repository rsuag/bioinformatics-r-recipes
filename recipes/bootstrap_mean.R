# Task: estimate uncertainty in a sample mean by nonparametric bootstrap.
# Input: finite numeric measurements, number of replicates, and random seed.
# Output: observed mean and percentile confidence interval.
# Resample at the independent biological-unit level; this is not a hierarchical bootstrap.

bootstrap_mean <- function(values, replicates = 1000L, level = 0.95, seed = 1L) {
  if (!is.numeric(values) || length(values) < 2L || anyNA(values) ||
      any(!is.finite(values)) || !is.numeric(replicates) || length(replicates) != 1L ||
      is.na(replicates) || !is.finite(replicates) || replicates < 2 ||
      replicates != floor(replicates) || !is.numeric(level) || length(level) != 1L ||
      is.na(level) || !is.finite(level) || level <= 0 || level >= 1 ||
      !is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != floor(seed)) {
    stop("values must be finite; use >=2 integer replicates, level in (0,1), and integer seed")
  }
  set.seed(seed)
  draws <- replicate(replicates, mean(sample(values, length(values), replace = TRUE)))
  interval <- quantile(draws, c((1 - level) / 2, (1 + level) / 2), names = FALSE)
  c(mean = mean(values), lower = interval[1], upper = interval[2])
}

if (sys.nframe() == 0L) print(bootstrap_mean(c(4.2, 5.1, 3.9, 4.8, 5.3)))
