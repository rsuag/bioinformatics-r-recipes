# Task: estimate a binomial proportion with a Beta prior.
# Input: integer successes and trials; positive Beta prior shape parameters.
# Output: posterior mean and equal-tail credible interval.
# Interpretation depends on the stated prior and a binomial sampling model.

beta_binomial <- function(successes, trials, prior_a = 1, prior_b = 1, level = 0.95) {
  if (length(successes) != 1L || length(trials) != 1L ||
      !is.numeric(successes) || !is.numeric(trials) || anyNA(c(successes, trials)) ||
      !all(is.finite(c(successes, trials))) || successes < 0 || trials < 1 ||
      successes > trials || successes != floor(successes) || trials != floor(trials)) {
    stop("successes and trials must be integers with 0 <= successes <= trials")
  }
  if (any(lengths(list(prior_a, prior_b, level)) != 1L) ||
      !all(vapply(list(prior_a, prior_b, level), is.numeric, logical(1))) ||
      anyNA(c(prior_a, prior_b, level)) || !all(is.finite(c(prior_a, prior_b, level))) ||
      prior_a <= 0 || prior_b <= 0 || level <= 0 || level >= 1) {
    stop("prior shapes must be positive and level must be between 0 and 1")
  }
  a <- prior_a + successes
  b <- prior_b + trials - successes
  tail <- (1 - level) / 2
  c(posterior_mean = a / (a + b), lower = qbeta(tail, a, b),
    upper = qbeta(1 - tail, a, b))
}

if (sys.nframe() == 0L) {
  print(beta_binomial(successes = 7, trials = 20))
}
