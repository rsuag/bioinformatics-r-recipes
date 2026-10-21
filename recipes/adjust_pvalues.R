# Task: control false discoveries across a family of tests.
# Input: named or unnamed numeric p-values in [0, 1].
# Output: original p-values and Benjamini-Hochberg adjusted q-values.

adjust_pvalues <- function(p_values, alpha = 0.05) {
  if (!is.numeric(p_values) || !length(p_values) || anyNA(p_values) ||
      any(!is.finite(p_values)) || any(p_values < 0 | p_values > 1)) {
    stop("p_values must be finite numbers between 0 and 1")
  }
  if (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number between 0 and 1")
  }
  labels <- names(p_values)
  if (is.null(labels)) labels <- paste0("test_", seq_along(p_values))
  q_values <- p.adjust(p_values, method = "BH")
  data.frame(test = labels, p_value = unname(p_values), q_value = unname(q_values),
             significant = unname(q_values <= alpha), row.names = NULL)
}

if (sys.nframe() == 0L) {
  print(adjust_pvalues(c(gene_a = 0.001, gene_b = 0.08, gene_c = 0.02, gene_d = 0.7)))
}
