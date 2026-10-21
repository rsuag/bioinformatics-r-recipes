# Task: compare normalized log2 expression between two biological groups.
# Input: genes x samples numeric matrix and a two-level sample group vector.
# Output: mean log2 difference, Welch p-value, and BH-adjusted q-value.
# Educational screen only; not a raw-count differential-expression model.

differential_expression <- function(log2_expression, group) {
  x <- log2_expression
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || ncol(x) < 4L ||
      anyNA(x) || any(!is.finite(x)) || length(group) != ncol(x) ||
      anyNA(group) || !is.factor(group) && !is.character(group)) {
    stop("supply a finite genes-by-samples matrix and a two-group label vector")
  }
  group <- factor(group)
  if (nlevels(group) != 2L || any(tabulate(group) < 2L)) {
    stop("each of exactly two groups needs at least two biological samples")
  }
  first <- which(group == levels(group)[1])
  second <- which(group == levels(group)[2])
  difference <- rowMeans(x[, second, drop = FALSE]) - rowMeans(x[, first, drop = FALSE])
  p <- apply(x, 1L, function(row) {
    a <- row[first]; b <- row[second]
    if (var(a) == 0 && var(b) == 0) return(if (a[1] == b[1]) 1 else NA_real_)
    t.test(b, a)$p.value
  })
  ids <- rownames(x)
  if (is.null(ids)) ids <- paste0("gene_", seq_len(nrow(x)))
  data.frame(gene = ids, log2_difference = difference, p_value = p,
             q_value = p.adjust(p, "BH"), row.names = NULL)
}

if (sys.nframe() == 0L) {
  x <- matrix(c(4, 6, 8, 4.2, 6.1, 8.1, 7, 6, 8, 7.2, 6.2, 8.2),
              nrow = 3, dimnames = list(c("gene_a", "gene_b", "gene_c"), NULL))
  print(differential_expression(x, c("control", "control", "treated", "treated")))
}
