# Create a small synthetic RNA-seq count matrix and sample sheet.
# Values are simulated for teaching and are not biological observations.

dir.create("workflows/rnaseq/example", recursive = TRUE, showWarnings = FALSE)
set.seed(42)
sample_ids <- paste0("sample_", 1:6)
condition <- rep(c("control", "treated"), each = 3)
counts <- matrix(rnbinom(120 * 6, mu = 80, size = 15), nrow = 120,
                 dimnames = list(paste0("gene_", 1:120), sample_ids))
counts[1:12, 4:6] <- matrix(rnbinom(12 * 3, mu = 220, size = 15), nrow = 12)
write.csv(data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
          "workflows/rnaseq/example/counts.csv", row.names = FALSE)
write.csv(data.frame(sample_id = sample_ids, condition = condition),
          "workflows/rnaseq/example/samples.csv", row.names = FALSE)
