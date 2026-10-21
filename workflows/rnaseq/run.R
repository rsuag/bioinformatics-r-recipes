# Run from the repository root:
# Rscript workflows/rnaseq/run.R counts.csv samples.csv output_dir control treated
# Counts: gene_id column followed by nonnegative integer sample counts.
# Samples: sample_id and condition columns; >=2 biological replicates per group.

read_rnaseq_inputs <- function(count_path, sample_path, reference, case) {
  counts_table <- read.csv(count_path, check.names = FALSE, stringsAsFactors = FALSE)
  samples <- read.csv(sample_path, stringsAsFactors = FALSE)
  if (ncol(counts_table) < 5L || !identical(names(counts_table)[1], "gene_id") ||
      nrow(counts_table) < 2L) {
    stop("counts must have gene_id and >=4 sample columns, with >=2 genes")
  }
  genes <- counts_table$gene_id
  ids <- names(counts_table)[-1L]
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes) ||
      any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("gene and sample IDs must be nonempty and unique")
  }
  if (!all(c("sample_id", "condition") %in% names(samples)) ||
      nrow(samples) != length(ids) || anyNA(samples$sample_id) ||
      anyDuplicated(samples$sample_id) || !setequal(samples$sample_id, ids)) {
    stop("samples must contain one sample_id and condition for every count column")
  }
  if (!is.character(reference) || length(reference) != 1L || !nzchar(reference) ||
      !is.character(case) || length(case) != 1L || !nzchar(case) ||
      identical(reference, case)) stop("supply distinct reference and case labels")
  samples <- samples[match(ids, samples$sample_id), , drop = FALSE]
  condition <- samples$condition
  if (anyNA(condition) || !setequal(unique(condition), c(reference, case)) ||
      any(table(factor(condition, levels = c(reference, case))) < 2L)) {
    stop("samples need exactly the reference and case groups, each with >=2 samples")
  }
  count_columns <- counts_table[-1L]
  if (!all(vapply(count_columns, is.numeric, logical(1)))) {
    stop("all sample count columns must be numeric")
  }
  counts <- as.matrix(count_columns)
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts)) || any(colSums(counts) == 0)) {
    stop("counts must be finite nonnegative integers and each sample needs counts")
  }
  rownames(counts) <- genes
  colnames(counts) <- ids
  rownames(samples) <- samples$sample_id
  samples$condition <- factor(condition, levels = c(reference, case))
  list(counts = counts, samples = samples)
}

write_qc <- function(counts, samples, out_dir) {
  library_sizes <- colSums(counts)
  write.csv(data.frame(sample_id = names(library_sizes),
                       condition = samples$condition, total_counts = as.numeric(library_sizes)),
            file.path(out_dir, "library_sizes.csv"), row.names = FALSE)
  log_cpm <- log2(t(t(counts) / library_sizes * 1e6) + 1)
  variable <- apply(log_cpm, 1L, var) > 0
  if (sum(variable) < 2L) stop("PCA needs at least two variable genes")
  pca <- prcomp(t(log_cpm[variable, , drop = FALSE]))
  coordinates <- data.frame(sample_id = rownames(pca$x),
                            condition = samples$condition,
                            PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  write.csv(coordinates, file.path(out_dir, "pca_coordinates.csv"), row.names = FALSE)
  png(file.path(out_dir, "pca.png"), width = 900, height = 700)
  on.exit(dev.off(), add = TRUE)
  colors <- c("#2369a8", "#d65f35")[as.integer(samples$condition)]
  plot(coordinates$PC1, coordinates$PC2, col = colors, pch = 19, cex = 1.5,
       xlab = "PC1", ylab = "PC2", main = "RNA-seq exploratory PCA (log2 CPM + 1)")
  text(coordinates$PC1, coordinates$PC2, labels = coordinates$sample_id,
       pos = 3, cex = 0.8)
  legend("topright", legend = levels(samples$condition),
         col = c("#2369a8", "#d65f35"), pch = 19, bty = "n")
  invisible(coordinates)
}

write_volcano <- function(results, out_dir) {
  valid <- is.finite(results$log2FoldChange) & is.finite(results$pvalue)
  if (!any(valid)) stop("DESeq2 produced no plottable genes")
  x <- results$log2FoldChange[valid]
  y <- -log10(pmax(results$pvalue[valid], .Machine$double.xmin))
  marked <- !is.na(results$padj[valid]) & results$padj[valid] < 0.05 & abs(x) >= 1
  png(file.path(out_dir, "volcano.png"), width = 900, height = 700)
  on.exit(dev.off(), add = TRUE)
  plot(x, y, col = ifelse(marked, "#d65f35", "#777777"), pch = 19,
       xlab = "log2 fold change (case / reference)", ylab = "-log10(p-value)",
       main = "Differential expression: DESeq2")
  abline(v = c(-1, 1), lty = 2, col = "#999999")
  legend("topright", legend = "padj < 0.05 and |log2FC| >= 1",
         col = "#d65f35", pch = 19, bty = "n")
}

run_rnaseq <- function(count_path, sample_path, out_dir, reference, case) {
  input <- read_rnaseq_inputs(count_path, sample_path, reference, case)
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required; see workflows/rnaseq/README.md for installation")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_qc(input$counts, input$samples, out_dir)
  keep <- rowSums(input$counts >= 10L) >= 2L
  if (sum(keep) < 2L) stop("fewer than two genes remain after count filtering")
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = input$counts[keep, , drop = FALSE],
    colData = input$samples, design = ~ condition)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  result <- as.data.frame(DESeq2::results(dds, contrast = c("condition", case, reference)))
  result <- data.frame(gene_id = rownames(result), result, row.names = NULL)
  write.csv(result, file.path(out_dir, "deseq2_results.csv"), row.names = FALSE)
  write_volcano(result, out_dir)
  writeLines(c(paste("reference:", reference), paste("case:", case),
               paste("genes tested:", sum(keep)), capture.output(sessionInfo())),
             file.path(out_dir, "session_info.txt"))
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5L) {
    stop("usage: Rscript workflows/rnaseq/run.R counts.csv samples.csv output_dir reference case")
  }
  run_rnaseq(args[1], args[2], args[3], args[4], args[5])
}
