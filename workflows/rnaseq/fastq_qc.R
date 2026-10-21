# Inspect paired FASTQ files without mapping or biological interpretation.
# Run: Rscript workflows/rnaseq/fastq_qc.R R1.fastq.gz R2.fastq.gz output_dir [max_pairs]

fastq_qc <- function(r1_path, r2_path, out_dir, max_pairs = 10000L) {
  if (length(max_pairs) != 1L || !is.numeric(max_pairs) || is.na(max_pairs) ||
      !is.finite(max_pairs) || max_pairs < 1 || max_pairs != floor(max_pairs)) {
    stop("max_pairs must be a positive integer")
  }
  r1 <- gzfile(r1_path, "rt")
  r2 <- gzfile(r2_path, "rt")
  on.exit({close(r1); close(r2)}, add = TRUE)
  lengths <- matrix(NA_integer_, nrow = max_pairs, ncol = 2L)
  qualities <- matrix(NA_real_, nrow = max_pairs, ncol = 2L)
  gc <- matrix(NA_real_, nrow = max_pairs, ncol = 2L)
  n <- 0L
  for (i in seq_len(max_pairs)) {
    a <- readLines(r1, n = 4L, warn = FALSE)
    b <- readLines(r2, n = 4L, warn = FALSE)
    if (!length(a) && !length(b)) break
    if (length(a) != 4L || length(b) != 4L ||
        !startsWith(a[1], "@") || !startsWith(b[1], "@") ||
        !startsWith(a[3], "+") || !startsWith(b[3], "+")) {
      stop("incomplete or malformed paired FASTQ record")
    }
    ids <- sub("/[12]$", "", sub(" .*", "", c(a[1], b[1])))
    if (ids[1] != ids[2]) stop("mate read IDs differ")
    for (mate in 1:2) {
      record <- list(a, b)[[mate]]
      sequence <- toupper(record[2])
      phred <- utf8ToInt(record[4]) - 33L
      if (!nzchar(sequence) || grepl("[^ACGTN]", sequence) ||
          nchar(sequence) != length(phred) || any(phred < 0L | phred > 93L)) {
        stop("invalid sequence or Phred+33 quality string")
      }
      bases <- strsplit(sequence, "", fixed = TRUE)[[1]]
      called <- sum(bases != "N")
      lengths[i, mate] <- length(bases)
      qualities[i, mate] <- mean(phred)
      gc[i, mate] <- if (called) sum(bases %in% c("G", "C")) / called else NA_real_
    }
    n <- i
  }
  if (!n) stop("FASTQ files contain no read pairs")
  summary <- data.frame(mate = c("R1", "R2"), pairs_checked = n,
                        median_length = apply(lengths[seq_len(n), , drop = FALSE], 2L, median),
                        mean_phred = colMeans(qualities[seq_len(n), , drop = FALSE]),
                        mean_gc = apply(gc[seq_len(n), , drop = FALSE], 2L, function(x) {
                          if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
                        }))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(summary, file.path(out_dir, "fastq_qc_summary.csv"), row.names = FALSE)
  png(file.path(out_dir, "read_lengths.png"), width = 900, height = 500)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2))
  hist(lengths[seq_len(n), 1], col = "#2369a8", border = "white",
       main = paste("R1: first", n, "pairs"), xlab = "Length (bases)")
  hist(lengths[seq_len(n), 2], col = "#d65f35", border = "white",
       main = paste("R2: first", n, "pairs"), xlab = "Length (bases)")
  summary
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3L || length(args) > 4L) {
    stop("usage: Rscript workflows/rnaseq/fastq_qc.R R1.fastq.gz R2.fastq.gz output_dir [max_pairs]")
  }
  limit <- if (length(args) == 4L) as.numeric(args[4]) else 10000L
  print(fastq_qc(args[1], args[2], args[3], limit))
}
