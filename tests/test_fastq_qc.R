source("workflows/rnaseq/fastq_qc.R")
root <- tempfile()
dir.create(root)
make_fastq <- function(path, mate, ids = 1:3) {
  connection <- gzfile(path, "wt")
  on.exit(close(connection))
  for (i in ids) writeLines(c(paste0("@read_", i, "/", mate), "ACGTN", "+", "IIIII"),
                           connection)
}
r1 <- file.path(root, "R1.fastq.gz")
r2 <- file.path(root, "R2.fastq.gz")
make_fastq(r1, 1)
make_fastq(r2, 2)
report <- fastq_qc(r1, r2, file.path(root, "output"))
stopifnot(all(report$pairs_checked == 3L), all(report$median_length == 5),
          all(report$mean_phred == 40), all(report$mean_gc == 0.5),
          file.exists(file.path(root, "output", "read_lengths.png")))
make_fastq(r2, 2, ids = 2:4)
stopifnot(inherits(try(fastq_qc(r1, r2, file.path(root, "bad")), silent = TRUE), "try-error"))
unlink(root, recursive = TRUE)
cat("FASTQ QC tests passed.\n")
