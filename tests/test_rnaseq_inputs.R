source("workflows/rnaseq/run.R")
base <- "workflows/rnaseq/example"
input <- read_rnaseq_inputs(file.path(base, "counts.csv"),
                            file.path(base, "samples.csv"), "control", "treated")
stopifnot(ncol(input$counts) == 6L, nrow(input$counts) == 120L,
          identical(colnames(input$counts), rownames(input$samples)),
          identical(levels(input$samples$condition), c("control", "treated")))

samples <- input$samples[6:1, c("sample_id", "condition")]
sample_path <- tempfile(fileext = ".csv")
write.csv(samples, sample_path, row.names = FALSE)
reordered <- read_rnaseq_inputs(file.path(base, "counts.csv"), sample_path,
                                "control", "treated")
stopifnot(identical(rownames(reordered$samples), colnames(reordered$counts)))

samples$sample_id[1] <- "wrong_id"
write.csv(samples, sample_path, row.names = FALSE)
stopifnot(inherits(try(read_rnaseq_inputs(file.path(base, "counts.csv"), sample_path,
                                         "control", "treated"), silent = TRUE), "try-error"))
unlink(sample_path)

qc_dir <- tempfile()
dir.create(qc_dir)
write_qc(input$counts, input$samples, qc_dir)
stopifnot(file.exists(file.path(qc_dir, "library_sizes.csv")),
          file.exists(file.path(qc_dir, "pca.png")),
          nrow(read.csv(file.path(qc_dir, "pca_coordinates.csv"))) == 6L)
unlink(qc_dir, recursive = TRUE)
cat("RNA-seq input and QC tests passed.\n")
