# Task: read a DNA FASTA file and calculate per-sequence GC fractions.
# Input: FASTA with unique headers and A/C/G/T/N bases (multi-line sequences allowed).
# Output: sequence ID, length, called bases, and GC fraction among called bases.

fasta_gc <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) stop("FASTA file is empty")
  headers <- grep("^>", lines)
  if (!length(headers) || headers[1] != 1L) stop("FASTA must start with a header")
  ids <- sub("^>([^[:space:]]+).*$", "\\1", lines[headers])
  if (any(ids == ">" | !nzchar(ids)) || anyDuplicated(ids)) stop("FASTA IDs must be nonempty and unique")
  ends <- c(headers[-1] - 1L, length(lines))
  sequences <- vapply(seq_along(headers), function(i) {
    if (ends[i] <= headers[i]) return("")
    paste0(gsub("[[:space:]]", "", lines[seq.int(headers[i] + 1L, ends[i])]), collapse = "")
  }, character(1))
  sequences <- toupper(sequences)
  if (any(!nzchar(sequences)) || any(grepl("[^ACGTN]", sequences))) {
    stop("each sequence must contain only A, C, G, T, or N")
  }
  bases <- strsplit(sequences, "", fixed = TRUE)
  called <- vapply(bases, function(x) sum(x != "N"), integer(1))
  gc <- vapply(bases, function(x) sum(x %in% c("G", "C")), integer(1))
  data.frame(id = ids, length = nchar(sequences), called_bases = called,
             gc_fraction = ifelse(called > 0, gc / called, NA_real_), row.names = NULL)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  path <- if (length(args)) args[1] else "data/toy_sequences.fasta"
  print(fasta_gc(path), row.names = FALSE)
}
