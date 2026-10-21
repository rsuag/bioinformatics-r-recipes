# Task: count overlapping DNA words in one sequence.
# Input: one A/C/G/T/N sequence and a positive integer k.
# Output: counts of observed k-mers; windows containing N are skipped.

kmer_frequencies <- function(sequence, k = 3L) {
  if (!is.character(sequence) || length(sequence) != 1L || is.na(sequence) ||
      !nzchar(sequence) || grepl("[^ACGTN]", toupper(sequence))) {
    stop("sequence must be one nonempty DNA string containing A, C, G, T, or N")
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || !is.finite(k) ||
      k < 1 || k != floor(k)) stop("k must be a positive integer")
  sequence <- toupper(sequence)
  if (k > nchar(sequence)) return(setNames(integer(), character()))
  starts <- seq_len(nchar(sequence) - k + 1L)
  words <- substring(sequence, starts, starts + k - 1L)
  words <- words[!grepl("N", words, fixed = TRUE)]
  sort(table(words), decreasing = TRUE)
}

if (sys.nframe() == 0L) print(kmer_frequencies("ACGTACGNAC", 2L))
