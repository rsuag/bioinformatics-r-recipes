# Task: estimate first-order transitions between discrete sequence states.
# Input: state vector and optional fixed state alphabet.
# Output: transition counts and row-normalized probabilities.
# A zero row means that state has no observed outgoing transition.

markov_transitions <- function(sequence, states = sort(unique(sequence))) {
  if (!is.character(sequence) || length(sequence) < 2L || anyNA(sequence) ||
      any(!nzchar(sequence)) || !is.character(states) || !length(states) ||
      anyNA(states) || any(!nzchar(states)) || anyDuplicated(states) ||
      !all(sequence %in% states)) {
    stop("sequence needs >=2 known nonempty states and states must be unique")
  }
  counts <- table(factor(head(sequence, -1L), levels = states),
                  factor(tail(sequence, -1L), levels = states))
  dimnames(counts) <- list(from = states, to = states)
  totals <- rowSums(counts)
  probabilities <- sweep(counts, 1L, ifelse(totals > 0, totals, 1), "/")
  list(counts = counts, probabilities = probabilities)
}

if (sys.nframe() == 0L) {
  print(markov_transitions(c("A", "C", "A", "A", "G"))$probabilities)
}
