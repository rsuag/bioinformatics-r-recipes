# Task: test gene-set overrepresentation in a selected gene list.
# Input: selected IDs, named list of gene sets, and measured-gene universe.
# Output: hypergeometric upper-tail p-values and BH q-values.
# The universe must be all genes eligible for selection, not an entire genome by default.

gene_set_enrichment <- function(selected, gene_sets, universe) {
  if (!is.character(universe) || !length(universe) || anyNA(universe) ||
      any(!nzchar(universe)) || anyDuplicated(universe) ||
      !is.character(selected) || anyNA(selected) || anyDuplicated(selected) ||
      !all(selected %in% universe) || !is.list(gene_sets) || !length(gene_sets) ||
      is.null(names(gene_sets)) || any(!nzchar(names(gene_sets))) ||
      anyDuplicated(names(gene_sets)) ||
      !all(vapply(gene_sets, is.character, logical(1)))) {
    stop("use unique gene IDs, named character gene sets, and selected IDs inside universe")
  }
  sets <- lapply(gene_sets, intersect, y = universe)
  overlap <- vapply(sets, function(s) sum(selected %in% s), integer(1))
  set_size <- lengths(sets)
  p <- phyper(overlap - 1L, set_size, length(universe) - set_size,
              length(selected), lower.tail = FALSE)
  data.frame(set = names(sets), set_size = set_size, overlap = overlap,
             p_value = p, q_value = p.adjust(p, "BH"), row.names = NULL)
}

if (sys.nframe() == 0L) {
  print(gene_set_enrichment(c("g1", "g2"),
        list(pathway_a = c("g1", "g2", "g3"), pathway_b = c("g4", "g5")),
        paste0("g", 1:10)))
}
