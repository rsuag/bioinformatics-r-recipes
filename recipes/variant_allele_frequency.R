# Task: summarize alternate allele frequency from genotype dosages.
# Input: numeric matrix with variants in rows, samples in columns; 0/1/2 or NA.
# Output: called samples, alternate alleles, frequency, and missing fraction.
# Assumes diploid, biallelic variants and alternate-allele dosages.

variant_allele_frequency <- function(genotypes) {
  if (!is.matrix(genotypes) || !is.numeric(genotypes) || !nrow(genotypes) ||
      !ncol(genotypes) || any(!is.na(genotypes) & !genotypes %in% 0:2)) {
    stop("genotypes must be a nonempty numeric matrix containing 0, 1, 2, or NA")
  }
  called <- rowSums(!is.na(genotypes))
  alternate <- rowSums(genotypes, na.rm = TRUE)
  ids <- rownames(genotypes)
  if (is.null(ids)) ids <- paste0("variant_", seq_len(nrow(genotypes)))
  data.frame(variant = ids, called_samples = called,
             alternate_alleles = alternate,
             alternate_frequency = ifelse(called > 0, alternate / (2 * called), NA_real_),
             missing_fraction = 1 - called / ncol(genotypes), row.names = NULL)
}

if (sys.nframe() == 0L) {
  genotypes <- matrix(c(0, 1, 2, NA, 1, 1, 2, 2), nrow = 2,
                      dimnames = list(c("v1", "v2"), NULL))
  print(variant_allele_frequency(genotypes))
}
