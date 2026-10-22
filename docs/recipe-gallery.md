# Recipe result gallery

Every standalone recipe has a result figure generated from synthetic or bundled example data. Rebuild all figures from the repository root with:

```bash
Rscript scripts/render_recipe_figures.R
```

The figures show the result type and expected shape; they are not biological findings.

| Recipe | Example result |
|---|---|
| [Adjust p-values](../recipes/adjust_pvalues.R) | ![P-values and adjusted q-values](figures/recipes/adjust_pvalues.png) |
| [Beta-binomial update](../recipes/beta_binomial.R) | ![Posterior density and credible interval](figures/recipes/beta_binomial.png) |
| [Bootstrap mean](../recipes/bootstrap_mean.R) | ![Observed measurements and bootstrap interval](figures/recipes/bootstrap_mean.png) |
| [Cluster samples](../recipes/cluster_samples.R) | ![Hierarchical sample clustering](figures/recipes/cluster_samples.png) |
| [Differential expression](../recipes/differential_expression.R) | ![Normalized-expression result plot](figures/recipes/differential_expression.png) |
| [Explore counts](../recipes/explore_counts.R) | ![PCA of count-matrix samples](figures/recipes/explore_counts.png) |
| [FASTA GC](../recipes/fasta_gc.R) | ![GC fraction by sequence](figures/recipes/fasta_gc.png) |
| [Gene-set enrichment](../recipes/gene_set_enrichment.R) | ![Gene-set enrichment scores](figures/recipes/gene_set_enrichment.png) |
| [K-mer frequencies](../recipes/kmer_frequencies.R) | ![DNA k-mer frequencies](figures/recipes/kmer_frequencies.png) |
| [Logistic classifier](../recipes/logistic_classifier.R) | ![Held-out prediction probabilities](figures/recipes/logistic_classifier.png) |
| [Markov transitions](../recipes/markov_transitions.R) | ![State transition probability matrix](figures/recipes/markov_transitions.png) |
| [Metabolomics fold change](../recipes/metabolomics_fold_change.R) | ![Metabolite log2 differences](figures/recipes/metabolomics_fold_change.png) |
| [Metabolomics missingness](../recipes/metabolomics_missingness.R) | ![Missing-value map](figures/recipes/metabolomics_missingness.png) |
| [Variant allele frequency](../recipes/variant_allele_frequency.R) | ![Allele frequency and missingness](figures/recipes/variant_allele_frequency.png) |
