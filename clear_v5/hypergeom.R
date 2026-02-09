hypergeom_test <- function(term_genes, all_pvals, threshold = 0.05) {
  total_genes <- length(all_pvals)

  significant_genes <- names(all_pvals)[all_pvals < threshold]
  num_significant <- length(significant_genes)

  genes_in_term <- intersect(term_genes, names(all_pvals))
  num_genes_in_term <- length(genes_in_term)

  genes_in_overlap <- intersect(genes_in_term, significant_genes)
  num_genes_in_overlap <- length(genes_in_overlap)

  # Perform the hypergeometric test
  p_value <- phyper(
    q = num_genes_in_overlap - 1,  # Number of successes in the sample - 1 (phyper is P(X <= q), so subtract 1 for P(X >= q))
    m = num_genes_in_term,        # Total successes in population (genes in term)
    n = total_genes - num_genes_in_term,  # Total failures in population
    k = num_significant,          # Sample size (significant genes)
    lower.tail = FALSE            # We want P(X >= q)
  )

  return(list(
    term_size = num_genes_in_term,
    overlap_size = num_genes_in_overlap,
    p_value = p_value
  ))
}
