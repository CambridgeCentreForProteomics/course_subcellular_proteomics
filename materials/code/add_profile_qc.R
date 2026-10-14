## addProfileQC()
##
## Records two indicators of protein profile quality in a QFeatures object
## and carries them through to the joined (concatenated) set:
##
##   - peptide support: the `.n` column added by aggregateFeatures()
##   - missingness: the proportion of precursor values that were missing
##     before imputation, calculated from the per-precursor `nNA` column
##
## Both are lost by joinAssays(), because they differ between replicates.
##
## Arguments
##   qf          a QFeatures object
##   prec_sets   named character vector of precursor-level sets that carry `nNA`
##   prot_sets   named character vector of protein-level sets that carry `.n`
##               (same names and order as prec_sets, e.g. rep1, rep2, rep3)
##   joined_set  name of the joined protein-level set
##   group_col   rowData column linking precursors to proteins
##
## Adds to each protein-level set:
##   prop_missing
## Adds to the joined set:
##   n_<rep>, n_min, n_mean, missing_<rep>, missing_max, missing_mean
##
## Returns the updated QFeatures object.

addProfileQC <- function(qf,
                         prec_sets,
                         prot_sets,
                         joined_set = "proteins_joined",
                         group_col = "Protein.Group") {

  stopifnot(identical(names(prec_sets), names(prot_sets)))
  reps <- names(prot_sets)
  joined_ids <- rownames(qf[[joined_set]])

  for (r in reps) {
    prec <- prec_sets[[r]]
    prot <- prot_sets[[r]]
    rd <- rowData(qf[[prec]])

    ## Proportion of missing values per protein:
    ## total missing values / (number of precursors x number of fractions)
    tot_na <- tapply(rd$nNA, rd[[group_col]], sum)
    tot_values <- tapply(rd$nNA, rd[[group_col]], length) * ncol(qf[[prec]])
    prop <- tot_na / tot_values

    ## Add to the protein-level set of this replicate
    rowData(qf[[prot]])$prop_missing <-
      as.numeric(prop[rownames(qf[[prot]])])
  }

  ## Extract a rowData column from each replicate, matched to the joined set
  per_rep <- function(col, prefix) {
    m <- sapply(prot_sets, function(s) {
      rd <- rowData(qf[[s]])
      rd[[col]][match(joined_ids, rownames(rd))]
    })
    m <- matrix(m, ncol = length(prot_sets))
    colnames(m) <- paste0(prefix, reps)
    m
  }

  n_mat <- per_rep(".n", "n_")
  miss_mat <- per_rep("prop_missing", "missing_")

  ## Per-replicate values plus summaries across replicates
  qc_df <- DataFrame(n_mat,
                     n_min = apply(n_mat, 1, min),
                     n_mean = rowMeans(n_mat),
                     miss_mat,
                     missing_max = apply(miss_mat, 1, max),
                     missing_mean = rowMeans(miss_mat))

  rowData(qf[[joined_set]]) <- cbind(rowData(qf[[joined_set]]), qc_df)

  qf
}
