#' Generate LD Matrix After Quality Control
#'
#' This function filters and reorders a reference LD correlation matrix to
#' include only SNPs that are present in both a GWAS summary and a QC-passed
#' SNP list. The output matrix will have rows and columns matching the order
#' of SNPs in the GWAS data.
#'
#' @param gwas A data frame of GWAS summary statistics. Must contain a column
#'   `snp` with SNP identifiers.
#' @param LD A square numeric matrix of pairwise LD correlations. Row names
#'   and column names must be SNP identifiers matching those in `gwas$snp`.
#' @param qc_snp_list A data frame or vector containing SNPs that passed QC.
#'   If a data frame, it must have a column `snp` with SNP identifiers.
#'
#' @return A numeric matrix of LD correlations filtered to include only
#'   SNPs common to `gwas$snp` and `qc_snp_list$snp`, and reordered to the
#'   same sequence as in `gwas$snp`.
#'
#' @export
#'
#' @examples
#' # Example GWAS summary
#' gwas_df <- data.frame(
#'   snp   = c("rs1", "rs2", "rs3", "rs4"),
#'   beta  = rnorm(4),
#'   P     = runif(4),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Example QC-passed SNP list
#' qc_list <- data.frame(
#'   snp = c("rs2", "rs4", "rs5"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Simulated LD matrix for SNPs rs1–rs5
#' snps_all <- paste0("rs", 1:5)
#' LD_mat <- diag(5)
#' rownames(LD_mat) <- snps_all
#' colnames(LD_mat) <- snps_all
#'
#' # Filter and reorder the LD matrix after QC
#' qc_LD_op(gwas_df, LD_mat, qc_list)
qc_LD_op <- function(gwas, LD, qc_snp_list) {
  # Identify SNPs common to GWAS and QC list
  snplist <- Reduce(intersect, list(gwas$snp,
                                     if (is.data.frame(qc_snp_list)) qc_snp_list$snp else qc_snp_list))
  
  # Determine the indices (column names) in LD matching the filtered SNP list
  LD_qced_snp_id <- colnames(LD)[match(snplist, colnames(LD))]
  
  # Subset and reorder the LD matrix to those SNPs
  LD_qced <- LD[LD_qced_snp_id, LD_qced_snp_id, drop = FALSE]
  
  return(LD_qced)
}