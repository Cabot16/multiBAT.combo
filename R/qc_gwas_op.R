#' Subset GWAS Summary Statistics by QC-Passed SNPs
#'
#' This function filters a GWAS summary statistics data frame to retain only SNPs
#' that passed quality control. The returned data frame preserves the order of
#' the original GWAS input.
#'
#' @param gwas A data frame of GWAS summary statistics. Must include a column
#'   `snp` with SNP identifiers.
#' @param qc_snp_list A data frame of QC-passed SNPs. Must include a column
#'   `snp` with SNP identifiers.
#'
#' @return A subset of `gwas` containing only rows where the `snp` appears in
#'   `qc_snp_list`. The row order matches the order of SNPs in `gwas`.
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
#' # QC-passed SNP list
#' qc_list <- data.frame(
#'   snp = c("rs2", "rs4"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Filter GWAS by QC-passed SNPs
#' qc_gwas_op(gwas_df, qc_list)
qc_gwas_op <- function(gwas, qc_snp_list) {
  # find common SNPs between GWAS and QC list
  snplist <- Reduce(intersect, list(gwas$snp, qc_snp_list$snp))
  
  # subset GWAS to those SNPs, preserving original order
  gwas_qced <- gwas[match(snplist, gwas$snp), ]
  
  return(gwas_qced)
}
