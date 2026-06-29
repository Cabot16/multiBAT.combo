#' Flip GWAS Alleles and Adjust Beta Values
#'
#' This function compares GWAS allele information with a reference LD panel,
#' flips mismatched alleles, adjusts beta values and allele frequencies,
#' and optionally filters SNPs based on allele frequency differences.
#'
#' @param gwas A data frame containing GWAS summary statistics. Must include columns:
#'   \code{snp}, \code{A1}, \code{A2}, \code{beta}, and \code{AF1}.
#' @param info_LD A data frame containing reference LD information. Must include columns:
#'   \code{ID}, \code{PhysPos}, \code{A1}, \code{A2}, and \code{A1Freq}.
#' @param no_AF Logical. If \code{TRUE}, skips filtering based on allele frequency difference.
#' @param af_threshold Numeric. Threshold for filtering SNPs based on absolute allele
#'   frequency difference. Default is \code{0.2}.
#'
#' @return A \code{data.table} with matched and flipped SNPs, adjusted beta values,
#'   allele frequencies, and LD-aligned alleles. Returns \code{NULL} if no SNPs
#'   are retained after allele matching or allele-frequency filtering.
#'
#' @examples
#' gwas <- data.frame(
#'   snp  = c("rs1", "rs2"),
#'   A1   = c("A", "G"),
#'   A2   = c("T", "C"),
#'   beta = c(0.2, -0.5),
#'   AF1  = c(0.1, 0.2),
#'   stringsAsFactors = FALSE
#' )
#'
#' ld_info <- data.frame(
#'   ID      = c("rs1", "rs2"),
#'   PhysPos = c(100, 200),
#'   A1      = c("A", "C"),
#'   A2      = c("T", "G"),
#'   A1Freq  = c(0.1, 0.8),
#'   stringsAsFactors = FALSE
#' )
#'
#' result <- allele_flip_gwas_beta(
#'   gwas         = gwas,
#'   info_LD      = ld_info,
#'   no_AF        = FALSE,
#'   af_threshold = 0.2
#' )
#' print(result)
#'
#' @import dplyr
#' @importFrom data.table as.data.table
#' @export
allele_flip_gwas_beta <- function(gwas, info_LD, no_AF, af_threshold = 0.2) {
  allele <- c("A1", "A2")

  gwas <- gwas %>% dplyr::mutate_at(allele, as.character)
  gwas$AF1 <- as.numeric(gwas$AF1)

  info_LD <- info_LD %>% dplyr::mutate_at(allele, as.character)
  info_LD$A1Freq <- as.numeric(info_LD$A1Freq)

  matched_id <- which(gwas$A1 == info_LD$A1 & gwas$A2 == info_LD$A2)
  gwas_tmp1 <- gwas[matched_id, ]

  mismatched_id <- which(gwas$A1 == info_LD$A2 & gwas$A2 == info_LD$A1)

  if (length(mismatched_id) > 0) {
    gwas$beta[mismatched_id] <- -gwas$beta[mismatched_id]
    gwas$AF1[mismatched_id] <- 1 - gwas$AF1[mismatched_id]
    gwas$A1[mismatched_id] <- info_LD$A1[mismatched_id]
    gwas$A2[mismatched_id] <- info_LD$A2[mismatched_id]
  }

  gwas_tmp2 <- gwas[mismatched_id, ]
  gwas_final <- rbind(gwas_tmp1, gwas_tmp2)

  if (nrow(gwas_final) == 0) {
    warning("No SNPs retained after allele matching.", call. = FALSE)
    return(NULL)
  }

  joined <- info_LD %>%
    dplyr::select(ID, PhysPos, A1Freq) %>%
    dplyr::rename(snp = ID)

  gwas_final <- dplyr::inner_join(gwas_final, joined, by = "snp") %>%
    dplyr::arrange(PhysPos)

  if (!no_AF) {
    gwas_final$dif_AF <- abs(gwas_final$AF1 - gwas_final$A1Freq)
    gwas_final <- dplyr::filter(gwas_final, dif_AF < af_threshold)
  }

  if (nrow(gwas_final) == 0) {
    warning("No SNPs retained after allele-frequency filtering.", call. = FALSE)
    return(NULL)
  }

  return(data.table::as.data.table(gwas_final))
}