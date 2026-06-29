#' Flip GWAS Alleles and Adjust Beta Values
#'
#' This function compares GWAS allele information with a reference LD panel,
#' flips mismatched alleles, adjusts beta values and allele frequencies,
#' and optionally filters SNPs based on allele frequency differences.
#'
#' The LD information table can contain either \code{A1Freq} or \code{A2Freq}.
#' If \code{A1Freq} is unavailable but \code{A2Freq} is available,
#' \code{A1Freq} is derived as \code{1 - A2Freq}. The SNP identifier column can
#' be either \code{ID} or \code{snp}. The physical position column
#' \code{PhysPos} is used for ordering when available.
#'
#' @param gwas A data frame containing GWAS summary statistics. Must include
#'   columns \code{snp}, \code{A1}, \code{A2}, \code{beta}, and \code{AF1}.
#' @param info_LD A data frame containing reference LD information. Must include
#'   allele columns \code{A1} and \code{A2}, a SNP identifier column
#'   \code{ID} or \code{snp}, and either \code{A1Freq} or \code{A2Freq}.
#' @param no_AF Logical. If \code{TRUE}, skips filtering based on allele
#'   frequency difference.
#' @param af_threshold Numeric. Threshold for filtering SNPs based on absolute
#'   allele frequency difference. Default is \code{0.2}.
#'
#' @return A \code{data.table} with matched and flipped SNPs, adjusted beta
#'   values, allele frequencies, and LD-aligned alleles. Returns \code{NULL}
#'   if no SNPs are retained after allele matching or allele-frequency filtering.
#'
#' @examples
#' gwas <- data.frame(
#'   snp  = c("rs1", "rs2"),
#'   A1   = c("A", "G"),
#'   A2   = c("T", "C"),
#'   beta = c(0.2, -0.5),
#'   AF1  = c(0.1, 0.2),
#'   P    = c(0.01, 0.20),
#'   stringsAsFactors = FALSE
#' )
#'
#' ld_info <- data.frame(
#'   ID      = c("rs1", "rs2"),
#'   PhysPos = c(100, 200),
#'   A1      = c("A", "C"),
#'   A2      = c("T", "G"),
#'   A2Freq  = c(0.9, 0.2),
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
#' @importFrom data.table as.data.table
#' @export
allele_flip_gwas_beta <- function(gwas, info_LD, no_AF, af_threshold = 0.2) {
  gwas <- as.data.frame(gwas, stringsAsFactors = FALSE)
  info_LD <- as.data.frame(info_LD, stringsAsFactors = FALSE)

  find_col <- function(dat, candidates) {
    idx <- match(tolower(candidates), tolower(names(dat)))
    idx <- idx[!is.na(idx)]

    if (length(idx) == 0) {
      return(NA_character_)
    }

    names(dat)[idx[1]]
  }

  gwas_snp_col <- find_col(gwas, c("snp", "SNP", "ID", "id"))
  gwas_a1_col <- find_col(gwas, c("A1", "a1"))
  gwas_a2_col <- find_col(gwas, c("A2", "a2"))
  gwas_af1_col <- find_col(gwas, c("AF1", "af1", "A1Freq", "a1freq"))
  gwas_beta_col <- find_col(gwas, c("beta", "BETA"))

  required_gwas <- c(gwas_snp_col, gwas_a1_col, gwas_a2_col, gwas_af1_col, gwas_beta_col)

  if (anyNA(required_gwas)) {
    warning(
      "GWAS data must contain SNP, A1, A2, AF1, and beta columns.",
      call. = FALSE
    )
    return(NULL)
  }

  names(gwas)[names(gwas) == gwas_snp_col] <- "snp"
  names(gwas)[names(gwas) == gwas_a1_col] <- "A1"
  names(gwas)[names(gwas) == gwas_a2_col] <- "A2"
  names(gwas)[names(gwas) == gwas_af1_col] <- "AF1"
  names(gwas)[names(gwas) == gwas_beta_col] <- "beta"

  ld_snp_col <- find_col(info_LD, c("snp", "SNP", "ID", "id"))
  ld_a1_col <- find_col(info_LD, c("A1", "a1"))
  ld_a2_col <- find_col(info_LD, c("A2", "a2"))
  ld_a1freq_col <- find_col(info_LD, c("A1Freq", "a1freq", "AF1", "af1"))
  ld_a2freq_col <- find_col(info_LD, c("A2Freq", "a2freq", "AF2", "af2"))
  ld_pos_col <- find_col(info_LD, c("PhysPos", "physpos", "POS", "pos", "BP", "bp"))

  required_ld <- c(ld_snp_col, ld_a1_col, ld_a2_col)

  if (anyNA(required_ld)) {
    warning(
      "LD info must contain SNP identifier, A1, and A2 columns.",
      call. = FALSE
    )
    return(NULL)
  }

  names(info_LD)[names(info_LD) == ld_snp_col] <- "snp"
  names(info_LD)[names(info_LD) == ld_a1_col] <- "A1"
  names(info_LD)[names(info_LD) == ld_a2_col] <- "A2"

  if (!is.na(ld_pos_col)) {
    names(info_LD)[names(info_LD) == ld_pos_col] <- "PhysPos"
  } else {
    info_LD$PhysPos <- seq_len(nrow(info_LD))
  }

  if (!is.na(ld_a1freq_col)) {
    info_LD$A1Freq <- suppressWarnings(as.numeric(info_LD[[ld_a1freq_col]]))
  } else if (!is.na(ld_a2freq_col)) {
    info_LD$A1Freq <- 1 - suppressWarnings(as.numeric(info_LD[[ld_a2freq_col]]))
  } else {
    if (!no_AF) {
      warning(
        "LD info must contain either A1Freq or A2Freq when allele-frequency filtering is enabled.",
        call. = FALSE
      )
      return(NULL)
    }

    info_LD$A1Freq <- NA_real_
  }

  gwas$A1 <- as.character(gwas$A1)
  gwas$A2 <- as.character(gwas$A2)
  gwas$AF1 <- suppressWarnings(as.numeric(gwas$AF1))
  gwas$beta <- suppressWarnings(as.numeric(gwas$beta))

  info_LD$A1 <- as.character(info_LD$A1)
  info_LD$A2 <- as.character(info_LD$A2)

  gwas$.original_order <- seq_len(nrow(gwas))

  ld_sub <- info_LD[, c("snp", "A1", "A2", "A1Freq", "PhysPos"), drop = FALSE]
  names(ld_sub) <- c("snp", "LD_A1", "LD_A2", "A1Freq", "PhysPos")

  merged <- merge(gwas, ld_sub, by = "snp", all = FALSE, sort = FALSE)
  merged <- merged[order(merged$.original_order), , drop = FALSE]

  if (nrow(merged) == 0) {
    warning("No SNPs retained after matching GWAS and LD information.", call. = FALSE)
    return(NULL)
  }

  matched_id <- which(merged$A1 == merged$LD_A1 & merged$A2 == merged$LD_A2)
  mismatched_id <- which(merged$A1 == merged$LD_A2 & merged$A2 == merged$LD_A1)

  keep_id <- sort(c(matched_id, mismatched_id))

  if (length(keep_id) == 0) {
    warning("No SNPs retained after allele matching.", call. = FALSE)
    return(NULL)
  }

  merged <- merged[keep_id, , drop = FALSE]

  reverse_now <- merged$A1 == merged$LD_A2 & merged$A2 == merged$LD_A1

  if (any(reverse_now, na.rm = TRUE)) {
    merged$beta[reverse_now] <- -merged$beta[reverse_now]
    merged$AF1[reverse_now] <- 1 - merged$AF1[reverse_now]
  }

  merged$A1 <- merged$LD_A1
  merged$A2 <- merged$LD_A2

  merged <- merged[order(merged$PhysPos), , drop = FALSE]

  if (!no_AF) {
    merged$dif_AF <- abs(merged$AF1 - merged$A1Freq)
    merged <- merged[!is.na(merged$dif_AF) & merged$dif_AF < af_threshold, , drop = FALSE]
  }

  if (nrow(merged) == 0) {
    warning("No SNPs retained after allele-frequency filtering.", call. = FALSE)
    return(NULL)
  }

  merged$.original_order <- NULL
  merged$LD_A1 <- NULL
  merged$LD_A2 <- NULL

  data.table::as.data.table(merged)
}
