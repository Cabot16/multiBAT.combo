#' Compute mBAT Statistic via SVD for Two Populations
#'
#' This function performs the multivariate Burden Association Test (mBAT)
#' by decomposing the SNP correlation matrix via singular value decomposition (SVD),
#' selecting principal components up to a proportion of explained variance,
#' and computing the mBAT chi-square statistic and its p-value.
#'
#' @param gwas_noqc_final
#'   A data frame of GWAS summary statistics for the gene region.
#'   Must contain at least columns:
#'     snp: SNP identifier
#'     P:   SNP-level p-values
#' @param eig
#'   An object returned by \code{eigen()}, with components:
#'     values:  numeric vector of eigenvalues
#'     vectors: matrix of eigenvectors
#' @param prop
#'   Numeric in (0,1]. Proportion of total variance to retain.
#' @param z
#'   Numeric vector of SNP z-scores in the same order as rows of \code{gwas_noqc_final}.
#' @param chr_name
#'   Character. Chromosome name or number.
#' @param ID
#'   Character. Gene identifier.
#' @param gene_annotate
#'   A data frame with gene annotations. Columns:
#'     gene, chr, start, end
#' @param top1
#'   Character. SNP identifier of the most significant SNP in population 1.
#' @param top1_p
#'   Numeric. P-value of \code{top1}.
#' @param top2
#'   Character. SNP identifier of the most significant SNP in population 2.
#' @param top2_p
#'   Numeric. P-value of \code{top2}.
#'
#' @return A one-row data.frame with columns:
#'   Gene, Chr, Start, End, No.SNPs,
#'   SNP_start, SNP_end,
#'   TopSNP_pop1, TopSNP_Pvalue_pop1,
#'   TopSNP_pop2, TopSNP_Pvalue_pop2,
#'   No.Eigenvalues, Chisq_mBAT, P_mBAT.
#'   \code{No.SNPs} reports the number of unique SNPs retained for the gene.
#'
#' @importFrom magrittr    %>%
#' @importFrom dplyr       filter pull
#' @importFrom stats       pchisq
#' @importFrom utils       tail
#' @export
#'
#' @examples
#' # Simulate two SNPs with correlated effect sizes
#' gwas_noqc_final <- data.frame(
#'   snp = c("rs1", "rs2"),
#'   P   = c(0.01, 0.05),
#'   stringsAsFactors = FALSE
#' )
#' cov_mat <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
#' eig     <- eigen(cov_mat)
#' z       <- c(qnorm(1 - 0.01 / 2), qnorm(1 - 0.05 / 2))
#' gene_annotate <- data.frame(
#'   gene  = "GENE1",
#'   chr   = "1",
#'   start = 1000,
#'   end   = 2000,
#'   stringsAsFactors = FALSE
#' )
#' top1   <- "rs1"
#' top1_p <- 0.01
#' top2   <- "rs2"
#' top2_p <- 0.05
#' mbat_SVD(
#'   gwas_noqc_final = gwas_noqc_final,
#'   eig             = eig,
#'   prop            = 0.9,
#'   z               = z,
#'   chr_name        = "1",
#'   ID              = "GENE1",
#'   gene_annotate   = gene_annotate,
#'   top1            = top1,
#'   top1_p          = top1_p,
#'   top2            = top2,
#'   top2_p          = top2_p
#' )
mbat_SVD <- function(gwas_noqc_final,
                     eig,
                     prop,
                     z,
                     chr_name,
                     ID,
                     gene_annotate,
                     top1,
                     top1_p,
                     top2,
                     top2_p) {
  # 1. Filter out near-zero eigenvalues
  posval <- eig$values[eig$values > 1e-10]

  # 2. Select number of components to reach prop variance
  if (prop == 1) {
    k_prop <- length(posval)
  } else {
    cum_var <- cumsum(posval / sum(posval))
    k_prop <- which(cum_var > as.numeric(prop))[1]
  }

  # 3. Subset eigenvalues/vectors
  sel <- seq_len(k_prop)
  lambda_prop <- eig$values[sel]
  U_prop <- eig$vectors[, sel]

  # 4. Project z-scores and compute chi-square
  UPz_prop <- crossprod(U_prop, z)
  chisq_prop <- as.numeric(crossprod(UPz_prop, (1 / lambda_prop) * UPz_prop))
  P_mBAT <- stats::pchisq(chisq_prop, df = k_prop, lower.tail = FALSE)

  # 5. Extract gene boundaries
  anno <- gene_annotate %>%
    dplyr::filter(gene == ID, chr == chr_name) %>%
    dplyr::slice(1)

  # 6. Build result row with both populations' Top SNPs
  res <- data.frame(
    Gene                = ID,
    Chr                 = chr_name,
    Start               = as.numeric(anno$start),
    End                 = as.numeric(anno$end),
    No.SNPs             = length(unique(gwas_noqc_final$snp)),
    SNP_start           = gwas_noqc_final$snp[1],
    SNP_end             = utils::tail(gwas_noqc_final$snp, 1),
    TopSNP_pop1         = top1,
    TopSNP_Pvalue_pop1  = top1_p,
    TopSNP_pop2         = top2,
    TopSNP_Pvalue_pop2  = top2_p,
    No.Eigenvalues      = k_prop,
    Chisq_mBAT          = chisq_prop,
    P_mBAT              = P_mBAT,
    stringsAsFactors    = FALSE
  )

  return(res)
}