#' Compute mBAT Statistic via SVD for Two Populations
#'
#' Computes the multivariate Burden Association Test (mBAT) statistic for a
#' joint two-population analysis. The function uses the eigendecomposition of
#' the joint SNP correlation matrix, retains principal components explaining
#' the specified proportion of variance, and calculates the corresponding
#' chi-square statistic and P value.
#'
#' In the multiBAT workflow, the joint correlation matrix is constructed as a
#' block-diagonal matrix from the population-specific LD matrices, and the
#' supplied Z-score vector contains the Z scores from both populations.
#'
#' @param gwas_noqc_final A data frame or data.table representing the final
#'   joint GWAS input used for the test. The number of rows must equal the
#'   length of \code{z}.
#' @param eig An object returned by \code{eigen()}, containing the eigenvalues
#'   and eigenvectors of the joint LD matrix.
#' @param prop Numeric value in \code{(0, 1]}. Proportion of eigen-variance to
#'   retain.
#' @param z Numeric vector containing the joint SNP Z scores.
#' @param chr_name Character value giving the chromosome name or number.
#' @param ID Character value giving the gene identifier.
#' @param gene_annotate A data frame containing gene annotations with columns
#'   \code{gene}, \code{chr}, \code{start}, and \code{end}.
#' @param top1 Character value giving the most significant SNP in population 1.
#' @param top1_p Numeric P value for \code{top1}.
#' @param top2 Character value giving the most significant SNP in population 2.
#' @param top2_p Numeric P value for \code{top2}.
#' @param no_snps_pop1 Number of SNPs retained for population 1 after
#'   population-specific allele alignment and allele-frequency filtering.
#' @param no_snps_pop2 Number of SNPs retained for population 2 after
#'   population-specific allele alignment and allele-frequency filtering.
#' @param no_snps_common Number of SNPs shared between the two populations and
#'   retained for the final joint multiBAT analysis.
#' @param common_snp_positions Numeric vector containing the genomic positions
#'   of the SNPs retained for the final joint multiBAT analysis.
#'
#' @return A one-row data.frame containing:
#'   \code{Gene}, \code{Chr}, \code{Start}, \code{End},
#'   \code{No.SNPs_pop1}, \code{No.SNPs_pop2}, \code{No.SNPs_common},
#'   \code{SNP_start}, \code{SNP_end},
#'   \code{TopSNP_pop1}, \code{TopSNP_Pvalue_pop1},
#'   \code{TopSNP_pop2}, \code{TopSNP_Pvalue_pop2},
#'   \code{No.Eigenvalues}, \code{Chisq_mBAT}, and \code{P_mBAT}.
#'
#' @importFrom stats pchisq
#' @export
#'
#' @examples
#' gwas_noqc_final <- data.frame(
#'   snp = c("rs1", "rs2", "rs1", "rs2"),
#'   P = c(0.01, 0.05, 0.02, 0.04),
#'   stringsAsFactors = FALSE
#' )
#'
#' cov_mat <- matrix(
#'   c(
#'     1.0, 0.3, 0.0, 0.0,
#'     0.3, 1.0, 0.0, 0.0,
#'     0.0, 0.0, 1.0, 0.2,
#'     0.0, 0.0, 0.2, 1.0
#'   ),
#'   nrow = 4,
#'   byrow = TRUE
#' )
#'
#' eig <- eigen(cov_mat)
#'
#' z <- c(
#'   2.0,
#'   1.5,
#'   1.8,
#'   1.6
#' )
#'
#' gene_annotate <- data.frame(
#'   gene = "GENE1",
#'   chr = "1",
#'   start = 1000,
#'   end = 2000,
#'   stringsAsFactors = FALSE
#' )
#'
#' mbat_SVD(
#'   gwas_noqc_final = gwas_noqc_final,
#'   eig = eig,
#'   prop = 0.9,
#'   z = z,
#'   chr_name = "1",
#'   ID = "GENE1",
#'   gene_annotate = gene_annotate,
#'   top1 = "rs1",
#'   top1_p = 0.01,
#'   top2 = "rs1",
#'   top2_p = 0.02,
#'   no_snps_pop1 = 2,
#'   no_snps_pop2 = 2,
#'   no_snps_common = 2,
#'   common_snp_positions = c(1200, 1800)
#' )
mbat_SVD <- function(
    gwas_noqc_final,
    eig,
    prop,
    z,
    chr_name,
    ID,
    gene_annotate,
    top1,
    top1_p,
    top2,
    top2_p,
    no_snps_pop1,
    no_snps_pop2,
    no_snps_common,
    common_snp_positions
) {

    # ------------------------------------------------------------
    # 1. Check consistency between the joint GWAS input and Z scores
    # ------------------------------------------------------------

    if (nrow(gwas_noqc_final) != length(z)) {
        stop(
            paste0(
                "The number of joint GWAS rows does not match the length ",
                "of the Z-score vector."
            ),
            call. = FALSE
        )
    }

    if (
        nrow(eig$vectors) != length(z) ||
        length(eig$values) != ncol(eig$vectors)
    ) {
        stop(
            paste0(
                "The eigendecomposition dimensions do not match the ",
                "joint Z-score vector."
            ),
            call. = FALSE
        )
    }


    # ------------------------------------------------------------
    # 2. Retain positive eigenvalues
    # ------------------------------------------------------------

    posval <- eig$values[
        eig$values > 1e-10
    ]

    if (length(posval) == 0L) {
        stop(
            "No positive LD eigenvalues remained after numerical filtering.",
            call. = FALSE
        )
    }


    # ------------------------------------------------------------
    # 3. Select components explaining the requested variance
    # ------------------------------------------------------------

    if (prop == 1) {

        k_prop <- length(
            posval
        )

    } else {

        cum_var <- cumsum(
            posval /
                sum(posval)
        )

        k_prop <- which(
            cum_var > as.numeric(prop)
        )[1]
    }


    # ------------------------------------------------------------
    # 4. Subset eigenvalues and eigenvectors
    # ------------------------------------------------------------

    sel <- seq_len(
        k_prop
    )

    lambda_prop <- eig$values[
        sel
    ]

    U_prop <- eig$vectors[
        ,
        sel,
        drop = FALSE
    ]


    # ------------------------------------------------------------
    # 5. Calculate the joint mBAT statistic
    # ------------------------------------------------------------

    UPz_prop <- crossprod(
        U_prop,
        z
    )

    chisq_prop <- as.numeric(
        crossprod(
            UPz_prop,
            (1 / lambda_prop) *
                UPz_prop
        )
    )

    P_mBAT <- stats::pchisq(
        chisq_prop,
        df = k_prop,
        lower.tail = FALSE
    )


    # ------------------------------------------------------------
    # 6. Extract gene annotation
    # ------------------------------------------------------------

    anno_index <- which(
        as.character(gene_annotate[["gene"]]) ==
            as.character(ID) &
        as.character(gene_annotate[["chr"]]) ==
            as.character(chr_name)
    )

    if (length(anno_index) == 0L) {
        stop(
            "No matching gene annotation record was found.",
            call. = FALSE
        )
    }

    anno <- gene_annotate[
        anno_index[1],
        ,
        drop = FALSE
    ]


    # ------------------------------------------------------------
    # 7. Determine the genomic range of the analysed SNPs
    # ------------------------------------------------------------

    common_snp_positions <- suppressWarnings(
        as.numeric(
            common_snp_positions
        )
    )

    valid_positions <- common_snp_positions[
        is.finite(common_snp_positions)
    ]

    if (length(valid_positions) == 0L) {
        stop(
            "No valid genomic positions were available for the analysed SNPs.",
            call. = FALSE
        )
    }

    snp_start <- min(
        valid_positions
    )

    snp_end <- max(
        valid_positions
    )


    # ------------------------------------------------------------
    # 8. Build the result row
    # ------------------------------------------------------------

    res <- data.frame(
        Gene = ID,
        Chr = chr_name,
        Start = as.numeric(
            anno$start
        ),
        End = as.numeric(
            anno$end
        ),
        No.SNPs_pop1 = as.integer(
            no_snps_pop1
        ),
        No.SNPs_pop2 = as.integer(
            no_snps_pop2
        ),
        No.SNPs_common = as.integer(
            no_snps_common
        ),
        SNP_start = snp_start,
        SNP_end = snp_end,
        TopSNP_pop1 = top1,
        TopSNP_Pvalue_pop1 = top1_p,
        TopSNP_pop2 = top2,
        TopSNP_Pvalue_pop2 = top2_p,
        No.Eigenvalues = k_prop,
        Chisq_mBAT = chisq_prop,
        P_mBAT = P_mBAT,
        stringsAsFactors = FALSE
    )


    # ------------------------------------------------------------
    # 9. Return the result
    # ------------------------------------------------------------

    return(res)
}