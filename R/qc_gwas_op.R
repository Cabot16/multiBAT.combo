#' Subset GWAS Summary Statistics by a SNP List
#'
#' Filters GWAS summary statistics to retain only SNPs present in a supplied
#' SNP list. The order of retained rows follows the original GWAS input.
#'
#' This function performs SNP subsetting only. It does not perform numeric
#' quality control of GWAS fields such as beta, standard error, P value,
#' allele frequency, or sample size. Numeric quality control is handled
#' separately before gene-level analysis.
#'
#' @param gwas A data.frame or data.table containing GWAS summary statistics.
#'   Must contain a column named \code{snp}.
#' @param qc_snp_list A data.frame or data.table containing SNPs to retain.
#'   Must contain a column named \code{snp}.
#'
#' @return A data.table containing GWAS rows whose SNP identifiers occur in
#'   \code{qc_snp_list}. The retained rows remain in the same order as in the
#'   original GWAS input.
#'
#' @importFrom data.table as.data.table copy
#' @export
#'
#' @examples
#' gwas_df <- data.frame(
#'   snp = c("rs1", "rs2", "rs3", "rs4"),
#'   beta = c(0.2, -0.1, -0.3, 0.5),
#'   se = c(0.1, 0.2, 0.2, 0.3),
#'   P = c(0.05, 0.01, 0.20, 0.30),
#'   stringsAsFactors = FALSE
#' )
#'
#' qc_list <- data.frame(
#'   snp = c("rs1", "rs3"),
#'   stringsAsFactors = FALSE
#' )
#'
#' qc_gwas_op(
#'   gwas = gwas_df,
#'   qc_snp_list = qc_list
#' )
qc_gwas_op <- function(
    gwas,
    qc_snp_list
) {

    # ------------------------------------------------------------
    # 1. Copy input objects
    # ------------------------------------------------------------

    gwas <- data.table::copy(
        data.table::as.data.table(gwas)
    )

    qc_snp_list <- data.table::copy(
        data.table::as.data.table(qc_snp_list)
    )


    # ------------------------------------------------------------
    # 2. Check required columns
    # ------------------------------------------------------------

    if (!"snp" %in% names(gwas)) {
        stop(
            "GWAS data must contain a 'snp' column.",
            call. = FALSE
        )
    }

    if (!"snp" %in% names(qc_snp_list)) {
        stop(
            "The SNP list must contain a 'snp' column.",
            call. = FALSE
        )
    }


    # ------------------------------------------------------------
    # 3. Standardise SNP identifiers as character values
    # ------------------------------------------------------------

    gwas[["snp"]] <- as.character(
        gwas[["snp"]]
    )

    qc_snp_list[["snp"]] <- as.character(
        qc_snp_list[["snp"]]
    )


    # ------------------------------------------------------------
    # 4. Extract valid SNP identifiers from the supplied SNP list
    # ------------------------------------------------------------

    qc_snps <- unique(
        qc_snp_list[["snp"]][
            !is.na(qc_snp_list[["snp"]]) &
            nzchar(qc_snp_list[["snp"]])
        ]
    )


    # ------------------------------------------------------------
    # 5. Return an empty GWAS table if the SNP list is empty
    # ------------------------------------------------------------

    if (length(qc_snps) == 0L) {
        return(
            gwas[0L]
        )
    }


    # ------------------------------------------------------------
    # 6. Retain GWAS rows matching the supplied SNP list
    # ------------------------------------------------------------

    keep_snp <- (
        !is.na(gwas[["snp"]]) &
        nzchar(gwas[["snp"]]) &
        gwas[["snp"]] %in% qc_snps
    )

    gwas_qced <- gwas[
        keep_snp
    ]


    # ------------------------------------------------------------
    # 7. Return the subsetted GWAS data
    # ------------------------------------------------------------

    gwas_qced
}