# ============================================================
# allele_flip_gwas_beta.R
#
# Harmonise GWAS alleles with LD reference information.
#
# GWAS numeric QC is handled separately.
# ============================================================


#' Harmonise GWAS Alleles with LD Reference Information
#'
#' Harmonises GWAS alleles with population-specific LD reference
#' information. When the GWAS alleles are reversed relative to the LD
#' reference, the effect estimate is multiplied by -1 and the effect-allele
#' frequency is transformed as \code{1 - AF1}.
#'
#' SNPs with incompatible allele pairs are removed. When allele-frequency
#' filtering is enabled, SNPs are additionally retained only when the absolute
#' difference between the GWAS and LD effect-allele frequencies is smaller
#' than \code{af_threshold}.
#'
#' The LD information may contain either \code{A1Freq} or \code{A2Freq}.
#' When only \code{A2Freq} is available, \code{A1Freq} is calculated as
#' \code{1 - A2Freq}.
#'
#' This function does not perform general GWAS numeric quality control.
#' Numeric quality control of beta, standard error, P value, allele frequency,
#' and sample size is handled separately before gene-based analysis.
#'
#' Duplicate SNP identifiers are not checked by this function because they
#' are handled separately as gene-level defensive checks.
#'
#' @param gwas A data frame containing GWAS summary statistics. It must
#'   contain SNP identifier, A1, A2, AF1, and beta columns.
#' @param info_LD A data frame containing population-specific LD information.
#'   It must contain a SNP identifier, A1, A2, and a physical-position column.
#'   When \code{no_AF = FALSE}, either A1Freq or A2Freq must also be available.
#' @param no_AF Logical. If \code{TRUE}, allele-frequency difference filtering
#'   is not applied.
#' @param af_threshold Numeric. Maximum allele-frequency difference threshold
#'   used when \code{no_AF = FALSE}. Must lie between 0 and 1, inclusive.
#'   Defaults to \code{0.2}.
#'
#' @return A \code{data.table} containing SNPs retained after allele
#'   harmonisation and, when requested, allele-frequency filtering.
#'   The output includes the LD physical position in \code{PhysPos}.
#'   Returns \code{NULL} when no SNPs can be retained.
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
#'   gwas = gwas,
#'   info_LD = ld_info,
#'   no_AF = FALSE,
#'   af_threshold = 0.2
#' )
#'
#' print(result)
#'
#' @importFrom data.table as.data.table
#' @export
allele_flip_gwas_beta <- function(
    gwas,
    info_LD,
    no_AF,
    af_threshold = 0.2
) {

    # ============================================================
    # 1. Validate function arguments
    # ============================================================

    if (!is.data.frame(gwas)) {
        stop(
            "'gwas' must be a data frame or data.table.",
            call. = FALSE
        )
    }

    if (!is.data.frame(info_LD)) {
        stop(
            "'info_LD' must be a data frame or data.table.",
            call. = FALSE
        )
    }

    if (
        !is.logical(no_AF) ||
        length(no_AF) != 1L ||
        is.na(no_AF)
    ) {
        stop(
            "'no_AF' must be TRUE or FALSE.",
            call. = FALSE
        )
    }

    if (
        !is.numeric(af_threshold) ||
        length(af_threshold) != 1L ||
        is.na(af_threshold) ||
        !is.finite(af_threshold) ||
        af_threshold < 0 ||
        af_threshold > 1
    ) {
        stop(
            "'af_threshold' must be a finite numeric value between 0 and 1.",
            call. = FALSE
        )
    }


    # ============================================================
    # 2. Work on local copies
    # ============================================================

    gwas <- as.data.frame(
        gwas,
        stringsAsFactors = FALSE
    )

    info_LD <- as.data.frame(
        info_LD,
        stringsAsFactors = FALSE
    )


    # ============================================================
    # 3. Internal column finder
    # ============================================================

    find_col <- function(
        dat,
        candidates
    ) {

        idx <- match(
            tolower(candidates),
            tolower(names(dat))
        )

        idx <- idx[
            !is.na(idx)
        ]

        if (length(idx) == 0L) {
            return(NA_character_)
        }

        names(dat)[idx[1L]]
    }


    # ============================================================
    # 4. Identify required GWAS columns
    # ============================================================

    gwas_snp_col <- find_col(
        gwas,
        c("snp", "ID")
    )

    gwas_a1_col <- find_col(
        gwas,
        c("A1")
    )

    gwas_a2_col <- find_col(
        gwas,
        c("A2")
    )

    gwas_af1_col <- find_col(
        gwas,
        c("AF1", "A1Freq")
    )

    gwas_beta_col <- find_col(
        gwas,
        c("beta")
    )


    required_gwas <- c(
        gwas_snp_col,
        gwas_a1_col,
        gwas_a2_col,
        gwas_af1_col,
        gwas_beta_col
    )


    if (anyNA(required_gwas)) {
        stop(
            paste0(
                "'gwas' must contain SNP identifier, A1, A2, AF1, ",
                "and beta columns."
            ),
            call. = FALSE
        )
    }


    # ============================================================
    # 5. Standardise GWAS column names
    # ============================================================

    names(gwas)[
        names(gwas) == gwas_snp_col
    ] <- "snp"

    names(gwas)[
        names(gwas) == gwas_a1_col
    ] <- "A1"

    names(gwas)[
        names(gwas) == gwas_a2_col
    ] <- "A2"

    names(gwas)[
        names(gwas) == gwas_af1_col
    ] <- "AF1"

    names(gwas)[
        names(gwas) == gwas_beta_col
    ] <- "beta"


    # ============================================================
    # 6. Identify required LD columns
    # ============================================================

    ld_snp_col <- find_col(
        info_LD,
        c("snp", "ID")
    )

    ld_a1_col <- find_col(
        info_LD,
        c("A1")
    )

    ld_a2_col <- find_col(
        info_LD,
        c("A2")
    )

    ld_a1freq_col <- find_col(
        info_LD,
        c("A1Freq", "AF1")
    )

    ld_a2freq_col <- find_col(
        info_LD,
        c("A2Freq", "AF2")
    )

    ld_pos_col <- find_col(
        info_LD,
        c("PhysPos", "POS", "BP")
    )

    if (
        is.na(ld_snp_col) ||
        is.na(ld_a1_col) ||
        is.na(ld_a2_col) ||
        is.na(ld_pos_col)
    ) {
        return(NULL)
    }


    if (
        !no_AF &&
        is.na(ld_a1freq_col) &&
        is.na(ld_a2freq_col)
    ) {
        return(NULL)
    }


    # ============================================================
    # 7. Standardise LD column names
    # ============================================================

    names(info_LD)[
        names(info_LD) == ld_snp_col
    ] <- "snp"

    names(info_LD)[
        names(info_LD) == ld_a1_col
    ] <- "A1"

    names(info_LD)[
        names(info_LD) == ld_a2_col
    ] <- "A2"

    names(info_LD)[
        names(info_LD) == ld_pos_col
    ] <- "PhysPos"


    # ============================================================
    # 8. Obtain the LD A1 frequency
    # ============================================================

    if (!is.na(ld_a1freq_col)) {

        info_LD$A1Freq <- suppressWarnings(
            as.numeric(
                info_LD[[ld_a1freq_col]]
            )
        )

    } else if (!is.na(ld_a2freq_col)) {

        info_LD$A1Freq <- 1 - suppressWarnings(
            as.numeric(
                info_LD[[ld_a2freq_col]]
            )
        )

    } else {

        info_LD$A1Freq <- NA_real_
    }


    # ============================================================
    # 9. Standardise fields required for allele harmonisation
    # ============================================================

    gwas$snp <- as.character(
        gwas$snp
    )

    gwas$A1 <- as.character(
        gwas$A1
    )

    gwas$A2 <- as.character(
        gwas$A2
    )

    gwas$AF1 <- suppressWarnings(
        as.numeric(
            gwas$AF1
        )
    )

    gwas$beta <- suppressWarnings(
        as.numeric(
            gwas$beta
        )
    )


    info_LD$snp <- as.character(
        info_LD$snp
    )

    info_LD$A1 <- as.character(
        info_LD$A1
    )

    info_LD$A2 <- as.character(
        info_LD$A2
    )

    info_LD$PhysPos <- suppressWarnings(
        as.numeric(
            info_LD$PhysPos
        )
    )

    if (
        any(
            !is.finite(
                info_LD$PhysPos
            )
        )
    ) {
        return(NULL)
    }


    # ============================================================
    # 10. Merge GWAS and LD information by SNP identifier
    # ============================================================

    gwas$.original_order <- seq_len(
        nrow(gwas)
    )


    ld_sub <- info_LD[
        ,
        c(
            "snp",
            "A1",
            "A2",
            "A1Freq",
            "PhysPos"
        ),
        drop = FALSE
    ]


    names(ld_sub) <- c(
        "snp",
        "LD_A1",
        "LD_A2",
        "A1Freq",
        "PhysPos"
    )


    merged <- merge(
        gwas,
        ld_sub,
        by = "snp",
        all = FALSE,
        sort = FALSE
    )


    merged <- merged[
        order(
            merged$.original_order
        ),
        ,
        drop = FALSE
    ]


    if (nrow(merged) == 0L) {
        return(NULL)
    }


    # ============================================================
    # 11. Identify directly matched and reversed alleles
    # ============================================================

    matched <- (
        !is.na(merged$A1) &
        !is.na(merged$A2) &
        !is.na(merged$LD_A1) &
        !is.na(merged$LD_A2) &
        merged$A1 == merged$LD_A1 &
        merged$A2 == merged$LD_A2
    )


    reversed <- (
        !is.na(merged$A1) &
        !is.na(merged$A2) &
        !is.na(merged$LD_A1) &
        !is.na(merged$LD_A2) &
        merged$A1 == merged$LD_A2 &
        merged$A2 == merged$LD_A1
    )


    keep <- matched | reversed


    if (!any(keep)) {
        return(NULL)
    }


    merged <- merged[
        keep,
        ,
        drop = FALSE
    ]


    # ============================================================
    # 12. Flip reversed alleles
    # ============================================================

    reversed_now <- (
        merged$A1 == merged$LD_A2 &
        merged$A2 == merged$LD_A1
    )


    if (
        any(
            reversed_now,
            na.rm = TRUE
        )
    ) {

        merged$beta[
            reversed_now
        ] <- -merged$beta[
            reversed_now
        ]


        merged$AF1[
            reversed_now
        ] <- 1 - merged$AF1[
            reversed_now
        ]
    }

    merged$A1 <- merged$LD_A1
    merged$A2 <- merged$LD_A2


    # ============================================================
    # 13. Apply allele-frequency difference filtering
    # ============================================================

    if (!no_AF) {

        merged$dif_AF <- abs(
            merged$AF1 -
            merged$A1Freq
        )

        merged <- merged[
            !is.na(merged$dif_AF) &
            merged$dif_AF < af_threshold,
            ,
            drop = FALSE
        ]
    }


    if (nrow(merged) == 0L) {
        return(NULL)
    }


    # ============================================================
    # 14. Order retained SNPs by genomic position
    # ============================================================

    merged <- merged[
        order(
            merged$PhysPos
        ),
        ,
        drop = FALSE
    ]


    # ============================================================
    # 15. Remove temporary columns
    # ============================================================

    merged$.original_order <- NULL
    merged$LD_A1 <- NULL
    merged$LD_A2 <- NULL


    # ============================================================
    # 16. Return retained SNPs
    # ============================================================

    data.table::as.data.table(
        merged
    )
}