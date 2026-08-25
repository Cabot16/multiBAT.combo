#' Quality Control of Numeric GWAS Summary Statistics
#'
#' Checks numeric fields in GWAS summary statistics and removes SNPs with
#' invalid values before downstream gene-based association analysis.
#'
#' The following checks are applied:
#'
#' \itemize{
#'   \item \code{beta} must be finite.
#'   \item \code{se} must be finite and greater than zero.
#'   \item \code{P} must be finite and between 0 and 1, inclusive.
#'   \item \code{AF1} must be finite and between 0 and 1, inclusive.
#'   \item \code{beta / se} must be finite.
#'   \item \code{N}, when checked, must be finite and greater than zero.
#' }
#'
#' A P value of zero is considered valid. SNPs failing one or more checks are
#' removed. A single summary warning is issued when any SNPs are removed.
#'
#' Sample-size validation can be disabled by setting \code{check_n = FALSE}.
#'
#' This function performs numeric quality control only. File existence,
#' file readability, required input structure, and gene-level defensive checks
#' are handled separately.
#'
#' @param gwas A data.frame or data.table containing GWAS summary statistics.
#' @param pop_label Character scalar used to identify the population.
#' @param check_n Logical. Whether to validate the \code{N} column.
#'
#' @return A data.table containing SNPs that passed numeric quality control.
#'
#' @export
gwas_numeric_qc <- function(
    gwas,
    pop_label = "GWAS",
    check_n = TRUE
) {

    # Validate arguments

    if (!is.data.frame(gwas)) {
        stop(
            "'gwas' must be a data frame or data.table.",
            call. = FALSE
        )
    }

    if (
        !is.character(pop_label) ||
        length(pop_label) != 1L ||
        is.na(pop_label) ||
        !nzchar(pop_label)
    ) {
        stop(
            "'pop_label' must be a non-empty character value of length 1.",
            call. = FALSE
        )
    }

    if (
        !is.logical(check_n) ||
        length(check_n) != 1L ||
        is.na(check_n)
    ) {
        stop(
            "'check_n' must be TRUE or FALSE.",
            call. = FALSE
        )
    }


    # Copy input

    gwas <- data.table::copy(
        data.table::as.data.table(gwas)
    )


    # Validate columns

    required_numeric_columns <- c(
        "beta",
        "se",
        "P",
        "AF1"
    )

    if (check_n) {
        required_numeric_columns <- c(
            required_numeric_columns,
            "N"
        )
    }

    missing_columns <- setdiff(
        required_numeric_columns,
        names(gwas)
    )

    if (length(missing_columns) > 0L) {
        stop(
            sprintf(
                paste0(
                    "%s GWAS summary statistics do not contain ",
                    "the required column(s) for numeric quality control: %s."
                ),
                pop_label,
                paste(missing_columns, collapse = ", ")
            ),
            call. = FALSE
        )
    }


    # Handle empty input

    if (nrow(gwas) == 0L) {
        return(gwas)
    }


    # Convert fields

    gwas[["beta"]] <- suppressWarnings(
        as.numeric(gwas[["beta"]])
    )

    gwas[["se"]] <- suppressWarnings(
        as.numeric(gwas[["se"]])
    )

    gwas[["P"]] <- suppressWarnings(
        as.numeric(gwas[["P"]])
    )

    gwas[["AF1"]] <- suppressWarnings(
        as.numeric(gwas[["AF1"]])
    )

    if (check_n) {
        gwas[["N"]] <- suppressWarnings(
            as.numeric(gwas[["N"]])
        )
    }


    # Identify invalid values

    invalid_beta <- !is.finite(
        gwas[["beta"]]
    )

    invalid_se <- !is.finite(
        gwas[["se"]]
    ) |
        gwas[["se"]] <= 0

    invalid_p <- !is.finite(
        gwas[["P"]]
    ) |
        gwas[["P"]] < 0 |
        gwas[["P"]] > 1

    invalid_af <- !is.finite(
        gwas[["AF1"]]
    ) |
        gwas[["AF1"]] < 0 |
        gwas[["AF1"]] > 1

    z_score <- gwas[["beta"]] / gwas[["se"]]

    invalid_z <- !is.finite(
        z_score
    )

    if (check_n) {

        invalid_n <- !is.finite(
            gwas[["N"]]
        ) |
            gwas[["N"]] <= 0

    } else {

        invalid_n <- rep(
            FALSE,
            nrow(gwas)
        )
    }


    # Combine filters

    invalid_snp <- (
        invalid_beta |
        invalid_se |
        invalid_p |
        invalid_af |
        invalid_z |
        invalid_n
    )

    n_input <- nrow(gwas)

    n_removed <- sum(
        invalid_snp
    )


    # Remove invalid SNPs

    if (n_removed > 0L) {

        gwas <- gwas[
            !invalid_snp
        ]


        # Report QC warning

        qc_details <- c(
            sprintf(
                "invalid beta=%d",
                sum(invalid_beta)
            ),
            sprintf(
                "invalid SE=%d",
                sum(invalid_se)
            ),
            sprintf(
                "invalid P value=%d",
                sum(invalid_p)
            ),
            sprintf(
                "invalid AF=%d",
                sum(invalid_af)
            ),
            sprintf(
                "non-finite beta/SE=%d",
                sum(invalid_z)
            )
        )

        if (check_n) {
            qc_details <- c(
                qc_details,
                sprintf(
                    "invalid N=%d",
                    sum(invalid_n)
                )
            )
        }

        warning(
            sprintf(
                paste0(
                    "%s GWAS numeric quality control removed %d of %d SNPs ",
                    "because they contained invalid numeric values (%s). ",
                    "%d SNPs were retained."
                ),
                pop_label,
                n_removed,
                n_input,
                paste(qc_details, collapse = "; "),
                nrow(gwas)
            ),
            call. = FALSE
        )
    }


    # Return result

    gwas
}
