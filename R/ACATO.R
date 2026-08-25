#' ACATO: Aggregated Cauchy Association Test
#'
#' Combines P values using the Aggregated Cauchy Association Test (ACAT)
#' method. Missing P values are removed automatically, and the corresponding
#' weights are removed before the remaining weights are normalised.
#'
#' If only one valid P value remains, that P value is returned directly.
#'
#' Extremely small P values are handled using an asymptotic approximation to
#' avoid numerical instability in the Cauchy transformation. A P value of zero
#' is permitted and results in a combined P value of zero when it contributes
#' to the test.
#'
#' @param p A numeric vector of P values.
#' @param w A numeric vector of weights corresponding to the P values.
#'   If \code{NULL}, equal weights are used. Defaults to \code{NULL}.
#'
#' @return A single combined P value. Returns \code{NA_real_} if all input
#'   P values are missing.
#'
#' @importFrom stats pcauchy
#' @export
ACATO <- function(
    p,
    w = NULL
) {

    # ------------------------------------------------------------
    # 1. Return NA when all P values are missing
    # ------------------------------------------------------------

    if (all(is.na(p))) {
        return(NA_real_)
    }


    # ------------------------------------------------------------
    # 2. Remove missing P values and their corresponding weights
    # ------------------------------------------------------------

    keep <- !is.na(p)

    p <- p[
        keep
    ]

    if (!is.null(w)) {
        w <- w[
            keep
        ]
    }


    # ------------------------------------------------------------
    # 3. Return the P value directly when only one value remains
    # ------------------------------------------------------------

    if (length(p) == 1L) {
        return(p)
    }


    # ------------------------------------------------------------
    # 4. Prepare and normalise weights
    # ------------------------------------------------------------

    if (!is.null(w)) {

        if (sum(w) <= 0) {

            w <- rep(
                1 / length(p),
                length(p)
            )

        } else {

            w <- w /
                sum(w)
        }

    } else {

        w <- rep(
            1 / length(p),
            length(p)
        )
    }


    # ------------------------------------------------------------
    # 5. Protect P values at the upper boundary
    # ------------------------------------------------------------

    p[p == 1] <- 1 - 1e-16


    # ------------------------------------------------------------
    # 6. Identify extremely small P values
    # ------------------------------------------------------------

    is.small <- (
        p < 1e-16
    )


    # ------------------------------------------------------------
    # 7. Calculate the Cauchy combination statistic
    # ------------------------------------------------------------

    if (sum(is.small) == 0L) {

        cct.stat <- sum(
            w *
                tan(
                    (0.5 - p) * pi
                )
        )

    } else {

        cct.stat <- sum(
            w[is.small] *
                (1 / p[is.small]) /
                pi
        )

        cct.stat <- cct.stat +
            sum(
                w[!is.small] *
                    tan(
                        (0.5 - p[!is.small]) *
                            pi
                    )
            )
    }


    # ------------------------------------------------------------
    # 8. Convert the Cauchy statistic to a combined P value
    # ------------------------------------------------------------

    if (cct.stat > 1e15) {

        pval <- (
            1 / cct.stat
        ) / pi

    } else {

        pval <- 1 -
            stats::pcauchy(
                cct.stat
            )
    }


    # ------------------------------------------------------------
    # 9. Return the combined P value
    # ------------------------------------------------------------

    return(pval)
}