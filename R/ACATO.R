#' ACATO: Aggregated Cauchy Association Test
#'
#' Combines P-values using the Aggregated Cauchy Association Test (ACAT) method.
#' Handles missing values (NA) automatically. If only one P-value is valid,
#' it returns that P-value.
#'
#' @param p A numeric vector of P-values.
#' @param w A numeric vector of weights corresponding to each p-value. Default is NULL (equal weights).
#'
#' @return A single combined P-value (numeric). Returns NA if all inputs are NA.
#'
#' @importFrom stats pcauchy
#' @export
ACATO <- function(p, w = NULL) {

  if (all(is.na(p))) return(NA_real_)

  keep <- !is.na(p)
  p <- p[keep]
  
  if (!is.null(w)) {
    w <- w[keep]

    if (sum(w) <= 0) {
      w <- rep(1/length(p), length(p))
    } else {
      w <- w / sum(w) 
    }
  } else {
    w <- rep(1/length(p), length(p))
  }
  

  p[p == 1] <- 1 - 1e-16
  is.small <- (p < 1e-16)
  
  if (sum(is.small) == 0) {
    cct.stat <- sum(w * tan((0.5 - p) * pi))
  } else {
    cct.stat <- sum(w[is.small] * (1 / p[is.small]) / pi)
    cct.stat <- cct.stat + sum(w[!is.small] * tan((0.5 - p[!is.small]) * pi))
  }
  
  if (cct.stat > 1e+15) {
    pval <- (1 / cct.stat) / pi
  } else {
    pval <- 1 - pcauchy(cct.stat)
  }
  return(pval)
}