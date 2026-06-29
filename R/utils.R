#' Internal helper utilities and global variable definitions
#'
#' A collection of miscellaneous helper functions and a set of
#' global variable declarations used internally in the mBAT.multipop package.
#'
#' @name utils
#' @keywords internal
NULL

#' Define Global Variables to Satisfy R CMD Check
#'
#' Prevents "no visible binding for global variable" notes during package
#' checks by declaring commonly used variable names as globals.
#'
#' @name utils
#' @rdname utils
#' @keywords internal
NULL

utils::globalVariables(
  c(".", "snp", "chr", "pos", "A1", "A2",
    "gene", "gene.start", "gene.end",
    "pvalue", "chisq",
    "P_multiBAT", "P_multiBAT.combo", 
    "Gene", "P_joint", "P_pop1", "P_pop2",
    "..req_cols" 
    )
)

#' Verify Required Columns in a Data Frame
#'
#' Ensures that the supplied data frame contains all of the specified
#' column names. If any columns are missing, the function throws an error
#' identifying which names are absent.
#'
#' @param df A data frame to be checked.
#' @param columns A character vector of column names that must be present
#'   in `df`.
#'
#' @return Invisibly returns `NULL` if all required columns are found.
#'   Otherwise stops execution with an error.
#'
#' @export
#'
#' @examples
#' # passes silently
#' columns_check(data.frame(a = 1, b = 2), c("a", "b"))
#'
#' # throws an error because "c" is missing
#' \dontrun{
#' columns_check(data.frame(a = 1), c("a", "c"))
#' }
columns_check <- function(df, columns) {
  df_name       <- deparse(substitute(df))
  df_columns    <- names(df)
  diff_columns  <- setdiff(columns, df_columns)
  if (length(diff_columns) > 0) {
    stop(
      "'", df_name, "' doesn't have column(s): ",
      paste(diff_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
}

#' Check That an Object Is a Data Frame
#'
#' Validates that the provided object is a data frame. If it is not,
#' execution is halted with an informative error.
#'
#' @param df The object to test.
#'
#' @return Invisibly returns `NULL` if `df` is a data frame.
#'   Otherwise, stops execution with an error.
#'
#' @export
#'
#' @examples
#' df_check(data.frame(x = 1:3))
#' \dontrun{
#' df_check(list(a = 1))
#' }
df_check <- function(df) {
  df_name <- deparse(substitute(df))
  if (!is.data.frame(df)) {
    stop(
      "'", df_name, "' is not a data frame.",
      call. = FALSE
    )
  }
}

#' Validate a Non-negative Numeric Argument
#'
#' Checks that an input is a single, non-negative numeric value.
#' Throws an error if the test fails.
#'
#' @param x A value to check.
#' @param arg A string naming the argument (used in the error message).
#'
#' @return Invisibly returns `TRUE` if `x` passes the checks.
#'   Otherwise, stops execution with an error.
#'
#' @export
#'
#' @examples
#' nonneg_num_check(0, "threshold")
#' \dontrun{
#' nonneg_num_check(-1, "threshold")
#' }
nonneg_num_check <- function(x, arg) {
  if (!is.numeric(x) || is.na(x) || x < 0L || length(x) != 1L) {
    stop(
      "'", arg, "' must be a non-negative number of length 1.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate a Numeric Argument Within a Range
#'
#' Checks that a numeric input lies between two specified bounds (inclusive).
#' If `x` falls outside the range or is not a single numeric value, an error
#' is thrown.
#'
#' @param x The value to check.
#' @param left The lower bound (inclusive).
#' @param right The upper bound (inclusive).
#' @param arg A string naming the argument (used in the error message).
#'
#' @return Invisibly returns `TRUE` if `x` is valid.
#'   Otherwise, stops execution with an error.
#'
#' @export
#'
#' @examples
#' num_between_check(0.5, 0, 1, "probability")
#' \dontrun{
#' num_between_check(1.5, 0, 1, "probability")
#' }
num_between_check <- function(x, left, right, arg) {
  if (!is.numeric(x) || is.na(x) || length(x) != 1L ||
      x < left || x > right) {
    stop(
      "'", arg, "' must be a number of length 1 between ",
      left, " and ", right, ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Check for Integer-Valued Numeric Vector
#'
#' Determines whether all elements of a numeric vector are effectively
#' integers (within a specified tolerance). Returns `TRUE` or `FALSE`.
#'
#' @param x Numeric vector to examine.
#' @param tol Numeric tolerance for detecting integer values. Defaults to
#'   machine double epsilon.
#'
#' @return A logical scalar: `TRUE` if every element of `x` is a whole
#'   number (within `tol`), otherwise `FALSE`.
#'
#' @export
#'
#' @examples
#' integer_vector_check(c(1, 2, 3))
#' integer_vector_check(c(1, 2.1, 3))
integer_vector_check <- function(x, tol = .Machine$double.eps) {
  is_wholenumber <- function(x, tol) {
    !is.na(x) && abs(x - round(x)) < tol
  }
  tryCatch(
    all(sapply(x, is_wholenumber, tol = tol)),
    error = function(e) FALSE
  )
}

#' Format Numbers with Thousands Separators
#'
#' A wrapper around `prettyNum()` that inserts commas as thousands
#' separators and disables scientific notation.
#'
#' @param x Numeric or character vector to format.
#' @param ... Additional arguments passed to `prettyNum()`.
#'
#' @return A character vector with formatted numbers.
#'
#' @export
#'
#' @examples
#' pretty_num(1234567)
#' pretty_num(c(1e3, 2e5), digits = 1)
pretty_num <- function(x, ...) {
  prettyNum(x, big.mark = ",", scientific = FALSE, ...)
}

#' Identify Indices of Duplicate Elements
#'
#' Returns the indices of all elements in `x` that appear more than once.
#'
#' @param x A vector to check for duplicates.
#'
#' @return An integer vector of positions corresponding to duplicated
#'   entries (including both the first and subsequent occurrences).
#'
#' @export
#'
#' @examples
#' if_dup_ind(c("a", "b", "a", "c", "b"))
if_dup_ind <- function(x) {
  which(duplicated(x) | duplicated(x, fromLast = TRUE))
}
