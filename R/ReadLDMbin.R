#' Read a GCTB-generated LD matrix from binary files
#'
#' This function reads a linkage disequilibrium (LD) matrix produced by GCTB.
#' It expects two files sharing a common prefix:
#' 1. \code{<prefix>.info}: a text file with SNP information and a header.
#' 2. \code{<prefix>.bin}: a binary file with the full LD matrix stored in
#'    row-major order as single-precision floats.
#'
#' If the input files are missing, inconsistent, or cannot be read, the function
#' returns \code{NULL} with a warning. This allows gene-level workflows to skip
#' the affected gene without stopping the full analysis.
#'
#' @param prefix Character string. File path prefix without extension for the
#'   GCTB LD files. The function appends \code{.info} and \code{.bin} to locate
#'   the input files.
#'
#' @return A numeric matrix of dimension \eqn{m \times m}, where \eqn{m} is the
#'   number of SNPs read from \code{<prefix>.info}. Returns \code{NULL} if the
#'   LD files are missing or inconsistent.
#'
#' @export
#'
#' @examples
#' td <- tempdir()
#' prefix <- file.path(td, "Chr1_gene1.ldm")
#'
#' info <- data.frame(
#'   ID = c("snp1", "snp2", "snp3"),
#'   stringsAsFactors = FALSE
#' )
#'
#' utils::write.table(
#'   info,
#'   file = paste0(prefix, ".info"),
#'   sep = "\t",
#'   row.names = FALSE,
#'   quote = FALSE
#' )
#'
#' mat <- matrix(
#'   c(1.0, 0.2, 0.1,
#'     0.2, 1.0, 0.3,
#'     0.1, 0.3, 1.0),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' bin_con <- file(paste0(prefix, ".bin"), "wb")
#' writeBin(as.numeric(mat), con = bin_con, size = 4, endian = "little")
#' close(bin_con)
#'
#' ld_matrix <- ReadLDMbin(prefix)
#' print(dim(ld_matrix))
ReadLDMbin <- function(prefix) {
  InfoFile <- paste0(prefix, ".info")
  BinFile <- paste0(prefix, ".bin")

  if (!file.exists(InfoFile)) {
    warning("LD info file not found: ", InfoFile, call. = FALSE)
    return(NULL)
  }

  if (!file.exists(BinFile)) {
    warning("LD binary file not found: ", BinFile, call. = FALSE)
    return(NULL)
  }

  info <- tryCatch(
    utils::read.table(InfoFile, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) {
      warning("Failed to read LD info file: ", InfoFile, ". ", e$message, call. = FALSE)
      return(NULL)
    }
  )

  if (is.null(info)) {
    return(NULL)
  }

  m <- nrow(info)

  if (m <= 0) {
    warning("LD info file contains no SNPs: ", InfoFile, call. = FALSE)
    return(NULL)
  }

  expected_size <- m * m * 4
  observed_size <- file.info(BinFile)$size

  if (is.na(observed_size) || observed_size != expected_size) {
    warning(
      "LD binary file size does not match the .info SNP count. Expected ",
      expected_size, " bytes for ", m, " SNPs, but observed ",
      observed_size, " bytes: ", BinFile,
      call. = FALSE
    )
    return(NULL)
  }

  bin_con <- file(BinFile, "rb")
  on.exit(close(bin_con), add = TRUE)

  vals <- readBin(
    bin_con,
    what = numeric(0),
    n = m * m,
    size = 4,
    endian = "little"
  )

  if (length(vals) != m * m) {
    warning(
      "Failed to read the expected number of LD values from binary file: ",
      BinFile,
      call. = FALSE
    )
    return(NULL)
  }

  matrix(vals, nrow = m, ncol = m, byrow = TRUE)
}