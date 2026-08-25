#' Run GCTA mBAT-combo Analysis
#'
#' Runs GCTA's gene-based mBAT-combo analysis for a single population.
#'
#' The function prepares temporary GCTA summary-statistics and gene-list
#' files, runs GCTA with \code{--mBAT-combo} and
#' \code{--mBAT-print-all-p}, and returns the resulting gene-based
#' association results.
#'
#' GWAS numeric quality control is not performed within this function.
#' Numeric quality control of beta, standard error, P value, allele frequency,
#' and sample size should be completed before this function is called.
#'
#' If \code{sample_size} is supplied, the supplied value is used as the
#' sample size for every SNP in the temporary GCTA summary-statistics file.
#' In this case, an \code{N} column is not required in \code{sumstats}.
#' If \code{sample_size = NULL}, the \code{N} column in \code{sumstats}
#' is used.
#'
#' GCTA is run with \code{--mBAT-print-all-p}, so the output is expected to
#' contain the population-specific \code{P_mBAT}, \code{P_fastBAT}, and
#' \code{P_mBATcombo} values. In the multiBAT-combo workflow,
#' \code{P_mBAT} is used as the single-population component of the final
#' weighted ACAT analysis.
#'
#' @param sumstats A data.frame or data.table containing GWAS summary
#'   statistics. Required columns are \code{snp}, \code{A1}, \code{A2},
#'   \code{AF1}, \code{beta}, \code{se}, and \code{P}. The \code{N}
#'   column is additionally required when \code{sample_size = NULL}.
#' @param gene_map A data.frame or data.table containing gene annotations.
#'   Required columns are \code{chr}, \code{start}, \code{end}, and either
#'   \code{gene.id} or \code{gene}.
#' @param bfile_path Character. Path prefix of the PLINK reference panel.
#'   The corresponding \code{.bed}, \code{.bim}, and \code{.fam} files
#'   must exist.
#' @param out_file Character. Full path and filename prefix for GCTA output.
#' @param gcta_bin Character. GCTA executable name available in PATH, or the
#'   full path to the executable. Defaults to \code{"gcta64"}.
#' @param af_threshold Numeric. Allele-frequency difference threshold supplied
#'   to GCTA through \code{--diff-freq}. Defaults to \code{0.2}.
#' @param gcta_wind_kb Numeric. Gene-region extension in kilobases supplied to
#'   GCTA through \code{--mBAT-wind}. Defaults to \code{50}.
#' @param threads Integer. Number of threads used by GCTA. If \code{NULL},
#'   one fewer than the number of detected physical CPU cores is used.
#' @param sample_size Optional externally supplied sample size. Must be a
#'   finite positive numeric scalar when provided. If supplied, this value
#'   replaces the GWAS \code{N} column in the temporary GCTA input.
#'
#' @return A data.table containing the GCTA gene-based association results.
#'   Returns \code{NULL} if GCTA fails or produces no readable result.
#'
#' @importFrom data.table fwrite fread setnames data.table as.data.table copy
#' @importFrom parallel detectCores
#' @export
run_gcta_mbat <- function(
    sumstats,
    gene_map,
    bfile_path,
    out_file,
    gcta_bin = "gcta64",
    af_threshold = 0.2,
    gcta_wind_kb = 50,
    threads = NULL,
    sample_size = NULL
) {

    # ------------------------------------------------------------
    # 1. Logging
    # ------------------------------------------------------------

    log_info <- function(msg) {
        message(
            sprintf(
                "[%s] [GCTA] [INFO] %s",
                Sys.time(),
                msg
            )
        )
    }

    log_warn <- function(msg) {
        warning(
            sprintf(
                "[%s] [GCTA] [WARN] %s",
                Sys.time(),
                msg
            ),
            call. = FALSE
        )
    }


    # ------------------------------------------------------------
    # 2. Validate basic parameters
    # ------------------------------------------------------------

    if (
        !is.character(bfile_path) ||
        length(bfile_path) != 1L ||
        is.na(bfile_path) ||
        !nzchar(bfile_path)
    ) {
        stop(
            "'bfile_path' must be a non-empty character value of length 1.",
            call. = FALSE
        )
    }

    if (
        !is.character(out_file) ||
        length(out_file) != 1L ||
        is.na(out_file) ||
        !nzchar(out_file)
    ) {
        stop(
            "'out_file' must be a non-empty character value of length 1.",
            call. = FALSE
        )
    }

    if (
        !is.character(gcta_bin) ||
        length(gcta_bin) != 1L ||
        is.na(gcta_bin) ||
        !nzchar(gcta_bin)
    ) {
        stop(
            "'gcta_bin' must be a non-empty character value of length 1.",
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

    if (
        !is.numeric(gcta_wind_kb) ||
        length(gcta_wind_kb) != 1L ||
        is.na(gcta_wind_kb) ||
        !is.finite(gcta_wind_kb) ||
        gcta_wind_kb < 0
    ) {
        stop(
            "'gcta_wind_kb' must be a finite non-negative numeric value.",
            call. = FALSE
        )
    }

    if (!is.null(sample_size)) {

        if (
            !is.numeric(sample_size) ||
            length(sample_size) != 1L ||
            is.na(sample_size) ||
            !is.finite(sample_size) ||
            sample_size <= 0
        ) {
            stop(
                paste0(
                    "'sample_size' must be a finite positive numeric value ",
                    "when supplied."
                ),
                call. = FALSE
            )
        }
    }


    # ------------------------------------------------------------
    # 3. Locate the GCTA executable
    # ------------------------------------------------------------

    gcta_exec <- Sys.which(
        gcta_bin
    )

    if (
        identical(gcta_exec, "") &&
        file.exists(gcta_bin)
    ) {
        gcta_exec <- normalizePath(
            gcta_bin,
            mustWork = TRUE
        )
    }

    if (identical(gcta_exec, "")) {
        stop(
            sprintf(
                paste0(
                    "GCTA executable was not found at '%s'. ",
                    "Provide an executable available in PATH or a full path."
                ),
                gcta_bin
            ),
            call. = FALSE
        )
    }


    # ------------------------------------------------------------
    # 4. Confirm that the PLINK reference files exist
    # ------------------------------------------------------------

    required_bfile_extensions <- c(
        ".bed",
        ".bim",
        ".fam"
    )

    required_bfile_files <- paste0(
        bfile_path,
        required_bfile_extensions
    )

    missing_bfile_files <- required_bfile_files[
        !file.exists(
            required_bfile_files
        )
    ]

    if (length(missing_bfile_files) > 0L) {
        stop(
            sprintf(
                "The PLINK reference is incomplete. Missing file(s): %s",
                paste(
                    missing_bfile_files,
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
    }


    # ------------------------------------------------------------
    # 5. Copy input objects
    # ------------------------------------------------------------

    sumstats <- data.table::copy(
        data.table::as.data.table(
            sumstats
        )
    )

    gene_map <- data.table::copy(
        data.table::as.data.table(
            gene_map
        )
    )


    # ------------------------------------------------------------
    # 6. Check required GWAS columns
    # ------------------------------------------------------------

    req_cols_base <- c(
        "snp",
        "A1",
        "A2",
        "AF1",
        "beta",
        "se",
        "P"
    )

    missing_sumstat_cols <- setdiff(
        req_cols_base,
        names(sumstats)
    )

    if (length(missing_sumstat_cols) > 0L) {
        stop(
            paste0(
                "The GWAS summary statistics are missing required column(s): ",
                paste(
                    missing_sumstat_cols,
                    collapse = ", "
                ),
                "."
            ),
            call. = FALSE
        )
    }


    # ------------------------------------------------------------
    # 7. Determine the sample-size source
    # ------------------------------------------------------------

    if (is.null(sample_size)) {

        if (!"N" %in% names(sumstats)) {
            stop(
                paste0(
                    "The GWAS summary statistics must contain an 'N' column ",
                    "when 'sample_size' is not supplied."
                ),
                call. = FALSE
            )
        }

    } else {

        sumstats[, N := sample_size]
    }


    # ------------------------------------------------------------
    # 8. Check required gene-map columns
    # ------------------------------------------------------------

    gene_col <- if ("gene.id" %in% names(gene_map)) {

        "gene.id"

    } else if ("gene" %in% names(gene_map)) {

        "gene"

    } else {

        stop(
            "The gene map must contain either a 'gene.id' or 'gene' column.",
            call. = FALSE
        )
    }

    required_gene_cols <- c(
        "chr",
        "start",
        "end",
        gene_col
    )

    missing_gene_cols <- setdiff(
        required_gene_cols,
        names(gene_map)
    )

    if (length(missing_gene_cols) > 0L) {
        stop(
            paste0(
                "The gene map is missing required column(s): ",
                paste(
                    missing_gene_cols,
                    collapse = ", "
                ),
                "."
            ),
            call. = FALSE
        )
    }


    # ------------------------------------------------------------
    # 9. Validate gene annotation values
    # ------------------------------------------------------------

    gene_chr <- as.character(
        gene_map[["chr"]]
    )

    gene_start <- suppressWarnings(
        as.numeric(
            gene_map[["start"]]
        )
    )

    gene_end <- suppressWarnings(
        as.numeric(
            gene_map[["end"]]
        )
    )

    gene_id <- as.character(
        gene_map[[gene_col]]
    )

    invalid_gene <- (
        is.na(gene_chr) |
        !nzchar(trimws(gene_chr)) |
        !is.finite(gene_start) |
        !is.finite(gene_end) |
        gene_start < 1 |
        gene_end < gene_start |
        is.na(gene_id) |
        !nzchar(trimws(gene_id))
    )

    if (any(invalid_gene)) {
        stop(
            sprintf(
                paste0(
                    "The gene map contains %d invalid gene annotation ",
                    "record(s)."
                ),
                sum(invalid_gene)
            ),
            call. = FALSE
        )
    }


    # ------------------------------------------------------------
    # 10. Determine the number of GCTA threads
    # ------------------------------------------------------------

    if (is.null(threads)) {

        sys_cores <- tryCatch(
            parallel::detectCores(
                logical = FALSE
            ),
            error = function(e) {
                2L
            }
        )

        if (
            length(sys_cores) != 1L ||
            is.na(sys_cores) ||
            !is.finite(sys_cores) ||
            sys_cores < 1L
        ) {
            sys_cores <- 2L
        }

        threads <- max(
            1L,
            as.integer(sys_cores) - 1L
        )

    } else {

        if (
            !is.numeric(threads) ||
            length(threads) != 1L ||
            is.na(threads) ||
            !is.finite(threads) ||
            threads < 1 ||
            abs(threads - round(threads)) >
                sqrt(.Machine$double.eps)
        ) {
            stop(
                "'threads' must be a positive whole number.",
                call. = FALSE
            )
        }

        threads <- as.integer(
            threads
        )
    }


    # ------------------------------------------------------------
    # 11. Prepare the output directory
    # ------------------------------------------------------------

    out_dir <- dirname(
        out_file
    )

    if (!dir.exists(out_dir)) {

        dir.create(
            out_dir,
            recursive = TRUE,
            showWarnings = FALSE
        )
    }

    if (!dir.exists(out_dir)) {
        stop(
            sprintf(
                "The GCTA output directory could not be created: %s",
                out_dir
            ),
            call. = FALSE
        )
    }

    log_info(
        sprintf(
            "Setting GCTA --mBAT-wind to %g kb.",
            gcta_wind_kb
        )
    )


    # ------------------------------------------------------------
    # 12. Prepare temporary file names
    # ------------------------------------------------------------

    tmp_ma_file <- paste0(
        out_file,
        ".temp.ma"
    )

    tmp_gene_file <- paste0(
        out_file,
        ".temp.gene"
    )

    on.exit(
        unlink(
            c(
                tmp_ma_file,
                tmp_gene_file
            )
        ),
        add = TRUE
    )


    # ------------------------------------------------------------
    # 13. Prepare the GCTA summary-statistics file
    # ------------------------------------------------------------

    log_info(
        "Preparing the temporary GCTA summary-statistics file."
    )

    req_cols <- c(
        "snp",
        "A1",
        "A2",
        "AF1",
        "beta",
        "se",
        "P",
        "N"
    )

    gcta_ma <- sumstats[
        ,
        ..req_cols
    ]

    data.table::setnames(
        gcta_ma,
        old = c(
            "snp",
            "A1",
            "A2",
            "AF1",
            "beta",
            "se",
            "P",
            "N"
        ),
        new = c(
            "SNP",
            "A1",
            "A2",
            "freq",
            "b",
            "se",
            "p",
            "N"
        )
    )

    data.table::fwrite(
        gcta_ma,
        file = tmp_ma_file,
        sep = "\t",
        quote = FALSE,
        na = "NA"
    )

    log_info(
        sprintf(
            "Created the temporary GCTA summary-statistics file with %d SNPs.",
            nrow(gcta_ma)
        )
    )


    # ------------------------------------------------------------
    # 14. Prepare the GCTA gene-list file
    # ------------------------------------------------------------

    log_info(
        "Preparing the temporary GCTA gene-list file."
    )

    gcta_gene <- gene_map[
        ,
        c(
            "chr",
            "start",
            "end",
            gene_col
        ),
        with = FALSE
    ]

    gcta_gene[, chr := as.character(chr)]

    gcta_gene[
        ,
        chr := gsub(
            "^chr",
            "",
            chr,
            ignore.case = TRUE
        )
    ]

    gcta_gene[, start := as.integer(start)]

    gcta_gene[, end := as.integer(end)]


    data.table::fwrite(
        gcta_gene,
        file = tmp_gene_file,
        sep = "\t",
        col.names = FALSE,
        quote = FALSE,
        na = "NA"
    )

    log_info(
        sprintf(
            "Created the temporary GCTA gene-list file with %d genes.",
            nrow(gcta_gene)
        )
    )


    # ------------------------------------------------------------
    # 15. Construct the GCTA command
    # ------------------------------------------------------------

    gcta_args <- c(
        "--bfile",
        bfile_path,
        "--mBAT-combo",
        tmp_ma_file,
        "--mBAT-gene-list",
        tmp_gene_file,
        "--mBAT-print-all-p",
        "--out",
        out_file,
        "--thread-num",
        as.character(threads),
        "--diff-freq",
        as.character(af_threshold),
        "--mBAT-wind",
        as.character(gcta_wind_kb)
    )

    log_info(
        paste(
            "Executing GCTA command:",
            shQuote(gcta_exec),
            paste(
                shQuote(gcta_args),
                collapse = " "
            )
        )
    )


    # ------------------------------------------------------------
    # 16. Run GCTA
    # ------------------------------------------------------------

    exit_code <- system2(
        command = gcta_exec,
        args = gcta_args,
        stdout = "",
        stderr = ""
    )

    log_info(
        sprintf(
            "GCTA command finished with exit code %d.",
            exit_code
        )
    )


    # ------------------------------------------------------------
    # 17. Read the GCTA result
    # ------------------------------------------------------------

    expected_out <- paste0(
        out_file,
        ".gene.assoc.mbat"
    )

    log_file <- paste0(
        out_file,
        ".log"
    )

    if (
        exit_code == 0L &&
        file.exists(expected_out)
    ) {

        res <- tryCatch(
            data.table::fread(
                expected_out
            ),
            error = function(e) {

                log_warn(
                    sprintf(
                        "The GCTA output could not be read: %s",
                        e$message
                    )
                )

                NULL
            }
        )

        if (
            !is.null(res) &&
            nrow(res) > 0L
        ) {

            log_info(
                sprintf(
                    paste0(
                        "GCTA mBAT-combo analysis completed successfully ",
                        "with %d gene results."
                    ),
                    nrow(res)
                )
            )

            return(
                res
            )
        }
    }


    # ------------------------------------------------------------
    # 18. Report GCTA failure
    # ------------------------------------------------------------

    log_warn(
        sprintf(
            paste0(
                "GCTA execution failed or produced no readable result ",
                "(exit code %d)."
            ),
            exit_code
        )
    )

    if (file.exists(log_file)) {

        log_warn(
            sprintf(
                "Check the GCTA log file for details: %s",
                log_file
            )
        )

        log_lines <- tryCatch(
            readLines(
                log_file,
                n = 10L
            ),
            error = function(e) {
                character()
            }
        )

        if (length(log_lines) > 0L) {

            message(
                paste0(
                    "--- GCTA log preview ---\n",
                    paste(
                        log_lines,
                        collapse = "\n"
                    ),
                    "\n--- End of GCTA log preview ---"
                )
            )
        }

    } else {

        log_warn(
            paste0(
                "GCTA did not generate a log file. ",
                "Check the PLINK reference and GCTA command inputs."
            )
        )
    }

    return(
        NULL
    )
}