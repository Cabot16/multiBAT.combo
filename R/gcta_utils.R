#' Run GCTA mBAT-combo Analysis
#'
#' A wrapper function to execute GCTA's gene-based association test
#' \code{mBAT-combo}.
#'
#' @param sumstats A data.table/data.frame containing GWAS summary statistics.
#' @param gene_map A data.table/data.frame containing gene annotations.
#' @param bfile_path Character. Path to the PLINK reference panel binary prefix.
#' @param out_file Character. Full path and filename prefix for the output.
#' @param gcta_bin Character. GCTA executable name available in PATH, or full
#'   path to the GCTA executable. Defaults to \code{"gcta64"}.
#' @param af_threshold Numeric. Threshold for allele frequency difference.
#' @param gcta_wind_kb Numeric. Gene region extension in kilobases for GCTA's
#'   \code{--mBAT-wind} parameter.
#' @param threads Integer. Number of threads to use.
#'
#' @return A data.table containing the GCTA results. Returns \code{NULL} if
#'   execution fails.
#'
#' @importFrom data.table fwrite fread setnames data.table as.data.table
#' @importFrom parallel detectCores
#' @export
run_gcta_mbat <- function(
    sumstats,
    gene_map,
    bfile_path,
    out_file,
    gcta_bin = "gcta64",
    af_threshold = 0.2,
    gcta_wind_kb = 50, # Renamed from 'extend' to explicitly take KB value
    threads = NULL
) {
    # --- 1. Log definitions (ensures visibility within function) ---
    log_info  <- function(msg) message(sprintf("[%s] [GCTA] [INFO]  %s",  Sys.time(), msg))
    log_warn  <- function(msg) message(sprintf("[%s] [GCTA] [WARN]  %s",  Sys.time(), msg))

    # 2. Environment Check
    gcta_exec <- Sys.which(gcta_bin)
    if (gcta_exec == "" && file.exists(gcta_bin)) {
        gcta_exec <- normalizePath(gcta_bin, mustWork = TRUE)
    }

    if (gcta_exec == "") {
        log_warn(sprintf("GCTA executable not found at '%s'. Analysis skipped.", gcta_bin))
        return(NULL)
    }

    sumstats <- data.table::as.data.table(sumstats)
    gene_map <- data.table::as.data.table(gene_map)

    # Thread setting
    if (is.null(threads)) {
        sys_cores <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) 2)
        threads <- max(1, sys_cores - 1)
    }

    out_dir <- dirname(out_file)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    ## GCTA Wind Parameter 📐
    # WARNING: Removed redundant /1000 calculation. The input gcta_wind_kb is assumed to be in KB.
    # We use the received gcta_wind_kb directly for the GCTA command.
    log_info(sprintf("Setting GCTA --mBAT-wind to %g kb.", gcta_wind_kb))

    # 3. Prepare Input Files
    log_info("Preparing temporary input files...")

    # GWAS Summary (.ma file)
    req_cols <- c("snp", "A1", "A2", "AF1", "beta", "se", "P", "N")
    tmp_ma_file <- paste0(out_file, ".temp.ma")
    # Use .SD and .SDcols if sumstats is a data.table and requires subsetting
    gcta_ma <- sumstats[, ..req_cols]
    setnames(gcta_ma,
             old = c("snp", "A1", "A2", "AF1",  "beta", "se", "P", "N"),
             new = c("SNP", "A1", "A2", "freq", "b",    "se", "p", "N"))

    fwrite(gcta_ma, tmp_ma_file, sep = "\t", quote = FALSE)

    # Prepare Gene List File
    log_info("Preparing Gene list file...")

    # Find gene column name
    gene_col <- if ("gene.id" %in% names(gene_map)) "gene.id" else "gene"

    # Ensure column order is as required by GCTA (CHR, Start, End, GeneID)
    gcta_gene <- gene_map[, c("chr", "start", "end", gene_col), with = FALSE]

    # Standardize chromosome format (ensure GCTA compatibility)
    gcta_gene[, chr := tolower(chr)]
    gcta_gene[, chr := gsub("^chr", "", chr)]
    gcta_gene[, chr := as.character(chr)]

    tmp_gene_file <- paste0(out_file, ".temp.gene")

    # GCTA gene list should not have a header or quotes around fields
    data.table::fwrite(gcta_gene,
                        tmp_gene_file,
                        sep = "\t",
                        col.names = FALSE,
                        quote = FALSE)
    log_info(sprintf("Created temporary Gene list file: %s", tmp_gene_file))

    # 4. Construct Command (using shQuote to handle spaces or special characters in paths)
    # **Using the received gcta_wind_kb value directly**
    cmd <- sprintf(
        "%s --bfile %s --mBAT-combo %s --mBAT-gene-list %s --out %s --thread-num %d --diff-freq %s --mBAT-wind %s",
        shQuote(gcta_exec),
        shQuote(bfile_path),
        shQuote(tmp_ma_file),
        shQuote(tmp_gene_file),
        shQuote(out_file), # GCTA output prefix
        threads,
        as.character(af_threshold),
        as.character(gcta_wind_kb) # Now uses the correct KB value (e.g., 50.0)
    )

    log_info(sprintf("Executing GCTA command: %s", cmd))

    # 5. Execute Command
    exit_code <- system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)
    log_info(sprintf("GCTA command finished with exit code: %d", exit_code))

    # 6. Clean temporary files
    unlink(c(tmp_ma_file, tmp_gene_file))

    # 7. Check Results (using the confirmed suffix)
    expected_out <- paste0(out_file, ".gene.assoc.mbat")
    log_file <- paste0(out_file, ".log")

    if (exit_code == 0 && file.exists(expected_out) && nrow(fread(expected_out)) > 0) {
        log_info("GCTA mBAT-combo analysis finished successfully.")
        res <- fread(expected_out)
        return(res)
    } else {
        # Diagnose GCTA failure reason
        log_warn(sprintf("GCTA execution failed (Exit Code %d).", exit_code))

        if (file.exists(log_file)) {
            log_warn(sprintf("Please check GCTA LOG file: %s", log_file))
            # Read and print the first few lines of the log file to aid debugging
            log_lines <- tryCatch(readLines(log_file, n = 5), error = function(e) "Could not read log file.")
            if (length(log_lines) > 0) {
                log_warn(sprintf("--- Log Preview ---\n%s\n--- End Preview ---",
                                 paste(log_lines, collapse = "\n")))
            }
        } else {
            log_warn("GCTA execution failed and NO log file was generated. Check if PLINK files are correctly referenced by GCTA.")
            log_warn(sprintf("Check if the --bfile parameter is correct: %s", bfile_path))
        }

        return(NULL)
    }
}