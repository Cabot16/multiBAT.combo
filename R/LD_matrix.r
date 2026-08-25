#' @title Calculate Gene-based LD Matrices using GCTB
#'
#' @description This function serves as a utility to preprocess data for
#' \code{multiBAT.combo} by calculating full LD matrices for specified gene
#' regions using PLINK to extract regions and GCTB's \code{--make-full-ldm}
#' command.
#'
#' @details The function processes two populations and creates separate
#' output directories named \code{pop1} and \code{pop2} within
#' \code{result_prefix}. For each gene, the GCTB output prefix is
#' \code{<result_prefix>/pop1/<Gene>} for population 1 and
#' \code{<result_prefix>/pop2/<Gene>} for population 2.
#'
#' The primary LD outputs for each gene are therefore:
#' \code{<Gene>.ldm.full.bin}, \code{<Gene>.ldm.full.info}, and
#' \code{<Gene>.ldm.full.snp}.
#'
#' Temporary PLINK files use the prefix \code{<Gene>.temp} and are removed
#' after each gene has been processed.
#'
#' The function relies on external system calls to \code{plink}, which must be
#' available in PATH, and \code{gctb}, which can be supplied either as a
#' command available in PATH or as a full executable path. Parallel execution
#' via \code{mclapply} is platform-dependent and is best supported on
#' Linux and macOS.
#'
#' Console output is limited to one final status line per gene, together with
#' an initial task summary and a final processing summary.
#'
#' @param gctb_bin GCTB executable name available in PATH, or full path to the
#'   GCTB executable.
#' @param map_file Path to the gene annotation file. The file should have no
#'   header and contain at least four tab-separated columns:
#'   \code{Chr}, \code{Start}, \code{End}, and \code{GeneName}.
#' @param plink_bfile_pop1 PLINK bfile prefix for population 1.
#' @param plink_bfile_pop2 PLINK bfile prefix for population 2.
#' @param extend Integer. Number of base pairs to extend the gene boundary.
#'   Defaults to \code{50000}.
#' @param result_prefix Root output directory. Population-specific outputs are
#'   written to \code{<result_prefix>/pop1} and
#'   \code{<result_prefix>/pop2}.
#' @param thread Integer. Maximum number of genes to process concurrently.
#'   Defaults to \code{8}.
#'
#' @return Invisibly returns a list of processing results with
#' \code{Success}, \code{Skipped}, or \code{Failed} messages. The primary
#' outputs are the LD matrix files \code{*.ldm.full.bin},
#' \code{*.ldm.full.info}, and \code{*.ldm.full.snp}, saved in
#' \code{<result_prefix>/pop1} and \code{<result_prefix>/pop2}.
#'
#' @importFrom data.table fread setnames
#' @importFrom parallel mclapply
#' @export
LD_matrix <- function(
    gctb_bin,
    map_file,
    plink_bfile_pop1,
    plink_bfile_pop2,
    extend = 50000,
    result_prefix,
    thread = 8
) {

    # ============================================================
    # 1. Locate required external software
    # ============================================================

    gctb_path <- Sys.which(
        gctb_bin
    )

    if (
        gctb_path == "" &&
        file.exists(gctb_bin)
    ) {
        gctb_path <- normalizePath(
            gctb_bin,
            mustWork = TRUE
        )
    }

    if (gctb_path == "") {
        stop(
            "GCTB executable not found: ",
            gctb_bin,
            ". Provide a command available in PATH or a full executable path.",
            call. = FALSE
        )
    }


    plink_path <- Sys.which(
        "plink"
    )

    if (plink_path == "") {
        stop(
            paste0(
                "PLINK command not found. Ensure that 'plink' is available ",
                "in the system PATH."
            ),
            call. = FALSE
        )
    }


    # ============================================================
    # 2. Check required inputs
    # ============================================================

    if (!file.exists(map_file)) {
        stop(
            "Gene annotation file not found: ",
            map_file,
            call. = FALSE
        )
    }


    if (
        is.null(result_prefix) ||
        result_prefix == ""
    ) {
        stop(
            "'result_prefix' must be specified.",
            call. = FALSE
        )
    }


    # ============================================================
    # 3. Prepare output directories
    # ============================================================

    pop1_dir <- file.path(
        result_prefix,
        "pop1"
    )

    pop2_dir <- file.path(
        result_prefix,
        "pop2"
    )


    dir.create(
        pop1_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    dir.create(
        pop2_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )


    # ============================================================
    # 4. Read gene annotation
    # ============================================================

    gene_map <- tryCatch(
        data.table::fread(
            map_file,
            header = FALSE,
            sep = "\t"
        ),
        error = function(e) {
            stop(
                paste0(
                    "Failed to read the annotation file. Check that the file ",
                    "is tab-separated and contains Chr, Start, End, and ",
                    "GeneName columns. ",
                    e$message
                ),
                call. = FALSE
            )
        }
    )


    if (ncol(gene_map) < 4) {
        stop(
            paste0(
                "The annotation file has insufficient columns. It must contain ",
                "at least four columns: Chr, Start, End, and GeneName."
            ),
            call. = FALSE
        )
    }


    data.table::setnames(
        gene_map,
        1:4,
        c(
            "Chr",
            "Start",
            "End",
            "Gene"
        )
    )


    gene_list <- split(
        gene_map,
        seq_len(
            nrow(gene_map)
        )
    )


    # ============================================================
    # 5. Report the analysis setup
    # ============================================================

    cat(
        sprintf(
            paste0(
                "Gene annotation loaded: %d genes.\n",
                "Output directories:\n",
                "  POP1: %s\n",
                "  POP2: %s\n",
                "Parallel workers: %d\n"
            ),
            nrow(gene_map),
            pop1_dir,
            pop2_dir,
            thread
        )
    )


    # ============================================================
    # 6. Helper for removing temporary PLINK files
    # ============================================================

    .cleanup_temp_files <- function(prefix) {

        temp_extensions <- c(
            ".bed",
            ".bim",
            ".fam",
            ".log",
            ".nosex",
            ".irem"
        )


        temp_files <- paste0(
            prefix,
            temp_extensions
        )


        for (f in temp_files) {

            if (file.exists(f)) {

                suppressWarnings(
                    file.remove(f)
                )
            }
        }
    }


    # ============================================================
    # 7. Helper for one-line gene status output
    # ============================================================

    .print_gene_status <- function(
        status,
        pop_name,
        gene,
        detail = NULL
    ) {

        if (
            is.null(detail) ||
            !nzchar(detail)
        ) {

            cat(
                sprintf(
                    "[%s] %s | %s\n",
                    status,
                    pop_name,
                    gene
                )
            )

        } else {

            cat(
                sprintf(
                    "[%s] %s | %s | %s\n",
                    status,
                    pop_name,
                    gene,
                    detail
                )
            )
        }
    }


    # ============================================================
    # 8. Process one gene for one population
    # ============================================================

    process_gene_pop <- function(
        gene_row,
        plink_prefix,
        output_dir,
        pop_name
    ) {

        chr <- gene_row$Chr

        start <- as.numeric(
            gene_row$Start
        )

        end <- as.numeric(
            gene_row$End
        )

        gene <- gene_row$Gene


        gene_output_prefix <- file.path(
            output_dir,
            gene
        )


        chr_clean <- gsub(
            "^chr",
            "",
            chr,
            ignore.case = TRUE
        )


        start_ext <- max(
            1,
            start - extend
        )

        end_ext <- end +
            extend


        # --------------------------------------------------------
        # Temporary PLINK prefix
        # --------------------------------------------------------

        plink_temp_prefix <- paste0(
            gene_output_prefix,
            ".temp"
        )


        # --------------------------------------------------------
        # Extract the gene region using PLINK
        # --------------------------------------------------------

        plink_args <- c(
            "--bfile",
            plink_prefix,
            "--chr",
            chr_clean,
            "--from-bp",
            start_ext,
            "--to-bp",
            end_ext,
            "--make-bed",
            "--out",
            plink_temp_prefix
        )


        plink_status <- system2(
            plink_path,
            args = plink_args,
            stdout = FALSE,
            stderr = FALSE,
            wait = TRUE
        )


        if (plink_status != 0) {

            .print_gene_status(
                status = "WARNING",
                pop_name = pop_name,
                gene = gene,
                detail = sprintf(
                    "PLINK failed with exit code %d.",
                    plink_status
                )
            )


            .cleanup_temp_files(
                plink_temp_prefix
            )


            return(
                paste0(
                    "Failed: ",
                    gene,
                    " (",
                    pop_name,
                    "; PLINK execution error)"
                )
            )
        }


        # --------------------------------------------------------
        # Check whether the region contains SNPs
        # --------------------------------------------------------

        bim_file_path <- paste0(
            plink_temp_prefix,
            ".bim"
        )


        if (
            !file.exists(bim_file_path) ||
            file.info(bim_file_path)$size < 10
        ) {

            .print_gene_status(
                status = "SKIPPED",
                pop_name = pop_name,
                gene = gene,
                detail = "No valid SNPs were found in the region."
            )


            .cleanup_temp_files(
                plink_temp_prefix
            )


            return(
                paste0(
                    "Skipped: ",
                    gene,
                    " (",
                    pop_name,
                    "; no valid SNPs in region)"
                )
            )
        }


        # --------------------------------------------------------
        # Calculate the full LD matrix using GCTB
        # --------------------------------------------------------

        gctb_args <- c(
            "--bfile",
            plink_temp_prefix,
            "--make-full-ldm",
            "--out",
            gene_output_prefix
        )


        gctb_status <- system2(
            gctb_path,
            args = gctb_args,
            stdout = FALSE,
            stderr = FALSE,
            wait = TRUE
        )


        ldm_bin_file <- paste0(
            gene_output_prefix,
            ".ldm.full.bin"
        )


        # --------------------------------------------------------
        # Successful gene
        # --------------------------------------------------------

        if (file.exists(ldm_bin_file)) {

            .cleanup_temp_files(
                plink_temp_prefix
            )


            .print_gene_status(
                status = "SUCCESS",
                pop_name = pop_name,
                gene = gene
            )


            return(
                paste0(
                    "Success: ",
                    gene,
                    " (",
                    pop_name,
                    ")"
                )
            )
        }


        # --------------------------------------------------------
        # GCTB failure
        # --------------------------------------------------------

        .print_gene_status(
            status = "WARNING",
            pop_name = pop_name,
            gene = gene,
            detail = sprintf(
                "GCTB failed with exit code %d.",
                gctb_status
            )
        )


        .cleanup_temp_files(
            plink_temp_prefix
        )


        return(
            paste0(
                "Failed: ",
                gene,
                " (",
                pop_name,
                "; GCTB execution error)"
            )
        )
    }


    # ============================================================
    # 9. Prepare population-specific tasks
    # ============================================================

    tasks_pop1 <- lapply(
        gene_list,
        function(g) {

            list(
                row = g,
                prefix = plink_bfile_pop1,
                dir = pop1_dir,
                name = "pop1"
            )
        }
    )


    tasks_pop2 <- lapply(
        gene_list,
        function(g) {

            list(
                row = g,
                prefix = plink_bfile_pop2,
                dir = pop2_dir,
                name = "pop2"
            )
        }
    )


    all_tasks <- c(
        tasks_pop1,
        tasks_pop2
    )


    cat(
        sprintf(
            "Starting LD-matrix calculation for %d population-gene tasks.\n",
            length(all_tasks)
        )
    )


    # ============================================================
    # 10. Run tasks in parallel
    # ============================================================

    results <- parallel::mclapply(
        all_tasks,
        function(task) {

            process_gene_pop(
                gene_row = task$row,
                plink_prefix = task$prefix,
                output_dir = task$dir,
                pop_name = task$name
            )
        },
        mc.cores = thread
    )


    # ============================================================
    # 11. Summarise processing results
    # ============================================================

    result_messages <- unlist(
        results,
        use.names = FALSE
    )


    n_success <- sum(
        grepl(
            "^Success:",
            result_messages
        )
    )


    n_skipped <- sum(
        grepl(
            "^Skipped:",
            result_messages
        )
    )


    n_failed <- sum(
        grepl(
            "^Failed:",
            result_messages
        )
    )


    cat(
        sprintf(
            paste0(
                "\nLD-matrix processing completed.\n",
                "Success: %d\n",
                "Skipped: %d\n",
                "Failed: %d\n",
                "Total tasks: %d\n"
            ),
            n_success,
            n_skipped,
            n_failed,
            length(result_messages)
        )
    )


    # ============================================================
    # 12. Return processing results
    # ============================================================

    invisible(
        results
    )
}