#' @title Calculate Gene-based LD Matrices using GCTB
#'
#' @description This function serves as a utility to preprocess data for
#' \code{multiBAT.combo} by calculating full LD matrices for specified gene
#' regions using PLINK to extract regions and GCTB's \code{--make-full-ldm}
#' command.
#'
#' @details The function handles two populations in parallel and creates separate
#' output folders (\code{pop1} and \code{pop2}) within the \code{result_prefix}
#' directory. It relies on external system calls to \code{plink}, which must be
#' available in PATH, and \code{gctb}, which can be supplied either as a command
#' available in PATH or as a full executable path. Parallel execution via
#' \code{mclapply} is platform-dependent and is best supported on Linux/macOS.
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
#' @param result_prefix Root directory where \code{pop1} and \code{pop2} output
#'   folders will be created.
#' @param thread Integer. Maximum number of genes to process concurrently.
#'   Defaults to \code{8}.
#'
#' @return Invisibly returns a list of processing results with
#' \code{Success}, \code{Skipped}, or \code{Failed} messages. The primary output
#' is the LD matrix files \code{*.ldm.full.bin}, \code{*.ldm.full.info}, and
#' \code{*.ldm.full.snp}, saved in \code{result_prefix/pop1} and
#' \code{result_prefix/pop2}.
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
    gctb_path <- Sys.which(gctb_bin)

    if (gctb_path == "" && file.exists(gctb_bin)) {
        gctb_path <- normalizePath(gctb_bin, mustWork = TRUE)
    }

    if (gctb_path == "") {
        stop(
            "GCTB executable not found: ", gctb_bin,
            ". Provide a command available in PATH or a full executable path.",
            call. = FALSE
        )
    }

    plink_path <- Sys.which("plink")
    if (plink_path == "") {
        stop(
            "PLINK command not found. Please ensure 'plink' is in the system PATH.",
            call. = FALSE
        )
    }

    if (!file.exists(map_file)) {
        stop("Gene annotation file not found: ", map_file, call. = FALSE)
    }

    if (is.null(result_prefix) || result_prefix == "") {
        stop("'result_prefix' must be specified.", call. = FALSE)
    }

    pop1_dir <- file.path(result_prefix, "pop1")
    pop2_dir <- file.path(result_prefix, "pop2")

    dir.create(pop1_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(pop2_dir, recursive = TRUE, showWarnings = FALSE)

    cat("Reading gene annotation file:", map_file, "\n")

    gene_map <- tryCatch(
        data.table::fread(map_file, header = FALSE, sep = "\t"),
        error = function(e) {
            stop(
                "Failed to read annotation file. Please check whether the file is tab-separated ",
                "and contains Chr, Start, End, and GeneName columns.\n",
                e$message,
                call. = FALSE
            )
        }
    )

    if (ncol(gene_map) < 4) {
        stop(
            "Annotation file has insufficient columns. It should contain at least four columns: ",
            "Chr, Start, End, and GeneName.",
            call. = FALSE
        )
    }

    data.table::setnames(gene_map, 1:4, c("Chr", "Start", "End", "Gene"))
    gene_list <- split(gene_map, seq(nrow(gene_map)))

    cat("Found", nrow(gene_map), "genes ready for processing.\n")

    .cleanup_temp_files <- function(prefix) {
        temp_extensions <- c(".bed", ".bim", ".fam", ".log", ".nosex", ".irem")
        temp_files <- paste0(prefix, temp_extensions)

        for (f in temp_files) {
            if (file.exists(f)) {
                suppressWarnings(file.remove(f))
            }
        }
    }

    process_gene_pop <- function(gene_row, plink_prefix, output_dir, pop_name) {
        chr <- gene_row$Chr
        start <- as.numeric(gene_row$Start)
        end <- as.numeric(gene_row$End)
        gene <- gene_row$Gene

        gene_output_prefix <- file.path(output_dir, gene)

        cat(
            paste0(
                "Processing gene for ", pop_name, " population: ",
                gene, " (", chr, ":", start, "-", end, ")\n"
            )
        )

        chr_clean <- gsub("^chr", "", chr, ignore.case = TRUE)
        start_ext <- max(1, start - extend)
        end_ext <- end + extend

        cat(paste0("Extended region: ", chr_clean, ":", start_ext, "-", end_ext, "\n"))

        plink_temp_prefix <- paste0(gene_output_prefix, ".temp")

        plink_args <- c(
            "--bfile", plink_prefix,
            "--chr", chr_clean,
            "--from-bp", start_ext,
            "--to-bp", end_ext,
            "--make-bed",
            "--out", plink_temp_prefix
        )

        plink_status <- system2(
            plink_path,
            args = plink_args,
            stdout = FALSE,
            stderr = FALSE,
            wait = TRUE
        )

        if (plink_status != 0) {
            cat(
                paste0(
                    "PLINK failed with status ", plink_status,
                    " for ", gene, " (", pop_name, "). Check PLINK logs.\n"
                )
            )
            .cleanup_temp_files(plink_temp_prefix)
            return(paste0("Failed: ", gene, " (PLINK execution error)"))
        }

        bim_file_path <- paste0(plink_temp_prefix, ".bim")

        if (!file.exists(bim_file_path) || file.info(bim_file_path)$size < 10) {
            cat(paste0("Skipping ", gene, " (", pop_name, "): No valid SNPs in the region.\n"))
            .cleanup_temp_files(plink_temp_prefix)
            return(paste0("Skipped: ", gene, " (No valid SNP in region)"))
        }

        gctb_args <- c(
            "--bfile", plink_temp_prefix,
            "--make-full-ldm",
            "--out", gene_output_prefix
        )

        cat("Calculating LD matrix using GCTB...\n")

        gctb_status <- system2(
            gctb_path,
            args = gctb_args,
            stdout = FALSE,
            stderr = FALSE,
            wait = TRUE
        )

        ldm_bin_file <- paste0(gene_output_prefix, ".ldm.full.bin")

        if (file.exists(ldm_bin_file)) {
            .cleanup_temp_files(plink_temp_prefix)
            cat(paste0("Completed ", gene, " (", pop_name, "). LD files saved in: ", output_dir, "\n"))
            return(paste0("Success: ", gene, " (", pop_name, ")"))
        }

        cat(
            paste0(
                "GCTB failed with status ", gctb_status,
                " for ", gene, " (", pop_name, "). Please check GCTB log for details.\n"
            )
        )

        .cleanup_temp_files(plink_temp_prefix)
        return(paste0("Failed: ", gene, " (", pop_name, ")"))
    }

    tasks_pop1 <- lapply(
        gene_list,
        function(g) list(row = g, prefix = plink_bfile_pop1, dir = pop1_dir, name = "pop1")
    )

    tasks_pop2 <- lapply(
        gene_list,
        function(g) list(row = g, prefix = plink_bfile_pop2, dir = pop2_dir, name = "pop2")
    )

    all_tasks <- c(tasks_pop1, tasks_pop2)

    cat(
        paste0(
            "Starting parallel processing (Threads: ",
            thread, "). Total tasks: ", length(all_tasks), "...\n"
        )
    )

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

    cat("\n--- Processing Complete ---\n")
    print(table(unlist(results)))

    invisible(results)
}