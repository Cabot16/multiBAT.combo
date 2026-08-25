# ============================================================
# input_validation.R
#
# Validate parameters, input structures, and LD matching.
# GWAS numeric QC is handled separately.
# ============================================================




# Validate logical scalar

.validate_logical_scalar <- function(x, arg) {

    if (
        !is.logical(x) ||
        length(x) != 1L ||
        is.na(x)
    ) {
        stop(
            sprintf(
                "'%s' must be TRUE or FALSE.",
                arg
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


# Validate character scalar

.validate_character_scalar <- function(x, arg) {

    if (
        !is.character(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !nzchar(trimws(x))
    ) {
        stop(
            sprintf(
                "'%s' must be a non-empty character value of length 1.",
                arg
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


# Validate whole number

.validate_nonnegative_whole_number <- function(x, arg) {

    if (
        !is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        x < 0 ||
        abs(x - round(x)) > sqrt(.Machine$double.eps)
    ) {
        stop(
            sprintf(
                "'%s' must be a finite, non-negative whole number.",
                arg
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


# Validate numeric range

.validate_numeric_range <- function(
    x,
    lower,
    upper,
    arg,
    lower_inclusive = TRUE,
    upper_inclusive = TRUE
) {

    if (
        !is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x)
    ) {
        stop(
            sprintf(
                "'%s' must be a finite numeric value of length 1.",
                arg
            ),
            call. = FALSE
        )
    }

    lower_invalid <- if (lower_inclusive) {
        x < lower
    } else {
        x <= lower
    }

    upper_invalid <- if (upper_inclusive) {
        x > upper
    } else {
        x >= upper
    }

    if (lower_invalid || upper_invalid) {

        lower_symbol <- if (lower_inclusive) "[" else "("
        upper_symbol <- if (upper_inclusive) "]" else ")"

        stop(
            sprintf(
                "'%s' must lie within %s%s, %s%s.",
                arg,
                lower_symbol,
                lower,
                upper,
                upper_symbol
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


# Validate sample size

.validate_external_sample_size <- function(x, arg) {

    if (is.null(x)) {
        return(
            invisible(TRUE)
        )
    }

    if (
        !is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        x <= 0
    ) {
        stop(
            sprintf(
                "'%s' must be a finite positive numeric value when supplied.",
                arg
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


# Validate input file

.validate_input_file <- function(path, arg) {

    .validate_character_scalar(
        path,
        arg
    )

    if (!file.exists(path)) {
        stop(
            sprintf(
                "Input file specified by '%s' was not found: %s",
                arg,
                path
            ),
            call. = FALSE
        )
    }

    if (dir.exists(path)) {
        stop(
            sprintf(
                "'%s' must refer to a file rather than a directory: %s",
                arg,
                path
            ),
            call. = FALSE
        )
    }

    if (file.access(path, mode = 4L) != 0L) {
        stop(
            sprintf(
                "Input file specified by '%s' is not readable: %s",
                arg,
                path
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


# Validate input directory

.validate_input_directory <- function(path, arg) {

    .validate_character_scalar(
        path,
        arg
    )

    if (!dir.exists(path)) {
        stop(
            sprintf(
                "Input directory specified by '%s' was not found: %s",
                arg,
                path
            ),
            call. = FALSE
        )
    }

    if (file.access(path, mode = 4L) != 0L) {
        stop(
            sprintf(
                "Input directory specified by '%s' is not readable: %s",
                arg,
                path
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


# Find input column

.find_input_column <- function(
    dat,
    candidates,
    field,
    source_label,
    required = TRUE
) {

    input_names <- names(
        dat
    )

    candidate_lower <- unique(
        tolower(
            candidates
        )
    )

    hit <- which(
        tolower(input_names) %in%
            candidate_lower
    )

    if (length(hit) == 0L) {

        if (!required) {
            return(
                NA_character_
            )
        }

        stop(
            sprintf(
                paste0(
                    "%s does not contain a column that can be recognised ",
                    "as '%s'. Accepted column name(s): %s."
                ),
                source_label,
                field,
                paste(
                    candidates,
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
    }

    if (length(hit) > 1L) {
        stop(
            sprintf(
                paste0(
                    "%s contains multiple columns that could be recognised ",
                    "as '%s': %s."
                ),
                source_label,
                field,
                paste(
                    input_names[hit],
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
    }

    input_names[
        hit
    ]
}


# Validate GWAS structure

.read_validate_gwas <- function(
    path,
    source_label,
    require_n = TRUE
) {

    dat <- tryCatch(
        data.table::fread(
            path,
            header = TRUE,
            check.names = FALSE,
            showProgress = FALSE
        ),
        error = function(e) {
            stop(
                sprintf(
                    "Failed to read %s: %s",
                    source_label,
                    e$message
                ),
                call. = FALSE
            )
        }
    )

    if (nrow(dat) == 0L) {
        stop(
            sprintf(
                "%s contains no data rows.",
                source_label
            ),
            call. = FALSE
        )
    }

    if (ncol(dat) == 0L) {
        stop(
            sprintf(
                "%s contains no readable columns.",
                source_label
            ),
            call. = FALSE
        )
    }

    col_snp <- .find_input_column(
        dat,
        c(
            "snp",
            "ID"
        ),
        "snp",
        source_label
    )

    col_a1 <- .find_input_column(
        dat,
        c(
            "A1"
        ),
        "A1",
        source_label
    )

    col_a2 <- .find_input_column(
        dat,
        c(
            "A2"
        ),
        "A2",
        source_label
    )

    col_af1 <- .find_input_column(
        dat,
        c(
            "AF1",
            "A1Freq"
        ),
        "AF1",
        source_label
    )

    col_beta <- .find_input_column(
        dat,
        c(
            "beta",
            "b"
        ),
        "beta",
        source_label
    )

    col_se <- .find_input_column(
        dat,
        c(
            "se"
        ),
        "se",
        source_label
    )

    col_p <- .find_input_column(
        dat,
        c(
            "P"
        ),
        "P",
        source_label
    )

    col_n <- .find_input_column(
        dat,
        c(
            "N"
        ),
        "N",
        source_label,
        required = require_n
    )

    if (
        !require_n &&
        is.na(col_n)
    ) {

        col_n_optional <- .find_input_column(
            dat,
            c(
                "N"
            ),
            "N",
            source_label,
            required = FALSE
        )

        col_n <- col_n_optional
    }


    old_names <- c(
        col_snp,
        col_a1,
        col_a2,
        col_af1,
        col_beta,
        col_se,
        col_p
    )

    new_names <- c(
        "snp",
        "A1",
        "A2",
        "AF1",
        "beta",
        "se",
        "P"
    )

    if (!is.na(col_n)) {

        old_names <- c(
            old_names,
            col_n
        )

        new_names <- c(
            new_names,
            "N"
        )
    }


    data.table::setnames(
        dat,
        old = old_names,
        new = new_names
    )

    dat
}


# Validate BIM structure

.read_validate_bim <- function(path) {

    bim <- tryCatch(
        data.table::fread(
            path,
            header = FALSE,
            check.names = FALSE,
            showProgress = FALSE
        ),
        error = function(e) {
            stop(
                sprintf(
                    "Failed to read the PLINK BIM file: %s",
                    e$message
                ),
                call. = FALSE
            )
        }
    )

    if (nrow(bim) == 0L) {
        stop(
            "The PLINK BIM file contains no SNP records.",
            call. = FALSE
        )
    }

    if (ncol(bim) != 6L) {
        stop(
            sprintf(
                paste0(
                    "The PLINK BIM file must contain exactly six columns, ",
                    "but %d columns were detected."
                ),
                ncol(bim)
            ),
            call. = FALSE
        )
    }

    data.table::setnames(
        bim,
        c(
            "chr",
            "snp",
            "bp",
            "pos",
            "A1",
            "A2"
        )
    )

    pos_numeric <- suppressWarnings(
        as.numeric(
            bim[["pos"]]
        )
    )

    invalid_pos <- (
        !is.finite(pos_numeric) |
        pos_numeric < 0
    )

    if (any(invalid_pos)) {
        stop(
            sprintf(
                paste0(
                    "The PLINK BIM file contains %d row(s) with invalid genomic ",
                    "positions. Confirm that the file is a headerless six-column ",
                    "PLINK BIM file."
                ),
                sum(invalid_pos)
            ),
            call. = FALSE
        )
    }

    bim[["pos"]] <- pos_numeric


    required_text_columns <- c(
        "chr",
        "snp",
        "A1",
        "A2"
    )

    for (column_name in required_text_columns) {

        values <- as.character(
            bim[[column_name]]
        )

        invalid_text <- (
            is.na(values) |
            !nzchar(
                trimws(values)
            )
        )

        if (any(invalid_text)) {
            stop(
                sprintf(
                    "The PLINK BIM file contains %d row(s) with missing '%s' values.",
                    sum(invalid_text),
                    column_name
                ),
                call. = FALSE
            )
        }
    }

    bim
}


# Validate gene annotation

.read_validate_gene_map <- function(path) {

    gene_map <- tryCatch(
        data.table::fread(
            path,
            header = FALSE,
            sep = "\t",
            check.names = FALSE,
            showProgress = FALSE
        ),
        error = function(e) {
            stop(
                sprintf(
                    "Failed to read the gene annotation file: %s",
                    e$message
                ),
                call. = FALSE
            )
        }
    )

    if (nrow(gene_map) == 0L) {
        stop(
            "The gene annotation file contains no gene records.",
            call. = FALSE
        )
    }

    if (ncol(gene_map) != 4L) {
        stop(
            sprintf(
                paste0(
                    "The gene annotation file must contain exactly four ",
                    "tab-separated columns: chr, start, end, and gene.id. ",
                    "%d columns were detected."
                ),
                ncol(gene_map)
            ),
            call. = FALSE
        )
    }

    data.table::setnames(
        gene_map,
        c(
            "chr",
            "start",
            "end",
            "gene.id"
        )
    )


    gene_map[["chr"]] <- as.character(
        gene_map[["chr"]]
    )

    gene_map[["gene.id"]] <- as.character(
        gene_map[["gene.id"]]
    )


    start_numeric <- suppressWarnings(
        as.numeric(
            gene_map[["start"]]
        )
    )

    end_numeric <- suppressWarnings(
        as.numeric(
            gene_map[["end"]]
        )
    )


    invalid_chr <- (
        is.na(gene_map[["chr"]]) |
        !nzchar(
            trimws(
                gene_map[["chr"]]
            )
        )
    )

    invalid_gene <- (
        is.na(gene_map[["gene.id"]]) |
        !nzchar(
            trimws(
                gene_map[["gene.id"]]
            )
        )
    )

    invalid_start <- (
        !is.finite(start_numeric) |
        start_numeric < 1
    )

    invalid_end <- (
        !is.finite(end_numeric) |
        end_numeric < 1
    )

    invalid_interval <- (
        !invalid_start &
        !invalid_end &
        end_numeric < start_numeric
    )


    if (
        any(invalid_chr) ||
        any(invalid_gene) ||
        any(invalid_start) ||
        any(invalid_end) ||
        any(invalid_interval)
    ) {

        stop(
            sprintf(
                paste0(
                    "The gene annotation file contains invalid records ",
                    "(missing chromosome=%d; missing gene identifier=%d; ",
                    "invalid start=%d; invalid end=%d; end before start=%d)."
                ),
                sum(invalid_chr),
                sum(invalid_gene),
                sum(invalid_start),
                sum(invalid_end),
                sum(invalid_interval)
            ),
            call. = FALSE
        )
    }


    gene_map[["start"]] <- start_numeric
    gene_map[["end"]] <- end_numeric

    gene_map
}


# Validate LD directory

.validate_ld_directory <- function(
    path,
    pop_label
) {

    bin_files <- list.files(
        path,
        pattern = "\\.ldm\\.full\\.bin$",
        full.names = TRUE
    )

    info_files <- list.files(
        path,
        pattern = "\\.ldm\\.full\\.info$",
        full.names = TRUE
    )


    if (length(bin_files) == 0L) {
        stop(
            sprintf(
                "%s LD directory contains no '.ldm.full.bin' files: %s",
                pop_label,
                path
            ),
            call. = FALSE
        )
    }


    if (length(info_files) == 0L) {
        stop(
            sprintf(
                "%s LD directory contains no '.ldm.full.info' files: %s",
                pop_label,
                path
            ),
            call. = FALSE
        )
    }


    # Strip the complete GCTB suffix so that:
    #
    # GENE1.test.ldm.full.bin
    #
    # becomes:
    #
    # GENE1.test

    bin_prefix <- sub(
        "\\.ldm\\.full\\.bin$",
        "",
        basename(
            bin_files
        )
    )


    info_prefix <- sub(
        "\\.ldm\\.full\\.info$",
        "",
        basename(
            info_files
        )
    )


    paired_prefix <- intersect(
        bin_prefix,
        info_prefix
    )


    paired_prefix <- unique(
        paired_prefix
    )


    if (length(paired_prefix) == 0L) {
        stop(
            sprintf(
                paste0(
                    "%s LD directory does not contain any matched ",
                    "'.ldm.full.bin' and '.ldm.full.info' file pairs: %s"
                ),
                pop_label,
                path
            ),
            call. = FALSE
        )
    }


    paired_prefix
}


# ------------------------------------------------------------
# Internal helper: resolve one gene to one unique LD prefix
#
# Matching rule:
#
# 1. A broad match accepts either:
#
#      GENE1
#
#    or any prefix beginning with:
#
#      GENE1.
#
# 2. If more than one complete LD pair matches, retain only
#    prefixes for which the first character immediately after
#    "GENE1." is non-numeric. An exact "GENE1" prefix is retained.
#
# 3. If exactly one candidate remains, it is selected.
#
# 4. If the match remains ambiguous, or if the secondary rule
#    removes all candidates, the conflicting files are reported
#    and analysis is stopped.
#
# 5. If no complete LD pair matches at the broad matching stage,
#    NA_character_ is returned. Missing LD for an individual gene
#    is therefore allowed and can be handled downstream.
# ------------------------------------------------------------

.resolve_ld_prefix <- function(
    ld_dir,
    gene,
    pop_label = "LD",
    paired_prefixes = NULL
) {

    .validate_character_scalar(
        gene,
        "gene"
    )


    if (is.null(paired_prefixes)) {

        .validate_input_directory(
            ld_dir,
            "ld_dir"
        )

        paired_prefixes <- .validate_ld_directory(
            ld_dir,
            pop_label
        )
    }


    paired_prefixes <- as.character(
        paired_prefixes
    )


    # ------------------------------------------------------------
    # Broad gene matching
    #
    # This allows:
    #
    #   GENE1
    #   GENE1.test
    #   GENE1.EUR
    #   GENE1.123
    #
    # but prevents:
    #
    #   GENE10
    #   GENE100
    #
    # from matching GENE1.
    # ------------------------------------------------------------

    gene_dot_prefix <- paste0(
        gene,
        "."
    )


    candidate_prefixes <- paired_prefixes[
        paired_prefixes == gene |
        startsWith(
            paired_prefixes,
            gene_dot_prefix
        )
    ]


    candidate_prefixes <- unique(
        candidate_prefixes
    )


    # ------------------------------------------------------------
    # No complete LD pair matched this gene
    # ------------------------------------------------------------

    if (length(candidate_prefixes) == 0L) {
        return(
            NA_character_
        )
    }


    # ------------------------------------------------------------
    # Exactly one complete LD pair matched
    # ------------------------------------------------------------

    if (length(candidate_prefixes) == 1L) {

        return(
            file.path(
                ld_dir,
                paste0(
                    candidate_prefixes,
                    ".ldm.full"
                )
            )
        )
    }


    # ------------------------------------------------------------
    # Multiple broad matches were found
    #
    # Apply the secondary rule:
    #
    # the first character immediately after "<gene>." must not
    # be numeric.
    #
    # An exact gene prefix is retained.
    # ------------------------------------------------------------

    original_candidates <- candidate_prefixes


    keep_secondary <- vapply(
        candidate_prefixes,
        function(candidate) {

            if (identical(
                candidate,
                gene
            )) {
                return(
                    TRUE
                )
            }


            if (!startsWith(
                candidate,
                gene_dot_prefix
            )) {
                return(
                    FALSE
                )
            }


            suffix <- substring(
                candidate,
                nchar(gene_dot_prefix) + 1L
            )


            if (!nzchar(suffix)) {
                return(
                    TRUE
                )
            }


            first_character <- substr(
                suffix,
                1L,
                1L
            )


            !grepl(
                "^[0-9]$",
                first_character
            )
        },
        logical(1)
    )


    secondary_candidates <- candidate_prefixes[
        keep_secondary
    ]


    # ------------------------------------------------------------
    # The secondary rule resolved the ambiguity
    # ------------------------------------------------------------

    if (length(secondary_candidates) == 1L) {

        return(
            file.path(
                ld_dir,
                paste0(
                    secondary_candidates,
                    ".ldm.full"
                )
            )
        )
    }


    # ------------------------------------------------------------
    # Matching remains ambiguous
    #
    # If the secondary rule retained multiple candidates, report
    # those remaining candidates.
    #
    # If it removed every candidate, report all original broad
    # matches so that the naming conflict can be inspected.
    # ------------------------------------------------------------

    candidates_to_report <- if (
        length(secondary_candidates) > 0L
    ) {
        secondary_candidates
    } else {
        original_candidates
    }


    conflicting_files <- unlist(
        lapply(
            candidates_to_report,
            function(candidate) {

                c(
                    file.path(
                        ld_dir,
                        paste0(
                            candidate,
                            ".ldm.full.info"
                        )
                    ),
                    file.path(
                        ld_dir,
                        paste0(
                            candidate,
                            ".ldm.full.bin"
                        )
                    )
                )
            }
        ),
        use.names = FALSE
    )


    stop(
        paste0(
            "Multiple LD file pairs matched gene '",
            gene,
            "' in ",
            pop_label,
            " and could not be resolved uniquely.\n",
            "Conflicting files:\n",
            paste(
                paste0(
                    "  ",
                    conflicting_files
                ),
                collapse = "\n"
            )
        ),
        call. = FALSE
    )
}


# Validate gene-LD matching

.validate_ld_gene_matching <- function(
    ld_dir,
    genes,
    pop_label,
    paired_prefixes
) {

    genes <- unique(
        as.character(
            genes
        )
    )


    genes <- genes[
        !is.na(genes) &
        nzchar(
            trimws(
                genes
            )
        )
    ]


    if (length(genes) == 0L) {
        stop(
            "No valid gene identifiers were available for LD-file validation.",
            call. = FALSE
        )
    }

    invisible(
        vapply(
            genes,
            function(gene) {

                .resolve_ld_prefix(
                    ld_dir = ld_dir,
                    gene = gene,
                    pop_label = pop_label,
                    paired_prefixes = paired_prefixes
                )
            },
            character(1)
        )
    )


    invisible(TRUE)
}


# Resolve GCTA executable

.resolve_gcta_executable <- function(gcta_bin) {

    .validate_character_scalar(
        gcta_bin,
        "gcta_bin"
    )

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

    gcta_exec
}


# Validate PLINK prefix

.validate_plink_bfile <- function(
    prefix,
    arg
) {

    .validate_character_scalar(
        prefix,
        arg
    )

    required_extensions <- c(
        ".bed",
        ".bim",
        ".fam"
    )

    required_files <- paste0(
        prefix,
        required_extensions
    )

    missing_files <- required_files[
        !file.exists(
            required_files
        )
    ]

    if (length(missing_files) > 0L) {
        stop(
            sprintf(
                paste0(
                    "PLINK reference specified by '%s' is incomplete. ",
                    "Missing file(s): %s"
                ),
                arg,
                paste(
                    missing_files,
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
    }


    unreadable_files <- required_files[
        file.access(
            required_files,
            mode = 4L
        ) != 0L
    ]

    if (length(unreadable_files) > 0L) {
        stop(
            sprintf(
                paste0(
                    "PLINK reference specified by '%s' contains unreadable ",
                    "file(s): %s"
                ),
                arg,
                paste(
                    unreadable_files,
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


# ============================================================
# Main validation function
# ============================================================


#' Validate Inputs for multiBAT.combo
#'
#' Performs parameter, file-existence, file-readability, structural, and LD
#' file-matching validation before multiBAT-combo analysis begins.
#'
#' GWAS numeric quality control is not performed by this function. Numeric
#' validation of beta, standard error, P value, allele frequency, and sample
#' size is handled separately by \code{gwas_numeric_qc()}.
#'
#' Duplicate SNP identifiers are not checked here because they are handled as
#' gene-level defensive checks during downstream analysis.
#'
#' For each gene in the supplied gene annotation, population-specific GCTB LD
#' files are checked using the package's gene-to-LD matching rules. Missing LD
#' files for individual genes are permitted. Ambiguous gene-to-LD matches stop
#' the analysis before gene-based testing begins.
#'
#' @param bim_file Path to a headerless six-column PLINK BIM file.
#' @param map_file Path to a headerless four-column gene annotation file.
#' @param assoc_pop1 Path to population 1 GWAS summary statistics.
#' @param assoc_pop2 Path to population 2 GWAS summary statistics.
#' @param LD_path_pop1 Directory containing population 1 GCTB LD files.
#' @param LD_path_pop2 Directory containing population 2 GCTB LD files.
#' @param result_prefix Output path and filename prefix.
#' @param extend Non-negative whole number specifying the gene-region extension
#'   in base pairs.
#' @param prop Proportion of eigen-variance to retain. Must be greater than zero
#'   and no greater than one.
#' @param af_threshold Allele-frequency difference threshold. Must lie between
#'   zero and one, inclusive.
#' @param run_gcta Logical. Whether GCTA is used for the single-population mBAT
#'   analyses.
#' @param print_all_p Logical. Whether intermediate multiBAT and
#'   single-population P values are included in the main output.
#' @param check_gwas_values Logical. Whether GWAS numeric quality control will
#'   be performed before analysis.
#' @param N_pop1 Optional externally supplied population 1 sample size.
#' @param N_pop2 Optional externally supplied population 2 sample size.
#' @param gcta_bin GCTA executable name or full path.
#' @param plink_bfile_pop1 PLINK reference prefix for population 1. Required
#'   when \code{run_gcta = TRUE}.
#' @param plink_bfile_pop2 PLINK reference prefix for population 2. Required
#'   when \code{run_gcta = TRUE}.
#'
#' @return Invisibly returns a list containing structurally validated input
#'   objects and information required by the main analysis.
#'
#' @keywords internal
#'
#' @importFrom data.table fread setnames
validate_multibat_inputs <- function(
    bim_file,
    map_file,
    assoc_pop1,
    assoc_pop2,
    LD_path_pop1,
    LD_path_pop2,
    result_prefix,
    extend = 50000,
    prop = 0.9,
    af_threshold = 0.2,
    run_gcta = FALSE,
    print_all_p = FALSE,
    check_gwas_values = TRUE,
    N_pop1 = NULL,
    N_pop2 = NULL,
    gcta_bin = "gcta64",
    plink_bfile_pop1 = NULL,
    plink_bfile_pop2 = NULL
) {

    # ============================================================
    # 1. Validate function parameters
    # ============================================================

    .validate_character_scalar(
        result_prefix,
        "result_prefix"
    )

    .validate_nonnegative_whole_number(
        extend,
        "extend"
    )

    .validate_numeric_range(
        prop,
        lower = 0,
        upper = 1,
        arg = "prop",
        lower_inclusive = FALSE,
        upper_inclusive = TRUE
    )

    .validate_numeric_range(
        af_threshold,
        lower = 0,
        upper = 1,
        arg = "af_threshold",
        lower_inclusive = TRUE,
        upper_inclusive = TRUE
    )

    .validate_logical_scalar(
        run_gcta,
        "run_gcta"
    )

    .validate_logical_scalar(
        print_all_p,
        "print_all_p"
    )

    .validate_logical_scalar(
        check_gwas_values,
        "check_gwas_values"
    )

    .validate_external_sample_size(
        N_pop1,
        "N_pop1"
    )

    .validate_external_sample_size(
        N_pop2,
        "N_pop2"
    )


    if (!check_gwas_values) {
        warning(
            paste0(
                "GWAS numeric quality control has been disabled by the user. ",
                "The analysis will proceed without validating numeric GWAS values."
            ),
            call. = FALSE
        )
    }


    # ============================================================
    # 2. Confirm that all explicitly required input files exist
    # ============================================================

    .validate_input_file(
        bim_file,
        "bim_file"
    )

    .validate_input_file(
        map_file,
        "map_file"
    )

    .validate_input_file(
        assoc_pop1,
        "assoc_pop1"
    )

    .validate_input_file(
        assoc_pop2,
        "assoc_pop2"
    )

    .validate_input_directory(
        LD_path_pop1,
        "LD_path_pop1"
    )

    .validate_input_directory(
        LD_path_pop2,
        "LD_path_pop2"
    )


    # ============================================================
    # 3. Validate GCTA-specific inputs when requested
    # ============================================================

    gcta_exec <- NULL

    if (run_gcta) {

        gcta_exec <- .resolve_gcta_executable(
            gcta_bin
        )

        if (is.null(plink_bfile_pop1)) {
            stop(
                "'plink_bfile_pop1' must be supplied when 'run_gcta = TRUE'.",
                call. = FALSE
            )
        }

        if (is.null(plink_bfile_pop2)) {
            stop(
                "'plink_bfile_pop2' must be supplied when 'run_gcta = TRUE'.",
                call. = FALSE
            )
        }

        .validate_plink_bfile(
            plink_bfile_pop1,
            "plink_bfile_pop1"
        )

        .validate_plink_bfile(
            plink_bfile_pop2,
            "plink_bfile_pop2"
        )
    }


    # ============================================================
    # 4. Validate LD directory structure
    # ============================================================

    paired_prefixes_pop1 <- .validate_ld_directory(
        LD_path_pop1,
        "POP1"
    )

    paired_prefixes_pop2 <- .validate_ld_directory(
        LD_path_pop2,
        "POP2"
    )


    # ============================================================
    # 5. Read and validate the structural format of all input files
    # ============================================================

    require_n_pop1 <- is.null(
        N_pop1
    )

    require_n_pop2 <- is.null(
        N_pop2
    )


    gwas1 <- .read_validate_gwas(
        path = assoc_pop1,
        source_label = "POP1 GWAS summary statistics",
        require_n = require_n_pop1
    )

    gwas2 <- .read_validate_gwas(
        path = assoc_pop2,
        source_label = "POP2 GWAS summary statistics",
        require_n = require_n_pop2
    )


    bim <- .read_validate_bim(
        bim_file
    )

    gene_map <- .read_validate_gene_map(
        map_file
    )


    # ============================================================
    # 6. Validate gene-to-LD matching before analysis
    #
    # Missing LD for an individual gene is permitted.
    #
    # Ambiguous LD matching is not permitted because the analysis
    # cannot determine which LD matrix should be used safely.
    # ============================================================

    genes_to_validate <- unique(
        as.character(
            gene_map[["gene.id"]]
        )
    )


    .validate_ld_gene_matching(
        ld_dir = LD_path_pop1,
        genes = genes_to_validate,
        pop_label = "POP1",
        paired_prefixes = paired_prefixes_pop1
    )


    .validate_ld_gene_matching(
        ld_dir = LD_path_pop2,
        genes = genes_to_validate,
        pop_label = "POP2",
        paired_prefixes = paired_prefixes_pop2
    )


    # ============================================================
    # 7. Return validated input objects
    # ============================================================

    invisible(
        list(
            gwas1 = gwas1,
            gwas2 = gwas2,
            bim = bim,
            gene_map = gene_map,
            require_n_pop1 = require_n_pop1,
            require_n_pop2 = require_n_pop2,
            gcta_exec = gcta_exec
        )
    )
}