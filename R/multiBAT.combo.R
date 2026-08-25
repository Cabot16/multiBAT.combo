#' multiBAT.combo: Joint Gene-based Association Test for Two Populations
#'
#' Performs multi-population gene-based association analysis using
#' population-specific GWAS summary statistics and LD matrices.
#'
#' Population-specific mBAT analyses and joint multiBAT analysis are attempted
#' independently for each gene. The full multiBAT-combo test is used only when
#' the joint multiBAT analysis and both population-specific mBAT analyses are
#' available.
#'
#' If the joint multiBAT analysis is unavailable, but both population-specific
#' mBAT analyses are available, the two population-specific P values are
#' combined using ACAT with sample-size-based weights.
#'
#' If either population-specific mBAT result is unavailable, the joint
#' multiBAT result is not used in the final combination. The analysis falls
#' back to the remaining available population-specific mBAT result. If only
#' one population-specific result is available, that P value is reported
#' directly as \code{P_multiBAT.combo}.
#'
#' Genes for which neither population-specific mBAT result is available are
#' not reported, irrespective of whether joint multiBAT was successful.
#'
#' When \code{run_gcta = FALSE}, population-specific mBAT statistics are
#' calculated internally. When \code{run_gcta = TRUE}, population-specific
#' mBAT statistics are obtained from GCTA. In both modes, GCTA
#' \code{P_mBATcombo} values are never used as inputs to the final
#' multiBAT-combo test.
#'
#' @param bim_file Path to a headerless six-column PLINK BIM file.
#' @param map_file Path to a headerless four-column gene annotation file
#'   containing chromosome, start, end, and gene identifier.
#' @param assoc_pop1 Path to population 1 GWAS summary statistics.
#' @param assoc_pop2 Path to population 2 GWAS summary statistics.
#' @param LD_path_pop1 Directory containing population 1 GCTB LD files.
#' @param LD_path_pop2 Directory containing population 2 GCTB LD files.
#' @param result_prefix Full path and filename prefix for output files.
#' @param extend Non-negative whole number specifying the upstream and
#'   downstream gene-region extension in base pairs. Defaults to
#'   \code{50000}.
#' @param prop Numeric value in \code{(0, 1]}. Proportion of eigen-variance
#'   retained in the mBAT calculation. Defaults to \code{0.9}.
#' @param af_threshold Numeric allele-frequency difference threshold in
#'   \code{[0, 1]}. Defaults to \code{0.2}.
#' @param run_gcta Logical. If \code{FALSE}, population-specific mBAT
#'   statistics are calculated internally. If \code{TRUE}, GCTA is used for
#'   the population-specific analyses. Defaults to \code{FALSE}.
#' @param print_all_p Logical. If \code{TRUE}, intermediate joint and
#'   population-specific test statistics and P values are included in the
#'   output. Defaults to \code{FALSE}. GCTA mBAT statistics and all three
#'   GCTA P values are always reported when \code{run_gcta = TRUE}.
#' @param check_gwas_values Logical. Whether numeric GWAS quality control is
#'   performed before analysis. Defaults to \code{TRUE}.
#' @param N_pop1 Optional externally supplied sample size for population 1.
#' @param N_pop2 Optional externally supplied sample size for population 2.
#' @param gcta_bin GCTA executable name available in PATH, or full path to
#'   the executable. Defaults to \code{"gcta64"}.
#' @param plink_bfile_pop1 PLINK reference prefix for population 1. Required
#'   when \code{run_gcta = TRUE}.
#' @param plink_bfile_pop2 PLINK reference prefix for population 2. Required
#'   when \code{run_gcta = TRUE}.
#'
#' @return Invisibly returns a list containing the final multipopulation
#'   results in element \code{multipop}.
#'
#' @importFrom data.table fread fwrite copy as.data.table data.table setnames setkey setorder rbindlist :=
#' @importFrom snpsettest map_snp_to_gene
#' @importFrom stats pnorm pchisq p.adjust setNames
#' @export
multiBAT.combo <- function(
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

    # Input validation

    validated <- validate_multibat_inputs(
        bim_file = bim_file,
        map_file = map_file,
        assoc_pop1 = assoc_pop1,
        assoc_pop2 = assoc_pop2,
        LD_path_pop1 = LD_path_pop1,
        LD_path_pop2 = LD_path_pop2,
        result_prefix = result_prefix,
        extend = extend,
        prop = prop,
        af_threshold = af_threshold,
        run_gcta = run_gcta,
        print_all_p = print_all_p,
        check_gwas_values = check_gwas_values,
        N_pop1 = N_pop1,
        N_pop2 = N_pop2,
        gcta_bin = gcta_bin,
        plink_bfile_pop1 = plink_bfile_pop1,
        plink_bfile_pop2 = plink_bfile_pop2
    )


    # Output and logging

    output_dir <- dirname(result_prefix)

    if (!dir.exists(output_dir)) {
        dir.create(
            output_dir,
            recursive = TRUE,
            showWarnings = FALSE
        )
    }

    if (!dir.exists(output_dir)) {
        stop(
            sprintf(
                "The output directory could not be created: %s",
                output_dir
            ),
            call. = FALSE
        )
    }


    log_file_path <- paste0(
        result_prefix,
        ".log"
    )

    log_conn <- file(
        log_file_path,
        "a"
    )


    log_to_file <- function(
        msg,
        level = "INFO",
        console_print = FALSE
    ) {

        line <- sprintf(
            "[%s] [%s] %s\n",
            Sys.time(),
            level,
            msg
        )

        cat(
            line,
            file = log_conn,
            append = TRUE
        )

        if (console_print) {
            cat(line)
        }
    }


    log_info <- function(
        msg,
        console_print = FALSE
    ) {
        log_to_file(
            msg,
            "INFO",
            console_print
        )
    }


    log_warn <- function(
        msg,
        console_print = FALSE
    ) {
        log_to_file(
            msg,
            "WARN",
            console_print
        )
    }


    on.exit(
        {
            if (isOpen(log_conn)) {
                close(log_conn)
            }
        },
        add = TRUE
    )


    log_info(
        "Starting multiBAT.combo analysis.",
        console_print = TRUE
    )

    log_info(
        "Input and parameter validation completed successfully."
    )


    if (!check_gwas_values) {
        log_warn(
            paste0(
                "GWAS numeric quality control was disabled by the user. ",
                "Numeric GWAS values were not validated before analysis."
            )
        )
    }


    # Load inputs

    gwas1 <- data.table::copy(
        validated$gwas1
    )

    gwas2 <- data.table::copy(
        validated$gwas2
    )

    bim <- data.table::copy(
        validated$bim
    )

    gene_annot_raw <- data.table::copy(
        validated$gene_map
    )


    bim[, chr := as.character(chr)]

    gene_annot_raw[
        ,
        chr := as.character(chr)
    ]


    # GWAS QC

    if (check_gwas_values) {

        n1_before <- nrow(gwas1)
        n2_before <- nrow(gwas2)


        gwas1 <- gwas_numeric_qc(
            gwas = gwas1,
            pop_label = "POP1",
            check_n = validated$require_n_pop1
        )


        gwas2 <- gwas_numeric_qc(
            gwas = gwas2,
            pop_label = "POP2",
            check_n = validated$require_n_pop2
        )


        log_info(
            sprintf(
                paste0(
                    "GWAS numeric quality control completed. ",
                    "POP1 retained %d of %d SNPs. ",
                    "POP2 retained %d of %d SNPs."
                ),
                nrow(gwas1),
                n1_before,
                nrow(gwas2),
                n2_before
            ),
            console_print = TRUE
        )
    }


    if (nrow(gwas1) == 0L) {
        stop(
            "No POP1 GWAS SNPs are available for analysis.",
            call. = FALSE
        )
    }

    if (nrow(gwas2) == 0L) {
        stop(
            "No POP2 GWAS SNPs are available for analysis.",
            call. = FALSE
        )
    }


    # Sample sizes

    get_population_sample_size <- function(
        gwas,
        external_n,
        pop_label
    ) {

        if (!is.null(external_n)) {
            return(
                as.numeric(external_n)
            )
        }


        if (!"N" %in% names(gwas)) {
            stop(
                sprintf(
                    paste0(
                        "%s GWAS summary statistics do not contain an N ",
                        "column and no external sample size was supplied."
                    ),
                    pop_label
                ),
                call. = FALSE
            )
        }


        n_values <- suppressWarnings(
            as.numeric(
                gwas[["N"]]
            )
        )


        n_value <- mean(
            n_values,
            na.rm = TRUE
        )


        if (
            !is.finite(n_value) ||
            n_value <= 0
        ) {
            stop(
                sprintf(
                    "A usable positive sample size could not be derived for %s.",
                    pop_label
                ),
                call. = FALSE
            )
        }


        n_value
    }


    n1 <- get_population_sample_size(
        gwas = gwas1,
        external_n = N_pop1,
        pop_label = "POP1"
    )

    n2 <- get_population_sample_size(
        gwas = gwas2,
        external_n = N_pop2,
        pop_label = "POP2"
    )


    source1 <- if (is.null(N_pop1)) {
        "GWAS N"
    } else {
        "externally supplied N"
    }

    source2 <- if (is.null(N_pop2)) {
        "GWAS N"
    } else {
        "externally supplied N"
    }


    log_info(
        sprintf(
            paste0(
                "Sample sizes used for weighting. ",
                "POP1=%.3f (%s). POP2=%.3f (%s)."
            ),
            n1,
            source1,
            n2,
            source2
        ),
        console_print = TRUE
    )


    # ACAT weights

    imbalance_rate <- abs(
        n1 - n2
    ) / (
        n1 + n2
    )


    w_joint <- 1 -
        imbalance_rate


    w_pop1 <- (
        1 - w_joint
    ) * (
        n1 /
            (
                n1 + n2
            )
    )


    w_pop2 <- (
        1 - w_joint
    ) * (
        n2 /
            (
                n1 + n2
            )
    )


    full_weights <- c(
        w_joint,
        w_pop1,
        w_pop2
    )



    singlepop_weights <- c(
        n1 / (n1 + n2),
        n2 / (n1 + n2)
    )


    log_info(
        sprintf(
            paste0(
                "Full multiBAT-combo weights. ",
                "Joint=%.6f. POP1=%.6f. POP2=%.6f."
            ),
            w_joint,
            w_pop1,
            w_pop2
        ),
        console_print = TRUE
    )


    # Helper functions

    find_col <- function(
        dat,
        candidates
    ) {

        idx <- match(
            tolower(candidates),
            tolower(names(dat))
        )

        idx <- idx[
            !is.na(idx)
        ]


        if (length(idx) == 0L) {
            return(
                NA_character_
            )
        }


        names(dat)[
            idx[1]
        ]
    }


    standardise_ld_info <- function(
        info,
        pop_label
    ) {

        info <- data.table::copy(
            data.table::as.data.table(
                info
            )
        )


        if (nrow(info) == 0L) {
            stop(
                sprintf(
                    "%s LD information contains no SNPs",
                    pop_label
                ),
                call. = FALSE
            )
        }


        snp_col <- find_col(
            info,
            c(
                "snp",
                "SNP",
                "ID",
                "id"
            )
        )

        a1_col <- find_col(
            info,
            c(
                "A1",
                "a1"
            )
        )

        a2_col <- find_col(
            info,
            c(
                "A2",
                "a2"
            )
        )

        pos_col <- find_col(
            info,
            c(
                "PhysPos",
                "physpos",
                "POS",
                "pos",
                "BP",
                "bp"
            )
        )

        a1freq_col <- find_col(
            info,
            c(
                "A1Freq",
                "a1freq",
                "AF1",
                "af1"
            )
        )

        a2freq_col <- find_col(
            info,
            c(
                "A2Freq",
                "a2freq",
                "AF2",
                "af2"
            )
        )


        if (
            is.na(snp_col) ||
            is.na(a1_col) ||
            is.na(a2_col)
        ) {
            stop(
                sprintf(
                    paste0(
                        "%s LD information must contain SNP identifier, ",
                        "A1, and A2 columns"
                    ),
                    pop_label
                ),
                call. = FALSE
            )
        }


        info[
            ,
            snp := as.character(
                get(snp_col)
            )
        ]

        info[
            ,
            A1 := as.character(
                get(a1_col)
            )
        ]

        info[
            ,
            A2 := as.character(
                get(a2_col)
            )
        ]


        if (!is.na(pos_col)) {

            info[
                ,
                PhysPos := suppressWarnings(
                    as.numeric(
                        get(pos_col)
                    )
                )
            ]

        } else {

            info[
                ,
                PhysPos := seq_len(.N)
            ]
        }


        if (!is.na(a1freq_col)) {

            info[
                ,
                A1Freq := suppressWarnings(
                    as.numeric(
                        get(a1freq_col)
                    )
                )
            ]

        } else if (!is.na(a2freq_col)) {

            info[
                ,
                A1Freq := 1 -
                    suppressWarnings(
                        as.numeric(
                            get(a2freq_col)
                        )
                    )
            ]

        } else {

            stop(
                sprintf(
                    "%s LD information must contain either A1Freq or A2Freq",
                    pop_label
                ),
                call. = FALSE
            )
        }


        invalid_snp <- (
            is.na(info[["snp"]]) |
            !nzchar(
                trimws(
                    info[["snp"]]
                )
            )
        )


        if (any(invalid_snp)) {
            stop(
                sprintf(
                    "%s LD information contains missing SNP identifiers",
                    pop_label
                ),
                call. = FALSE
            )
        }


        if (anyDuplicated(info[["snp"]]) > 0L) {
            stop(
                sprintf(
                    "%s LD information contains duplicated SNP identifiers",
                    pop_label
                ),
                call. = FALSE
            )
        }


        info
    }


    calculate_internal_mbat <- function(
        z,
        R
    ) {

        if (
            length(z) == 0L ||
            nrow(R) == 0L ||
            ncol(R) == 0L
        ) {
            stop(
                "no SNPs were available for mBAT",
                call. = FALSE
            )
        }


        if (
            nrow(R) != length(z) ||
            ncol(R) != length(z)
        ) {
            stop(
                "the LD matrix dimensions do not match the Z-score vector",
                call. = FALSE
            )
        }


        if (any(!is.finite(z))) {
            stop(
                "non-finite Z scores were detected",
                call. = FALSE
            )
        }


        eig <- eigen(
            R,
            symmetric = TRUE
        )


        posval <- eig$values[
            eig$values > 1e-10
        ]


        if (length(posval) == 0L) {
            stop(
                "no positive LD eigenvalues remained",
                call. = FALSE
            )
        }


        if (prop == 1) {

            k <- length(
                posval
            )

        } else {

            cum_var <- cumsum(
                posval /
                    sum(posval)
            )

            k <- which(
                cum_var >
                    as.numeric(prop)
            )[1]
        }


        sel <- seq_len(k)

        lambda_prop <- eig$values[
            sel
        ]

        U_prop <- eig$vectors[
            ,
            sel,
            drop = FALSE
        ]


        UPz_prop <- crossprod(
            U_prop,
            z
        )


        T_mBAT <- as.numeric(
            crossprod(
                UPz_prop,
                (
                    1 /
                        lambda_prop
                ) *
                    UPz_prop
            )
        )


        P_mBAT <- stats::pchisq(
            T_mBAT,
            df = k,
            lower.tail = FALSE
        )


        list(
            T = T_mBAT,
            P = P_mBAT,
            k = k
        )
    }


    attempt <- function(fun) {

        tryCatch(
            list(
                ok = TRUE,
                value = fun(),
                reason = NA_character_
            ),
            error = function(e) {

                list(
                    ok = FALSE,
                    value = NULL,
                    reason = conditionMessage(e)
                )
            }
        )
    }


    valid_p <- function(x) {

        length(x) == 1L &&
            !is.na(x) &&
            is.finite(x) &&
            x >= 0 &&
            x <= 1
    }


    pick_first_existing_col <- function(
        dt,
        candidates
    ) {

        hit <- candidates[
            candidates %in%
                names(dt)
        ]


        if (length(hit) == 0L) {
            return(
                NA_character_
            )
        }


        hit[1]
    }


    # SNP-gene mapping

    gene_annot_R <- data.table::copy(
        gene_annot_raw
    )


    gene_annot_R[
        ,
        start := pmax(
            1,
            start - extend
        )
    ]

    gene_annot_R[
        ,
        end := end + extend
    ]


    bim_map <- bim[
        ,
        .(
            id = snp,
            chr,
            pos
        )
    ]


    snp2gene_df <- snpsettest::map_snp_to_gene(
        info_snp = bim_map,
        info_gene = data.frame(
            chr = gene_annot_R$chr,
            start = gene_annot_R$start,
            end = gene_annot_R$end,
            gene.id = gene_annot_R$gene.id,
            stringsAsFactors = FALSE
        ),
        only_sets = FALSE
    )$map


    snp2gene <- data.table::as.data.table(
        snp2gene_df
    )[
        ,
        .(
            id,
            chr,
            pos,
            gene.id
        )
    ]


    snp_gene <- merge(
        bim[
            ,
            .(
                snp,
                chr,
                pos
            )
        ],
        snp2gene,
        by.x = c(
            "snp",
            "chr",
            "pos"
        ),
        by.y = c(
            "id",
            "chr",
            "pos"
        )
    )


    data.table::setnames(
        snp_gene,
        "gene.id",
        "gene"
    )


    snp_gene <- snp_gene[
        !is.na(gene) &
            !is.na(snp)
    ]


    if (nrow(snp_gene) == 0L) {
        stop(
            "No SNPs could be mapped to the supplied gene annotation.",
            call. = FALSE
        )
    }


    # Merge GWAS data

    m1 <- merge(
        snp_gene,
        gwas1,
        by = "snp"
    )

    m2 <- merge(
        snp_gene,
        gwas2,
        by = "snp"
    )


    data.table::setkey(
        m1,
        gene
    )

    data.table::setkey(
        m2,
        gene
    )


    m1_list <- split(
        m1,
        by = "gene",
        keep.by = FALSE
    )

    m2_list <- split(
        m2,
        by = "gene",
        keep.by = FALSE
    )


    log_info(
        sprintf(
            paste0(
                "Mapped GWAS genes. ",
                "POP1=%d. POP2=%d. ",
                "Mapped SNP-gene rows. POP1=%d. POP2=%d."
            ),
            length(m1_list),
            length(m2_list),
            nrow(m1),
            nrow(m2)
        ),
        console_print = TRUE
    )


    # Gene list

    genes_all <- union(
        names(m1_list),
        names(m2_list)
    )


    if (length(genes_all) == 0L) {
        stop(
            "No genes contained mapped GWAS data.",
            call. = FALSE
        )
    }


    # Gene-level preparation

    prepare_population_gene <- function(
        g,
        pop_label,
        g_list,
        ld_dir
    ) {

        if (!g %in% names(g_list)) {
            stop(
                sprintf(
                    "no %s mapped GWAS data were available",
                    pop_label
                ),
                call. = FALSE
            )
        }


        gdat <- g_list[[g]]


        if (nrow(gdat) == 0L) {
            stop(
                sprintf(
                    "no %s GWAS SNPs were available",
                    pop_label
                ),
                call. = FALSE
            )
        }


        if (anyDuplicated(gdat[["snp"]]) > 0L) {
            stop(
                sprintf(
                    paste0(
                        "duplicated SNP identifiers were detected in the ",
                        "%s gene-level GWAS data"
                    ),
                    pop_label
                ),
                call. = FALSE
            )
        }


        ld_obj <- ReadLDMbin(
            ld_dir = ld_dir,
            gene = g,
            pop_label = pop_label
        )


        if (is.null(ld_obj)) {
            stop(
                sprintf(
                    "the %s LD matrix was not available or could not be read",
                    pop_label
                ),
                call. = FALSE
            )
        }


        LD_raw <- ld_obj$LD

        info <- standardise_ld_info(
            ld_obj$info,
            pop_label
        )


        common <- intersect(
            gdat[["snp"]],
            info[["snp"]]
        )


        if (length(common) == 0L) {
            stop(
                sprintf(
                    paste0(
                        "no %s SNPs were shared between the GWAS ",
                        "and LD information"
                    ),
                    pop_label
                ),
                call. = FALSE
            )
        }


        if (
            nrow(LD_raw) != nrow(info) ||
            ncol(LD_raw) != nrow(info)
        ) {
            stop(
                sprintf(
                    paste0(
                        "the %s LD matrix dimensions did not match ",
                        "the LD information"
                    ),
                    pop_label
                ),
                call. = FALSE
            )
        }


        rownames(LD_raw) <- info[["snp"]]
        colnames(LD_raw) <- info[["snp"]]


        g_qc <- qc_gwas_op(
            gdat,
            data.frame(
                snp = common,
                stringsAsFactors = FALSE
            )
        )


        R_qc <- qc_LD_op(
            g_qc,
            LD_raw,
            data.frame(
                snp = common,
                stringsAsFactors = FALSE
            )
        )


        if (
            nrow(g_qc) == 0L ||
            nrow(R_qc) == 0L
        ) {
            stop(
                sprintf(
                    "no %s SNPs remained after GWAS-LD matching",
                    pop_label
                ),
                call. = FALSE
            )
        }


        pos_dt <- unique(
            snp_gene[
                gene == g &
                    snp %in%
                    g_qc[["snp"]],
                .(
                    snp,
                    pos
                )
            ]
        )


        if (anyDuplicated(pos_dt[["snp"]]) > 0L) {
            stop(
                sprintf(
                    paste0(
                        "duplicated SNP identifiers were detected in the ",
                        "%s gene-level SNP-position mapping"
                    ),
                    pop_label
                ),
                call. = FALSE
            )
        }


        data.table::setorder(
            pos_dt,
            pos
        )


        ordered_snps <- pos_dt[["snp"]]


        if (length(ordered_snps) == 0L) {
            stop(
                sprintf(
                    "no %s genomic positions were available",
                    pop_label
                ),
                call. = FALSE
            )
        }


        g_qc <- g_qc[
            match(
                ordered_snps,
                g_qc[["snp"]]
            )
        ]


        R_qc <- R_qc[
            ordered_snps,
            ordered_snps,
            drop = FALSE
        ]


        info_sub <- info[
            match(
                ordered_snps,
                info[["snp"]]
            )
        ]


        if (anyNA(info_sub[["snp"]])) {
            stop(
                sprintf(
                    "%s LD information could not be aligned to GWAS SNPs",
                    pop_label
                ),
                call. = FALSE
            )
        }


        df_b <- allele_flip_gwas_beta(
            g_qc,
            info_sub,
            no_AF = FALSE,
            af_threshold = af_threshold
        )


        if (
            is.null(df_b) ||
            nrow(df_b) == 0L
        ) {
            stop(
                sprintf(
                    paste0(
                        "no %s SNPs remained after allele alignment ",
                        "or allele-frequency filtering"
                    ),
                    pop_label
                ),
                call. = FALSE
            )
        }


        df_b <- data.table::as.data.table(
            df_b
        )


        kept_snps <- as.character(
            df_b[["snp"]]
        )


        if (anyDuplicated(kept_snps) > 0L) {
            stop(
                sprintf(
                    paste0(
                        "duplicated SNP identifiers remained in %s ",
                        "after allele alignment"
                    ),
                    pop_label
                ),
                call. = FALSE
            )
        }


        R_final <- R_qc[
            kept_snps,
            kept_snps,
            drop = FALSE
        ]


        z <- suppressWarnings(
            as.numeric(
                df_b[["beta"]]
            ) /
                as.numeric(
                    df_b[["se"]]
                )
        )


        if (any(!is.finite(z))) {
            stop(
                sprintf(
                    "non-finite %s Z scores were detected",
                    pop_label
                ),
                call. = FALSE
            )
        }


        p_values <- suppressWarnings(
            as.numeric(
                df_b[["P"]]
            )
        )


        finite_p <- is.finite(p_values) &
            p_values >= 0 &
            p_values <= 1


        if (!any(finite_p)) {
            stop(
                sprintf(
                    "no valid %s SNP P values were available",
                    pop_label
                ),
                call. = FALSE
            )
        }


        ranking_p <- ifelse(
            finite_p,
            p_values,
            Inf
        )


        top_index <- which.min(
            ranking_p
        )


        list(
            gwas = df_b,
            R = R_final,
            z = z,
            snps = kept_snps,
            no_snps = length(kept_snps),
            top_snp = df_b[["snp"]][top_index],
            top_p = p_values[top_index]
        )
    }


    # Gene-level analysis

    analyse_gene <- function(g) {

        pop1_prepare <- attempt(
            function() {
                prepare_population_gene(
                    g = g,
                    pop_label = "POP1",
                    g_list = m1_list,
                    ld_dir = LD_path_pop1
                )
            }
        )


        pop2_prepare <- attempt(
            function() {
                prepare_population_gene(
                    g = g,
                    pop_label = "POP2",
                    g_list = m2_list,
                    ld_dir = LD_path_pop2
                )
            }
        )


        T_mBAT_pop1 <- NA_real_
        P_mBAT_pop1 <- NA_real_
        No.Eigenvalues_pop1 <- NA_integer_

        T_mBAT_pop2 <- NA_real_
        P_mBAT_pop2 <- NA_real_
        No.Eigenvalues_pop2 <- NA_integer_


        pop1_test_reason <- if (pop1_prepare$ok) {
            NA_character_
        } else {
            pop1_prepare$reason
        }


        pop2_test_reason <- if (pop2_prepare$ok) {
            NA_character_
        } else {
            pop2_prepare$reason
        }



        if (
            !run_gcta &&
            pop1_prepare$ok
        ) {

            pop1_test <- attempt(
                function() {
                    calculate_internal_mbat(
                        z = pop1_prepare$value$z,
                        R = pop1_prepare$value$R
                    )
                }
            )


            if (pop1_test$ok) {

                T_mBAT_pop1 <- pop1_test$value$T
                P_mBAT_pop1 <- pop1_test$value$P
                No.Eigenvalues_pop1 <- pop1_test$value$k

            } else {

                pop1_test_reason <- pop1_test$reason
            }
        }


        if (
            !run_gcta &&
            pop2_prepare$ok
        ) {

            pop2_test <- attempt(
                function() {
                    calculate_internal_mbat(
                        z = pop2_prepare$value$z,
                        R = pop2_prepare$value$R
                    )
                }
            )


            if (pop2_test$ok) {

                T_mBAT_pop2 <- pop2_test$value$T
                P_mBAT_pop2 <- pop2_test$value$P
                No.Eigenvalues_pop2 <- pop2_test$value$k

            } else {

                pop2_test_reason <- pop2_test$reason
            }
        }



        No.SNPs_pop1 <- if (pop1_prepare$ok) {
            pop1_prepare$value$no_snps
        } else {
            NA_integer_
        }


        No.SNPs_pop2 <- if (pop2_prepare$ok) {
            pop2_prepare$value$no_snps
        } else {
            NA_integer_
        }


        TopSNP_pop1 <- if (pop1_prepare$ok) {
            pop1_prepare$value$top_snp
        } else {
            NA_character_
        }


        TopSNP_Pvalue_pop1 <- if (pop1_prepare$ok) {
            pop1_prepare$value$top_p
        } else {
            NA_real_
        }


        TopSNP_pop2 <- if (pop2_prepare$ok) {
            pop2_prepare$value$top_snp
        } else {
            NA_character_
        }


        TopSNP_Pvalue_pop2 <- if (pop2_prepare$ok) {
            pop2_prepare$value$top_p
        } else {
            NA_real_
        }



        common_snps <- character()

        if (
            pop1_prepare$ok &&
            pop2_prepare$ok
        ) {

            common_snps <- intersect(
                pop1_prepare$value$snps,
                pop2_prepare$value$snps
            )
        }


        No.SNPs_common <- if (
            pop1_prepare$ok &&
            pop2_prepare$ok
        ) {
            length(common_snps)
        } else {
            NA_integer_
        }



        T_multiBAT <- NA_real_
        P_multiBAT <- NA_real_

        SNP_start <- NA_real_
        SNP_end <- NA_real_

        No.Eigenvalues <- NA_integer_

        joint_reason <- NA_character_


        if (
            !pop1_prepare$ok ||
            !pop2_prepare$ok
        ) {

            joint_reason <- paste0(
                "population-specific data preparation was unavailable ",
                "for at least one population"
            )

        } else if (length(common_snps) == 0L) {

            joint_reason <- paste0(
                "no SNPs were shared between populations after allele ",
                "alignment and allele-frequency filtering"
            )

        } else {

            joint_attempt <- attempt(
                function() {

                    common_pos <- unique(
                        snp_gene[
                            gene == g &
                                snp %in%
                                common_snps,
                            .(
                                snp,
                                pos
                            )
                        ]
                    )


                    if (anyDuplicated(common_pos[["snp"]]) > 0L) {
                        stop(
                            paste0(
                                "duplicated SNP identifiers were detected ",
                                "in the final joint SNP set"
                            ),
                            call. = FALSE
                        )
                    }


                    data.table::setorder(
                        common_pos,
                        pos
                    )


                    ordered_common <- common_pos[["snp"]]


                    if (length(ordered_common) == 0L) {
                        stop(
                            paste0(
                                "no valid genomic positions were available ",
                                "for the final joint SNP set"
                            ),
                            call. = FALSE
                        )
                    }


                    common_positions <- suppressWarnings(
                        as.numeric(
                            common_pos[["pos"]]
                        )
                    )


                    if (any(!is.finite(common_positions))) {
                        stop(
                            paste0(
                                "invalid genomic positions were detected ",
                                "in the final joint SNP set"
                            ),
                            call. = FALSE
                        )
                    }


                    df1_joint <- pop1_prepare$value$gwas[
                        match(
                            ordered_common,
                            pop1_prepare$value$gwas[["snp"]]
                        )
                    ]


                    df2_joint <- pop2_prepare$value$gwas[
                        match(
                            ordered_common,
                            pop2_prepare$value$gwas[["snp"]]
                        )
                    ]


                    R1_joint <- pop1_prepare$value$R[
                        ordered_common,
                        ordered_common,
                        drop = FALSE
                    ]


                    R2_joint <- pop2_prepare$value$R[
                        ordered_common,
                        ordered_common,
                        drop = FALSE
                    ]


                    z1 <- suppressWarnings(
                        as.numeric(
                            df1_joint[["beta"]]
                        ) /
                            as.numeric(
                                df1_joint[["se"]]
                            )
                    )


                    z2 <- suppressWarnings(
                        as.numeric(
                            df2_joint[["beta"]]
                        ) /
                            as.numeric(
                                df2_joint[["se"]]
                            )
                    )


                    if (
                        any(!is.finite(z1)) ||
                        any(!is.finite(z2))
                    ) {
                        stop(
                            "non-finite Z scores were detected in the joint analysis",
                            call. = FALSE
                        )
                    }


                    m <- length(
                        ordered_common
                    )


                    zero_block <- matrix(
                        0,
                        nrow = m,
                        ncol = m
                    )


                    R_joint <- rbind(
                        cbind(
                            R1_joint,
                            zero_block
                        ),
                        cbind(
                            zero_block,
                            R2_joint
                        )
                    )


                    joint_gwas <- data.table::rbindlist(
                        list(
                            data.table::data.table(
                                snp = df1_joint[["snp"]],
                                P = 2 *
                                    stats::pnorm(
                                        -abs(z1)
                                    )
                            ),
                            data.table::data.table(
                                snp = df2_joint[["snp"]],
                                P = 2 *
                                    stats::pnorm(
                                        -abs(z2)
                                    )
                            )
                        )
                    )


                    eig <- eigen(
                        R_joint,
                        symmetric = TRUE
                    )


                    chr_g <- unique(
                        as.character(
                            snp_gene[
                                gene == g,
                                chr
                            ]
                        )
                    )


                    if (length(chr_g) != 1L) {
                        stop(
                            "the gene could not be assigned to one unique chromosome",
                            call. = FALSE
                        )
                    }


                    this_annot <- gene_annot_raw[
                        gene.id == g &
                            chr == chr_g
                    ]


                    if (nrow(this_annot) == 0L) {
                        stop(
                            "the gene annotation record could not be found",
                            call. = FALSE
                        )
                    }


                    data.table::setnames(
                        this_annot,
                        "gene.id",
                        "gene"
                    )


                    result <- mbat_SVD(
                        gwas_noqc_final = joint_gwas,
                        eig = eig,
                        prop = prop,
                        z = c(
                            z1,
                            z2
                        ),
                        chr_name = chr_g,
                        ID = g,
                        gene_annotate = this_annot,
                        top1 = pop1_prepare$value$top_snp,
                        top1_p = pop1_prepare$value$top_p,
                        top2 = pop2_prepare$value$top_snp,
                        top2_p = pop2_prepare$value$top_p,
                        no_snps_pop1 = pop1_prepare$value$no_snps,
                        no_snps_pop2 = pop2_prepare$value$no_snps,
                        no_snps_common = length(ordered_common),
                        common_snp_positions = common_positions
                    )


                    data.table::as.data.table(
                        result
                    )
                }
            )


            if (joint_attempt$ok) {

                joint_result <- joint_attempt$value


                T_multiBAT <- suppressWarnings(
                    as.numeric(
                        joint_result[["Chisq_mBAT"]][1]
                    )
                )


                P_multiBAT <- suppressWarnings(
                    as.numeric(
                        joint_result[["P_mBAT"]][1]
                    )
                )


                SNP_start <- suppressWarnings(
                    as.numeric(
                        joint_result[["SNP_start"]][1]
                    )
                )


                SNP_end <- suppressWarnings(
                    as.numeric(
                        joint_result[["SNP_end"]][1]
                    )
                )


                No.Eigenvalues <- suppressWarnings(
                    as.integer(
                        joint_result[["No.Eigenvalues"]][1]
                    )
                )

            } else {

                joint_reason <- joint_attempt$reason
            }
        }


        data.table::data.table(
            Gene = g,

            No.SNPs_pop1 = No.SNPs_pop1,
            No.SNPs_pop2 = No.SNPs_pop2,
            No.SNPs_common = No.SNPs_common,

            SNP_start = SNP_start,
            SNP_end = SNP_end,

            TopSNP_pop1 = TopSNP_pop1,
            TopSNP_Pvalue_pop1 = TopSNP_Pvalue_pop1,

            TopSNP_pop2 = TopSNP_pop2,
            TopSNP_Pvalue_pop2 = TopSNP_Pvalue_pop2,

            No.Eigenvalues = No.Eigenvalues,
            No.Eigenvalues_pop1 = No.Eigenvalues_pop1,
            No.Eigenvalues_pop2 = No.Eigenvalues_pop2,

            T_multiBAT = T_multiBAT,
            P_multiBAT = P_multiBAT,

            T_mBAT_pop1 = T_mBAT_pop1,
            P_mBAT_pop1 = P_mBAT_pop1,

            T_mBAT_pop2 = T_mBAT_pop2,
            P_mBAT_pop2 = P_mBAT_pop2,

            reason_pop1 = pop1_test_reason,
            reason_pop2 = pop2_test_reason,
            reason_joint = joint_reason
        )
    }


    log_info(
        sprintf(
            "Preparing population-specific and joint data for %d genes.",
            length(genes_all)
        ),
        console_print = TRUE
    )


   # Run gene analysis

total_genes <- length(genes_all)

start_time <- Sys.time()

last_report_time <- start_time


gene_records <- vector(
    "list",
    total_genes
)


for (i in seq_along(genes_all)) {

    gene_records[[i]] <- analyse_gene(
        genes_all[[i]]
    )


    current_time <- Sys.time()


    elapsed_seconds <- as.numeric(
        difftime(
            current_time,
            start_time,
            units = "secs"
        )
    )


    report_interval <- as.numeric(
        difftime(
            current_time,
            last_report_time,
            units = "secs"
        )
    )


    if (
        report_interval >= 1800 ||
        i == total_genes
    ) {

        completed_pct <- i /
            total_genes *
            100


        elapsed_time <- sprintf(
            "%02d:%02d:%02d",
            floor(elapsed_seconds / 3600),
            floor(
                elapsed_seconds %% 3600 / 60
            ),
            floor(
                elapsed_seconds %% 60
            )
        )


        log_info(
            sprintf(
                paste0(
                    "Gene processing progress: ",
                    "%d/%d (%.2f%%). ",
                    "Elapsed time=%s."
                ),
                i,
                total_genes,
                completed_pct,
                elapsed_time
            ),
            console_print = TRUE
        )


        last_report_time <- current_time
    }
}


    result_dt <- data.table::rbindlist(
        gene_records,
        fill = TRUE
    )


    # GCTA analyses

    prepare_gcta_result <- function(
        x,
        pop_label
    ) {

        if (is.null(x)) {
            return(
                data.table::data.table()
            )
        }


        dt <- data.table::copy(
            data.table::as.data.table(
                x
            )
        )


        if (nrow(dt) == 0L) {
            return(
                data.table::data.table()
            )
        }


        gene_col <- pick_first_existing_col(
            dt,
            c(
                "Gene",
                "gene",
                "GENE"
            )
        )


        t_mbat_col <- pick_first_existing_col(
            dt,
            c(
                "Chisq_mBAT",
                "Chisq.mBAT",
                "Chisq_mbat",
                "T_mBAT"
            )
        )


        p_mbat_col <- pick_first_existing_col(
            dt,
            c(
                "P_mBAT",
                "P_mbat"
            )
        )


        p_fastbat_col <- pick_first_existing_col(
            dt,
            c(
                "P_fastBAT",
                "P_fastbat",
                "P_FastBAT"
            )
        )


        p_combo_col <- pick_first_existing_col(
            dt,
            c(
                "P_mBATcombo",
                "P_mBAT_combo",
                "P_mBAT-combo",
                "P_mbatcombo"
            )
        )


        required_columns <- c(
            gene_col,
            t_mbat_col,
            p_mbat_col,
            p_fastbat_col,
            p_combo_col
        )


        if (anyNA(required_columns)) {
            stop(
                sprintf(
                    paste0(
                        "GCTA %s output does not contain all required ",
                        "mBAT result columns. Available columns are %s."
                    ),
                    pop_label,
                    paste(
                        names(dt),
                        collapse = ", "
                    )
                ),
                call. = FALSE
            )
        }


        nsnp_col <- pick_first_existing_col(
            dt,
            c(
                "No.SNPs",
                "No.SNP",
                "NSNP",
                "NSnps",
                "nSNP",
                "nSNPs"
            )
        )


        eigen_col <- pick_first_existing_col(
            dt,
            c(
                "No.Eigenvalues",
                "No_Eigenvalues",
                "NoEigenvalues"
            )
        )


        top_col <- pick_first_existing_col(
            dt,
            c(
                "TopSNP",
                "topSNP",
                "Top_SNP"
            )
        )


        top_p_col <- pick_first_existing_col(
            dt,
            c(
                "TopSNP_Pvalue",
                "TopSNP_P",
                "TopSNP_Pval",
                "topSNP_P"
            )
        )


        get_numeric <- function(column_name) {

            if (is.na(column_name)) {
                return(
                    rep(
                        NA_real_,
                        nrow(dt)
                    )
                )
            }


            suppressWarnings(
                as.numeric(
                    dt[[column_name]]
                )
            )
        }


        get_character <- function(column_name) {

            if (is.na(column_name)) {
                return(
                    rep(
                        NA_character_,
                        nrow(dt)
                    )
                )
            }


            as.character(
                dt[[column_name]]
            )
        }


        out <- data.table::data.table(
            Gene = as.character(
                dt[[gene_col]]
            ),

            No.SNPs = get_numeric(
                nsnp_col
            ),

            No.Eigenvalues = get_numeric(
                eigen_col
            ),

            TopSNP = get_character(
                top_col
            ),

            TopSNP_Pvalue = get_numeric(
                top_p_col
            ),

            T_mBAT = get_numeric(
                t_mbat_col
            ),

            P_mBAT = get_numeric(
                p_mbat_col
            ),

            P_fastBAT = get_numeric(
                p_fastbat_col
            ),

            P_mBATcombo = get_numeric(
                p_combo_col
            )
        )


        out <- out[
            !is.na(Gene) &
                nzchar(
                    trimws(Gene)
                )
        ]


        if (anyDuplicated(out[["Gene"]]) > 0L) {

            log_warn(
                sprintf(
                    paste0(
                        "GCTA %s output contains duplicated Gene entries. ",
                        "The row with the minimum finite P_mBAT was retained."
                    ),
                    pop_label
                )
            )


            rank_p <- ifelse(
                is.finite(
                    out[["P_mBAT"]]
                ) &
                    out[["P_mBAT"]] >= 0 &
                    out[["P_mBAT"]] <= 1,
                out[["P_mBAT"]],
                Inf
            )


            gene_indices <- split(
                seq_len(
                    nrow(out)
                ),
                out[["Gene"]]
            )


            keep_indices <- vapply(
                gene_indices,
                function(idx) {
                    idx[
                        which.min(
                            rank_p[idx]
                        )
                    ]
                },
                integer(1)
            )


            out <- out[
                sort(
                    keep_indices
                )
            ]
        }


        out
    }


    gcta_pop1 <- data.table::data.table()
    gcta_pop2 <- data.table::data.table()


    if (run_gcta) {

        gcta_wind_kb <- extend /
            1000


        log_info(
            "Running GCTA mBAT-combo for POP1.",
            console_print = TRUE
        )


        gcta_pop1_attempt <- attempt(
            function() {

                raw <- run_gcta_mbat(
                    sumstats = gwas1,
                    gene_map = gene_annot_raw,
                    bfile_path = plink_bfile_pop1,
                    out_file = paste0(
                        result_prefix,
                        ".gcta.pop1"
                    ),
                    gcta_bin = gcta_bin,
                    af_threshold = af_threshold,
                    gcta_wind_kb = gcta_wind_kb,
                    sample_size = N_pop1
                )


                prepare_gcta_result(
                    raw,
                    "POP1"
                )
            }
        )


        if (gcta_pop1_attempt$ok) {

            gcta_pop1 <- gcta_pop1_attempt$value

        } else {

            log_warn(
                sprintf(
                    "POP1 GCTA mBAT analysis failed because %s.",
                    gcta_pop1_attempt$reason
                ),
                console_print = TRUE
            )
        }


        log_info(
            "Running GCTA mBAT-combo for POP2.",
            console_print = TRUE
        )


        gcta_pop2_attempt <- attempt(
            function() {

                raw <- run_gcta_mbat(
                    sumstats = gwas2,
                    gene_map = gene_annot_raw,
                    bfile_path = plink_bfile_pop2,
                    out_file = paste0(
                        result_prefix,
                        ".gcta.pop2"
                    ),
                    gcta_bin = gcta_bin,
                    af_threshold = af_threshold,
                    gcta_wind_kb = gcta_wind_kb,
                    sample_size = N_pop2
                )


                prepare_gcta_result(
                    raw,
                    "POP2"
                )
            }
        )


        if (gcta_pop2_attempt$ok) {

            gcta_pop2 <- gcta_pop2_attempt$value

        } else {

            log_warn(
                sprintf(
                    "POP2 GCTA mBAT analysis failed because %s.",
                    gcta_pop2_attempt$reason
                ),
                console_print = TRUE
            )
        }



        if (nrow(gcta_pop1) > 0L) {

            data.table::setnames(
                gcta_pop1,
                old = c(
                    "No.SNPs",
                    "No.Eigenvalues",
                    "TopSNP",
                    "TopSNP_Pvalue",
                    "T_mBAT",
                    "P_mBAT",
                    "P_fastBAT",
                    "P_mBATcombo"
                ),
                new = c(
                    "GCTA_No.SNPs_pop1",
                    "No.Eigenvalues_pop1",
                    "GCTA_TopSNP_pop1",
                    "GCTA_TopSNP_Pvalue_pop1",
                    "GCTA_T_mBAT_pop1",
                    "GCTA_P_mBAT_pop1",
                    "P_fastBAT_pop1",
                    "P_mBATcombo_pop1"
                )
            )
        }


        if (nrow(gcta_pop2) > 0L) {

            data.table::setnames(
                gcta_pop2,
                old = c(
                    "No.SNPs",
                    "No.Eigenvalues",
                    "TopSNP",
                    "TopSNP_Pvalue",
                    "T_mBAT",
                    "P_mBAT",
                    "P_fastBAT",
                    "P_mBATcombo"
                ),
                new = c(
                    "GCTA_No.SNPs_pop2",
                    "No.Eigenvalues_pop2",
                    "GCTA_TopSNP_pop2",
                    "GCTA_TopSNP_Pvalue_pop2",
                    "GCTA_T_mBAT_pop2",
                    "GCTA_P_mBAT_pop2",
                    "P_fastBAT_pop2",
                    "P_mBATcombo_pop2"
                )
            )
        }



        final_genes <- union(
            result_dt[["Gene"]],
            union(
                if (nrow(gcta_pop1) > 0L) {
                    gcta_pop1[["Gene"]]
                } else {
                    character()
                },
                if (nrow(gcta_pop2) > 0L) {
                    gcta_pop2[["Gene"]]
                } else {
                    character()
                }
            )
        )


        combined <- data.table::data.table(
            Gene = final_genes
        )


        combined <- merge(
            combined,
            result_dt,
            by = "Gene",
            all.x = TRUE,
            sort = FALSE
        )


        if (nrow(gcta_pop1) > 0L) {

            combined <- merge(
                combined,
                gcta_pop1,
                by = "Gene",
                all.x = TRUE,
                sort = FALSE
            )
        }


        if (nrow(gcta_pop2) > 0L) {

            combined <- merge(
                combined,
                gcta_pop2,
                by = "Gene",
                all.x = TRUE,
                sort = FALSE
            )
        }



        if (!"GCTA_T_mBAT_pop1" %in% names(combined)) {
            combined[
                ,
                GCTA_T_mBAT_pop1 := NA_real_
            ]
        }

        if (!"GCTA_P_mBAT_pop1" %in% names(combined)) {
            combined[
                ,
                GCTA_P_mBAT_pop1 := NA_real_
            ]
        }

        if (!"GCTA_T_mBAT_pop2" %in% names(combined)) {
            combined[
                ,
                GCTA_T_mBAT_pop2 := NA_real_
            ]
        }

        if (!"GCTA_P_mBAT_pop2" %in% names(combined)) {
            combined[
                ,
                GCTA_P_mBAT_pop2 := NA_real_
            ]
        }


        if (!"No.Eigenvalues_pop1" %in% names(combined)) {
            combined[
                ,
                No.Eigenvalues_pop1 := NA_integer_
            ]
        }


        if (!"No.Eigenvalues_pop2" %in% names(combined)) {
            combined[
                ,
                No.Eigenvalues_pop2 := NA_integer_
            ]
        }


        combined[
            ,
            T_mBAT_pop1 :=
                GCTA_T_mBAT_pop1
        ]

        combined[
            ,
            P_mBAT_pop1 :=
                GCTA_P_mBAT_pop1
        ]

        combined[
            ,
            T_mBAT_pop2 :=
                GCTA_T_mBAT_pop2
        ]

        combined[
            ,
            P_mBAT_pop2 :=
                GCTA_P_mBAT_pop2
        ]



        if ("GCTA_TopSNP_pop1" %in% names(combined)) {

            use_gcta <- !is.na(
                combined[["GCTA_TopSNP_pop1"]]
            )

            combined[
                use_gcta,
                TopSNP_pop1 :=
                    GCTA_TopSNP_pop1
            ]

            combined[
                use_gcta,
                TopSNP_Pvalue_pop1 :=
                    GCTA_TopSNP_Pvalue_pop1
            ]
        }


        if ("GCTA_TopSNP_pop2" %in% names(combined)) {

            use_gcta <- !is.na(
                combined[["GCTA_TopSNP_pop2"]]
            )

            combined[
                use_gcta,
                TopSNP_pop2 :=
                    GCTA_TopSNP_pop2
            ]

            combined[
                use_gcta,
                TopSNP_Pvalue_pop2 :=
                    GCTA_TopSNP_Pvalue_pop2
            ]
        }


        result_dt <- combined
    }


    # Result columns

    required_result_columns <- c(
        "Gene",
        "No.SNPs_pop1",
        "No.SNPs_pop2",
        "No.SNPs_common",
        "SNP_start",
        "SNP_end",
        "TopSNP_pop1",
        "TopSNP_Pvalue_pop1",
        "TopSNP_pop2",
        "TopSNP_Pvalue_pop2",
        "No.Eigenvalues",
        "No.Eigenvalues_pop1",
        "No.Eigenvalues_pop2",
        "T_multiBAT",
        "P_multiBAT",
        "T_mBAT_pop1",
        "P_mBAT_pop1",
        "T_mBAT_pop2",
        "P_mBAT_pop2"
    )


    character_columns <- c(
        "Gene",
        "TopSNP_pop1",
        "TopSNP_pop2"
    )


    for (column_name in required_result_columns) {

        if (!column_name %in% names(result_dt)) {

            if (column_name %in% character_columns) {

                result_dt[
                    ,
                    (column_name) :=
                        NA_character_
                ]

            } else {

                result_dt[
                    ,
                    (column_name) :=
                        NA_real_
                ]
            }
        }
    }


    # Combination fallback

    combine_gene <- function(
        p_joint,
        p_pop1,
        p_pop2
    ) {

        has_joint <- valid_p(
            p_joint
        )

        has_pop1 <- valid_p(
            p_pop1
        )

        has_pop2 <- valid_p(
            p_pop2
        )


        if (
            has_joint &&
            has_pop1 &&
            has_pop2
        ) {

            p_values <- c(
                p_joint,
                p_pop1,
                p_pop2
            )

            weights <- full_weights


            keep_weight <- weights > 0


            return(
                ACATO(
                    p = p_values[
                        keep_weight
                    ],
                    w = weights[
                        keep_weight
                    ]
                )
            )
        }


        if (
            has_pop1 &&
            has_pop2
        ) {

            return(
                ACATO(
                    p = c(
                        p_pop1,
                        p_pop2
                    ),
                    w = singlepop_weights
                )
            )
        }


        if (has_pop1) {
            return(
                p_pop1
            )
        }


        if (has_pop2) {
            return(
                p_pop2
            )
        }


        NA_real_
    }


    result_dt[
        ,
        P_multiBAT.combo :=
            vapply(
                seq_len(.N),
                function(i) {

                    combine_gene(
                        p_joint = P_multiBAT[i],
                        p_pop1 = P_mBAT_pop1[i],
                        p_pop2 = P_mBAT_pop2[i]
                    )
                },
                numeric(1)
            )
    ]


    # Record status

    result_dt[
        ,
        analysis_status :=
            vapply(
                seq_len(.N),
                function(i) {

                    has_joint <- valid_p(
                        P_multiBAT[i]
                    )

                    has_pop1 <- valid_p(
                        P_mBAT_pop1[i]
                    )

                    has_pop2 <- valid_p(
                        P_mBAT_pop2[i]
                    )


                    if (
                        has_joint &&
                        has_pop1 &&
                        has_pop2
                    ) {
                        return(
                            "full"
                        )
                    }


                    if (
                        has_pop1 &&
                        has_pop2
                    ) {
                        return(
                            "two_population_fallback"
                        )
                    }


                    if (has_pop1) {
                        return(
                            "pop1_fallback"
                        )
                    }


                    if (has_pop2) {
                        return(
                            "pop2_fallback"
                        )
                    }


                    "no_result"
                },
                character(1)
            )
    ]


    # Fallback summary

    fallback_two <- result_dt[
        analysis_status ==
            "two_population_fallback",
        Gene
    ]


    fallback_pop1 <- result_dt[
        analysis_status ==
            "pop1_fallback",
        Gene
    ]


    fallback_pop2 <- result_dt[
        analysis_status ==
            "pop2_fallback",
        Gene
    ]


    no_result_genes <- result_dt[
        analysis_status ==
            "no_result",
        Gene
    ]


    if (length(fallback_two) > 0L) {

        for (g in fallback_two) {

            log_warn(
                sprintf(
                    paste0(
                        "Gene %s used the two-population mBAT fallback ",
                        "because the full multiBAT-combo analysis was unavailable."
                    ),
                    g
                )
            )
        }
    }


    if (length(fallback_pop1) > 0L) {

        for (g in fallback_pop1) {

            log_warn(
                sprintf(
                    paste0(
                        "Gene %s used the POP1 mBAT fallback because the ",
                        "POP2 mBAT result was unavailable."
                    ),
                    g
                )
            )
        }
    }


    if (length(fallback_pop2) > 0L) {

        for (g in fallback_pop2) {

            log_warn(
                sprintf(
                    paste0(
                        "Gene %s used the POP2 mBAT fallback because the ",
                        "POP1 mBAT result was unavailable."
                    ),
                    g
                )
            )
        }
    }


    if (length(no_result_genes) > 0L) {

        for (g in no_result_genes) {

            log_warn(
                sprintf(
                    paste0(
                        "Gene %s was skipped because no population-specific ",
                        "mBAT result was available."
                    ),
                    g
                )
            )
        }
    }


    # Remove unusable genes

    result_dt <- result_dt[
        analysis_status !=
            "no_result"
    ]


    if (nrow(result_dt) == 0L) {
        stop(
            paste0(
                "No genes had a usable population-specific mBAT result. ",
                "No multiBAT-combo results could be generated."
            ),
            call. = FALSE
        )
    }


    # Multiple testing correction

    result_dt[
        ,
        P_multiBAT.combo_BH :=
            stats::p.adjust(
                P_multiBAT.combo,
                method = "BH"
            )
    ]


    # Add annotation

    gene_meta <- gene_annot_raw[
        ,
        .(
            Chr = as.character(
                chr[1]
            ),
            Start = start[1],
            End = end[1]
        ),
        by = gene.id
    ]


    data.table::setnames(
        gene_meta,
        "gene.id",
        "Gene"
    )


    result_dt <- merge(
        result_dt,
        gene_meta,
        by = "Gene",
        all.x = TRUE,
        sort = FALSE
    )


    # Output structure

    base_columns <- c(
        "Gene",
        "Chr",
        "Start",
        "End",

        "No.SNPs_pop1",
        "No.SNPs_pop2",
        "No.SNPs_common",

        "SNP_start",
        "SNP_end",

        "TopSNP_pop1",
        "TopSNP_Pvalue_pop1",

        "TopSNP_pop2",
        "TopSNP_Pvalue_pop2",

        "No.Eigenvalues"
    )


    intermediate_columns <- character()


    if (!run_gcta) {

        if (print_all_p) {

            intermediate_columns <- c(
                "T_multiBAT",
                "P_multiBAT",

                "T_mBAT_pop1",
                "P_mBAT_pop1",

                "T_mBAT_pop2",
                "P_mBAT_pop2",

                "No.Eigenvalues_pop1",
                "No.Eigenvalues_pop2"
            )
        }

    } else {

        if (print_all_p) {

            intermediate_columns <- c(
                "T_multiBAT",
                "P_multiBAT",

                "No.Eigenvalues_pop1",
                "No.Eigenvalues_pop2"
            )
        }


        gcta_columns <- c(
            "T_mBAT_pop1",
            "P_mBAT_pop1",
            "P_fastBAT_pop1",
            "P_mBATcombo_pop1",

            "T_mBAT_pop2",
            "P_mBAT_pop2",
            "P_fastBAT_pop2",
            "P_mBATcombo_pop2"
        )


        for (column_name in gcta_columns) {

            if (!column_name %in% names(result_dt)) {

                result_dt[
                    ,
                    (column_name) :=
                        NA_real_
                ]
            }
        }


        intermediate_columns <- c(
            intermediate_columns,
            gcta_columns
        )
    }


    final_columns <- c(
        base_columns,
        intermediate_columns,
        "P_multiBAT.combo",
        "P_multiBAT.combo_BH"
    )


    # Remove internal columns

    output_dt <- result_dt[
        ,
        final_columns,
        with = FALSE
    ]


    # Write output

    multipop_file <- paste0(
        result_prefix,
        ".multipop.txt"
    )


    data.table::fwrite(
        output_dt,
        file = multipop_file,
        sep = "\t",
        quote = FALSE,
        na = "NA"
    )


    # Summary

    n_full <- sum(
        result_dt[["analysis_status"]] ==
            "full"
    )

    n_two_fallback <- sum(
        result_dt[["analysis_status"]] ==
            "two_population_fallback"
    )

    n_pop1_fallback <- sum(
        result_dt[["analysis_status"]] ==
            "pop1_fallback"
    )

    n_pop2_fallback <- sum(
        result_dt[["analysis_status"]] ==
            "pop2_fallback"
    )


    log_info(
        sprintf(
            paste0(
                "Analysis completed. ",
                "Full multiBAT-combo=%d. ",
                "Two-population fallback=%d. ",
                "POP1-only fallback=%d. ",
                "POP2-only fallback=%d. ",
                "Total reported genes=%d."
            ),
            n_full,
            n_two_fallback,
            n_pop1_fallback,
            n_pop2_fallback,
            nrow(output_dt)
        ),
        console_print = TRUE
    )


    log_info(
        sprintf(
            "Final results were written to %s.",
            multipop_file
        ),
        console_print = TRUE
    )


    invisible(
        list(
            multipop = output_dt
        )
    )
}