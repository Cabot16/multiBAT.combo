#' multiBAT.combo: Joint gene-based association test for two populations
#'
#' Reads GWAS summaries and per-gene LD matrices, applies allele flipping, then
#' wraps up with \code{mbat_SVD()} for the chi-square statistic and p-value.
#' Includes optional GCTA mBAT-combo integration and multiBAT-combo ACAT analysis.
#'
#' @param bim_file Path to PLINK *.bim file with columns chr, snp, bp, pos, A1, A2.
#' @param map_file Gene annotation file without header: chr, start, end, gene.id.
#' @param assoc_pop1 GWAS summary for pop1 with columns snp, A1, A2, AF1, beta, se, P, N.
#' @param assoc_pop2 GWAS summary for pop2 with the same format.
#' @param LD_path_pop1 Directory of GCTB LD files for pop1, including .bin and .info files.
#' @param LD_path_pop2 Directory of GCTB LD files for pop2, including .bin and .info files.
#' @param result_prefix Full path and base name for output files.
#' @param extend Integer. Gene region extension length in base pairs applied upstream
#'   and downstream for both R-internal mapping and GCTA calls. Defaults to
#'   \code{50000}, i.e. 50 kb.
#' @param prop Proportion of eigen-variance to retain. Defaults to \code{0.9}.
#' @param af_threshold Allele frequency difference threshold for filtering. Defaults to \code{0.2}.
#' @param run_gcta Logical. If \code{TRUE}, runs GCTA mBAT-combo for each population.
#'   Defaults to \code{FALSE}.
#' @param gcta_bin Path to GCTA executable, or executable name available in PATH.
#'   Defaults to \code{"gcta64"}.
#' @param plink_bfile_pop1 Path prefix to PLINK reference panel for pop1.
#' @param plink_bfile_pop2 Path prefix to PLINK reference panel for pop2.
#'
#' @return A list containing multipop and singlepop results.
#'
#' @importFrom data.table fread setnames fwrite rbindlist copy dcast merge.data.table setcolorder := as.data.table setkey setorder data.table
#' @importFrom snpsettest map_snp_to_gene
#' @importFrom stats na.omit pnorm pchisq setNames
#' @importFrom utils tail
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
    gcta_bin = "gcta64",
    plink_bfile_pop1 = NULL,
    plink_bfile_pop2 = NULL
) {

    # --- Log setup and functions ---
    log_file_path <- paste0(result_prefix, ".log")
    log_conn <- file(log_file_path, "a")

    log_to_file <- function(msg, level = "INFO", console_print = FALSE) {
        line <- sprintf("[%s] [%s] %s\n", Sys.time(), level, msg)
        cat(line, file = log_conn, append = TRUE)
        if (console_print) {
            cat(line)
        }
    }

    log_info <- function(msg, console_print = FALSE) log_to_file(msg, "INFO", console_print)
    log_warn <- function(msg, console_print = FALSE) log_to_file(msg, "WARN", console_print)
    log_error <- function(msg, console_print = FALSE) log_to_file(msg, "ERROR", console_print)

    log_info("Starting multiBAT.combo", console_print = TRUE)

    # Initialize counters
    stats <- list(
        multipop_success = 0,
        multipop_ld_missing = 0,
        multipop_qc_fail = 0,
        singlepop_success = 0,
        singlepop_ld_missing = 0,
        singlepop_qc_fail = 0
    )

    on.exit({
        if (isOpen(log_conn)) close(log_conn)
    })

    ## 1. Load gene annotation & GWAS
    log_info(sprintf("Reading map_file: %s", map_file))
    gene_annot_raw <- fread(map_file, header = FALSE)
    setnames(gene_annot_raw, c("chr", "start", "end", "gene.id"))

    gene_annot_R <- copy(gene_annot_raw)
    gene_annot_R[, start := pmax(1, start - extend)]
    gene_annot_R[, end := end + extend]

    log_info("Reading GWAS summaries...")
    gwas1 <- fread(assoc_pop1)
    setnames(gwas1, c("snp", "A1", "A2", "AF1", "beta", "se", "P", "N"))

    gwas2 <- fread(assoc_pop2)
    setnames(gwas2, names(gwas1))

    bim <- fread(bim_file)
    setnames(bim, c("chr", "snp", "bp", "pos", "A1", "A2"))
    bim[, chr := as.character(chr)]

    pop_weights <- NULL
    if (run_gcta) {
        n1 <- if ("N" %in% names(gwas1)) mean(gwas1$N, na.rm = TRUE) else NA_real_
        n2 <- if ("N" %in% names(gwas2)) mean(gwas2$N, na.rm = TRUE) else NA_real_

        if (!is.na(n1) && !is.na(n2) && n1 > 0 && n2 > 0) {
            ir <- abs(n1 - n2) / (n1 + n2)
            w_joint <- 1 - ir
            w_pop1 <- (1 - w_joint) * (n1 / (n1 + n2))
            w_pop2 <- (1 - w_joint) * (n2 / (n1 + n2))
            pop_weights <- c(w_joint, w_pop1, w_pop2)
            log_info(sprintf("Weights: Joint=%.3f, Pop1=%.3f, Pop2=%.3f", w_joint, w_pop1, w_pop2))
        }
    }

    ## 2. Map SNP -> gene
    bim_map <- bim[, .(id = snp, chr, pos)]
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

    snp2gene <- as.data.table(snp2gene_df)[, .(id, chr, pos, gene.id)]
    snp_gene <- merge(
        bim[, .(snp, chr, pos)],
        snp2gene,
        by.x = c("snp", "chr", "pos"),
        by.y = c("id", "chr", "pos")
    )
    setnames(snp_gene, "gene.id", "gene")
    snp_gene <- na.omit(snp_gene)

    ## 3. Pre-merge GWAS & split by gene
    m1 <- merge(snp_gene, gwas1, by = "snp")
    m2 <- merge(snp_gene, gwas2, by = "snp")
    setkey(m1, gene)
    setkey(m2, gene)

    m1_list <- split(m1, by = "gene", keep.by = FALSE)
    m2_list <- split(m2, by = "gene", keep.by = FALSE)

    ## 4. Pre-scan LD .bin files
    bins1_all <- list.files(LD_path_pop1, pattern = "\\.bin$", full.names = TRUE)
    bins2_all <- list.files(LD_path_pop2, pattern = "\\.bin$", full.names = TRUE)

    ## 5. Identify valid genes (Multipop)
    genes1 <- intersect(names(m1_list), names(m2_list))

    gene_prefix1 <- setNames(sapply(genes1, function(g) {
        f <- bins1_all[grep(paste0("^", g, "\\."), basename(bins1_all), perl = TRUE)]
        if (length(f)) sub("\\.bin$", "", f[1]) else NA_character_
    }, USE.NAMES = TRUE), genes1)

    gene_prefix2 <- setNames(sapply(genes1, function(g) {
        f <- bins2_all[grep(paste0("^", g, "\\."), basename(bins2_all), perl = TRUE)]
        if (length(f)) sub("\\.bin$", "", f[1]) else NA_character_
    }, USE.NAMES = TRUE), genes1)

    valid_genes <- genes1[!is.na(gene_prefix1) & !is.na(gene_prefix2)]
    stats$multipop_ld_missing <- length(genes1[is.na(gene_prefix1) | is.na(gene_prefix2)])

    ## 6. Worker: process one gene (multipop)
    process_gene <- function(g) {
        prefix1 <- gene_prefix1[[g]]
        prefix2 <- gene_prefix2[[g]]
        g1 <- m1_list[[g]]
        g2 <- m2_list[[g]]

        skip_gene <- function(reason) {
            log_warn(sprintf("Skipping gene %s: %s", g, reason))
            stats$multipop_qc_fail <<- stats$multipop_qc_fail + 1
            return(NULL)
        }

        if (!nrow(g1) || !nrow(g2)) {
            return(skip_gene("Insufficient GWAS data"))
        }

        info1 <- fread(paste0(prefix1, ".info"))
        common1 <- intersect(g1$snp, info1$ID)
        if (!length(common1)) {
            return(skip_gene("POP1: No common SNPs"))
        }

        LD1_raw <- ReadLDMbin(prefix1)
        if (is.null(LD1_raw)) {
            return(skip_gene("POP1: Failed to read LD matrix"))
        }

        rownames(LD1_raw) <- colnames(LD1_raw) <- info1$ID
        g1_qc <- qc_gwas_op(g1, data.frame(snp = common1, stringsAsFactors = FALSE))
        R1_qc <- qc_LD_op(g1_qc, LD1_raw, data.frame(snp = common1, stringsAsFactors = FALSE))

        info2 <- fread(paste0(prefix2, ".info"))
        common2 <- intersect(g2$snp, info2$ID)
        if (!length(common2)) {
            return(skip_gene("POP2: No common SNPs"))
        }

        LD2_raw <- ReadLDMbin(prefix2)
        if (is.null(LD2_raw)) {
            return(skip_gene("POP2: Failed to read LD matrix"))
        }

        rownames(LD2_raw) <- colnames(LD2_raw) <- info2$ID
        g2_qc <- qc_gwas_op(g2, data.frame(snp = common2, stringsAsFactors = FALSE))
        R2_qc <- qc_LD_op(g2_qc, LD2_raw, data.frame(snp = common2, stringsAsFactors = FALSE))

        common_snps <- intersect(g1_qc$snp, g2_qc$snp)
        if (!length(common_snps)) {
            return(skip_gene("No common SNPs between populations"))
        }

        pos_dt <- snp_gene[gene == g & snp %in% common_snps, .(snp, pos)]
        setorder(pos_dt, pos)
        common_snps <- pos_dt$snp

        df1_b <- allele_flip_gwas_beta(
            g1_qc[match(common_snps, snp)],
            info1[match(common_snps, ID), ],
            no_AF = FALSE,
            af_threshold = af_threshold
        )
        if (is.null(df1_b)) {
            return(skip_gene("POP1: No SNPs retained after allele alignment or AF filtering"))
        }

        df2_b <- allele_flip_gwas_beta(
            g2_qc[match(common_snps, snp)],
            info2[match(common_snps, ID), ],
            no_AF = FALSE,
            af_threshold = af_threshold
        )
        if (is.null(df2_b)) {
            return(skip_gene("POP2: No SNPs retained after allele alignment or AF filtering"))
        }

        common_snps2 <- intersect(df1_b$snp, df2_b$snp)
        if (!length(common_snps2)) {
            return(skip_gene("No common SNPs retained after allele alignment or AF filtering"))
        }

        df1_b <- df1_b[match(common_snps2, snp)]
        df2_b <- df2_b[match(common_snps2, snp)]

        R1_fin <- R1_qc[common_snps2, common_snps2, drop = FALSE]
        R2_fin <- R2_qc[common_snps2, common_snps2, drop = FALSE]

        z1 <- df1_b$beta / df1_b$se
        z2 <- df2_b$beta / df2_b$se

        m <- length(common_snps2)
        Rj <- rbind(
            cbind(R1_fin, matrix(0, m, m)),
            cbind(matrix(0, m, m), R2_fin)
        )

        gwas_fin <- rbind(
            data.table(snp = df1_b$snp, P = 2 * pnorm(-abs(z1))),
            data.table(snp = df2_b$snp, P = 2 * pnorm(-abs(z2)))
        )

        eig <- eigen(Rj, symmetric = TRUE)

        chr_g <- unique(snp_gene[gene == g, chr])
        this_annot <- gene_annot_raw[gene.id == g & chr == chr_g]
        setnames(this_annot, "gene.id", "gene")

        res <- tryCatch({
            mbat_SVD(
                gwas_noqc_final = gwas_fin,
                eig = eig,
                prop = prop,
                z = c(z1, z2),
                chr_name = chr_g,
                ID = g,
                gene_annotate = this_annot,
                top1 = g1_qc[which.min(P), snp],
                top1_p = min(g1_qc$P),
                top2 = g2_qc[which.min(P), snp],
                top2_p = min(g2_qc$P)
            )
        }, error = function(e) {
            log_warn(sprintf("Gene %s: mbat_SVD error: %s", g, e$message))
            return(NULL)
        })

        if (!is.null(res)) {
            stats$multipop_success <<- stats$multipop_success + 1
            res[1, ]
        } else {
            skip_gene("mbat_SVD processing failed")
        }
    }

    ## 7. Multipop execution
    log_info("Executing Multipop analysis...")
    results <- lapply(valid_genes, function(g) tryCatch(process_gene(g), error = function(e) {
        log_error(sprintf("Gene %s (Joint) error: %s", g, e$message))
        return(NULL)
    }))
    out_df <- rbindlist(Filter(Negate(is.null), results))

    ## 8. Single-pop execution
    log_info("Executing Single-pop analysis...")
    genes_pop1_only <- setdiff(names(m1_list), names(m2_list))
    genes_pop2_only <- setdiff(names(m2_list), names(m1_list))

    gene_prefix1_single <- setNames(sapply(genes_pop1_only, function(g) {
        f <- bins1_all[grep(paste0("^", g, "\\."), basename(bins1_all), perl = TRUE)]
        if (length(f)) sub("\\.bin$", "", f[1]) else NA_character_
    }, USE.NAMES = TRUE), genes_pop1_only)

    gene_prefix2_single <- setNames(sapply(genes_pop2_only, function(g) {
        f <- bins2_all[grep(paste0("^", g, "\\."), basename(bins2_all), perl = TRUE)]
        if (length(f)) sub("\\.bin$", "", f[1]) else NA_character_
    }, USE.NAMES = TRUE), genes_pop2_only)

    stats$singlepop_ld_missing <- sum(is.na(gene_prefix1_single)) + sum(is.na(gene_prefix2_single))

    process_gene_single <- function(g, pop, g_list, prefix_map) {
        prefix <- prefix_map[[g]]

        skip_gene_single <- function(reason) {
            log_warn(sprintf("Skipping gene %s (%s): %s", g, pop, reason))
            stats$singlepop_qc_fail <<- stats$singlepop_qc_fail + 1
            return(NULL)
        }

        if (is.na(prefix)) {
            return(NULL)
        }

        gdat <- g_list[[g]]
        if (!nrow(gdat)) {
            return(skip_gene_single("No GWAS data"))
        }

        info <- fread(paste0(prefix, ".info"))
        common <- intersect(gdat$snp, info$ID)
        if (!length(common)) {
            return(skip_gene_single("No common SNPs"))
        }

        LD_raw <- ReadLDMbin(prefix)
        if (is.null(LD_raw)) {
            return(skip_gene_single("Failed to read LD matrix"))
        }

        rownames(LD_raw) <- colnames(LD_raw) <- info$ID
        g_qc <- qc_gwas_op(gdat, data.frame(snp = common, stringsAsFactors = FALSE))
        R_qc <- qc_LD_op(g_qc, LD_raw, data.frame(snp = common, stringsAsFactors = FALSE))

        if (!nrow(g_qc) || !ncol(R_qc)) {
            return(skip_gene_single("QC failed"))
        }

        pos_dt <- snp_gene[gene == g & snp %in% g_qc$snp, .(snp, pos)]
        setorder(pos_dt, pos)
        ordered_snps <- pos_dt$snp

        g_qc <- g_qc[match(ordered_snps, snp)]
        R_qc <- R_qc[ordered_snps, ordered_snps, drop = FALSE]

        df_b <- allele_flip_gwas_beta(
            g_qc,
            info[match(ordered_snps, ID), ],
            no_AF = FALSE,
            af_threshold = af_threshold
        )
        if (is.null(df_b)) {
            return(skip_gene_single("No SNPs retained after allele alignment or AF filtering"))
        }

        kept_snps <- df_b$snp
        R_qc <- R_qc[kept_snps, kept_snps, drop = FALSE]

        z <- df_b$beta / df_b$se
        eig <- eigen(R_qc, symmetric = TRUE)

        posval <- eig$values[eig$values > 1e-10]

        if (!length(posval)) {
            return(skip_gene_single("No positive LD eigenvalues"))
        }

        if (prop == 1) {
            k <- length(posval)
        } else {
            cum_var <- cumsum(posval / sum(posval))
            k <- which(cum_var > as.numeric(prop))[1]
        }

        sel <- seq_len(k)
        lambda_prop <- eig$values[sel]
        U_prop <- eig$vectors[, sel, drop = FALSE]

        UPz_prop <- crossprod(U_prop, z)
        chisq_mBAT <- as.numeric(crossprod(UPz_prop, (1 / lambda_prop) * UPz_prop))
        P_mBATcombo <- pchisq(chisq_mBAT, df = k, lower.tail = FALSE)

        chr_g <- unique(snp_gene[gene == g, chr])
        this_annot <- gene_annot_raw[gene.id == g & chr == chr_g]

        stats$singlepop_success <<- stats$singlepop_success + 1

        data.table(
            Gene = g,
            POP = pop,
            Chr = chr_g,
            Start = this_annot$start,
            End = this_annot$end,
            No.SNPs = length(kept_snps),
            SNP_start = min(pos_dt$pos),
            SNP_end = max(pos_dt$pos),
            TopSNP = df_b$snp[which.min(df_b$P)],
            TopSNP_Pvalue = min(df_b$P),
            No.Eigenvalues = k,
            Chisq_mBAT = chisq_mBAT,
            P_mBATcombo = P_mBATcombo
        )
    }

    single_results <- list()
    if (length(genes_pop1_only)) {
        single_results <- c(
            single_results,
            lapply(genes_pop1_only, function(g) process_gene_single(g, "POP1", m1_list, gene_prefix1_single))
        )
    }
    if (length(genes_pop2_only)) {
        single_results <- c(
            single_results,
            lapply(genes_pop2_only, function(g) process_gene_single(g, "POP2", m2_list, gene_prefix2_single))
        )
    }
    single_out <- rbindlist(Filter(Negate(is.null), single_results))

    ## 9. Run GCTA Analysis
    if (run_gcta) {
        if (!is.null(plink_bfile_pop1) && !is.null(plink_bfile_pop2)) {
            gcta_wind_kb <- extend / 1000
            run_gcta_mbat(
                gwas1,
                gene_annot_raw,
                plink_bfile_pop1,
                paste0(result_prefix, ".gcta.pop1"),
                gcta_bin,
                af_threshold,
                gcta_wind_kb
            )
            run_gcta_mbat(
                gwas2,
                gene_annot_raw,
                plink_bfile_pop2,
                paste0(result_prefix, ".gcta.pop2"),
                gcta_bin,
                af_threshold,
                gcta_wind_kb
            )
        }
    }

    ## 10. Run multiBAT-combo Analysis (Weighted ACATO)
    log_info("Starting final multiBAT-combo P-value merging...")
    out_df_final <- copy(out_df)

    if ("P_mBAT" %in% names(out_df_final)) {
        setnames(out_df_final, "P_mBAT", "P_multiBAT")
    }

    if (run_gcta && nrow(out_df_final) > 0) {
        f1 <- paste0(result_prefix, ".gcta.pop1.gene.assoc.mbat")
        f2 <- paste0(result_prefix, ".gcta.pop2.gene.assoc.mbat")

        if (file.exists(f1) && file.exists(f2)) {
            dt_gcta1 <- fread(f1, select = c("Gene", "P_mBATcombo"))
            setnames(dt_gcta1, "P_mBATcombo", "P_pop1")

            dt_gcta2 <- fread(f2, select = c("Gene", "P_mBATcombo"))
            setnames(dt_gcta2, "P_mBATcombo", "P_pop2")

            dt_acato <- merge(
                out_df_final[, .(Gene, P_joint = P_multiBAT)],
                dt_gcta1,
                by = "Gene",
                all.x = TRUE
            )
            dt_acato <- merge(dt_acato, dt_gcta2, by = "Gene", all.x = TRUE)

            if (exists("ACATO")) {
                out_df_final[, P_multiBAT.combo := apply(
                    dt_acato[, .(P_joint, P_pop1, P_pop2)],
                    1,
                    function(p) ACATO(as.numeric(p), w = pop_weights)
                )]
            } else {
                out_df_final[, P_multiBAT.combo := P_multiBAT]
            }
        } else {
            out_df_final[, P_multiBAT.combo := P_multiBAT]
        }
    } else {
        if (nrow(out_df_final) > 0) {
            out_df_final[, P_multiBAT.combo := P_multiBAT]
        }
    }

    ## 11. Write outputs
    if (nrow(out_df_final) > 0) {
        cols <- names(out_df_final)
        main_cols <- setdiff(cols, c("No.Eigenvalues", "Chisq_mBAT", "P_multiBAT", "P_multiBAT.combo"))
        setcolorder(out_df_final, c(main_cols, "No.Eigenvalues", "Chisq_mBAT", "P_multiBAT", "P_multiBAT.combo"))
        fwrite(out_df_final, file = paste0(result_prefix, ".multipop.txt"), sep = "\t", quote = FALSE)
    }

    if (nrow(single_out) > 0) {
        fwrite(single_out, file = paste0(result_prefix, ".singlepop.txt"), sep = "\t", quote = FALSE)
    }

    log_info("Analysis complete.", console_print = TRUE)

    invisible(list(multipop = out_df_final, singlepop = single_out))
}