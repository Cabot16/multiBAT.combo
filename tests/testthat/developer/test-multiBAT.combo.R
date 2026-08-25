library(testthat)


# ============================================================
# Helpers
# ============================================================


get_test_inputs <- function() {

    data_dir <- system.file(
        "extdata",
        package = "multiBAT.combo"
    )


    list(
        data_dir = data_dir,

        bim_file = file.path(
            data_dir,
            "1000G_eur_unrel_ukbcom22.bim"
        ),

        bfile_prefix = file.path(
            data_dir,
            "1000G_eur_unrel_ukbcom22"
        ),

        map_file = file.path(
            data_dir,
            "hg19_v40_glist_Pat_v1.txt"
        ),

        assoc1 = file.path(
            data_dir,
            "fallback_pop1.ma"
        ),

        assoc2 = file.path(
            data_dir,
            "fallback_pop2.ma"
        ),

        LD_path1 = file.path(
            data_dir,
            "test_LD_Chr22",
            "pop1"
        ),

        LD_path2 = file.path(
            data_dir,
            "test_LD_Chr22",
            "pop2"
        )
    )
}



find_gcta <- function() {

    configured_gcta <- Sys.getenv(
        "MULTIBAT_GCTA_BIN",
        unset = ""
    )


    if (nzchar(configured_gcta)) {

        if (file.exists(configured_gcta)) {

            return(
                normalizePath(
                    configured_gcta,
                    mustWork = TRUE
                )
            )
        }


        configured_exec <- Sys.which(
            configured_gcta
        )


        if (nzchar(configured_exec)) {

            return(
                unname(configured_exec)
            )
        }
    }


    gcta_exec <- Sys.which(
        "gcta64"
    )

    if (nzchar(gcta_exec)) {
        return(
            unname(gcta_exec)
        )
    }


    gcta_exec <- Sys.which(
        "gcta"
    )

    if (nzchar(gcta_exec)) {
        return(
            unname(gcta_exec)
        )
    }


    ""
}



check_test_inputs <- function(
    inputs,
    require_gcta_reference = FALSE
) {


    expect_true(
        dir.exists(inputs$data_dir),
        info = "Package extdata directory was not found."
    )


    required_files <- c(
        inputs$bim_file,
        inputs$map_file,
        inputs$assoc1,
        inputs$assoc2
    )


    expect_true(
        all(file.exists(required_files)),
        info = paste(
            "One or more required multiBAT.combo test files were missing:",
            paste(
                required_files[
                    !file.exists(required_files)
                ],
                collapse = ", "
            )
        )
    )


    expect_true(
        dir.exists(inputs$LD_path1),
        info = "POP1 LD test directory was not found."
    )


    expect_true(
        dir.exists(inputs$LD_path2),
        info = "POP2 LD test directory was not found."
    )


    if (require_gcta_reference) {


        required_plink_files <- paste0(
            inputs$bfile_prefix,
            c(
                ".bed",
                ".bim",
                ".fam"
            )
        )


        expect_true(
            all(file.exists(required_plink_files)),
            info = paste(
                "Required PLINK reference files were missing:",
                paste(
                    required_plink_files[
                        !file.exists(required_plink_files)
                    ],
                    collapse = ", "
                )
            )
        )
    }


    invisible(TRUE)
}



expect_required_columns <- function(
    x,
    expected_columns,
    output_name = "multipop output"
) {


    missing_columns <- setdiff(
        expected_columns,
        names(x)
    )


    expect_true(
        length(missing_columns) == 0L,
        info = paste0(
            "Missing columns in ",
            output_name,
            ": ",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    )
}



expect_absent_columns <- function(
    x,
    absent_columns,
    output_name = "multipop output"
) {


    unexpected_columns <- intersect(
        absent_columns,
        names(x)
    )


    expect_true(
        length(unexpected_columns) == 0L,
        info = paste0(
            "Unexpected columns in ",
            output_name,
            ": ",
            paste(
                unexpected_columns,
                collapse = ", "
            )
        )
    )
}
# ============================================================
# Common result checks
# ============================================================


expect_common_result <- function(
    res,
    result_prefix
) {


    expect_true(
        is.list(res),
        info = "multiBAT.combo did not return a list."
    )


    expect_true(
        "multipop" %in% names(res),
        info = "Returned object does not contain 'multipop'."
    )


    expect_true(
        data.table::is.data.table(
            res$multipop
        ),
        info = "res$multipop is not a data.table."
    )


    expect_true(
        nrow(res$multipop) > 0L,
        info = "No multiBAT.combo results were produced."
    )


    required_columns <- c(

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

        "No.Eigenvalues",

        "P_multiBAT.combo",
        "P_multiBAT.combo_BH"
    )


    expect_required_columns(
        res$multipop,
        required_columns
    )


    # --------------------------------------------------------
    # Fixed developer test data:
    #
    # CCT8L2:
    # two-population analysis
    #
    # BID:
    # POP1-only fallback
    # --------------------------------------------------------


    expect_equal(
        nrow(res$multipop),
        2L,
        info = paste0(
            "Developer test data should return exactly two genes."
        )
    )


    fallback_rows <- (
        !is.na(
            res$multipop$No.SNPs_pop1
        ) &
        is.na(
            res$multipop$No.SNPs_pop2
        )
    )


    regular_rows <- (
        !is.na(
            res$multipop$No.SNPs_pop1
        ) &
        !is.na(
            res$multipop$No.SNPs_pop2
        )
    )


    expect_equal(
        sum(fallback_rows),
        1L,
        info = "Expected exactly one POP1-only fallback gene."
    )


    expect_equal(
        sum(regular_rows),
        1L,
        info = "Expected exactly one two-population gene."
    )


    expect_true(
        all(
            is.finite(
                res$multipop[
                    ,
                    P_multiBAT.combo
                ]
            )
        ),
        info = "Invalid P_multiBAT.combo values were produced."
    )


    expect_true(
        any(
            is.finite(
                res$multipop$P_multiBAT.combo_BH
            )
        ),
        info = "No valid BH-adjusted P values were produced."
    )


    # --------------------------------------------------------
    # Output file check
    # --------------------------------------------------------


    output_file <- paste0(
        result_prefix,
        ".multipop.txt"
    )


    cat("\n============================\n")
    cat("Checking multipop output file:\n")
    cat(output_file, "\n")

    cat("\nFiles in output directory:\n")
    print(
        list.files(
            dirname(output_file),
            full.names = TRUE
        )
    )

    cat("============================\n")


    expect_true(
        file.exists(output_file),
        info = paste0(
            "The multipop output file was not created: ",
            output_file
        )
    )


    expect_true(
        file.info(output_file)$size > 0,
        info = paste0(
            "The multipop output file is empty: ",
            output_file
        )
    )


    written <- data.table::fread(
        output_file,
        showProgress = FALSE
    )


    expect_equal(
        nrow(written),
        2L,
        info = "Written output does not contain two genes."
    )


    expect_true(
    	setequal(
        written$Gene,
        res$multipop$Gene
        ),
        info = "Written genes do not match returned genes."
    )

}



# ============================================================
# Internal mBAT mode
# ============================================================


run_internal_mode <- function(
    print_all_p
) {


    inputs <- get_test_inputs()


    check_test_inputs(
        inputs,
        require_gcta_reference = FALSE
    )


    out_dir <- tempfile(
        "multiBAT_internal_"
    )


    dir.create(
        out_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )



    result_prefix <- file.path(
        out_dir,
        ifelse(
            print_all_p,
            "internal_all_p",
            "internal_default"
        )
    )


    res <- multiBAT.combo::multiBAT.combo(
        
        bim_file = inputs$bim_file,

        map_file = inputs$map_file,

        assoc_pop1 = inputs$assoc1,

        assoc_pop2 = inputs$assoc2,


        LD_path_pop1 = inputs$LD_path1,

        LD_path_pop2 = inputs$LD_path2,


        result_prefix = result_prefix,


        extend = 50000,

        prop = 0.9,

        af_threshold = 1,


        run_gcta = FALSE,

        print_all_p = print_all_p
    )


    list(
        res = res,
        result_prefix = result_prefix
    )
}



# ============================================================
# GCTA mode
# ============================================================


run_gcta_mode <- function(
    print_all_p
) {


    gcta_exec <- find_gcta()


    expect_true(
        nzchar(gcta_exec),
        info = paste0(
            "GCTA executable was not detected. ",
            "Please set MULTIBAT_GCTA_BIN."
        )
    )


    inputs <- get_test_inputs()


    check_test_inputs(
        inputs,
        require_gcta_reference = TRUE
    )


    out_dir <- tempfile(
        "multiBAT_gcta_"
    )


    dir.create(
        out_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )



    result_prefix <- file.path(
        out_dir,
        ifelse(
            print_all_p,
            "gcta_all_p",
            "gcta_default"
        )
    )


    res <- multiBAT.combo::multiBAT.combo(

        bim_file = inputs$bim_file,

        map_file = inputs$map_file,

        assoc_pop1 = inputs$assoc1,

        assoc_pop2 = inputs$assoc2,


        LD_path_pop1 = inputs$LD_path1,

        LD_path_pop2 = inputs$LD_path2,


        result_prefix = result_prefix,


        extend = 50000,

        prop = 0.9,

        af_threshold = 1,


        run_gcta = TRUE,

        print_all_p = print_all_p,


        gcta_bin = gcta_exec,


        plink_bfile_pop1 = inputs$bfile_prefix,

        plink_bfile_pop2 = inputs$bfile_prefix
    )


    list(
        res = res,
        result_prefix = result_prefix
    )
}
# ============================================================
# Mode 1
#
# run_gcta = FALSE
# print_all_p = FALSE
# ============================================================

test_that(
    "multiBAT.combo runs with internal mBAT and default P-value output",
    {

        test_result <- run_internal_mode(
            print_all_p = FALSE
        )


        res <- test_result$res


        expect_common_result(
            res,
            test_result$result_prefix
        )


        expect_absent_columns(
            res$multipop,
            c(
                "T_multiBAT",
                "P_multiBAT",

                "T_mBAT_pop1",
                "P_mBAT_pop1",
                "No.Eigenvalues_pop1",

                "T_mBAT_pop2",
                "P_mBAT_pop2",
                "No.Eigenvalues_pop2"
            )
        )
    }
)



# ============================================================
# Mode 2
#
# run_gcta = FALSE
# print_all_p = TRUE
# ============================================================

test_that(
    "multiBAT.combo runs with internal mBAT and all P values",
    {

        test_result <- run_internal_mode(
            print_all_p = TRUE
        )


        res <- test_result$res


        expect_common_result(
            res,
            test_result$result_prefix
        )


        expect_required_columns(
            res$multipop,
            c(
                "T_multiBAT",
                "P_multiBAT",

                "T_mBAT_pop1",
                "P_mBAT_pop1",
                "No.Eigenvalues_pop1",

                "T_mBAT_pop2",
                "P_mBAT_pop2",
                "No.Eigenvalues_pop2"
            )
        )


        expect_required_columns(
            res$multipop,
            c(
                "No.Eigenvalues_pop1",
                "No.Eigenvalues_pop2"
            )
        )


        expect_true(
            any(
                is.finite(
                    res$multipop$P_multiBAT
                )
            ),
            info = "No valid P_multiBAT values were produced."
        )


        expect_true(
            any(
                is.finite(
                    res$multipop$P_mBAT_pop1
                )
            ),
            info = "No valid P_mBAT_pop1 values were produced."
        )


        expect_true(
            any(
                is.finite(
                    res$multipop$P_mBAT_pop2
                )
            ),
            info = "No valid P_mBAT_pop2 values were produced."
        )


        expect_true(
            any(
                is.finite(
                    res$multipop$No.Eigenvalues_pop1
                )
            ),
            info = "No valid No.Eigenvalues_pop1 values were produced."
        )


        expect_true(
            any(
                is.finite(
                    res$multipop$No.Eigenvalues_pop2
                )
            ),
            info = "No valid No.Eigenvalues_pop2 values were produced."
        )

    }
)



# ============================================================
# Mode 3
#
# run_gcta = TRUE
# print_all_p = FALSE
# ============================================================

test_that(
    "multiBAT.combo runs with GCTA and default P-value output",
    {

        test_result <- run_gcta_mode(
            print_all_p = FALSE
        )


        res <- test_result$res


        expect_common_result(
            res,
            test_result$result_prefix
        )


        # ----------------------------------------------------
        # GCTA output is always reported when run_gcta = TRUE
        # ----------------------------------------------------

        expect_required_columns(
            res$multipop,
            c(

                "P_mBAT_pop1",
                "P_fastBAT_pop1",
                "P_mBATcombo_pop1",

                "P_mBAT_pop2",
                "P_fastBAT_pop2",
                "P_mBATcombo_pop2"

            )
        )


        # Joint statistics are controlled by print_all_p

        expect_absent_columns(
            res$multipop,
            c(
                "T_multiBAT",
                "P_multiBAT"
            )
        )


        expect_true(
            any(
                is.finite(
                    res$multipop$P_mBAT_pop1
                )
            ),
            info = "No valid GCTA P_mBAT_pop1 values were produced."
        )


        expect_true(
            any(
                is.finite(
                    res$multipop$P_mBAT_pop2
                )
            ),
            info = "No valid GCTA P_mBAT_pop2 values were produced."
        )


        # ----------------------------------------------------
        # GCTA files
        # ----------------------------------------------------

        pop1_file <- paste0(
            test_result$result_prefix,
            ".gcta.pop1.gene.assoc.mbat"
        )


        pop2_file <- paste0(
            test_result$result_prefix,
            ".gcta.pop2.gene.assoc.mbat"
        )


        expect_true(
            file.exists(pop1_file),
            info = "POP1 GCTA output file was not created."
        )


        expect_true(
            file.exists(pop2_file),
            info = "POP2 GCTA output file was not created."
        )


        expect_true(
            file.info(pop1_file)$size > 0,
            info = "POP1 GCTA output file is empty."
        )


        expect_true(
            file.info(pop2_file)$size > 0,
            info = "POP2 GCTA output file is empty."
        )
        pop1_gcta <- data.table::fread(
            pop1_file,
            showProgress = FALSE
        )


        pop2_gcta <- data.table::fread(
            pop2_file,
            showProgress = FALSE
        )


        expect_true(
            "No.Eigenvalues" %in% names(pop1_gcta),
            info = "POP1 GCTA output does not contain No.Eigenvalues."
        )


        expect_true(
            "No.Eigenvalues" %in% names(pop2_gcta),
            info = "POP2 GCTA output does not contain No.Eigenvalues."
        )


    }
)



# ============================================================
# Mode 4
#
# run_gcta = TRUE
# print_all_p = TRUE
# ============================================================

test_that(
    "multiBAT.combo runs with GCTA and all P values",
    {

        test_result <- run_gcta_mode(
            print_all_p = TRUE
        )


        res <- test_result$res


        expect_common_result(
            res,
            test_result$result_prefix
        )


        expect_required_columns(
            res$multipop,
            c(

                # joint multiBAT

                "T_multiBAT",
                "P_multiBAT",


                # POP1 GCTA outputs

                "P_mBAT_pop1",
                "P_fastBAT_pop1",
                "P_mBATcombo_pop1",
                "No.Eigenvalues_pop1",


                # POP2 GCTA outputs

                "P_mBAT_pop2",
                "P_fastBAT_pop2",
                "P_mBATcombo_pop2",
                "No.Eigenvalues_pop2"

            )
        )


        expect_true(
            any(
                is.finite(
                    res$multipop$P_multiBAT
                )
            ),
            info = "No valid P_multiBAT values were produced."
        )


        expect_true(
            any(
                is.finite(
                    res$multipop$P_mBAT_pop1
                )
            ),
            info = "No valid GCTA P_mBAT_pop1 values were produced."
        )


        expect_true(
            any(
                is.finite(
                    res$multipop$P_mBAT_pop2
                )
            ),
            info = "No valid GCTA P_mBAT_pop2 values were produced."
        )


        pop1_file <- paste0(
            test_result$result_prefix,
            ".gcta.pop1.gene.assoc.mbat"
        )


        pop2_file <- paste0(
            test_result$result_prefix,
            ".gcta.pop2.gene.assoc.mbat"
        )


        expect_true(
            file.exists(pop1_file),
            info = "POP1 GCTA output file was not created."
        )


        expect_true(
            file.exists(pop2_file),
            info = "POP2 GCTA output file was not created."
        )


        expect_true(
            file.info(pop1_file)$size > 0,
            info = "POP1 GCTA output file is empty."
        )


        expect_true(
            file.info(pop2_file)$size > 0,
            info = "POP2 GCTA output file is empty."
        )
        pop1_gcta <- data.table::fread(
            pop1_file,
            showProgress = FALSE
        )


        pop2_gcta <- data.table::fread(
            pop2_file,
            showProgress = FALSE
        )


        expect_true(
            "No.Eigenvalues" %in% names(pop1_gcta),
            info = "POP1 GCTA output does not contain No.Eigenvalues."
        )


        expect_true(
            "No.Eigenvalues" %in% names(pop2_gcta),
            info = "POP2 GCTA output does not contain No.Eigenvalues."
        )


    }
)
