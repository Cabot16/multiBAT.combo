library(testthat)


# ============================================================
# Lightweight user test for multiBAT.combo
#
# This test only verifies that the main workflow can run using
# bundled example data.
#
# No external dependency:
#   - GCTA
#   - GCTB
#   - PLINK
#
# ============================================================


test_that(
    "multiBAT.combo runs with bundled example data",
    {


        # ----------------------------------------------------
        # Locate package example data
        # ----------------------------------------------------

        data_dir <- system.file(
            "extdata",
            package = "multiBAT.combo"
        )


        expect_true(
            nzchar(data_dir),
            info = "Package extdata directory was not found."
        )


        # ----------------------------------------------------
        # Input files
        # ----------------------------------------------------

        bim_file <- file.path(
            data_dir,
            "1000G_eur_unrel_ukbcom22.bim"
        )


        map_file <- file.path(
            data_dir,
            "hg19_v40_glist_Pat_v1.txt"
        )


        assoc_pop1 <- file.path(
            data_dir,
            "fallback_pop1.ma"
        )


        assoc_pop2 <- file.path(
            data_dir,
            "fallback_pop2.ma"
        )


        LD_path_pop1 <- file.path(
            data_dir,
            "test_LD_Chr22",
            "pop1"
        )


        LD_path_pop2 <- file.path(
            data_dir,
            "test_LD_Chr22",
            "pop2"
        )


        required_files <- c(
            bim_file,
            map_file,
            assoc_pop1,
            assoc_pop2
        )


        expect_true(
            all(file.exists(required_files)),
            info = paste(
                "Missing bundled test files:",
                paste(
                    required_files[
                        !file.exists(required_files)
                    ],
                    collapse = ", "
                )
            )
        )


        expect_true(
            dir.exists(LD_path_pop1),
            info = "POP1 LD directory missing."
        )


        expect_true(
            dir.exists(LD_path_pop2),
            info = "POP2 LD directory missing."
        )



        # ----------------------------------------------------
        # Run lightweight workflow
        # ----------------------------------------------------

        out_prefix <- file.path(
            tempdir(),
            "multiBAT_combo_user_test"
        )


        result <- multiBAT.combo::multiBAT.combo(

            bim_file = bim_file,

            map_file = map_file,

            assoc_pop1 = assoc_pop1,

            assoc_pop2 = assoc_pop2,

            LD_path_pop1 = LD_path_pop1,

            LD_path_pop2 = LD_path_pop2,

            result_prefix = out_prefix,

            extend = 50000,

            prop = 0.9,

            af_threshold = 1,

            run_gcta = FALSE,

            print_all_p = FALSE
        )



        # ----------------------------------------------------
        # Check returned object
        # ----------------------------------------------------

        expect_true(
            is.list(result),
            info = "multiBAT.combo did not return a list."
        )


        expect_true(
            "multipop" %in% names(result),
            info = "Returned object does not contain multipop."
        )


        expect_true(
            data.table::is.data.table(
                result$multipop
            ),
            info = "multipop output is not a data.table."
        )


        expect_gt(
            nrow(result$multipop),
            0
        )



        # ----------------------------------------------------
        # Check essential output columns
        # ----------------------------------------------------

        expected_columns <- c(
            "Gene",
            "Chr",
            "Start",
            "End",
            "No.SNPs_pop1",
            "No.SNPs_pop2",
            "P_multiBAT.combo",
            "P_multiBAT.combo_BH"
        )


        expect_true(
            all(
                expected_columns %in%
                    names(result$multipop)
            ),
            info = paste(
                "Missing output columns:",
                paste(
                    setdiff(
                        expected_columns,
                        names(result$multipop)
                    ),
                    collapse = ", "
                )
            )
        )



        # ----------------------------------------------------
        # Check written result file
        # ----------------------------------------------------

        output_file <- paste0(
            out_prefix,
            ".multipop.txt"
        )


        expect_true(
            file.exists(output_file),
            info = "Output multipop file was not created."
        )


        expect_gt(
            file.info(output_file)$size,
            0
        )


        written <- data.table::fread(
            output_file,
            showProgress = FALSE
        )


        expect_gt(
            nrow(written),
            0
        )


        expect_true(
            "Gene" %in% names(written)
        )

    }
)
