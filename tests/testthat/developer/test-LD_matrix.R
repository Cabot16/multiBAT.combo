library(testthat)


test_that(
    "LD_matrix runs and creates valid LD files",
    {

        # ========================================================
        # 1. Testing environment
        # ========================================================

        withr::local_envvar(
            R_TESTS = ""
        )


        # ========================================================
        # 2. Locate package test data
        # ========================================================

        data_dir <- system.file(
            "extdata",
            package = "multiBAT.combo"
        )


        skip_if_not(
            dir.exists(data_dir),
            message = "Package extdata directory was not found."
        )


        # ========================================================
        # 3. Check PLINK reference files
        #
        # POP1 and POP2 intentionally use the same PLINK panel
        # to minimise package size.
        # ========================================================

        plink_base <- "1000G_eur_unrel_ukbcom22"


        plink_path_prefix <- file.path(
            data_dir,
            plink_base
        )


        required_plink_files <- paste0(
            plink_path_prefix,
            c(
                ".bed",
                ".bim",
                ".fam"
            )
        )


        if (
            !all(
                file.exists(required_plink_files)
            )
        ) {

            skip(
                "Required PLINK reference files were unavailable."
            )
        }


        # ========================================================
        # 4. Check external software
        # ========================================================

        sys_gctb <- Sys.which(
            "gctb"
        )

        sys_plink <- Sys.which(
            "plink"
        )


        if (
            !nzchar(sys_gctb) ||
            !nzchar(sys_plink)
        ) {

            skip(
                "GCTB or PLINK executable was unavailable."
            )
        }


        # ========================================================
        # 5. Prepare test gene annotation
        # ========================================================

        map_file <- file.path(
            data_dir,
            "hg19_v40_glist_Pat_v1.txt"
        )


        skip_if_not(
            file.exists(map_file),
            message = "Gene annotation file was not found."
        )


        gene_map <- data.table::fread(
            map_file,
            header = FALSE,
            sep = "\t",
            showProgress = FALSE
        )


        skip_if(
            nrow(gene_map) == 0L,
            "Gene annotation file is empty."
        )


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


        genes_to_check <- head(
            unique(
                as.character(
                    gene_map$Gene
                )
            ),
            2L
        )


        selected_rows <- match(
            genes_to_check,
            gene_map$Gene
        )


        test_gene_map <- gene_map[
            selected_rows,
            .(
                Chr,
                Start,
                End,
                Gene
            )
        ]


        test_map_file <- file.path(
            withr::local_tempdir(),
            "test_gene_map.txt"
        )


        data.table::fwrite(
            test_gene_map,
            test_map_file,
            sep = "\t",
            col.names = FALSE,
            quote = FALSE
        )


        # ========================================================
        # 6. Output directory
        # ========================================================

        out_root_dir <- withr::local_tempdir()


        pop1_ld_dir <- file.path(
            out_root_dir,
            "pop1"
        )


        pop2_ld_dir <- file.path(
            out_root_dir,
            "pop2"
        )


        message(
            "\nStarting LD_matrix test for ",
            length(genes_to_check),
            " genes across two populations..."
        )


        # ========================================================
        # 7. Run LD_matrix
        # ========================================================

        ldm_results <- multiBAT.combo::LD_matrix(
            gctb_bin = sys_gctb,
            map_file = test_map_file,
            plink_bfile_pop1 = plink_path_prefix,
            plink_bfile_pop2 = plink_path_prefix,
            extend = 50000,
            result_prefix = out_root_dir,
            thread = 1
        )


        status <- unlist(
            ldm_results,
            use.names = FALSE
        )


        expect_true(
            length(status) > 0L,
            info = "LD_matrix returned no status information."
        )


        expect_false(
            any(
                grepl(
                    "^Failed",
                    status
                )
            ),
            info = paste(
                "LD_matrix failure:",
                paste(status, collapse = "\n")
            )
        )


        expect_true(
            any(
                grepl(
                    "^Success",
                    status
                )
            ),
            info = "No LD_matrix task succeeded."
        )


        # ========================================================
        # 8. Check generated LD files
        # ========================================================

        check_population_ld_files <- function(
            ld_dir,
            pop_label
        ) {


            expect_true(
                dir.exists(ld_dir),
                info = paste(
                    pop_label,
                    "LD directory missing."
                )
            )


            bin_files <- list.files(
                ld_dir,
                pattern = "\\.ldm\\.full\\.bin$",
                full.names = TRUE
            )


            info_files <- list.files(
                ld_dir,
                pattern = "\\.ldm\\.full\\.info$",
                full.names = TRUE
            )


            expect_true(
                length(bin_files) > 0L,
                info = paste(
                    pop_label,
                    "no LD binary files generated."
                )
            )


            expect_true(
                length(info_files) > 0L,
                info = paste(
                    pop_label,
                    "no LD info files generated."
                )
            )


            for (bin_file in bin_files) {


                prefix <- sub(
                    "\\.bin$",
                    "",
                    bin_file
                )


                info_file <- paste0(
                    prefix,
                    ".info"
                )


                expect_true(
                    file.exists(info_file),
                    info = paste(
                        pop_label,
                        "missing info file:",
                        basename(info_file)
                    )
                )


                expect_gt(
                    file.info(bin_file)$size,
                    0
                )


                expect_gt(
                    file.info(info_file)$size,
                    0
                )
            }


            invisible(TRUE)
        }


        # ========================================================
        # 9. Check both populations
        # ========================================================

        check_population_ld_files(
            pop1_ld_dir,
            "POP1"
        )


        check_population_ld_files(
            pop2_ld_dir,
            "POP2"
        )

    }
)