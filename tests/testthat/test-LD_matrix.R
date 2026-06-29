library(testthat)

# Test the LD_matrix utility function
test_that("LD_matrix utility runs and creates necessary LD files", {
    
    # Set R_TESTS environment variable to avoid mclapply issues during testing
    withr::local_envvar(R_TESTS = "")

    # 1. Locate inst/extdata
    data_dir <- system.file("extdata", package = "multiBAT.combo")
    skip_if_not(dir.exists(data_dir), message = "inst/extdata not found.")

    # 2. GCTB/PLINK Environment Setup & Logic
    # ---------------------------------------------------------
    # A. Check PLINK file integrity
    plink_base <- "1000G_eur_unrel_ukbcom22"
    plink_path_prefix <- file.path(data_dir, plink_base)
    
    has_plink_files <- file.exists(paste0(plink_path_prefix, ".bed")) && 
                       file.exists(paste0(plink_path_prefix, ".fam")) &&
                       file.exists(paste0(plink_path_prefix, ".bim"))

    # B. Check if gctb and plink are available on the system
    sys_gctb <- Sys.which("gctb")
    sys_plink <- Sys.which("plink")
    has_gctb_bin <- sys_gctb != ""
    has_plink_bin <- sys_plink != ""

    # C. Decide whether to enable the test
    do_ldm_test <- has_plink_files && has_gctb_bin && has_plink_bin

    if (!do_ldm_test) {
        skip("LD_matrix test skipped: GCTB binary ('gctb') or PLINK binary/files missing.")
    }
    # ---------------------------------------------------------

    # 3. Define Inputs and Output
    map_file <- file.path(data_dir, "hg19_v40_glist_Pat_v1.txt")
    out_root_dir <- withr::local_tempdir()
    
    # Define expected output directories
    pop1_ld_dir <- file.path(out_root_dir, "pop1")
    pop2_ld_dir <- file.path(out_root_dir, "pop2")

    # Load gene map to determine which genes to check
    gene_map <- data.table::fread(map_file, header = FALSE, sep = '\t')
    data.table::setnames(gene_map, 1:4, c("Chr", "Start", "End", "Gene"))
    
    # Select the first two genes to assert existence (to keep test fast)
    genes_to_check <- head(gene_map$Gene, 2)
    skip_if(length(genes_to_check) == 0, "Gene map file is empty or cannot be read.")

    message("\nStarting LD_matrix test, processing 2 populations for 2 genes...")

    # 4. Call LD_matrix (Accessing via internal namespace ::: for robust testing)
    ldm_results <- multiBAT.combo:::LD_matrix(
        gctb_bin = sys_gctb,
        map_file = map_file,
        plink_bfile_pop1 = plink_path_prefix,
        plink_bfile_pop2 = plink_path_prefix, # Use the same reference panel for testing
        extend = 50000,
        result_prefix = out_root_dir,
        thread = 1 # Use 1 thread in tests for stability
    )

    # 5. Assert LD Matrix File Generation and Content
    
    # Create a mapping from result index to gene/pop name for dynamic checking
    gene_pop_names <- c(
        paste0(genes_to_check, "_pop1"),
        paste0(genes_to_check, "_pop2")
    )
    
    # Iterate through the results and check files ONLY for successful runs
    for (i in seq_along(ldm_results)) {
        result_msg <- ldm_results[[i]]
        
        # Determine the gene and population for the current result
        gene_pop_parts <- strsplit(gene_pop_names[i], "_")[[1]]
        gene <- gene_pop_parts[1]
        pop_name <- gene_pop_parts[2] # "pop1" or "pop2"
        pop_ld_dir <- ifelse(pop_name == "pop1", pop1_ld_dir, pop2_ld_dir)
        
        if (grepl("^Success", result_msg)) {
            
            # --- File Checks (Only executed if the function reported success) ---
            pop_bin_file <- file.path(pop_ld_dir, paste0(gene, ".ldm.full.bin"))
            pop_info_file <- file.path(pop_ld_dir, paste0(gene, ".ldm.full.info"))
            
            # Binary file existence and non-zero size
            expect_true(file.exists(pop_bin_file), info = paste(pop_name, "LDM binary file missing for gene", gene))
            expect_true(file.info(pop_bin_file)$size > 0, info = paste(pop_name, "LDM binary file is empty for gene", gene))
            
            # Info file existence and non-zero size
            expect_true(file.exists(pop_info_file), info = paste(pop_name, "LDM info file missing for gene", gene))
            expect_true(file.info(pop_info_file)$size > 0, info = paste(pop_name, "LDM info file is empty for gene", gene))

        } else if (grepl("^Failed", result_msg)) {
             # If it explicitly failed (e.g., GCTB error), fail the test immediately
             fail(paste("Task explicitly failed for gene", gene, "(", pop_name, ") with message:", result_msg))
        }
        # If the result is "Skipped", we do nothing. The test passes that task as 'Skipped' is a valid outcome.
    }
    
    # Final assertion: Ensure at least one task succeeded.
    # If the manual run succeeded, we expect success here too.
    expect_true(any(grepl("Success", unlist(ldm_results))), 
                info = "Expected at least one task to succeed based on manual check, but all tasks were Skipped or Failed.")
    
    # Final assertion 2: Ensure no tasks explicitly failed (Skipped is OK, Failed is not)
    expect_false(any(grepl("Failed", unlist(ldm_results))), 
                 info = "One or more tasks explicitly failed (Status != Skipped/Success). Check logs for GCTB/PLINK error messages.")

})