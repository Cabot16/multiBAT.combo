library(testthat)

test_that("multiBAT.combo runs on example data with optional GCTA", {
  # 1. Locate inst/extdata
  data_dir <- system.file("extdata", package = "multiBAT.combo")
  expect_true(dir.exists(data_dir), info = "inst/extdata not found")

  # 2. Construct input paths
  bim_file <- file.path(data_dir, "1000G_eur_unrel_ukbcom22.bim")
  map_file <- file.path(data_dir, "hg19_v40_glist_Pat_v1.txt")
  assoc1   <- file.path(data_dir, "test.ma")
  assoc2   <- file.path(data_dir, "test2.ma")
  LD_path1 <- file.path(data_dir, "POP1_LD_Chr22")
  LD_path2 <- file.path(data_dir, "POP1_LD_Chr22")
  out_dir  <- withr::local_tempdir()

  # 3. Ensure basic files exist
  expect_true(file.exists(bim_file))
  expect_true(file.exists(map_file))
  
  # 4. GCTA Environment Setup
  plink_base <- "1000G_eur_unrel_ukbcom22"
  plink_path_prefix <- file.path(data_dir, plink_base)
  
  has_plink_files <- file.exists(paste0(plink_path_prefix, ".bed")) && 
                     file.exists(paste0(plink_path_prefix, ".fam"))
  
  sys_gcta <- Sys.which("gcta64")
  has_gcta_bin <- sys_gcta != ""
  do_gcta_test <- has_plink_files && has_gcta_bin

  # 5. Call multiBAT.combo
  result_prefix <- file.path(out_dir, "test")
  
  res <- multiBAT.combo (
    bim_file         = bim_file,
    map_file         = map_file,
    assoc_pop1       = assoc1,
    assoc_pop2       = assoc2,
    LD_path_pop1     = LD_path1,
    LD_path_pop2     = LD_path2,
    result_prefix    = result_prefix,
    prop             = 0.9,
    af_threshold     = 1,
    run_gcta         = do_gcta_test,
    gcta_bin         = if(do_gcta_test) sys_gcta else "gcta64",
    plink_bfile_pop1 = if(do_gcta_test) plink_path_prefix else NULL,
    plink_bfile_pop2 = if(do_gcta_test) plink_path_prefix else NULL
  )

  # 6. Assert R Analysis Structure
  expect_true(is.list(res))

  expected_cols_multi <- c(
    "Gene","Chr","Start","End","No.SNPs","SNP_start","SNP_end",
    "TopSNP_pop1","TopSNP_Pvalue_pop1","TopSNP_pop2","TopSNP_Pvalue_pop2",
    "No.Eigenvalues","Chisq_mBAT","P_multiBAT",
    "P_multiBAT.combo"
  )
  
  if (nrow(res$multipop) > 0) {
    expect_true(all(expected_cols_multi %in% names(res$multipop)))
  }

  expected_cols_single <- c(
    "Gene","POP","Chr","Start","End","No.SNPs","SNP_start","SNP_end",
    "TopSNP","TopSNP_Pvalue","No.Eigenvalues","Chisq_mBAT","P_mBATcombo"
  )
  if (nrow(res$singlepop) > 0) {
    expect_true(all(expected_cols_single %in% names(res$singlepop)))
  }

  # 7. Assert Standard Output Files
  out_file_multi  <- paste0(result_prefix, ".multipop.txt")
  out_file_single <- paste0(result_prefix, ".singlepop.txt")
  
  if (nrow(res$multipop) > 0) expect_true(file.exists(out_file_multi))
  if (nrow(res$singlepop) > 0) expect_true(file.exists(out_file_single))

  # 8. Assert GCTA Output Files (Conditional)
  if (do_gcta_test) {
    gcta_out1 <- paste0(result_prefix, ".gcta.pop1.gene.assoc.mbat")
    gcta_out2 <- paste0(result_prefix, ".gcta.pop2.gene.assoc.mbat")
    expect_true(file.exists(gcta_out1))
    expect_true(file.exists(gcta_out2))
  }
})