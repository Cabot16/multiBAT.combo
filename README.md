
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mBAT.multipop

<!-- badges: start -->

<!-- badges: end -->

mBAT.multipop is developed based on the mBAT-combo package
(<https://github.com/Share-AL-work/mBAT.combo>) to enhance the ability
to identify function-related genes by increasing the cross-ethnic sample
size. We further proposed a combined method, MBAT-Global, by integrating
the single-population mBAT analysis results achieved by gcta with the
MBAT.multipop results, to reduce the impact of excessive differences in
available sample sizes among different ethnic groups on the results. In
addition, we offer an extra LD_matrix feature to quickly build available
LD matrices through gctb software and reference plink files

## Installation

You can install the development version of mBAT.multipop like so:

``` r
install.packages("mBAT.multipop_0.0.0.9000.tar.gz", repos = NULL, type = "source")
```

## Example

This section provides basic examples of the `mBAT.multipop` package’s
core functionalities: running the LD matrix preparation, performing the
core gene-based association test (mBAT-Global), and demonstrating the
optional GCTA integration. 

Sample files can be found in the “/inst/extdata/” directory 
### 1. LD Matrix Preparation (Utility) Use
the LD_matrix utility to calculate population-specific LD matrices using
GCTB prior to running the core analysis.

``` r
# Assuming GCTB and PLINK are installed and in PATH
# Set appropriate paths for your reference data
ldm_results <- mBAT.multipop::LD_matrix(
  gctb_bin         = "gctb", 
  map_file         = file.path(data_dir, "your_gene_map.txt"),
  plink_bfile_pop1 = "path/to/pop1/reference_panel",
  plink_bfile_pop2 = "path/to/pop2/reference_panel",
  extend = 50000,
  result_prefix    = file.path(tempdir(), "ld_prep_output")
)

# Check the status of LD calculation tasks
print(ldm_results)
```

### 2. Gene-Based Association Testing (mBAT-Global & Single-pop)

The primary function for joint analysis and single-population mBAT-combo
tests. This requires pre-calculated LD matrix files.

``` r
library(mBAT.multipop)

# Replace these placeholders with actual paths to your example files
data_dir <- system.file("extdata", package = "mBAT.multipop")

results <- mBAT_multipop(
  bim_file          = file.path(data_dir, "your_bfile.bim"),
  map_file          = file.path(data_dir, "your_gene_map.txt"),
  assoc_pop1        = file.path(data_dir, "pop1_gwas.ma"),
  assoc_pop2        = file.path(data_dir, "pop2_gwas.ma"),
  LD_path_pop1      = file.path(data_dir, "pop1_ld_matrices"), # Directory containing .ldm files
  LD_path_pop2      = file.path(data_dir, "pop2_ld_matrices"),
  result_prefix     = file.path(tempdir(), "analysis_results"),
  prop         = 0.9,
  af_threshold = 0.2,
  run_gcta          = FALSE # Run core functionality only
)

# View the joint (multipop) results
head(results$multipop)
```

### 3. Running with GCTA Integration (Optional)

This demonstrates running GCTA’s mBAT-combo test alongside the core
analysis, using PLINK files from the same reference panel.

``` r
# This requires GCTA (gcta64) to be available in your system's PATH.
plink_prefix <- file.path(data_dir, "your_bfile_prefix")

gcta_integrated_results <- mBAT_multipop(
  # ... (all input parameters from Example 1) ...
  run_gcta          = TRUE,
  gcta_bin          = "gcta64",
  plink_bfile_pop1  = plink_prefix, 
  plink_bfile_pop2  = plink_prefix
)

# GCTA results are written to files named: [result_prefix].gcta.pop[X].gene.assoc.mbat
# Check the generated GCTA output files
print(list.files(tempdir(), pattern = "gcta.*mbat", full.names = TRUE))
```
