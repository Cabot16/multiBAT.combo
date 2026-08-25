# multiBAT.combo

<!-- badges: start -->
<!-- badges: end -->

`multiBAT.combo` provides an R implementation of `multiBAT-combo`, a multi-ancestry gene-based association testing framework for complex trait discovery. The method integrates GWAS summary statistics and population-specific linkage disequilibrium (LD) information from multiple populations to improve gene-level discovery.

## Installation

Install the source package locally with:

```r
install.packages("multiBAT.combo_0.0.1.3000.tar.gz", repos = NULL, type = "source")
```

Alternatively:

```r
devtools::install()
```

## External software requirements

`multiBAT.combo()` requires pre-computed population-specific LD matrices for gene-level association testing. These LD matrices should be generated using `LD_matrix()`, which relies on external software:

- `LD_matrix()` requires GCTB and PLINK.

After LD matrices have been prepared, the main gene-based association testing workflow can be performed entirely within R.

The recommended usage is:

```r
run_gcta = FALSE
```

In this mode, `multiBAT.combo()` calculates population-specific mBAT statistics internally in R and combines them with the joint multiBAT statistic using the multiBAT-combo framework. This mode provides a fully integrated workflow without requiring GCTA.

Alternatively, users can specify:

```r
run_gcta = TRUE
```

to enable optional GCTA mBAT-combo integration. In this mode, population-specific mBAT-combo analyses are performed using the external GCTA implementation before being combined with the joint multiBAT result.

Because GCTA is implemented in C and may provide faster computation in some large-scale analyses, this option can be useful when computational efficiency is a priority. However, we recommend the internal R implementation (`run_gcta = FALSE`) as the standard workflow because it avoids additional software dependencies and provides a unified analysis pipeline within `multiBAT.combo()`.

## Input files

### PLINK BIM file

The PLINK `.bim` file should contain six columns:

```text
chr snp bp pos A1 A2
```

### Gene annotation file

The annotation file should contain four tab-separated columns without a header:

```text
chr start end gene.id
```

The `extend` argument expands gene regions upstream and downstream in base pairs. A value of 50000 is recommended.

### GWAS summary statistics

Each population-specific GWAS summary statistics file should contain:

```text
snp A1 A2 AF1 beta se P N
```

The sample size column `N` is required. Accurate population-specific sample sizes should be provided because sample sizes are used to determine the contribution of population-specific association signals in the weighted ACAT combination.

Although `multiBAT.combo()` performs internal numerical validation of GWAS fields, users should perform standard GWAS quality control before analysis. Recommended preprocessing includes:

- removal of SNPs with missing or invalid effect sizes;
- removal of SNPs with invalid standard errors;
- allele frequency filtering;
- imputation quality filtering where applicable;
- removal of duplicated SNP identifiers.

The internal checks are designed to prevent computational errors and do not replace standard GWAS quality control procedures.

## LD matrix preparation

`LD_matrix()` generates population-specific LD matrices using GCTB.

Example:

```r
ldm_results <- LD_matrix(
    gctb_bin = "gctb",
    map_file = "path/to/gene_annotation.txt",
    plink_bfile_pop1 = "path/to/pop1_reference",
    plink_bfile_pop2 = "path/to/pop2_reference",
    extend = 50000,
    result_prefix = "path/to/output",
    thread = 1
)
```

## Gene-based association testing

The main function is:

```r
multiBAT.combo()
```

Example:

```r
res <- multiBAT.combo(
    bim_file = "path/to/reference.bim",
    map_file = "path/to/gene_annotation.txt",
    assoc_pop1 = "path/to/pop1_gwas.txt",
    assoc_pop2 = "path/to/pop2_gwas.txt",
    LD_path_pop1 = "path/to/pop1_ld",
    LD_path_pop2 = "path/to/pop2_ld",
    result_prefix = "path/to/output/result",
    extend = 50000,
    prop = 0.9,
    af_threshold = 0.2,
    run_gcta = FALSE
)
```

When `run_gcta = FALSE`, the package calculates population-specific mBAT results internally in R.

The final association combines:

```text
P_multiBAT
P_mBAT_pop1
P_mBAT_pop2
```

using weighted ACAT.

## Print all P values

By default:

```r
print_all_p = FALSE
```

Only the final multiBAT-combo association statistics are written.

When:

```r
print_all_p = TRUE
```

additional gene-level component statistics are reported, including intermediate population-specific and joint association P values.

This option is recommended for users who want to inspect the contribution of each component test.

## Sample size information

The GWAS sample size column:

```text
N
```

is required.

The sample size information is used in the weighted ACAT procedure to account for sample-size imbalance between populations. Users should provide accurate population-specific sample sizes.

## Optional GCTA mBAT-combo integration

When:

```r
run_gcta = TRUE
```

the package runs GCTA separately for the two populations and combines the resulting population-specific GCTA mBAT-combo P values with the joint multiBAT result.

Required arguments include:

```r
gcta_bin
plink_bfile_pop1
plink_bfile_pop2
```

This mode requires GCTA installation and is intended for advanced users.

## Output files

### Default mode

With:

```r
run_gcta = FALSE
print_all_p = FALSE
```

the main outputs are:

```text
<result_prefix>.multipop.txt
<result_prefix>.log
```

### Printing all P values

With:

```r
print_all_p = TRUE
```

the file:

```text
<result_prefix>.multipop.txt
```

contains additional component-level association statistics.

### GCTA mode

With:

```r
run_gcta = TRUE
```

additional files are generated:

```text
<result_prefix>.gcta.pop1.gene.assoc.mbat
<result_prefix>.gcta.pop2.gene.assoc.mbat
```

## Testing the installation

Users can verify installation using the bundled example workflow:

```r
library(multiBAT.combo)

data_dir <- system.file("extdata", package = "multiBAT.combo")
result_prefix <- file.path(tempdir(), "multiBAT_combo_test")

res <- multiBAT.combo(
    bim_file = file.path(data_dir, "1000G_eur_unrel_ukbcom22.bim"),
    map_file = file.path(data_dir, "hg19_v40_glist_Pat_v1.txt"),
    assoc_pop1 = file.path(data_dir, "test.ma"),
    assoc_pop2 = file.path(data_dir, "test2.ma"),
    LD_path_pop1 = file.path(data_dir, "POP1_LD_Chr22"),
    LD_path_pop2 = file.path(data_dir, "POP1_LD_Chr22"),
    result_prefix = result_prefix,
    extend = 50000,
    prop = 0.9,
    af_threshold = 1,
    run_gcta = FALSE
)

head(res$multipop)
```

This test uses the internal R implementation and does not require GCTA.

## Notes

- Alleles are harmonised against the LD reference panel before gene-level testing.
- SNPs with large allele-frequency differences can be filtered using `af_threshold`.
- Gene-specific failures are recorded in the log.
- `run_gcta = FALSE` does not require GCTA.

## Acknowledgements

This implementation builds on the original `mBAT-combo` framework developed by Ang Li.
