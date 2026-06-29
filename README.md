
# multiBAT.combo

<!-- badges: start -->

<!-- badges: end -->

`multiBAT.combo` provides an R implementation of `multiBAT-combo`, a
multi-ancestry gene-based association testing framework for complex
trait discovery. The method integrates GWAS summary statistics and
population-specific LD information from multiple populations to improve
gene-level discovery while avoiding reliance on a single LD reference
for multi-ancestry data.

The package supports allele harmonisation, gene-level SNP mapping,
per-gene LD matrix processing, block-diagonal joint LD construction,
SVD-based multiBAT testing, and optional integration of
population-specific GCTA `mBAT-combo` results through a weighted
ACAT-based combination framework.

## Naming convention

This package uses the following naming convention:

- `multiBAT` refers to the joint multi-population gene-based association
  method.
- `multiBAT.combo` refers to the package name, main function name, and
  code-level P-value combination workflow.
- `multiBAT-combo` is used for output and display labels.
- `mBAT-combo` refers to the original single-population method
  implemented in GCTA.

## Installation

Install the source package locally with:

``` r
install.packages("multiBAT.combo_0.0.1.2000.tar.gz", repos = NULL, type = "source")
```

Or install from a local development directory with:

``` r
devtools::install()
```

## External software requirements

Some workflows require external command-line tools:

- `LD_matrix()` requires GCTB and PLINK.
- `multiBAT.combo(run_gcta = TRUE)` requires GCTA.

These tools should either be available in `PATH`, or their executable
paths should be provided explicitly.

## Input files

### PLINK BIM file

The PLINK `.bim` file should contain six columns:

``` text
chr snp bp pos A1 A2
```

### Gene annotation file

The gene annotation file should not contain a header and should include
four columns:

``` text
chr start end gene.id
```

The `extend` argument expands each gene region upstream and downstream
in base pairs.

### GWAS summary statistics

Each population-specific GWAS summary file should contain the following
columns:

``` text
snp A1 A2 AF1 beta se P N
```

The two GWAS files should use the same column order.

### LD matrix directories

`LD_path_pop1` and `LD_path_pop2` should point to directories containing
per-gene GCTB LD files:

``` text
<gene>.ldm.full.bin
<gene>.ldm.full.info
```

The `.info` file should contain a header and include SNP-level LD
reference information.

## LD matrix preparation

`LD_matrix()` can be used to generate per-gene LD matrices from PLINK
reference panels using GCTB.

``` r
library(multiBAT.combo)

ldm_results <- LD_matrix(
    gctb_bin = "gctb",
    map_file = "path/to/gene_annotation.txt",
    plink_bfile_pop1 = "path/to/pop1_reference_prefix",
    plink_bfile_pop2 = "path/to/pop2_reference_prefix",
    extend = 50000,
    result_prefix = "path/to/ld_output",
    thread = 1
)

print(ldm_results)
```

This creates population-specific LD directories:

``` text
<result_prefix>/pop1/
<result_prefix>/pop2/
```

Each directory contains per-gene `.ldm.full.bin` and `.ldm.full.info`
files.

## Gene-based association testing

The main function is `multiBAT.combo()`. It performs joint gene-based
association testing across two populations using pre-computed
population-specific LD matrices.

``` r
library(multiBAT.combo)

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

head(res$multipop)
```

When `run_gcta = FALSE`, the final `P_multiBAT.combo` column is set
equal to `P_multiBAT`.

## Optional GCTA mBAT-combo integration

When `run_gcta = TRUE`, the package runs population-specific GCTA
`mBAT-combo` analyses and combines the joint multiBAT P value with
population-specific GCTA P values using weighted ACAT.

``` r
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
    run_gcta = TRUE,
    gcta_bin = "gcta64",
    plink_bfile_pop1 = "path/to/pop1_reference_prefix",
    plink_bfile_pop2 = "path/to/pop2_reference_prefix"
)
```

The final multiBAT-combo P value combines:

``` text
P_multiBAT
P_pop1
P_pop2
```

using sample-size-imbalance-aware weights.

## Output files

The main workflow writes:

``` text
<result_prefix>.multipop.txt
<result_prefix>.singlepop.txt
<result_prefix>.log
```

If `run_gcta = TRUE`, it also writes:

``` text
<result_prefix>.gcta.pop1.gene.assoc.mbat
<result_prefix>.gcta.pop2.gene.assoc.mbat
```

### Multi-population output

The multi-population output contains:

``` text
Gene
Chr
Start
End
No.SNPs
SNP_start
SNP_end
TopSNP_pop1
TopSNP_Pvalue_pop1
TopSNP_pop2
TopSNP_Pvalue_pop2
No.Eigenvalues
Chisq_mBAT
P_multiBAT
P_multiBAT.combo
```

### Single-population fallback output

The single-population fallback output contains genes that are only
available in one population-specific GWAS input:

``` text
Gene
POP
Chr
Start
End
No.SNPs
SNP_start
SNP_end
TopSNP
TopSNP_Pvalue
No.Eigenvalues
Chisq_mBAT
P_mBATcombo
```

This fallback is separate from GCTA `mBAT-combo`. When
`run_gcta = TRUE`, the population-specific evidence used for final ACAT
combination is taken from the GCTA outputs.

## Testing the installation

After installing `multiBAT.combo`, users can test whether the package
was installed successfully by running the example workflow included with
the package.

``` r
library(multiBAT.combo)

data_dir <- system.file("extdata", package = "multiBAT.combo")
stopifnot(dir.exists(data_dir))

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

str(res)
head(res$multipop)

file.exists(paste0(result_prefix, ".multipop.txt"))
```

A successful installation should return a list containing `multipop` and
`singlepop` results. If at least one multi-population gene is
successfully analysed, the following output file should also be
generated:

``` text
<result_prefix>.multipop.txt
```

## Optional GCTA test

If GCTA and the required PLINK reference files are available, users can
also test the optional GCTA integration:

``` r
data_dir <- system.file("extdata", package = "multiBAT.combo")
plink_prefix <- file.path(data_dir, "1000G_eur_unrel_ukbcom22")
result_prefix <- file.path(tempdir(), "multiBAT_combo_gcta_test")

res_gcta <- multiBAT.combo(
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
    run_gcta = TRUE,
    gcta_bin = "gcta64",
    plink_bfile_pop1 = plink_prefix,
    plink_bfile_pop2 = plink_prefix
)

head(res_gcta$multipop)
file.exists(paste0(result_prefix, ".gcta.pop1.gene.assoc.mbat"))
file.exists(paste0(result_prefix, ".gcta.pop2.gene.assoc.mbat"))
```

## Notes

- Alleles are harmonised against the LD reference panel before
  gene-level testing.
- SNPs with large allele-frequency differences can be filtered using
  `af_threshold`.
- Gene-specific file or QC failures are skipped and recorded in the log,
  rather than stopping the full analysis.
- Strand-complement allele matching is not performed and should be
  handled before running this package if required.

## Acknowledgements

The implementation builds on the original `mBAT-combo` framework
developed by Ang Li. Parts of the testing workflow were adapted from the
test scripts of the original `mBAT.combo` R package and modified for the
`multiBAT.combo` workflow.
