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

By default, `multiBAT.combo()` obtains population-specific sample sizes from the GWAS summary statistics files using the sample size column:

```text
N
```

The sample size information is used in the weighted ACAT procedure to account for sample-size imbalance between populations and determine the contribution of population-specific association signals.

Alternatively, users can directly specify population-specific sample sizes through:

```r
N_pop1
N_pop2
```

in the `multiBAT.combo()` function. When provided, these user-specified values are used instead of the sample size information extracted from the GWAS summary statistics.

Accurate population-specific sample sizes should be provided because they affect the weighting of population-specific components in the multiBAT-combo association framework.

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
## Output files

`multiBAT.combo()` generates output files according to the selected analysis mode and output options.

The output prefix is controlled by:

```r
result_prefix
```

All output files are named using this prefix.

---

### Default mode

With:

```r
run_gcta = FALSE
print_all_p = FALSE
```

the main output files are:

```text
<result_prefix>.multipop.txt
<result_prefix>.log
```

#### `<result_prefix>.multipop.txt`

This file contains the final multiBAT-combo gene-based association results.

| Column | Description |
|---|---|
| `Gene` | Gene identifier |
| `Chr` | Chromosome |
| `Start` | Gene start position |
| `End` | Gene end position |
| `No.SNPs_pop1` | Number of SNPs available for population 1 |
| `No.SNPs_pop2` | Number of SNPs available for population 2 |
| `No.SNPs_common` | Number of SNPs shared between populations |
| `SNP_start` | Start position of SNPs included in the analysis |
| `SNP_end` | End position of SNPs included in the analysis |
| `TopSNP_pop1` | Most significant SNP in population 1 |
| `TopSNP_Pvalue_pop1` | P value of the top SNP in population 1 |
| `TopSNP_pop2` | Most significant SNP in population 2 |
| `TopSNP_Pvalue_pop2` | P value of the top SNP in population 2 |
| `No.Eigenvalues` | Number of retained eigenvalues in SVD-based multiBAT calculation |
| `P_multiBAT.combo` | Final multiBAT-combo association P value |
| `P_multiBAT.combo_BH` | Benjamini-Hochberg adjusted P value |

#### `<result_prefix>.log`

The log file records analysis parameters, quality control information, fallback events, and analysis summaries.

---

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

Additional columns include:

| Column | Description |
|---|---|
| `P_joint` | P value from the joint multi-ancestry multiBAT analysis |
| `P_pop1` | Population 1-specific mBAT P value |
| `P_pop2` | Population 2-specific mBAT P value |

These columns allow users to evaluate the contribution of joint and population-specific association signals to the final ACAT-combined association.

---

### GCTA mode

With:

```r
run_gcta = TRUE
```

additional GCTA-derived population-specific mBAT results are generated:

```text
<result_prefix>.gcta.pop1.gene.assoc.mbat
<result_prefix>.gcta.pop2.gene.assoc.mbat
```

These files contain population-specific gene-based association results generated by GCTA.

| `No.Eigenvalues` | Number of eigenvalues retained in the mBAT calculation |
| `Chisq_mBAT` | mBAT test statistic |
| `P_mBATcombo` | mBAT-combo association P value generated by GCTA |
| `P_mBAT` | mBAT association P value |
| `Chisq_fastBAT` | fastBAT test statistic |
| `P_fastBAT` | fastBAT association P value |

The final `P_multiBAT.combo` value reported by `multiBAT.combo()` is constructed based on the mBAT-derived association statistics.

Detailed descriptions of the GCTA mBAT-combo method, command options, and output file formats are available in the official GCTA documentation:

https://yanglab.westlake.edu.cn/software/gcta/#mBAT-combo

When:

```r
run_gcta = TRUE
print_all_p = TRUE
```

the final `.multipop.txt` file additionally reports the GCTA-derived population-specific association statistics used in the multiBAT-combo framework.


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

## Acknowledgements

This implementation builds on the original `mBAT-combo` framework developed by Ang Li.
