# multiBAT.combo

<!-- badges: start -->

<!-- badges: end -->

`multiBAT.combo` provides an R implementation of `multiBAT-combo`, a multi-ancestry gene-based association testing framework for complex trait discovery. The method integrates ancestry-specific GWAS summary statistics and linkage disequilibrium (LD) information to jointly evaluate gene-level association signals across populations.

The standard workflow uses two ancestry-specific GWAS summary statistics files together with ancestry-matched LD matrices. `multiBAT.combo()` combines the joint multi-ancestry multiBAT result with ancestry-specific results using a weighted aggregated Cauchy association test (ACAT) to generate the final multiBAT-combo gene-level association P value.

## Workflow overview

A typical analysis consists of the following steps:

1. Prepare population-specific GWAS summary statistics.
2. Prepare ancestry-matched genotype reference panels.
3. Generate population-specific gene-level LD matrices using `LD_matrix()`.
4. Run `multiBAT.combo()` using the GWAS summary statistics and corresponding LD matrices.
5. Examine the final gene-level association results in `<result_prefix>.multipop.txt`.

## Installation

Before installing `multiBAT.combo`, first install the required R package dependencies by running:

```bash
Rscript multiBAT.combo-main/install_dependencies.R
```

The script checks whether the required R packages are available in the current R environment and automatically installs any missing dependencies.

After all dependencies have been installed successfully, install the source package locally with:

```r
install.packages("multiBAT.combo-main", repos = NULL, type = "source")
```

Alternatively:

```r
devtools::install()
```

## External software requirements

`multiBAT.combo()` requires population-specific LD matrices in GCTB format for gene-based association testing.

For users starting from PLINK genotype reference data, the package provides `LD_matrix()` as a convenient way to generate the required LD matrices. `LD_matrix()` uses predefined settings to extract gene-region genotype data from PLINK files and construct gene-level LD matrices using GCTB.

`LD_matrix()` relies on:

* **GCTB**, for constructing the LD matrices;
* **PLINK**, for extracting genotype data for individual gene regions.

Users who already have population-specific LD matrices in the required GCTB format do not need to run `LD_matrix()` and can provide the existing LD matrices directly to `multiBAT.combo()` through `LD_path_pop1` and `LD_path_pop2`.

## Input files

### PLINK BIM file

The PLINK `.bim` file is used to provide the genomic position and allele information for each SNP during gene annotation and allele matching.

The file should follow the standard PLINK BIM format with six columns:

```text
CHR SNP CM BP A1 A2
```

where `CHR` is the chromosome, `SNP` is the SNP identifier, `CM` is the genetic distance, `BP` is the base-pair position, and `A1` and `A2` are the two alleles.

We recommend using the same `.bim` file as that used for the LD reference data. When population-specific LD matrices are constructed from two ancestry-specific reference panels, the `.bim` file from either population can be used, provided that the SNP identifiers, genomic build, and genomic positions are consistent with the input GWAS summary statistics and LD data.

### Gene annotation file

The gene annotation file is used to define the genomic coordinates of each gene for SNP-to-gene mapping.

The file should contain four tab-separated columns without a header:

```text
chr start end gene.id
```

For example:

```text
1	65419	71585	OR4F5
1	134901	139379	AL627309.1
1	367640	368634	OR4F29
```
The genes listed in this file define the set of genes evaluated throughout the downstream analysis. Therefore, a simple way to restrict or customise the genes included in the analysis is to modify the gene list in this annotation file.

The chromosome and genomic coordinates should use the same genome build as the PLINK BIM file, GWAS summary statistics, and LD reference data.

### GWAS summary statistics

Each population-specific GWAS summary statistics file should contain the following columns:

```text
snp A1 A2 AF1 beta se P N
```

The sample size column `N` is required. Accurate population-specific sample sizes (or effect sample size for binary traits) should be provided because sample sizes are used to determine the contribution of population-specific association signals in the weighted ACAT combination.

Although `multiBAT.combo()` performs internal numerical validation of GWAS fields, users should perform standard GWAS quality control before analysis. Recommended preprocessing includes: 
    removal of SNPs with missing or invalid effect sizes; 
    removal of SNPs with invalid standard errors; 
    allele frequency filtering; 
    imputation quality filtering where applicable; 
    removal of duplicated SNP identifiers.

For more detailed guidance on GWAS summary statistic formatting and preparation, users may refer to the GCTA documentation (https://yanglab.westlake.edu.cn/software/gcta/#DataResource).

## LD matrix preparation

For users starting from PLINK genotype reference data, `LD_matrix()` can be used to generate the population-specific gene-level LD matrices required by `multiBAT.combo()`.

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

The `extend` argument specifies the number of base pairs added to both the upstream and downstream boundaries of each gene. For example:

```r
extend = 50000
```

includes SNPs located within 50 kb upstream and 50 kb downstream of the annotated gene region.

For each gene and population, `LD_matrix()` first extracts the corresponding genomic region from the PLINK reference panel and then uses GCTB to construct a full LD matrix. The underlying procedure is equivalent to:

```bash
plink --bfile <PLINK_PREFIX> \
      --chr <CHR> \
      --from-bp <START_EXTENDED> \
      --to-bp <END_EXTENDED> \
      --make-bed \
      --out <GENE>.<POPULATION>

gctb --bfile <GENE>.<POPULATION> \
     --make-full-ldm \
     --out <GENE>.<POPULATION>
```

where `<START_EXTENDED>` and `<END_EXTENDED>` are determined from the gene coordinates and the specified `extend` value.

If no SNPs are available within a gene region after extraction, that gene is skipped for the corresponding population.

`LD_matrix()` applies this procedure automatically to all genes in the annotation file and to both population-specific reference panels.

As a reference for computational requirements, using 1000 Genomes Project Phase 3 reference data and approximately 20,000 genes, LD matrix construction required approximately 90 minutes per ancestry when four parallel tasks were used. 

Users who already have compatible population-specific GCTB full LD matrices do not need to run `LD_matrix()` and can provide the existing LD matrix directories directly to `multiBAT.combo()`.

LD matrices for different populations should be stored in separate directories. For each gene, the LD files should follow the naming pattern:

```
<gene>.*ldm.full.info
<gene>.*ldm.full.bin
```

where <gene> must exactly match the gene identifier used in the gene annotation file.

## Gene-based association testing

The main analysis is performed using:

```r
multiBAT.combo()
```

A standard analysis using the internal population-specific mBAT implementation can be run as follows:

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
    run_gcta = FALSE,
    print_all_p = FALSE
)
```

### Arguments

| Argument            | Description                                                                                                                                                                                                                                                                           |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `bim_file`          | Path to the PLINK BIM file used to provide SNP genomic positions and allele information.                                                                                                                                                                                              |
| `map_file`          | Path to the headerless four-column gene annotation file containing chromosome, gene start position, gene end position, and gene identifier.                                                                                                                                           |
| `assoc_pop1`        | Path to the GWAS summary statistics file for population 1.                                                                                                                                                                                                                            |
| `assoc_pop2`        | Path to the GWAS summary statistics file for population 2.                                                                                                                                                                                                                            |
| `LD_path_pop1`      | Directory containing the population 1 GCTB full LD matrices and corresponding `.info` files.                                                                                                                                                                                          |
| `LD_path_pop2`      | Directory containing the population 2 GCTB full LD matrices and corresponding `.info` files.                                                                                                                                                                                          |
| `result_prefix`     | Full output path and filename prefix. For example, `result_prefix = "results/BMI"` generates files beginning with `results/BMI`.                                                                                                                                                      |
| `extend`            | Number of base pairs added to both the upstream and downstream boundaries of each annotated gene when assigning SNPs to genes. The default is `50000`, corresponding to 50 kb on each side. The same gene-region definition should be used when preparing the LD matrices.            |
| `prop`              | Proportion of eigen-variance retained in the SVD-based calculation. The value must be greater than 0 and no greater than 1. The default is `0.9`.                                                                                                                                |
| `af_threshold`      | Maximum allowed absolute difference between the effect-allele frequency in the GWAS summary statistics and the corresponding allele frequency in the LD reference after allele harmonisation. SNPs are retained when the difference is smaller than this value. The default is `0.2`. |
| `run_gcta`          | Controls how the population-specific mBAT analyses are performed. With `FALSE`, population-specific mBAT statistics are calculated internally in R. With `TRUE`, GCTA is used for the two population-specific mBAT analyses. The default is `FALSE`.                                  |
| `print_all_p`       | Controls whether intermediate joint and population-specific test statistics and P values are included in the main output file. The default is `FALSE`.                                                                                                                                |
| `check_gwas_values` | Controls whether numeric GWAS quality control is performed before gene-based analysis. The default is `TRUE`.                                                                                                                                                                         |
| `N_pop1`            | Optional population 1 sample size supplied directly by the user. When provided, this value is used instead of the sample size obtained from the `N` column of the population 1 GWAS summary statistics.                                                                               |
| `N_pop2`            | Optional population 2 sample size supplied directly by the user. When provided, this value is used instead of the sample size obtained from the `N` column of the population 2 GWAS summary statistics.                                                                               |
| `gcta_bin`          | GCTA executable name available in the system `PATH`, or the full path to the GCTA executable. Used only when `run_gcta = TRUE`. The default is `"gcta64"`.                                                                                                                            |
| `plink_bfile_pop1`  | PLINK binary reference panel prefix for population 1. Required only when `run_gcta = TRUE`. The prefix should be supplied without the `.bed`, `.bim`, or `.fam` extension.                                                                                                            |
| `plink_bfile_pop2`  | PLINK binary reference panel prefix for population 2. Required only when `run_gcta = TRUE`. The prefix should be supplied without the `.bed`, `.bim`, or `.fam` extension.                                                                                                            |

### Population-specific mBAT implementation

The recommended setting is:

```r
run_gcta = FALSE
```

With this setting, population-specific mBAT statistics are calculated internally in R using the corresponding population-specific LD matrices.

Alternatively:

```r
run_gcta = TRUE
```

uses GCTA to calculate the population-specific mBAT statistics. In addition to mBAT statistics, GCTA also calculates mBAT-combo and fastBAT results for each population. In this mode, the following additional arguments must be provided:

```r
gcta_bin
plink_bfile_pop1
plink_bfile_pop2
```

The joint multi-ancestry multiBAT statistic is calculated internally in both modes.

The final multiBAT-combo association combines the following three gene-level P values:

```text
P_multiBAT
P_mBAT_pop1
P_mBAT_pop2
```

using a weighted aggregated Cauchy association test (ACAT).

When run_gcta = TRUE, the additional GCTA-derived mBAT-combo and fastBAT results are retained as supplementary population-specific outputs but are not used in the final multiBAT-combo weighted ACAT combination.

### Sample size information

Population-specific sample sizes are used to determine the contribution of the population-specific components in the weighted ACAT combination.

By default, `multiBAT.combo()` obtains the sample size for each population from the `N` column of the corresponding GWAS summary statistics file.

Alternatively, population-level sample sizes can be supplied directly using:

```r
N_pop1
N_pop2
```

When `N_pop1` and/or `N_pop2` are provided, the supplied values are used instead of the sample size information obtained from the corresponding GWAS summary statistics file.

Accurate population-specific sample sizes are important because the multiBAT-combo weighting scheme accounts for sample-size imbalance between populations. Greater sample-size imbalance reduces the relative contribution of the joint multi-ancestry component and increases the contribution of the population-specific components.


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
| `SNP_start` | Start position of SNPs included in the joint analysis |
| `SNP_end` | End position of SNPs included in the joint analysis |
| `TopSNP_pop1` | Most significant SNP in population 1 |
| `TopSNP_Pvalue_pop1` | P value of the top SNP in population 1 |
| `TopSNP_pop2` | Most significant SNP in population 2 |
| `TopSNP_Pvalue_pop2` | P value of the top SNP in population 2 |
| `No.Eigenvalues` | Number of eigenvalues retained for the joint multi-ancestry multiBAT calculation |
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
| `T_mBAT_pop1` | mBAT test statistic for population 1 |
| `T_mBAT_pop2` | mBAT test statistic for population 2 |
| `T_multiBAT` | Joint multiBAT test statistic |
| `P_multiBAT` | P value from the joint multi-ancestry multiBAT analysis |
| `P_mBAT_pop1` | Population 1-specific mBAT P value |
| `P_mBAT_pop2` | Population 2-specific mBAT P value |
| `No.Eigenvalues_pop1` | Number of eigenvalues retained for population 1 mBAT calculation |
| `No.Eigenvalues_pop2` | Number of eigenvalues retained for population 2 mBAT calculation |

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

These files are directly generated by GCTA and follow the standard mBAT-combo output format.

| Column | Description |
|---|---|
| `Gene` | Gene identifier |
| `Chr` | Chromosome |
| `Start` | Gene start position |
| `End` | Gene end position |
| `No.SNPs` | Number of SNPs included in the gene-based analysis |
| `SNP_start` | First SNP included in the analysis window |
| `SNP_end` | Last SNP included in the analysis window |
| `TopSNP` | Most significant SNP within the gene region |
| `TopSNP_Pvalue` | P value of the top SNP |
| `No.Eigenvalues` | Number of eigenvalues retained in the mBAT calculation |
| `Chisq_mBAT` | mBAT test statistic |
| `P_mBATcombo` | mBAT-combo association P value generated by GCTA |
| `P_mBAT` | mBAT association P value |
| `Chisq_fastBAT` | fastBAT test statistic |
| `P_fastBAT` | fastBAT association P value |

Detailed descriptions of the GCTA mBAT-combo method, command options, and output file formats are available in the official GCTA documentation:

https://yanglab.westlake.edu.cn/software/gcta/#mBAT-combo

In addition to the original GCTA output files, the final `<result_prefix>.multipop.txt` file also reports the population-specific mBAT, fastBAT, and mBAT-combo results generated by GCTA.

These include:

```text
T_mBAT_pop1
P_mBAT_pop1
P_fastBAT_pop1
P_mBATcombo_pop1
T_mBAT_pop2
P_mBAT_pop2
P_fastBAT_pop2
P_mBATcombo_pop2

The final `P_multiBAT.combo` value reported by `multiBAT.combo()` is constructed based on mBAT-derived association statistics and combined with the joint multi-ancestry component using the weighted ACAT procedure.

## Testing the installation

Users can verify that the package has been installed correctly by running the bundled minimal example workflow:

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
This example uses the test data bundled with the package to verify that multiBAT.combo() can run successfully after installation.

## Acknowledgements

This implementation builds on the original `mBAT-combo` framework developed by Ang Li (https://github.com/Share-AL-work/mBAT.combo).
