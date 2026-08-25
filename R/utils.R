#' Internal Helper Utilities and Global Variable Definitions
#'
#' A collection of miscellaneous helper functions and global variable
#' declarations used internally in the multiBAT.combo package.
#'
#' @name utils
#' @keywords internal
NULL


#' Define Global Variables to Satisfy R CMD Check
#'
#' Prevents "no visible binding for global variable" notes during package
#' checks by declaring variable names used in data.table expressions.
#'
#' @name utils
#' @rdname utils
#' @keywords internal
NULL


utils::globalVariables(
    c(
        ".",
        "snp",
        "chr",
        "pos",
        "A1",
        "A2",
        "gene",
        "gene.start",
        "gene.end",
        "pvalue",
        "chisq",

        "ID",
        "id",
        "gene.id",

        "PhysPos",
        "A1Freq",
        "dif_AF",

        "info1",
        "info2",

        "start",
        "end",
        "P",

        "P_multiBAT",
        "P_multiBAT.combo",
        "P_multiBAT.combo_BH",

        "Gene",

        "P_joint",
        "P_pop1",
        "P_pop2",
        "P_value",

        "P_mBAT",
        "P_mBATcombo",

        "No.SNPs",
        "No.SNPs_pop1",
        "No.SNPs_pop2",
        "No.SNPs_common",

        "SNP_start",
        "SNP_end",

        "TopSNP",
        "TopSNP_Pvalue",

        "POP",

        "Chr",
        "Chr_annot",

        "Start",
        "Start_annot",

        "End",
        "End_annot",

        "..req_cols",
            ".N",

    "GCTA_T_mBAT_pop1",
    "GCTA_P_mBAT_pop1",
    "GCTA_T_mBAT_pop2",
    "GCTA_P_mBAT_pop2",

    "T_mBAT_pop1",
    "P_mBAT_pop1",
    "T_mBAT_pop2",
    "P_mBAT_pop2",
    
    "No.Eigenvalues_pop1",
    "No.Eigenvalues_pop2",

    "TopSNP_pop1",
    "TopSNP_Pvalue_pop1",
    "TopSNP_pop2",
    "TopSNP_Pvalue_pop2",

    "GCTA_TopSNP_pop1",
    "GCTA_TopSNP_Pvalue_pop1",
    "GCTA_TopSNP_pop2",
    "GCTA_TopSNP_Pvalue_pop2",

    "N",

    "analysis_status"
    )
)