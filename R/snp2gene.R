#' Map SNPs to Genes with Windowed Annotation
#'
#' Annotate single‐nucleotide polymorphisms (SNPs) to genes based on their genomic
#' positions. This function extends each gene’s boundaries upstream and downstream
#' by specified kilobase windows, then assigns SNPs falling within those extended
#' regions. It returns a list of SNP sets per gene, and—unless only the sets are
#' requested—a mapping data frame with adjusted positions.
#'
#' @param info_snp Data frame of SNP annotations. Must contain columns:
#'   \code{snp} (SNP ID), \code{chr} (chromosome), and \code{pos} (base-pair position).
#' @param info_gene Data frame of gene boundaries. Must contain columns:
#'   \code{gene} (gene ID), \code{chr} (chromosome), \code{start} (gene start in bp),
#'   and \code{end} (gene end in bp).
#' @param window_start Integer. Upstream extension length in kilobases to apply to each gene’s
#'   start coordinate. Defaults to \code{50L}.
#' @param window_end Integer. Downstream extension length in kilobases to apply to each gene’s
#'   end coordinate. Defaults to \code{50L}.
#' @param only_sets Logical. If \code{TRUE}, return only the list of SNP sets per gene.
#'   If \code{FALSE}, also return a data frame mapping SNPs to genes with adjusted positions.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{sets}{Named list of character vectors. Each element holds the SNP IDs mapped to that gene.}
#'   \item{map}{Data frame of SNP‐to‐gene mappings (only if \code{only_sets = FALSE}). Columns:
#'     \code{snp}, \code{chr}, \code{pos} (SNP position),
#'     \code{gene}, \code{gene.start}, \code{gene.end},
#'     \code{gene.adj.start} (start minus \code{window_start} kb),
#'     \code{gene.adj.end}   (end plus \code{window_end} kb).}
#' }
#'
#' @export
#'
#' @examples
#' # Example SNP annotation
#' snps <- data.frame(
#'   snp = c("rs1", "rs2", "rs3"),
#'   chr = c("1", "1", "2"),
#'   pos = c(100000, 150000, 200000),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Example gene annotation
#' genes <- data.frame(
#'   gene  = c("GENE1", "GENE2"),
#'   chr   = c("1", "2"),
#'   start = c(120000, 190000),
#'   end   = c(140000, 210000),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Map SNPs using a 20 kb upstream and downstream window
#' snp2gene(
#'   info_snp      = snps,
#'   info_gene     = genes,
#'   window_start  = 20L,
#'   window_end    = 20L,
#'   only_sets     = FALSE
#' )
snp2gene <- function(info_snp, info_gene,
                     window_start = 50L,
                     window_end   = 50L,
                     only_sets    = FALSE) {

  info_snp_name <- deparse(substitute(info_snp))

  df_check(info_snp)
  df_check(info_gene)
  columns_check(info_snp,  c("snp",  "chr",  "pos"))
  columns_check(info_gene, c("gene", "chr",  "start", "end"))
  nonneg_num_check(window_start, "window_start")
  nonneg_num_check(window_end,   "window_end")

  if (anyDuplicated(info_snp$snp) > 0L) {
    stop("SNP IDs in '", info_snp_name, "' must be unique.", call. = FALSE)
  }

  snpdat <- data.table::data.table(
    snp   = as.character(info_snp$snp),
    chr   = as.character(info_snp$chr),
    start = as.integer(info_snp$pos),
    end   = as.integer(info_snp$pos),
    key   = c("chr", "start", "end")
  )

  genedat <- data.table::data.table(
    gene       = as.character(info_gene$gene),
    chr        = as.character(info_gene$chr),
    start      = as.integer(
                   pmax(1L, info_gene$start - window_start * 1000L)
                 ),
    end        = as.integer(info_gene$end + window_end * 1000L),
    gene.start = as.integer(info_gene$start),
    gene.end   = as.integer(info_gene$end),
    key        = c("chr", "start", "end")
  )

  mapped <- data.table::foverlaps(snpdat, genedat, type = "within")

  snp_sets <- lapply(
    split(mapped[!is.na(gene), .(snp, gene)], by = "gene", keep.by = FALSE),
    function(df) unique(df$snp)
  )

  if (only_sets) {
    return(list(sets = snp_sets))
  }

  data.table::setnames(
    mapped,
    old = c("start", "end", "i.start"),
    new = c("gene.adj.start", "gene.adj.end", "pos")
  )

  vkeep <- c("snp", "chr", "pos",
             "gene", "gene.start", "gene.end",
             "gene.adj.start", "gene.adj.end")

  map_df <- data.table::setDF(
    data.table::setorder(mapped[, vkeep, with = FALSE],
                         chr, gene.start, gene.end, pos, na.last = TRUE)
  )

  return(list(sets = snp_sets, map = map_df))
}
