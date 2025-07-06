#' Normalize RNA-Seq Read Count Data
#'
#' Normalizes RNA-Seq read count data.
#'
#' This function first calculates gene expression by summing homeolog expression
#' across all subgenomes. Then, it normalizes the gene expression values using
#' the specified method. After normalization, it redistributes the normalized
#' gene-level expression back to homeologs according to the original homeolog
#' expression ratios.
#'
#' @param x An \linkS4class{ExpMX} class object.
#' @param method Character. Specifies the normalization method to use.  
#'          Supported methods are
#'          `'tmm'` (TMM: trimmed mean of M-values; used by \pkg{edgeR}),
#'          `'med'` (median ratio method; used by \pkg{DESeq2}),
#'          and `'cpm'` (CPM: counts per million).
#'          Note that CPM is generally used for visualizations and
#'          is **NOT RECOMMENDED** for statistical analysis,
#'          as it may distort variance structure.
#' @param round Logical. If `TRUE`, rounds normalized expression values
#'          to integers using the `round()` function.
#'
#' @return An \linkS4class{ExpMX} object with normalized read counts.
#'
#' @examples
#' x <- sim_homeolog_counts(n_genes = 100)
#' x <- norm_counts(x)
#'
#' @export
norm_counts <- function(x, method = c('tmm', 'med', 'cpm'), round = FALSE) {
    method <- match.arg(method)
    
    gexp <- .calc_gexp(x)
    
    her <- vector('list', length(x@data))
    for (i in seq_along(x@data)) {
        her[[i]] <- x@data[[i]] / gexp
    }

    gexp_norm <- switch(method,
        'cpm' = .norm_cpm(gexp),
        'tmm' = .norm_tmm(gexp),
        'med' = .norm_med(gexp, x@exp_design$group)
    )
    for (i in seq_along(x@data)) {
        hexp <- gexp_norm * her[[i]]
        if (round)
            hexp <- round(hexp)
        x@data[[i]] <- hexp
        x@data[[i]][is.na(x@data[[i]])] <- 0
    }
    x
}


.norm_cpm <- function(gexp) {
    sweep(gexp, 2, 1e6 / colSums(gexp), '*')
}


#' @importFrom edgeR DGEList normLibSizes cpm
.norm_tmm <- function(gexp) {
    y <- DGEList(counts = gexp)
    y <- normLibSizes(y)
    ef_libsizes <- y$samples$norm.factors * colSums(gexp)
    sweep(gexp, 2, mean(ef_libsizes) / ef_libsizes, "*")
}


#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts
.norm_med <- function(gexp, group) {
    y <- suppressMessages(DESeqDataSetFromMatrix(
                                countData = gexp,
                                colData = data.frame(condition = group),
                                design = ~ condition))
    y <- estimateSizeFactors(y)
    counts(y, normalized = TRUE)
}







