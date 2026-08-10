# Select user-facing output columns unless diagnostic output is requested.
.hobit.select_output <- function(output, .debug = FALSE) {
    if (isTRUE(.debug)) {
        return(output)
    }

    output[, c(
        intersect(c("gene", "pvalue", "qvalue", "Dmax", "LRmax", "RRmax"), names(output)),
        grep("^theta0__", names(output), value = TRUE)
    ), drop = FALSE]
}


#' HOBIT: Detect Shifts in Homeolog Expression Ratios
#'
#' Tests whether homeolog expression ratios differ among experimental
#' conditions in allopolyploid RNA-seq data.
#'
#' `hobit()` is the main entry point for HOBIT analysis. It supports maximum
#' likelihood and MCMC-based inference through the `method` argument, with MLE
#' as the default.
#'
#' Detailed methodological background and usage guidance are available at
#' \url{https://bitdessin.github.io/hespresso/}.
#'
#' @param x An \linkS4class{SeqCountData} object.
#' @param method Character string specifying the inference method. Available
#'   values are `"mcmc"` and `"mle"`. The default is `"mle"`.
#' @param .debug Logical. If `TRUE`, return all statistical test results for
#'   detailed analysis or debugging. If `FALSE`, return only the main
#'   statistical test results.
#' @param ... Additional arguments passed to [hobit_mcmc()] or [hobit_mle()],
#'   according to `method`.
#'
#' @return A data frame with one row per homeolog tuple. By default, the result
#'   contains `gene`, `pvalue`, `qvalue`, `Dmax`, `LRmax`, `RRmax`, and
#'   `theta0__*` columns. Set `.debug = TRUE` to return all method-specific
#'   columns.
#'
#' @seealso [hobit_mcmc()], [hobit_mle()], \linkS4class{SeqCountData}
#'
#' @references
#' Sun J, Sese J, Shimizu KK.
#' A likelihood ratio test for detecting shifts in homeolog expression ratios
#' in allopolyploids.
#' \emph{New Phytologist}. 2026;249(1):603--617.
#' \doi{10.1111/nph.70648}
#'
#' @examples
#' x <- sim_homeolog_counts(n_genes = 10)
#'
#' x_output_mle_1 <- hobit(x)
#' # x_output_mle_2 <- hobit(x, method = "mle")
#'
#' x_output_mcmc <- hobit(x, method = "mcmc")
#'
#' @export
hobit <- function(x, method = c("mle", "mcmc"), .debug = FALSE, ...) {
    method <- match.arg(method)
    switch(
        method,
        mcmc = hobit_mcmc(x, .debug = .debug, ...),
        mle = hobit_mle(x, .debug = .debug, ...)
    )
}
