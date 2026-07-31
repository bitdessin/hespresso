#' HOBIT: Detect Shifts in Homeolog Expression Ratios
#'
#' Tests whether homeolog expression ratios differ among experimental
#' conditions in allopolyploid RNA-seq data.
#'
#' `hobit()` is the main entry point for HOBIT analysis. It dispatches the
#' analysis to either the MCMC implementation or the maximum-likelihood
#' implementation according to `method`.
#'
#' Two methods are available:
#'
#' - `"mcmc"` uses the original Stan-based MCMC implementation evaluated in
#'   Sun et al. (2026).
#' - `"mle"` fits the null and alternative models by maximum likelihood,
#'   and then performs a likelihood ratio test.
#'
#' Method-specific arguments supplied through `...` are passed directly to
#' [hobit_mcmc()] or [hobit_mle()].
#'
#' @param x An \linkS4class{ExpMX} object.
#' @param method Character string specifying the inference method. Available
#'   values are `"mcmc"` and `"mle"`. The default is `"mcmc"`.
#' @param ... Additional arguments passed to [hobit_mcmc()] or [hobit_mle()],
#'   according to `method`.
#'
#' @return A data frame with one row per homeolog tuple. The returned columns
#'   depend on the selected inference method. See [hobit_mcmc()] and
#'   [hobit_mle()] for details.
#'
#' @seealso [hobit_mcmc()], [hobit_mle()], \linkS4class{ExpMX}
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
#' x_output_mcmc_1 <- hobit(x)
#' # x_output_mcmc_2 <- hobit(x, method = "mcmc")
#'
#' x_output_mle <- hobit(x, method = "mle")
#'
#' @export
hobit <- function(x, method = c("mcmc", "mle"), ...) {
    method <- match.arg(method)
    switch(
        method,
        mcmc = hobit_mcmc(x, ...),
        mle = hobit_mle(x, ...)
    )
}
