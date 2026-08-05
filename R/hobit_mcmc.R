#' HOBIT Using MCMC
#'
#' Detects changes in homeolog expression ratios using an MCMC-based
#' implementation of HOBIT.
#'
#' HOBIT uses a likelihood ratio test (LRT) to compare two hierarchical models:
#' a full model that allows homeolog expression ratios to vary among conditions
#' and a reduced model that assumes constant ratios across all conditions.
#' Two fitting modes are available, depending on how the reduced model is handled:
#'
#' \itemize{
#'   \item `"fast"` uses the original approximate workflow, in which the reduced
#'     model likelihood is constructed from posterior samples of the
#'     full model. This reduces computational time but may lead to biased
#'     p-values and q-values.
#'     This strategy is described and evaluated in Sun et al. (2026).
#'   \item `"strict"` samples the null and alternative models independently.
#' }
#'
#' @param x An \linkS4class{ExpMX} object containing homeolog RNA-seq count
#'   data. The required count representation may depend on `mode`.
#' @param mode Character. MCMC fitting mode. Either `"fast"` or `"strict"`.
#' @param use_Dirichlet Logical. Whether to use data-derived Dirichlet priors
#'   for homeolog expression ratios.
#' @param no_replicate Logical. Whether to use the dispersion-estimation
#'   procedure intended for data without biological replicates.
#' @param eps Positive numeric value used as a lower bound for expected
#'   expression values.
#' @param dist Character. Distribution used for count modeling. Available
#'   choices are `"NB"` and `"ZINB"`. Currently, `"strict"` mode supports
#'   only `"NB"`.
#' @param n_threads Positive integer specifying the number of parallel R
#'   workers used across homeolog tuples.
#' @param parallel_chains Positive integer specifying the number of Stan chains
#'   run in parallel for each model.
#' @param .debug Logical. If `TRUE`, return all statistical test results for
#'   detailed analysis or debugging. If `FALSE`, return only the main
#'   statistical test results.
#' @param ... Additional arguments passed to
#'   \code{\link[cmdstanr]{model-method-sample}}, such as `chains`,
#'   `iter_warmup`, `iter_sampling`, `thin`, and `seed`.
#'
#' @seealso [hobit()], [hobit_mle()]
#'
#' @return By default, a \code{data.frame} with one row per homeolog containing
#'   `gene`, `pvalue`, `qvalue`, `Dmax`, `LRmax`, `RRmax`, and `theta0__*`
#'   columns. When `.debug = TRUE`, the data frame additionally contains:
#'      \itemize{
#'          \item `pvalue`: p-value from the MCMC likelihood ratio test (LRT).
#'          \item `qvalue`: Benjamini–Hochberg adjusted `pvalue`.
#'          \item `lrt`: Likelihood ratio statistic.
#'          \item `df`: Degrees of freedom used for the LRT.
#'          \item `status`: Fitting status for the homeolog.
#'          \item `logLik_H0`: Log-likelihood of the reduced model.
#'          \item `logLik_H1`: Log-likelihood of the full model.
#'          \item `lrt_negative_prob`: Proportion of posterior LRT draws below zero.
#'          \item `D__$__*`: Difference in expression ratios for subgenome `$`
#'                between the two groups represented by `*`.
#'          \item `LR__$__*`: Log ratio of expression ratios between the
#'                subgenomes and groups represented by `$` and `*`.
#'          \item `RR__$__*`: Ratio of expression ratios between the
#'                subgenomes and groups represented by `$` and `*`.
#'          \item `Dmax`: Maximum absolute difference in expression ratios (`D__$__*`)
#'                observed across all subgenomes and group comparisons.
#'          \item `LRmax`: Maximum absolute log ratio in expression ratios (`LR__$__*`)
#'                observed across all subgenome-pair and group comparisons.
#'          \item `RRmax`: Maximum ratio in expression ratios (`RR__$__*`)
#'                observed across all subgenomes and group comparisons.
#'          \item `theta0__$`: Posterior expression ratio estimates shared
#'                across all conditions, where `$` denotes the subgenome name.
#'          \item `theta1__$__*`: Posterior expression ratio estimates specific
#'                to each condition,
#'                where `$` denotes the subgenome name and `*` denotes the group name.
#'      }
#'      All statistics are derived from MCMC samples and may differ from
#'      statistics calculated directly from RNA-seq read counts. An `NA` value
#'      indicates that the homeolog is expressed in only one condition, making
#'      ratio comparisons infeasible.
#'
#' @references Sun J, Sese J, Shimizu KK.
#'      A likelihood ratio test for detecting shifts in homeolog expression ratios in allopolyploids.
#'      New Phytol., 249(1):603-617, 2026.
#'      \doi{10.1111/nph.70648}
#'
#' @examples
#' x <- sim_homeolog_counts(10)
#' x_output <- hobit_mcmc(x)
#'
#' # fast mode
#' x_output <- hobit_mcmc(x, chains = 2, iter_warmup = 100, iter_sampling = 100)
#'
#' # strict mode
#' x_output <- hobit_mcmc(x, mode = "strict", chains = 2, iter_warmup = 100, iter_sampling = 100)
#'
#' # parallel processing
#' x_output <- hobit_mcmc(x, n_threads = 1, parallel_chains = 8,
#'                        iter_warmup = 100, iter_sampling = 100)
#' x_output <- hobit_mcmc(x, n_threads = 8, parallel_chains = 1,
#'                        iter_warmup = 100, iter_sampling = 100)
#'
#' @importFrom stats p.adjust
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_vapply
#' @importFrom progressr progressor with_progress handlers
#' @importFrom progress progress_bar
#' @importFrom cmdstanr cmdstan_model
#' @export
hobit_mcmc <- function(x,
    mode = c("fast", "strict"),
    use_Dirichlet = FALSE,
    no_replicate = FALSE,
    eps = 1e-3,
    dist = c("NB", "ZINB"),
    n_threads = getOption("mc.cores", 1L),
    parallel_chains = 1L,
    .debug = FALSE,
    ...
) {
    mode <- match.arg(mode)
    dist <- match.arg(dist)

    output <- if (mode == "fast") {
        .hobit_mcmc.fast(x = x,
            use_Dirichlet = use_Dirichlet, no_replicate = no_replicate, eps = eps, dist = dist,
            n_threads = n_threads, parallel_chains = parallel_chains, ...)
    } else {
        .hobit_mcmc.strict(x = x,
            use_Dirichlet = use_Dirichlet, no_replicate = no_replicate, eps = eps, dist = dist,
            n_threads = n_threads, parallel_chains = parallel_chains, ...)
    }

    .hobit.select_output(output, .debug)
}
