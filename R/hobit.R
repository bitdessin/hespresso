# Calculate Prior Parameters from Observed Count Data
# 
# @details
# Homeolog expression ratio (HER) is calculated using three methods:  
#  
# - **sum**: Calculates the total counts of all replicates, then computes HER.  
# - **mean**: Calculates HER for each replicate, then computes the mean.  
# - **median**: Calculates HER for each replicate, then computes the median.  
#  
# Simulation evaluations indicate that using the `sum` method produces stable
# and accurate results.
#' @importFrom stats median
.hb.init_params <- function(x, eps, prior_func = 'sum') {
    # gene expression
    gene_exp_mu <- rep(0, length = x$N_CONDITIONS)
    replicate_id <- 1
    for (i in seq_len(x$N_CONDITIONS)) {
        for (j in seq_len(x$N_REPLICATES[i])) {
            gene_exp_mu[i] <-
                gene_exp_mu[i] + (x$GENE_EXP[replicate_id] / x$N_REPLICATES[i])
            replicate_id <- replicate_id + 1
        }
    }
    gene_exp_mu[gene_exp_mu == 0] <- eps
    
    # homeolog expression ratios
    hexp_ratios <- .hb.init_params.hexp_ratios(x)
    alpha1 <- switch(prior_func,
        'sum' = sweep(hexp_ratios$sum, 1, x$N_SUBGENOMES, '*'),
        'median' = sweep(hexp_ratios$median, 1, x$N_SUBGENOMES, '*'),
        'mean' = sweep(hexp_ratios$mean, 1, x$N_SUBGENOMES, '*'),
        stop('The supported functions are `sum`, `median`, and `mean`.')
    )
    alpha1[alpha1 == 0] <- 0.01
    alpha0 <- colMeans(alpha1)
    
    list(gene_exp_mu = gene_exp_mu,
         alpha0 = alpha0, alpha1 = alpha1,
         ratios = hexp_ratios)
}


.hb.init_params.hexp_ratios <- function(x) {
    med_p1 <- avg_p1 <- matrix(1 / x$N_SUBGENOMES,
                               nrow = x$N_CONDITIONS, ncol = x$N_SUBGENOMES)
    raw_p1 <- x$HOMEOLOG_EXP / rowSums(x$HOMEOLOG_EXP)
    raw_p1[is.na(raw_p1)] <- 1 / x$N_SUBGENOMES
    replicate_id <- 1
    for (i in seq_len(x$N_CONDITIONS)) {
        raw_p1_subset <-
            raw_p1[seq(replicate_id, replicate_id + x$N_REPLICATES[i] - 1), ]
        if (is.null(nrow(raw_p1_subset))) {
            med_p1[i, ] <- avg_p1[i, ] <- raw_p1_subset
        } else if (nrow(raw_p1_subset) > 1) {
            med_p1[i, ] <- apply(raw_p1_subset, 2, median)
            avg_p1[i, ] <- apply(raw_p1_subset, 2, mean)
        }
        replicate_id <- replicate_id + x$N_REPLICATES[i]
    }
    
    sum_p1 <- matrix(1 / x$N_SUBGENOMES,
                     nrow = x$N_CONDITIONS, ncol = x$N_SUBGENOMES)
    replicate_id <- 1
    for (i in seq_len(x$N_CONDITIONS)) {
        exp_sum <- rep(0, length = x$N_SUBGENOME)
        for (j in seq_len(x$N_REPLICATES[i])) {
            exp_sum <- exp_sum + x$HOMEOLOG_EXP[replicate_id, ]
            replicate_id <- replicate_id + 1
        }
        if (sum(exp_sum) > 0)
            sum_p1[i, ] <- exp_sum / sum(exp_sum)
    }
    
    list(sum = sum_p1, median = med_p1, mean = avg_p1)
}


# Estimate Dispersion from Observed Count Data Using edgeR
# 
# @details
# Dispersion is estimated for each group
# if the data consists of multiple biological replicates.
# If `no_replicate` is set to `TRUE`, dispersion is estimated from  
# all counts without using group information. 
#
#' @importFrom edgeR DGEList estimateDisp
.hb.est_dispersion_mx <- function(x, group, no_replicate) {
    n_genes <- dim(x)
    group_names <- base::unique(group)
    
    gene_disp <- matrix(0, nrow = nrow(x), ncol = length(group_names))
    if (no_replicate) {
        y <- edgeR::DGEList(counts = x)
        y <- suppressMessages(edgeR::estimateDisp(y))
        for (group_id in seq_along(group_names)) {
            gene_disp[, group_id] <- y$tagwise.dispersion
        }
    } else {
        for (group_id in seq_along(group_names)) {
            y <- edgeR::DGEList(counts = x[, group == group_names[group_id]])
            y <- suppressMessages(edgeR::estimateDisp(y))
            gene_disp[, group_id] <- y$tagwise.dispersion
        }
    }
    1 / gene_disp
}


# Estimate Dispersion for Genes and Homeologs
# 
# @details
# This function estimates dispersions for genes and homeologs.  
# To estimate gene dispersions, the function sums all homeolog expressions.  
# It then calculates homeolog dispersions for each subgenome.
.hb.est_dispersion <- function(x, no_replicate) {
    gene_exp <- 0
    for (i in seq_along(x@data)) {
        gene_exp <- gene_exp + x@data[[i]]
    }
    gene_disp <- .hb.est_dispersion_mx(gene_exp, x@exp_design$group, no_replicate)
    
    homeolog_disp <- vector('list', length(x@data))
    for (i in seq_along(x@data)) {
        homeolog_disp[[i]] <- .hb.est_dispersion_mx(x@data[[i]],
                                    x@exp_design$group, no_replicate)
    }
    
    list(gene = gene_disp, homeolog = homeolog_disp)
}


# Convert Matrix Data to List
# 
# Converts a list of homeolog expression matrices into a list containing  
# homeolog expression data, prior parameters,
# and other information for modeling. 
.hb.format_data <- function(x, use_Dirichlet, eps, no_replicate) {
    x_list <- vector('list', length = nrow(x@data[[1]]))
    
    exp_group <- x@exp_design$group
    .groups <- base::table(exp_group)
    group_names <- names(.groups)
    n_replicates <- as.numeric(.groups)
    
    disp <- .hb.est_dispersion(x, no_replicate)
   
    for (i in seq_along(x_list)) {
        idx <- 1
        hexp <- matrix(NA, nrow = sum(n_replicates), ncol = length(x@data))
        for (g in seq_along(group_names)) {
            for (j in seq_along(x@data)) {
                hexp[seq(idx, (idx - 1) + n_replicates[g]), j] <-
                    x@data[[j]][i, exp_group == group_names[g]]
            }
            idx <- idx + n_replicates[g]
        }
        
        gexp_upper <- max(max(rowSums(hexp)) * 10, 100)
        gexp_lower <- max(rowSums(hexp)) * 0.1
        if (gexp_lower < 100) gexp_lower <- 0
        
        x_fmt <- list(
            N_SUBGENOMES = ncol(hexp),
            N_CONDITIONS = length(group_names),
            N_REPLICATES = n_replicates,
            GENE_EXP = rowSums(hexp),
            GENE_EXP_PHI = disp$gene[i, ],
            GENE_EXP_UPPER = gexp_upper,
            GENE_EXP_LOWER = gexp_lower,
            HOMEOLOG_EXP = hexp,
            HOMEOLOG_EXP_PHI = vapply(disp$homeolog, function(x, i) {
                x[i, ]}, numeric(ncol(disp$homeolog[[1]])), i),
            
            USE_DIRICHLET = ifelse(use_Dirichlet, 1, 0),
            EPS = eps,
            .META = list(condition_names = group_names,
                         subgenome_names = names(x@data),
                         init_params = NULL))
        x_fmt$.META$init_params <- .hb.init_params(x_fmt, eps)
        x_fmt$PRIOR_ALPHA0 <- x_fmt$.META$init_params$alpha0
        x_fmt$PRIOR_ALPHA1 <- x_fmt$.META$init_params$alpha1
        x_fmt$GENE_EXP_MU <- x_fmt$.META$init_params$gene_exp_mu
        
        x_list[[i]] <- x_fmt
    }
    x_list
}

#' @importFrom stats pchisq
.hb.calc_p <- function(log_lik, inputs, norm_log_lik = FALSE) {
    if (norm_log_lik) {
        n_replicates <- mean(inputs$N_REPLICATES)
        k <- (log(n_replicates) + 1) / n_replicates
        log_lik <- log_lik * k
    }
    
    df <- inputs$N_CONDITIONS * (inputs$N_SUBGENOMES - 1) - (inputs$N_SUBGENOMES - 1)
    lrt_lambda <- - 2 * (log_lik[, 1] - log_lik[, 2])
    median(pchisq(lrt_lambda, df = df, lower.tail = FALSE))
}


.hb.calc_shifts <- function(theta) {
    d <- or <- d_names <- or_names <- NULL
    for (s in seq_len(ncol(theta[[1]]))) {
        for (i in seq(1, length(theta) - 1)) {
            for (j in seq(i + 1, length(theta))) {
                i_tags <- strsplit(colnames(theta[[i]])[s], '__..__')[[1]]
                j_tags <- strsplit(colnames(theta[[j]])[s], '__..__')[[1]]
                d <- cbind(d, theta[[i]][, s] - theta[[j]][, s])
                or <- cbind(or, (theta[[i]][, s] / (1 - theta[[i]][, s])) / (theta[[j]][, s] / (1 - theta[[j]][, s])))
                d_names <- c(d_names,
                             paste0('D__', i_tags[2], '__(', i_tags[3], '-', j_tags[3], ')'))
                or_names <- c(or_names,
                              paste0('OR__', i_tags[2], '__(', i_tags[3], '/', j_tags[3], ')'))
            }
        }
        
    }
    d <- apply(d, 2, median)
    or <- apply(or, 2, median)
    d_stats <- c(d, or, max(abs(d), na.rm = TRUE), max(c(or, 1 / or), na.rm = TRUE))
    names(d_stats) <- c(d_names, or_names, 'Dmax', 'ORmax')
    d_stats
}


# Parse STAN Output  
#  
# Parses STAN output and converts it into a data frame containing
# posterior probabilities of homeolog expression ratios,
# log likelihoods of the full and reduced models, p-Value from LRT,
# Dmax, and ORmax.
.hb.format_draws <- function(inputs, dmat) {
    # theta for null hypothesis
    theta0 <- dmat[, base::match(paste0('theta0[',
                        seq(1, inputs$N_SUBGENOMES), ']'), colnames(dmat))]
    colnames(theta0) <- paste0('theta0__..__', inputs$.META$subgenome_names)
    
    # theta for alternative hypothesis
    theta <- vector('list', length = inputs$N_CONDITIONS)
    for (i in seq_len(inputs$N_CONDITIONS)) {
        theta[[i]] <- dmat[, match(paste0('theta1[', i,',',
                        seq(1, inputs$N_SUBGENOMES), ']'), colnames(dmat))]
        colnames(theta[[i]]) <- paste0('theta1__..__',
                                       inputs$.META$subgenome_names, '__..__',
                                       inputs$.META$condition_names[i])
    }

    # log-likelihood
    log_lik <- dmat[, match(c('log_lik[1]', 'log_lik[2]'), colnames(dmat))]
    colnames(log_lik) <- c('logLik_H0', 'logLik_H1')
    
    # format
    v <- apply(theta0, 2, median)
    for (i in seq_along(theta)) {
        v <- c(v, apply(theta[[i]], 2, median))
    }
    v <- c(v, apply(log_lik, 2, median))
    v <- c(v,
           .hb.calc_shifts(theta),
           raw_pvalue = .hb.calc_p(log_lik, inputs, FALSE),
           raw_qvalue = NA,
           pvalue = .hb.calc_p(log_lik, inputs, TRUE),
           qvalue = NA)
    
    
    #v_mat <- matrix(v, nrow = 1)
    #colnames(v_mat) <- names(v)
    #v_mat
    v
}


#' HOBIT: Detecting Shifts in Homeolog Expression Ratios
#' 
#' A statistical framework for identifying changes in homeolog expression ratios
#' across experimental conditions in allopolyploid species.
#' 
#' HOBIT is a statistical test for detecting homeologs with differential
#' homeolog expression ratios across multiple experimental conditions
#' in allopolyploid species, using RNA-Seq read count data.
#' It employs a likelihood ratio test (LRT) to compare two hierarchical models:
#' a full model that allows HERs to vary among conditions,
#' and a reduced model that assumes constant ratios across all conditions.
#' 
#' Expression counts are modeled using a negative binomial (NB) distribution.
#' Specifically, the expression level of a homeolog from the \eqn{i}-th subgenome,
#' denoted \eqn{x_{i}}, is assumed to follow an NB distribution with mean
#' \eqn{\mu \theta_{i}} and dispersion \eqn{\phi}, where \eqn{\mu} is the total
#' gene expression across subgenomes, \eqn{\theta_{i}} is the proportion of
#' expression attributable to the \eqn{i}-th subgenome, and \eqn{\phi} is the
#' dispersion parameter:
#' 
#' \deqn{
#' x_{i} \sim NB(\mu \theta_{i}, \phi)
#' }
#' 
#' By default, \eqn{\mu} and \eqn{\boldsymbol{\theta}} are sampled from
#' uninformative priors, while dispersion \eqn{\phi} is estimated using
#' \code{\link[edgeR]{estimateDisp}} from the \pkg{edgeR} package.
#' 
#' HOBIT uses Markov chain Monte Carlo (MCMC) sampling to fit both the full and
#' reduced models and compute likelihoods for the LRT. The main entry point is
#' the `hobit()` function, which performs parameter estimation, model fitting,
#' likelihood computation, and hypothesis testing.
#'
#' @param x An \linkS4class{ExpMX} object containing normalized homeolog
#'      expression data (RNA-Seq read counts).
#' @param use_Dirichlet Logical. Whether to apply a Dirichlet prior distribution
#'      when sampling expression ratios.
#' @param no_replicate Logical. If `TRUE`, all replicates are treated as a
#'      single condition when estimating dispersion. This avoids errors in cases
#'      where no biological replicates are available.
#' @param eps Numeric. Minimum threshold for homeolog expression. Values below
#'      this threshold are replaced with 0 to avoid instability during MCMC sampling,
#'      especially when fitting NB models with near-zero means.
#' @param dist Character. Distribution used to model expression data.
#'      Options: `'NB'` (negative binomial) or `'ZINB'` (zero-inflated negative
#'      binomial).
#' @param n_threads Integer. Number of threads for parallelizing computation
#'      across homeologs. This parallelization is independent of MCMC sampling,
#'      which is controlled by `parallel_chains`. Defaults to the global option
#'      `'mc.cores'`, or 1 if not set. In practice, using a larger `n_threads`
#'      with `parallel_chains = 1` typically reduces overall execution time 
#'      compared to allocating more threads to MCMC chains.
#' @param parallel_chains Integer. Number of parallel MCMC chains. Passed
#'      directly to \code{\link[cmdstanr]{sample}} in the \pkg{cmdstanr} package.
#'      If both `n_threads` and `parallel_chains` exceed 1, the total number of
#'      threads used is `n_threads × parallel_chains`. Ensure adequate resources
#'      are available before increasing both parameters.
#' @param ... Additional arguments passed to
#'      \code{\link[cmdstanr]{sample}} (e.g., `chains`, `iter_warmup`,
#'      `iter_sampling`, `thin`).
#'
#' @return A \code{data.frame} with one row per homeolog, containing:
#'      \itemize{
#'          \item `pvalue`: p-value from the LRT using normalized likelihoods.
#'          \item `qvalue`: Benjamini–Hochberg adjusted `pvalue`.
#'          \item `raw_pvalue`: p-value from the LRT using raw likelihoods.
#'          \item `raw_qvalue`: Benjamini–Hochberg adjusted `raw_pvalue`.
#'          \item `D__$__*`: Difference in expression ratios for subgenome `$`
#'                between the two groups represented by `*`.
#'          \item `OR__$__*`: Odds ratio of expression ratios for subgenome `$`
#'                between the two groups represented by `*`.
#'          \item `Dmax`: Maximum absolute difference in expression ratios (`D__$__*`)
#'                observed across all subgenomes and group comparisons.
#'          \item `ORmax`: Maximum odds ratio in expression ratios (`OR__$__*`)
#'                observed across all subgenomes and group comparisons.
#'          \item `theta0__$`: Posterior expression ratio estimates shared
#'                across all conditions, where `$` denotes the subgenome name.
#'          \item `theta1__$__*`: Posterior expression ratio estimates specific
#'                to each condition,
#'                where `$` denotes the subgenome name and `*` denotes the group name.
#'          \item `logLik_H0`: Log-likelihood of the reduced model.
#'          \item `logLik_H1`: Log-likelihood of the full model.
#'      }
#'      All statistics are derived from MCMC samples
#'      and may differ from those calculated directly from the RNA-Seq read counts.
#'
#' @references Sun J, Sese J, and Shimizu KK.
#'      A moderated statistical test for detecting shifts in homeolog expression
#'      ratios in allopolyploids.
#'      bioRxiv 2025;2025.07.01.660977. \doi{10.1101/2025.07.01.660977}
#'      
#' @examples
#' x <- sim_homeolog_counts(10)
#' x_output <- hobit(x)
#' 
#' # MCMC sampling options
#' x_output <- hobit(x, chains = 2, iter_warmup = 100, iter_sampling = 100)
#' 
#' # parallel processing
#' x_output <- hobit(x, n_threads = 1, parallel_chains = 8, iter_warmup = 100, iter_sampling = 100)
#' x_output <- hobit(x, n_threads = 8, parallel_chains = 1, iter_warmup = 100, iter_sampling = 100)
#' 
#' @importFrom stats p.adjust
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_vapply
#' @importFrom progressr progressor with_progress handlers
#' @importFrom progress progress_bar
#' @importFrom cmdstanr cmdstan_model 
#' @export
hobit <- function(x,
                  use_Dirichlet = FALSE,
                  no_replicate = FALSE,
                  eps = 1e-3,
                  dist = c('NB', 'ZINB'),
                  n_threads = getOption('mc.cores', 1),
                  parallel_chains = 1,
                  ...) {
    dist <- match.arg(dist)
    if (any(grepl("__\\.\\.__", x@exp_design$group))) {
        warning('One or more group names contain the reserved string "__..__".',
                'Please rename the group(s) and try again.')
        user_ans <- readline("Would you like to rename them now? (y/n): ")
        if (tolower(user_ans) == "y") {
            message('Please rename the group name(s) and try again.')
            return(NULL)
        } else if (tolower(user_ans) == "n") {
            message('Continuing the program despite reserved group names.')
        } else {
            message('Invalid input. Continuing the program by default.')
        }
    }
    
    input_params <- list(...)
    input_params$parallel_chains <- parallel_chains
    
    data <- .hb.format_data(x, use_Dirichlet, eps, no_replicate)
    
    stan_code_fpath <- system.file(package = 'hespresso', 'extdata',
                                   paste0('HOBIT.', dist, '.stan'))
    m <- cmdstan_model(stan_file = stan_code_fpath,
                       dir = tempdir(), quiet = TRUE, compile = TRUE)
    
    # calculate n of outputs
    .input_params <- input_params
    .input_params$data <- data[[1]]
    .input_params$data$.META <- NULL
    .outputs <- do.call(m$sample, .input_params)
    .outputs_fmt <- .hb.format_draws(data[[1]],
                        .outputs$draws(inc_warmup = FALSE, format = 'matrix'))
    n_stats <- length(.outputs_fmt)
    
    # HOBIT
    plan(multisession, workers = n_threads)
    handlers("progress")
    if (FALSE) progress_bar
    with_progress({
        pb <- progressor(length(data))
        stats <- future_vapply(seq_along(data), function(i) {
            input_params$data <- data[[i]]
            input_params$data$.META <- NULL
            outputs <- do.call(m$sample, input_params)
            outputs_fmt <- .hb.format_draws(data[[i]],
                                    outputs$draws(inc_warmup = FALSE, format = 'matrix'))
            pb()
            outputs_fmt
        },
        FUN.VALUE = numeric(n_stats),
        future.seed = TRUE,
        future.packages = c('cmdstanr'))
    })
    plan(sequential)

    stats_names <- rownames(stats)
    rownames(stats) <- gsub("__\\.\\.__", "__", stats_names)
    stats <- data.frame(t(stats), check.names = FALSE)
    stats$qvalue <- p.adjust(stats$pvalue, method = 'BH')
    stats$raw_qvalue <- p.adjust(stats$raw_pvalue, method = 'BH')
    data.frame(gene = x@gene_names,
               stats[, c('pvalue', 'qvalue', 'raw_pvalue', 'raw_qvalue')],
               stats[, grep('^D__', colnames(stats))],
               stats[, grep('^OR__', colnames(stats))],
               stats[, c('Dmax', 'ORmax')],
               stats[, grep('^theta0__', colnames(stats))],
               stats[, grep('^theta1__', colnames(stats))],
               stats[, c('logLik_H0', 'logLik_H1')], 
               check.names = FALSE)
}
