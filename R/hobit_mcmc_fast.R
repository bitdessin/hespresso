# Calculate prior parameters from observed count data.
# Use total counts (`sum`) to compute H1 prior parameters because simulation
# evaluations indicate that total counts produce stable and accurate results.
#' @importFrom stats median
.hb.mcmc.fast.init_params <- function(x, prior_func = 'sum') {
    # Gene expression.
    gene_exp_mu <- rep(0, length = x$N_CONDITIONS)
    replicate_id <- 1
    for (i in seq_len(x$N_CONDITIONS)) {
        for (j in seq_len(x$N_REPLICATES[i])) {
            gene_exp_mu[i] <-
                gene_exp_mu[i] + (x$GENE_EXP[replicate_id] / x$N_REPLICATES[i])
            replicate_id <- replicate_id + 1
        }
    }

    # Homeolog expression ratios.
    hexp_ratios <- .hb.mcmc.fast.init_params.hexp_ratios(x)
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


.hb.mcmc.fast.init_params.hexp_ratios <- function(x) {
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
        exp_sum <- rep(0, length = x$N_SUBGENOMES)
        for (j in seq_len(x$N_REPLICATES[i])) {
            exp_sum <- exp_sum + x$HOMEOLOG_EXP[replicate_id, ]
            replicate_id <- replicate_id + 1
        }
        if (sum(exp_sum) > 0)
            sum_p1[i, ] <- exp_sum / sum(exp_sum)
    }
    
    list(sum = sum_p1, median = med_p1, mean = avg_p1)
}


# Estimate dispersion for each group with edgeR.
# If `no_replicate = TRUE`, dispersion is estimated by treating all replicates
# as one group.
#' @importFrom edgeR DGEList estimateDisp
.hb.mcmc.fast.est_dispersion_mx <- function(x, group, group_names, no_replicate) {
    gene_disp <- matrix(0, nrow = nrow(x), ncol = length(group_names))
    if (no_replicate) {
        y <- edgeR::DGEList(counts = x)
        y <- suppressMessages(edgeR::estimateDisp(y))
        for (group_id in seq_along(group_names)) {
            gene_disp[, group_id] <- y$tagwise.dispersion
        }
    } else {
        for (group_id in seq_along(group_names)) {
            idx <- group == group_names[group_id]
            y <- edgeR::DGEList(counts = x[, idx, drop = FALSE])
            y <- suppressMessages(edgeR::estimateDisp(y))
            gene_disp[, group_id] <- y$tagwise.dispersion
        }
    }
    colnames(gene_disp) <- group_names
    1 / gene_disp
}


# Estimate tuple-wise and subgenome-wise dispersion for homeologs.
.hb.mcmc.fast.est_dispersion <- function(x, exp_group, group_names, no_replicate) {
    gene_exp <- 0
    for (i in seq_along(x@data)) {
        gene_exp <- gene_exp + x@data[[i]]
    }
    gene_disp <- .hb.mcmc.fast.est_dispersion_mx(gene_exp, exp_group, group_names, no_replicate)
    homeolog_disp <- vector("list", length(x@data))
    for (i in seq_along(x@data)) {
        homeolog_disp[[i]] <- .hb.mcmc.fast.est_dispersion_mx(x@data[[i]],
            exp_group, group_names, no_replicate)
    }
    list(gene = gene_disp, homeolog = homeolog_disp)
}


# Format count data for Stan fitting.
.hb.mcmc.fast.format_data <- function(x, use_Dirichlet, eps, no_replicate) {
    x_list <- vector('list', length = nrow(x@data[[1]]))

    exp_group <- droplevels(x@exp_design$group)
    if (!is.factor(exp_group)) {
        exp_group <- factor(exp_group, levels = unique(exp_group))
    }
    group_names <- levels(exp_group)
    n_replicates <- as.numeric(table(exp_group))

    disp <- .hb.mcmc.fast.est_dispersion(x, exp_group, group_names, no_replicate)

    for (i in seq_along(x_list)) {
        idx <- 1
        hexp <- matrix(NA_real_, nrow = sum(n_replicates), ncol = length(x@data))
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
        x_fmt$.META$init_params <- .hb.mcmc.fast.init_params(x_fmt)
        x_fmt$PRIOR_ALPHA0 <- x_fmt$.META$init_params$alpha0
        x_fmt$PRIOR_ALPHA1 <- x_fmt$.META$init_params$alpha1
        x_fmt$GENE_EXP_MU <- x_fmt$.META$init_params$gene_exp_mu
        
        x_list[[i]] <- x_fmt
    }
    x_list
}

#' @importFrom stats pchisq
.hb.mcmc.fast.calc_p <- function(log_lik, inputs, norm_log_lik = FALSE) {
    if (norm_log_lik) {
        n_replicates <- mean(inputs$N_REPLICATES)
        k <- (log(n_replicates) + 1) / n_replicates
        log_lik <- log_lik * k
    }
    
    df <- inputs$N_CONDITIONS * (inputs$N_SUBGENOMES - 1) - (inputs$N_SUBGENOMES - 1)
    lrt_lambda <- - 2 * (log_lik[, 1] - log_lik[, 2])
    median(pchisq(lrt_lambda, df = df, lower.tail = FALSE))
}


.hb.mcmc.fast.calc_shifts <- function(theta) {
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


# Parse Stan output into the common MCMC result format.
.hb.mcmc.fast.format_draws <- function(inputs, dmat) {
    # Theta for the null hypothesis.
    theta0 <- dmat[, base::match(paste0('theta0[',
                        seq(1, inputs$N_SUBGENOMES), ']'), colnames(dmat))]
    colnames(theta0) <- paste0('theta0__..__', inputs$.META$subgenome_names)
    
    # Theta for the alternative hypothesis.
    theta <- vector('list', length = inputs$N_CONDITIONS)
    for (i in seq_len(inputs$N_CONDITIONS)) {
        theta[[i]] <- dmat[, match(paste0('theta1[', i,',',
                        seq(1, inputs$N_SUBGENOMES), ']'), colnames(dmat))]
        colnames(theta[[i]]) <- paste0('theta1__..__',
                                       inputs$.META$subgenome_names, '__..__',
                                       inputs$.META$condition_names[i])
    }

    # Log-likelihood.
    log_lik <- dmat[, match(c('log_lik[1]', 'log_lik[2]'), colnames(dmat))]
    colnames(log_lik) <- c('logLik_H0', 'logLik_H1')
    lrt_draws <- -2 * (log_lik[, 'logLik_H0'] - log_lik[, 'logLik_H1'])
    df <- (inputs$N_CONDITIONS - 1L) * (inputs$N_SUBGENOMES - 1L)
    
    # Format output.
    v <- apply(theta0, 2, median)
    for (i in seq_along(theta)) {
        v <- c(v, apply(theta[[i]], 2, median))
    }
    v <- c(v, apply(log_lik, 2, median))
    v <- c(v,
           pvalue = .hb.mcmc.fast.calc_p(log_lik, inputs, TRUE),
           lrt = median(lrt_draws),
           df = df,
           lrt_negative_prob = mean(lrt_draws < 0),
           status_code = 0,
           .hb.mcmc.strict.calc_effects(theta, inputs),
           qvalue = NA_real_)
    v
}


.hb.mcmc.fast.fit_gene <- function(model, inputs, sample_args, stats_template) {
    if (any(inputs$GENE_EXP_MU == 0)) {
        out <- stats_template * NA_real_
        out['df'] <- (inputs$N_CONDITIONS - 1L) * (inputs$N_SUBGENOMES - 1L)
        out['status_code'] <- 1
        return(out)
    }

    sample_args <- sample_args
    sample_args$data <- inputs
    sample_args$data$.META <- NULL

    fit <- do.call(model$sample, sample_args)
    .hb.mcmc.fast.format_draws(inputs, fit$draws(inc_warmup = FALSE, format = 'matrix'))
}


# HOBIT MCMC implementation (fast mode).
#
# Fits the null and alternative models using MCMC and performs a likelihood ratio test.
# The null model likelihood is constructed from posterior samples of the full model.
# This implementation is reported and evaluated in Sun et al.,
# New Phytol., 249(1):603-617, 2026.
#
#' @importFrom stats p.adjust
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_vapply
#' @importFrom progressr progressor with_progress handlers
#' @importFrom progress progress_bar
.hobit_mcmc.fast <- function(x,
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

    x <- .norm_counts(x)
    
    input_params <- list(...)
    input_params$parallel_chains <- parallel_chains
    
    data <- .hb.mcmc.fast.format_data(x, use_Dirichlet, eps, no_replicate)
    
    stan_code_fpath <- system.file(package = 'hespresso', 'extdata',
                                   paste0('HOBIT.', dist, '.stan'))
    m <- cmdstanr::cmdstan_model(stan_file = stan_code_fpath,
                       dir = tempdir(), quiet = TRUE, compile = TRUE)
    
    # Calculate the number of outputs.
    .input_params <- input_params
    .input_params$data <- data[[1]]
    .input_params$data$.META <- NULL
    .input_params$refresh <- 0
    .input_params$show_messages <- FALSE
    .outputs <- do.call(m$sample, .input_params)
    .outputs_fmt <- .hb.mcmc.fast.format_draws(data[[1]],
                        .outputs$draws(inc_warmup = FALSE, format = 'matrix'))
    stats_template <- .outputs_fmt
    stats_names <- names(stats_template)

    future_plan_bkup <- future::plan()
    on.exit(future::plan(future_plan_bkup), add = TRUE)

    if (n_threads > 1L) {
        future::plan(future::multisession, workers = n_threads)
    } else {
        future::plan(future::sequential)
    }

    # HOBIT
    stats <- progressr::with_progress({
        pb <- progressr::progressor(length(data))
        future.apply::future_vapply(seq_along(data), function(i) {
            on.exit(pb(message = 'HOBIT (MCMC fast)'), add = TRUE)
            tryCatch(
                .hb.mcmc.fast.fit_gene(m, data[[i]], input_params, stats_template),
                error = function(e) {
                    warning(sprintf("HOBIT MCMC (fast) fitting failed for gene `%s`: %s",
                            x@gene_names[i], conditionMessage(e)),
                        call. = FALSE)
                    fit_error <- stats_template * NA_real_
                    fit_error['status_code'] <- 6
                    fit_error
                }
            )
        },
        FUN.VALUE = stats_template,
        future.seed = TRUE,
        future.packages = c('hespresso', 'cmdstanr'))
        },
        handlers = progressr::handler_progress(format = ":message :percent [:bar] :current/:total :eta")
    )


    rownames(stats) <- gsub("__\\.\\.__", "__", stats_names)
    stats <- data.frame(t(stats), check.names = FALSE)
    stats$qvalue <- p.adjust(stats$pvalue, method = 'BH')

    status_labels <- c(
        `0` = 'success', `1` = 'zero_expression', `6` = 'internal_error'
    )

    status <- unname(status_labels[as.character(stats$status_code)])
    data.frame(
        gene = x@gene_names,
        pvalue = stats$pvalue,
        qvalue = stats$qvalue,
        lrt = stats$lrt,
        df = stats$df,
        status = status,
        logLik_H0 = stats$logLik_H0,
        logLik_H1 = stats$logLik_H1,
        lrt_negative_prob = stats$lrt_negative_prob,
        stats[, grep('^D__', names(stats), value = TRUE), drop = FALSE],
        stats[, grep('^LR__', names(stats), value = TRUE), drop = FALSE],
        stats[, grep('^RR__', names(stats), value = TRUE), drop = FALSE],
        Dmax = stats$Dmax,
        LRmax = stats$LRmax,
        RRmax = stats$RRmax,
        stats[, grep('^theta0__', names(stats), value = TRUE), drop = FALSE],
        stats[, grep('^theta1__', names(stats), value = TRUE), drop = FALSE],
        check.names = FALSE)
}
