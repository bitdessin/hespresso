.hb.mcmc.strict.make_prior <- function(init_params, n_subgenomes, concentration = n_subgenomes, min_alpha = 0.01) {
    alpha0 <- init_params$theta0 * concentration
    alpha1 <- init_params$theta1 * concentration
    alpha0 <- pmax(alpha0, min_alpha)
    alpha1 <- pmax(alpha1, min_alpha)
    alpha0 <- alpha0 / rowSums(alpha0) * concentration
    alpha1 <- alpha1 / rowSums(alpha1) * concentration
    list(H0 = unname(alpha0), H1 = unname(alpha1))
}

.hb.mcmc.strict.format_data <- function(x, use_Dirichlet, eps, no_replicate) {
    gdata <- .parse_group_names(x@exp_design$group)
    nf <- .hb.mle.calc_nf(x)
    disp <- .hb.mle.est_dispersion(x, gdata, nf, no_replicate)
    log_offset <- nf$log_offset[gdata$sample_order]

    if (length(log_offset) != sum(gdata$n_replicates) || any(!is.finite(log_offset))) {
        stop("Invalid normalization offsets.", call. = FALSE)
    }

    n_genes <- nrow(x@data[[1L]])
    n_subgenomes <- length(x@data)
    n_conditions <- gdata$n_conditions
    x_list <- vector(mode = "list", length = n_genes)

    for (i in seq_len(n_genes)) {
        homeolog_exp <- do.call(cbind, lapply(x@data, function(m) { as.integer(m[i, ]) }))
        homeolog_exp <- homeolog_exp[gdata$sample_order, , drop = FALSE]
        colnames(homeolog_exp) <- names(x@data)

        homeolog_dispersion <- vapply(disp$subgenome, function(d) { d[i, ] }, numeric(n_conditions))

        dimnames(homeolog_dispersion) <- list(gdata$names, names(x@data))
        init_params <- .hb.mle.init_params(homeolog_exp, gdata, log_offset)
        prior <- .hb.mcmc.strict.make_prior(init_params, n_subgenomes)

        x_list[[i]] <- list(
            USE_DIRICHLET = as.integer(use_Dirichlet),
            EPS = as.numeric(eps),
            N_SUBGENOMES = as.integer(n_subgenomes),
            N_CONDITIONS = as.integer(n_conditions),
            N_REPLICATES = as.integer(gdata$n_replicates),
            HOMEOLOG_EXP = homeolog_exp,
            HOMEOLOG_EXP_DISPERSION = homeolog_dispersion,
            LOG_OFFSET = as.numeric(log_offset),
            PRIOR_ALPHA0 = prior$H0,
            PRIOR_ALPHA1 = prior$H1,
            .META = list(
                condition_names = gdata$names,
                subgenome_names = names(x@data),
                init_params = init_params,
                valid_disp = all(is.finite(homeolog_dispersion) & homeolog_dispersion > 0)
            )
        )
    }
    x_list
}

.hb.mcmc.strict.as_stan_data <- function(inputs, hypothesis = c("H0", "H1")) {
    hypothesis <- match.arg(hypothesis)

    if (hypothesis == "H0") {
        n_theta <- 1L
        theta_id <- rep.int(1L, inputs$N_CONDITIONS)
        prior_alpha <- inputs$PRIOR_ALPHA0
    } else {
        n_theta <- inputs$N_CONDITIONS
        theta_id <- seq_len(inputs$N_CONDITIONS)
        prior_alpha <- inputs$PRIOR_ALPHA1
    }

    list(
        USE_DIRICHLET = as.integer(inputs$USE_DIRICHLET),
        EPS = as.numeric(inputs$EPS),
        N_SUBGENOMES = as.integer(inputs$N_SUBGENOMES),
        N_CONDITIONS = as.integer(inputs$N_CONDITIONS),
        N_REPLICATES = as.integer(inputs$N_REPLICATES),
        N_THETA = as.integer(n_theta),
        THETA_ID = as.integer(theta_id),
        HOMEOLOG_EXP = inputs$HOMEOLOG_EXP,
        HOMEOLOG_EXP_DISPERSION = unname(inputs$HOMEOLOG_EXP_DISPERSION),
        LOG_OFFSET = as.numeric(inputs$LOG_OFFSET),
        PRIOR_ALPHA = unname(prior_alpha)
    )
}

.hb.mcmc.strict.make_init <- function(inputs, hypothesis = c("H0", "H1")) {
    hypothesis <- match.arg(hypothesis)
    init_params <- inputs$.META$init_params
    gene_exp_mu <- pmax(as.numeric(init_params$gene_exp_mu), inputs$EPS)

    theta <- if (hypothesis == "H0") {
        init_params$theta0
    } else {
        init_params$theta1
    }

    theta <- pmax(as.matrix(theta), 1e-10)
    theta <- theta / rowSums(theta)
    list(log_gene_exp_mu = log(gene_exp_mu), theta = unname(theta))
}

.hb.mcmc.strict.make_init_function <- function(inputs, hypothesis, jitter_sd = 0.02) {
    base_init <- .hb.mcmc.strict.make_init(inputs, hypothesis)
    function(chain_id) {
        theta <- base_init$theta
        theta <- theta * exp(matrix(stats::rnorm(length(theta), sd = jitter_sd),
                                    nrow = nrow(theta), ncol = ncol(theta)))
        theta <- theta / rowSums(theta)
        list(log_gene_exp_mu = base_init$log_gene_exp_mu + stats::rnorm(length(base_init$log_gene_exp_mu), sd = jitter_sd),
             theta = unname(theta))
    }
}


.hb.mcmc.strict.seed <- function(base_seed, gene_id, hypothesis = c("H0", "H1")) {
    hypothesis <- match.arg(hypothesis)
    model_offset <- ifelse(hypothesis == "H0", 0, 1)
    if (is.null(base_seed)) {
        return(sample.int(2147483646L, size = 1L))
    }
    raw_seed <- as.double(base_seed) + 2 * (gene_id - 1L) + model_offset
    as.integer(((raw_seed - 1) %% 2147483646) + 1)
}


.hb.mcmc.strict.sample_one <- function(model, inputs, hypothesis, sample_args, seed) {
    args <- sample_args
    args$data <- .hb.mcmc.strict.as_stan_data(inputs, hypothesis)
    init_values <- .hb.mcmc.strict.make_init(inputs, hypothesis)
    args$init <- .hb.mcmc.strict.make_init_function(inputs, hypothesis)
    args$seed <- seed
    do.call(model$sample, args)
}

.hb.mcmc.strict.extract_draws <- function(fit, n_theta, n_subgenomes) {
    draws <- fit$draws(
        variables = c("theta", "gene_exp_mu", "log_lik_total"),
        inc_warmup = FALSE,
        format = "matrix"
    )
    if (nrow(draws) < 1L) {
        stop("Stan returned no post-warmup draws.", call. = FALSE)
    }
    theta <- lapply(seq_len(n_theta), function(k) {
        theta_names <- paste0("theta[", k, ",", seq_len(n_subgenomes), "]")
        if (!all(theta_names %in% colnames(draws))) {
            stop("Required `theta` draws were not returned.", call. = FALSE)
        }
        as.matrix(draws[, theta_names, drop = FALSE])
    })
    if (!("log_lik_total" %in% colnames(draws))) {
        stop("`log_lik_total` was not returned by Stan.", call. = FALSE)
    }
    log_lik <- as.numeric(draws[, "log_lik_total"])
    if (any(!is.finite(log_lik))) {
        stop("Stan returned non-finite log-likelihood values.", call. = FALSE)
    }
    list(theta = theta, log_lik = log_lik)
}

.hb.mcmc.strict.calc_effects <- function(theta_draws, inputs, eps = 1e-12) {
    condition_names <- inputs$.META$condition_names
    subgenome_names <- inputs$.META$subgenome_names

    theta_draws <- lapply(theta_draws, function(theta) {
        theta <- pmax(theta, eps)
        theta / rowSums(theta)
    })

    d <- numeric()
    d_names <- character()

    for (s in seq_len(inputs$N_SUBGENOMES)) {
        for (i in seq_len(inputs$N_CONDITIONS - 1L)) {
            for (j in seq.int(i + 1L, inputs$N_CONDITIONS)) {
                effect <- theta_draws[[i]][, s] - theta_draws[[j]][, s]
                d <- c(d, stats::median(effect))
                d_names <- c(d_names,
                    paste0("D__", subgenome_names[s], "__(", condition_names[i], "-", condition_names[j], ")"))
            }
        }
    }

    lr <- numeric()
    rr <- numeric()
    lr_names <- character()
    rr_names <- character()

    subgenome_pairs <- utils::combn(seq_len(inputs$N_SUBGENOMES), 2L)

    for (pair_id in seq_len(ncol(subgenome_pairs))) {
        a <- subgenome_pairs[1L, pair_id]
        b <- subgenome_pairs[2L, pair_id]
        pair_name <- paste0(subgenome_names[a], "/", subgenome_names[b])

        for (i in seq_len(inputs$N_CONDITIONS - 1L)) {
            for (j in seq.int(i + 1L, inputs$N_CONDITIONS)) {
                lr_draw <- log(theta_draws[[i]][, a] / theta_draws[[i]][, b]) - log(theta_draws[[j]][, a] / theta_draws[[j]][, b])
                lr_median <- stats::median(lr_draw)
                lr <- c(lr, lr_median)
                rr <- c(rr, exp(lr_median))
                lr_names <- c(lr_names,
                    paste0("LR__", pair_name, "__(", condition_names[i], "-", condition_names[j], ")"))
                rr_names <- c(rr_names,
                    paste0("RR__", pair_name, "__(", condition_names[i], "/", condition_names[j], ")"))
            }
        }
    }

    c(stats::setNames(d, d_names), stats::setNames(lr, lr_names), stats::setNames(rr, rr_names),
      Dmax = max(abs(d)), LRmax = max(abs(lr)), RRmax = exp(max(abs(lr))))
}


.hb.mcmc.strict.init_stats_tmpl <- function(inputs) {
    dummy_theta <- lapply(seq_len(inputs$N_CONDITIONS), function(c) {
        matrix(inputs$.META$init_params$theta1[c, ], nrow = 1L)
    })

    effect_stats <- .hb.mcmc.strict.calc_effects(dummy_theta, inputs)
    theta0_names <- paste0("theta0__", inputs$.META$subgenome_names)

    theta1_names <- unlist(lapply(seq_len(inputs$N_CONDITIONS), function(c) {
            paste0("theta1__", inputs$.META$subgenome_names, "__", inputs$.META$condition_names[c])
        }), use.names = FALSE)

    stats_names <- c(
        "pvalue", "lrt", "df", "logLik_H0", "logLik_H1",
        "lrt_negative_prob", "status_code", names(effect_stats), theta0_names, theta1_names
    )

    stats::setNames(rep(NA_real_, length(stats_names)), stats_names)
}

.hb.mcmc.strict.format_pair <- function(inputs, fit_h0, fit_h1) {
    draws_h0 <- .hb.mcmc.strict.extract_draws(
        fit_h0, n_theta = 1L, n_subgenomes = inputs$N_SUBGENOMES
    )
    draws_h1 <- .hb.mcmc.strict.extract_draws(
        fit_h1, n_theta = inputs$N_CONDITIONS, n_subgenomes = inputs$N_SUBGENOMES
    )
    n_draws <- min(length(draws_h0$log_lik), length(draws_h1$log_lik))
    if (n_draws < 1L) {
        stop("No paired H0/H1 draws are available.", call. = FALSE)
    }

    log_lik_h0 <- draws_h0$log_lik[seq_len(n_draws)]
    log_lik_h1 <- draws_h1$log_lik[seq_len(n_draws)]
    lrt_draws <- 2 * (log_lik_h1 - log_lik_h0)
    finite <- is.finite(lrt_draws)

    if (!any(finite)) {
        stop("All likelihood-ratio draws are non-finite.", call. = FALSE)
    }

    lrt_draws <- lrt_draws[finite]
    df <- (inputs$N_CONDITIONS - 1L) * (inputs$N_SUBGENOMES - 1L)

    pvalue_draws <- stats::pchisq(pmax(lrt_draws, 0), df = df, lower.tail = FALSE)
    theta0 <- apply(draws_h0$theta[[1L]], 2L, stats::median)

    theta1 <- unlist(
        lapply(draws_h1$theta, function(theta) { apply(theta, 2L, stats::median) }),
        use.names = FALSE
    )

    effects <- .hb.mcmc.strict.calc_effects(draws_h1$theta, inputs)
    out <- .hb.mcmc.strict.init_stats_tmpl(inputs)
    out["pvalue"] <- stats::median(pvalue_draws)
    out["lrt"] <- stats::median(lrt_draws)
    out["df"] <- df
    out["logLik_H0"] <- stats::median(log_lik_h0)
    out["logLik_H1"] <- stats::median(log_lik_h1)
    out["lrt_negative_prob"] <- mean(lrt_draws < 0)
    out["status_code"] <- 0
    out[names(effects)] <- effects

    theta0_names <- paste0("theta0__", inputs$.META$subgenome_names)
    out[theta0_names] <- theta0

    theta1_names <- unlist(lapply(seq_len(inputs$N_CONDITIONS), function(c) {
            paste0("theta1__", inputs$.META$subgenome_names, "__", inputs$.META$condition_names[c])
        }), use.names = FALSE)

    out[theta1_names] <- theta1
    out
}

.hb.mcmc.strict.fit_gene <- function(model, inputs, sample_args, base_seed, gene_id) {
    out <- .hb.mcmc.strict.init_stats_tmpl(inputs)
    out["df"] <- (inputs$N_CONDITIONS - 1L) * (inputs$N_SUBGENOMES - 1L)

    if (isTRUE(inputs$.META$init_params$zero_condition)) {
        out["status_code"] <- 1
        return(out)
    }

    if (!isTRUE(inputs$.META$valid_disp)) {
        out["status_code"] <- 2
        return(out)
    }

    seed_h0 <- .hb.mcmc.strict.seed(base_seed, gene_id, "H0")
    seed_h1 <- .hb.mcmc.strict.seed(base_seed, gene_id, "H1")

    fit_h0 <- tryCatch(
        .hb.mcmc.strict.sample_one(model, inputs, "H0", sample_args, seed_h0),
        error = function(e) e
    )

    if (inherits(fit_h0, "error")) {
        out["status_code"] <- 3
        return(out)
    }

    fit_h1 <- tryCatch(
        .hb.mcmc.strict.sample_one(model, inputs, "H1", sample_args, seed_h1),
        error = function(e) e
    )

    if (inherits(fit_h1, "error")) {
        out["status_code"] <- 4
        return(out)
    }

    result <- tryCatch(
        .hb.mcmc.strict.format_pair(inputs, fit_h0, fit_h1),
        error = function(e) e
    )

    if (inherits(result, "error")) {
        out["status_code"] <- 5
        return(out)
    }

    result
}

# HOBIT MCMC implementation (strict mode).
#
# Fits the null and alternative HOBIT models in separate Stan sampling
# runs using raw integer counts, TMM normalization offsets, and fixed
# edgeR dispersions.
#
#' @importFrom stats p.adjust
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_vapply
#' @importFrom progressr progressor with_progress handlers
#' @importFrom progress progress_bar
.hobit_mcmc.strict <- function(x, use_Dirichlet = FALSE, no_replicate = FALSE, eps = 1e-3, dist = 'NB', n_threads = getOption("mc.cores", 1L), parallel_chains = 1L, ...) {
    dist <- match.arg(dist)
    eps <- .as_positive_float(eps, "eps")
    n_threads <- .as_positive_int(n_threads, "n_threads")
    parallel_chains <- .as_positive_int(parallel_chains, "parallel_chains")

    sample_args <- list(...)
    reserved_args <- intersect(names(sample_args), c("data", "init", "parallel_chains"))

    if (length(reserved_args) > 0L) {
        stop(paste0(
            "The following arguments are controlled internally: ",
            paste(reserved_args, collapse = ", "), "."
        ), call. = FALSE)
    }

    base_seed <- sample_args$seed
    sample_args$seed <- NULL

    if (!is.null(base_seed)) {
        base_seed <- .as_positive_int(base_seed, "seed")
    }
    sample_args$parallel_chains <- parallel_chains

    data <- .hb.mcmc.strict.format_data(x, use_Dirichlet, eps, no_replicate)

    stan_file <- system.file("extdata", "HOBIT.NB.strict.stan", package = "hespresso")
    model <- cmdstanr::cmdstan_model(
        stan_file = stan_file, dir = tempdir(), quiet = TRUE, compile = TRUE
    )

    stats_template <- .hb.mcmc.strict.init_stats_tmpl(data[[1L]])

    future_plan_bkup <- future::plan()
    on.exit(future::plan(future_plan_bkup), add = TRUE)

    if (n_threads > 1L) {
        future::plan(future::multisession, workers = n_threads)
    } else {
        future::plan(future::sequential)
    }

    stats <- progressr::with_progress({
        pb <- progressr::progressor(length(data))
        future.apply::future_vapply(seq_along(data), function(i) {
                on.exit(pb(message = 'HOBIT (MCMC strict)'), add = TRUE)
                tryCatch(
                    .hb.mcmc.strict.fit_gene(
                        model = model, inputs = data[[i]], sample_args = sample_args,
                        base_seed = base_seed, gene_id = i
                    ),
                    error = function(e) {
                        warning(sprintf("HOBIT MCMC (strict) fitting failed for gene `%s`: %s", x@gene_names[i], conditionMessage(e)),
                            call. = FALSE)
                        fit_error <- stats_template
                        fit_error["status_code"] <- 6
                        fit_error
                    }
                )
            },
            FUN.VALUE = stats_template,
            future.seed = TRUE,
            future.packages = c("hespresso", "cmdstanr"))
        },
        handlers = progressr::handler_progress(format = ":message :percent [:bar] :current/:total :eta")
    )

    stats <- data.frame(t(stats), check.names = FALSE)
    stats$qvalue <- stats::p.adjust(stats$pvalue, method = "BH")

    status_labels <- c(
        `0` = "success", `1` = "zero_expression", `2` = "dispersion_failed",
        `3` = "h0_sampling_failed", `4` = "h1_sampling_failed",
        `5` = "result_extraction_failed", `6` = "internal_error"
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
        stats[, grep("^D__", names(stats), value = TRUE), drop = FALSE],
        stats[, grep("^LR__", names(stats), value = TRUE), drop = FALSE],
        stats[, grep("^RR__", names(stats), value = TRUE), drop = FALSE],
        Dmax = stats$Dmax,
        LRmax = stats$LRmax,
        RRmax = stats$RRmax,
        stats[, grep("^theta0__", names(stats), value = TRUE), drop = FALSE],
        stats[, grep("^theta1__", names(stats), value = TRUE), drop = FALSE],
        check.names = FALSE)
}
