# Load the separately compiled TMB model library.
#' @useDynLib hespresso, .registration = TRUE
NULL


.parse_group_names <- function(group) {
    if (!is.factor(group)) {
        group <- factor(group, levels = unique(group))
        #stop("`group` must be a factor.", call. = FALSE)
    }
    if (anyNA(group)) {
        stop("`group` must not contain NA.", call. = FALSE)
    }

    group <- droplevels(group)
    group_names <- levels(group)

    if (length(group_names) < 2L) {
        stop("At least two conditions are required.", call. = FALSE)
    }

    sample_order <- unlist(
        lapply(group_names, function(g) which(group == g)),
        use.names = FALSE)

    ordered_group <- factor(
        group[sample_order],
        levels = group_names
    )

    list(
        group = group,
        ordered_group = ordered_group,
        sample_order = sample_order,
        names = group_names,
        n_replicates = as.integer(table(group)),
        n_conditions = length(group_names)
    )
}

.hb.mle.calc_nf <- function(x) {
    h_exp <- .calc_hexp(x)
    y <- edgeR::DGEList(counts = h_exp)
    y <- edgeR::calcNormFactors(y)
    effective_lib_size <- y$samples$lib.size * y$samples$norm.factors
    log_offset <- log(effective_lib_size)
    log_offset <- log_offset - mean(log_offset)
    list(nf = y$samples$norm.factors,
         libsize = y$samples$lib.size,
         effective_libsize = effective_lib_size,
         log_offset = log_offset)
}


.hb.mle.init_nlminb_control <- function(dots) {
    if (length(dots) > 0L && (is.null(names(dots)) || any(names(dots) == ""))) {
        stop("All arguments passed through `...` must be named.", call. = FALSE)
    }
    defaults <- list(
        eval.max = 2000L,
        iter.max = 2000L,
        rel.tol = 1e-10,
        x.tol = 1e-8,
        trace = 0L
    )
    utils::modifyList(defaults, dots)
}


# Estimate dispersion from a matrix using edgeR.
#
#' @importFrom edgeR DGEList estimateDisp
.hb.mle.est_dispersion_from_matrix <- function(x, gdata, nf, no_replicate) {
    if (no_replicate) {
        y <- edgeR::DGEList(counts = x, group = rep(1, length(gdata$group)),
                            lib.size = nf$libsize, norm.factors = nf$nf)
        y <- suppressMessages(edgeR::estimateDisp(y))
    } else {
        y <- edgeR::DGEList(counts = x, group = gdata$group, lib.size = nf$libsize, norm.factors = nf$nf)
        design <- stats::model.matrix(~ 0 + group, data = data.frame(group = gdata$group))
        y <- suppressMessages(edgeR::estimateDisp(y, design = design, robust = TRUE))
    }

    y_disp <- y$tagwise.dispersion
    if (is.null(y_disp) || length(y_disp) != nrow(x)) {
        stop("edgeR failed to return tagwise dispersions.", call. = FALSE)
    }
    disp <- matrix(rep(y_disp, times = gdata$n_conditions), nrow = nrow(x), ncol = gdata$n_conditions,
                   dimnames = list(NULL, gdata$names))
    disp
}




# Estimate tuple-wise and subgenome-wise dispersion from an SeqCountData object.
#
# Tuple-wise dispersion is estimated from counts summed across all subgenomes
# for each homeolog tuple. In contrast, subgenome-wise dispersion is estimated
# independently for each subgenome.
.hb.mle.est_dispersion <- function(x, gdata, nf, no_replicate) {
    # tuple-wise
    h_disp <- .hb.mle.est_dispersion_from_matrix(.calc_hexp(x), gdata, nf, no_replicate)

    # subgenome-wise
    s_disp <- lapply(x@data, function(s_exp) {
        .hb.mle.est_dispersion_from_matrix(s_exp, gdata, nf, no_replicate)
    })
    names(s_disp) <- names(x@data)

    list(homeolog = h_disp, subgenome = s_disp)
}


.hb.mle.norm_ratio <- function(x, esp = 1e-8) {
    x <- as.numeric(x)
    n_subgenomes <- length(x)
    total_counts <- sum(x)
    if (!is.finite(total_counts) || total_counts <= 0) {
        return(rep(1 / n_subgenomes, n_subgenomes))
    }
    ratio <- x / total_counts
    ratio <- pmax(ratio, esp)
    ratio / sum(ratio)
}


.hb.mle.init_params <- function(x, gdata, log_offset) {
    gene_exp_mu <- numeric(gdata$n_conditions)
    names(gene_exp_mu) <- gdata$names

    theta1 <- matrix(
        1 / ncol(x),
        nrow = gdata$n_conditions,
        ncol = ncol(x)
    )
    rownames(theta1) <- gdata$names

    relative_library_size <- exp(log_offset)
    adjusted_x <- sweep(
        x,
        MARGIN = 1L,
        STATS = relative_library_size,
        FUN = "/"
    )

    zero_condition <- FALSE

    for (g_name in gdata$names) {
        idx <- gdata$ordered_group == g_name
        v <- x[idx, , drop = FALSE]
        adjusted_v <- adjusted_x[idx, , drop = FALSE]

        if (sum(v) == 0L) {
            zero_condition <- TRUE
            next
        }

        gene_exp_mu[g_name] <-
            mean(rowSums(adjusted_v))

        theta1[g_name, ] <-
            .hb.mle.norm_ratio(colSums(adjusted_v))
    }

    theta0 <- matrix(
        .hb.mle.norm_ratio(colSums(adjusted_x)),
        nrow = 1L
    )

    list(
        gene_exp_mu = gene_exp_mu,
        theta0 = theta0,
        theta1 = theta1,
        zero_condition = zero_condition
    )
}


.hb.mle.format_data <- function(x, eps, no_replicate) {
    gdata <- .parse_group_names(x@exp_design$group)
    nf <- .hb.mle.calc_nf(x)
    nf_log_offset <- nf$log_offset[gdata$sample_order]
    disp <- .hb.mle.est_dispersion(x, gdata, nf, no_replicate)

    n_genes <- nrow(x@data[[1L]])
    n_subgenomes <- length(x@data)
    n_conditions <- length(gdata$names)

    x_list <- vector("list", n_genes)
    for (i in seq_len(n_genes)) {
        s_expr <- do.call(
            cbind,
            lapply(x@data, function(m) {
                as.integer(m[i, ])
            })
        )
        s_expr <- s_expr[gdata$sample_order, , drop = FALSE]
        colnames(s_expr) <- names(x@data)

        s_disp <- vapply(disp$subgenome, function(phi) phi[i, ], numeric(n_conditions))
        dimnames(s_disp) <- list(gdata$names, names(x@data))

        x_list[[i]] <- list(
            N_SUBGENOMES = n_subgenomes,
            N_CONDITIONS = n_conditions,
            N_REPLICATES = gdata$n_replicates,
            HOMEOLOG_EXP = s_expr,
            HOMEOLOG_EXP_DISPERSION = s_disp,
            LOG_OFFSET = nf_log_offset,
            EPS = eps,
            .META = list(
                condition_names = gdata$names,
                subgenome_names = names(x@data),
                init_params = .hb.mle.init_params(s_expr, gdata, nf_log_offset),
                valid_disp = all(is.finite(s_disp) & s_disp > 0)
            )
        )
    }

    n_invalid_disp <- sum(!vapply(x_list, function(z) isTRUE(z$.META$valid_disp), logical(1L)))
    if (n_invalid_disp > 0L) {
        warning(sprintf("%d genes have invalid dispersion (infinite or negative values).", n_invalid_disp),
                call. = FALSE)
    }

    x_list
}


# Generate an NA placeholder for estimation results.
#' @importFrom utils combn
.hb.mle.init_stats_tmpl <- function(inputs) {
    # Get all combinations of D_*, LR_*, and RR_* statistics with dummy values.
    effect_stats <- .hb.mle.calc_effects(inputs$.META$init_params$theta1, inputs)

    # Theta parameters.
    theta0_names <- paste0("theta0__", inputs$.META$subgenome_names)
    theta1_names <- unlist(lapply(seq_len(inputs$N_CONDITIONS), function(c) {
            paste0("theta1__", inputs$.META$subgenome_names, "__", inputs$.META$condition_names[c])
        }),
        use.names = FALSE)

    stats_names <- c(
        "pvalue",
        "lrt",
        "df",
        "logLik_H0",
        "logLik_H1",
        "converged_H0",
        "converged_H1",
        "optimizer_code_H0",
        "optimizer_code_H1",
        "max_gradient_H0",
        "max_gradient_H1",
        "lrt_negative",
        "status_code",
        names(effect_stats),
        theta0_names,
        theta1_names
    )
    stats::setNames(rep(NA_real_, length(stats_names)), stats_names)
}


.hb.mle.theta_to_logit <- function(theta, esp = 1e-12) {
    theta <- as.matrix(theta)
    n_theta <- nrow(theta)
    n_subgenomes <- ncol(theta)

    logit_mx <- matrix(NA_real_, nrow = n_theta, ncol = n_subgenomes - 1L)

    for (k in seq_len(n_theta)) {
        p <- pmax(theta[k, ], esp)
        p <- p / sum(p)
        logit_mx[k, ] <- log(p[seq_len(n_subgenomes - 1L)] / p[n_subgenomes])
    }

    logit_mx
}

# Convert data for TMB fitting.
.hb.mle.as_tmb_data <- function(inputs, model = c("H0", "H1"), theta = NULL, gene_exp_mu = NULL) {
    model <- match.arg(model)

    # TMB data
    if (model == "H0") {
        n_theta <- 1L
        theta_id <- rep.int(0L, inputs$N_CONDITIONS)
    } else {
        n_theta <- inputs$N_CONDITIONS
        theta_id <- seq_len(inputs$N_CONDITIONS) - 1L
    }

    tmb_data <- list(
        N_SUBGENOMES = as.integer(inputs$N_SUBGENOMES),
        N_CONDITIONS = as.integer(inputs$N_CONDITIONS),
        N_REPLICATES = as.integer(inputs$N_REPLICATES),
        N_THETA = as.integer(n_theta),
        THETA_ID = as.integer(theta_id),
        HOMEOLOG_EXP = inputs$HOMEOLOG_EXP,
        HOMEOLOG_EXP_DISPERSION = inputs$HOMEOLOG_EXP_DISPERSION,
        LOG_OFFSET = as.numeric(inputs$LOG_OFFSET),
        EPS = as.numeric(inputs$EPS)
    )

    # TMB parameters.
    if (is.null(gene_exp_mu)) {
        gene_exp_mu <- inputs$.META$init_params$gene_exp_mu
    }
    gene_exp_mu <- pmax(as.numeric(gene_exp_mu), inputs$EPS)
    if (is.null(theta)) {
        theta <- if (model == "H0") {
            inputs$.META$init_params$theta0
        } else {
            inputs$.META$init_params$theta1
        }
    }

    tmb_params <- list(log_gene_exp_mu = log(gene_exp_mu),
                       theta_logit = .hb.mle.theta_to_logit(theta))

    list(data = tmb_data, params = tmb_params)
}


# Return the model with the maximum log-likelihood.
.hb.mle.best_fit <- function(fit_results) {
    valid <- vapply(fit_results, function(z) isTRUE(z$success), logical(1L))
    if (any(valid)) {
        valid_results <- fit_results[valid]
        log_lik <- vapply(valid_results, `[[`, numeric(1L), "log_lik")
        return(valid_results[[which.max(log_lik)]])
    }

    finite <- vapply(fit_results, function(z) length(z$log_lik) == 1L && is.finite(z$log_lik), logical(1L))
    if (any(finite)) {
        finite_results <- fit_results[finite]
        log_lik <- vapply(finite_results, `[[`, numeric(1L), "log_lik")
        return(finite_results[[which.max(log_lik)]])
    }

    fit_results[[1L]]
}


.hb.mle.update_stats <- function(input, out, fit_result, model) {
    out[paste0("converged_", model)] <- as.numeric(fit_result$success)
    out[paste0("optimizer_code_", model)] <- fit_result$convergence
    out[paste0("max_gradient_", model)] <- fit_result$max_gradient
    out[paste0("logLik_", model)] <- fit_result$log_lik
    # Store fitted theta values.
    if (fit_result$success) {
        if (model == 'H0') {
            theta0_names <- paste0("theta0__", input$.META$subgenome_names)
            out[theta0_names] <- fit_result$theta[1L, ]
        } else if (model == 'H1') {
            theta1_names <- unlist(lapply(seq_len(input$N_CONDITIONS), function(c) {
                    paste0("theta1__", input$.META$subgenome_names, "__", input$.META$condition_names[c])
                }), use.names = FALSE)
            out[theta1_names] <- as.vector(t(fit_result$theta))
        }
    }
    # Store fitting status codes.
    if (!fit_result$success) {
        if (model == 'H0') {
            out["status_code"] <- 2
        } else if (model == 'H1') {
            out["status_code"] <- 3
        }
    }
    out
}

.hb.mle.as_report_matrix <- function(x, nrow, ncol, name) {
    if (is.null(x) || length(x) != nrow * ncol || any(!is.finite(x))) {
        stop(sprintf("Invalid `%s` value returned by TMB.", name), call. = FALSE)
    }
    matrix(as.numeric(x), nrow = nrow, ncol = ncol)
}

#' @importFrom TMB MakeADFun FreeADFun
#' @importFrom stats nlminb
.hb.mle.tmbfit <- function(data, nlminb_control, gradient_tolerance, dll = "hespresso") {
    # Generate the objective function.
    obj <- tryCatch(
        TMB::MakeADFun(data = data$data, parameters = data$params, DLL = dll, silent = TRUE),
        error = function(e) e
    )
    if (inherits(obj, "error")) {
        return(list(
            success = FALSE,
            log_lik = NA_real_,
            theta = NULL,
            gene_exp_mu = NULL,
            convergence = NA_integer_,
            max_gradient = NA_real_,
            message = conditionMessage(obj)
        ))
    }

    on.exit(try(TMB::FreeADFun(obj), silent = TRUE), add = TRUE)

    # Fit the model.
    fit <- tryCatch(
        stats::nlminb(start = obj$par,  objective = obj$fn, gradient = obj$gr, control = nlminb_control),
        error = function(e) e
    )
    if (inherits(fit, "error")) {
        return(list(
            success = FALSE,
            log_lik = NA_real_,
            theta = NULL,
            gene_exp_mu = NULL,
            convergence = NA_integer_,
            max_gradient = NA_real_,
            message = conditionMessage(fit)
        ))
    }

    objective <- tryCatch(obj$fn(fit$par), error = function(e) NA_real_)
    gradient <- tryCatch(obj$gr(fit$par), error = function(e) rep(NA_real_, length(fit$par)))
    report <- tryCatch(obj$report(fit$par), error = function(e) NULL)
    max_gradient <- if (length(gradient) > 0L && all(is.finite(gradient))) {
        max(abs(gradient))
    } else {
        NA_real_
    }

    basic_valid <- isTRUE(fit$convergence == 0L) && length(objective) == 1L && is.finite(objective) &&
        is.finite(max_gradient) && max_gradient <= gradient_tolerance && !is.null(report)
    if (!basic_valid) {
        return(list(
            success = FALSE,
            log_lik = if (is.finite(objective)) -objective else NA_real_,
            theta = NULL,
            gene_exp_mu = NULL,
            convergence = as.integer(fit$convergence),
            max_gradient = max_gradient,
            message = fit$message
        ))
    }

    theta <- tryCatch({
        if (is.null(report$theta) || length(report$theta) != data$data$N_THETA * data$data$N_SUBGENOMES || any(!is.finite(report$theta))) {
            stop("Invalid theta value returned by TMB.", call. = FALSE)
        }
        matrix(as.numeric(report$theta), nrow = data$data$N_THETA, ncol = data$data$N_SUBGENOMES)
    }, error = function(e) NULL)

    gene_exp_mu <- tryCatch(as.numeric(report$gene_exp_mu), error = function(e) NULL)

    report_valid <-
        !is.null(theta) &&
        identical(
            dim(theta),
            c(
                data$data$N_THETA,
                data$data$N_SUBGENOMES
            )
        ) &&
        length(gene_exp_mu) == data$data$N_CONDITIONS &&
        all(is.finite(gene_exp_mu)) &&
        all(gene_exp_mu > 0)

    if (!isTRUE(report_valid)) {
        theta <- NULL
        gene_exp_mu <- NULL
    }

    list(
        success = isTRUE(report_valid),
        log_lik = -objective,
        theta = theta,
        gene_exp_mu = gene_exp_mu,
        convergence = as.integer(fit$convergence),
        max_gradient = max_gradient,
        message = fit$message
    )
}



# Compute homeolog expression ratio differences from observed counts.
#
#' @importFrom utils combn
#' @importFrom stats setNames
.hb.mle.calc_effects <- function(theta1, inputs, eps = 1e-12) {
    c_names <- inputs$.META$condition_names
    s_names <- inputs$.META$subgenome_names

    theta1 <- pmax(theta1, eps)
    theta1 <- theta1 / rowSums(theta1)

    # Ratio differences.
    d <- numeric()
    d_names <- character()
    for (s in seq_len(inputs$N_SUBGENOMES)) {
        for (i in seq_len(inputs$N_CONDITIONS - 1L)) {
            for (j in seq.int(i + 1L, inputs$N_CONDITIONS)) {
                d <- c(d, theta1[i, s] - theta1[j, s])
                d_names <- c(d_names, paste0("D__", s_names[s], "__(", c_names[i], "-", c_names[j], ")"))
            }
        }
    }

    # Ratio fold changes.
    log_ratio_shift <- numeric()
    ratio_of_ratios <- numeric()
    lr_names <- character()
    rr_names <- character()
    subgenome_pairs <- utils::combn(seq_len(inputs$N_SUBGENOMES), 2L)
    for (pair_id in seq_len(ncol(subgenome_pairs))) {
        a <- subgenome_pairs[1L, pair_id]
        b <- subgenome_pairs[2L, pair_id]
        pair_name <- paste0(s_names[a], "/", s_names[b])
        for (i in seq_len(inputs$N_CONDITIONS - 1L)) {
            for (j in seq.int(i + 1L, inputs$N_CONDITIONS)) {
                lr <- log(theta1[i, a] / theta1[i, b]) - log(theta1[j, a] / theta1[j, b])
                log_ratio_shift <- c(log_ratio_shift, lr)
                ratio_of_ratios <- c(ratio_of_ratios, exp(lr))
                lr_names <- c(lr_names, paste0("LR__", pair_name, "__(", c_names[i], "-", c_names[j], ")"))
                rr_names <- c(rr_names, paste0("RR__", pair_name, "__(", c_names[i], "/", c_names[j], ")"))
            }
        }
    }

    c(stats::setNames(d, d_names),
      stats::setNames(log_ratio_shift, lr_names),
      stats::setNames(ratio_of_ratios, rr_names),
      Dmax = max(abs(d), na.rm = TRUE),
      LRmax = max(abs(log_ratio_shift), na.rm = TRUE),
      RRmax = exp(max(abs(log_ratio_shift), na.rm = TRUE)))
}




#' @importFrom stats pchisq
.hb.mle.fit_gene <- function(inputs, nlminb_control, gradient_tolerance, lrt_tolerance) {
    out <- .hb.mle.init_stats_tmpl(inputs)
    out["df"] <- (inputs$N_CONDITIONS - 1L) * (inputs$N_SUBGENOMES - 1L)
    out["lrt_negative"] <- 0

    if (isTRUE(inputs$.META$init_params$zero_condition)) {
        out["converged_H0"] <- 0
        out["converged_H1"] <- 0
        out["status_code"] <- 1
        return(out)
    }

    if (!isTRUE(inputs$.META$valid_disp)) {
        out["converged_H0"] <- 0
        out["converged_H1"] <- 0
        out["status_code"] <- 5
        return(out)
    }

    # H0
    uniform_theta <- matrix(1 / inputs$N_SUBGENOMES, nrow = 1L, ncol = inputs$N_SUBGENOMES)
    data_h0_uniform <- .hb.mle.as_tmb_data(inputs, "H0", theta = uniform_theta)
    data_h0_prior <- .hb.mle.as_tmb_data(inputs, "H0", theta = NULL)

    fit_h0 <- list(
        .hb.mle.tmbfit(data_h0_prior, nlminb_control, gradient_tolerance),
        .hb.mle.tmbfit(data_h0_uniform, nlminb_control, gradient_tolerance))
    fit_h0 <- .hb.mle.best_fit(fit_h0)
    out <- .hb.mle.update_stats(inputs, out, fit_h0, 'H0')

    if (!fit_h0$success) {
        return(out)
    }

    # H1
    nested_theta <- matrix(rep(fit_h0$theta[1L, ], times = inputs$N_CONDITIONS),
                           nrow = inputs$N_CONDITIONS, ncol = inputs$N_SUBGENOMES, byrow = TRUE)
    data_h1_nested <- .hb.mle.as_tmb_data(inputs, 'H1', theta = nested_theta, gene_exp_mu = fit_h0$gene_exp_mu)
    data_h1_prior <- .hb.mle.as_tmb_data(inputs, "H1", theta = NULL)

    fit_h1 <- list(
        .hb.mle.tmbfit(data_h1_nested, nlminb_control, gradient_tolerance),
        .hb.mle.tmbfit(data_h1_prior, nlminb_control, gradient_tolerance))
    fit_h1 <- .hb.mle.best_fit(fit_h1)
    out <- .hb.mle.update_stats(inputs, out, fit_h1, 'H1')

    if (!fit_h1$success) {
        return(out)
    }

    # Calculate ratio differences.
    r_effects <- .hb.mle.calc_effects(fit_h1$theta, inputs)
    out[names(r_effects)] <- r_effects

    # LRT
    lrt <- 2 * (fit_h1$log_lik - fit_h0$log_lik)

    if (!is.finite(lrt)) {
        out["status_code"] <- 4
        return(out)
    }

    lrt_negative <- (is.finite(lrt)) && (lrt < -lrt_tolerance)
    out["lrt_negative"] <- as.numeric(lrt_negative)
    if (lrt_negative) {
        out["lrt"] <- lrt
        out["status_code"] <- 4
        return(out)
    }

    lrt <- max(lrt, 0)
    out["lrt"] <- lrt
    out["pvalue"] <- stats::pchisq(lrt, df = out["df"], lower.tail = FALSE)
    out["status_code"] <- 0

    out
}


#' HOBIT Using MLE
#'
#' Tests whether homeolog expression ratios differ among experimental
#' conditions by fitting null and alternative negative-binomial models with
#' maximum likelihood estimation.
#'
#' HOBIT MLE fits the null and alternative models independently and compares
#' them with a likelihood ratio test. Detailed methodological background and
#' usage guidance are available at
#' \url{https://bitdessin.github.io/hespresso/}.
#'
#' @param x An \linkS4class{SeqCountData} object containing raw integer RNA-seq counts
#'   for corresponding homeologs.
#' @param no_replicate Logical. If `FALSE`, a common tagwise dispersion
#'   is estimated for each homeolog using all samples and a design matrix
#'   containing the experimental group. The estimated dispersion is shared
#'   across conditions. If `TRUE`, a pooled dispersion estimate is used as
#'   an approximation for unreplicated data.
#' @param eps Positive numeric lower bound applied to expected homeolog counts
#'   for numerical stability.
#' @param n_threads Positive integer specifying the number of parallel R workers
#'   used across homeolog tuples.
#' @param gradient_tolerance Positive numeric threshold for the maximum absolute
#'   gradient component accepted after optimization.
#' @param lrt_tolerance Non-negative numeric tolerance for a slightly negative
#'   likelihood ratio statistic caused by numerical error. Values between
#'   `-lrt_tolerance` and zero are set to zero. More negative values are treated
#'   as optimization failures.
#' @param .debug Logical. If `TRUE`, return all statistical test results for
#'   detailed analysis or debugging. If `FALSE`, return only the main
#'   statistical test results.
#' @param ... Named arguments added to the `control` argument of
#'   \code{\link[stats]{nlminb}}, such as `iter.max`, `eval.max`, `rel.tol`,
#'   and `x.tol`.
#'
#' @return By default, a data frame with one row per homeolog tuple containing
#'   `gene`, `pvalue`, `qvalue`, `Dmax`, `LRmax`, `RRmax`, and `theta0__*`
#'   columns. When `.debug = TRUE`, the data frame additionally contains:
#'
#'   \itemize{
#'     \item `pvalue`: p-value from the likelihood ratio test;
#'     \item `qvalue`: Benjamini--Hochberg adjusted p-value;
#'     \item `lrt`: likelihood ratio statistic;
#'     \item `df`: asymptotic degrees of freedom;
#'     \item `theta0__*`: null model maximum-likelihood ratio estimates;
#'     \item `theta1__*__*`: condition-specific alternative model estimates;
#'     \item `D__*`: differences in individual subgenome proportions;
#'     \item `LR__*`: changes in pairwise subgenome log-ratios;
#'     \item `RR__*`: exponentiated pairwise log-ratio changes;
#'     \item `logLik_H0` and `logLik_H1`: maximized log-likelihoods;
#'     \item `status`: fitting and testing status.
#'   }
#'
#'   Possible values of `status` include:
#'
#'   \itemize{
#'     \item `success`: both models converged and the test was calculated;
#'     \item `zero_expression`: at least one condition had zero total expression;
#'     \item `dispersion_failed`: a non-finite or non-positive negative-binomial dispersion was obtained;
#'     \item `h0_fit_failed`: null model optimization failed;
#'     \item `h1_fit_failed`: alternative model optimization failed;
#'     \item `lrt_failed`: the optimized alternative model likelihood was
#'       materially smaller than the null model likelihood.
#'   }
#'
#' @seealso [hobit()], [hobit_mcmc()], \code{\link[TMB]{MakeADFun}},
#'   \code{\link[stats]{nlminb}}
#'
#' @examples
#' x <- sim_homeolog_counts(n_genes = 10)
#'
#' output <- hobit_mle(x, n_threads = 1)
#'
#' output <- hobit_mle(x, n_threads = 2, iter.max = 5000, eval.max = 5000)
#'
#' @importFrom stats nlminb p.adjust pchisq setNames
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_vapply
#' @importFrom progressr progressor with_progress
#' @export
hobit_mle <- function(x,
                      no_replicate = FALSE,
                      eps = 1e-3,
                      n_threads = 1L,
                      gradient_tolerance = 1e-2,
                      lrt_tolerance = 1e-6,
                      .debug = FALSE,
                      ...) {
    # Validate parameters.
    n_threads <- .as_positive_int(n_threads, 'n_threads')
    eps <- .as_positive_float(eps, "eps")
    gradient_tolerance <- .as_positive_float(gradient_tolerance, 'gradient_tolerance')
    lrt_tolerance <- .as_nonnegative_float(lrt_tolerance, 'lrt_tolerance')
    nlminb_control <- .hb.mle.init_nlminb_control(list(...))

    # Format data.
    data <- .hb.mle.format_data(x, eps, no_replicate)

    # Set output names.
    stats_names <- names(.hb.mle.init_stats_tmpl(data[[1L]]))
    fun_value <- stats::setNames(numeric(length(stats_names)), stats_names)

    # Run HOBIT estimation across multiple processes.
    future_plan_bkup <- future::plan()
    on.exit(future::plan(future_plan_bkup), add = TRUE)

    if (n_threads > 1L) {
        future::plan(future::multisession, workers = n_threads)
    } else {
        future::plan(future::sequential)
    }

    stats <- progressr::with_progress({
        pb <- progressr::progressor(length(data))
        future.apply::future_vapply(
            seq_along(data),
            function(i) {
                on.exit(pb(message = 'HOBIT (MEL)'), add = TRUE)
                fit_output <- tryCatch(
                    .hb.mle.fit_gene(data[[i]], nlminb_control, gradient_tolerance, lrt_tolerance),
                    error = function(e) {
                        warning(sprintf("MLE fitting failed for gene `%s`: %s",
                                x@gene_names[i], conditionMessage(e)),
                            call. = FALSE)
                        fit_errout <- .hb.mle.init_stats_tmpl(data[[i]])
                        fit_errout["converged_H0"] <- 0
                        fit_errout["converged_H1"] <- 0
                        fit_errout["status_code"] <- 6
                        fit_errout
                    })
                fit_output
            },
            FUN.VALUE = fun_value,
            future.seed = FALSE,
            future.packages = c('hespresso', 'TMB')
        )
        },
        handlers = progressr::handler_progress(format = ":message :percent [:bar] :current/:total :eta")
    )

    stats <- data.frame(t(stats), check.names = FALSE)
    stats$qvalue <- stats::p.adjust(stats$pvalue, method = "BH")

    status_labels <- c(`0` = "success",
                       `1` = "zero_expression",
                       `2` = "h0_fit_failed", `3` = "h1_fit_failed",
                       `4` = "lrt_failed",
                       `5` = "dispersion_failed",
                       `6` = "internal_error")
    status <- unname(status_labels[as.character(stats$status_code)])

    stats$converged_H0 <- as.logical(stats$converged_H0)
    stats$converged_H1 <- as.logical(stats$converged_H1)
    stats$lrt_negative <- as.logical(stats$lrt_negative)

    output <- data.frame(gene = x@gene_names,
        pvalue = stats$pvalue, qvalue = stats$qvalue,
        lrt = stats$lrt, df = stats$df,
        status = status,
        converged_H0 = stats$converged_H0, converged_H1 = stats$converged_H1,
        optimizer_code_H0 = stats$optimizer_code_H0, optimizer_code_H1 = stats$optimizer_code_H1,
        max_gradient_H0 = stats$max_gradient_H0, max_gradient_H1 = stats$max_gradient_H1,
        lrt_negative = stats$lrt_negative,
        stats[, grep("^D__", names(stats), value = TRUE), drop = FALSE],
        stats[, grep("^LR__", names(stats), value = TRUE), drop = FALSE],
        stats[, grep("^RR__", names(stats), value = TRUE), drop = FALSE],
        Dmax = stats$Dmax,
        LRmax = stats$LRmax,
        RRmax = stats$RRmax,
        stats[, grep("^theta0__", names(stats), value = TRUE), drop = FALSE],
        stats[, grep("^theta1__", names(stats), value = TRUE), drop = FALSE],
        logLik_H0 = stats$logLik_H0,
        logLik_H1 = stats$logLik_H1,
        check.names = FALSE
    )

    .hobit.select_output(output, .debug)
}
