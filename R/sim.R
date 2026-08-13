# Load the empirical matrix used to construct the mean-dispersion population.
.init_seed_counts <- function(expmx, n_subgenomes) {
    if (is.null(expmx)) {
        seed_file <- if (n_subgenomes == 2L) {
            "seed_matrix.C_flexuosa.tsv.gz"
        } else if (n_subgenomes == 3L) {
            "seed_matrix.T_aestivum.tsv.gz"
        } else {
            stop("For more than three subgenomes, provide `seed_counts` explicitly.", call. = FALSE)
        }

        expmx <- utils::read.table(system.file("extdata", seed_file, package = "hespresso"), check.names = FALSE)
    } else if (!is.matrix(expmx) && !is.data.frame(expmx)) {
        stop("`seed_counts` must be `NULL`, a matrix, or a data.frame.", call. = FALSE)
    }

    expmx <- as.matrix(expmx)
    storage.mode(expmx) <- "double"

    if (nrow(expmx) < 3L || ncol(expmx) < 2L) {
        stop("`seed_counts` must contain at least three rows and two samples.", call. = FALSE)
    }
    if (anyNA(expmx) || any(!is.finite(expmx)) || any(expmx < 0)) {
        stop("`seed_counts` must contain finite, non-negative values.", call. = FALSE)
    }
    if (any(colSums(expmx) <= 0)) {
        stop("Every column of `seed_counts` must have a positive sum.", call. = FALSE)
    }

    expmx
}


# Construct an empirical population of total-expression means and NB
# dispersions, then sample baseline means for the requested number of genes.
.get_seed_params <- function(n_genes, seed_matrix) {
    lib_size <- colSums(seed_matrix)
    target_lib_size <- exp(mean(log(lib_size)))
    scaled <- sweep(seed_matrix, 2L, target_lib_size / lib_size, "*")

    mean_value <- rowMeans(scaled)
    variance_value <- apply(scaled, 1L, stats::var)
    dispersion_value <- (variance_value - mean_value) / pmax(mean_value^2, 1e-12)

    valid <- is.finite(mean_value) & mean_value >= 1 &
        is.finite(variance_value) & variance_value > 0 &
        is.finite(dispersion_value) & dispersion_value > 0

    population <- data.frame(
        mean = mean_value[valid],
        variance = variance_value[valid],
        dispersion = dispersion_value[valid]
    )

    if (nrow(population) < 20L) {
        stop("The seed matrix did not provide enough genes with positive empirical dispersion.", call. = FALSE)
    }

    sampled_id <- sample.int(nrow(population), size = n_genes, replace = TRUE)
    sampled <- population[sampled_id, , drop = FALSE]
    rownames(sampled) <- NULL

    list(population = population, sample = sampled)
}


# Fit a robust log-dispersion trend as a quadratic function of log mean.
.fit_mean_dispersion <- function(population) {
    model_data <- population[
        is.finite(population$mean) & population$mean > 0 &
            is.finite(population$dispersion) & population$dispersion > 0,
        ,
        drop = FALSE
    ]

    model_data$log_mean <- log10(model_data$mean)
    model_data$log_dispersion <- log10(model_data$dispersion)

    lower <- stats::quantile(model_data$log_dispersion, 0.005, na.rm = TRUE, names = FALSE)
    upper <- stats::quantile(model_data$log_dispersion, 0.995, na.rm = TRUE, names = FALSE)
    model_data <- model_data[
        model_data$log_dispersion >= lower & model_data$log_dispersion <= upper, ,
        drop = FALSE]

    if (nrow(model_data) >= 20L && length(unique(model_data$log_mean)) >= 3L) {
        fit <- stats::lm(log_dispersion ~ log_mean + I(log_mean^2), data = model_data)
    } else {
        fit <- stats::lm(log_dispersion ~ 1, data = model_data)
    }

    residual_sd <- stats::mad(stats::residuals(fit), constant = 1.4826, na.rm = TRUE)
    if (!is.finite(residual_sd) || residual_sd <= 0) {
        residual_sd <- stats::sd(stats::residuals(fit), na.rm = TRUE)
    }
    if (!is.finite(residual_sd) || residual_sd <= 0) {
        residual_sd <- 0.25
    }

    dispersion_quantiles <- stats::quantile(
        model_data$dispersion, probs = c(0.005, 0.995), na.rm = TRUE, names = FALSE)

    dispersion_floor <- max(1e-8, dispersion_quantiles[1L] / 2)
    dispersion_ceiling <- min(100, max(
        dispersion_quantiles[2L] * 2,
        dispersion_floor * 10
    ))

    list(fit = fit,
         residual_sd = residual_sd,
         log_mean_range = range(model_data$log_mean),
         dispersion_floor = dispersion_floor,
         dispersion_ceiling = dispersion_ceiling)
}


# Predict positive NB dispersions from expected homeolog means.
.predict_dispersion <- function(model, mean_value) {
    mean_value <- as.numeric(mean_value)
    if (any(!is.finite(mean_value)) || any(mean_value <= 0)) {
        stop("Dispersion prediction requires finite, positive means.", call. = FALSE)
    }

    log_mean <- log10(mean_value)
    log_mean <- pmin(pmax(log_mean, model$log_mean_range[1L]), model$log_mean_range[2L])

    predicted <- as.numeric(stats::predict(model$fit, newdata = data.frame(log_mean = log_mean)))

    log_dispersion <- stats::rnorm(length(mean_value), mean = predicted, sd = model$residual_sd)
    dispersion <- 10^log_dispersion
    dispersion <- pmin(pmax(dispersion, model$dispersion_floor), model$dispersion_ceiling)
    dispersion
}


# Sample exactly round(prop * n) true entries.
.sample_indicator <- function(n, prop) {
    n_true <- as.integer(round(n * prop))
    out <- rep(FALSE, n)
    if (n_true > 0L) {
        out[sample.int(n, size = n_true, replace = FALSE)] <- TRUE
    }
    out
}


# Independently generate total-expression DEG truth and HER-shift truth.
.sample_truth <- function(n_genes, group_names, settings) {
    is_shift <- .sample_indicator(n_genes, settings$p_shift)
    is_deg <- .sample_indicator(n_genes, settings$p_deg)

    shift_group <- rep("", n_genes)
    deg_group <- rep("", n_genes)
    log2fc <- numeric(n_genes)

    if (any(is_shift)) {
        shift_group[is_shift] <- sample(group_names, sum(is_shift), replace = TRUE)
    }
    if (any(is_deg)) {
        deg_group[is_deg] <- sample(group_names, sum(is_deg), replace = TRUE)
        magnitude <- pmax(
            stats::rnorm(sum(is_deg), mean = settings$deg_log2fc_mean, sd = settings$deg_log2fc_sd),
            settings$deg_log2fc_min)
        direction <- sample(c(-1, 1), sum(is_deg), replace = TRUE)
        log2fc[is_deg] <- direction * magnitude
    }

    list(is_shift = is_shift, is_deg = is_deg, shift_group = shift_group, deg_group = deg_group, log2fc = log2fc)
}


# Stable softmax for one numeric vector.
.softmax <- function(x, floor = 1e-8) {
    x <- as.numeric(x)
    x <- x - max(x)
    p <- exp(x)
    p <- p / sum(p)
    p <- pmax(p, floor)
    p / sum(p)
}


# Draw baseline HER vectors from a symmetric Dirichlet distribution whose
# concentration varies among genes.
.sample_baseline_her <- function(n_genes, n_subgenomes, gene_names, subgenome_names, settings) {
    total_concentration <- exp(stats::rnorm(n_genes,
        mean = log(settings$baseline_her_concentration),
        sd = settings$baseline_her_concentration_sd))

    theta <- matrix(NA_real_, nrow = n_genes, ncol = n_subgenomes,
                    dimnames = list(gene_names, subgenome_names))

    for (i in seq_len(n_genes)) {
        alpha <- rep(total_concentration[i] / n_subgenomes, n_subgenomes)
        draw <- stats::rgamma(n_subgenomes, shape = alpha, rate = 1)
        theta[i, ] <- draw / sum(draw)
    }

    theta
}


# Create exact-H0 and explicit-H1 HER matrices.  H1 is introduced as a
# zero-sum perturbation on log-ratio coordinates, avoiding a privileged
# reference subgenome.
.generate_her <- function(baseline_her, truth, group_names, settings) {
    n_genes <- nrow(baseline_her)
    n_subgenomes <- ncol(baseline_her)

    her <- lapply(group_names, function(g) baseline_her)
    names(her) <- group_names

    shift_strength <- numeric(n_genes)
    structural_zero <- rep(FALSE, n_genes)
    structural_zero_group <- rep(NA_character_, n_genes)

    shifted_ids <- which(truth$is_shift)
    for (i in shifted_ids) {
        delta <- stats::rnorm(n_subgenomes)
        delta <- delta - mean(delta)
        norm_delta <- sqrt(sum(delta^2))
        if (!is.finite(norm_delta) || norm_delta <= 0) {
            delta <- c(1, rep(-1 / (n_subgenomes - 1), n_subgenomes - 1))
            norm_delta <- sqrt(sum(delta^2))
        }
        delta <- delta / norm_delta

        magnitude <- settings$shift_magnitude_min + stats::rgamma(1L,
            shape = settings$shift_magnitude_shape,
            scale = settings$shift_magnitude_scale
        )
        shift_strength[i] <- magnitude

        shifted_theta <- .softmax(
            log(pmax(baseline_her[i, ], settings$theta_floor)) +
                magnitude * delta,
            floor = settings$theta_floor
        )

        her[[truth$shift_group[i]]][i, ] <- shifted_theta
    }

    if (settings$p_zero > 0) {
        for (g in seq_along(group_names)) {
            boundary_ids <- which(stats::rbinom(
                n_genes, size = 1L, prob = settings$p_zero
            ) == 1L)
            for (i in boundary_ids) {
                absent <- sample.int(n_subgenomes, size = 1L)
                theta <- her[[g]][i, ]
                theta[absent] <- 0
                theta <- theta / sum(theta)
                her[[g]][i, ] <- theta
                structural_zero[i] <- TRUE
                if (is.na(structural_zero_group[i])) {
                    structural_zero_group[i] <- group_names[g]
                }
            }
        }
    }

    list(
        her = her,
        shift_strength = shift_strength,
        structural_zero = structural_zero,
        structural_zero_group = structural_zero_group
    )
}


# Generate group-specific total expression means independently of HER.
.generate_total_mu <- function(baseline_total_mean, truth, group_names, gene_names) {
    total_mu <- lapply(group_names, function(g) {
        stats::setNames(as.numeric(baseline_total_mean), gene_names)
    })
    names(total_mu) <- group_names
    deg_ids <- which(truth$is_deg)
    for (i in deg_ids) {
        g <- truth$deg_group[i]
        total_mu[[g]][i] <- total_mu[[g]][i] * 2^truth$log2fc[i]
    }

    total_mu
}


# Combine total-expression means and HERs:
#   mu[g, s] = total_mu[g] * theta[g, s].
.compose_homeolog_mu <- function(total_mu, her) {
    out <- vector("list", length(her))
    names(out) <- names(her)

    for (g in seq_along(her)) {
        out[[g]] <- sweep(her[[g]], 1L, total_mu[[g]], "*")
    }

    out
}


# Generate centered sample-specific log library-size offsets.
.generate_log_offset <- function(n_replicates, group_names, settings) {
    n_samples <- sum(n_replicates)
    log_offset <- stats::rnorm(n_samples, mean = 0, sd = settings$offset_sd)
    log_offset <- log_offset - mean(log_offset)
    sample_names <- paste(rep(group_names, times = n_replicates), sequence(n_replicates), sep = "__")
    names(log_offset) <- sample_names
    log_offset
}


# Sample one shared dispersion per gene and subgenome, using the average
# expected homeolog expression across conditions as the trend covariate.
.generate_dispersion <- function(mu, dispersion_model, subgenome_names) {
    n_groups <- length(mu)
    n_subgenomes <- ncol(mu[[1L]])

    dispersion <- vector("list", n_subgenomes)
    names(dispersion) <- subgenome_names
    for (s in seq_len(n_subgenomes)) {
        mean_value <- Reduce("+", lapply(mu, function(m) m[, s])) / n_groups
        dispersion[[s]] <- .predict_dispersion(dispersion_model, pmax(mean_value, 1e-8))
        names(dispersion[[s]]) <- rownames(mu[[1L]])
    }

    dispersion
}


# Draw raw integer counts exactly once from the model-matched NB distribution.
.sample_counts <- function(mu, dispersion, log_offset, n_replicates, group_names, gene_names, subgenome_names) {
    n_genes <- length(gene_names)
    n_subgenomes <- length(subgenome_names)
    n_samples <- sum(n_replicates)
    sample_names <- names(log_offset)

    counts <- lapply(seq_len(n_subgenomes), function(s) {
        matrix(0, nrow = n_genes, ncol = n_samples,
               dimnames = list(gene_names, sample_names))
    })
    names(counts) <- subgenome_names

    r <- 1L
    for (g in seq_along(group_names)) {
        for (j in seq_len(n_replicates[g])) {
            offset_factor <- exp(log_offset[r])
            for (s in seq_len(n_subgenomes)) {
                expected <- mu[[g]][, s] * offset_factor
                if (any(!is.finite(expected)) || any(expected < 0)) {
                    stop("The simulated expected counts became invalid.", call. = FALSE)
                }
                counts[[s]][, r] <- stats::rnbinom(n = n_genes, mu = expected, size = 1 / dispersion[[s]])
            }
            r <- r + 1L
        }
    }

    counts
}


# Calculate true Dmax, pairwise log-ratio/RR maxima, and one-vs-rest ORmax.
.calc_true_effects <- function(her, eps = 1e-12) {
    n_groups <- length(her)
    n_genes <- nrow(her[[1L]])
    n_subgenomes <- ncol(her[[1L]])

    dmax <- numeric(n_genes)
    lrmax <- numeric(n_genes)
    ormax <- rep(1, n_genes)

    subgenome_pairs <- utils::combn(seq_len(n_subgenomes), 2L)

    for (g1 in seq_len(n_groups - 1L)) {
        for (g2 in seq.int(g1 + 1L, n_groups)) {
            theta1 <- pmax(her[[g1]], eps)
            theta2 <- pmax(her[[g2]], eps)
            theta1 <- theta1 / rowSums(theta1)
            theta2 <- theta2 / rowSums(theta2)

            dmax <- pmax(dmax, apply(abs(theta1 - theta2), 1L, max))

            odds1 <- theta1 / pmax(1 - theta1, eps)
            odds2 <- theta2 / pmax(1 - theta2, eps)
            one_vs_rest_or <- odds1 / odds2
            ormax <- pmax(
                ormax,
                apply(pmax(one_vs_rest_or, 1 / one_vs_rest_or), 1L, max)
            )

            for (pair_id in seq_len(ncol(subgenome_pairs))) {
                a <- subgenome_pairs[1L, pair_id]
                b <- subgenome_pairs[2L, pair_id]
                lr <- log(theta1[, a] / theta1[, b]) -
                    log(theta2[, a] / theta2[, b])
                lrmax <- pmax(lrmax, abs(lr))
            }
        }
    }

    list(dmax = dmax, lrmax = lrmax, rrmax = exp(lrmax), ormax = ormax)
}


#' Generate Artificial RNA-seq Read Counts for Homeolog-Expression Simulation
#'
#' Generates model-matched raw integer RNA-seq counts for evaluating HOBIT.
#' The simulator separates total expression from homeolog expression ratios,
#' creates exact null genes and explicit ratio-shift genes, uses one shared
#' dispersion per gene and subgenome across conditions, adds sample-specific
#' library-size offsets, and draws each observation once from a negative-
#' binomial distribution.
#'
#' The returned `SeqCountData` object contains the simulated count matrices and an
#' associated `SimParams` object in `x@meta`. The latter records the generating
#' parameters, exact H0/H1 labels, and true effect-size summaries.
#'
#' @param n_genes Integer. Number of homeolog tuples to simulate.
#' @param n_replicates Integer vector. Number of replicates for each condition.
#' @param n_subgenomes Integer. Number of subgenomes to simulate.
#' @param p_shift Numeric between 0 and 1. Proportion of genes assigned to the
#'   exact homeolog-expression-ratio shift alternative. The default is 0.10.
#' @param p_deg Numeric between 0 and 1. Proportion of genes assigned a
#'   total-expression change independently of homeolog-expression-ratio status.
#'   The default is 0.10.
#' @param p_zero Numeric between 0 and 1. Proportion of genes
#'   assigned a condition-by-subgenome structural zero, producing homeolog
#'   expression ratios near 0 or 1.
#' @param group_names Character vector of condition names.
#' @param subgenome_names Character vector of subgenome names.
#' @param seed_counts `NULL`, matrix, or data.frame used to construct the
#'   empirical mean-dispersion population. If `NULL`, bundled reference data
#'   are used for two or three subgenomes.
#' @return An `SeqCountData` object containing raw integer counts. Its `meta` slot is a
#'   `SimParams` object containing the exact homeolog-expression-ratio labels
#'   (`is_shift`), independent total-expression labels (`is_deg`), true means,
#'   ratios, dispersions, library-size offsets, and effect-size summaries.
#'
#' @examples
#' x <- sim_homeolog_counts()
#' x <- sim_homeolog_counts(
#'     n_genes = 100,
#'     n_replicates = c(3, 3),
#'     n_subgenomes = 2
#' )
#'
#' table(x@meta@is_shift)
#' table(x@meta@is_deg, x@meta@is_shift)
#'
#' # Simulate a 1% homeolog-expression-ratio shift prevalence.
#' x_rare <- sim_homeolog_counts(n_genes = 10000, p_shift = 0.01)
#'
#' @importFrom methods new
#' @export
sim_homeolog_counts <- function(n_genes = 10000,
                                n_replicates = c(3, 3),
                                n_subgenomes = 2,
                                p_shift = 0.10,
                                p_deg = 0.10,
                                p_zero = 0.05,
                                group_names = NULL,
                                subgenome_names = NULL,
                                seed_counts = NULL) {
    # validate the design before deriving labels or sampling parameters.
    n_genes <- .as_positive_int(n_genes, "n_genes")
    n_replicates <- vapply(seq_along(n_replicates),
        function(ri) .as_positive_int(n_replicates[[ri]], paste0("n_replicates[", ri, "]")), integer(1L))
    n_groups <- .as_positive_int(length(n_replicates), "n_replicates", min_threshold = 2L)
    n_subgenomes <- .as_positive_int(n_subgenomes, "n_subgenomes", min_threshold = 2L)
    p_shift <- .as_prob(p_shift, "p_shift")
    p_deg <- .as_prob(p_deg, "p_deg")
    p_zero <- .as_prob(p_zero, "p_zero")

    # init labels
    if (is.null(group_names)) {
        group_names <- paste0("group_", seq_len(n_groups))
    }
    if (is.null(subgenome_names)) {
        if (n_subgenomes <= length(LETTERS)) {
            subgenome_names <- LETTERS[seq_len(n_subgenomes)]
        } else {
            subgenome_names <- paste0("S", seq_len(n_subgenomes))
        }
    }
    gene_names <- paste0("gene_", seq_len(n_genes))

    seed_counts <- .init_seed_counts(seed_counts, n_subgenomes)
    seed_params <- .get_seed_params(n_genes, seed_counts)
    rownames(seed_params$sample) <- gene_names

    dispersion_model <- .fit_mean_dispersion(seed_params$population)
    settings <- list(
        p_shift = p_shift,
        p_deg = p_deg,
        baseline_her_concentration = 12,
        baseline_her_concentration_sd = 0.45,
        shift_magnitude_min = 0.45,
        shift_magnitude_shape = 2,
        shift_magnitude_scale = 0.35,
        deg_log2fc_mean = 1,
        deg_log2fc_sd = 0.35,
        deg_log2fc_min = 0.50,
        offset_sd = 0.35,
        theta_floor = 1e-8,
        dispersion_shared_across_conditions = TRUE,
        p_zero = p_zero,
        n_subgenomes = n_subgenomes
    )
    truth <- .sample_truth(n_genes, group_names, settings)

    baseline_her <- .sample_baseline_her(
        n_genes = n_genes,
        n_subgenomes = n_subgenomes,
        gene_names = gene_names,
        subgenome_names = subgenome_names,
        settings = settings
    )

    her_result <- .generate_her(
        baseline_her = baseline_her,
        truth = truth,
        group_names = group_names,
        settings = settings
    )
    her <- her_result$her
    boundary_ids <- her_result$structural_zero & !truth$is_shift
    truth$is_shift[boundary_ids] <- TRUE
    truth$shift_group[boundary_ids] <- her_result$structural_zero_group[boundary_ids]

    total_mu <- .generate_total_mu(
        baseline_total_mean = seed_params$sample$mean,
        truth = truth,
        group_names = group_names,
        gene_names = gene_names
    )

    mu <- .compose_homeolog_mu(total_mu, her)

    dispersion <- .generate_dispersion(
        mu = mu,
        dispersion_model = dispersion_model,
        subgenome_names = subgenome_names
    )

    log_offset <- .generate_log_offset(
        n_replicates = n_replicates,
        group_names = group_names,
        settings = settings
    )

    counts <- .sample_counts(
        mu = mu,
        dispersion = dispersion,
        log_offset = log_offset,
        n_replicates = n_replicates,
        group_names = group_names,
        gene_names = gene_names,
        subgenome_names = subgenome_names
    )

    true_effects <- .calc_true_effects(her)

    group_vec <- rep(group_names, times = n_replicates)
    group_vec <- factor(group_vec, levels = group_names)

    gene_level_vectors <- list(
        is_shift = truth$is_shift,
        is_deg = truth$is_deg,
        shift_group = truth$shift_group,
        deg_group = truth$deg_group,
        log2fc = truth$log2fc,
        shift_strength = her_result$shift_strength,
        dmax_true = true_effects$dmax,
        lrmax_true = true_effects$lrmax,
        rrmax_true = true_effects$rrmax,
        ormax_true = true_effects$ormax
    )
    gene_level_vectors <- lapply(gene_level_vectors, function(v) {
        names(v) <- gene_names
        v
    })

    methods::new(
        "SeqCountData",
        data = counts,
        gene_names = gene_names,
        exp_design = data.frame(
            group = group_vec,
            replicate = sequence(n_replicates)
        ),
        meta = methods::new(
            "SimParams",
            n_subgenomes = n_subgenomes,
            n_genes = n_genes,
            n_groups = n_groups,
            n_replicates = n_replicates,

            params_population = seed_params$population,
            params = seed_params$sample,
            nls = dispersion_model,

            mu = mu,
            her = her,
            total_mu = total_mu,
            dispersion = dispersion,
            log_offset = log_offset,

            is_shift = gene_level_vectors$is_shift,
            is_deg = gene_level_vectors$is_deg,
            shift_group = gene_level_vectors$shift_group,
            deg_group = gene_level_vectors$deg_group,
            log2fc = gene_level_vectors$log2fc,
            shift_strength = gene_level_vectors$shift_strength,

            dmax_true = gene_level_vectors$dmax_true,
            lrmax_true = gene_level_vectors$lrmax_true,
            rrmax_true = gene_level_vectors$rrmax_true,
            ormax_true = gene_level_vectors$ormax_true,

            settings = settings
        )
    )
}


#' Define Ground-Truth Homeolog-Expression-Ratio Shifts
#'
#' Returns the exact simulated homeolog-expression-ratio shift labels when
#' `Dmax` and `ORmax` are both `NULL`. If `base` is specified, the global truth
#' labels are restricted to genes whose selected subgenome's ratio changes
#' between conditions. When either threshold is supplied, the function
#' calculates ratio differences from the true homeolog expression ratios and
#' filters genes using the supplied threshold or thresholds.
#'
#' @param x An `SeqCountData` object generated by [sim_homeolog_counts()].
#' @param base Integer subgenome index used for thresholding. If `0`, the
#'   maximum value across all subgenomes is used.
#' @param groups Optional character vector of length two specifying the
#'   conditions to compare. If `NULL`, all condition pairs are considered.
#' @param Dmax `NULL` or a numeric threshold for the maximum absolute
#'   homeolog-expression-ratio difference.
#' @param ORmax `NULL` or a numeric threshold for the maximum one-versus-rest
#'   odds ratio of homeolog-expression-ratio changes.
#' @param operator Character string specifying how `Dmax` and `ORmax` criteria
#'   are combined when both are supplied: either `"OR"` or `"AND"`.
#'
#' @return A named logical vector. If both thresholds are `NULL` and `base` is
#'   `0`, it is exactly `x@meta@is_shift`. For a specified base, it contains
#'   the global truth labels restricted to genes whose selected subgenome's
#'   ratio changes between conditions. Otherwise, it indicates genes passing
#'   the supplied threshold or thresholds.
#'
#' @examples
#' x <- sim_homeolog_counts(n_genes = 100)
#'
#' # Exact simulation truth.
#' is_shift <- def_sigShift(x)
#'
#' # Effect-size filtering.
#' is_large_shift <- def_sigShift(x, Dmax = 0.15, ORmax = 2)
#'
#' @seealso [sim_homeolog_counts()]
#' @export
def_sigShift <- function(x,
    base = 0, groups = NULL,
    Dmax = NULL, ORmax = NULL, operator = c('OR', 'AND')) {
    if (is.null(Dmax) && is.null(ORmax)) {
        if (base == 0) {
            return(x@meta@is_shift)
        }

        her <- if (is.null(groups)) {
            x@meta@her
        } else {
            x@meta@her[groups]
        }
        base_shift <- matrix(FALSE, nrow = nrow(her[[1L]]), ncol = 1L)
        if (length(her) > 1L) {
            for (g1 in seq_len(length(her) - 1L)) {
                for (g2 in (g1 + 1L):length(her)) {
                    base_shift[, 1L] <- base_shift[, 1L] |
                        (her[[g1]][, base] != her[[g2]][, base])
                }
            }
        }
        gt <- x@meta@is_shift & base_shift[, 1L]
        names(gt) <- x@gene_names
        return(gt)
    }

    operator <- match.arg(operator)
    ormx <- dmx <- matrix(0, nrow = nrow(x@meta@her[[1]]),
                             ncol = ncol(x@meta@her[[1]]))

    if (!is.null(groups)) {
        dmx <- .calc_exp_dist(x@meta@her[[groups[1]]], x@meta@her[[groups[2]]])
        ormx <- .calc_exp_oddsratio(x@meta@her[[groups[1]]],
                                    x@meta@her[[groups[2]]])
    } else {
        for (g1 in seq_along(x@meta@her)) {
            for (g2 in seq_along(x@meta@her)) {
                .dmx <- .calc_exp_dist(x@meta@her[[g1]], x@meta@her[[g2]])
                .ormx <- .calc_exp_oddsratio(x@meta@her[[g1]], x@meta@her[[g2]])
                for (i in seq_len(ncol(x@meta@her[[1]]))) {
                    dmx[, i] <- apply(cbind(dmx[, i], .dmx[, i]), 1, max, na.rm = TRUE)
                    ormx[, i] <- apply(cbind(ormx[, i], .ormx[, i]), 1, max, na.rm = TRUE)
                }
            }
        }
    }

    criteria <- list()
    if (!is.null(Dmax)) {
        criteria$Dmax <- if (base > 0) {
            dmx[, base] > Dmax
        } else {
            apply(dmx > Dmax, 1, any)
        }
    }
    if (!is.null(ORmax)) {
        criteria$ORmax <- if (base > 0) {
            ormx[, base] > ORmax
        } else {
            apply(ormx > ORmax, 1, any)
        }
    }

    gt <- if (length(criteria) == 1L) {
        criteria[[1L]]
    } else if (operator == 'OR') {
        criteria$Dmax | criteria$ORmax
    } else {
        criteria$Dmax & criteria$ORmax
    }

    names(gt) <- x@gene_names
    gt
}
