# Initialize Seed Expression Matrix
# 
# Initializes a seed expression matrix for creating a population of means and  
# dispersions for simulating count data.
#
# This function checks the provided `seed_matrix`.  
# If a `data.frame` or `matrix` is supplied,
# the input object is returned as is.  
# If `NULL` is provided, the function loads a preset seed expression matrix  
# based on the specified number of subgenomes:
# 
# - If `n_subgenomes` is set to 2, real data quantified
#   from the allotetraploid *Cardamine flexuosa* will be loaded.  
# - If `n_subgenomes` is set to 3, real data quantified
#   from the allohexaploid *Triticum aestivum* (wheat) will be loaded.  
# 
# To use other real data as a population for sampling,
# the user should provide an expression matrix via `seed_matrix`.  
# 
# @param exp_mx `NULL` or an expression matrix in `data.frame` or `matrix`
#                  format used to create a population for sampling.  
# @param n_subgenomes An integer specifying the number of subgenomes.  
# 
# @return An expression matrix in `data.frame` or `matrix` format.
# 
#' @importFrom utils read.table
.init_seed_expmx <- function(expmx, n_subgenomes) {
    if (is.data.frame(expmx) || is.matrix(expmx)) return(expmx)

    if (is.null(expmx)) {
        if (n_subgenomes == 2) {
            expmx <- read.table(system.file(package = 'hespresso', 'extdata',
                                            'seed_matrix.C_flexuosa.tsv.gz'))
        } else if (n_subgenomes == 3) {
            expmx <- read.table(system.file(package = 'hespresso', 'extdata',
                                            'seed_matrix.T_aestivum.tsv.gz'))
        } else {
            stop("Only 2 or 3 subgenomes are supported in current version.",
                 "For other cases,",
                 "please provide a data.frame or matrix object.")
        }
    } else {
        stop("The `seed_expmat`` must be a data.frame or matrix object.")
    }
    expmx
}



# Create Population for Sampling
#
# Creates a population of means and dispersions from a seed matrix.
#
# This function calculates a sampling population of means and dispersions  
# from the expression matrix provided in `seed_matrix`.
# Additionally, it creates a regression model
# to describe the relationship between mean and dispersion  
# using the `nls` function by fitting the calculated population.
#
# @param n_genes An integer specifying the number of genes.  
# @param seed_matrix An expression matrix in `data.frame` or `matrix` format.  
# @param random_seed An integer or `NULL` to set the seed for random sampling.  
#
# @return A list containing two elements:
#  - `population` — the sampling population calculated from `seed_matrix`.  
#  - `sample` — the sampled parameters.
#  
#' @importFrom stats var
.get_seed_params <- function(n_genes, seed_matrix, seed = NULL) {
    # get valid data
    col_sums <- colSums(seed_matrix)
    seed_matrix <- sweep(seed_matrix, 2, mean(col_sums) / col_sums, "*")
    seed_matrix <- seed_matrix[rowMeans(seed_matrix) > 1, ]
    
    # calculate params population
    m <- apply(seed_matrix, 1, mean)
    v <- apply(seed_matrix, 1, var)
    p_pop <- data.frame(mean = m, variance = v, dispersion = (v - m) / (m * m))
    
    list(population = p_pop,
         sample = p_pop[sample(seq_len(nrow(p_pop)),
                               size = n_genes, replace = TRUE), ])
}



# Generate model to predict variance from mean
# 
# @param p_pop A data.frame consisting of mean and variance columns.
#' @importFrom stats nls coef residuals
.nls_mean_dispersion <- function(p_pop) {
    nls_object <- nls(y ~ a * x * x + b * x + c,
                      start = c(a = 1, b = 1, c = 0),
                      data = data.frame(x = log10(p_pop$mean),
                                        y = log10(p_pop$variance)))
    
    nls_variance <- sum(residuals(nls_object) ^ 2) / (nrow(p_pop) - length(coef(nls_object)))
    
    list(nls = nls_object, nls_variance = nls_variance)
}


# Predict dispersion with fitted-model
#' @importFrom stats predict
.nls_mean_dispersion_preidct <- function(object, newx,
                                         type = c("dispersion", "variance")) {
    type <- match.arg(type)
    newx <- as.numeric(newx)
    
    is_0 <- (newx < 1)
    newx[is_0] <- 1
    
    v_m <- predict(object$nls, data.frame(x = log10(newx)))
    v_m[v_m < 0] <- 1e-4
    v_v <- sqrt(object$nls_variance)
    newv <- 10 ^ (rnorm(length(newx), v_m, v_v))

    newd <- (newv - newx) / (newx * newx)
    newd[newd < 0] <- runif(sum(newd < 0), min = 1e-6, max = 1e-2)

    newd[is_0] <- NA
    newv[is_0] <- NA
    
    if (type == "variance") {
        newv
    } else {
        newd
    }
}



# Introduce DEG by multityping randomly sampled fold-changes.
#' @importFrom stats runif rgamma
.sim_counts_DEG <- function(gene_names, subgenome_names, group_names,
                            seed_params) {
    mu_list <- vector('list', length = length(group_names))
    names(mu_list) <- group_names
    for (g in seq_along(group_names)) {
        mu_list[[g]] <- .init_matrix(NA, gene_names, subgenome_names)
    }
    
    fc <- 1 + rgamma(length(gene_names),
                     shape = runif(1, 0.05, 0.5),
                     scale = runif(length(gene_names),
                                   min = runif(1, 0, 0.1) , max = 1))
    
    # select group to apply fold changes
    fc_sides_rand <- rbeta(length(gene_names), 2, runif(1, 2, 6))
    fc_sides <- rep(1, length = length(fc_sides_rand))
    for (g in seq_along(group_names)) {
        fc_sides[fc_sides_rand > ((g - 1) / length(group_names))] <- g
    }
    
    for (g in seq_along(group_names)) {
        for (s in seq_along(subgenome_names)) {
            mu_list[[g]][, s] <-
                seed_params$sample$mean * ifelse(fc_sides == g, fc, 1)
        }
    }

    mu_list
}


# Generate Base HER and Adding Noises to Introduce BEH
#' @importFrom stats rnorm runif rgamma
.sim_counts_HER <- function(gene_names, subgenome_names, group_names) {
    n_genes <- length(gene_names)
    n_genes_2x <- 2 * n_genes
    n_subgenomes <- length(subgenome_names)
    n_groups <- length(group_names)
    
    p_list <- vector('list', length = length(group_names))
    names(p_list) <- group_names
    for (g in seq_along(group_names)) {
        p_list[[g]] <- .init_matrix(NA,
                                    paste0('tmp_', seq_len(n_genes_2x)),
                                    subgenome_names)
    }
    
    # select group to apply ratio changes
    p_sides <- sample(seq_len(n_groups), n_genes_2x, replace = TRUE)
    
    # sampling ratios
    for (s in seq_along(subgenome_names)) {
        # base ratios for s-th subgenome
        p0 <- rnorm(n_genes_2x,
                    mean = 1 / n_subgenomes, sd = (1 / n_subgenomes) ^ 2)
        for (g in seq_along(group_names)) {
            # changes between conditions
            pd <- rgamma(n_genes_2x, shape = 1, scale = 0.1)
            pd <- pd / max(pd) / n_subgenomes * ifelse(runif(n_genes_2x) > 0.5,
                                                       1, -1)
            # base + changes
            p_list[[g]][, s] <- p0 + pd
            # set base as the current condition for generating next condition
            p0 <- p_list[[g]][, s]
        }
    }
    
    # select the valid ratios (0 <= p <= 1)
    is_valid <- rep(TRUE, length = n_genes_2x)
    for (g in seq_along(group_names)) {
        for (s in seq_along(subgenome_names)) {
            is_valid <- ((is_valid) & (0 <= p_list[[g]][, s]))
        }
    }
    idx <- sample(seq_len(n_genes_2x)[is_valid], n_genes, replace = FALSE)
    
    for (g in seq_along(group_names)) {
        p_list[[g]] <- p_list[[g]][idx, ]
        p_list[[g]] <- p_list[[g]] / rowSums(p_list[[g]])
        rownames(p_list[[g]]) <- gene_names
    }
    p_list
}


# Remove outliers
.sim_counts_rmOutliers <- function(mh) {
    n_groups <- length(mh$m)
    n_subgenomes <- ncol(mh$m[[1]])
    n_genes <- nrow(mh$m[[1]])
    
    mu_mx_fc <- 0
    for (g1 in seq_len(n_groups)) {
        for (g2 in seq_len(n_groups)) {
            mu_mx_fc <- max(mu_mx_fc,
                abs(log2(rowSums(mh$m[[g1]]) / rowSums(mh$m[[g2]]))))
        }
    }

    # adjust ground truth means after 2nd subgenomes to reflect HER 
    # (do not change the 1st mean but adjust means after 2nd subgenomes)
    for (g in seq_len(n_groups)) {
        for (s in seq_len(n_subgenomes)) {
            mh$m[[g]][, s] <- mh$m[[g]][, s] * (mh$h[[g]][, s] / mh$h[[g]][, 1])
        }
    }
    
    # resampling mean and HER if outliers
    is_outliers <- rep(FALSE, length = n_genes)
    for (g1 in seq_len(n_groups)) {
        for (g2 in seq_len(n_groups)) {
            is_outliers <- ((is_outliers) | (abs(log2(rowSums(mh$m[[g1]]) / rowSums(mh$m[[g2]]))) > (mu_mx_fc * 1.1)))
        }
    }
    
    if (sum(is_outliers) > 0) {
        mu_sampled <- sample(as.numeric(unlist(mh$m)),
                             size = sum(is_outliers), replace = TRUE)
        for (g in seq_len(n_groups)) {
            mh$m[[g]][is_outliers, ] <- mu_sampled
            mh$h[[g]][is_outliers, ] <- 1 / n_subgenomes
        }
    }
    
    .correct_HER(mh)
}


# Shuffle parameters
.sim_counts_shuffleParams <- function(mh) {
    n_groups <- length(mh$m)
    n_subgenomes <- ncol(mh$m[[1]])
    n_genes <- nrow(mh$m[[1]])
    
    for (i in seq_len(n_genes)) {
        mmx <- hmx <- matrix(NA, ncol = n_subgenomes, nrow = n_groups)
        for (g in seq_len(n_groups)) {
            mmx[g, ] <- mh$m[[g]][i, ]
            hmx[g, ] <- mh$h[[g]][i, ]
        }
        
        c_rnd <- sample(seq_len(ncol(mmx)), replace = FALSE)
        mmx <- mmx[, c_rnd]
        hmx <- hmx[, c_rnd]
        
        r_rnd <- sample(seq_len(nrow(mmx)), replace = FALSE)
        mmx <- mmx[r_rnd, ]
        hmx <- hmx[r_rnd, ]
        
        for (g in seq_len(n_groups)) {
            mh$m[[g]][i, ] <- mmx[g, ]
            mh$h[[g]][i, ] <- hmx[g, ]
        }
    }
    
    .correct_HER(mh)
}


#' @importFrom stats rbeta
.sim_counts_addZeros <- function(mh) {
    n_groups <- length(mh$m)
    n_subgenomes <- ncol(mh$m[[1]])
    n_genes <- nrow(mh$m[[1]])
    
    for (s in seq_len(n_subgenomes)) {
        set_to_0 <- (rbeta(n_genes, 10, runif(1, 2, 3)) < 0.5)
        z_sides <- sample(seq_len(n_groups), n_genes, replace = TRUE)
        for (g in seq_len(n_groups)) {
            mh$m[[g]][set_to_0 & (z_sides == g), s] <- 0
        }
    }
    
    .correct_HER(mh)
}


# Recalculate HER based on mean
.correct_HER <- function(mh) {
    for (g in seq_len(length(mh$m))) {
        mh$h[[g]] <- mh$m[[g]] / rowSums(mh$m[[g]])
        mh$h[[g]][is.nan(mh$h[[g]])] <- NA
    }
    mh
}


#' @importFrom stats rnbinom
.sim_counts_sample1s <- function(n, mu, disp) {
    nb_counts <- rep(0, length = n)
    if (mu > 0) {
        while (all(nb_counts == nb_counts[1])) {
            nb_counts <- rnbinom(n = n, mu = mu, size = 1 / disp)
            if (sum(nb_counts) == 0) break
            if (n == 1) break
        }
    }
    nb_counts
}


.sim_counts <- function(mh, n_replicates, seed_params, nls_object) {
    n_groups <- length(mh$m)
    n_subgenomes <- ncol(mh$m[[1]])
    
    exp_mx_list <- vector('list', n_subgenomes)
    for (s in seq_len(n_subgenomes)) {
        exp_mx_list[[s]] <-
            matrix(NA, ncol = sum(n_replicates), nrow = nrow(mh$m[[1]]))
    }
    
    for (s in seq_len(n_subgenomes)) {
        r <- 1
        for (g in seq_len(n_groups)) {
            mh$m[[g]][, s][mh$m[[g]][, s] < 1] <-
                sample(c(1, 0), sum(mh$m[[g]][, s] < 1), replace = TRUE)
            disp <- .nls_mean_dispersion_preidct(nls_object, mh$m[[g]][, s])
            for (i in seq_len(nrow(mh$m[[1]]))) {
                exp_mx_list[[s]][i, r:(r + n_replicates[g] - 1)] <-
                .sim_counts_sample1s(n_replicates[g], mh$m[[g]][i, s], disp[i])
            }
            r <- r + n_replicates[g]
        }
    }
    
    mh <- .correct_HER(mh)
    list(mh = mh, exp = exp_mx_list)
}



#' Generate Artificial RNA-Seq Read Counts Simulating Homeolog Expression
#'
#' Simulates artificial RNA-Seq read counts for allopolyploids based on a
#' negative binomial distribution, mimicking homeolog expression patterns
#' across subgenomes.
#'
#' This function performs the following steps to generate synthetic read counts:
#' (1) prepares a sampling population,  
#' (2) modifies mean expression values to simulate differentially expressed genes,  
#' (3) introduces variation in homeolog expression ratios,  
#' (4) estimates dispersion parameters based on adjusted means, and  
#' (5) samples counts from a negative binomial distribution using the calculated
#' means and dispersions.
#'
#' By default, the function uses real RNA-Seq datasets from *Cardamine flexuosa*
#' and wheat to define the sampling population for allopolyploids consiting of
#' two or three subgenomes, respectively, in step (1). Alternatively, users can
#' provide their own normalized expression matrix via the `seed_expmx` argument.
#' The matrix should have genes as rows and replicates (from the same condition)
#' as columns.
#'
#' @param n_genes Integer. Number of genes to simulate.
#' @param n_replicates Integer vector. Number of replicates per group.
#' @param n_subgenomes Integer. Number of subgenomes.
#' @param group_names Character vector. Names of the condition groups.
#' @param subgenome_names Character vector. Names of the subgenomes.
#' @param seed_expmx A data.frame or matrix representing a normalized expression
#'      matrix to be used for defining the sampling population.
#'
#' @return An object of class \linkS4class{ExpMX},
#'      containing simulated expression data.
#' @examples
#' x <- sim_homeolog_counts()
#' 
#' # simulating allotetraploids with 3 replicates in each group
#' x <- sim_homeolog_counts(n_genes = 100,
#'                          n_replicates = c(3, 3), n_subgenomes = 2)
#' 
#' # simulating allohexaploids with 5 replicats in each group
#' x <- sim_homeolog_counts(n_genes = 100,
#'                          n_replicates = c(5, 5), n_subgenomes =  3)
#' @importFrom methods new
#' @export
sim_homeolog_counts <- function(n_genes = 10000, n_replicates = c(3, 3),
                                n_subgenomes = 2,
                                group_names = NULL,
                                subgenome_names = NULL,
                                seed_expmx = NULL) {
    # annotations
    gene_names <- paste0("gene_", seq_len(n_genes))
    if (is.null(group_names))
        group_names <- c(paste0('group_', seq_len(length(n_replicates))))
    if (is.null(subgenome_names))
        subgenome_names <- LETTERS[seq_len(n_subgenomes)]
    # prepare mean and dispersion
    seed_expmx <- .init_seed_expmx(seed_expmx, n_subgenomes)
    seed_params <- .get_seed_params(n_genes, seed_expmx)
    nls_object <- .nls_mean_dispersion(seed_params$population)
    # initialize means and ratios of homeolog expression (inducing DEG and BEH)
    mh <- list(m = .sim_counts_DEG(gene_names, subgenome_names, group_names,
                                   seed_params),
               h = .sim_counts_HER(gene_names, subgenome_names, group_names))
    # adjust to accommodate each other and add some unexpressed homeologs
    mh <- .sim_counts_rmOutliers(mh)
    mh <- .sim_counts_shuffleParams(mh)
    mh <- .sim_counts_addZeros(mh)
    # generate artificial counts
    counts_obj <- .sim_counts(mh, n_replicates, seed_params, nls_object)
    mh <- counts_obj$mh
    counts <- counts_obj$exp
    # clean up
    names(counts) <- subgenome_names
    for (i in seq_along(counts)) {
        rownames(counts[[i]]) <- gene_names
        colnames(counts[[i]]) <- paste(rep(group_names, times = n_replicates),
                                       sequence(n_replicates), sep = '__')
    }
    new("ExpMX",
        data = counts,
        gene_names = gene_names,
        exp_design = data.frame(group = rep(group_names, times = n_replicates),
                                replicate = sequence(n_replicates)),
        meta = new("SimParams",
                   n_subgenomes = n_subgenomes,
                   n_genes = n_genes,
                   n_groups = length(n_replicates),
                   n_replicates = n_replicates,
                   params_population = seed_params$population,
                   params = seed_params$sample,
                   nls = nls_object,
                   mu = mh$m,
                   her = mh$h))
}


#' Define Ground-truth Homeologs with Differential Expression Ratios
#'
#' Defines ground-truth homeologs with differential expression ratios among
#' subgenomes across multiple conditions in allopolyploid species,
#' from simulated count data stored in an \linkS4class{ExpMX} object.
#'
#' An \linkS4class{ExpMX} object generated by `sim_homeolog_counts()`
#' contains the true homeolog expression ratios (HERs), derived from the mean
#' expression values used to simulate RNA-Seq read counts under a negative
#' binomial distribution.
#' This function calculates two summary statistics from these HERs:
#' **Dmax**, the maximum absolute difference in HERs between two conditions,
#' and **ORmax**, the maximum odds ratio of HERs across all subgenomes.
#' Homeologs are classified as significantly HER-shifting if their
#' Dmax and/or ORmax values exceed the specified thresholds.
#'
#' @param x An object of class \linkS4class{ExpMX}.
#' @param base Either `0`, a positive integer, or a character string.
#'      Values to specify the subgenome to use as the reference for calculating
#'      Dmax and ORmax. If `0`, all subgenomes are used in the calculation.
#' @param groups Vector. A character vector of length two specifying the
#'      condition groups to compare.
#' @param Dmax Numeric. Homeologs with Dmax values greater than this
#'      threshold are flagged as HER-shifting.
#' @param ORmax Numeric. Homeologs with ORmax values greater than this
#'      threshold are flagged as HER-shifting.
#' @param operator Character. Either `'AND'` or `'OR'`, specifying how
#'      the `Dmax` and `ORmax` criteria are combined to ground-truth homeologs
#'      with differential expression ratios between conditions.
#'
#' @return A logical vector indicating whether each homeolog is classified
#'     as a ground-truth homeolog with differential expression ratios
#'     that satisfy the given thresholds (`TRUE`) or not (`FALSE`).
#' 
#' @examples
#' x <- sim_homeolog_counts()
#' is_sig <- def_sigShift(x)
#' table(is_sig)
#'
#' @seealso [sim_homeolog_counts()], [def_sigShift()]
#' @export
def_sigShift <- function(x,
                         base = 0, groups = NULL,
                         Dmax = 0.15, ORmax = 2, operator = c('OR', 'AND')) {
    operator <- match.arg(operator)
    ormx <- dmx <- matrix(0, nrow = nrow(x@meta@her[[1]]),
                             ncol = ncol(x@meta@her[[1]]))
    
    if (!is.null(groups)) {
        # find BEHs in the target groups
        dmx <- .calc_exp_dist(x@meta@her[[groups[1]]], x@meta@her[[groups[2]]])
        ormx <- .calc_exp_oddsratio(x@meta@her[[groups[1]]],
                                    x@meta@her[[groups[2]]])
    } else {
        # find BEHs among whole dataset
        for (g1 in seq_along(x@meta@her)) {
            for (g2 in seq_along(x@meta@her)) {
                .dmx <- .calc_exp_dist(x@meta@her[[g1]], x@meta@her[[g2]])
                .ormx <- .calc_exp_oddsratio(x@meta@her[[g1]], x@meta@her[[g2]])
                for (i in seq_len(ncol(x@meta@her[[1]]))) {
                    dmx[, i] <- apply(cbind(dmx[, i], .dmx[, i]), 1, max,
                                      na.rm = TRUE)
                    ormx[, i] <- apply(cbind(ormx[, i], .ormx[, i]), 1, max,
                                       na.rm = TRUE)
                }
            }
        }
    }
    
    if (base > 0) {
        gt_d <- (dmx[, base] > Dmax)
        gt_or <- (ormx[, base] > ORmax)
    } else {
        gt_d <- apply(dmx > Dmax, 1, any)
        gt_or <- apply(ormx > ORmax, 1, any)
    }
    
    if (operator == 'OR') {
        gt <- (gt_d | gt_or)
    } else if (operator == 'AND') {
        gt <- (gt_d & gt_or)
    }
    
    names(gt) <- x@gene_names
    gt
}



#' Retrieve Simulation Parameters
#'
#' Retrieves the mean and dispersion values used to simulate read count data
#' from a negative binomial distribution.
#'
#' The `sim_homeolog_counts()` function generates artificial read counts
#' for each homeolog using a negative binomial distribution.
#' The mean and dispersion values used for this simulation are derived
#' from real RNA-Seq datasets and stored in a \linkS4class{SimParams} object,
#' which is embedded within the returned \linkS4class{ExpMX} object.
#'
#' This function provides an interface to extract specific simulation parameters
#' from the \linkS4class{SimParams} object.
#'
#' @param x An \linkS4class{ExpMX} class object
#'      generated by `sim_homeolog_counts()`.
#' @param param Character. A character string specifying which parameter
#'      to retrieve. Valid options include `"mu"` (true means),
#'      `"her"` (homeolog expression ratios),
#'      `"params"` (sampled means and dispersions),
#'      and `"params_population"` (population-level values).
#'
#' @return A list of data.frame objects
#'      containing the requested parameter values.
#'
#' @examples
#' x <- sim_homeolog_counts()
#' param_mu <- get_sim_params(x, 'mu')
#' head(param_mu)
#'
#' @seealso [sim_homeolog_counts()], \linkS4class{ExpMX},
#'     \linkS4class{SimParams}
#' 
#' @export
get_sim_params <- function(x, param = NULL) {
    if (is.null(param)) return(NULL)
    
    output <- NULL
    
    if (param == 'mu' || param == 'mean') {
        output <- x@meta@mu
    } else if (param == 'disp' || param == 'dispersion') {
        output <- NULL
    } else if (param == 'HER' || param == 'her') {
        output <- x@meta@her
    }
    
    for (i in seq_along(output)) {
        output[[i]] <- data.frame(output[[i]])
    }
    
    output
}

