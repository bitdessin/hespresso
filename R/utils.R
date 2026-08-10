# Generate a matrix with the given row and column names.
.init_matrix <- function(x = NA, row_names = NULL, col_names = NULL) {
    if (is.null(row_names ) || is.null(col_names))
        stop('The row and column names are required to generaet a matrix.')
    x <- matrix(x, nrow = length(row_names), ncol = length(col_names))
    rownames(x) <- row_names
    colnames(x) <- col_names
    x
}


# Set default group name
#
# @param groups A vector of characters indicating group names.
# @param sel_groups NULL or a vector containing two elements indicating the
#        selected group names. If NULL is given, the first two group names
#        will be selected.
# @return A vector containing the two selected group names.
.set_default_groups <- function(group_names = NULL, sel_group_names = NULL) {
    if (is.null(group_names))
        stop('The candidate group names should be given. ')
    
    if (is.null(sel_group_names)) {
        sel_group_names <- unique(group_names)[seq_len(2)]
        message('The group names were set to "',
                sel_group_names[1], '" and "', sel_group_names[2], '".')
    } else {
        if (length(sel_group_names) != 2)
            stop('The number of selected group name should be 2, but ',
                 length(sel_group_names), ' was given.')
    }
    sel_group_names
}


# Calculate homeolog tuple-wise expression by summing up all homeolog expression
#' @importFrom methods is
.calc_gexp <- function(x) {
    if (!is(x, 'SeqCountData'))
        stop('The input data should be stored as SeqCountData class.')

    h_exp <- 0
    for (i in seq_along(x@data)) {
        h_exp <- h_exp + x@data[[i]]
    }
    h_exp
}
.calc_hexp <- function(x) {
    .calc_gexp(x)
}

# Calculate homeolog expression ratios from expression matrix
#
# @param x A list with each element of homeolog expression matrix, where
#          each row of matrix indicates homeolog and each column represents
#          biological replicates. 
# @return a matrix consists of homeolog expression ratios.
.calc_hexp_ratios <- function(x) {
    n_subgenomes <- length(x)

    gene_exp <- 0
    for (i in seq_len(n_subgenomes)) {
        gene_exp <- gene_exp + x[[i]]
    }
    
    hexp_ratios <- matrix(NA, ncol = n_subgenomes, nrow = nrow(gene_exp))
    for (i in seq_len(n_subgenomes)) {
        ratios <- as.matrix(x[[i]] / gene_exp)
        ratios[!is.finite(gene_exp) | gene_exp == 0] <- NA_real_
        hexp_ratios[, i] <- rowMeans(ratios, na.rm = TRUE)
    }
    hexp_ratios[is.nan(hexp_ratios)] <- NA_real_

    colnames(hexp_ratios) <- names(x)
    if (!is.null(rownames(x[[1]])))
        rownames(hexp_ratios) <- rownames(x[[1]])
    
    hexp_ratios
}



# Calculate absolute distances of expression between two subgenomes
#
# @param exp1 An expression matrix of the first subgenome.
# @param exp2 An expression matrix of the second subgenome. The number of rows
#               and columns should equal to that of `exp1`.
.calc_exp_dist <-function(exp1, exp2) {
    if ((nrow(exp1) != nrow(exp2)) || (ncol(exp1) != ncol(exp2)))
        stop('The dimensions of the two provided expression matrices',
             'do not match.')
    abs(exp1 - exp2)
}


# Calculate odds-ratio of expression between two subgenomes
#
# @param exp1 An expression matrix of the first subgenome.
# @param exp2 An expression matrix of the second subgenome. The number of rows
#               and columns should equal to that of `exp1`.
.calc_exp_oddsratio <-function(exp1, exp2) {
    if ((nrow(exp1) != nrow(exp2)) || (ncol(exp1) != ncol(exp2)))
        stop('The dimensions of the two provided expression matrices',
             'do not match.')
    if (any(!is.finite(exp1)) || any(!is.finite(exp2)) ||
        any(exp1 < 0 | exp1 > 1) || any(exp2 < 0 | exp2 > 1)) {
        stop('The provided expression ratios must be finite values between 0 and 1.', call. = FALSE)
    }

    log_odds1 <- log(exp1) - log1p(-exp1)
    log_odds2 <- log(exp2) - log1p(-exp2)
    log_or <- abs(log_odds1 - log_odds2)
    log_or[is.nan(log_or)] <- 0
    exp(log_or)
}


# Convert indexes to a boolean vector
#
# @param i A vector of indexes.
# @param n A integer to specify the length of the boolean vector.
.int2logicalvec <- function(i, n) {
    if (is.null(i)) {
        boolvec <- rep(TRUE, length = n)
    } else {
        if (is.logical(i)) {
            boolvec <- i
        }
        if (is.numeric(i)) {
            if ((max(i) > n) || (min(i) < 1))
                stop('The given indexes out of range.')
            boolvec <- rep(FALSE, length = n)
            boolvec[i] <- TRUE
        }
    }
    boolvec
}


.as_positive_int <- function(x, arg_name, min_threshold = NULL) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 ||
        x > .Machine$integer.max || x != floor(x)) {
        stop(paste0("The `", arg_name, "` must be one positive integer."), call. = FALSE)
    }
    if (!is.null(min_threshold)) {
        if (!is.numeric(min_threshold) || length(min_threshold) != 1L ||
            !is.finite(min_threshold) || min_threshold < 1 ||
            min_threshold > .Machine$integer.max || min_threshold != floor(min_threshold)) {
            stop("`min_threshold` must be one positive integer.", call. = FALSE)
        }
        if (x < min_threshold) {
            stop(paste0("The `", arg_name, "` must be at least ", min_threshold, "."), call. = FALSE)
        }
    }

    as.integer(x)
}


.as_positive_float <- function(x, arg_name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
        stop(paste0("`", arg_name, "` must be one finite positive number."), call. = FALSE)
    }
    as.numeric(x)
}

.as_nonnegative_float <- function(x, arg_name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
        stop(paste0("`", arg_name, "` must be one non-negative number."), call. = FALSE)
    }
    as.numeric(x)
}

.as_prob <- function(x, arg_name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
        stop(paste0("`", arg_name, "` must be one finite number between 0 and 1."), call. = FALSE)
    }
    as.numeric(x)
}
