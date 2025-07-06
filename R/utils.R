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


# Calculate gene expression by aggregating all homeolog expression
#' @importFrom methods is
.calc_gexp <- function(x) {
    if (!is(x, 'ExpMX'))
        stop('The input data should be stored as ExpMX class.')
    
    gexp <- 0
    for (i in seq_along(x@data)) {
        gexp <- gexp + x@data[[i]]
    }
    gexp
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
        hexp_ratios[, i] <- rowMeans(as.matrix(x[[i]] / gene_exp), na.rm = TRUE)
    }

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
    
    ormax <- matrix(NA, ncol = ncol(exp1), nrow = nrow(exp1))
    or_12 <- (exp1 / (1 - exp1)) / (exp2 / (1 - exp2))
    or_21 <- (exp2 / (1 - exp2)) / (exp1 / (1 - exp1))
    for (i in seq_len(ncol(exp1))) {
        ormax[, i] <- apply(cbind(or_12[, i], or_21[, i]), 1, max, na.rm = TRUE)
    }
    ormax
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

