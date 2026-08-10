#' Class to Store RNA-seq Read Count Data
#'
#' An S4 class for storing RNA-seq read count data, including gene or homeolog
#' expression, experimental design, and associated metadata.
#'
#' @slot data A list of read count matrices, with one matrix per subgenome.
#' @slot gene_names A character vector containing gene names.
#' @slot exp_design A data frame describing the experimental design.
#'      It must include a column named `group` indicating sample groupings.
#' @slot meta A list containing additional metadata related to the experiment
#'      or processing steps.
#' 
#' @return An object of class \linkS4class{SeqCountData}.
#' @seealso [newSeqCountData()]
#' @exportClass SeqCountData
setClass("SeqCountData",
         slots = c(
             data = "list",
             gene_names = "character",
             exp_design = "data.frame",
             meta = "ANY"
         ))

setValidity("SeqCountData", function(object) {
    if (!('group' %in% colnames(object@exp_design))) {
        stop('The data.frame `exp_design` should contain a column named `group`.',
             call. = FALSE)
    }

    if (any(grepl("__\\.\\.__", object@exp_design$group))) {
        stop("The group name can not contain `__..__` as the `__..__` is reserved for treating data labels during the analysis progress.",
             call. = FALSE)
    }

    if (length(unique(object@exp_design$group)) < 2) {
        stop("The number of group should be equal or larger than 2.",
             call. = FALSE)
    }

    if (anyNA(object@exp_design$group)) {
        stop("The `group` column must not contain NA.", call. = FALSE)
    }

    if (length(object@data) < 2) {
        stop("The number of subgenome should be equal or larger than 2.")
    }

    ref_dim <- dim(object@data[[1L]])
    valid_dims <- vapply(object@data, function(m) identical(dim(m), ref_dim), logical(1L))
    if (!all(valid_dims)) {
        stop("All subgenome count matrices must have identical dimensions.", call. = FALSE)
    }

    if (ncol(object@data[[1L]]) != nrow(object@exp_design)) {
        stop("The number of samples in `data` must match the number of rows in `exp_design`.",
        call. = FALSE)
    }

    TRUE
})


#' Print Contents of an `SeqCountData` Object
#'
#' Prints the contents of an `SeqCountData` object to the console.
#'
#' @param object An \linkS4class{SeqCountData} object.
#'
#' @return The input \linkS4class{SeqCountData} object, returned invisibly.
#'
#' @examples
#' x <- sim_homeolog_counts(100)
#' x
#' 
#' @seealso \linkS4class{SeqCountData}
#' @importFrom utils head
#' @importFrom methods show
#' @export
setMethod(
    f = "show",
    signature = "SeqCountData",
    definition = function(object) {
        msg <- c('# ', length(object@data), ' subgenome sets', ' (', paste(names(object@data), collapse = ', '),')', '\n',
                 '# ', nrow(object@data[[1]]), ' homeolog tuples', '\n')
        cat(paste(msg, collapse = ''))
        cat('---------------------\n')
        cat('Experiment Design:\n')
        print(object@exp_design)
        cat('---------------------\n')
        for (i in seq_along(object@data)) {
            cat('> subgenome:', names(object@data)[i], '\n')
            print(head(object@data[[i]]))
            if (i != length(object@data)) {
                cat('+++++++++++++++++++++\n')
            } else {
                cat('---------------------\n')
            }
        }
        invisible(object)
})


#' Extract a Subset of an SeqCountData Object
#'
#' Extracts homeolog tuples or samples from an `SeqCountData` object.
#'
#' @name subset-SeqCountData
#' @rdname subset-SeqCountData
#' @docType methods
#' @aliases [,SeqCountData-method
#' @aliases [,SeqCountData,ANY,ANY,ANY-method
#'
#' @param x An \linkS4class{SeqCountData} object.
#' @param i A logical or integer vector specifying the homeolog tuples
#'      to retain.
#' @param j A logical or integer vector specifying the samples to retain.
#' @param ... Unused. Included for consistency with the `[` generic.
#' @param drop Unused. Included for consistency with the `[` generic.
#'
#' @return A subsetted \linkS4class{SeqCountData} object.
#'
#' @examples
#' x <- sim_homeolog_counts(n_genes = 100)
#' x_10 <- x[seq_len(10), ]
#'
#' @seealso \linkS4class{SeqCountData}
#' @export
setMethod(
    f = "[",
    signature = signature(
        x = "SeqCountData",
        i = "ANY",
        j = "ANY",
        drop = "ANY"
    ),
    definition = function(x, i, j, ..., drop) {
        if (missing(i)) i <- NULL
        if (missing(j)) j <- NULL
        i <- .int2logicalvec(i, nrow(x@data[[1]]))
        j <- .int2logicalvec(j, ncol(x@data[[1]]))
        for (s in seq_along(x@data)) {
            x@data[[s]] <- x@data[[s]][i, j, drop = FALSE]
        }
        x@gene_names <- x@gene_names[i]
        x@exp_design <- x@exp_design[j, , drop = FALSE]
        x
    }
)


#' Combine Homeolog Expression Across Subgenomes
#'
#' Sums homeolog expression values from multiple subgenomes into a single
#' subgenome.
#'
#' An `SeqCountData` object stores homeolog expression data for each subgenome as a
#' list of matrices in the `data` slot. This function aggregates selected
#' subgenome expression matrices into one combined matrix and updates the
#' relevant fields in the `SeqCountData` object accordingly.
#'
#' This function is particularly useful for analyzing homeolog expression in
#' species like wheat, which has three subgenomes (commonly named A, B, and D).
#' For example, to analyze the expression ratio between the combined expression
#' of subgenomes A and B versus D, users can apply this function to merge
#' A and B into a single subgenome.
#'
#' @param x An \linkS4class{SeqCountData} object.
#' @param subgenomes A vector of indices specifying which subgenomes to combine.
#' @param name_to A character string specifying the name of the new combined
#'      subgenome. If not provided, the name will be generated by concatenating
#'      the names of the combined subgenomes.
#'
#' @return An updated \linkS4class{SeqCountData} object with the specified subgenomes
#'      combined.
#'
#' @examples
#' x <- sim_homeolog_counts(n_genes = 100, n_subgenomes = 3)
#' x_combined <- combine_hexp(x, subgenomes = c(1, 2))
#' x_combined
#'
#' @export
combine_hexp <-function(x, subgenomes, name_to = NULL) {
    combined_data <- vector('list', length = length(x@data) - length(subgenomes) + 1)
    combined_name <- rep('', length = length(combined_data))
    i <- 1
    last_i <- length(combined_data)
    for (s in seq_along(x@data)) {
        if (s %in% subgenomes) {
            if (is.null(combined_data[[last_i]])) {
                combined_data[[last_i]] <- x@data[[s]]
            } else {
                combined_data[[last_i]] <- combined_data[[last_i]] + x@data[[s]]
            }
        } else {
            combined_data[[i]] <- x@data[[s]]
            combined_name[i] <- names(x@data)[s]
            i <- i + 1
        }
    }
    if (is.null(name_to)) {
        combined_name[last_i] <- paste(names(x@data)[subgenomes], collapse = '+')
    } else {
        combined_name[last_i] <- name_to
    }

    names(combined_data) <- combined_name
    x@data <- combined_data
    x
}


#' Load Expression Dataset
#'
#' Organize gene expression data into homeolog expression matrices and store
#' them as an object of class `SeqCountData`.
#'
#' This function arranges the given gene expression matrix
#' into homeolog expression matrices using a mapping table (`mapping_table`),
#' and stores the result as an \linkS4class{SeqCountData} object for downstream analyses.
#'
#' @param x A matrix or data frame of gene expression data,
#'      where each column corresponds to a sample
#'      and each row corresponds to a gene.
#'      The row names should be gene names that match those in the `mapping_table`.
#' @param group A vector or a data frame describing the experimental design.
#'      If a data frame, it must include a column named `group`
#'      indicating sample groupings.
#' @param mapping_table A data frame
#'      representing the mapping table of homeolog names.
#'      Each column should correspond to a progenitor genome.
#'  
#' @return An \linkS4class{SeqCountData} object containing homeolog expression data.
#'
#' @examples
#' gexp <- read.table(system.file(package = 'hespresso', 'extdata', 'C_flexuosa.tsv.gz'),
#'                    header = TRUE, sep = '\t', row.names = 1)
#' group <- c('wet', 'wet', 'wet', 'dry', 'dry', 'dry')
#' hnames <- read.table(system.file(package = 'hespresso', 'extdata', 'C_flexuosa.homeolog.tsv.gz'),
#'                      header = TRUE, sep = '\t')
#' expmx <- newSeqCountData(gexp, group, hnames)
#'
#' @seealso \linkS4class{SeqCountData}
#' @export
newSeqCountData <- function(x, group, mapping_table) {
    if (is.vector(group)) {
        group <- data.frame(group = group)
    } else {
        if (!('group' %in% colnames(group))) {
            stop('The experimental design must has a column named "group".')
        }
    }
    if (!is.factor(group$group)) {
        message('The `group` has been converted into factors.')
        group$group <- factor(group$group, levels = unique(group$group))
    }

    n_subgenomes <- ncol(mapping_table)
    x_counts <- vector('list', length = n_subgenomes)
    names(x_counts) <- colnames(mapping_table)

    # remove undef genes
    undef_genes <- setdiff(rownames(x), unlist(mapping_table))
    if (length(undef_genes) > 0) {
        warning(length(undef_genes),
                ' genes in expression data are not defined in mapping table. They will be ignored.',
                '\n',
                'Undefined genes: ',
                paste(head(undef_genes, 10), collapse = ', '),
                ifelse(length(undef_genes) > 10, ', ...', ''),
                '\n')
        x <- x[setdiff(rownames(x), undef_genes), , drop = FALSE]
    }

    # find homeolog tuples in expression data
    valid_homeologs <- NULL
    for ( i in seq_len(n_subgenomes)) {
        ids <- seq(1, nrow(mapping_table))
        if (is.null(valid_homeologs)) {
            valid_homeologs <- ids[(mapping_table[, i] %in% rownames(x))]
        } else {
            valid_homeologs <- intersect(valid_homeologs, ids[(mapping_table[, i] %in% rownames(x))])
        }
    }
    valid_homeologs <- sort(valid_homeologs)
    unvalid_homeologs <- mapping_table[!(seq_len(nrow(mapping_table)) %in% valid_homeologs), ]
    if (nrow(unvalid_homeologs) > 0) {
        warning(nrow(unvalid_homeologs),
                ' homeolog tuples are only defined in one of subgenomes. They will be ignored.',
                '\n',
                'Examples of ignored tuples:\n',
                paste(apply(head(unvalid_homeologs, 10), 1, paste, collapse = ', '), collapse = '\n'),
                ifelse(nrow(unvalid_homeologs) > 10, '\n...', ''),
                '\n')
    }

    # format
    for (i in seq_len(n_subgenomes)) {
        x_counts[[i]] <- as.matrix(x[mapping_table[valid_homeologs, i], ])
        rownames(x_counts[[i]]) <- NULL
        colnames(x_counts[[i]]) <- NULL
    }

    new("SeqCountData",
        data = x_counts,
        gene_names = mapping_table[valid_homeologs, 1],
        exp_design = group,
        meta = NULL)
}


