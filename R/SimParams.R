#' Simulation Parameters and Ground Truth
#'
#' An S4 class that stores the design, generating parameters, and ground truth
#' for RNA-seq count data simulated by [sim_homeolog_counts()]. The object is
#' stored in the `meta` slot of the returned `SeqCountData` object.
#'
#' @slot n_subgenomes Number of simulated subgenomes.
#' @slot n_genes Number of simulated homeolog tuples.
#' @slot n_groups Number of simulated conditions.
#' @slot n_replicates Integer vector giving the number of replicates per
#'   condition.
#' @slot params_population Data frame containing the empirical mean, variance,
#'   and dispersion population used to define the simulation model.
#' @slot params Data frame containing the sampled baseline parameters for the
#'   simulated homeolog tuples.
#' @slot nls Fitted mean-dispersion trend model used to generate dispersions.
#' @slot mu Named list of condition-specific expected mean-expression matrices
#'   for homeologs before library-size offsets are applied.
#' @slot her Named list of condition-specific true homeolog expression-ratio
#'   matrices.
#' @slot total_mu Named list of condition-specific total-expression means.
#' @slot dispersion Named list of per-gene negative-binomial dispersions, one
#'   vector for each subgenome.
#' @slot log_offset Centered log library-size offsets, one per simulated sample.
#' @slot is_shift Logical vector indicating genes with exact
#'   homeolog-expression-ratio shifts.
#' @slot is_deg Logical vector indicating genes with an independent
#'   total-expression change.
#' @slot shift_group Character vector giving the condition receiving each
#'   homeolog-expression-ratio shift.
#' @slot deg_group Character vector giving the condition receiving each
#'   total-expression change.
#' @slot log2fc Numeric vector of total-expression log2 fold changes.
#' @slot shift_strength Numeric vector of homeolog-expression-ratio shift
#'   magnitudes.
#' @slot dmax_true Numeric vector of true maximum absolute ratio differences.
#' @slot lrmax_true Numeric vector of true maximum absolute pairwise log-ratio
#'   changes.
#' @slot rrmax_true Numeric vector of true maximum pairwise ratios of ratios.
#' @slot ormax_true Numeric vector of true maximum one-versus-rest odds ratios.
#' @slot settings List of simulation settings used to generate the object.
#'
#' @return An object of class `SimParams`.
#' @seealso [sim_homeolog_counts()]
#' @exportClass SimParams
setClass(
    "SimParams",
    slots = c(
        n_subgenomes = "numeric",
        n_genes = "numeric",
        n_groups = "numeric",
        n_replicates = "numeric",

        params_population = "data.frame",
        params = "data.frame",
        nls = "ANY",

        mu = "list",
        her = "list",
        total_mu = "list",
        dispersion = "list",
        log_offset = "numeric",

        is_shift = "logical",
        is_deg = "logical",
        shift_group = "character",
        deg_group = "character",
        log2fc = "numeric",
        shift_strength = "numeric",

        dmax_true = "numeric",
        lrmax_true = "numeric",
        rrmax_true = "numeric",
        ormax_true = "numeric",

        settings = "list"
    ),
    prototype = list(
        total_mu = list(),
        dispersion = list(),
        log_offset = numeric(),
        is_shift = logical(),
        is_deg = logical(),
        shift_group = character(),
        deg_group = character(),
        log2fc = numeric(),
        shift_strength = numeric(),
        dmax_true = numeric(),
        lrmax_true = numeric(),
        rrmax_true = numeric(),
        ormax_true = numeric(),
        settings = list()
    )
)

setValidity("SimParams", function(object) {
    is_count <- function(value, minimum) {
        is.numeric(value) && length(value) == 1L && is.finite(value) &&
            value >= minimum && value <= .Machine$integer.max && value == floor(value)
    }

    if (!is_count(object@n_genes, 1L)) {
        return("`n_genes` must be one positive integer.")
    }
    if (!is_count(object@n_subgenomes, 2L)) {
        return("`n_subgenomes` must be an integer of at least 2.")
    }
    if (!is_count(object@n_groups, 2L)) {
        return("`n_groups` must be an integer of at least 2.")
    }

    n_genes <- as.integer(object@n_genes)
    n_subgenomes <- as.integer(object@n_subgenomes)
    n_groups <- as.integer(object@n_groups)
    if (!is.numeric(object@n_replicates) ||
        length(object@n_replicates) != n_groups ||
        any(!is.finite(object@n_replicates)) || any(object@n_replicates < 1) ||
        any(object@n_replicates > .Machine$integer.max) ||
        any(object@n_replicates != floor(object@n_replicates))) {
        return("`n_replicates` must contain one positive integer per condition.")
    }

    if (length(object@is_shift) != n_genes ||
        length(object@is_deg) != n_genes ||
        length(object@shift_group) != n_genes ||
        length(object@deg_group) != n_genes ||
        length(object@log2fc) != n_genes ||
        length(object@shift_strength) != n_genes ||
        length(object@dmax_true) != n_genes ||
        length(object@lrmax_true) != n_genes ||
        length(object@rrmax_true) != n_genes ||
        length(object@ormax_true) != n_genes) {
        return("Gene-level ground-truth vectors must have length `n_genes`.")
    }

    if (length(object@mu) != n_groups ||
        length(object@her) != n_groups ||
        length(object@total_mu) != n_groups) {
        return("Condition-level parameter lists must have length `n_groups`.")
    }

    if (length(object@dispersion) != n_subgenomes) {
        return("`dispersion` must contain one vector per subgenome.")
    }

    if (length(object@log_offset) != sum(object@n_replicates)) {
        return("`log_offset` must contain one value per simulated sample.")
    }
    if (any(!is.finite(object@log_offset))) {
        return("`log_offset` must contain finite values.")
    }
    if (nrow(object@params) != n_genes) {
        return("`params` must contain one row per simulated gene.")
    }

    is_parameter_matrix <- function(value) {
        is.matrix(value) && is.numeric(value) && nrow(value) == n_genes &&
            ncol(value) == n_subgenomes
    }
    if (!all(vapply(object@mu, is_parameter_matrix, logical(1L)))) {
        return("`mu` must contain one numeric gene-by-subgenome matrix per condition.")
    }
    if (!all(vapply(object@her, is_parameter_matrix, logical(1L)))) {
        return("`her` must contain one numeric gene-by-subgenome matrix per condition.")
    }

    is_gene_vector <- function(value) is.numeric(value) && length(value) == n_genes
    if (!all(vapply(object@total_mu, is_gene_vector, logical(1L)))) {
        return("`total_mu` must contain one numeric vector per condition with length `n_genes`.")
    }
    if (!all(vapply(object@dispersion, is_gene_vector, logical(1L)))) {
        return("`dispersion` must contain one numeric vector per subgenome with length `n_genes`.")
    }

    TRUE
})
