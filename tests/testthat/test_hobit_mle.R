local_hobit_mle_2x_input <- function() {
    set.seed(1)
    sim_homeolog_counts(n_genes = 100, n_replicates = c(3, 3), n_subgenomes = 2)
}

local_hobit_mle_3x_input <- function() {
    set.seed(1)
    sim_homeolog_counts(n_genes = 10, n_replicates = c(3, 3, 3), n_subgenomes = 3)
}

local_expect_mle_output <- function(result, x) {
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), length(x@gene_names))
    expect_equal(result$gene, x@gene_names)
    expect_named(result, c(
        "gene", "pvalue", "qvalue", "Dmax", "LRmax", "RRmax",
        paste0("theta0__", names(x@data))
    ))
    expect_type(result$pvalue, "double")
    expect_type(result$qvalue, "double")
}

local_expect_mle_debug_output <- function(result, x, required_columns) {
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), length(x@gene_names))
    expect_equal(result$gene, x@gene_names)
    expect_true(all(required_columns %in% names(result)))
}

test_that("hobit_mle() returns end-user MLE results for 2x", {
    x <- local_hobit_mle_2x_input()

    result <- suppressWarnings(hobit_mle(x, n_threads = 1, iter.max = 20, eval.max = 20))

    local_expect_mle_output(result, x)
})

test_that("hobit_mle(.debug = TRUE) returns detailed MLE results for 3x", {
    x <- local_hobit_mle_3x_input()

    result <- suppressWarnings(hobit_mle(x, n_threads = 1, iter.max = 20, eval.max = 20, .debug = TRUE))

    local_expect_mle_debug_output(result, x, c(
        "lrt",
        "df",
        "status",
        "converged_H0",
        "converged_H1",
        "theta0__C",
        "theta1__C__group_3",
        "D__C__(group_2-group_3)",
        "LR__B/C__(group_2-group_3)",
        "RR__B/C__(group_2/group_3)"
    ))
    expect_true(all(result$df == 4))
})
