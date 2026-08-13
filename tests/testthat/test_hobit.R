local_hobit_2x_input <- function() {
    set.seed(1)
    sim_homeolog_counts(n_genes = 100, n_replicates = c(3, 3), n_subgenomes = 2)
}

local_hobit_3x_input <- function() {
    set.seed(1)
    sim_homeolog_counts(n_genes = 100, n_replicates = c(3, 3, 3), n_subgenomes = 3)
}

local_mcmc_args <- list(
    n_threads = 1,
    parallel_chains = 1,
    chains = 1,
    iter_warmup = 10,
    iter_sampling = 10,
    refresh = 0,
    show_messages = FALSE
)

local_expect_hobit_output <- function(result, x) {
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

local_expect_hobit_debug_output <- function(result, x, required_columns) {
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), length(x@gene_names))
    expect_equal(result$gene, x@gene_names)
    expect_true(all(required_columns %in% names(result)))
}

local_skip_if_no_cmdstan <- function() {
    skip_if_not_installed("cmdstanr")
    skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed.")
}

test_that("hobit() returns end-user MLE results for 2x", {
    x <- local_hobit_2x_input()

    result <- suppressWarnings(hobit(x, n_threads = 1, iter.max = 20, eval.max = 20))

    local_expect_hobit_output(result, x)
})

test_that("hobit(method = 'mle') returns end-user MLE results for 2x", {
    x <- local_hobit_2x_input()

    result <- suppressWarnings(hobit(x, method = 'mle', n_threads = 1, iter.max = 20, eval.max = 20))

    local_expect_hobit_output(result, x)
})

test_that("hobit(method = 'mcmc', mode = 'fast') returns end-user MCMC results for 2x", {
    local_skip_if_no_cmdstan()
    x <- local_hobit_2x_input()

    result <- suppressWarnings(do.call(
        hobit,
        c(list(x = x, method = 'mcmc', mode = 'fast'), local_mcmc_args)
    ))

    local_expect_hobit_output(result, x)
})

test_that("hobit(method = 'mcmc', mode = 'strict') returns end-user MCMC results for 2x", {
    local_skip_if_no_cmdstan()
    x <- local_hobit_2x_input()

    result <- suppressWarnings(do.call(
        hobit,
        c(list(x = x, method = 'mcmc', mode = 'strict'), local_mcmc_args)
    ))

    local_expect_hobit_output(result, x)
})

test_that("hobit(method = 'mcmc') rejects misspelled modes for 2x", {
    x <- local_hobit_2x_input()

    expect_error(
        hobit(x, method = 'mcmc', mode = 'strcit'),
        "'arg' should be one of"
    )
})

test_that("hobit(.debug = TRUE) supports 3x in fast MCMC", {
    local_skip_if_no_cmdstan()
    x <- local_hobit_3x_input()

    result <- suppressWarnings(do.call(
        hobit,
        c(list(x = x, method = 'mcmc', mode = 'fast', .debug = TRUE), local_mcmc_args)
    ))

    local_expect_hobit_debug_output(result, x, c(
        "lrt",
        "df",
        "status",
        "theta0__C",
        "theta1__C__group_3",
        "D__C__(group_2-group_3)",
        "LR__B/C__(group_2-group_3)",
        "RR__B/C__(group_2/group_3)"
    ))
    expect_true(all(result$df == 4))
})

test_that("hobit(.debug = TRUE) supports 3x in strict MCMC", {
    local_skip_if_no_cmdstan()
    x <- local_hobit_3x_input()

    result <- suppressWarnings(do.call(
        hobit,
        c(list(x = x, method = 'mcmc', mode = 'strict', .debug = TRUE), local_mcmc_args)
    ))

    local_expect_hobit_debug_output(result, x, c(
        "lrt",
        "df",
        "status",
        "theta0__C",
        "theta1__C__group_3",
        "D__C__(group_2-group_3)",
        "LR__B/C__(group_2-group_3)",
        "RR__B/C__(group_2/group_3)"
    ))
    expect_true(all(result$df == 4))
})

test_that("hobit(.debug = TRUE) supports 3x in MLE", {
    x <- local_hobit_3x_input()

    result <- suppressWarnings(hobit(x, method = 'mle', n_threads = 1, iter.max = 20, eval.max = 20, .debug = TRUE))

    local_expect_hobit_debug_output(result, x, c(
        "lrt",
        "df",
        "status",
        "theta0__C",
        "theta1__C__group_3",
        "D__C__(group_2-group_3)",
        "LR__B/C__(group_2-group_3)",
        "RR__B/C__(group_2/group_3)"
    ))
    expect_true(all(result$df == 4))
})
