local_hobit_mcmc_2x_input <- function() {
    set.seed(1)
    sim_homeolog_counts(n_genes = 100, n_replicates = c(3, 3), n_subgenomes = 2)
}

local_hobit_mcmc_3x_input <- function() {
    set.seed(1)
    sim_homeolog_counts(n_genes = 10, n_replicates = c(3, 3, 3), n_subgenomes = 3)
}

local_hobit_mcmc_args <- list(
    n_threads = 1,
    parallel_chains = 1,
    chains = 1,
    iter_warmup = 10,
    iter_sampling = 10,
    refresh = 0,
    show_messages = FALSE
)

local_expect_mcmc_output <- function(result, x) {
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

local_expect_mcmc_debug_output <- function(result, x, required_columns) {
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), length(x@gene_names))
    expect_equal(result$gene, x@gene_names)
    expect_true(all(required_columns %in% names(result)))
}

local_skip_if_no_cmdstan <- function() {
    skip_if_not_installed("cmdstanr")
    skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed.")
}

test_that("hobit_mcmc(mode = 'fast') returns end-user MCMC results for 2x", {
    local_skip_if_no_cmdstan()
    x <- local_hobit_mcmc_2x_input()

    result <- suppressWarnings(do.call(
        hobit_mcmc,
        c(list(x = x, mode = 'fast'), local_hobit_mcmc_args)
    ))

    local_expect_mcmc_output(result, x)
})

test_that("hobit_mcmc() returns end-user strict MCMC results for 2x", {
    local_skip_if_no_cmdstan()
    x <- local_hobit_mcmc_2x_input()

    result <- suppressWarnings(do.call(
        hobit_mcmc,
        c(list(x = x), local_hobit_mcmc_args)
    ))

    local_expect_mcmc_output(result, x)
})

test_that("hobit_mcmc(mode = 'fast', .debug = TRUE) supports 3x", {
    local_skip_if_no_cmdstan()
    x <- local_hobit_mcmc_3x_input()

    result <- suppressWarnings(do.call(
        hobit_mcmc,
        c(list(x = x, mode = 'fast', .debug = TRUE), local_hobit_mcmc_args)
    ))

    local_expect_mcmc_debug_output(result, x, c(
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

test_that("hobit_mcmc(mode = 'strict', .debug = TRUE) supports 3x", {
    local_skip_if_no_cmdstan()
    x <- local_hobit_mcmc_3x_input()

    result <- suppressWarnings(do.call(
        hobit_mcmc,
        c(list(x = x, mode = 'strict', .debug = TRUE), local_hobit_mcmc_args)
    ))

    local_expect_mcmc_debug_output(result, x, c(
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
