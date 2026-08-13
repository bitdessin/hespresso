test_that('Test function to generate artificial dataset.', {
    x <- sim_homeolog_counts()

    x_100 <- sim_homeolog_counts(n_genes = 100)

    x_100_4x <- sim_homeolog_counts(n_genes = 100, n_subgenomes = 2)

    x_100_6x <- sim_homeolog_counts(n_genes = 100, n_subgenomes = 3)

    x_100_5r <- sim_homeolog_counts(n_genes = 100, n_replicates = c(5, 5))

    x_100_pshift <- sim_homeolog_counts(n_genes = 100, prop_shift = 0.20)

    x_100_pdeg <- sim_homeolog_counts(n_genes = 100, prop_deg = 0.20)

    x_100_pzero <- sim_homeolog_counts(n_genes = 100, prop_zero = 0)
})


test_that('Generate artificial dataset using custom seed.', {
    seed_counts <- read.table(system.file(package = 'hespresso', 'extdata', 'seed_matrix.C_flexuosa.tsv.gz'))
    x <- sim_homeolog_counts(n_genes = 100, seed_counts = seed_counts)
})


test_that('Test function to define RCHs.', {
    set.seed(1)
    x <- sim_homeolog_counts()

    expect_identical(def_sigShift(x), x@meta@is_shift)

    expect_true(any(def_sigShift(x)))

    expect_true(any(def_sigShift(x, base = 1)))

    expect_true(any(def_sigShift(x, Dmax = 0.2, ORmax = 2.0, operator = 'AND')))

    # following condition should be FALSE with Dmax>Inf or ORmax>Inf.

    expect_false(any(def_sigShift(x, Dmax = Inf)))

    expect_false(any(def_sigShift(x, ORmax = Inf)))

    expect_identical(
        def_sigShift(x, Dmax = 0),
        def_sigShift(x, Dmax = 0, ORmax = Inf, operator = 'OR')
    )

    expect_false(any(def_sigShift(x, Dmax = 0, ORmax = Inf, operator = 'AND')))
})


test_that('ratio helpers handle zero counts and boundary ratios.', {
    zero_counts <- list(A = matrix(c(0, 1), ncol = 1), B = matrix(c(0, 1), ncol = 1))
    ratios <- hespresso:::.calc_hexp_ratios(zero_counts)
    expect_true(all(is.na(ratios[1, ])))
    expect_equal(unname(ratios[2, ]), c(0.5, 0.5))

    boundary <- matrix(c(0, 1, 0.5), ncol = 1)
    expect_equal(hespresso:::.calc_exp_oddsratio(boundary, boundary), matrix(1, nrow = 3))
    expect_true(is.infinite(hespresso:::.calc_exp_oddsratio(matrix(0, 1, 1), matrix(1, 1, 1))[1, 1]))
})
