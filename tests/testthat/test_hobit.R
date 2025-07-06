set.seed(1)
x2 <- sim_homeolog_counts(n_genes = 10, n_subgenomes = 2)
x3 <- sim_homeolog_counts(n_genes = 10, n_subgenomes = 3)
x2nr <- sim_homeolog_counts(n_genes = 10, n_replicates = c(1, 1, 1, 1, 1), n_subgenomes = 3)

test_that('Test HOBIT for 4x.', {
    x_output <- hobit(x2, chains = 1, iter_warmup = 100, iter_sampling = 100)
})

test_that('Test HOBIT for 4x with multiple threads.', {
    x11 <- hobit(x2, n_threads = 1, parallel_chains = 1,
                 chains = 4, iter_warmup = 1000, iter_sampling = 2000)
    
    x14 <- hobit(x2, n_threads = 1, parallel_chains = 4,
                 chains = 4, iter_warmup = 1000, iter_sampling = 2000)
    
    x41 <- hobit(x2, n_threads = 4, parallel_chains = 1,
                 chains = 4, iter_warmup = 1000, iter_sampling = 2000)
    
    d1 <- x11$pvalue - x14$pvalue
    d2 <- x11$pvalue - x41$pvalue
    
    expect_setequal(d1 < 0.05, rep(TRUE, length(d1)))
    expect_setequal(d2 < 0.05, rep(TRUE, length(d2)))
})

test_that('Test HOBIT for 6x.', {
    x_output <- hobit(x3, chains = 1, iter_warmup = 100, iter_sampling = 100)
})

test_that('Test HOBIT for 4x using ZINB.', {
    x_output <- hobit(x2,
                      use_prior = TRUE,
                      dist = 'ZINB',
                      chains = 1, iter_warmup = 100, iter_sampling = 100)
})

test_that('Test HOBIT for 4x without replicates.', {
    x_output <- hobit(x2nr, no_replicate = TRUE,
                      chains = 1, iter_warmup = 100, iter_sampling = 100)
})
