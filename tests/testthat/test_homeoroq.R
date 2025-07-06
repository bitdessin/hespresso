set.seed(1)
x <- sim_homeolog_counts(n_genes = 1000, n_subgenomes = 2)

test_that('Test HomeoRoq.', {
    x_output <- homeoroq(x, 2, 100)
})

test_that('Test HomeoRoq with multiple threads.', {
    x_output <- homeoroq(x, 4, 100, n_threads = 4)
})

test_that('Test HomeoRoq with unexpected conditions.', {
    set.seed(1)
    x <- sim_homeolog_counts(n_genes = 10, n_subgenomes = 2)
    expect_error(homeoroq(x))
    
    x <- sim_homeolog_counts(n_genes = 100, n_subgenomes = 3)
    expect_error(homeoroq(x))
})



