set.seed(1)
x <- sim_homeolog_counts(n_genes = 100, n_subgenomes = 2)
is_sig <- def_sigShift(x)


test_that('Test HER visualization function.', {
    plot_HER_distr(x)
    
    plot_HER(x)
    
    plot_HER(x, label = is_sig)
    
    plot_HER(x, groups = unique(x@exp_design$group))
})
