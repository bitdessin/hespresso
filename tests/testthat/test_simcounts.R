test_that('Test function to generate artificial dataset.', {
    x <- sim_homeolog_counts()
    
    x_100 <- sim_homeolog_counts(n_genes = 100)
    
    x_100_4x <- sim_homeolog_counts(n_genes = 100, n_subgenomes = 2)
    
    x_100_6x <- sim_homeolog_counts(n_genes = 100, n_subgenomes = 3)
    
    x_100_5 <- sim_homeolog_counts(n_genes = 100, n_replicates = c(3, 3), n_subgenomes = 2)
})


test_that('Generate artificial dataset using custom seed.', {
    expmx <- read.table(system.file(package = 'hespresso', 'extdata', 'seed_matrix.C_flexuosa.tsv.gz'))
    head(expmx)
    
    x <- sim_homeolog_counts(n_genes = 100, seed_expmx = expmx)
})


test_that('Test function to define DEGs and RCHs.', {
    x <- sim_homeolog_counts()
    
    is_sig_0 <- def_sigShift(x)
    
    is_sig_1 <- def_sigShift(x, base = 1)
    
    is_sig_1a <- def_sigShift(x, operator = 'AND')
})
