test_that('Test ExpMX.', {
    gexp <- read.table(system.file(package = 'hespresso', 'extdata', 'C_flexuosa.tsv.gz'),
                       header = TRUE, sep = '\t', row.names = 1)
    group <- c('wet', 'wet', 'wet', 'dry', 'dry', 'dry')
    hnames <- read.table(system.file(package = 'hespresso', 'extdata', 'C_flexuosa.homeolog.tsv.gz'),
                         header = TRUE, sep = '\t')
    
    expmx <- newExpMX(gexp, group, hnames)
    show(expmx)
})


test_that('Test ExpMX combine.', {
    x <- sim_homeolog_counts(n_subgenomes = 3)
    x_combined <- combine_hexp(x, subgenomes = c(1, 2))
    
    x_combined <- combine_hexp(x, subgenomes = c(1, 2), name_to = '1+2')
})


test_that('Test ExpMX subset', {
    x <- sim_homeolog_counts(n_genes = 100)
    
    boolvec <- rep(FALSE, length = 100)
    idx <- sample(seq_along(boolvec), 10, FALSE)
    boolvec[idx] <- TRUE
        
    x1 <- x[idx, ]
    x2 <- x[boolvec, ]
    expect_setequal(x1@data[[1]], x2@data[[1]])
    expect_setequal(x1@data[[2]], x2@data[[2]])
    
    
    boolvec <- c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE)
    idx <- seq_along(boolvec)[boolvec]
    x1 <- x[, idx]
    x2 <- x[, boolvec]
    expect_setequal(x1@data[[1]], x2@data[[1]])
    expect_setequal(x1@data[[2]], x2@data[[2]])
    
    
    idx <- c(1, 2, 3, 101)
    expect_error(x[idx, ])
    
    idx <- c(0, 2, 3, 4)
    expect_error(x[, idx])
})

