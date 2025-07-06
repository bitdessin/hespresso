x <- sim_homeolog_counts()

.calc_gexp_libsize <- function(x) {
    gexp <- NULL
    for (i in seq_along(x@data)) {
        if (is.null(gexp)) {
            gexp <- x@data[[i]]
        } else {
            gexp <- gexp + x@data[[i]]
        }
    }
    colSums(gexp)
}


test_that('Test TMM.', {
    x_norm <- norm_counts(x, method = 'tmm')
    x_norm_round <- norm_counts(x, method = 'tmm', round = TRUE)
})


test_that('Test MED.', {
    x_norm <- norm_counts(x, method = 'med')
    x_norm_round <- norm_counts(x, method = 'med', round = TRUE)
})



test_that('Test CPM.', {
    x_norm <- norm_counts(x, method = 'cpm')
    
    err <- .calc_gexp_libsize(x_norm) - rep(1e6, length = ncol(x@data[[1]]))
    expect_setequal(err < 1.0, rep(TRUE, length(err)))
})

