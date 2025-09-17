rm(list = ls(all = TRUE))

library(hespresso)
options(mc.cores = 8)
options(width = 1024)
options(ggplot2.discrete.colour = c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
knitr::opts_chunk$set(tidy = FALSE,
                      cache = TRUE,
                      error = TRUE,
                      warning = FALSE,
                      message = FALSE,
                      dev = "png")
set.seed(1)


pkg <- function(pkg_name, repos = NULL) {
    pkg_html <- ''
    if (is.null(repos)) {
        if (pkg_name == 'hespresso') {
            pkg_html <- paste0('<span class="pkg-name">',
                               '<a href="https://bitdessin.github.io/hespresso" target="_blank">',
                               pkg_name, '</a></span>')
        } else if (pkg_name == 'cmdstanr') {
            pkg_html <- paste0('<span class="pkg-name">',
                               '<a href="https://mc-stan.org/cmdstanr" target="_blank">',
                               pkg_name, '</a></span>')
        } else {
            pkg_html <- paste0('<span class="pkg-name">', pkg_name, '</span>')
        }
    } else if (tolower(repos) == 'cran') {
        pkg_html <- paste0('<span class="pkg-name">',
                           '<a href="https://CRAN.R-project.org/package=',
                           pkg_name, '" target="_blank">', pkg_name, '</a></span>')
    } else if (tolower(repos) == 'bioc') {
        pkg_html <- paste0('<span class="pkg-name">',
                           '<a href="https://bioconductor.org/packages/',
                           pkg_name, '" target="_blank">', pkg_name, '</a></span>')
    }
    pkg_html
}
