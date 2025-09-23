#' Visualize distribution of homeolog expression ratios
#' 
#' Visualizes distribution of homeolog expression ratio using histogram.
#'  
#' @param x An \linkS4class{ExpMX} class object.
#' @param base Integer or character. An index or name specifying the subgenome
#'      used to calculate homeolog expression ratios.
#' @return A list of ggplot histogram objects, one for each group.
#' 
#' @examples
#' x <- sim_homeolog_counts(100)
#' plot_HER_distr(x)
#' 
#' @importFrom ggplot2 ggplot geom_point geom_histogram aes_string
#' @export
plot_HER_distr <- function(x, base = 1) {
    n_subgenomes <- length(x@data)
    
    her_figs <- vector('list', length = length(unique(x@exp_design$group)))
    names(her_figs) <- unique(x@exp_design$group)
    
    for (group_name in unique(x@exp_design$group)) {
        x_counts <- lapply(x@data, function(y, cl) {
            as.matrix(y[, cl])
        }, cl = (x@exp_design$group == group_name))
        hexp_ratios <- .calc_hexp_ratios(x_counts)
        
        if (base > 0) {
            hexp_ratios <- data.frame(HER = hexp_ratios[, base])
            her_figs[[group_name]] <- ggplot(data = hexp_ratios, 
                                             aes_string(x = 'HER')) +
                geom_histogram()
        }
    }
    her_figs
}




#' Visualize Homeolog Expression Ratios between Two Groups
#' 
#' Creates a scatter plot to visualize homeolog expression ratios between two
#' groups.
#'  
#' @param x An \linkS4class{ExpMX} class object.
#' @param base Integer or character. An index or name specifying the subgenome
#'      used for calculating homeolog expression ratios.
#' @param groups Character vector of length two specifying the group names
#'      to compare.
#' @param label Optional character vector or factor indicating categories
#'      or types of each homeolog. Points will be colored based on this grouping
#'      if provided.
#' @param size Numeric. A value specifying the point size in the scatter plot.
#' @param alpha Numeric. A decimal number between 0 and 1 specifying the
#'      transparency level of points in the scatter plot.
#' 
#' @return A ggplot scatter chart object.
#' 
#' @examples
#' x <- sim_homeolog_counts(100)
#' plot_HER(x)
#' 
#' @importFrom ggplot2 ggplot geom_point geom_hline geom_vline geom_abline
#' @importFrom ggplot2 aes_string xlab ylab scale_colour_discrete
#' @export
plot_HER <- function(x,
                     base = 1, groups = NULL, label = NULL,
                     size = 3, alpha = 0.8) {
    groups <- .set_default_groups(x@exp_design$group, groups)
    
    hexp_1 <- lapply(x@data, function(y, cl) {y[, cl]},
                     cl = (x@exp_design$group == groups[1]))
    hexp_2 <- lapply(x@data, function(y, cl) {y[, cl]},
                     cl = (x@exp_design$group == groups[2]))
    her_1 <- .calc_hexp_ratios(hexp_1)
    her_2 <- .calc_hexp_ratios(hexp_2)
    
    her_df <- data.frame(x = her_1[, base], y = her_2[, base])
    if (!is.null(label)) her_df$label <- label
    
    g <- ggplot()
    if (is.null(label)) {
        g <- g + geom_point(aes_string(x = 'x', y = 'y'),
                            alpha = alpha, size = size, data = her_df)
    } else {
        if (!is.factor(label)) {
            tt <- table(label)
            od <- order(tt, decreasing = TRUE)
            her_df$label <- factor(her_df$label, levels = names(tt)[od])
        }
        for (cl in levels(her_df$label)) {
            g <- g + geom_point(aes_string(x = 'x', y = 'y', color = 'label'),
                                alpha = alpha, size = size,
                                data  = her_df[her_df$label == cl, ])
        }
        g <- g + scale_colour_discrete(drop = FALSE)
    }
    
    g + geom_hline(yintercept = 0.5, linetype = "dashed") +
        geom_vline(xintercept = 0.5, linetype = "dashed") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        xlab(groups[1]) +
        ylab(groups[2])
}
