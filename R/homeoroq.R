# Format an ExpMX Class Object to Input HomeoRoq
.hq.data.format <- function(x) {
    c1 <- x@data[[1]]
    c2 <- x@data[[2]]
    c1[c1 < c2 * 0.01] <- (c2 * 0.01)[c1 < c2 * 0.01]
    c2[c2 < c1 * 0.01] <- (c1 * 0.01)[c2 < c1 * 0.01]
    
    # set column names
    glabel <- NULL
    gdict <- list()
    for (gn in x@exp_design$group) {
        if (!(gn %in% names(gdict)))
            gdict[[gn]] <- 0
        gdict[[gn]] <- gdict[[gn]] + 1
        glabel <- c(glabel, paste0(gn, '__', gdict[[gn]]))
    }
    colnames(c1) <- paste(names(x@data)[1], glabel, sep = '__')
    colnames(c2) <- paste(names(x@data)[2], glabel, sep = '__')
    counts <- cbind(c1, c2)
    
    gn <- unique(x@exp_design$group)
    
    list(DATA = counts,
         OC = colnames(c1)[x@exp_design$group == gn[1]],
         AC = colnames(c2)[x@exp_design$group == gn[1]],
         OD = colnames(c1)[x@exp_design$group == gn[2]],
         AD = colnames(c2)[x@exp_design$group == gn[2]])
}


# Calculate Summary Statistics from Counts Data of Group 2
.hq.data.OD <- function(d) {
    counts <- d$DATA[, d$OD] + d$DATA[, d$AD]
    logmean <- rowMeans(log(counts))
    logmean_std <- exp(log(d$DATA[, d$OD]) - (log(counts) - logmean))
    v <- apply(logmean_std, 1, var) * (length(d$OD) - 1) / length(d$OD) 
    ratios <- rowMeans(d$DATA[, d$OD] / counts)
    ratios_v <- apply(d$DATA[, d$OD] / counts,
                      1, var) * (length(d$OD) - 1) / length(d$OD) 
    
    list(logmean = logmean, v = v, r = ratios, rv = ratios_v)
}


# Calculate Summary Statistics from Counts Data of Group 2
.hq.data.OC <- function(d) {
    counts <- d$DATA[, d$OC] + d$DATA[, d$AC]
    logmean <- rowMeans(log(counts))
    logmean_std <- exp(log(d$DATA[, d$OC]) - (log(counts) - logmean))
    v <- apply(logmean_std, 1, var) * (length(d$OC) - 1) / length(d$OC) 
    ratios <- rowMeans(d$DATA[, d$OC] / counts)
    ratios_v <- apply(d$DATA[, d$OC] / counts,
                      1, var) * (length(d$OC) - 1) / length(d$OC)
    
    list(logmean = logmean, v = v, r = ratios, rv = ratios_v)
}


# Calculate Summary Statistics from Counts Data of All Groups
.hq.data.OR <- function(d) {
    ORG <- c(d$OC, d$OD)
    totalCounts <- d$DATA[, ORG] + d$DATA[, c(d$AC, d$AD)]
    logmean <- rowMeans(log(totalCounts))
    ORGscaledCounts <- exp(log(d$DATA[, ORG]) - (log(totalCounts) - logmean))
    v <- apply(ORGscaledCounts, 1, var) * (length(ORG) - 1) / length(ORG) 
    ratios <- rowMeans(d$DATA[, ORG] / totalCounts)
    ORGRatioVarList <- apply(d$DATA[, ORG] / totalCounts,
                             1, var) * (length(ORG) - 1) / length(ORG)
    ORGRatioSD <- sqrt(ORGRatioVarList)
    
    list(logmean = logmean, v = v, r = ratios, ratioSD = ORGRatioSD)
}


# Calculate P-values from Counts and Given Parameters
#' @importFrom stats dnorm
.hq.calc_pval <- function(i, j, k,
                          ctrlRatio, ratioSD, ctrlExpMean, objExpMean,
                          ctrlSD, objSD) {
    stdRPval <- dnorm(i, mean = ctrlRatio, sd = ratioSD)
    ctrlPval <- dnorm(j, mean = ctrlExpMean * i, sd = ctrlSD)
    objPval <- dnorm(k, mean = objExpMean * i, sd = objSD)
    stdRPval * ctrlPval * objPval
}


# Eestimate Variance
#
# Estimate Variance with `fv` object. If `residue` is given, then correct
# the predicted values with the residues. If `obs_v` is given, then correct
# the predicted values with the observed counts.
#
#' @importFrom stats predict
.hq.estimate_var <- function(fv, logmean, i, residue = NULL, obs_v = NULL) {
    logmean[logmean < 0] <- 0
    pred_v <- predict(fv, lp(logmean, i))
    
    if ((!is.null(residue)) && (!is.null(obs_v)))
        stop('The `residue` and `obs_v` should NOT be non NULL values',
              'at the same time.')
    
    if (!is.null(residue)) {
        pred_v <- exp(pred_v) + residue
        pred_v <- vapply(seq_along(pred_v),
                         function(x) max(pred_v[x], 1.0e-5),
                         numeric(1))
    }
    
    if (!is.null(obs_v)) {
        pred_v <- exp(pred_v)
        for(i in seq_along(pred_v)) {
            obs_v[i] <- ifelse(obs_v[i] > pred_v[i] * 2, pred_v[i] * 2, 
                               ifelse(obs_v[i] < pred_v[i] / 2,
                                      pred_v[i] / 2, obs_v[i]))
        }
        pred_v <- obs_v
    }
    
    pred_v
}





# Calculate Statistics for Output from Observed Counts
.hq.calc_obscounts_stats <- function(x) {
    dOR <- .hq.data.OR(.hq.data.format(x))
    
    gn <- unique(x@exp_design$group)
    sn <- names(x@data)
    g1_sub1 <- rowSums(x@data[[1]][, x@exp_design$group == gn[1]])
    g1_sub2 <- rowSums(x@data[[2]][, x@exp_design$group == gn[1]])
    g2_sub1 <- rowSums(x@data[[1]][, x@exp_design$group == gn[2]])
    g2_sub2 <- rowSums(x@data[[2]][, x@exp_design$group == gn[2]])
    df <- data.frame(ctrlFirst = g1_sub1,
                     ctrlSecond = g1_sub2,
                     objFirst = g2_sub1,
                     objSecond  = g2_sub2,
                     ctrlRatio = g1_sub1 / (g1_sub1 + g1_sub2),
                     objRatio = g2_sub1 / (g2_sub1 + g1_sub2),
                     ratioSD = dOR$ratioSD)
    
    colnames(df) <- c(
        paste('sumexp__',
              paste(rep(gn, each = length(sn)),
                    rep(sn, times = length(gn)),
                    sep = '__'),
              sep = ''),
        paste('ratio__', paste(gn, sn[1], sep = '__'), sep = ''),
        'ratio_sd')
    df
}



# Sample Values from Truncated Normal Distribution with Given Parameters
#' @importFrom stats rnorm
.hq.rnorm <- function(n, m, s, cutoff, meta_is_valid) {
    if (length(cutoff) == 1)
        cutoff <- rep(cutoff, n)
    y <- rep(NA_real_, n)
    invalid_idx <- rep(TRUE, n)
    n_tries <- 0
    while(sum(invalid_idx) > 0) {
        n_tries <- n_tries + 1
        y[invalid_idx] <- rnorm(sum(invalid_idx),
                                mean = m[invalid_idx], sd = s[invalid_idx])
        invalid_idx[invalid_idx] <- (y[invalid_idx] < 0) | (y[invalid_idx] > cutoff[invalid_idx])
        
        if (n_tries > 1000) {
            gene_ids <- seq_along(meta_is_valid)
            err_ids <- gene_ids[meta_is_valid][invalid_idx]
            stop('HomeoRoq attempted to sample ratios/counts ',
                 'from a truncated normal distribution up to 1000 times. ',
                 'However, it failed to sample valid samples. ',
                 'Try removing the homeologs on rows ',
                 paste0(err_ids, sep = ', ', collapse = ''),
                 ', then re-run `homeoroq()`')
        }
    }
    y
}


# Repeat Sampling to Compute P-values
.hq.sigChangeMH <- function(d, dOR, dOC, dOD, fv, iter_sampling) {
    n_genes <- seq_len(nrow(d$DATA))
    
    is_valid <- is.finite(dOC$r) & is.finite(dOC$rv) & is.finite(dOC$logmean) & is.finite(dOD$r)
    valid_genes <- n_genes[is_valid]
    n_valid_genes <- length(valid_genes)
    
    # pre-processing
    ctrlRatio <- dOC$r[valid_genes]
    ctrlLogMean <- dOC$logmean[valid_genes]
    ctrlExpMean <- exp(ctrlLogMean)
    ctrlVar   <- .hq.estimate_var(fv, ctrlLogMean, ctrlRatio, obs_v = dOC$v[valid_genes])
    ctrlResidue <- ctrlVar - exp(.hq.estimate_var(fv, ctrlLogMean, ctrlRatio))
    objRatio  <- dOD$r[valid_genes]
    objLogMean <- dOD$logmean[valid_genes]
    objExpMean <- exp(objLogMean)
    objVar    <- .hq.estimate_var(fv, objLogMean, objRatio, obs_v = dOD$v[valid_genes])
    objResidue <- objVar - exp(.hq.estimate_var(fv, objLogMean, objRatio))
    ratioSD <- sqrt(vapply(dOC$rv[valid_genes], function(i) max(i, 1.0e-2), numeric(1)))
    
    overallRatio <- dOR$r[valid_genes]
    overallRatioSD <- dOR$ratioSD[valid_genes]
    
    obs_prob <- .hq.calc_pval(dOR$r[valid_genes],
                        exp(dOC$logmean[valid_genes]) * dOC$r[valid_genes],
                        exp(dOD$logmean[valid_genes]) * dOD$r[valid_genes],
                        overallRatio, overallRatioSD, ctrlExpMean, objExpMean, 
                        sqrt(ctrlVar), sqrt(objVar))
    
    # sampling
    sig_counts <- rep(0, n_valid_genes)
    for (nidx in seq_len(iter_sampling)) {
        ix <- .hq.rnorm(n_valid_genes, overallRatio, overallRatioSD, 1, is_valid)
        ctrlSD <- sqrt(.hq.estimate_var(fv, ctrlLogMean, ix, residue = ctrlResidue))
        objSD <- sqrt(.hq.estimate_var(fv, objLogMean, ix, residue = objResidue))
        jx <- .hq.rnorm(n_valid_genes, ctrlExpMean * ix, ctrlSD, ctrlExpMean, is_valid)
        kx <- .hq.rnorm(n_valid_genes, objExpMean * ix, objSD, objExpMean, is_valid)
        est_probs <- .hq.calc_pval(ix, jx, kx,
                                   overallRatio, overallRatioSD, ctrlExpMean, objExpMean, 
                                   ctrlSD, objSD)
        sig_counts <- sig_counts + ifelse(est_probs <= obs_prob, 1, 0)
    }
    
    pvalues <- rep(NA_real_, length(n_genes))
    pvalues[is_valid] <- sig_counts / iter_sampling
    pvalues
}


# A Single Run of HomeoRoq test
#' @importFrom locfit locfit lp
.hq.homeoroq <- function(x, iter_sampling) {
    d <- .hq.data.format(x)
    
    # calculating parameters from observed rations
    dOR <- .hq.data.OR(d) # observed ratios
    dOC <- .hq.data.OC(d) # group 1
    dOD <- .hq.data.OD(d) # group 2
    
    # variance estimation model
    is_valid <- is.finite(dOR$v) & dOR$logmean > 0 & dOR$v > 1.0e-10
    fv <- locfit(log(dOR$v[is_valid]) ~ lp(dOR$logmean[is_valid],
                                           dOR$r[is_valid], scale = TRUE),
                 family = 'gaussian', maxk = 1000)
    
    # run sampling-resampling to calculate p-values
    .hq.sigChangeMH(d, dOR, dOC, dOD, fv, iter_sampling)
}


#' HomeoRoq: Detecting Shifts in Homeolog Expression Ratios between Two Groups
#'
#' A statistical method for detecting shifts in homeolog expression ratios
#' among two subgenomes in allopolyploids across two groups.
#' 
#' HomeoRoq estimates the probability that homeolog expression ratios
#' remain constant between two groups using RNA-Seq read count data for
#' each homeolog. Since this probability cannot be calculated analytically,
#' HomeoRoq uses a Bayesian framework and performs Markov chain Monte Carlo (MCMC)
#' sampling to estimate the distribution of homeolog expression ratios
#' under the null hypothesis of no change between conditions.
#' 
#' The `homeoroq()` function provides a comprehensive interface
#' for detecting shifts in homeolog expression ratios across conditions.
#' By default, it performs 10,000 sampling iterations
#' per MCMC chain (`iter_sampling = 1e4`) across 4 chains (`chians = 4`).
#' However, for robust and reproducible estimation of statistical significance
#' (i.e., _p_-values), it is recommended to conduct multiple independent sampling runs,
#' each with a sufficiently large number of iterations.
#' To replicate the settings used in the original HomeoRoq implementation
#' (Akama et al., NAR, 2014), set `chains = 10` and `iter_sampling = 1e4`.
#'
#' @param x An \linkS4class{ExpMX} class object containing normalized homeolog
#'      expression data (i.e., RNA-Seq read counts).
#' @param chains Integer. Specifies the number of independent sampling
#'      iterations used in the statistical test.
#' @param iter_sampling Integer. Specifies the number of samples to draw per
#'      iteration for estimating p-values.
#' @param n_threads Integer. Specifies the number of threads to use for parallel
#'      processing across sampling iterations.
#'      The default value is taken from the global option `'mc.cores'`,
#'      which can be set at the beginning of an R session as
#'      `options(mc.cores = 4)`. If `'mc.cores'` is not set, the default is 1.
#'
#' @return A data.frame with one row per homeolog,
#'      containing the following columns:
#'      \itemize{
#'          \item `pvalue`: p-value from the statistical test.
#'          \item `qvalue`: Adjusted p-value using the Benjamini-Hochberg method.
#'          \item `sumexp__*__$`: Total read counts for subgenome `$` under condition `*`.
#'          \item `ratio__*__$`: Homeolog expression ratio for subgenome `$` under condition `*`.
#'          \item `ratio_sd`: Standard deviation of homeolog expression ratios calculated from observed counts.
#'      }
#'      All statistics are computed directly from the observed read counts.
#'      Additionally, `NA` in the output indicates
#'      that the homeolog is expressed in only one condition,
#'      making ratio comparisons infeasible.
#'
#' @references Akama S, Shimizu-Inatsugi R, Shimizu KK, Sese J.
#'      Genome-wide quantification of homeolog expression ratio revealed
#'      nonstochastic gene regulation in synthetic allopolyploid Arabidopsis.
#'      Nucleic Acids Res. 2014;42(6):e46. \doi{10.1093/nar/gkt1376}
#'
#' @examples
#' x <- sim_homeolog_counts(100)
#' x_output <- homeoroq(x)
#' 
#' # sampling options
#' x_output <- homeoroq(x, chains = 4, iter_sampling = 1000, n_threads = 4)
#' 
#' @importFrom stats p.adjust
#' @importFrom locfit locfit lp
#' @importFrom future.apply future_vapply
#' @importFrom future plan multisession sequential
#' @importFrom progressr progressor with_progress handlers
#' @importFrom progress progress_bar
#' @export
homeoroq <- function(x,
                     chains = 4,
                     iter_sampling = 1e4,
                     n_threads = getOption('mc.cores', 1)) {
    if ((length(x@data) != 2) || (length(unique(x@exp_design$group)) != 2))
        stop('HomeoRoq only supports testing allopolyploids',
             'with two subgenomes under two conditions.')
    
    plan(multisession, workers = n_threads)
    handlers("progress")
    with_progress({
        pb <- progressor(chains)
        pvalues <- future_vapply(seq_len(chains), function(i) {
                outputs <- .hq.homeoroq(x, iter_sampling)
                pb()
                outputs
            },
            FUN.VALUE = numeric(nrow(x@data[[1]])),
            future.seed = TRUE,
            future.packages = c('locfit'))
    })
    plan(sequential)
    
    pvalues <- rowMeans(pvalues)
    data.frame(gene = x@gene_names,
               pvalue = pvalues,
               qvalue = p.adjust(pvalues, 'BH'),
               .hq.calc_obscounts_stats(x))
}
