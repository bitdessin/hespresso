library(parallel)
library(tools)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(pROC)
library(cmdstanr)
library(hespresso)

cmdstan_version()
options(mc.cores = 4)
options(dplyr.summarise.inform = FALSE)

MODEL_NAMES <- c('HOBIT_MLE' = 'MLE', 'HOBIT_MCMCSTRICT' = 'MCMC (strict)', 'HOBIT_MCMCFAST' = 'MCMC (fast)', 'HOMEOROQ' = 'HOMEOROQ')


parse_simdata_fpath <- function(data_fpath) {
    fname <- basename(data_fpath)
    fname <- strsplit(fname, '\\.')[[1]][1]
    fname <- strsplit(fname, '_')[[1]]
    n_subgenome <- as.integer(fname[2])
    n_rep <- as.integer(fname[3])
    seed <- as.integer(fname[4])
    list(n_subgenome = n_subgenome, n_rep = n_rep, seed = seed)
}


create_simdata <- function(data_fpath, seed, n_subgenomes, n_reps, overwrite = TRUE) {
    if (!file.exists(data_fpath) || overwrite) {
        set.seed(seed)
        x <- sim_homeolog_counts(n_genes = 100,
                                 n_subgenomes = n_subgenomes,
                                 n_replicates = c(n_reps, n_reps),
                                 group_names = c('ctrl', 'test'),
                                 prop_shift = runif(1, 0.03, 0.10),
                                 prop_deg = runif(1, 0.05, 0.25))
        saveRDS(x, data_fpath)
    }
}


run_test <- function(data_fpath, overwrite = FALSE) {
    output_fpath <- file.path(dirname(data_fpath),
                              paste0(file_path_sans_ext(basename(data_fpath)), '.MODELS.RData'))

    # load dataset
    x <- readRDS(data_fpath)

    # modeling
    MODELS <- list()
    PROC_TIMES <- list()
    if (file.exists(output_fpath))
        load(output_fpath)

    # HOBIT (MCMC fast)
    if (is.null(MODELS$HOBIT_MCMCFAST) || overwrite) {
        tryCatch({
            t0 <- proc.time()
            MODELS$HOBIT_MCMCFAST <- hobit(x, method = "mcmc", mode = 'fast', show_messages = FALSE)
            PROC_TIMES$HOBIT_MCMCFAST <- proc.time() - t0
            save(MODELS, PROC_TIMES, file = output_fpath)
        }, error = function(e) {
            MODELS$HOBIT_MCMCFAST <- NULL
            PROC_TIMES$HOBIT_MCMCFAST <- NULL
        })
    }

    # HOBIT (MCMC strict)
    if (is.null(MODELS$HOBIT_MCMCSTRICT) || overwrite) {
        tryCatch({
            t0 <- proc.time()
            MODELS$HOBIT_MCMCSTRICT <- hobit(x, method = "mcmc", mode = 'strict', show_messages = FALSE)
            PROC_TIMES$HOBIT_MCMCSTRICT <- proc.time() - t0
            save(MODELS, PROC_TIMES, file = output_fpath)
        }, error = function(e) {
            MODELS$HOBIT_MCMCSTRICT <- NULL
            PROC_TIMES$HOBIT_MCMCSTRICT <- NULL
        })
    }

    # HOBIT (MLE)
    if (is.null(MODELS$HOBIT_MLE) || overwrite || TRUE) {
        tryCatch({
            t0 <- proc.time()
            MODELS$HOBIT_MLE <- hobit(x, method = "mle", show_messages = FALSE)
            PROC_TIMES$HOBIT_MLE <- proc.time() - t0
            save(MODELS, PROC_TIMES, file = output_fpath)
        }, error = function(e) {
            MODELS$HOBIT_MLE <- NULL
            PROC_TIMES$HOBIT_MLE <- NULL
        })
    }

    # HOMEOROQ
    if (is.null(MODELS$HOMEOROQ) || overwrite) {
        tryCatch({
            t0 <- proc.time()
            MODELS$HOMEOROQ <- homeoroq(x, chains = 10, iter_sampling = 1e4, n_threads = 8)
            PROC_TIMES$HOMEOROQ <- proc.time() - t0
            save(MODELS, PROC_TIMES, file = output_fpath)
        }, error = function(e) {
            MODELS$HOMEOROQ <- NULL
            PROC_TIMES$HOMEOROQ <- NULL
        })
    }

    rm(x, MODELS, PROC_TIMES)
}


run_tests <- function(dpath, seeds, n_subgenomes, n_reps, overwrite = FALSE) {
    if (!dir.exists(dpath)) {
        dir.create(dpath, recursive = TRUE)
    }

    for (n_subgenome in n_subgenomes) {
        for (seed in seeds) {
            for (n_rep in n_reps) {
                data_fpath <- file.path(dpath, paste0('simdata_', n_subgenome, '_', n_rep, '_', seed, '.rds'))
                create_simdata(data_fpath, seed, n_subgenome, n_rep, overwrite = overwrite)
                run_test(data_fpath, overwrite = overwrite)
            }
        }
    }
}


calc_metrics <- function(gt, model_output, q_cutoff = 0.05) {
    model_output$pvalue[is.na(model_output$pvalue)] <- 1
    model_output$qvalue[is.na(model_output$qvalue)] <- 1

    tp <- sum(model_output$qvalue < q_cutoff & gt)
    fp <- sum(model_output$qvalue < q_cutoff & !gt)
    fn <- sum(model_output$qvalue >= q_cutoff & gt)
    tn <- sum(model_output$qvalue >= q_cutoff & !gt)
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    accuracy <- (tp + tn) / (tp + fp + fn + tn)
    f1 <- 2 * precision * recall / (precision + recall)
    auc <- roc(gt, 1 - model_output$pvalue, quiet = TRUE)$auc[1]

    data.frame(tp = tp, fp = fp, fn = fn, tn = tn,
               precision = precision,
               recall = recall,
               accuracy = accuracy,
               f1 = f1, auc = auc)
}


summarise_test <- function(data_fpath) {
    tags <- parse_simdata_fpath(data_fpath)
    model_fpath <- file.path(dirname(data_fpath),
                              paste0(file_path_sans_ext(basename(data_fpath)), '.MODELS.RData'))
    if (!file.exists(model_fpath)) return(NULL)

    x <- readRDS(data_fpath)
    load(model_fpath)

    gt <- def_sigShift(x)

    stats <- NULL
    for (model in c('HOBIT_MCMCFAST', 'HOBIT_MCMCSTRICT', 'HOBIT_MLE', 'HOMEOROQ')) {
        if (!is.null(MODELS[[model]])) {
            stats <- rbind(stats, data.frame(
                                model = MODEL_NAMES[model],
                                exe_time = as.numeric(PROC_TIMES[[model]][3]),
                                calc_metrics(gt, MODELS[[model]])))
        }
    }

    stats$model <- factor(stats$model, levels = MODEL_NAMES)
    stats
}


summarise_tests <- function(dpath) {
    stats <- NULL
    for (data_fpath in list.files(dpath, pattern = 'simdata_.*\\.rds$', full.names = TRUE)) {
        tags <- parse_simdata_fpath(data_fpath)
        stats <- rbind(stats, data.frame(
                    n_subgenomes = tags$n_subgenome,
                    n_replicates = tags$n_rep,
                    seed = tags$seed,
                    summarise_test(data_fpath)))
    }

    stats
}





if (FALSE) {
    simdata_dpath <- file.path('data', 'simdata')
    run_tests(simdata_dpath,
              seeds = 1:20,
              n_subgenomes = c(2, 3),
              n_reps = c(3, 5),
              overwrite = TRUE)

    stats <- summarise_tests(simdata_dpath)
    write.table(stats, file = file.path('data', 'hobit_metrics.tsv.gz'),
                sep = '\t', row.names = FALSE, quote = FALSE)

    stats %>% group_by(n_subgenomes, n_replicates, model) %>%
        summarise(exe_time = mean(exe_time),
                  tp = mean(tp), fp = mean(fp), fn = mean(fn), tn = mean(tn),
                  precision = mean(precision), recall = mean(recall),
                  accuracy = mean(accuracy), f1 = mean(f1), auc = mean(auc)) %>%
        ungroup()

    stats %>%
        pivot_longer(cols = c(precision, recall, f1, auc),
                     names_to = 'metric', values_to = 'value') %>%
        mutate(metric = factor(metric, levels = c('auc', 'f1', 'precision', 'recall'))) %>%
        ggplot(aes(x = model, y = value, group = model)) +
        geom_violin() +
        geom_jitter(width = 0.2) +
        facet_grid(metric ~ n_subgenomes + n_replicates, scales = 'free_y',
                   labeller = labeller(.rows = label_value, .cols = label_both)) +
        theme(axis.text.x = element_text(angle = 30, hjust = 0.5))

    stats %>%
        ggplot(aes(x = model, y = exe_time, group = model)) +
        geom_jitter(width = 0.2) +
        facet_grid( ~ n_subgenomes + n_replicates, scales = 'free_y',
                   labeller = labeller(.rows = label_value, .cols = label_both)) +
        theme(axis.text.x = element_text(angle = 30, hjust = 0.5)) 
}
