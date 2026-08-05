data {
    int<lower=0, upper=1> USE_DIRICHLET;
    real<lower=1e-12> EPS;

    int<lower=2> N_SUBGENOMES;
    int<lower=2> N_CONDITIONS;

    array[N_CONDITIONS]
        int<lower=1> N_REPLICATES;

    // H0:
    //   N_THETA = 1
    //   THETA_ID[c] = 1
    //
    // H1:
    //   N_THETA = N_CONDITIONS
    //   THETA_ID[c] = c
    int<lower=1, upper=N_CONDITIONS> N_THETA;

    array[N_CONDITIONS] int<lower=1, upper=N_THETA> THETA_ID;
    array[sum(N_REPLICATES), N_SUBGENOMES] int<lower=0> HOMEOLOG_EXP;

    // edgeR dispersion:
    //
    //   Var(Y) = mu + dispersion * mu^2
    matrix<lower=0>[N_CONDITIONS, N_SUBGENOMES] HOMEOLOG_EXP_DISPERSION;

    // Centered log effective library sizes.
    vector[sum(N_REPLICATES)] LOG_OFFSET;

    // One row for H0 and N_CONDITIONS rows for H1.
    matrix<lower=0>[N_THETA, N_SUBGENOMES] PRIOR_ALPHA;
}

parameters {
    vector[N_CONDITIONS] log_gene_exp_mu;

    // H0 has one shared simplex.
    // H1 has one simplex per condition.
    array[N_THETA] simplex[N_SUBGENOMES] theta;
}

model {
    if (USE_DIRICHLET == 1) {
        for (k in 1:N_THETA) {
            theta[k] ~ dirichlet(PRIOR_ALPHA[k]');
        }
    }

    {
        int r = 1;

        for (c in 1:N_CONDITIONS) {
            int theta_id = THETA_ID[c];

            for (i in 1:N_REPLICATES[c]) {
                for (s in 1:N_SUBGENOMES) {
                    real log_mu_raw = log_gene_exp_mu[c] + LOG_OFFSET[r] + log(theta[theta_id][s]);
                    real log_mu = log_sum_exp(log_mu_raw, log(EPS));
                    real precision = 1.0 /  HOMEOLOG_EXP_DISPERSION[c, s];
                    target += neg_binomial_2_log_lpmf(HOMEOLOG_EXP[r, s]| log_mu, precision);
                }
                r += 1;
            }
        }
    }
}

generated quantities {
    vector[N_CONDITIONS] gene_exp_mu = exp(log_gene_exp_mu);
    real log_lik_total = 0;

    {
        int r = 1;

        for (c in 1:N_CONDITIONS) {
            int theta_id = THETA_ID[c];

            for (i in 1:N_REPLICATES[c]) {
                for (s in 1:N_SUBGENOMES) {
                    real log_mu = log_sum_exp(log_gene_exp_mu[c] + LOG_OFFSET[r] + log(theta[theta_id][s]), log(EPS));
                    real precision = 1.0 / HOMEOLOG_EXP_DISPERSION[c, s];
                    log_lik_total += neg_binomial_2_log_lpmf(HOMEOLOG_EXP[r, s] | log_mu, precision);
                }
                r += 1;
            }
        }
    }
}
