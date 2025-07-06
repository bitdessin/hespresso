data {
    int USE_PRIOR;
    real EPS;
    real GENE_EXP_UPPER;

    // experiment design
    int<lower=1> N_SUBGENOMES;
    int<lower=1> N_CONDITIONS;
    array[N_CONDITIONS] int<lower=1> N_REPLICATES;
    
    // observed expression data
    array[sum(N_REPLICATES), N_SUBGENOMES] int<lower=0> HOMEOLOG_EXP;
    array[N_CONDITIONS, N_SUBGENOMES] real<lower=0> HOMEOLOG_EXP_PHI;
    
    // prior params
    vector[N_SUBGENOMES] PRIOR_ALPHA0;
    array[N_CONDITIONS] vector[N_SUBGENOMES] PRIOR_ALPHA1;
}

parameters {
    array[N_CONDITIONS] real<lower=0, upper=GENE_EXP_UPPER> gene_exp_mu;
    array[N_CONDITIONS] simplex[N_SUBGENOMES] theta1;
}

transformed parameters {
    array[N_CONDITIONS] vector<lower=0>[N_SUBGENOMES] homeolog_exp_mu;
    for (c in 1:N_CONDITIONS) {
        homeolog_exp_mu[c] = fmax(gene_exp_mu[c] * theta1[c], EPS);
    }
}

model {
    if (USE_PRIOR == 1) {
        for (c in 1:N_CONDITIONS) {
            theta1[c] ~ dirichlet(PRIOR_ALPHA1[c]);
        }
    }
    
    int r = 1;
    for (c in 1:N_CONDITIONS) {
        for (i in 1:N_REPLICATES[c]) {
            for (s in 1:N_SUBGENOMES) {
                target += neg_binomial_2_lpmf(HOMEOLOG_EXP[r, s] | homeolog_exp_mu[c, s], HOMEOLOG_EXP_PHI[c, s]);
            }
            r += 1;
        }
    }
}


generated quantities {
    vector[2] log_lik = rep_vector(0, 2);
    
    //array[sum(N_REPLICATES), N_SUBGENOMES] int hexp0;
    //array[sum(N_REPLICATES), N_SUBGENOMES] int hexp1;
    
    vector[N_SUBGENOMES] theta0 = rep_vector(0, N_SUBGENOMES);
    for (c in 1:N_CONDITIONS) {
        theta0 += theta1[c];
    }
    theta0 /= N_CONDITIONS;
    
    
    int r = 1;
    for (c in 1:N_CONDITIONS) {
        vector[N_SUBGENOMES] homeolog_exp_mu_0 = fmax(gene_exp_mu[c] * theta0, EPS);
        for (i in 1:N_REPLICATES[c]) {
            for (s in 1:N_SUBGENOMES) {
                log_lik[1] += neg_binomial_2_lpmf(HOMEOLOG_EXP[r, s] | homeolog_exp_mu_0[s], HOMEOLOG_EXP_PHI[c, s]);
                log_lik[2] += neg_binomial_2_lpmf(HOMEOLOG_EXP[r, s] | homeolog_exp_mu[c, s], HOMEOLOG_EXP_PHI[c, s]);
                //hexp0[r, s] = neg_binomial_2_rng(fmax(homeolog_exp_mu_0[s], EPS), HOMEOLOG_EXP_PHI[c, s]);
                //hexp1[r, s] = neg_binomial_2_rng(fmax(homeolog_exp_mu[c, s], EPS), HOMEOLOG_EXP_PHI[c, s]);
            }
            r += 1;
        }
    }
}
