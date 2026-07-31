#ifndef HESPRESSO_HOBIT_HPP
#define HESPRESSO_HOBIT_HPP


// Maximum-likelihood negative-binomial model for HOBIT.
//
// H0:
//   N_THETA = 1
//   THETA_ID[c] = 0 for all conditions.
//
// H1:
//   N_THETA = N_CONDITIONS
//   THETA_ID[c] = c.
//
// HOMEOLOG_EXP_DISPERSION follows edgeR's dispersion convention:
//
//   Var(Y) = mu + dispersion * mu^2.
//
// For TMB::dnbinom_robust:
//
//   log(var - mu)
//     = log(dispersion * mu^2)
//     = 2 * log(mu) + log(dispersion).

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type hobit_mle(objective_function<Type>* obj) {
    DATA_INTEGER(N_SUBGENOMES);
    DATA_INTEGER(N_CONDITIONS);
    DATA_IVECTOR(N_REPLICATES);

    DATA_INTEGER(N_THETA);
    DATA_IVECTOR(THETA_ID);

    DATA_IMATRIX(HOMEOLOG_EXP);
    DATA_MATRIX(HOMEOLOG_EXP_DISPERSION);
    DATA_VECTOR(LOG_OFFSET);
    DATA_SCALAR(EPS);

    PARAMETER_VECTOR(log_gene_exp_mu);
    PARAMETER_MATRIX(theta_logit);

    matrix<Type> log_theta(
        N_THETA,
        N_SUBGENOMES
    );

    matrix<Type> theta(
        N_THETA,
        N_SUBGENOMES
    );

    for (int k = 0; k < N_THETA; ++k) {
        // log(exp(0)) for the reference subgenome.
        Type log_denom = Type(0.0);
        for (int s = 0; s < N_SUBGENOMES - 1; ++s) {
            log_denom = logspace_add(log_denom, theta_logit(k, s));
        }
        for (int s = 0; s < N_SUBGENOMES - 1; ++s) {
            log_theta(k, s) = theta_logit(k, s) - log_denom;
            theta(k, s) = exp(log_theta(k, s));
        }
        log_theta(k, N_SUBGENOMES - 1) = -log_denom;
        theta(k, N_SUBGENOMES - 1) = exp(-log_denom);
    }

    Type nll = Type(0.0);
    const Type log_eps = log(EPS);

    int r = 0;

    for (int c = 0; c < N_CONDITIONS; ++c) {
        const int theta_id = THETA_ID[c];

        for (int i = 0; i < N_REPLICATES[c]; ++i) {
            for (int s = 0; s < N_SUBGENOMES; ++s) {
                const Type log_mu_raw = log_gene_exp_mu[c] + LOG_OFFSET[r] + log_theta(theta_id, s);
                // const Type log_mu = CppAD::CondExpGt(log_mu_raw, log_eps, log_mu_raw, log_eps);
                const Type log_mu = logspace_add(log_mu_raw, log_eps);

                const Type log_dispersion = log(HOMEOLOG_EXP_DISPERSION(c, s));
                const Type log_var_minus_mu = Type(2.0) * log_mu + log_dispersion;

                nll -= dnbinom_robust(Type(HOMEOLOG_EXP(r, s)), log_mu, log_var_minus_mu, true);
            }
            ++r;
        }
    }

    vector<Type> gene_exp_mu(N_CONDITIONS);

    for (int c = 0; c < N_CONDITIONS; ++c) {
        gene_exp_mu[c] = exp(log_gene_exp_mu[c]);
    }

    const Type log_lik = -nll;

    REPORT(gene_exp_mu);
    REPORT(theta);
    REPORT(log_lik);

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
