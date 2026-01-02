hurdle_cratio <- custom_family(
    "hurdle_cratio",
    dpars       = c("mu", "hu", "disc"),
    links       = c("logit", "logit", "log"),
    specials    = c("ordinal", "extra_cat"),
    type        = "int",
    threshold   = "flexible"
)

stan_funs <- "
    // cratio
        real cratio_logit_lpmf(int y, real mu, real disc, vector thres) {
            int nthres = num_elements(thres);
            real mu_logit = logit(mu);
            vector[nthres + 1] p;
            vector[nthres] q;
            int k = 1;
            while (k <= min(y, nthres)) {
                q[k] = log_inv_logit(disc * (mu_logit - thres[k]));
                p[k] = log1m_exp(q[k]);
                for (kk in 1:(k - 1)) p[k] = p[k] + q[kk];
                k += 1;
            }
            if (y == nthres + 1) {
                p[nthres + 1] = sum(q);
            }
            return p[y];
        }

    // hurdle cratio
    real hurdle_cratio_lpmf(int y, real mu, real hu, real disc, vector thres) {
        if (y == 0){
            return bernoulli_lpmf(1 | hu);
        } else {
            return bernoulli_lpmf(0 | hu) +
                cratio_logit_lpmf(y | mu, disc, thres);
        }
    }

"
stan_vars <- stanvar(scode = stan_funs, block = "functions")

log_lik_hurdle_cratio <- function(i, prep) {
    hu <- get_dpar(prep, "hu", i = i)
    mu <- get_dpar(prep, "mu", i = i)
    mu_logit <- logit(mu)
    disc <- get_dpar(prep, "disc", i = i)
    thres <- subset_thres(prep, i)
    nthres <- NCOL(thres)
    eta <- disc * (mu_logit - thres)
    y <- prep$data$Y[i]
    q <- sapply(
        seq_len(min(y, nthres)),
        function(k) log_cdf(eta[, k], prep$family$link)
    )
    if (y == 0) {
        out <- out <- dbinom(1, size = 1, prob = hu, log = TRUE)
    } else if (y == 1L) {
        out <- log1m_exp(q[, 1L]) +
            dbinom(0, size = 1, prob = hu, log = TRUE)
    } else if (y == 2L) {
        out <- log1m_exp(q[, 2L]) + q[, 1L] +
            dbinom(0, size = 1, prob = hu, log = TRUE)
    } else if (y == nthres + 1L) {
        out <- rowSums(q) +
            dbinom(0, size = 1, prob = hu, log = TRUE)
    } else {
        out <- log1m_exp(q[, y]) + rowSums(q[, 1L:(y - 1L)]) +
            dbinom(0, size = 1, prob = hu, log = TRUE)
    }
    log_lik_weight(out, i = i, prep = prep)
}
posterior_predict_hurdle_cratio <- function(i, prep, ...) {
    mu <- get_dpar(prep, "mu", i = i)
    hu <- get_dpar(prep, "hu", i = i)
    disc <- get_dpar(prep, "disc", i = i)
    thres <- subset_thres(prep, i)
    nthres <- NCOL(thres)
    ndraws <- prep$ndraws
    p <- pordinal(
        seq_len(nthres + 1),
        eta = logit(mu),
        disc = disc,
        thres = thres,
        family = "cratio",
        link = prep$family$link
    )
    tmp <- runif(ndraws, 0, 1)
    draws <- ifelse(tmp < hu, 0,
        first_greater(p, target = runif(prep$ndraws, min = 0, max = 1))
    )
    return(draws)
}
posterior_epred_hurdle_cratio <- function(prep) {
    adjust <- ifelse(prep$family$link == "identity", 0, 1)
    ncat_max <- max(prep$data$nthres) + adjust
    nact_min <- min(prep$data$nthres) + adjust
    init_mat <- matrix(
        ifelse(prep$family$link == "identity", NA, 0),
        nrow = prep$ndraws, ncol = ncat_max - nact_min
    )
    args <- list(link = prep$family$link)
    out <- vector("list", prep$nobs)

    for (i in seq_along(out)) {
        args_i <- args
        args_i$eta <- logit(slice_col(get_dpar(prep, "mu", i)))
        args_i$disc <- slice_col(get_dpar(prep, "disc", i))
        args_i$thres <- subset_thres(prep, i)
        ncat_i <- NCOL(args_i$thres) + adjust
        args_i$x <- seq_len(ncat_i)
        out[[i]] <- do_call(dcratio, args_i)

        if (ncat_i < ncat_max) {
            sel <- seq_len(ncat_max - ncat_i)
            out[[i]] <- cbind(out[[i]], init_mat[, sel])
        }

        hu <- get_dpar(prep, "hu", i)
        out[[i]] <- cbind(hu, out[[i]] * (1 - hu))
    }

    out <- abind::abind(out, along = 3)
    out <- aperm(out, perm = c(1, 3, 2))
    dimnames(out)[[3]] <- c(paste0(0), seq_len(ncat_max))
    return(out)
}
