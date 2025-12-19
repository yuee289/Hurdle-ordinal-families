

hurdle_acat <- custom_family(
    "hurdle_acat",
    dpars       = c("mu", "hu", "disc"),
    links       = c("logit", "logit", "log"),
    specials    = c("ordinal", "extra_cat"),
    type        = "int",
    threshold   = "flexible"
)

# Define stan function
stan_funs <- "
    // acat
    real acat_logit_lpmf(int y, real mu, real disc, vector thres) {
        int nthres = num_elements(thres);
        real mu_logit = logit(mu);
        vector[nthres + 1] p = append_row(0, cumulative_sum(disc * (mu_logit - thres)));
        return p[y] - log_sum_exp(p);
    }

    // hurdle acat
    real hurdle_acat_lpmf(int y, real mu, real hu, real disc, vector thres) {
        if (y == 0) {
            return bernoulli_lpmf(1 | hu);
        } else {
            return bernoulli_lpmf(0 | hu) +
                acat_logit_lpmf(y | mu, disc, thres);
        }
    }

"
stan_vars <- stanvar(scode = stan_funs, block = "functions")

log_lik_hurdle_acat <- function(i, prep) {
    hu <- brms::get_dpar(prep, "hu", i = i)    
    mu <- brms::get_dpar(prep, "mu", i = i)
    mu_logit <- logit(mu)
    disc <- brms::get_dpar(prep, "disc", i = i)
    thres <- subset_thres(prep, i)
    nthres <- NCOL(thres)
    eta <- disc * (mu_logit - thres)
    y <- prep$data$Y[i]
    if (y == 0L) {
        out <- dbinom(1, size = 1, prob = hu, log = TRUE)
    } else {
        q <- sapply(1:nthres, function(k) eta[, k])
        p <- cbind(rep(0, nrow(eta)), q[, 1],
                matrix(0, nrow = nrow(eta), ncol = nthres - 1))
        if (nthres > 1L) {
            p[, 3:(nthres + 1)] <-
                sapply(3:(nthres + 1), function(k) rowSums(q[, 1:(k - 1)]))
        }
        out <- p[, y] - log(rowSums(exp(p))) + 
            dbinom(0, size = 1, prob = hu, log = TRUE)
    }
    log_lik_weight(out, i = i, prep = prep)
}
posterior_predict_hurdle_acat <- function(i, prep, ...) {
    mu <- get_dpar(prep, "mu", i = i)
    hu <- get_dpar(prep, "hu", i = i)
    disc <- get_dpar(prep, "disc", i = i)
    thres <- subset_thres(prep, i)
    nthres <- NCOL(thres)
    ndraws <- prep$ndraws
    p <- pordinal(
        seq_len(nthres + 1L),
        eta = logit(mu),
        disc = disc,
        thres = thres,
        family = "acat",
        link = prep$family$link
    )
    tmp <- runif(ndraws, 0, 1)
    draws <- ifelse(
        tmp < hu, 0L,
        first_greater(p, target = runif(prep$ndraws, min = 0, max = 1))
    )
    return(draws)
}
posterior_epred_hurdle_acat <- function(prep) {
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
        out[[i]] <- do_call(dacat, args_i)

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





