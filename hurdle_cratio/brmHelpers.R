collapse <- function(..., sep = "") {
    paste(..., sep = sep, collapse = "")
}
extract <- function(x, ..., drop = FALSE, drop_dim = NULL) {
    if (!length(dim(x))) {
        return(x[...])
    }
    if (length(drop_dim)) {
        drop <- FALSE
    } else {
        drop <- as_one_logical(drop)
    }
    out <- x[..., drop = drop]
    if (drop || !length(drop_dim) || any(dim(out) == 0L)) {
        return(out)
    }
    if (is.numeric(drop_dim)) {
        drop_dim <- seq_along(dim(x)) %in% drop_dim
    }
    if (!is.logical(drop_dim)) {
        stop2("'drop_dim' needs to be logical or numeric.")
    }
    keep <- dim(out) > 1L | !drop_dim
    new_dim <- dim(out)[keep]
    if (length(new_dim) <= 1L) {
        # use vectors instead of 1D arrays
        new_dim <- NULL
    }
    dim(out) <- new_dim
    out
}
pordinal <- function(q, eta, thres, disc = 1, family = NULL, link = "logit") {
    family <- as_one_character(family)
    link <- as_one_character(link)
    args <- nlist(x = seq_len(max(q)), eta, thres, disc, link)
    p <- do_call(paste0("d", family), args)
    .fun <- function(j) rowSums(as.matrix(p[, 1:j, drop = FALSE]))
    cblapply(q, .fun)
}
cblapply <- function(X, FUN, ...) {
    do.call(cbind, lapply(X, FUN, ...))
}
first_greater <- function(A, target, i = 1) {
    ifelse(target <= A[, i] | ncol(A) == i, i, first_greater(A, target, i + 1))
}
slice_col <- function(x, i) {
    if (length(dim(x)) < 2L) {
        return(x)
    }
    slice(x, 2, i)
}
slice <- function(x, dim, i, drop = TRUE) {
    ndim <- length(dim(x))
    commas1 <- collapse(rep(", ", dim - 1))
    commas2 <- collapse(rep(", ", ndim - dim))
    drop_dim <- ifelse(drop, ", drop_dim = dim", "")
    expr <- paste0("extract(x, ", commas1, "i", commas2, drop_dim, ")")
    eval2(expr)
}
subset_thres <- function(prep, i) {
    thres <- prep$thres$thres
    Jthres <- prep$thres$Jthres
    if (!is.null(Jthres)) {
        thres <- thres[, Jthres[i, 1]:Jthres[i, 2], drop = FALSE]
    }
    thres
}
log_lik_weight <- function(x, i, prep) {
    weight <- prep$data$weights[i]
    if (!is.null(weight)) {
        x <- x * weight
    }
    x
}
as_one_character <- function(x, allow_na = FALSE) {
    s <- substitute(x)
    x <- as.character(x)
    if (length(x) != 1L || anyNA(x) && !allow_na) {
        s <- deparse0(s, max_char = 100L)
        stop2("Cannot coerce '", s, "' to a single character value.")
    }
    x
}
as_one_logical <- function(x, allow_na = FALSE) {
    s <- substitute(x)
    x <- as.logical(x)
    if (length(x) != 1L || anyNA(x) && !allow_na) {
        s <- deparse0(s, max_char = 100L)
        stop2("Cannot coerce '", s, "' to a single logical value.")
    }
    x
}
inv_link_cratio <- function(x, link) {
    x <- inv_link(x, link)
    ndim <- length(dim(x))
    dim_noncat <- dim(x)[-ndim]
    nthres <- dim(x)[ndim]
    marg_noncat <- seq_along(dim(x))[-ndim]
    ones_arr <- array(1, dim = c(dim_noncat, 1))
    dim_t <- c(nthres, dim_noncat)
    x_cumprod <- aperm(
        array(apply(x, marg_noncat, cumprod), dim = dim_t),
        perm = c(marg_noncat + 1, 1)
    )
    abind::abind(1 - x, ones_arr) * abind::abind(ones_arr, x_cumprod)
}
log_cdf <- function(x, link) {
    switch(link,
        logit = log_inv_logit(x),
        probit = pnorm(x, log.p = TRUE),
        cauchit = pcauchy(x, log.p = TRUE),
        cloglog = log1m_exp(-exp(x)),
        probit_approx = pnorm(x, log.p = TRUE),
        softit = log_inv_softit(x),
        stop2("Link '", link, "' is not supported.")
    )
}
log_inv_logit <- function(x) {
    log(inv_logit(x))
}
log1m_exp <- function(x) {
    ifelse(x < 0, log1p(-exp(x)), NaN)
}
dcratio <- function(x, eta, thres, disc = 1, link = "logit") {
    eta <- disc * (eta - thres)
    if (link == "identity") {
        out <- eta
    } else {
        out <- inv_link_cratio(eta, link = link)
    }
    out[, x, drop = FALSE]
}
inv_link <- function(x, link) {
    switch(link,
        identity = x,
        log = exp(x),
        logm1 = expp1(x),
        log1p = expm1(x),
        inverse = 1 / x,
        sqrt = x^2,
        "1/mu^2" = 1 / sqrt(x),
        tan_half = 2 * atan(x),
        logit = inv_logit(x),
        probit = pnorm(x),
        cauchit = pcauchy(x),
        cloglog = inv_cloglog(x),
        probit_approx = pnorm(x),
        softplus = log1p_exp(x),
        squareplus = (x + sqrt(x^2 + 4)) / 2,
        softit = inv_softit(x),
        stop2("Link '", link, "' is not supported.")
    )
}
logit <- function(p) {
    log(p) - log1p(-p)
}
inv_logit <- function(x) {
    1 / (1 + exp(-x))
}
inv_cloglog <- function(x) {
    1 - exp(-exp(x))
}
log1p_exp <- function(x) {
    out <- log1p(exp(x))
    ifelse(out < Inf, out, x)
}
inv_softit <- function(x) {
    y <- log1p_exp(x)
    y / (1 + y)
}
log_inv_softit <- function(x) {
    y <- log1p_exp(x)
    log(y) - log1p(y)
}
expp1 <- function(x) {
    exp(x) + 1
}
nlist <- function(...) {
    m <- match.call()
    dots <- list(...)
    no_names <- is.null(names(dots))
    has_name <- if (no_names) FALSE else nzchar(names(dots))
    if (all(has_name)) {
        return(dots)
    }
    nms <- as.character(m)[-1]
    if (no_names) {
        names(dots) <- nms
    } else {
        names(dots)[!has_name] <- nms[!has_name]
    }
    dots
}
do_call <- function(what, args, pkg = NULL, envir = parent.frame()) {
    call <- ""
    if (length(args)) {
        if (!is.list(args)) {
            stop2("'args' must be a list.")
        }
        fun_args <- names(args)
        if (is.null(fun_args)) {
            fun_args <- rep("", length(args))
        } else {
            nzc <- nzchar(fun_args)
            fun_args[nzc] <- paste0("`", fun_args[nzc], "` = ")
        }
        names(args) <- paste0(".x", seq_along(args))
        call <- paste0(fun_args, names(args), collapse = ",")
    } else {
        args <- list()
    }
    if (is.function(what)) {
        args$.fun <- what
        what <- ".fun"
    } else {
        what <- paste0("`", as_one_character(what), "`")
        if (!is.null(pkg)) {
            what <- paste0(as_one_character(pkg), "::", what)
        }
    }
    call <- paste0(what, "(", call, ")")
    eval2(call, envir = args, enclos = envir)
}
eval2 <- function(expr, envir = parent.frame(), ...) {
    if (is.character(expr)) {
        expr <- str2expression(expr)
    }
    eval(expr, envir, ...)
}
deparse0 <- function(x, max_char = NULL, ...) {
    out <- collapse(deparse(x, ...))
    if (isTRUE(max_char > 0)) {
        out <- substr(out, 1L, max_char)
    }
    out
}
stop2 <- function(message = "", ..., .subclass = NULL,
                  call = NULL, .envir = parent.frame()) {
    if (is.null(call)) {
        call <- rlang::caller_call()
    }

    rlang::abort(
        message   = glue::glue(message, ..., .envir = .envir),
        .subclass = c(.subclass, "brms_error"),
        call      = call
    )
}
