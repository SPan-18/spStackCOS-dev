#' Fit a Bayesian linear model
#' @param y A numeric vector of response values.
#' @param X A numeric matrix of predictors (design matrix).
#' @param Z A vector consisting of a covariate.
#' @param wts A numeric vector of weights for the observations.
#' @param priors A list containing prior distributions for the model parameters.
#' @param n.samples Number of posterior samples to draw.
#' @return An object of class \code{fitLM} containing the fitted model and posterior samples.
#' @export
fitLM <- function(y, X, Z, wts, priors, n.samples = 1000){

    n <- length(y)
    p <- dim(X)[2]
    p1 <- p + 1
    # Check dimensions
    if (dim(X)[1] != n) {
        stop("The number of rows in X must be equal to the length of y")
    }
    if (length(Z) != n) {
        stop("The number of rows in Z must be equal to the length of y")
    }
    if (length(wts) != n) {
        stop("The length of wts must be equal to the length of y")
    }

    # Check priors
    if(missing(priors)){
        beta.Norm <- list(rep(0.0, p1), diag(1000.0, p1))
        sigma.sq.IG <- c(0.01, 0.01)
    }else{
        if (length(priors$beta.Norm) != 2) {
            stop("The length of beta.Norm must be equal to two")
        }
        if (length(priors$beta.Norm[[1]]) != p1) {
            stop("The first element of beta.Norm must have length equal to p + 1")
        }
        if (length(priors$beta.Norm[[2]]) != p1 * p1) {
            stop("The second element of beta.Norm must be a square matrix of size p + 1")
        }
        if (length(priors$sigma.sq.IG) != 2) {
            stop("The length of sigma.sq.IG must be equal to 2")
        }
        beta.Norm <- priors$beta.Norm
        sigma.sq.IG <- priors$sigma.sq.IG
    }

    beta.prior.mean <- beta.Norm[[1]]
    beta.prior.cov <- as.matrix(beta.Norm[[2]])
    prior.shape <- sigma.sq.IG[1]
    prior.rate <- sigma.sq.IG[2]

    V.beta.inv <- solve(beta.prior.cov)
    mu <- beta.prior.mean

    tX1X1 <- crossprod(X / sqrt(wts))
    tX1Z <- crossprod(X / wts, Z)
    tZZ <- crossprod(Z / sqrt(wts))
    tXX <- rbind(cbind(tX1X1, tX1Z), cbind(t(tX1Z), tZZ))

    tX1y <- crossprod(X / wts, y)
    tZy <- crossprod(Z / wts, y)
    tXy <- c(tX1y, tZy)

    yty <- sum(y * y / wts)

    muTVBetaInvmu <- t(mu) %*% V.beta.inv %*% mu

    posterior.precision <- V.beta.inv + tXX
    posterior.variance <- chol2inv(chol(posterior.precision))
    posterior.mean <- posterior.variance%*%(V.beta.inv%*%mu + tXy)

    posterior.shape <- prior.shape + n/2
    posterior.rate <- prior.rate + 0.5*(muTVBetaInvmu + yty - t(posterior.mean) %*% posterior.precision %*% posterior.mean)

    ptm <- proc.time()

    posterior.samples <- as.matrix(normalIGammaSampler(n.samples, posterior.mean, posterior.variance, posterior.shape, posterior.rate))

    run.time <- proc.time() - ptm

    # Extract the samples
    out <- list()
    out$y <- y
    out$X <- X
    out$Z <- Z
    out$wts <- wts
    out$priors <- list(beta.Norm = beta.Norm, sigma.sq.IG = sigma.sq.IG)
    out$samples <- list(beta = posterior.samples[, grep("beta", colnames(posterior.samples))],
                        sigmaSq = posterior.samples[, grep("sigmaSq", colnames(posterior.samples))])
    out$run.time <- run.time
    class(out) <- "fitLM"
    return(out)

}

#' Fit a Bayesian linear model under modular framework for COS
#' @param y A numeric vector of response values.
#' @param X A numeric matrix of predictors (design matrix).
#' @param Z A numeric matrix containg posterior samples of a covariate (e.g., spatial-temporal process).
#' Each column corresponds to a different sample, and number of rows must match the length of \code{y}.
#' @param wts A numeric vector of weights for the observations.
#' @param priors A list containing prior distributions for the model parameters.
#' @return An object of class \code{modularLM} containing the fitted model and posterior samples.
#' @export
modularLM <- function(y, X, Z, wts, priors){

    n <- length(y)
    p <- dim(X)[2]
    p1 <- p + 1
    # Check dimensions
    if (dim(X)[1] != n) {
        stop("The number of rows in X must be equal to the length of y")
    }
    if (dim(Z)[1] != n) {
        stop("The number of rows in Z must be equal to the length of y")
    }
    if (length(wts) != n) {
        stop("The length of wts must be equal to the length of y")
    }
    n.samples <- dim(Z)[2]

    # Check priors
    if(missing(priors)){
        beta.Norm <- list(rep(0.0, p1), diag(1000.0, p1))
        sigma.sq.IG <- c(0.01, 0.01)
    }else{
        if (length(priors$beta.Norm) != 2) {
            stop("The length of beta.Norm must be equal to two")
        }
        if (length(priors$beta.Norm[[1]]) != p1) {
            stop("The first element of beta.Norm must have length equal to p + 1")
        }
        if (length(priors$beta.Norm[[2]]) != p1 * p1) {
            stop("The second element of beta.Norm must be a square matrix of size p + 1")
        }
        if (length(priors$sigma.sq.IG) != 2) {
            stop("The length of sigma.sq.IG must be equal to 2")
        }
        beta.Norm <- priors$beta.Norm
        sigma.sq.IG <- priors$sigma.sq.IG
    }

    beta.prior.mean <- beta.Norm[[1]]
    beta.prior.cov <- as.matrix(beta.Norm[[2]])
    prior.shape <- sigma.sq.IG[1]
    prior.rate <- sigma.sq.IG[2]

    V.beta.inv <- solve(beta.prior.cov)
    mu <- beta.prior.mean

    # Compute the necessary quantities
    tX1X1 <- crossprod(X / sqrt(wts))
    tX1y <- crossprod(X / wts, y)
    yty <- sum(y * y / wts)
    muTVBetaInvmu <- t(mu) %*% V.beta.inv %*% mu
    posterior.shape <- prior.shape + n/2

    # storage for posterior samples
    posterior.samples <- array(0.0, dim = c(n.samples, p1 + 1))
    colnames(posterior.samples) <- c(paste0("beta", 1:p1), "sigmaSq")

    ptm <- proc.time()

    for(i in 1:n.samples){

        tX1Z <- crossprod(X / wts, Z[, i])
        tZZ <- crossprod(Z[, i] / sqrt(wts))
        tXX <- rbind(cbind(tX1X1, tX1Z), cbind(t(tX1Z), tZZ))
        tZy <- crossprod(Z[, i] / wts, y)
        tXy <- c(tX1y, tZy)

        # compute the posterior parameters
        posterior.precision <- V.beta.inv + tXX
        posterior.variance <- chol2inv(chol(posterior.precision))
        posterior.mean <- posterior.variance%*%(V.beta.inv%*%mu + tXy)
        posterior.rate <- prior.rate + 0.5*(muTVBetaInvmu + yty - t(posterior.mean) %*% posterior.precision %*% posterior.mean)

        posterior.samples[i, ] <- as.matrix(normalIGammaSampler(1, posterior.mean, posterior.variance, posterior.shape, posterior.rate))

    }

    run.time <- proc.time() - ptm

    # Extract the samples
    out <- list()
    out$y <- y
    out$X <- X
    out$Z <- Z
    out$wts <- wts
    out$priors <- list(beta.Norm = beta.Norm, sigma.sq.IG = sigma.sq.IG)
    out$samples <- list(beta = posterior.samples[, grep("beta", colnames(posterior.samples))],
                        sigmaSq = posterior.samples[, grep("sigmaSq", colnames(posterior.samples))])
    out$run.time <- run.time
    class(out) <- "modularLM"
    return(out)

}

#' Find ELPD for a fitted model using WAIC
#' @param fit An object of class \code{fitLM} or \code{modularLM} containing the fitted model.
#' @return A numeric value representing the WAIC (Widely Applicable Information Criterion) for the fitted model.
#' @importFrom stats dnorm
#' @importFrom stats var
#' @export
waicLM <- function(fit){

    if(inherits(fit, "fitLM")){

        # Extract the samples
        beta_samples <- fit$samples$beta
        sigma_sq_samples <- fit$samples$sigmaSq

        # Compute the log-likelihood for each sample
        log_lik <- matrix(0, nrow = nrow(beta_samples), ncol = length(fit$y))
        for (i in 1:nrow(beta_samples)) {
            mu <- cbind(fit$X, fit$Z) %*% beta_samples[i, ]
            log_lik[i, ] <- dnorm(fit$y, mean = mu, sd = sqrt(fit$wts * sigma_sq_samples[i]), log = TRUE)
        }

        # Compute WAIC
        lppd <- sum(apply(log_lik, 2, logSumExp) - log(nrow(log_lik)))
        p_waic2 <- sum(apply(log_lik, 2, var))
        waic <- -2 * (lppd - p_waic2)
        return(waic)

    }else if(inherits(fit, "modularLM")){

        # Extract the samples
        beta_samples <- fit$samples$beta
        sigma_sq_samples <- fit$samples$sigmaSq

        # Compute the log-likelihood for each sample
        log_lik <- matrix(0, nrow = nrow(beta_samples), ncol = length(fit$y))
        for (i in 1:nrow(beta_samples)) {
            mu <- cbind(fit$X, fit$Z[, i]) %*% beta_samples[i, ]
            log_lik[i, ] <- dnorm(fit$y, mean = mu, sd = sqrt(fit$wts * sigma_sq_samples[i]), log = TRUE)
        }

        # Compute WAIC
        lppd <- sum(apply(log_lik, 2, logSumExp) - log(nrow(log_lik)))
        p_waic2 <- sum(apply(log_lik, 2, var))
        waic <- -2 * (lppd - p_waic2)
        return(waic)

    }
}