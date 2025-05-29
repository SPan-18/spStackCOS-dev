#' Univariate Bayesian spatial-temporal linear model for misaligned data
#' @param X A numeric vector of values (e.g., exposure).
#' @param Xcov A numeric matrix of covariates (e.g., basis functions).
#' @param spcoords A numeric matrix of spatial coordinates (longitude, latitude).
#' @param timecoords A numeric matrix of time intervals (e.g., year, month).
#' @param phi_s A numeric value for the spatial decay parameter.
#' @param phi_t A numeric value for the temporal decay parameter.
#' @param nu_s A numeric value for the smoothness parameter.
#' @param noise_sp_ratio A numeric value for the noise to spatial-temporal variance ratio.
#' @param n.samples An integer specifying the number of posterior samples to draw.
#' @param verbose A logical value indicating whether to print progress messages.
#' @return An object of class \code{sptLMexactCOSTimeAgg} containing the fitted model and posterior samples.
#' @export
sptLMexactTimeAgg <- function(X, Xcov, spcoords, timecoords,
                              phi_s, phi_t, nu_s, noise_sp_ratio,
                              n.samples = 100, verbose = TRUE){

  N <- length(X)
  if(dim(Xcov)[1] != N){
    stop("The number of rows in Xcov must be equal to the length of X")
  }

  r <- dim(Xcov)[2]

  if(dim(spcoords)[1] != N){
      stop("The number of rows in spcoords must be equal to the length of X")
  }

  if(dim(spcoords)[2] != 2){
      stop("spcoords must have two columns")
  }

  if(dim(timecoords)[1] != N){
      stop("The number of rows in timecoords must be equal to the length of X")
  }

  if(dim(timecoords)[2] != 2){
      stop("Xcoords must have two columns")
  }

  # prior setup
  beta.Norm <- list(rep(0.0, r), diag(1000.0, r))
  sigma.sq.IG <- c(2, 0.1)

  storage.mode(X) <- "double"
  storage.mode(Xcov) <- "double"
  storage.mode(spcoords) <- "double"
  storage.mode(timecoords) <- "double"
  storage.mode(phi_s) <- "double"
  storage.mode(phi_t) <- "double"
  storage.mode(nu_s) <- "double"
  storage.mode(noise_sp_ratio) <- "double"
  storage.mode(n.samples) <- "integer"
  storage.mode(N) <- "integer"
  storage.mode(r) <- "integer"

  ##### main function call #####
  ptm <- proc.time()

  samps <- .Call("sptLMexactTimeAgg",
                 X, Xcov, N, r, spcoords, timecoords,
                 beta.Norm, sigma.sq.IG,
                 phi_s, phi_t, nu_s, noise_sp_ratio, n.samples)

  run.time <- proc.time() - ptm

  # Extract the samples
  out <- list()
  out$X <- X
  out$Xcov <- Xcov
  out$spcoords <- spcoords
  out$timecoords <- timecoords
  out$phi_s <- phi_s
  out$phi_t <- phi_t
  out$nu_s <- nu_s
  out$noise_sp_ratio <- noise_sp_ratio
  out$n.samples <- n.samples
  out$samples <- samps[c("beta", "sigmaSq", "z")]
  out$run.time <- run.time

  class(out) <- "sptLMexactCOSTimeAgg"

  return(out)

}