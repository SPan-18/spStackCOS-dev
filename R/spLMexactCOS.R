#' Univariate Bayesian spatial-temporal linear model for misaligned data
#' @param X A numeric vector (e.g., exposure).
#' @param Xcov A numeric basis matrix, used to model the mean function.
#' @param X_spcoords A numeric matrix of spatial coordinates (longitude, latitude).
#' @param X_timecoords A numeric matrix of time intervals (e.g., year, month).
#' @param phi_s A numeric value for the spatial decay parameter.
#' @param phi_t A numeric value for the temporal decay parameter.
#' @param nu_s A numeric value for the smoothness parameter.
#' @param noise_sp_ratio A numeric value for the noise to spatial-temporal variance ratio.
#' @param Xpredcov1 A numeric matrix of covariates for the downscale prediction.
#' @param Xpred_spcoords1 A numeric matrix of spatial coordinates for the downscale prediction.
#' @param Xpred_timecoords1 A numeric vector of temporal coordinates for the downscale prediction.
#' @param Xpred_polycov A numeric matrix of covariates for upscale polygon-level prediction.
#' @param Xpred_polycoords A data frame with two columns: "ID" and "geometry", where "geometry" is an \code{sf} object.
#' @param Xpred_polytimecoords A numeric matrix of time intervals for the upscale polygon-level prediction.
#' @param nMCpoly_coords An integer specifying the number of Monte Carlo points to sample inside each polygon.
#' @param n.samples An integer specifying the number of posterior samples to draw.
#' @param verbose A logical value indicating whether to print progress messages.
#' @return An object of class \code{spLMexactCOS} containing the fitted model and posterior samples.
#' @importFrom dplyr distinct
#' @importFrom magrittr %>%
#' @importFrom sf st_coordinates
#' @importFrom sf st_sample
#' @importFrom sf st_area
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @export
spLMexactCOS <- function(X, Xcov, X_spcoords, X_timecoords,
                         phi_s, phi_t, nu_s, noise_sp_ratio,
                         Xpredcov1 = NULL, Xpred_spcoords1 = NULL, Xpred_timecoords1 = NULL,
                         Xpred_polycov = NULL, Xpred_polycoords = NULL, Xpred_polytimecoords = NULL,
                         nMCpoly_coords = 100, n.samples = 100, verbose = TRUE){

  N <- length(X)
  if(dim(Xcov)[1] != N){
    stop("The number of rows in Xcov must be equal to the length of X")
  }

  r <- dim(Xcov)[2]

  if(dim(X_spcoords)[1] != N){
      stop("The number of rows in X_spcoords must be equal to the length of X")
  }

  if(dim(X_spcoords)[2] != 2){
      stop("X_spcoords must have two columns")
  }

  if(dim(X_timecoords)[1] != N){
      stop("The number of rows in X_timecoords must be equal to the length of X")
  }

  if(dim(X_timecoords)[2] != 2){
      stop("Xcoords must have two columns")
  }

  # prior setup
  beta.Norm <- list(rep(0.0, r), diag(1000.0, r))
  sigma.sq.IG <- c(2, 0.1)

  # posterior predictive setup
  if(is.null(Xpredcov1)){
    pred1 <- FALSE
    N_pred1 <- 0
    Xpredcov1 <- 0
    Xpred_spcoords1 <- 0
    Xpred_timecoords1 <- 0
  }else{
    pred1 <- TRUE
    N_pred1 <- nrow(Xpredcov1)
    if(dim(Xpredcov1)[2] != r){
        stop("Xpredcov1 must have the same number of columns as Xcov")
    }
    if(dim(Xpred_spcoords1)[1] != N_pred1){
        stop("The number of rows in Xpred_spcoords1 must be equal to the number of rows in Xpredcov1")
    }
    if(dim(Xpred_spcoords1)[2] != 2){
        stop("Xpred_spcoords1 must have two columns")
    }
    Xpred_timecoords1 <- matrix(Xpred_timecoords1, nrow = N_pred1, ncol = 1)
    if(dim(Xpred_timecoords1)[1] != N_pred1){
        stop("The number of rows in Xpred_timecoords1 must be equal to the number of rows in Xpredcov1")
    }
    if(dim(Xpred_timecoords1)[2] != 1){
        stop("Xpred_timecoords1 must have two columns")
    }
  }

  if(is.null(Xpred_polycov)){
    pred2 <- FALSE
    N_pred2 <- 0
    Xpred_polycov <- 0
    Xpred_polycoords <- 0
    Xpred_polytimecoords <- 0
    poly_areas <- 0
    id_map <- 0
    montecarlo_coords <- 0
    nMCpoly_coords <- 0
  }else{
    pred2 <- TRUE
    N_pred2 <- nrow(Xpred_polycov)
    if(dim(Xpred_polycov)[2] != r){
        stop("Xpred_polycov must have the same number of columns as Xcov")
    }
    if(!inherits(Xpred_polycoords, "sf") || !is.data.frame(Xpred_polycoords)){
        stop("Xpred_polycoords must be both an sf object and a dataframe")
    }
    if(!all(c("ID", "geometry") %in% colnames(Xpred_polycoords))){
        stop("Xpred_polycoords must have columns 'ID' and 'geometry'")
    }
    if(dim(Xpred_polytimecoords)[1] != N_pred2){
        stop("The number of rows in Xpred_polytimecoords must be equal to the number of rows in Xpred_polycoords")
    }
    if(dim(Xpred_polytimecoords)[2] != 2){
        stop("Xpred_polytimecoords must have two columns")
    }

    # read unique polygons from Xpred_polycoords
    unique_polysf <- Xpred_polycoords %>%
      distinct(across(all_of("ID")), .keep_all = TRUE)

    id_map <- match(Xpred_polycoords$ID, unique(unique_polysf$ID))
    unique_ids <- unique_polysf$ID
    unique_polygons <- unique_polysf$geometry

    if(verbose){
        cat("Sampling Monte Carlo points inside polygons...")
    }
    poly.sample.time.start <- proc.time()
    montecarlo_coords <- vector(mode = "list", length = length(unique_ids))
    poly_areas <- array(dim = length(unique_ids))
    for(i in 1:length(unique_ids)){
        MC_points <- st_sample(unique_polygons[[i]], size = 2 * nMCpoly_coords, type = "random")
        montecarlo_coords[[i]] <- st_coordinates(MC_points)
        poly_areas[i] <- st_area(unique_polygons[[i]])
    }
    poly.sample.time.end <- proc.time()
    poly.sample.time <- poly.sample.time.end - poly.sample.time.start
    if(verbose){
      cat("Elapsed ", format_elapsed_compact(poly.sample.time), ".\n")
    }

  }

  storage.mode(X) <- "double"
  storage.mode(Xcov) <- "double"
  storage.mode(X_spcoords) <- "double"
  storage.mode(X_timecoords) <- "double"
  storage.mode(phi_s) <- "double"
  storage.mode(phi_t) <- "double"
  storage.mode(noise_sp_ratio) <- "double"
  storage.mode(n.samples) <- "integer"
  storage.mode(N) <- "integer"
  storage.mode(r) <- "integer"
  storage.mode(nu_s) <- "double"

  storage.mode(Xpredcov1) <- "double"
  storage.mode(Xpred_spcoords1) <- "double"
  storage.mode(Xpred_timecoords1) <- "double"
  storage.mode(N_pred1) <- "integer"

  storage.mode(Xpred_polycov) <- "double"
  storage.mode(Xpred_polytimecoords) <- "double"
  storage.mode(N_pred2) <- "integer"
  storage.mode(id_map) <- "integer"
  storage.mode(nMCpoly_coords) <- "integer"
  storage.mode(poly_areas) <- "double"

  ##### main function call #####
  ptm <- proc.time()

  samps <- .Call("sptLMexactCOS",
                 X, Xcov, N, r, X_spcoords, X_timecoords,
                 pred1, N_pred1, Xpredcov1, Xpred_spcoords1, Xpred_timecoords1,
                 pred2, N_pred2, Xpred_polycov, Xpred_polytimecoords,
                 id_map, montecarlo_coords, nMCpoly_coords, poly_areas,
                 beta.Norm, sigma.sq.IG,
                 phi_s, phi_t, nu_s, noise_sp_ratio, n.samples)

  run.time <- proc.time() - ptm

  # Extract the samples
  out <- list()
  out$X <- X
  out$Xcov <- Xcov
  out$X_spcoords <- X_spcoords
  out$X_timecoords <- X_timecoords
  out$phi_s <- phi_s
  out$phi_t <- phi_t
  out$noise_sp_ratio <- noise_sp_ratio
  out$n.samples <- n.samples

  if(pred1){
    out$Xpredcov1 <- Xpredcov1
    out$Xpred_spcoords1 <- Xpred_spcoords1
    out$Xpred_timecoords1 <- Xpred_timecoords1
    out$samples <- samps[c("beta", "sigmaSq", "z", "zpred1", "ypred1")]
  }else if(pred2){
    out$Xpred_polycov <- Xpred_polycov
    out$Xpred_polytimecoords <- Xpred_polycoords
    out$Xpred_polytimecoords <- Xpred_polytimecoords
    out$samples <- samps[c("beta", "sigmaSq", "z", "zpred2", "mupred2")]
  }else{
    out$samples <- samps[c("beta", "sigmaSq", "z")]
  }
  out$run.time <- run.time

  class(out) <- "spLMexactCOS"

  return(out)



}