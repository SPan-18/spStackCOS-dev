#' Point-level spatial-temporal prediction
#' @param mod.out An object of class \code{sptLMexactCOSTimeAgg} containing the fitted model.
#' @param Xcov_pred A numeric matrix of basis functions for the prediction.
#' @param spcoords_pred A numeric matrix of spatial coordinates for downscale prediction (longitude, latitude).
#' @param timecoords_pred A numeric matrix of temporal coordinates for downscale prediction.
#' @return An object of class \code{sptPointPredictTimeAgg} containing the posterior predictive samples.
#' @export
sptPointPredictTimeAgg <- function(mod.out, Xcov_pred, spcoords_pred, timecoords_pred){

    # Check if the model output is valid
    if (is.null(mod.out) || !inherits(mod.out, "sptLMexactCOSTimeAgg")) {
        stop("Invalid model output. Please provide a valid sptLMexactCOSTimeAgg object.")
    }

    # Read model output
    post_sigmasq <- mod.out$samples[["sigmaSq"]]
    post_z <- mod.out$samples[["z"]]
    post_beta <- mod.out$samples[["beta"]]

    phi_s <- mod.out$phi_s
    phi_t <- mod.out$phi_t
    nu_s <- mod.out$nu_s
    noise_sp_ratio <- mod.out$noise_sp_ratio
    n.samples <- mod.out$n.samples

    N <- length(mod.out$X)
    r <- dim(mod.out$Xcov)[2]
    spcoords <- mod.out$spcoords
    timecoords <- mod.out$timecoords

    N_pred <- nrow(Xcov_pred)
    if(dim(Xcov_pred)[2] != r){
        stop("Xcov_pred must have the same number of columns as Xcov")
    }
    if(dim(spcoords_pred)[1] != N_pred){
        stop("The number of rows in spcoords_pred must be equal to the number of rows in Xcov_pred")
    }
    if(dim(spcoords_pred)[2] != 2){
        stop("spcoords_pred must have two columns")
    }
    timecoords_pred <- matrix(timecoords_pred, nrow = N_pred, ncol = 1)
    if(dim(timecoords_pred)[1] != N_pred){
        stop("The number of rows in timecoords_pred must be equal to the number of rows in Xcov_pred")
    }
    if(dim(timecoords_pred)[2] != 1){
        stop("timecoords_pred must have one column")
    }

    storage.mode(N) <- "integer"
    storage.mode(r) <- "integer"
    storage.mode(spcoords) <- "double"
    storage.mode(timecoords) <- "double"
    storage.mode(phi_s) <- "double"
    storage.mode(phi_t) <- "double"
    storage.mode(noise_sp_ratio) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(nu_s) <- "double"

    storage.mode(post_sigmasq) <- "double"
    storage.mode(post_z) <- "double"
    storage.mode(post_beta) <- "double"

    storage.mode(N_pred) <- "integer"
    storage.mode(Xcov_pred) <- "double"
    storage.mode(spcoords_pred) <- "double"
    storage.mode(timecoords_pred) <- "double"

    ##### main function call #####
    ptm <- proc.time()

    samps <- .Call("sptPointPredictTimeAgg",
                   N, r, N_pred, spcoords, timecoords,
                   Xcov_pred, spcoords_pred, timecoords_pred,
                   post_sigmasq, post_z, post_beta,
                   phi_s, phi_t, nu_s, noise_sp_ratio, n.samples)

    run.time <- proc.time() - ptm

    # Extract the samples
    out <- list()
    out$Xcov_pred <- Xcov_pred
    out$spcoords_pred <- spcoords_pred
    out$timecoords_pred <- timecoords_pred
    out$samples <- samps[c("zpred", "xpred")]
    out$run.time <- run.time

    class(out) <- "sptPointPredictTimeAgg"
    return(out)

}

#' Block-level spatial-temporal prediction
#' @param mod.out An object of class \code{sptLMexactCOSTimeAgg} containing the fitted model.
#' @param Xcov_pred A numeric matrix of basis functions for the prediction.
#' @param polycoords_pred An \code{sf} object containing the polygons for upscale prediction.
#' @param timecoords_pred A numeric matrix of time intervals for upscale prediction.
#' @param nMCpoly_coords An integer specifying the number of Monte Carlo points to sample inside each polygon.
#' @param verbose A logical value indicating whether to print progress messages.
#' @return An object of class \code{sptBlockPredictTimeAgg} containing the posterior predictive samples.
#' @export
sptBlockPredictTimeAgg <- function(mod.out, Xcov_pred, polycoords_pred, timecoords_pred,
                                   nMCpoly_coords = 500, verbose = TRUE){

    # Check if the model output is valid
    if (is.null(mod.out) || !inherits(mod.out, "sptLMexactCOSTimeAgg")) {
        stop("Invalid model output. Please provide a valid sptLMexactCOSTimeAgg object.")
    }

    # Read model output
    post_sigmasq <- mod.out$samples[["sigmaSq"]]
    post_z <- mod.out$samples[["z"]]
    post_beta <- mod.out$samples[["beta"]]

    phi_s <- mod.out$phi_s
    phi_t <- mod.out$phi_t
    nu_s <- mod.out$nu_s
    noise_sp_ratio <- mod.out$noise_sp_ratio
    n.samples <- mod.out$n.samples

    N <- length(mod.out$X)
    r <- dim(mod.out$Xcov)[2]
    spcoords <- mod.out$spcoords
    timecoords <- mod.out$timecoords

    N_pred <- nrow(Xcov_pred)
    if(dim(Xcov_pred)[2] != r){
        stop("Xcov_pred must have the same number of columns as Xcov")
    }
    if(!inherits(polycoords_pred, "sf") || !is.data.frame(polycoords_pred)){
        stop("polycoords_pred must be both an sf object and a dataframe")
    }
    if(!all(c("ID", "geometry") %in% colnames(polycoords_pred))){
        stop("polycoords_pred must have columns 'ID' and 'geometry'")
    }
    if(dim(timecoords_pred)[1] != N_pred){
        stop("The number of rows in timecoords_pred must be equal to the number of rows in polycoords_pred")
    }
    if(dim(timecoords_pred)[2] != 2){
        stop("timecoords_pred must have two columns")
    }

    # read unique polygons from polycoords_pred
    unique_polysf <- polycoords_pred %>%
      distinct(across(all_of("ID")), .keep_all = TRUE)

    id_map <- match(polycoords_pred$ID, unique(unique_polysf$ID))
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
    poly.sample.time <- proc.time() - poly.sample.time.start

    if(verbose){
      cat("Elapsed ", format_elapsed_compact(poly.sample.time), ".\n")
    }

    storage.mode(N) <- "integer"
    storage.mode(r) <- "integer"
    storage.mode(spcoords) <- "double"
    storage.mode(timecoords) <- "double"
    storage.mode(phi_s) <- "double"
    storage.mode(phi_t) <- "double"
    storage.mode(noise_sp_ratio) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(nu_s) <- "double"

    storage.mode(post_sigmasq) <- "double"
    storage.mode(post_z) <- "double"
    storage.mode(post_beta) <- "double"

    storage.mode(N_pred) <- "integer"
    storage.mode(Xcov_pred) <- "double"
    storage.mode(timecoords_pred) <- "double"
    storage.mode(poly_areas) <- "double"
    storage.mode(id_map) <- "integer"
    storage.mode(nMCpoly_coords) <- "integer"

    # main function call
    ptm <- proc.time()

    samps <- .Call("sptBlockPredictTimeAgg",
                   N, r, N_pred, spcoords, timecoords,
                   Xcov_pred, montecarlo_coords, nMCpoly_coords, poly_areas, id_map,
                   timecoords_pred,
                   post_sigmasq, post_z, post_beta,
                   phi_s, phi_t, nu_s, noise_sp_ratio, n.samples)

    run.time <- proc.time() - ptm

    # Extract the samples
    out <- list()
    out$Xcov_pred <- Xcov_pred
    out$polycoords_pred <- polycoords_pred
    out$timecoords_pred <- timecoords_pred
    out$samples <- samps[c("zpred", "processPred")]
    out$run.time <- run.time

    class(out) <- "sptBlockPredictTimeAgg"

    return(out)

}