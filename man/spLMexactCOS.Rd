% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spLMexactCOS.R
\name{spLMexactCOS}
\alias{spLMexactCOS}
\title{Univariate Bayesian spatial-temporal linear model for misaligned data}
\usage{
spLMexactCOS(
  X,
  Xcov,
  X_spcoords,
  X_timecoords,
  phi_s,
  phi_t,
  nu_s,
  noise_sp_ratio,
  Xpredcov1 = NULL,
  Xpred_spcoords1 = NULL,
  Xpred_timecoords1 = NULL,
  Xpred_polycov = NULL,
  Xpred_polycoords = NULL,
  Xpred_polytimecoords = NULL,
  nMCpoly_coords = 100,
  n.samples = 100,
  verbose = TRUE
)
}
\arguments{
\item{X}{A numeric vector (e.g., exposure).}

\item{Xcov}{A numeric basis matrix, used to model the mean function.}

\item{X_spcoords}{A numeric matrix of spatial coordinates (longitude, latitude).}

\item{X_timecoords}{A numeric matrix of time intervals (e.g., year, month).}

\item{phi_s}{A numeric value for the spatial decay parameter.}

\item{phi_t}{A numeric value for the temporal decay parameter.}

\item{nu_s}{A numeric value for the smoothness parameter.}

\item{noise_sp_ratio}{A numeric value for the noise to spatial-temporal variance ratio.}

\item{Xpredcov1}{A numeric matrix of covariates for the downscale prediction.}

\item{Xpred_spcoords1}{A numeric matrix of spatial coordinates for the downscale prediction.}

\item{Xpred_timecoords1}{A numeric vector of temporal coordinates for the downscale prediction.}

\item{Xpred_polycov}{A numeric matrix of covariates for upscale polygon-level prediction.}

\item{Xpred_polycoords}{A data frame with two columns: "ID" and "geometry", where "geometry" is an \code{sf} object.}

\item{Xpred_polytimecoords}{A numeric matrix of time intervals for the upscale polygon-level prediction.}

\item{nMCpoly_coords}{An integer specifying the number of Monte Carlo points to sample inside each polygon.}

\item{n.samples}{An integer specifying the number of posterior samples to draw.}

\item{verbose}{A logical value indicating whether to print progress messages.}
}
\value{
An object of class \code{spLMexactCOS} containing the fitted model and posterior samples.
}
\description{
Univariate Bayesian spatial-temporal linear model for misaligned data
}
