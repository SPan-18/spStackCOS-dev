#define USE_FC_LEN_T
#include <algorithm>
#include <string>
#include "util.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Memory.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

extern "C" {

    SEXP sptBlockPredictTimeAgg(SEXP n_r, SEXP p_r, SEXP Npred_r, SEXP spcoords_r, SEXP timecoords_r,
                                SEXP Xpredcov_r, SEXP MCpolycoords_r, SEXP nMCpolycoords_r, SEXP poly_areas_r,
                                SEXP polyid_map_r, SEXP timecoords_pred_r,
                                SEXP post_sigmasq_r, SEXP post_z_r, SEXP post_beta_r,
                                SEXP phi_s_r, SEXP phi_t_r, SEXP nu_r, SEXP deltasq_r, SEXP nSamples_r){

    /*****************************************
     Common variables
    *****************************************/
    int i, s, info, nProtect = 0;
    char const *lower = "L";
    char const *nUnit = "N";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *lside = "L";
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    const int incOne = 1;

    /*****************************************
     Set-up
    *****************************************/
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    int nn = n * n;

    double *coords_sp = REAL(spcoords_r);
    double *coords_tm = REAL(timecoords_r);

    std::string corfn = "matern-exponential";
    double phi_s = REAL(phi_s_r)[0];
    double phi_t = REAL(phi_t_r)[0];
    double nu = REAL(nu_r)[0];
    // double deltasq = REAL(deltasq_r)[0];
    // const double delta = sqrt(deltasq);
    int nSamples = INTEGER(nSamples_r)[0];

    int Npred2 = INTEGER(Npred_r)[0];
    int NpredNpred2 = Npred2 * Npred2;
    double *Xpredcov2 = REAL(Xpredcov_r);
    double *Xpred_timecoords2 = REAL(timecoords_pred_r);
    MatrixList MCpolycoords(MCpolycoords_r);
    int *polyID_map = INTEGER(polyid_map_r);
    int nMCpolycoords = INTEGER(nMCpolycoords_r)[0];
    // double *poly_areas = REAL(poly_areas_r);
    double jitter = 1e-1;

    // Read posterior samples
    double *post_sigmaSq = REAL(post_sigmasq_r);
    double *post_z = REAL(post_z_r);
    double *post_beta = REAL(post_beta_r);

    /*****************************************
     Set-up posterior sample vector/matrices etc.
    *****************************************/
    double *Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);              // correlation matrix
    double *thetasp = (double *) R_alloc(3, sizeof(double));                         // spatial-temporal process parameters

    //construct covariance matrix (full)
    thetasp[0] = phi_s;
    thetasp[1] = nu;
    thetasp[2] = phi_t;
    sptCorCOSFull(n, 2, coords_sp, coords_tm, thetasp, corfn, Vz);

    double *Cz = (double *) R_chk_calloc(Npred2 * n, sizeof(double)); zeros(Cz, Npred2 * n);                    // allocate memory for n x Npred2 matrix
    double *Vzpred = (double *) R_chk_calloc(NpredNpred2, sizeof(double)); zeros(Vzpred, NpredNpred2);          // allocate memory for Npred2 x Npred2 matrix
    double *cholVz = (double *) R_chk_calloc(nn, sizeof(double)); zeros(cholVz, nn);                            // allocate memory for n x n matrix
    double *tmp_NpredNpred2 = (double *) R_chk_calloc(NpredNpred2, sizeof(double)); zeros(tmp_NpredNpred2, NpredNpred2);  // allocate memory for Npred2 x Npred2 matrix

    F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
    F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
    mkLT(cholVz, n);

    sptCorCrossCOS2(n, Npred2, 2, coords_sp, coords_tm, MCpolycoords, nMCpolycoords, polyID_map, Xpred_timecoords2, thetasp, corfn, Cz);
    sptCorCOSFull2(Npred2, 2, MCpolycoords, nMCpolycoords, polyID_map, Xpred_timecoords2, thetasp, corfn, Vzpred);

    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &Npred2, &one, cholVz, &n, Cz, &n FCONE FCONE FCONE FCONE);             // Cz = cholinv(Vz)*Cz
    F77_NAME(dgemm)(ytran, ntran, &Npred2, &Npred2, &n, &one, Cz, &n, Cz, &n, &zero, tmp_NpredNpred2, &Npred2 FCONE FCONE); // tmp_NpredNpred2 = t(Cz)*inv(Vz)*Cz
    F77_NAME(daxpy)(&NpredNpred2, &negOne, tmp_NpredNpred2, &incOne, Vzpred, &incOne);                                      // Vzpred = Vzpred - Cz*inv(Vz)*t(Cz)
    for(i = 0; i < Npred2; i++){
        Vzpred[i*Npred2 + i] += jitter; // add jitter to diagonal
    }
    F77_NAME(dpotrf)(lower, &Npred2, Vzpred, &Npred2, &info FCONE); if(info != 0){perror("c++ error: Vzpred dpotrf failed\n");}
    mkLT(Vzpred, Npred2);

    // memory allocations
    double sigmasq2 = 0;
    double *z2 = (double *) R_chk_calloc(n, sizeof(double)); zeros(z2, n);                          // read n x 1 samples of z
    double *zpred2 = (double *) R_chk_calloc(Npred2, sizeof(double)); zeros(zpred2, Npred2);        // Npred2 x 1 posterior predictive z
    double *zpred2_mu = (double *) R_chk_calloc(Npred2, sizeof(double)); zeros(zpred2_mu, Npred2);  // Npred2 x 1 posterior predictive mean of z
    double *beta2 = (double *) R_chk_calloc(p, sizeof(double)); zeros(beta2, p);                    // read p x 1 samples of beta

    SEXP samples_zpred2_r = PROTECT(Rf_allocMatrix(REALSXP, Npred2, nSamples)); nProtect++;
    SEXP samples_mupred2_r = PROTECT(Rf_allocMatrix(REALSXP, Npred2, nSamples)); nProtect++;

    GetRNGstate();

    for(s = 0; s < nSamples; s++){

        sigmasq2 = post_sigmaSq[s];
        F77_NAME(dcopy)(&n, &post_z[s*n], &incOne, z2, &incOne);                                                          // z2 = samples_z_r[s]
        F77_NAME(dtrsv)(lower, ntran, nUnit, &n, cholVz, &n, z2, &incOne FCONE FCONE FCONE);                              // z2 = cholinv(Vz)*z2
        F77_NAME(dgemv)(ytran, &n, &Npred2, &one, Cz, &n, z2, &incOne, &zero, zpred2_mu, &incOne FCONE);                  // zpred2_mu = t(Cz)*inv(Vz)*z2
        for(i = 0; i < Npred2; i++){
            zpred2[i] = rnorm(0, sqrt(sigmasq2));                                                                         // zpred1 = rnorm(0, sigmasq1*I)
        }
        F77_NAME(dgemv)(ntran, &Npred2, &Npred2, &one, Vzpred, &Npred2, zpred2, &incOne, &one, zpred2_mu, &incOne FCONE); // zpred2_mu = zpred2_mu + chol(Vzpred)*zpred1
        F77_NAME(dcopy)(&Npred2, zpred2_mu, &incOne, &REAL(samples_zpred2_r)[s*Npred2], &incOne);                         // samples_zpred2_r[, s] = zpred2_mu

        F77_NAME(dcopy)(&p, &post_beta[s*p], &incOne, beta2, &incOne);                                                    // beta2 = samples_beta_r[s]
        F77_NAME(dgemv)(ntran, &Npred2, &p, &one, Xpredcov2, &Npred2, beta2, &incOne, &zero, zpred2, &incOne FCONE);      // zpred2 = Xpredcov2*beta2
        F77_NAME(daxpy)(&Npred2, &one, zpred2, &incOne, zpred2_mu, &incOne);                                              // zpred2_mu = zpred2_mu + Xpredcov2*beta2
        F77_NAME(dcopy)(&Npred2, zpred2_mu, &incOne, &REAL(samples_mupred2_r)[s*Npred2], &incOne);

    }

    PutRNGstate();

    R_chk_free(Cz);
    R_chk_free(Vzpred);
    R_chk_free(cholVz);
    R_chk_free(tmp_NpredNpred2);
    R_chk_free(z2);
    R_chk_free(zpred2);
    R_chk_free(zpred2_mu);
    R_chk_free(beta2);

    int nResultListObjs = 2;

    SEXP result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    SEXP resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // samples of zpred1
    SET_VECTOR_ELT(result_r, 0, samples_zpred2_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("zpred"));

    // samples of ypred1
    SET_VECTOR_ELT(result_r, 1, samples_mupred2_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("processPred"));

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return result_r;

    }
}