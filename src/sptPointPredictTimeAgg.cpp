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

    SEXP sptPointPredictTimeAgg(SEXP n_r, SEXP p_r, SEXP Npred_r, SEXP X_spcoords_r, SEXP X_timecoords_r,
                                SEXP Xpredcov_r, SEXP Xpred_spcoords_r, SEXP Xpred_timecoords_r,
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

    double *coords_sp = REAL(X_spcoords_r);
    double *coords_tm = REAL(X_timecoords_r);

    std::string corfn = "matern-exponential";
    double phi_s = REAL(phi_s_r)[0];
    double phi_t = REAL(phi_t_r)[0];
    double nu = REAL(nu_r)[0];
    double deltasq = REAL(deltasq_r)[0];
    int nSamples = INTEGER(nSamples_r)[0];

    // Handle inputs for prediction
    int Npred1 = INTEGER(Npred_r)[0];
    int NpredNpred1 = Npred1 * Npred1;

    double *Xpredcov1 = REAL(Xpredcov_r);
    double *Xpred_spcoords1 = REAL(Xpred_spcoords_r);
    double *Xpred_timecoords1 = REAL(Xpred_timecoords_r);

    // Read posterior samples
    double *post_sigmaSq = REAL(post_sigmasq_r);
    double *post_z = REAL(post_z_r);
    double *post_beta = REAL(post_beta_r);

    /*****************************************
     Set-up posterior predictive sample vector/matrices etc.
    *****************************************/
    double *Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);              // correlation matrix
    double *thetasp = (double *) R_alloc(3, sizeof(double));                         // spatial-temporal process parameters

    //construct covariance matrix (full)
    thetasp[0] = phi_s;
    thetasp[1] = nu;
    thetasp[2] = phi_t;
    sptCorCOSFull(n, 2, coords_sp, coords_tm, thetasp, corfn, Vz);

    double *Cz = (double *) R_chk_calloc(Npred1 * n, sizeof(double)); zeros(Cz, Npred1 * n);                    // allocate memory for n x Npred1 matrix
    double *Vzpred = (double *) R_chk_calloc(NpredNpred1, sizeof(double)); zeros(Vzpred, NpredNpred1);          // allocate memory for Npred1 x Npred1 matrix
    double *cholVz = (double *) R_chk_calloc(nn, sizeof(double)); zeros(cholVz, nn);                            // allocate memory for n x n matrix
    double *tmp_NpredNpred1 = (double *) R_chk_calloc(NpredNpred1, sizeof(double)); zeros(tmp_NpredNpred1, NpredNpred1);  // allocate memory for Npred1 x Npred1 matrix

    F77_NAME(dcopy)(&nn, Vz, &incOne, cholVz, &incOne);
    F77_NAME(dpotrf)(lower, &n, cholVz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}
    mkLT(cholVz, n);

    sptCorCrossCOS1(n, Npred1, 2, coords_sp, coords_tm, Xpred_spcoords1, Xpred_timecoords1, thetasp, corfn, Cz);
    sptCorFull(Npred1, 2, Xpred_spcoords1, Xpred_timecoords1, thetasp, corfn, Vzpred);

    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &Npred1, &one, cholVz, &n, Cz, &n FCONE FCONE FCONE FCONE);                 // Cz = cholinv(Vz)*Cz
    F77_NAME(dgemm)(ytran, ntran, &Npred1, &Npred1, &n, &one, Cz, &n, Cz, &n, &zero, tmp_NpredNpred1, &Npred1 FCONE FCONE);     // tmp_NpredNpred1 = t(Cz)*inv(Vz)*Cz
    F77_NAME(daxpy)(&NpredNpred1, &negOne, tmp_NpredNpred1, &incOne, Vzpred, &incOne);                                          // Vzpred = Vzpred - Cz*inv(Vz)*t(Cz)
    F77_NAME(dpotrf)(lower, &Npred1, Vzpred, &Npred1, &info FCONE); if(info != 0){perror("c++ error: Vzpred dpotrf failed\n");}
    mkLT(Vzpred, Npred1);

    // memory allocations
    double sigmasq1 = 0;
    double *z1 = (double *) R_chk_calloc(n, sizeof(double)); zeros(z1, n);                          // read n x 1 samples of z
    double *zpred1 = (double *) R_chk_calloc(Npred1, sizeof(double)); zeros(zpred1, Npred1);        // Npred1 x 1 posterior predictive z
    double *zpred1_mu = (double *) R_chk_calloc(Npred1, sizeof(double)); zeros(zpred1_mu, Npred1);  // Npred1 x 1 posterior predictive mean of z
    double *beta1 = (double *) R_chk_calloc(p, sizeof(double)); zeros(beta1, p);                    // read p x 1 samples of beta

    SEXP samples_zpred1_r = PROTECT(Rf_allocMatrix(REALSXP, Npred1, nSamples)); nProtect++;
    SEXP samples_ypred1_r = PROTECT(Rf_allocMatrix(REALSXP, Npred1, nSamples)); nProtect++;

    GetRNGstate();

    for(s = 0; s < nSamples; s++){

        sigmasq1 = post_sigmaSq[s];
        F77_NAME(dcopy)(&n, &post_z[s*n], &incOne, z1, &incOne);                                                          // z1 = samples_z_r[s]
        F77_NAME(dtrsv)(lower, ntran, nUnit, &n, cholVz, &n, z1, &incOne FCONE FCONE FCONE);                              // z1 = cholinv(Vz)*z1
        F77_NAME(dgemv)(ytran, &n, &Npred1, &one, Cz, &n, z1, &incOne, &zero, zpred1_mu, &incOne FCONE);                  // zpred1_mu = t(Cz)*inv(Vz)*z1
        for(i = 0; i < Npred1; i++){
            zpred1[i] = rnorm(0, sqrt(sigmasq1));                                                                         // zpred1 = rnorm(0, sigmasq1*I)
        }
        F77_NAME(dgemv)(ntran, &Npred1, &Npred1, &one, Vzpred, &Npred1, zpred1, &incOne, &one, zpred1_mu, &incOne FCONE); // zpred1_mu = zpred1_mu + chol(Vzpred)*zpred1
        F77_NAME(dcopy)(&Npred1, zpred1_mu, &incOne, &REAL(samples_zpred1_r)[s*Npred1], &incOne);                         // samples_zpred1_r[, s] = zpred1_mu

        F77_NAME(dcopy)(&p, &post_beta[s*p], &incOne, beta1, &incOne);                                                    // beta1 = samples_beta_r[s]
        F77_NAME(dgemv)(ntran, &Npred1, &p, &one, Xpredcov1, &Npred1, beta1, &incOne, &zero, zpred1, &incOne FCONE);      // zpred1 = Xpredcov1*beta1
        F77_NAME(daxpy)(&Npred1, &one, zpred1, &incOne, zpred1_mu, &incOne);                                              // zpred1_mu = zpred1_mu + Xpredcov1*beta1

        for(i = 0; i < Npred1; i++){
            REAL(samples_ypred1_r)[s*Npred1 + i] = zpred1_mu[i];                                                           // samples_zpred1_r[, s] = zpred1_mu
            // REAL(samples_ypred1_r)[s*Npred1 + i] = rnorm(zpred1_mu[i], delta * sqrt(sigmasq1));                           // samples_ypred1_r[j, s] = N(zpred1_mu, delta^2*sigmasq1)
        }

    }

    PutRNGstate();

    R_chk_free(Cz);
    R_chk_free(Vzpred);
    R_chk_free(cholVz);
    R_chk_free(tmp_NpredNpred1);
    R_chk_free(z1);
    R_chk_free(zpred1);
    R_chk_free(zpred1_mu);
    R_chk_free(beta1);


    int nResultListObjs = 2;

    SEXP result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    SEXP resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // samples of zpred1
    SET_VECTOR_ELT(result_r, 0, samples_zpred1_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("zpred"));

    // samples of ypred1
    SET_VECTOR_ELT(result_r, 1, samples_ypred1_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("xpred"));

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return result_r;

    }
}