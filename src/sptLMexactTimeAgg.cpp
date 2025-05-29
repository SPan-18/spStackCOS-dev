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

    SEXP sptLMexactTimeAgg(SEXP X_r, SEXP Xcov_r, SEXP n_r, SEXP p_r,
                           SEXP X_spcoords_r, SEXP X_timecoords_r,
                           SEXP betaNorm_r, SEXP sigmaSqIG_r,
                           SEXP phi_s_r, SEXP phi_t_r, SEXP nu_r, SEXP deltasq_r,
                           SEXP nSamples_r){

    /*****************************************
     Common variables
     *****************************************/
    int i, j, s, info, nProtect = 0;
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
    double *Y = REAL(X_r);
    double *X = REAL(Xcov_r);
    int p = INTEGER(p_r)[0];
    int pp = p * p;
    int n = INTEGER(n_r)[0];
    int nn = n * n;
    int np = n * p;

    double *coords_sp = REAL(X_spcoords_r);
    double *coords_tm = REAL(X_timecoords_r);

    //priors
    double *betaMu = NULL;
    double *betaV = NULL;

    betaMu = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, REAL(VECTOR_ELT(betaNorm_r, 0)), &incOne, betaMu, &incOne);

    betaV = (double *) R_alloc(pp, sizeof(double));
    F77_NAME(dcopy)(&pp, REAL(VECTOR_ELT(betaNorm_r, 1)), &incOne, betaV, &incOne);

    double sigmaSqIGa = REAL(sigmaSqIG_r)[0];
    double sigmaSqIGb = REAL(sigmaSqIG_r)[1];

    std::string corfn = "matern-exponential";
    double phi_s = REAL(phi_s_r)[0];
    double phi_t = REAL(phi_t_r)[0];
    double nu = REAL(nu_r)[0];
    double deltasq = REAL(deltasq_r)[0];
    int nSamples = INTEGER(nSamples_r)[0];

    /*****************************************
     Set-up posterior sample vector/matrices etc.
     *****************************************/
    double sigmaSqIGaPost = 0, sigmaSqIGbPost = 0;
    double sse = 0;
    double dtemp = 0;
    double muBetatVbetaInvmuBeta = 0;

    // const double deltasqInv = 1.0 / deltasq;
    const double delta = sqrt(deltasq);

    double *Vz = (double *) R_alloc(nn, sizeof(double)); zeros(Vz, nn);              // correlation matrix
    double *cholVy = (double *) R_alloc(nn, sizeof(double)); zeros(cholVy, nn);      // allocate memory for n x n matrix
    double *thetasp = (double *) R_alloc(2, sizeof(double));                         // spatial process parameters

    double *tmp_n = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n, n);          // allocate memory for n x 1 vector

    double *tmp_p1 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p1, p);        // allocate memory for p x 1 vector
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double)); zeros(tmp_p2, p);        // allocate memory for p x 1 vector

    double *VbetaInv = (double *) R_alloc(pp, sizeof(double)); zeros(VbetaInv, pp);  // allocate VbetaInv
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp, pp);      // allocate memory for p x p matrix
    double *tmp_pp2 = (double *) R_alloc(pp, sizeof(double)); zeros(tmp_pp2, pp);    // allocate memory for p x p matrix

    /*****************************************
     Construct covariance matrix (full)
     *****************************************/
    //construct covariance matrix (full)
    thetasp[0] = phi_s;
    thetasp[1] = nu;
    thetasp[2] = phi_t;

    sptCorCOSFull(n, 2, coords_sp, coords_tm, thetasp, corfn, Vz);
    // F77_NAME(dpotrf)(lower, &n, Vz, &n, &info FCONE); if(info != 0){perror("c++ error: Vz dpotrf failed\n");}

    // construct marginal covariance matrix (Vz+deltasq*I)
    F77_NAME(dcopy)(&nn, Vz, &incOne, cholVy, &incOne);
    for(i = 0; i < n; i++){
      cholVy[i*n + i] += deltasq / fabs(coords_tm[i] - coords_tm[n + i]);
    }

    // find sse to sample sigmaSq
    // chol(Vy)
    F77_NAME(dpotrf)(lower, &n, cholVy, &n, &info FCONE); if(info != 0){perror("c++ error: Vy dpotrf failed\n");}

    // find YtVyInvY
    F77_NAME(dcopy)(&n, Y, &incOne, tmp_n, &incOne);                                         // tmp_n = Y
    F77_NAME(dtrsv)(lower, ntran, nUnit, &n, cholVy, &n, tmp_n, &incOne FCONE FCONE FCONE);  // tmp_n = cholinv(Vy)*Y
    dtemp = pow(F77_NAME(dnrm2)(&n, tmp_n, &incOne), 2);                                     // dtemp = t(Y)*VyInv*Y
    sse += dtemp;                                                                            // sse = YtVyinvY

    // find VbetaInvmuBeta
    F77_NAME(dcopy)(&pp, betaV, &incOne, VbetaInv, &incOne);                                                     // VbetaInv = Vbeta
    F77_NAME(dpotrf)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: VBetaInv dpotrf failed\n");} // VbetaInv = chol(Vbeta)
    F77_NAME(dpotri)(lower, &p, VbetaInv, &p, &info FCONE); if(info != 0){perror("c++ error: VBetaInv dpotri failed\n");} // VbetaInv = chol2inv(Vbeta)
    F77_NAME(dsymv)(lower, &p, &one, VbetaInv, &p, betaMu, &incOne, &zero, tmp_p2, &incOne FCONE);               // tmp_p2 = VbetaInv*muBeta

    // find muBetatVbetaInvmuBeta
    muBetatVbetaInvmuBeta = F77_CALL(ddot)(&p, betaMu, &incOne, tmp_p2, &incOne);                               // t(muBeta)*VbetaInv*muBeta
    sse += muBetatVbetaInvmuBeta;                                                                               // sse = YtVyinvY + muBetatVbetaInvmuBeta

    //  find XtVyInvY
    double *tmp_np = (double *) R_chk_calloc(np, sizeof(double)); zeros(tmp_np, np);                            // allocate temporary memory for n x p matrix
    F77_NAME(dcopy)(&np, X, &incOne, tmp_np, &incOne);                                                          // tmp_np = X
    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &p, &one, cholVy, &n, tmp_np, &n FCONE FCONE FCONE FCONE);  // tmp_np = cholinv(Vy)*X
    F77_NAME(dgemv)(ytran, &n, &p, &one, tmp_np, &n, tmp_n, &incOne, &zero, tmp_p1, &incOne FCONE);             // tmp_p1 = t(X)*VyInv*Y

    // find betahat = inv(XtVyInvX + VbetaInv)(XtVyInvY + VbetaInvmuBeta)
    F77_NAME(daxpy)(&p, &one, tmp_p2, &incOne, tmp_p1, &incOne);                                                // tmp_p1 = XtVyInvY + VbetaInvmuBeta
    F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, tmp_np, &n, tmp_np, &n, &zero, tmp_pp, &p FCONE FCONE);     // tmp_pp = t(X)*VyInv*X

    // deallocate tmp_np
    R_chk_free(tmp_np);

    F77_NAME(daxpy)(&pp, &one, VbetaInv, &incOne, tmp_pp, &incOne);                                             // tmp_pp = t(X)*VyInv*X + VbetaInv
    F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){perror("c++ error: tmp_pp dpotrf failed\n");}  // tmp_pp = chol(XtVyInvX + VbetaInv)
    F77_NAME(dtrsv)(lower, ntran, nUnit, &p, tmp_pp, &p, tmp_p1, &incOne FCONE FCONE FCONE);                    // tmp_p1 = cholinv(XtVyInvX + VbetaInv)*tmp_p1
    dtemp = pow(F77_NAME(dnrm2)(&p, tmp_p1, &incOne), 2);                                                       // dtemp = t(m)*M*m
    sse -= dtemp;                                                                                               // sse = YtVyinvY + muBetatVbetaInvmuBeta - mtMm

    // set-up for sampling spatial random effects
    double *tmp_nn2 = (double *) R_chk_calloc(nn, sizeof(double)); zeros(tmp_nn2, nn);                           // calloc n x n matrix
    F77_NAME(dcopy)(&nn, Vz, &incOne, tmp_nn2, &incOne);                                                         // tmp_nn2 = Vz
    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &n, &one, cholVy, &n, tmp_nn2, &n FCONE FCONE FCONE FCONE);  // tmp_nn2 = cholinv(Vy)*Vz
    F77_NAME(dtrsm)(lside, lower, ytran, nUnit, &n, &n, &one, cholVy, &n, tmp_nn2, &n FCONE FCONE FCONE FCONE);  // tmp_nn2 = inv(Vy)*Vz
    // tmp_nn2 = D*inv(Vy)*Vz
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            tmp_nn2[i*n + j] = tmp_nn2[i*n + j] / fabs(coords_tm[j] - coords_tm[n + j]);
        }
    }
    F77_NAME(dpotrf)(lower, &n, tmp_nn2, &n, &info FCONE); if(info != 0){perror("c++ error: inv(Vy)*Vz dpotrf failed\n");}  // tmp_nn2 = chol(inv(Vy)*Vz)
    mkLT(tmp_nn2, n);                                                                                            // make cholDinv lower-triangular

    // posterior parameters of sigmaSq
    sigmaSqIGaPost += sigmaSqIGa;
    sigmaSqIGaPost += 0.5 * n;

    sigmaSqIGbPost += sigmaSqIGb;
    sigmaSqIGbPost += 0.5 * sse;

    // posterior samples of sigma-sq and beta
    SEXP samples_sigmaSq_r = PROTECT(Rf_allocVector(REALSXP, nSamples)); nProtect++;
    SEXP samples_beta_r = PROTECT(Rf_allocMatrix(REALSXP, p, nSamples)); nProtect++;
    SEXP samples_z_r = PROTECT(Rf_allocMatrix(REALSXP, n, nSamples)); nProtect++;

    // sample storage at s-th iteration temporary allocation
    double sigmaSq = 0;
    double *beta = (double *) R_chk_calloc(p, sizeof(double)); zeros(beta, p);
    double *z = (double *) R_chk_calloc(n, sizeof(double)); zeros(z, n);

    GetRNGstate();

    for(s = 0; s < nSamples; s++){
      // sample sigmaSq from its marginal posterior
      dtemp = 1.0 / sigmaSqIGbPost;
      dtemp = rgamma(sigmaSqIGaPost, dtemp);
      sigmaSq = 1.0 / dtemp;
      REAL(samples_sigmaSq_r)[s] = sigmaSq;

      // sample fixed effects by composition sampling
      dtemp = sqrt(sigmaSq);
      for(j = 0; j < p; j++){
        beta[j] = rnorm(tmp_p1[j], dtemp);                                                   // beta ~ N(tmp_p1, sigmaSq*I)
      }
      F77_NAME(dtrsv)(lower, ytran, nUnit, &p, tmp_pp, &p, beta, &incOne FCONE FCONE FCONE); // beta = t(cholinv(tmp_pp))*beta

      dtemp = dtemp * delta;                                                                 // dtemp = sqrt(sigmaSq*deltasq)
      // sample spatial effects by composition sampling
      for(i = 0; i < n; i++){
        tmp_n[i] = rnorm(0.0, dtemp);                                                               // tmp_n ~ N(0, sigmaSq*I)
      }
      F77_NAME(dcopy)(&n, Y, &incOne, z, &incOne);                                                  // z = Y
      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &one, z, &incOne FCONE);        // z = Y-X*beta
      F77_NAME(dgemv)(ytran, &n, &n, &one, tmp_nn2, &n, z, &incOne, &one, tmp_n, &incOne FCONE);    // tmp_n = tmp_n + t(chol(tmp_nn2))*(Y-X*beta)/deltasq
      F77_NAME(dgemv)(ntran, &n, &n, &one, tmp_nn2, &n, tmp_n, &incOne, &zero, z, &incOne FCONE);   // z = chol(tmp_nn2)*tmp_n

      // copy samples into SEXP return object
      F77_NAME(dcopy)(&p, &beta[0], &incOne, &REAL(samples_beta_r)[s*p], &incOne);
      F77_NAME(dcopy)(&n, &z[0], &incOne, &REAL(samples_z_r)[s*n], &incOne);

    }

    PutRNGstate();

    // Free stuff
    R_chk_free(tmp_nn2);
    R_chk_free(beta);
    R_chk_free(z);

    // make return object for posterior samples
    SEXP result_r, resultName_r;
    int nResultListObjs = 3;

    result_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    resultName_r = PROTECT(Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    // samples of beta
    SET_VECTOR_ELT(result_r, 0, samples_beta_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("beta"));

    // samples of sigma-sq
    SET_VECTOR_ELT(result_r, 1, samples_sigmaSq_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("sigmaSq"));

    // samples of z
    SET_VECTOR_ELT(result_r, 2, samples_z_r);
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("z"));

    Rf_namesgets(result_r, resultName_r);

    UNPROTECT(nProtect);

    return result_r;

    }

}