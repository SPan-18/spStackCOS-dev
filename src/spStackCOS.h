#include <R.h>
#include <Rinternals.h>

extern "C" {

  SEXP idist(SEXP coords1_r, SEXP n1_r, SEXP coords2_r, SEXP n2_r, SEXP p_r, SEXP D_r);

  SEXP R_cholRankOneUpdate(SEXP L_r, SEXP n_r, SEXP v_r, SEXP alpha_r, SEXP beta_r, SEXP lower_r);

  SEXP R_cholRowDelUpdate(SEXP L_r, SEXP n_r, SEXP row_r, SEXP lower_r);

  SEXP R_cholRowBlockDelUpdate(SEXP L_r, SEXP n_r, SEXP start_r, SEXP end_r, SEXP lower_r);

  SEXP spLMexact(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
                 SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                 SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                 SEXP nSamples_r, SEXP verbose_r);

  SEXP spLMexactLOO(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
                    SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r,
                    SEXP phi_r, SEXP nu_r, SEXP deltasq_r, SEXP corfn_r,
                    SEXP nSamples_r, SEXP loopd_r, SEXP loopd_method_r,
                    SEXP verbose_r);

  SEXP sptLMexactCOS(SEXP X_r, SEXP Xcov_r, SEXP N_r, SEXP r_r,
                     SEXP X_spcoords_r, SEXP X_timecoords_r,
                     SEXP pred1_r, SEXP Npred1_r, SEXP Xpredcov1_r, SEXP Xpred_spcoords1_r, SEXP Xpred_timecoords1_r,
                     SEXP pred2_r, SEXP Npred2_r, SEXP Xpredcov2_r, SEXP Xpred_timecoords2_r,
                     SEXP polyid_map_r, SEXP MCpolycoords_r, SEXP nMCpolycoords_r, SEXP poly_areas_r,
                     SEXP betaNorm_r, SEXP sigmaSqIG_r,
                     SEXP phi_s_r, SEXP phi_t_r, SEXP nu_r, SEXP deltasq_r,
                     SEXP nSamples_r);

  SEXP sptLMexactTimeAgg(SEXP X_r, SEXP Xcov_r, SEXP n_r, SEXP p_r,
                         SEXP X_spcoords_r, SEXP X_timecoords_r,
                         SEXP betaNorm_r, SEXP sigmaSqIG_r,
                         SEXP phi_s_r, SEXP phi_t_r, SEXP nu_r, SEXP deltasq_r,
                         SEXP nSamples_r);

  SEXP sptPointPredictTimeAgg(SEXP n_r, SEXP p_r, SEXP Npred_r, SEXP X_spcoords_r, SEXP X_timecoords_r,
                              SEXP Xpredcov_r, SEXP Xpred_spcoords_r, SEXP Xpred_timecoords_r,
                              SEXP post_sigmasq_r, SEXP post_z_r, SEXP post_beta_r,
                              SEXP phi_s_r, SEXP phi_t_r, SEXP nu_r, SEXP deltasq_r, SEXP nSamples_r);

  SEXP sptBlockPredictTimeAgg(SEXP n_r, SEXP p_r, SEXP Npred_r, SEXP spcoords_r, SEXP timecoords_r,
                              SEXP Xpredcov_r, SEXP MCpolycoords_r, SEXP nMCpolycoords_r, SEXP poly_areas_r,
                              SEXP polyid_map_r, SEXP timecoords_pred_r,
                              SEXP post_sigmasq_r, SEXP post_z_r, SEXP post_beta_r,
                              SEXP phi_s_r, SEXP phi_t_r, SEXP nu_r, SEXP deltasq_r, SEXP nSamples_r);
}
