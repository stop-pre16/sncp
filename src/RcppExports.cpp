// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bd_process_test
Rcpp::List bd_process_test(arma::mat obs_points, double mean_mu_alpha, double sd_log_alpha, double sd_prop_alpha, double beta, int n_it, double window_hw, int df_iw_prior, int df_iw_prop, arma::mat sigma_prior, arma::mat lung_data, double var_mu_alpha, double pen_dist, double pen_val, int n_cent_init, double prior_n_cent, int max_bd_events, double max_bd_vt);
RcppExport SEXP _sncp_bd_process_test(SEXP obs_pointsSEXP, SEXP mean_mu_alphaSEXP, SEXP sd_log_alphaSEXP, SEXP sd_prop_alphaSEXP, SEXP betaSEXP, SEXP n_itSEXP, SEXP window_hwSEXP, SEXP df_iw_priorSEXP, SEXP df_iw_propSEXP, SEXP sigma_priorSEXP, SEXP lung_dataSEXP, SEXP var_mu_alphaSEXP, SEXP pen_distSEXP, SEXP pen_valSEXP, SEXP n_cent_initSEXP, SEXP prior_n_centSEXP, SEXP max_bd_eventsSEXP, SEXP max_bd_vtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type obs_points(obs_pointsSEXP);
    Rcpp::traits::input_parameter< double >::type mean_mu_alpha(mean_mu_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sd_log_alpha(sd_log_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sd_prop_alpha(sd_prop_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_it(n_itSEXP);
    Rcpp::traits::input_parameter< double >::type window_hw(window_hwSEXP);
    Rcpp::traits::input_parameter< int >::type df_iw_prior(df_iw_priorSEXP);
    Rcpp::traits::input_parameter< int >::type df_iw_prop(df_iw_propSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_prior(sigma_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lung_data(lung_dataSEXP);
    Rcpp::traits::input_parameter< double >::type var_mu_alpha(var_mu_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type pen_dist(pen_distSEXP);
    Rcpp::traits::input_parameter< double >::type pen_val(pen_valSEXP);
    Rcpp::traits::input_parameter< int >::type n_cent_init(n_cent_initSEXP);
    Rcpp::traits::input_parameter< double >::type prior_n_cent(prior_n_centSEXP);
    Rcpp::traits::input_parameter< int >::type max_bd_events(max_bd_eventsSEXP);
    Rcpp::traits::input_parameter< double >::type max_bd_vt(max_bd_vtSEXP);
    rcpp_result_gen = Rcpp::wrap(bd_process_test(obs_points, mean_mu_alpha, sd_log_alpha, sd_prop_alpha, beta, n_it, window_hw, df_iw_prior, df_iw_prop, sigma_prior, lung_data, var_mu_alpha, pen_dist, pen_val, n_cent_init, prior_n_cent, max_bd_events, max_bd_vt));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sncp_bd_process_test", (DL_FUNC) &_sncp_bd_process_test, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_sncp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
