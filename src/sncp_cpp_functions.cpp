////    Functions for sncp package    ////

#include <iostream>
//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins("cpp11")]]

//using namespace std;

const double log2pi = std::log(2.0 * M_PI);

////////////////////////////////////////////////////////////////////////////////
/*
 *Helper functions for doing multivariate normal density evaluation and random draws
 */
////////////////////////////////////////////////////////////////////////////////

//  Function to evaluate density of Wishart distribution
double dwish_arma(const arma::mat &W, const double &v, const arma::mat &S){
  int i, k = S.n_rows;
  double ret = 0;
  double gammapart = 1, num, denom, detS, detW, tracehold, dk = k;
  arma::mat hold;

  // denominator
  for(i = 1; i <= k; i++){
    gammapart = gammapart * R::gammafn((v + 1.0 - i) / 2.0);
  }
  denom = gammapart *  pow(2.0 ,(v * dk / 2.0)) * pow(M_PI, (dk*(dk-1.0)/4.0));

  //numerator
  detS = det(S);
  detW = det(W);
  hold = S.i() * W;

  tracehold = arma::trace(hold);
  num = pow(detS, (-v/2.0)) * pow(detW, ((v - dk - 1.0)/2.0)) * exp((-1.0 * tracehold) / 2.0);
  ret = num / denom;
  return(ret);
}

//  Function to evaluate density of inverse Wishart distribution
double diwish_arma(const arma::mat &W, const double &v, const arma::mat &S){
  int i, k = S.n_rows;
  double ret = 0;
  double gammapart = 1, num, denom, detS, detW, tracehold, dk = k;
  arma::mat hold;

  // denominator
  for(i = 1; i <= k; i++){
    gammapart = gammapart * R::gammafn((v + 1.0 - i) / 2.0);
  }
  denom = gammapart *  pow(2.0 ,(v * dk / 2.0)) * pow(M_PI, (dk*(dk-1.0)/4.0));

  //numerator
  detS = det(S);
  detW = det(W);
  hold = S * W.i();
  tracehold = arma::trace(hold);
  num = pow(detS, (v/2.0)) * pow(detW, (-(v + dk + 1.0)/2.0)) * exp((-1.0 * tracehold) / 2.0);
  ret = num / denom;
  return(ret);
}

//  Function to generate random draws from a Wishart distribution
arma::mat rwish_arma(const double &v, const arma::mat &S){
  int p = S.n_rows, i, j;
  arma::mat Z(p, p), CC, ret(p, p);
  arma::vec chi_vec(p);
  Z.zeros();
  CC = arma::chol(S);
  for(j = 0; j < p;j++){
    chi_vec(j) = R::rchisq(v-j);
  }
  Z.diag() = sqrt(chi_vec);
  for(i=1; i<p; i++){
    Z.diag(i) = arma::randn(p-i);
  }
  ret = arma::trans(Z*CC) * (Z*CC);
  return(ret);
}

//  Function to generate random draws from an inverse Wishart distribution
arma::mat riwish_arma(const double &v, const arma::mat &S){
  return(arma::inv(rwish_arma(v, arma::inv(S))));
}

//  Function that computes multivariate normal pdf for a single point
arma::vec dmvnrm_vec_arma_1f(const arma::mat &x,
                             const arma::vec &mean,
                             const arma::mat &sigma
) {
  unsigned int i;
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::mat rooti;
  arma::vec z, out;
  out.set_size(n);
  double rootisum, constants;
  constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  rootisum = arma::sum(log(rooti.diag()));

  for (i=0; i < n; i++) {
    z = rooti * arma::trans(x.row(i) - mean.t()) ;
    out(i) = constants - 0.5 * arma::sum(z%z) + rootisum;
  }
  out = exp(out);
  return(out);
}

//   Function that determines if any pairs of points within a matrix are
//   within a given (Euclidian) distance of each other (pen_dist)
arma::umat P_mat_gen1(const arma::mat &A, const double &pen_dist) {

  arma::colvec An =  arma::sum(square(A),1);
  arma::colvec Bn =  arma::sum(square(A),1);

  arma::mat C = -2 * (A * A.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  C = sqrt(C);

  return 1 * (C < pen_dist);
}

//   Function that determines if any points within a matrix are
//   within a given (Euclidian) distance of a separate point (pen_dist)
arma::uvec P_mat_gen2(const arma::mat &A, const arma::rowvec &B, const double &pen_dist) {

  arma::colvec An =  arma::sum(square(A),1);
  arma::colvec Bn =  arma::sum(square(B),1);

  arma::mat C = -2 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  C = sqrt(C);

  return 1 * (C < pen_dist);
}

arma::vec fastDMVNorm_diag_norm_sum(const arma::mat &A, const arma::mat &B, const double &sigma2) {

  arma::colvec An =  arma::sum(square(A),1);
  arma::colvec Bn =  arma::sum(square(B),1);
  arma::vec ret;
  int n = A.n_cols;
  double norm_const = pow(2.0 * M_PI * sigma2, n / 2.0);
  arma::mat C = -2.0 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  ret = arma::sum(exp(-(C) / (2.0 * sigma2)), 1);
  return (ret / norm_const);
}

////////////////////////////////////////////////////////////////////////////////
/*
 * Functions for dealing with linked list elements that represent clusters in the SNCP
 */
////////////////////////////////////////////////////////////////////////////////

//   Structure declaration for nodes in the linked list
struct Node {
  arma::vec mvn_dens;
  double log_alpha;
  double alpha;
  double log_bd_surf;
  arma::mat sigma;
  arma::rowvec loc;
  int n_cent_close;
  Node* next;
};

// Constructor, only for the 1st Node
void initNode(struct Node *head,
              const arma::rowvec &loc,
              const double &mu_alpha,
              const double &sd_log_alpha,
              const arma::mat &sigma_prior,
              const int &df_iw_prior,
              const arma::mat &obs_points,
              const double &log_bd_surf){
  head -> next = NULL;
  head -> loc = loc;
  head -> log_alpha = R::rnorm(mu_alpha, sd_log_alpha);
  head -> alpha = exp(head -> log_alpha);
  head -> sigma = riwish_arma(df_iw_prior, (df_iw_prior-3) * sigma_prior);
  head -> mvn_dens = dmvnrm_vec_arma_1f(obs_points, head -> loc.t(), head -> sigma);
  head -> n_cent_close = 0;
  head -> log_bd_surf = log_bd_surf;
}

// Apending nodes to end of list
void addNode(struct Node *head,
             const arma::rowvec &loc,
             const double &mu_alpha,
             const double &sd_log_alpha,
             const arma::mat &sigma_prior,
             const int &df_iw_prior,
             const arma::mat &obs_points,
             const double &log_bd_surf){
  Node *newNode = new Node;
  newNode -> loc = loc;
  newNode -> log_alpha = R::rnorm(mu_alpha, sd_log_alpha);
  newNode -> alpha = exp(newNode -> log_alpha);
  newNode -> sigma = riwish_arma(df_iw_prior, (df_iw_prior-3) * sigma_prior);
  newNode -> mvn_dens = dmvnrm_vec_arma_1f(obs_points, newNode -> loc.t(), newNode -> sigma);
  newNode -> n_cent_close = 0;
  newNode -> log_bd_surf = log_bd_surf;
  newNode->next = NULL;

  Node *cur = head;
  while(cur) {
    if(cur->next == NULL) {
      cur->next = newNode;
      return;
    }
    cur = cur->next;
  }

}

//   Adding new node to the front of the LL
void insertFront(struct Node **head,
                 const arma::rowvec &loc,
                 const double &mu_alpha,
                 const double &sd_log_alpha,
                 const arma::mat &sigma_prior,
                 const int &df_iw_prior,
                 const arma::mat &obs_points,
                 const double &log_bd_surf){
  Node *newNode = new Node;
  newNode -> loc = loc;
  newNode -> log_alpha = R::rnorm(mu_alpha, sd_log_alpha);
  newNode -> alpha = exp(newNode -> log_alpha);
  newNode -> sigma = riwish_arma(df_iw_prior, (df_iw_prior-3) * sigma_prior);
  newNode -> mvn_dens = dmvnrm_vec_arma_1f(obs_points, newNode -> loc.t(), newNode -> sigma);
  newNode -> n_cent_close = 0;
  newNode -> log_bd_surf = log_bd_surf;
  newNode->next = *head;

  *head = newNode;
}

//   Display location attribute for each node (cluster center) in the list
void display(struct Node *head) {
  Node *list = head;
  while(list) {
    Rcpp::Rcout << list -> loc << std::endl;
    list = list->next;
  }
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << std::endl;
}

//   Remove a node from the list
bool deleteNode(struct Node **head, Node *ptrDel) {
  Node *cur = *head;
  if(ptrDel == *head) {
    *head = cur->next;
    delete ptrDel;
    return true;
  }

  while(cur) {
    if(cur->next == ptrDel) {
      cur->next = ptrDel->next;
      delete ptrDel;
      return true;
    }
    cur = cur->next;
  }
  return false;
}

//   Search for a node at a specific position
struct Node *searchNode(struct Node *head, const int &n) {
  Node *cur = head;
  int cur_pos = 0;
  while(cur) {
    if(cur_pos == n) return cur;
    cur_pos++;
    cur = cur->next;
  }
  Rcpp::Rcout << "No Node " << n << " in list.\n";
  return NULL;
}

//   Calculate the Poisson process likelihood for all observed points
arma::vec calc_PP_lik_vec(Node *head){
  arma::vec ret(head -> mvn_dens.n_elem);
  ret.zeros();
  Node *cur = head;
  while(cur){
    ret += (cur -> alpha) * (cur -> mvn_dens);
    cur = cur -> next;
  }
  return ret;
}

//   Return a vector with log_alpha values in the LL
arma::vec get_log_alphas(Node *head, const int &n_cent){
  arma::vec ret(n_cent);
  ret.zeros();
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    ret(idx_cur) = (cur -> log_alpha);
    cur = cur -> next;
    idx_cur++;
  }
  return ret;
}

//   Return matrix of cluster center locations in the LL
arma::mat get_centers(Node *head, const int &n_cent, const int &xdim){
  arma::mat ret(n_cent, xdim);
  ret.zeros();
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    ret.row(idx_cur) = (cur -> loc);
    cur = cur -> next;
    idx_cur++;
  }
  return ret;
}

//   Return matrix of likelihood values for each observed point and cluster center pair
arma::mat get_L_mat(Node *head, const int &n_cent, const int &n_points){
  arma::mat ret(n_points, n_cent);
  ret.zeros();
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    ret.col(idx_cur) = (cur -> mvn_dens);
    cur = cur -> next;
    idx_cur++;
  }
  return ret;
}

//   Return cube of dispersion matrices from LL
arma::cube get_sigmas(Node *head, const int &n_cent, const int &xdim){
  arma::cube ret( xdim,  xdim, n_cent);
  ret.zeros();
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    ret.slice(idx_cur) = (cur -> sigma);
    cur = cur -> next;
    idx_cur++;
  }
  return ret;
}

//   Function to update penalty parameter for Strauss process repulsion
void update_Ps(Node *head, const int &n_cent, const double &pen_dist, const int &xdim){
  arma::mat center_sub(n_cent, xdim);
  center_sub.zeros();
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    center_sub.row(idx_cur) = (cur -> loc);
    cur = cur -> next;
    idx_cur++;
  }
  cur = head;
  while(cur){
    cur -> n_cent_close = arma::sum(P_mat_gen2(center_sub, cur -> loc, pen_dist)) - 1;
    cur = cur -> next;
  }
}

//   Function to calculated BD-MCMC death rates for each node in LL
arma::vec calc_DR_vec(Node *head,
                      const int &n_cent,
                      const arma::vec &PP_lik_cur,
                      const double &ll_cur,
                      const double &beta_it,
                      const int &n_points,
                      const double &log_pen_val,
                      const double &prior_n_cent_log){
  arma::vec DR_vec(n_cent);
  DR_vec.zeros();
  Node *cur = head;
  int idx_cur = 0;
  double bd_ll_tmp;

  while(cur){
    bd_ll_tmp = arma::sum(arma::log(PP_lik_cur - ((cur -> alpha) * (cur -> mvn_dens)) + beta_it));
    DR_vec(idx_cur) = exp(bd_ll_tmp + (cur -> log_bd_surf) - ll_cur + (cur -> alpha) - prior_n_cent_log -
      ((cur -> n_cent_close) * log_pen_val));
    cur = cur -> next;
    idx_cur++;
  }

  return DR_vec;
}

//   Function to calculated BD-MCMC death rates for each node in LL
arma::vec calc_DR_vec2(Node *head,
                       const int &n_cent,
                       const arma::vec &PP_lik_cur,
                       const double &ll_cur,
                       const double &beta_it,
                       const int &n_points,
                       const double &log_pen_val,
                       const double &log_bd_surf,
                       const double &prior_n_cent_log){
  arma::vec DR_vec(n_cent);
  DR_vec.zeros();
  Node *cur = head;
  int idx_cur = 0;
  double bd_ll_tmp;

  while(cur){
    bd_ll_tmp = arma::sum(arma::log(PP_lik_cur - ((cur -> alpha) * (cur -> mvn_dens)) + beta_it));
    DR_vec(idx_cur) = exp(bd_ll_tmp + (log_bd_surf) - ll_cur + (cur -> alpha) - prior_n_cent_log -
      ((cur -> n_cent_close) * log_pen_val));
    cur = cur -> next;
    idx_cur++;
  }

  return DR_vec;
}

//   Update log_alpha values in list based on vector of same length
void update_log_alphas(Node *head, const arma::vec &log_alpha_vec){
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    cur -> log_alpha = log_alpha_vec(idx_cur);
    cur -> alpha = exp(cur -> log_alpha);
    cur = cur -> next;
    idx_cur++;
  }
}

//   Update center locations in LL based on matrix of values
void update_centers(Node *head, const arma::mat &centers){
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    cur -> loc = centers.row(idx_cur);
    cur = cur -> next;
    idx_cur++;
  }
}

//    Update dispersion matrices in LL based on cube of new values
void update_sigmas(Node *head, const arma::cube &sigmas){
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    cur -> sigma = sigmas.slice(idx_cur);
    cur = cur -> next;
    idx_cur++;
  }
}

//   Update mvn density values based on matrix of updated values
void update_mvn_dens(Node *head, const arma::mat &L_mat){
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    cur -> mvn_dens = L_mat.col(idx_cur);
    cur = cur -> next;
    idx_cur++;
  }
}

//   Update log_alpha values in list based on vector of same length
void update_log_bd_surf(Node *head, const arma::vec &log_bd_vec){
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    cur -> log_bd_surf = log_bd_vec(idx_cur);
    cur = cur -> next;
    idx_cur++;
  }
}

//   Update all node values
void update_node_values(Node *head,
                        const arma::mat &centers,
                        const arma::vec &log_alpha_vec,
                        const arma::cube &sigmas,
                        const arma::mat &L_mat){
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    cur -> mvn_dens = L_mat.col(idx_cur);
    cur -> sigma = sigmas.slice(idx_cur);
    cur -> loc = centers.row(idx_cur);
    cur -> log_alpha = log_alpha_vec(idx_cur);
    cur -> alpha = exp(cur -> log_alpha);
    cur = cur -> next;
    idx_cur++;
  }
}

//   Update all node values
void update_node_values2(Node *head,
                         const arma::mat &centers,
                         const arma::vec &log_alpha_vec,
                         const arma::cube &sigmas,
                         const arma::mat &L_mat,
                         const arma::vec &log_bd_vec){
  Node *cur = head;
  int idx_cur = 0;
  while(cur){
    cur -> mvn_dens = L_mat.col(idx_cur);
    cur -> sigma = sigmas.slice(idx_cur);
    cur -> loc = centers.row(idx_cur);
    cur -> log_alpha = log_alpha_vec(idx_cur);
    cur -> alpha = exp(cur -> log_alpha);
    cur -> log_bd_surf = log_bd_vec(idx_cur);
    cur = cur -> next;
    idx_cur++;
  }
}

/*
 * Functions to do version with uniform dispersion density
 */

arma::rowvec fast_eig2d(const arma::mat &A) {
  arma::rowvec ret(2);
  double T = arma::trace(A);
  double B = arma::det(A);
  double C = sqrt((pow(T, 2) / 4.0) - B);
  ret(0) = (T / 2.0) + C;
  ret(1) = (T / 2.0) - C;
  return (ret);
}

arma::vec mahaInt_bin(const arma::mat & X,
                      const arma::vec & mu,
                      const arma::mat & sigma,
                      const unsigned int ncores,
                      const bool isChol = false){
  using namespace arma;

  // Some sanity checks
  //if(ncores < 0) Rcpp::stop("ncores has to be positive.");
  if(mu.n_elem != sigma.n_cols) Rcpp::stop("The mean vector has a different dimensions from the covariance matrix.");
  if(X.n_cols != sigma.n_cols)  Rcpp::stop("The number of columns of X is different from the dimension of the covariance matrix.");

  // Calculate transposed cholesky dec. unless sigma is alread a cholesky dec.
  mat cholDec;
  if( isChol == false ) {
    cholDec = trimatl(chol(sigma).t());
  }
  else{
    cholDec = trimatl(sigma.t());
    if(any(cholDec.diag() <= 0.0))  Rcpp::stop("The supplied cholesky decomposition has values <= 0.0 on the main diagonal.");
  }

  vec D = cholDec.diag();

  vec out(X.n_rows);
  uvec out2(X.n_rows);

#ifdef _OPENMP
#pragma omp parallel num_threads(ncores) if(ncores > 1)
{
#endif

  // Declaring some private variables
  uint32_t d = X.n_cols;
  uint32_t n = X.n_rows;

  vec tmp(d);

  double acc;
  uint32_t icol, irow, ii;

  // For each of the "n" random vectors, forwardsolve the corresponding linear system.
  // Forwardsolve because I'm using the lower triangle Cholesky.
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for(icol = 0; icol < n; icol++)
  {

    for(irow = 0; irow < d; irow++)
    {
      acc = 0.0;

      for(ii = 0; ii < irow; ii++) acc += tmp.at(ii) * cholDec.at(irow, ii);

      tmp.at(irow) = ( X.at(icol, irow) - mu.at(irow) - acc ) / D.at(irow);
    }

    out.at(icol) = sum(square(tmp));
  }

#ifdef _OPENMP
}
#endif
double area = arma::prod(arma::sqrt(fast_eig2d(sigma))) * M_1_PI;
out2 = out <= 1;
return ((1.0 / area) * (arma::conv_to<arma::vec>::from(out2)));
}


// Constructor, only for the 1st Node
void initNode_u(struct Node *head,
                const arma::rowvec &loc,
                const double &mu_alpha,
                const double &sd_log_alpha,
                const arma::mat &sigma_prior,
                const int &df_iw_prior,
                const arma::mat &obs_points){
  //head->data = n;
  head -> next = NULL;
  head -> loc = loc;
  head -> log_alpha = R::rnorm(mu_alpha, sd_log_alpha);
  head -> alpha = exp(head -> log_alpha);
  head -> sigma = riwish_arma(df_iw_prior, (df_iw_prior-3) * sigma_prior);
  //head -> mvn_dens = dmvnrm_vec_arma_1f(obs_points, head -> loc.t(), head -> sigma);
  head -> mvn_dens = mahaInt_bin(obs_points, head -> loc.t(), head -> sigma, 1, false);
  head -> n_cent_close = 0;
}

// Apending nodes to end of list
void addNode_u(struct Node *head,
               const arma::rowvec &loc,
               const double &mu_alpha,
               const double &sd_log_alpha,
               const arma::mat &sigma_prior,
               const int &df_iw_prior,
               const arma::mat &obs_points){
  Node *newNode = new Node;
  //newNode->data = n;
  newNode -> loc = loc;
  newNode -> log_alpha = R::rnorm(mu_alpha, sd_log_alpha);
  newNode -> alpha = exp(newNode -> log_alpha);
  newNode -> sigma = riwish_arma(df_iw_prior, (df_iw_prior-3) * sigma_prior);
  //newNode -> mvn_dens = dmvnrm_vec_arma_1f(obs_points, newNode -> loc.t(), newNode -> sigma);
  newNode -> mvn_dens = mahaInt_bin(obs_points, newNode -> loc.t(), newNode -> sigma, 1, false);
  newNode -> n_cent_close = 0;
  newNode->next = NULL;

  Node *cur = head;
  while(cur) {
    if(cur->next == NULL) {
      cur->next = newNode;
      return;
    }
    cur = cur->next;
  }

}

//   Adding new node to the front of the LL
void insertFront_u(struct Node **head,
                   const arma::rowvec &loc,
                   const double &mu_alpha,
                   const double &sd_log_alpha,
                   const arma::mat &sigma_prior,
                   const int &df_iw_prior,
                   const arma::mat &obs_points){
  Node *newNode = new Node;
  newNode -> loc = loc;
  newNode -> log_alpha = R::rnorm(mu_alpha, sd_log_alpha);
  newNode -> alpha = exp(newNode -> log_alpha);
  newNode -> sigma = riwish_arma(df_iw_prior, (df_iw_prior-3) * sigma_prior);
  //newNode -> mvn_dens = dmvnrm_vec_arma_1f(obs_points, newNode -> loc.t(), newNode -> sigma);
  newNode -> mvn_dens = mahaInt_bin(obs_points, newNode -> loc.t(), newNode -> sigma, 1, false);
  newNode -> n_cent_close = 0;
  newNode->next = *head;

  *head = newNode;
}




////////////////////////////////////////////////////////////////////////////////
/*
 * Wrapper functions to do BD-MCMC for whole data set that are exported to R
 */
////////////////////////////////////////////////////////////////////////////////
using namespace Rcpp;

//' Bayesian SNCP fit using BD-MCMC
//'
//' Run a BD-MCMC chain for a SNCP (mvn dispersion density) on lung ct data using a uniform proposal surface for the BD process
//'
//'
//'
//'
//'
//' @param obs_points matrix of coordinates for observed points
//' @param mean_mu_alpha prior mean for mu_alpha
//' @param sd_log_alpha prior std. dev. for log_alphas
//' @param sd_prop_alpha std. dev. for random walk proposal for log_alphas
//' @param beta prior estimate of homogeneous Poisson intensity function for cluster centers
//' @param n_it number of iterations to run MCMC
//' @param window_hw half-width of random walk cube for proposing new cluster center locations
//' @param df_iw_prior degrees of freedom for inverse-wishart prior on dispersion matrices
//' @param df_iw_prop degrees of freedom for inverse-wishart proposal on dispersion matrices
//' @param sigma_prior prior for dispersion matrices sigma
//' @param lung_data matrix that contains coordinates of all points within 2d slice of lung CT that is lung tissue
//' @param var_mu_alpha prior variance on mu_alpha
//' @param pen_dist distance for strauss process repulsion
//' @param pen_val vaule between 0 and 1 that is the penalty for 2 cluster centers being within pen_dist of each other
//' @param n_cent_init number of clusters to start with
//' @param prior_n_cent prior value for number of clusters
//' @param max_bd_events max events (births + deaths) to allow at each iteration of BD-MCMC
//' @param max_bd_vt max ammount of virtual time to spend in BD process at each BD_MCMC iteration
//'
//' @author Brian Vestal
//'
//' @return
//' Returns a named list containing posterior samples of model parameters with the following elements:
//' \itemize{
//' \item log_alpha_sample = list whose elements are vectors containing the estimated log-alpha estimates for the clusters present in that iteration
//' \item centers_sample = list whose elements matrices containing the estimated center locations for the clusters present in that iteration
//' \item sigmas_sample = list whose elements are lists containing the estimated covariance matrix estimates for the clusters present in that iteration
//' \item mu_alpha_sample = vector of mu_alpha samples
//' \item eta_sample = vector of eta samples,
//' \item n_centers_sample = vector that contains the number of estimated clusters at each iteration,
//' \item cumulative_intensity_sample = vector that contains the cumulative intensity estimate for the SNCP at each iteration
//' }
//'
//' @export
// [[Rcpp::export]]

Rcpp::List sncp_bdmcmc(arma::mat obs_points,
                       double mean_mu_alpha,
                       double sd_log_alpha,
                       double sd_prop_alpha,
                       double beta,
                       int n_it,
                       double window_hw,
                       int df_iw_prior,
                       int df_iw_prop,
                       arma::mat sigma_prior,
                       arma::mat lung_data,
                       double var_mu_alpha,
                       double pen_dist,
                       double pen_val,
                       int n_cent_init,
                       double prior_n_cent,
                       int max_bd_events,
                       double max_bd_vt){
  int i, j, xdim = obs_points.n_cols, n_cent_it, idx_birth, n_points = obs_points.n_rows;
  double ll_cur, ll_prop, mh_log, log_alpha_prop, beta_it = beta, mu_alpha_it, log_pen_val = log(pen_val);
  arma::vec log_alpha_it, beta_sample(n_it), mu_alpha_sample(n_it), xwin(2), ywin(2);
  //arma::mat log_alpha_sample(n_it, n_cent);
  arma::umat P_mat;
  arma::uvec P_vec_tmp, bd_idx_inf;
  arma::rowvec center_prop;
  arma::mat centers_it(n_cent_init, xdim), L_mat, L_mat_tmp, centers_tmp, sigma_tmp, L_mat_tmp_BD;
  arma::field<arma::mat> centers_sample(n_it, 1), log_alpha_sample(n_it, 1);
  arma::cube sigma_cube_it(xdim, xdim, n_cent_init);
  arma::field<arma::cube> sigma_sample(n_it, 1);
  double beta_tmp, var_log_alpha = pow(sd_log_alpha, 2), mean_mu_post, var_mu_post, LM_slice = lung_data.n_rows;
  arma::vec lung_data_idx = arma::linspace(0, LM_slice - 1, LM_slice), BD_probs = arma::ones(LM_slice) / LM_slice;
  arma::vec idx_centers_init, DR_vec(n_cent_init), log_alpha_tmp_BD, PP_lik_cur, PP_lik_prop, L_vec_prop, bd_prob_vec;
  double log_bd_surf = log(BD_probs(0)), DR_tot, prior_n_cent_log = log(prior_n_cent / LM_slice);
  bool bd_flag;
  double bd_vt, BR_tot = 1;
  int n_bd_events, bd_kill_idx;
  struct Node *head = new Node;
  Node *kill_node;
  arma::vec sample_n_cent(n_it), sample_cum_int(n_it);
  //Rcpp::Rcout << "Initialization of variables OK" << std::endl;

  xwin(0) = lung_data.col(0).min();
  xwin(1) = lung_data.col(0).max();
  ywin(0) = lung_data.col(1).min();
  ywin(1) = lung_data.col(1).max();

  idx_centers_init = RcppArmadillo::sample(lung_data_idx, n_cent_init, false, BD_probs);
  //Rcpp::Rcout << "idx_centers_init = " << idx_centers_init << std::endl;
  //Rcpp::Rcout << "Sample for initial center_idx OK" << std::endl;
  mu_alpha_it = R::rnorm(mean_mu_alpha, sqrt(var_mu_alpha));
  initNode(head, lung_data.row(idx_centers_init(0)), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  for(i = 1; i < n_cent_init; i++){
    addNode(head, lung_data.row(idx_centers_init(i)), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  }
  update_Ps(head, n_cent_init, pen_dist, xdim);
  //log_alpha_it = arma::randn(n_cent_init) * sd_log_alpha + mean_mu_alpha;


  sigma_cube_it = get_sigmas(head, n_cent_init, xdim);
  log_alpha_it = get_log_alphas(head, n_cent_init);
  centers_it = get_centers(head, n_cent_init, xdim);

  //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
  //Rcpp::Rcout << "log_alphas_it = " << log_alpha_it << std::endl;
  //Rcpp::Rcout << "sigma_cube_it = " << sigma_cube_it << std::endl;

  sigma_sample(0, 0) = sigma_cube_it;
  log_alpha_sample(0, 0) = log_alpha_it;
  centers_sample(0, 0) = centers_it;
  mu_alpha_sample(0) = mu_alpha_it;
  beta_sample(0) = beta_it;
  sample_n_cent(0) = n_cent_init;
  sample_cum_int(0) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  //Rcpp::Rcout << "Storage of initials OK" << std::endl;
  //L_mat = dmvnrm_vec_arma_cube(obs_points, centers_it, sigma_cube_it);
  //Rcpp::Rcout << "L_mat creation OK" << std::endl;
  PP_lik_cur = calc_PP_lik_vec(head);
  ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
  //Rcpp::Rcout << "ll_cur OK" << std::endl;
  //P_mat = P_mat_gen1(centers_it, pen_dist);
  //Rcpp::Rcout << "P_mat OK" << std::endl;

  //  Initial BD death rates
  DR_vec = calc_DR_vec(head, n_cent_init, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);

  //Rcpp::Rcout << "Initial death rates OK" << std::endl;
  //Rcpp::Rcout << "DR_vec = " << DR_vec << std::endl;

  n_cent_it = n_cent_init;

  //   Birth-death process:
  for(i = 1; i < n_it; ++i){
    bd_flag = true;
    bd_vt = 0;
    n_bd_events = 0;
    while(bd_flag & (n_bd_events < max_bd_events)){
      //Rcpp::Rcout << "n_bd_events = " << n_bd_events << std::endl;
      //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
      //Rcpp::Rcout << "log_alpha_it = " << log_alpha_it << std::endl;
      //display(head);
      DR_tot = arma::sum(DR_vec);

      //Rcpp::Rcout << "DR_tot = " << DR_tot << std::endl;
      //Rcpp::Rcout << "DR_vec = " << DR_tot << std::endl;
      if(DR_vec.has_inf()){
        bd_idx_inf = arma::find_nonfinite(DR_vec);
        //Rcpp::Rcout << "Infinite DR, bd_idx_inf = " << bd_idx_inf << std::endl;
        if(bd_idx_inf.n_elem < 2){
          bd_kill_idx = bd_idx_inf(0);
        }
        else{
          bd_kill_idx = RcppArmadillo::sample(bd_idx_inf, 1, false)(0);
        }

        kill_node = searchNode(head, bd_kill_idx);
        if(n_cent_it < 2){
          PP_lik_cur.zeros();
        }
        else{
          PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
        }
        ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
        deleteNode(&head, kill_node);
        n_bd_events++;
        n_cent_it--;
        //Rcpp::Rcout << "Node killed succesfully (inf)" << std::endl;
      }
      else{
        bd_vt += R::rexp(1.0 / (BR_tot + DR_tot));
        //Rcpp::Rcout << "bd_vt = " << bd_vt << std::endl;
        if(bd_vt < max_bd_vt){
          if(R::runif(0, BR_tot + DR_tot) < BR_tot){
            idx_birth = RcppArmadillo::sample(lung_data_idx, 1, false)(0);
            insertFront(&head, lung_data.row(idx_birth), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
            PP_lik_cur += (head -> alpha) * (head -> mvn_dens);
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            n_bd_events++;
            n_cent_it++;
            //Rcpp::Rcout << "Node birthed succesfully" << std::endl;
          }
          else{
            bd_prob_vec = DR_vec / DR_tot;
            bd_kill_idx = RcppArmadillo::sample(arma::linspace(0, n_cent_it - 1, n_cent_it), 1, false, bd_prob_vec)(0);
            kill_node = searchNode(head, bd_kill_idx);
            if(n_cent_it < 2){
              PP_lik_cur.zeros();
            }
            else{
              PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
            }
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            deleteNode(&head, kill_node);
            n_bd_events++;
            n_cent_it--;
            //Rcpp::Rcout << "Node killed succesfully" << std::endl;
          }
        }
        else{
          bd_flag = false;
        }
      }
      if(bd_flag){
        if(n_cent_it < 1){
          DR_vec.set_size(1);
          DR_vec(0) = 0;
        }
        else{
          update_Ps(head, n_cent_it, pen_dist, xdim);
          //Rcpp::Rcout << "Ps updated" << std::endl;
          DR_vec.set_size(n_cent_it);
          DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
          //Rcpp::Rcout << "DRs updated succesfully" << std::endl;
        }
      }
    }

    //Rcpp::Rcout << "Left BD Process" << std::endl;

    L_mat.set_size(n_points, n_cent_it);
    L_mat = get_L_mat(head, n_cent_it, n_points);

    log_alpha_it.set_size(n_cent_it);
    log_alpha_it = get_log_alphas(head, n_cent_it);

    centers_it.set_size(n_cent_it, xdim);
    centers_it = get_centers(head, n_cent_it, xdim);
    P_mat = P_mat_gen1(centers_it, pen_dist);

    sigma_cube_it.set_size(xdim, xdim, n_cent_it);
    sigma_cube_it = get_sigmas(head, n_cent_it, xdim);

    //PP_lik_cur = calc_PP_lik_vec(head);
    //ll_cur = arma::sum(arma::log(PP_lik_cur + beta_it));
    for(j = 0; j < n_cent_it; j++){
      //log_alpha_temp = log_alpha_it;
      log_alpha_prop = R::rnorm(log_alpha_it(j), sd_prop_alpha);
      //log_alpha_temp(j) = log_alpha_prop;
      PP_lik_prop = PP_lik_cur + (exp(log_alpha_prop) - exp(log_alpha_it(j))) * L_mat.col(j);
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop - ll_cur -
        exp(log_alpha_prop) + exp(log_alpha_it(j)) +
        (1.0/(2.0 * pow(sd_log_alpha, 2))) * (pow(log_alpha_it(j) - mean_mu_alpha, 2) - pow(log_alpha_prop - mean_mu_alpha, 2));
      if(R::runif(0, 1) < exp(mh_log)){
        //Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        log_alpha_it(j) = log_alpha_prop;
      }
    }

    // Doing center location MH steps
    for(j = 0; j < n_cent_it; j++){
      center_prop = centers_it.row(j);
      center_prop(0) += R::runif(-window_hw, window_hw);
      center_prop(1) += R::runif(-window_hw, window_hw);
      if(center_prop(0) < xwin(1) & center_prop(0) > xwin(0) &
         center_prop(1) < ywin(1) & center_prop(1) > ywin(0)){
        //L_mat_tmp = L_mat;
        //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        L_vec_prop = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        P_vec_tmp = P_mat_gen2(centers_it, center_prop, pen_dist);
        //P_vec_tmp.shed_row(j);
        PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
        ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
        mh_log = ll_prop - ll_cur + (arma::sum(P_vec_tmp) - P_vec_tmp(j) - (arma::sum(P_mat.col(j)) - 1)) * log_pen_val;
        if(R::runif(0, 1) < exp(mh_log)){
          PP_lik_cur = PP_lik_prop;
          ll_cur = ll_prop;
          centers_it.row(j) = center_prop;
          L_mat.col(j) = L_vec_prop;
          P_mat = P_mat_gen1(centers_it, pen_dist);
        }
      }
    }

    // Doing sigma matrix MH steps
    for(j = 0; j < n_cent_it; j++){
      sigma_tmp = riwish_arma(df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j));
      /*
       log(diwish(W = sigma_mat_it[[j]], v = df_inv_wish_prior,
       S = (df_inv_wish_prior - 3) * S_mat_prior)) + ## Prior
       log(diwish(W = sigma_prop_j, v = df_inv_wish_prop,
       S = (df_inv_wish_prop - 3) * sigma_mat_it[[j]]))
       */
      //L_mat_tmp = L_mat;
      //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      L_vec_prop = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop +
        log(diwish_arma(sigma_tmp, df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) +
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prop, (df_iw_prop - 3.0) * sigma_tmp)) -
        ll_cur -
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) -
        log(diwish_arma(sigma_tmp, df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j)));
      if(R::runif(0, 1) < exp(mh_log)){
        ////Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        L_mat.col(j) = L_vec_prop;
        sigma_cube_it.slice(j) = sigma_tmp;
        //L_mat = L_mat_tmp;
      }
    }

    //  Updating beta
    beta_tmp =  std::fabs(beta_it + R::runif(-.0005, .0005));
    ll_prop = arma::sum((arma::log(PP_lik_cur + beta_tmp)));
    mh_log = ll_prop - beta_tmp * LM_slice + R::dgamma(beta_tmp, 0.01, 1.0 / 0.01, 1) -
      (ll_cur - beta_it * LM_slice + R::dgamma(beta_it, 0.01, 1.0 / 0.01, 1));
    if(R::runif(0, 1) < exp(mh_log)){
      ////Rcpp::Rcout <<"MH accept" << std::endl;
      ll_cur = ll_prop;
      beta_it = beta_tmp;
    }

    //  Updating mu_alpha

    mean_mu_post = ((mean_mu_alpha / var_mu_alpha) + (arma::sum(log_alpha_it) / var_log_alpha)) /
      ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    var_mu_post = 1.0 / ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    mu_alpha_it = R::rnorm(mean_mu_post, sqrt(var_mu_post));

    //Putting new values in the linked list
    update_node_values(head, centers_it, log_alpha_it, sigma_cube_it, L_mat);

    //  Updating death rates
    if(n_cent_it > 0){
      DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
    }

    // Storing MCMC samples
    centers_sample(i, 0) = centers_it;
    log_alpha_sample(i, 0) = log_alpha_it.t();
    sigma_sample(i, 0) = sigma_cube_it;
    beta_sample(i) = beta_it;
    mu_alpha_sample(i) = mu_alpha_it;
    sample_n_cent(i) = n_cent_it;
    sample_cum_int(i) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  }
  return Rcpp::List::create(Rcpp::Named("log_alpha_sample") = log_alpha_sample,
                            Rcpp::Named("centers_sample") = centers_sample,
                            Rcpp::Named("sigmas_sample") = sigma_sample,
                            Rcpp::Named("mu_alpha_sample") = mu_alpha_sample,
                            Rcpp::Named("beta_sample") = beta_sample,
                            Rcpp::Named("n_centers_sample") = sample_n_cent,
                            Rcpp::Named("cumulative_intensity_sample") = sample_cum_int);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Doing a version for continuous points (i.e. not pixelated)
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//' Bayesian SNCP fit using BD-MCMC
//'
//' Run a BD-MCMC chain for a SNCP (mvn dispersion density) on lung ct data using a uniform proposal surface for the BD process
//'
//'
//'
//'
//'
//' @param obs_points matrix of coordinates for observed points
//' @param mean_mu_alpha prior mean for mu_alpha
//' @param sd_log_alpha prior std. dev. for log_alphas
//' @param sd_prop_alpha std. dev. for random walk proposal for log_alphas
//' @param beta prior estimate of homogeneous Poisson intensity function for cluster centers
//' @param n_it number of iterations to run MCMC
//' @param window_hw half-width of random walk cube for proposing new cluster center locations
//' @param df_iw_prior degrees of freedom for inverse-wishart prior on dispersion matrices
//' @param df_iw_prop degrees of freedom for inverse-wishart proposal on dispersion matrices
//' @param sigma_prior prior for dispersion matrices sigma
//' @param xwin vector that has the min and max x values of the observation window
//' @param ywin vector that has the min and max y values of the observation window
//' @param var_mu_alpha prior variance on mu_alpha
//' @param pen_dist distance for strauss process repulsion
//' @param pen_val vaule between 0 and 1 that is the penalty for 2 cluster centers being within pen_dist of each other
//' @param n_cent_init number of clusters to start with
//' @param prior_n_cent prior value for number of clusters
//' @param max_bd_events max events (births + deaths) to allow at each iteration of BD-MCMC
//' @param max_bd_vt max ammount of virtual time to spend in BD process at each BD_MCMC iteration
//'
//' @author Brian Vestal
//'
//' @return
//' Returns a named list containing posterior samples of model parameters with the following elements:
//' \itemize{
//' \item log_alpha_sample = list whose elements are vectors containing the estimated log-alpha estimates for the clusters present in that iteration
//' \item centers_sample = list whose elements matrices containing the estimated center locations for the clusters present in that iteration
//' \item sigmas_sample = list whose elements are lists containing the estimated covariance matrix estimates for the clusters present in that iteration
//' \item mu_alpha_sample = vector of mu_alpha samples
//' \item eta_sample = vector of eta samples,
//' \item n_centers_sample = vector that contains the number of estimated clusters at each iteration,
//' \item cumulative_intensity_sample = vector that contains the cumulative intensity estimate for the SNCP at each iteration
//' }
//'
//' @export
// [[Rcpp::export]]

Rcpp::List sncp_bdmcmc_cont(arma::mat obs_points,
                            double mean_mu_alpha,
                            double sd_log_alpha,
                            double sd_prop_alpha,
                            double beta,
                            int n_it,
                            double window_hw,
                            int df_iw_prior,
                            int df_iw_prop,
                            arma::mat sigma_prior,
                            arma::mat obs_window,
                            double LM,
                            double var_mu_alpha,
                            double pen_dist,
                            double pen_val,
                            int n_cent_init,
                            double prior_n_cent,
                            int max_bd_events,
                            double max_bd_vt){
  int i, j, xdim = obs_points.n_cols, n_cent_it, n_points = obs_points.n_rows;
  double ll_cur, ll_prop, mh_log, log_alpha_prop, beta_it = beta, mu_alpha_it, log_pen_val = log(pen_val);
  arma::vec log_alpha_it, beta_sample(n_it), mu_alpha_sample(n_it);
  //arma::mat log_alpha_sample(n_it, n_cent);
  arma::umat P_mat;
  arma::uvec P_vec_tmp, bd_idx_inf;
  arma::rowvec center_prop;
  arma::mat centers_it(n_cent_init, xdim), L_mat, L_mat_tmp, centers_tmp, sigma_tmp, L_mat_tmp_BD;
  arma::field<arma::mat> centers_sample(n_it, 1), log_alpha_sample(n_it, 1);
  arma::cube sigma_cube_it(xdim, xdim, n_cent_init);
  arma::field<arma::cube> sigma_sample(n_it, 1);
  double beta_tmp, var_log_alpha = pow(sd_log_alpha, 2), mean_mu_post, var_mu_post, LM_slice = LM;
  arma::vec idx_centers_init, DR_vec(n_cent_init), log_alpha_tmp_BD, PP_lik_cur, PP_lik_prop, L_vec_prop, bd_prob_vec;
  double log_bd_surf = log(1.0 / LM_slice), DR_tot, prior_n_cent_log = log(prior_n_cent / LM_slice);
  bool bd_flag;
  double bd_vt, BR_tot = 1;
  int n_bd_events, bd_kill_idx;
  struct Node *head = new Node;
  Node *kill_node;
  arma::vec sample_n_cent(n_it), sample_cum_int(n_it);
  arma::rowvec cent_add(xdim);
  arma::vec bd_event_vec(n_it), vt_vec(n_it);
  int kk;
  // Rcpp::Rcout << "Initialization of variables OK" << std::endl;

  // xwin(0) = lung_data.col(0).min();
  // xwin(1) = lung_data.col(0).max();
  // ywin(0) = lung_data.col(1).min();
  // ywin(1) = lung_data.col(1).max();
  // double x_range = xwin(1) - xwin(0), y_range = ywin(1) - ywin(0);
  arma::mat centers_init(n_cent_init, xdim);
  for(int jj = 0; jj < n_cent_init; jj++){
    for(kk = 0; kk < xdim; kk++){
      centers_init(jj, kk) = R::runif(obs_window(kk, 0), obs_window(kk, 1));
    }
    // centers_init(jj, 1) = R::runif(ywin(0), ywin(1));
  }
  // idx_centers_init = RcppArmadillo::sample(lung_data_idx, n_cent_init, false, BD_probs);
  //Rcpp::Rcout << "idx_centers_init = " << idx_centers_init << std::endl;
  // Rcpp::Rcout << "Sample for initial center_idx OK" << std::endl;
  mu_alpha_it = R::rnorm(mean_mu_alpha, sqrt(var_mu_alpha));
  initNode(head, centers_init.row(0), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  for(i = 1; i < n_cent_init; i++){
    addNode(head, centers_init.row(i), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  }
  update_Ps(head, n_cent_init, pen_dist, xdim);
  //log_alpha_it = arma::randn(n_cent_init) * sd_log_alpha + mean_mu_alpha;


  sigma_cube_it = get_sigmas(head, n_cent_init, xdim);
  log_alpha_it = get_log_alphas(head, n_cent_init);
  centers_it = get_centers(head, n_cent_init, xdim);

  //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
  //Rcpp::Rcout << "log_alphas_it = " << log_alpha_it << std::endl;
  //Rcpp::Rcout << "sigma_cube_it = " << sigma_cube_it << std::endl;

  sigma_sample(0, 0) = sigma_cube_it;
  log_alpha_sample(0, 0) = log_alpha_it;
  centers_sample(0, 0) = centers_it;
  mu_alpha_sample(0) = mu_alpha_it;
  beta_sample(0) = beta_it;
  sample_n_cent(0) = n_cent_init;
  sample_cum_int(0) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  // Rcpp::Rcout << "Storage of initials OK" << std::endl;
  //L_mat = dmvnrm_vec_arma_cube(obs_points, centers_it, sigma_cube_it);
  //Rcpp::Rcout << "L_mat creation OK" << std::endl;
  PP_lik_cur = calc_PP_lik_vec(head);
  ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
  //Rcpp::Rcout << "ll_cur OK" << std::endl;
  //P_mat = P_mat_gen1(centers_it, pen_dist);
  //Rcpp::Rcout << "P_mat OK" << std::endl;

  //  Initial BD death rates
  DR_vec = calc_DR_vec(head, n_cent_init, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);

  // Rcpp::Rcout << "Initial death rates OK" << std::endl;
  //Rcpp::Rcout << "DR_vec = " << DR_vec << std::endl;

  n_cent_it = n_cent_init;

  //   Birth-death process:
  for(i = 1; i < n_it; ++i){
    bd_flag = true;
    bd_vt = 0;
    n_bd_events = 0;
    while(bd_flag & (n_bd_events < max_bd_events)){
      //Rcpp::Rcout << "n_bd_events = " << n_bd_events << std::endl;
      //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
      //Rcpp::Rcout << "log_alpha_it = " << log_alpha_it << std::endl;
      //display(head);
      DR_tot = arma::sum(DR_vec);

      //Rcpp::Rcout << "DR_tot = " << DR_tot << std::endl;
      //Rcpp::Rcout << "DR_vec = " << DR_tot << std::endl;
      if(DR_vec.has_inf()){
        bd_idx_inf = arma::find_nonfinite(DR_vec);
        //Rcpp::Rcout << "Infinite DR, bd_idx_inf = " << bd_idx_inf << std::endl;
        if(bd_idx_inf.n_elem < 2){
          bd_kill_idx = bd_idx_inf(0);
        }
        else{
          bd_kill_idx = RcppArmadillo::sample(bd_idx_inf, 1, false)(0);
        }

        kill_node = searchNode(head, bd_kill_idx);
        if(n_cent_it < 2){
          PP_lik_cur.zeros();
        }
        else{
          PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
        }
        ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
        deleteNode(&head, kill_node);
        n_bd_events++;
        n_cent_it--;
        //Rcpp::Rcout << "Node killed succesfully (inf)" << std::endl;
      }
      else{
        bd_vt += R::rexp(1.0 / (BR_tot + DR_tot));
        //Rcpp::Rcout << "bd_vt = " << bd_vt << std::endl;
        if(bd_vt < max_bd_vt){
          if(R::runif(0, BR_tot + DR_tot) < BR_tot){
            // idx_birth = RcppArmadillo::sample(lung_data_idx, 1, false)(0);
            for(kk = 0; kk < xdim; kk++){
              cent_add(kk) = R::runif(obs_window(kk, 0), obs_window(kk, 1));
            }
            // cent_add(0) = R::runif(xwin(0), xwin(1));
            // cent_add(1) = R::runif(ywin(0), ywin(1));
            insertFront(&head, cent_add, mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
            PP_lik_cur += (head -> alpha) * (head -> mvn_dens);
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            n_bd_events++;
            n_cent_it++;
            //Rcpp::Rcout << "Node birthed succesfully" << std::endl;
          }
          else{
            bd_prob_vec = DR_vec / DR_tot;
            bd_kill_idx = RcppArmadillo::sample(arma::linspace(0, n_cent_it - 1, n_cent_it), 1, false, bd_prob_vec)(0);
            kill_node = searchNode(head, bd_kill_idx);
            if(n_cent_it < 2){
              PP_lik_cur.zeros();
            }
            else{
              PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
            }
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            deleteNode(&head, kill_node);
            n_bd_events++;
            n_cent_it--;
            //Rcpp::Rcout << "Node killed succesfully" << std::endl;
          }
        }
        else{
          bd_flag = false;
        }
      }
      if(bd_flag){
        if(n_cent_it < 1){
          DR_vec.set_size(1);
          DR_vec(0) = 0;
        }
        else{
          update_Ps(head, n_cent_it, pen_dist, xdim);
          //Rcpp::Rcout << "Ps updated" << std::endl;
          DR_vec.set_size(n_cent_it);
          DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
          //Rcpp::Rcout << "DRs updated succesfully" << std::endl;
        }
      }
    }

    // Rcpp::Rcout << "Left BD Process" << std::endl;
    bd_event_vec(i) = n_bd_events;
    vt_vec(i) = bd_vt;
    L_mat.set_size(n_points, n_cent_it);
    L_mat = get_L_mat(head, n_cent_it, n_points);

    log_alpha_it.set_size(n_cent_it);
    log_alpha_it = get_log_alphas(head, n_cent_it);

    centers_it.set_size(n_cent_it, xdim);
    centers_it = get_centers(head, n_cent_it, xdim);
    P_mat = P_mat_gen1(centers_it, pen_dist);

    sigma_cube_it.set_size(xdim, xdim, n_cent_it);
    sigma_cube_it = get_sigmas(head, n_cent_it, xdim);

    //PP_lik_cur = calc_PP_lik_vec(head);
    //ll_cur = arma::sum(arma::log(PP_lik_cur + beta_it));
    for(j = 0; j < n_cent_it; j++){
      //log_alpha_temp = log_alpha_it;
      log_alpha_prop = R::rnorm(log_alpha_it(j), sd_prop_alpha);
      //log_alpha_temp(j) = log_alpha_prop;
      PP_lik_prop = PP_lik_cur + (exp(log_alpha_prop) - exp(log_alpha_it(j))) * L_mat.col(j);
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop - ll_cur -
        exp(log_alpha_prop) + exp(log_alpha_it(j)) +
        (1.0/(2.0 * pow(sd_log_alpha, 2))) * (pow(log_alpha_it(j) - mean_mu_alpha, 2) - pow(log_alpha_prop - mean_mu_alpha, 2));
      if(R::runif(0, 1) < exp(mh_log)){
        //Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        log_alpha_it(j) = log_alpha_prop;
      }
    }

    // Doing center location MH steps
    for(j = 0; j < n_cent_it; j++){
      center_prop = centers_it.row(j);
      for(kk = 0; kk < xdim; kk++){
        center_prop(kk) += R::runif(-window_hw, window_hw);
      }
      // center_prop(0) += R::runif(-window_hw, window_hw);
      // center_prop(1) += R::runif(-window_hw, window_hw);
      // if(center_prop(0) < xwin(1) & center_prop(0) > xwin(0) &
      //    center_prop(1) < ywin(1) & center_prop(1) > ywin(0)){
      //   //L_mat_tmp = L_mat;
      //   //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
      //   L_vec_prop = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
      //   P_vec_tmp = P_mat_gen2(centers_it, center_prop, pen_dist);
      //   //P_vec_tmp.shed_row(j);
      //   PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
      //   ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      //   mh_log = ll_prop - ll_cur + (arma::sum(P_vec_tmp) - P_vec_tmp(j) - (arma::sum(P_mat.col(j)) - 1)) * log_pen_val;
      //   if(R::runif(0, 1) < exp(mh_log)){
      //     PP_lik_cur = PP_lik_prop;
      //     ll_cur = ll_prop;
      //     centers_it.row(j) = center_prop;
      //     L_mat.col(j) = L_vec_prop;
      //     P_mat = P_mat_gen1(centers_it, pen_dist);
      //   }
      // }

      //L_mat_tmp = L_mat;
      //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
      L_vec_prop = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
      P_vec_tmp = P_mat_gen2(centers_it, center_prop, pen_dist);
      //P_vec_tmp.shed_row(j);
      PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop - ll_cur + (arma::sum(P_vec_tmp) - P_vec_tmp(j) - (arma::sum(P_mat.col(j)) - 1)) * log_pen_val;
      if(R::runif(0, 1) < exp(mh_log)){
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        centers_it.row(j) = center_prop;
        L_mat.col(j) = L_vec_prop;
        P_mat = P_mat_gen1(centers_it, pen_dist);
      }

    }

    // Doing sigma matrix MH steps
    for(j = 0; j < n_cent_it; j++){
      sigma_tmp = riwish_arma(df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j));
      /*
       log(diwish(W = sigma_mat_it[[j]], v = df_inv_wish_prior,
       S = (df_inv_wish_prior - 3) * S_mat_prior)) + ## Prior
       log(diwish(W = sigma_prop_j, v = df_inv_wish_prop,
       S = (df_inv_wish_prop - 3) * sigma_mat_it[[j]]))
       */
      //L_mat_tmp = L_mat;
      //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      L_vec_prop = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop +
        log(diwish_arma(sigma_tmp, df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) +
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prop, (df_iw_prop - 3.0) * sigma_tmp)) -
        ll_cur -
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) -
        log(diwish_arma(sigma_tmp, df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j)));
      if(R::runif(0, 1) < exp(mh_log)){
        ////Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        L_mat.col(j) = L_vec_prop;
        sigma_cube_it.slice(j) = sigma_tmp;
        //L_mat = L_mat_tmp;
      }
    }

    //  Updating beta
    beta_tmp =  std::fabs(beta_it + R::runif(-.0005, .0005));
    ll_prop = arma::sum((arma::log(PP_lik_cur + beta_tmp)));
    mh_log = ll_prop - beta_tmp * LM_slice + R::dgamma(beta_tmp, 0.01, 1.0 / 0.01, 1) -
      (ll_cur - beta_it * LM_slice + R::dgamma(beta_it, 0.01, 1.0 / 0.01, 1));
    if(R::runif(0, 1) < exp(mh_log)){
      ////Rcpp::Rcout <<"MH accept" << std::endl;
      ll_cur = ll_prop;
      beta_it = beta_tmp;
    }

    //  Updating mu_alpha

    mean_mu_post = ((mean_mu_alpha / var_mu_alpha) + (arma::sum(log_alpha_it) / var_log_alpha)) /
      ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    var_mu_post = 1.0 / ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    mu_alpha_it = R::rnorm(mean_mu_post, sqrt(var_mu_post));

    //Putting new values in the linked list
    update_node_values(head, centers_it, log_alpha_it, sigma_cube_it, L_mat);

    //  Updating death rates
    if(n_cent_it > 0){
      DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
    }

    // Storing MCMC samples
    centers_sample(i, 0) = centers_it;
    log_alpha_sample(i, 0) = log_alpha_it.t();
    sigma_sample(i, 0) = sigma_cube_it;
    beta_sample(i) = beta_it;
    mu_alpha_sample(i) = mu_alpha_it;
    sample_n_cent(i) = n_cent_it;
    sample_cum_int(i) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  }
  return Rcpp::List::create(Rcpp::Named("log_alpha_sample") = log_alpha_sample,
                            Rcpp::Named("centers_sample") = centers_sample,
                            Rcpp::Named("sigmas_sample") = sigma_sample,
                            Rcpp::Named("mu_alpha_sample") = mu_alpha_sample,
                            Rcpp::Named("beta_sample") = beta_sample,
                            Rcpp::Named("n_centers_sample") = sample_n_cent,
                            Rcpp::Named("cumulative_intensity_sample") = sample_cum_int,
                            Rcpp::Named("n_bd_events") = bd_event_vec,
                            Rcpp::Named("bd_vt") = vt_vec);
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Doing a version for continuous points (i.e. not pixelated)
//  Now with no noise process (assuming it has already been filtered)
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//' Bayesian SNCP fit using BD-MCMC
//'
//' Run a BD-MCMC chain for a SNCP (mvn dispersion density) on lung ct data using a uniform proposal surface for the BD process
//'
//'
//'
//'
//'
//' @param obs_points matrix of coordinates for observed points
//' @param mean_mu_alpha prior mean for mu_alpha
//' @param sd_log_alpha prior std. dev. for log_alphas
//' @param sd_prop_alpha std. dev. for random walk proposal for log_alphas
//' @param n_it number of iterations to run MCMC
//' @param window_hw half-width of random walk cube for proposing new cluster center locations
//' @param df_iw_prior degrees of freedom for inverse-wishart prior on dispersion matrices
//' @param df_iw_prop degrees of freedom for inverse-wishart proposal on dispersion matrices
//' @param sigma_prior prior for dispersion matrices sigma
//' @param xwin vector that has the min and max x values of the observation window
//' @param ywin vector that has the min and max y values of the observation window
//' @param var_mu_alpha prior variance on mu_alpha
//' @param pen_dist distance for strauss process repulsion
//' @param pen_val vaule between 0 and 1 that is the penalty for 2 cluster centers being within pen_dist of each other
//' @param n_cent_init number of clusters to start with
//' @param prior_n_cent prior value for number of clusters
//' @param max_bd_events max events (births + deaths) to allow at each iteration of BD-MCMC
//' @param max_bd_vt max ammount of virtual time to spend in BD process at each BD_MCMC iteration
//'
//' @author Brian Vestal
//'
//' @return
//' Returns a named list containing posterior samples of model parameters with the following elements:
//' \itemize{
//' \item log_alpha_sample = list whose elements are vectors containing the estimated log-alpha estimates for the clusters present in that iteration
//' \item centers_sample = list whose elements matrices containing the estimated center locations for the clusters present in that iteration
//' \item sigmas_sample = list whose elements are lists containing the estimated covariance matrix estimates for the clusters present in that iteration
//' \item mu_alpha_sample = vector of mu_alpha samples
//' \item n_centers_sample = vector that contains the number of estimated clusters at each iteration,
//' \item cumulative_intensity_sample = vector that contains the cumulative intensity estimate for the SNCP at each iteration
//' }
//'
//' @export
// [[Rcpp::export]]

Rcpp::List sncp_bdmcmc_cont_no_noise(arma::mat obs_points,
                                     double mean_mu_alpha,
                                     double sd_log_alpha,
                                     double sd_prop_alpha,
                                     int n_it,
                                     double window_hw,
                                     int df_iw_prior,
                                     int df_iw_prop,
                                     arma::mat sigma_prior,
                                     arma::vec xwin,
                                     arma::vec ywin,
                                     double var_mu_alpha,
                                     double pen_dist,
                                     double pen_val,
                                     int n_cent_init,
                                     double prior_n_cent,
                                     int max_bd_events,
                                     double max_bd_vt){
  int i, j, xdim = obs_points.n_cols, n_cent_it, n_points = obs_points.n_rows;
  double ll_cur, ll_prop, mh_log, log_alpha_prop, beta_it = 0.000000001, mu_alpha_it, log_pen_val = log(pen_val);
  arma::vec log_alpha_it, mu_alpha_sample(n_it);
  //arma::mat log_alpha_sample(n_it, n_cent);
  arma::umat P_mat;
  arma::uvec P_vec_tmp, bd_idx_inf;
  arma::rowvec center_prop;
  arma::mat centers_it(n_cent_init, xdim), L_mat, L_mat_tmp, centers_tmp, sigma_tmp, L_mat_tmp_BD;
  arma::field<arma::mat> centers_sample(n_it, 1), log_alpha_sample(n_it, 1);
  arma::cube sigma_cube_it(xdim, xdim, n_cent_init);
  arma::field<arma::cube> sigma_sample(n_it, 1);
  double var_log_alpha = pow(sd_log_alpha, 2), mean_mu_post, var_mu_post, LM_slice = (xwin(1) - xwin(0)) * (ywin(1) - ywin(0));
  arma::vec idx_centers_init, DR_vec(n_cent_init), log_alpha_tmp_BD, PP_lik_cur, PP_lik_prop, L_vec_prop, bd_prob_vec;
  double log_bd_surf = log(1.0 / LM_slice), DR_tot, prior_n_cent_log = log(prior_n_cent / LM_slice);
  bool bd_flag;
  double bd_vt, BR_tot = 1;
  int n_bd_events, bd_kill_idx;
  struct Node *head = new Node;
  Node *kill_node;
  arma::vec sample_n_cent(n_it), sample_cum_int(n_it);
  arma::rowvec cent_add(2);
  arma::vec bd_event_vec(n_it), vt_vec(n_it);
  //Rcpp::Rcout << "Initialization of variables OK" << std::endl;

  // xwin(0) = lung_data.col(0).min();
  // xwin(1) = lung_data.col(0).max();
  // ywin(0) = lung_data.col(1).min();
  // ywin(1) = lung_data.col(1).max();
  // double x_range = xwin(1) - xwin(0), y_range = ywin(1) - ywin(0);
  arma::mat centers_init(n_cent_init, xdim);
  for(int jj = 0; jj < n_cent_init; jj++){
    centers_init(jj, 0) = R::runif(xwin(0), xwin(1));
    centers_init(jj, 1) = R::runif(ywin(0), ywin(1));
  }
  // idx_centers_init = RcppArmadillo::sample(lung_data_idx, n_cent_init, false, BD_probs);
  //Rcpp::Rcout << "idx_centers_init = " << idx_centers_init << std::endl;
  //Rcpp::Rcout << "Sample for initial center_idx OK" << std::endl;
  mu_alpha_it = R::rnorm(mean_mu_alpha, sqrt(var_mu_alpha));
  initNode(head, centers_init.row(0), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  for(i = 1; i < n_cent_init; i++){
    addNode(head, centers_init.row(i), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  }
  update_Ps(head, n_cent_init, pen_dist, xdim);
  //log_alpha_it = arma::randn(n_cent_init) * sd_log_alpha + mean_mu_alpha;


  sigma_cube_it = get_sigmas(head, n_cent_init, xdim);
  log_alpha_it = get_log_alphas(head, n_cent_init);
  centers_it = get_centers(head, n_cent_init, xdim);

  //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
  //Rcpp::Rcout << "log_alphas_it = " << log_alpha_it << std::endl;
  //Rcpp::Rcout << "sigma_cube_it = " << sigma_cube_it << std::endl;

  sigma_sample(0, 0) = sigma_cube_it;
  log_alpha_sample(0, 0) = log_alpha_it;
  centers_sample(0, 0) = centers_it;
  mu_alpha_sample(0) = mu_alpha_it;
  sample_n_cent(0) = n_cent_init;
  sample_cum_int(0) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  //Rcpp::Rcout << "Storage of initials OK" << std::endl;
  //L_mat = dmvnrm_vec_arma_cube(obs_points, centers_it, sigma_cube_it);
  //Rcpp::Rcout << "L_mat creation OK" << std::endl;
  PP_lik_cur = calc_PP_lik_vec(head);
  ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
  //Rcpp::Rcout << "ll_cur OK" << std::endl;
  //P_mat = P_mat_gen1(centers_it, pen_dist);
  //Rcpp::Rcout << "P_mat OK" << std::endl;

  //  Initial BD death rates
  DR_vec = calc_DR_vec(head, n_cent_init, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);

  //Rcpp::Rcout << "Initial death rates OK" << std::endl;
  //Rcpp::Rcout << "DR_vec = " << DR_vec << std::endl;

  n_cent_it = n_cent_init;

  //   Birth-death process:
  for(i = 1; i < n_it; ++i){
    bd_flag = true;
    bd_vt = 0;
    n_bd_events = 0;
    while(bd_flag & (n_bd_events < max_bd_events)){
      //Rcpp::Rcout << "n_bd_events = " << n_bd_events << std::endl;
      //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
      //Rcpp::Rcout << "log_alpha_it = " << log_alpha_it << std::endl;
      //display(head);
      DR_tot = arma::sum(DR_vec);

      //Rcpp::Rcout << "DR_tot = " << DR_tot << std::endl;
      //Rcpp::Rcout << "DR_vec = " << DR_tot << std::endl;
      if(DR_vec.has_inf()){
        bd_idx_inf = arma::find_nonfinite(DR_vec);
        //Rcpp::Rcout << "Infinite DR, bd_idx_inf = " << bd_idx_inf << std::endl;
        if(bd_idx_inf.n_elem < 2){
          bd_kill_idx = bd_idx_inf(0);
        }
        else{
          bd_kill_idx = RcppArmadillo::sample(bd_idx_inf, 1, false)(0);
        }

        kill_node = searchNode(head, bd_kill_idx);
        if(n_cent_it < 2){
          PP_lik_cur.zeros();
        }
        else{
          PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
        }
        ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
        deleteNode(&head, kill_node);
        n_bd_events++;
        n_cent_it--;
        //Rcpp::Rcout << "Node killed succesfully (inf)" << std::endl;
      }
      else{
        bd_vt += R::rexp(1.0 / (BR_tot + DR_tot));
        //Rcpp::Rcout << "bd_vt = " << bd_vt << std::endl;
        if(bd_vt < max_bd_vt){
          if(R::runif(0, BR_tot + DR_tot) < BR_tot){
            // idx_birth = RcppArmadillo::sample(lung_data_idx, 1, false)(0);
            cent_add(0) = R::runif(xwin(0), xwin(1));
            cent_add(1) = R::runif(ywin(0), ywin(1));
            insertFront(&head, cent_add, mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
            PP_lik_cur += (head -> alpha) * (head -> mvn_dens);
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            n_bd_events++;
            n_cent_it++;
            //Rcpp::Rcout << "Node birthed succesfully" << std::endl;
          }
          else{
            bd_prob_vec = DR_vec / DR_tot;
            bd_kill_idx = RcppArmadillo::sample(arma::linspace(0, n_cent_it - 1, n_cent_it), 1, false, bd_prob_vec)(0);
            kill_node = searchNode(head, bd_kill_idx);
            if(n_cent_it < 2){
              PP_lik_cur.zeros();
            }
            else{
              PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
            }
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            deleteNode(&head, kill_node);
            n_bd_events++;
            n_cent_it--;
            //Rcpp::Rcout << "Node killed succesfully" << std::endl;
          }
        }
        else{
          bd_flag = false;
        }
      }
      if(bd_flag){
        if(n_cent_it < 1){
          DR_vec.set_size(1);
          DR_vec(0) = 0;
        }
        else{
          update_Ps(head, n_cent_it, pen_dist, xdim);
          //Rcpp::Rcout << "Ps updated" << std::endl;
          DR_vec.set_size(n_cent_it);
          DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
          //Rcpp::Rcout << "DRs updated succesfully" << std::endl;
        }
      }
    }
    bd_event_vec(i) = n_bd_events;
    vt_vec(i) = bd_vt;
    //Rcpp::Rcout << "Left BD Process" << std::endl;

    L_mat.set_size(n_points, n_cent_it);
    L_mat = get_L_mat(head, n_cent_it, n_points);

    log_alpha_it.set_size(n_cent_it);
    log_alpha_it = get_log_alphas(head, n_cent_it);

    centers_it.set_size(n_cent_it, xdim);
    centers_it = get_centers(head, n_cent_it, xdim);
    P_mat = P_mat_gen1(centers_it, pen_dist);

    sigma_cube_it.set_size(xdim, xdim, n_cent_it);
    sigma_cube_it = get_sigmas(head, n_cent_it, xdim);

    //PP_lik_cur = calc_PP_lik_vec(head);
    //ll_cur = arma::sum(arma::log(PP_lik_cur + beta_it));
    for(j = 0; j < n_cent_it; j++){
      //log_alpha_temp = log_alpha_it;
      log_alpha_prop = R::rnorm(log_alpha_it(j), sd_prop_alpha);
      //log_alpha_temp(j) = log_alpha_prop;
      PP_lik_prop = PP_lik_cur + (exp(log_alpha_prop) - exp(log_alpha_it(j))) * L_mat.col(j);
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop - ll_cur -
        exp(log_alpha_prop) + exp(log_alpha_it(j)) +
        (1.0/(2.0 * pow(sd_log_alpha, 2))) * (pow(log_alpha_it(j) - mean_mu_alpha, 2) - pow(log_alpha_prop - mean_mu_alpha, 2));
      if(R::runif(0, 1) < exp(mh_log)){
        //Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        log_alpha_it(j) = log_alpha_prop;
      }
    }

    // Doing center location MH steps
    for(j = 0; j < n_cent_it; j++){
      center_prop = centers_it.row(j);
      center_prop(0) += R::runif(-window_hw, window_hw);
      center_prop(1) += R::runif(-window_hw, window_hw);
      if(center_prop(0) < xwin(1) & center_prop(0) > xwin(0) &
         center_prop(1) < ywin(1) & center_prop(1) > ywin(0)){
        //L_mat_tmp = L_mat;
        //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        L_vec_prop = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        P_vec_tmp = P_mat_gen2(centers_it, center_prop, pen_dist);
        //P_vec_tmp.shed_row(j);
        PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
        ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
        mh_log = ll_prop - ll_cur + (arma::sum(P_vec_tmp) - P_vec_tmp(j) - (arma::sum(P_mat.col(j)) - 1)) * log_pen_val;
        if(R::runif(0, 1) < exp(mh_log)){
          PP_lik_cur = PP_lik_prop;
          ll_cur = ll_prop;
          centers_it.row(j) = center_prop;
          L_mat.col(j) = L_vec_prop;
          P_mat = P_mat_gen1(centers_it, pen_dist);
        }
      }
    }

    // Doing sigma matrix MH steps
    for(j = 0; j < n_cent_it; j++){
      sigma_tmp = riwish_arma(df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j));
      /*
       log(diwish(W = sigma_mat_it[[j]], v = df_inv_wish_prior,
       S = (df_inv_wish_prior - 3) * S_mat_prior)) + ## Prior
       log(diwish(W = sigma_prop_j, v = df_inv_wish_prop,
       S = (df_inv_wish_prop - 3) * sigma_mat_it[[j]]))
       */
      //L_mat_tmp = L_mat;
      //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      L_vec_prop = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop +
        log(diwish_arma(sigma_tmp, df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) +
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prop, (df_iw_prop - 3.0) * sigma_tmp)) -
        ll_cur -
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) -
        log(diwish_arma(sigma_tmp, df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j)));
      if(R::runif(0, 1) < exp(mh_log)){
        ////Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        L_mat.col(j) = L_vec_prop;
        sigma_cube_it.slice(j) = sigma_tmp;
        //L_mat = L_mat_tmp;
      }
    }

    //  Updating beta
    // beta_tmp =  std::fabs(beta_it + R::runif(-.0005, .0005));
    // ll_prop = arma::sum((arma::log(PP_lik_cur + beta_tmp)));
    // mh_log = ll_prop - beta_tmp * LM_slice + R::dgamma(beta_tmp, 0.01, 1.0 / 0.01, 1) -
    //   (ll_cur - beta_it * LM_slice + R::dgamma(beta_it, 0.01, 1.0 / 0.01, 1));
    // if(R::runif(0, 1) < exp(mh_log)){
    //   ////Rcpp::Rcout <<"MH accept" << std::endl;
    //   ll_cur = ll_prop;
    //   beta_it = beta_tmp;
    // }

    //  Updating mu_alpha

    mean_mu_post = ((mean_mu_alpha / var_mu_alpha) + (arma::sum(log_alpha_it) / var_log_alpha)) /
      ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    var_mu_post = 1.0 / ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    mu_alpha_it = R::rnorm(mean_mu_post, sqrt(var_mu_post));

    //Putting new values in the linked list
    update_node_values(head, centers_it, log_alpha_it, sigma_cube_it, L_mat);

    //  Updating death rates
    if(n_cent_it > 0){
      DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
    }

    // Storing MCMC samples
    centers_sample(i, 0) = centers_it;
    log_alpha_sample(i, 0) = log_alpha_it.t();
    sigma_sample(i, 0) = sigma_cube_it;
    // beta_sample(i) = beta_it;
    mu_alpha_sample(i) = mu_alpha_it;
    sample_n_cent(i) = n_cent_it;
    sample_cum_int(i) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  }
  return Rcpp::List::create(Rcpp::Named("log_alpha_sample") = log_alpha_sample,
                            Rcpp::Named("centers_sample") = centers_sample,
                            Rcpp::Named("sigmas_sample") = sigma_sample,
                            Rcpp::Named("mu_alpha_sample") = mu_alpha_sample,
                            Rcpp::Named("n_centers_sample") = sample_n_cent,
                            Rcpp::Named("cumulative_intensity_sample") = sample_cum_int,
                            Rcpp::Named("n_bd_events") = bd_event_vec,
                            Rcpp::Named("bd_vt") = vt_vec);
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Doing a version for continuous points (i.e. not pixelated)
//  Now with fixed noise process (assuming it has already been estimated)
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//' Bayesian SNCP fit using BD-MCMC
//'
//' Run a BD-MCMC chain for a SNCP (mvn dispersion density) on lung ct data using a uniform proposal surface for the BD process
//'
//'
//'
//'
//'
//' @param obs_points matrix of coordinates for observed points
//' @param mean_mu_alpha prior mean for mu_alpha
//' @param sd_log_alpha prior std. dev. for log_alphas
//' @param sd_prop_alpha std. dev. for random walk proposal for log_alphas
//' @param n_it number of iterations to run MCMC
//' @param window_hw half-width of random walk cube for proposing new cluster center locations
//' @param df_iw_prior degrees of freedom for inverse-wishart prior on dispersion matrices
//' @param df_iw_prop degrees of freedom for inverse-wishart proposal on dispersion matrices
//' @param sigma_prior prior for dispersion matrices sigma
//' @param xwin vector that has the min and max x values of the observation window
//' @param ywin vector that has the min and max y values of the observation window
//' @param var_mu_alpha prior variance on mu_alpha
//' @param pen_dist distance for strauss process repulsion
//' @param pen_val vaule between 0 and 1 that is the penalty for 2 cluster centers being within pen_dist of each other
//' @param n_cent_init number of clusters to start with
//' @param prior_n_cent prior value for number of clusters
//' @param max_bd_events max events (births + deaths) to allow at each iteration of BD-MCMC
//' @param max_bd_vt max ammount of virtual time to spend in BD process at each BD_MCMC iteration
//'
//' @author Brian Vestal
//'
//' @return
//' Returns a named list containing posterior samples of model parameters with the following elements:
//' \itemize{
//' \item log_alpha_sample = list whose elements are vectors containing the estimated log-alpha estimates for the clusters present in that iteration
//' \item centers_sample = list whose elements matrices containing the estimated center locations for the clusters present in that iteration
//' \item sigmas_sample = list whose elements are lists containing the estimated covariance matrix estimates for the clusters present in that iteration
//' \item mu_alpha_sample = vector of mu_alpha samples
//' \item n_centers_sample = vector that contains the number of estimated clusters at each iteration,
//' \item cumulative_intensity_sample = vector that contains the cumulative intensity estimate for the SNCP at each iteration
//' }
//'
//' @export
// [[Rcpp::export]]

Rcpp::List sncp_bdmcmc_cont_fixed_noise(arma::mat obs_points,
                                        double mean_mu_alpha,
                                        double sd_log_alpha,
                                        double sd_prop_alpha,
                                        double eta,
                                        int n_it,
                                        double window_hw,
                                        int df_iw_prior,
                                        int df_iw_prop,
                                        arma::mat sigma_prior,
                                        arma::vec xwin,
                                        arma::vec ywin,
                                        double var_mu_alpha,
                                        double pen_dist,
                                        double pen_val,
                                        int n_cent_init,
                                        double prior_n_cent,
                                        int max_bd_events,
                                        double max_bd_vt){
  int i, j, xdim = obs_points.n_cols, n_cent_it, n_points = obs_points.n_rows;
  double ll_cur, ll_prop, mh_log, log_alpha_prop, beta_it = eta, mu_alpha_it, log_pen_val = log(pen_val);
  arma::vec log_alpha_it, mu_alpha_sample(n_it);
  //arma::mat log_alpha_sample(n_it, n_cent);
  arma::umat P_mat;
  arma::uvec P_vec_tmp, bd_idx_inf;
  arma::rowvec center_prop;
  arma::mat centers_it(n_cent_init, xdim), L_mat, L_mat_tmp, centers_tmp, sigma_tmp, L_mat_tmp_BD;
  arma::field<arma::mat> centers_sample(n_it, 1), log_alpha_sample(n_it, 1);
  arma::cube sigma_cube_it(xdim, xdim, n_cent_init);
  arma::field<arma::cube> sigma_sample(n_it, 1);
  double var_log_alpha = pow(sd_log_alpha, 2), mean_mu_post, var_mu_post, LM_slice = (xwin(1) - xwin(0)) * (ywin(1) - ywin(0));
  arma::vec idx_centers_init, DR_vec(n_cent_init), log_alpha_tmp_BD, PP_lik_cur, PP_lik_prop, L_vec_prop, bd_prob_vec;
  double log_bd_surf = log(1.0 / LM_slice), DR_tot, prior_n_cent_log = log(prior_n_cent / LM_slice);
  bool bd_flag;
  double bd_vt, BR_tot = 1;
  int n_bd_events, bd_kill_idx;
  struct Node *head = new Node;
  Node *kill_node;
  arma::vec sample_n_cent(n_it), sample_cum_int(n_it);
  arma::rowvec cent_add(2);
  arma::vec bd_event_vec(n_it), vt_vec(n_it);
  //Rcpp::Rcout << "Initialization of variables OK" << std::endl;

  // xwin(0) = lung_data.col(0).min();
  // xwin(1) = lung_data.col(0).max();
  // ywin(0) = lung_data.col(1).min();
  // ywin(1) = lung_data.col(1).max();
  // double x_range = xwin(1) - xwin(0), y_range = ywin(1) - ywin(0);
  arma::mat centers_init(n_cent_init, xdim);
  for(int jj = 0; jj < n_cent_init; jj++){
    centers_init(jj, 0) = R::runif(xwin(0), xwin(1));
    centers_init(jj, 1) = R::runif(ywin(0), ywin(1));
  }
  // idx_centers_init = RcppArmadillo::sample(lung_data_idx, n_cent_init, false, BD_probs);
  //Rcpp::Rcout << "idx_centers_init = " << idx_centers_init << std::endl;
  //Rcpp::Rcout << "Sample for initial center_idx OK" << std::endl;
  mu_alpha_it = R::rnorm(mean_mu_alpha, sqrt(var_mu_alpha));
  initNode(head, centers_init.row(0), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  for(i = 1; i < n_cent_init; i++){
    addNode(head, centers_init.row(i), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  }
  update_Ps(head, n_cent_init, pen_dist, xdim);
  //log_alpha_it = arma::randn(n_cent_init) * sd_log_alpha + mean_mu_alpha;


  sigma_cube_it = get_sigmas(head, n_cent_init, xdim);
  log_alpha_it = get_log_alphas(head, n_cent_init);
  centers_it = get_centers(head, n_cent_init, xdim);

  //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
  //Rcpp::Rcout << "log_alphas_it = " << log_alpha_it << std::endl;
  //Rcpp::Rcout << "sigma_cube_it = " << sigma_cube_it << std::endl;

  sigma_sample(0, 0) = sigma_cube_it;
  log_alpha_sample(0, 0) = log_alpha_it;
  centers_sample(0, 0) = centers_it;
  mu_alpha_sample(0) = mu_alpha_it;
  sample_n_cent(0) = n_cent_init;
  sample_cum_int(0) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  //Rcpp::Rcout << "Storage of initials OK" << std::endl;
  //L_mat = dmvnrm_vec_arma_cube(obs_points, centers_it, sigma_cube_it);
  //Rcpp::Rcout << "L_mat creation OK" << std::endl;
  PP_lik_cur = calc_PP_lik_vec(head);
  ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
  //Rcpp::Rcout << "ll_cur OK" << std::endl;
  //P_mat = P_mat_gen1(centers_it, pen_dist);
  //Rcpp::Rcout << "P_mat OK" << std::endl;

  //  Initial BD death rates
  DR_vec = calc_DR_vec(head, n_cent_init, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);

  //Rcpp::Rcout << "Initial death rates OK" << std::endl;
  //Rcpp::Rcout << "DR_vec = " << DR_vec << std::endl;

  n_cent_it = n_cent_init;

  //   Birth-death process:
  for(i = 1; i < n_it; ++i){
    bd_flag = true;
    bd_vt = 0;
    n_bd_events = 0;
    while(bd_flag & (n_bd_events < max_bd_events)){
      //Rcpp::Rcout << "n_bd_events = " << n_bd_events << std::endl;
      //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
      //Rcpp::Rcout << "log_alpha_it = " << log_alpha_it << std::endl;
      //display(head);
      DR_tot = arma::sum(DR_vec);

      //Rcpp::Rcout << "DR_tot = " << DR_tot << std::endl;
      //Rcpp::Rcout << "DR_vec = " << DR_tot << std::endl;
      if(DR_vec.has_inf()){
        bd_idx_inf = arma::find_nonfinite(DR_vec);
        //Rcpp::Rcout << "Infinite DR, bd_idx_inf = " << bd_idx_inf << std::endl;
        if(bd_idx_inf.n_elem < 2){
          bd_kill_idx = bd_idx_inf(0);
        }
        else{
          bd_kill_idx = RcppArmadillo::sample(bd_idx_inf, 1, false)(0);
        }

        kill_node = searchNode(head, bd_kill_idx);
        if(n_cent_it < 2){
          PP_lik_cur.zeros();
        }
        else{
          PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
        }
        ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
        deleteNode(&head, kill_node);
        n_bd_events++;
        n_cent_it--;
        //Rcpp::Rcout << "Node killed succesfully (inf)" << std::endl;
      }
      else{
        bd_vt += R::rexp(1.0 / (BR_tot + DR_tot));
        //Rcpp::Rcout << "bd_vt = " << bd_vt << std::endl;
        if(bd_vt < max_bd_vt){
          if(R::runif(0, BR_tot + DR_tot) < BR_tot){
            // idx_birth = RcppArmadillo::sample(lung_data_idx, 1, false)(0);
            cent_add(0) = R::runif(xwin(0), xwin(1));
            cent_add(1) = R::runif(ywin(0), ywin(1));
            insertFront(&head, cent_add, mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
            PP_lik_cur += (head -> alpha) * (head -> mvn_dens);
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            n_bd_events++;
            n_cent_it++;
            //Rcpp::Rcout << "Node birthed succesfully" << std::endl;
          }
          else{
            bd_prob_vec = DR_vec / DR_tot;
            bd_kill_idx = RcppArmadillo::sample(arma::linspace(0, n_cent_it - 1, n_cent_it), 1, false, bd_prob_vec)(0);
            kill_node = searchNode(head, bd_kill_idx);
            if(n_cent_it < 2){
              PP_lik_cur.zeros();
            }
            else{
              PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
            }
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            deleteNode(&head, kill_node);
            n_bd_events++;
            n_cent_it--;
            //Rcpp::Rcout << "Node killed succesfully" << std::endl;
          }
        }
        else{
          bd_flag = false;
        }
      }
      if(bd_flag){
        if(n_cent_it < 1){
          DR_vec.set_size(1);
          DR_vec(0) = 0;
        }
        else{
          update_Ps(head, n_cent_it, pen_dist, xdim);
          //Rcpp::Rcout << "Ps updated" << std::endl;
          DR_vec.set_size(n_cent_it);
          DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
          //Rcpp::Rcout << "DRs updated succesfully" << std::endl;
        }
      }
    }
    bd_event_vec(i) = n_bd_events;
    vt_vec(i) = bd_vt;
    //Rcpp::Rcout << "Left BD Process" << std::endl;

    L_mat.set_size(n_points, n_cent_it);
    L_mat = get_L_mat(head, n_cent_it, n_points);

    log_alpha_it.set_size(n_cent_it);
    log_alpha_it = get_log_alphas(head, n_cent_it);

    centers_it.set_size(n_cent_it, xdim);
    centers_it = get_centers(head, n_cent_it, xdim);
    P_mat = P_mat_gen1(centers_it, pen_dist);

    sigma_cube_it.set_size(xdim, xdim, n_cent_it);
    sigma_cube_it = get_sigmas(head, n_cent_it, xdim);

    //PP_lik_cur = calc_PP_lik_vec(head);
    //ll_cur = arma::sum(arma::log(PP_lik_cur + beta_it));
    for(j = 0; j < n_cent_it; j++){
      //log_alpha_temp = log_alpha_it;
      log_alpha_prop = R::rnorm(log_alpha_it(j), sd_prop_alpha);
      //log_alpha_temp(j) = log_alpha_prop;
      PP_lik_prop = PP_lik_cur + (exp(log_alpha_prop) - exp(log_alpha_it(j))) * L_mat.col(j);
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop - ll_cur -
        exp(log_alpha_prop) + exp(log_alpha_it(j)) +
        (1.0/(2.0 * pow(sd_log_alpha, 2))) * (pow(log_alpha_it(j) - mean_mu_alpha, 2) - pow(log_alpha_prop - mean_mu_alpha, 2));
      if(R::runif(0, 1) < exp(mh_log)){
        //Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        log_alpha_it(j) = log_alpha_prop;
      }
    }

    // Doing center location MH steps
    for(j = 0; j < n_cent_it; j++){
      center_prop = centers_it.row(j);
      center_prop(0) += R::runif(-window_hw, window_hw);
      center_prop(1) += R::runif(-window_hw, window_hw);
      if(center_prop(0) < xwin(1) & center_prop(0) > xwin(0) &
         center_prop(1) < ywin(1) & center_prop(1) > ywin(0)){
        //L_mat_tmp = L_mat;
        //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        L_vec_prop = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        P_vec_tmp = P_mat_gen2(centers_it, center_prop, pen_dist);
        //P_vec_tmp.shed_row(j);
        PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
        ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
        mh_log = ll_prop - ll_cur + (arma::sum(P_vec_tmp) - P_vec_tmp(j) - (arma::sum(P_mat.col(j)) - 1)) * log_pen_val;
        if(R::runif(0, 1) < exp(mh_log)){
          PP_lik_cur = PP_lik_prop;
          ll_cur = ll_prop;
          centers_it.row(j) = center_prop;
          L_mat.col(j) = L_vec_prop;
          P_mat = P_mat_gen1(centers_it, pen_dist);
        }
      }
    }

    // Doing sigma matrix MH steps
    for(j = 0; j < n_cent_it; j++){
      sigma_tmp = riwish_arma(df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j));
      /*
       log(diwish(W = sigma_mat_it[[j]], v = df_inv_wish_prior,
       S = (df_inv_wish_prior - 3) * S_mat_prior)) + ## Prior
       log(diwish(W = sigma_prop_j, v = df_inv_wish_prop,
       S = (df_inv_wish_prop - 3) * sigma_mat_it[[j]]))
       */
      //L_mat_tmp = L_mat;
      //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      L_vec_prop = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop +
        log(diwish_arma(sigma_tmp, df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) +
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prop, (df_iw_prop - 3.0) * sigma_tmp)) -
        ll_cur -
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) -
        log(diwish_arma(sigma_tmp, df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j)));
      if(R::runif(0, 1) < exp(mh_log)){
        ////Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        L_mat.col(j) = L_vec_prop;
        sigma_cube_it.slice(j) = sigma_tmp;
        //L_mat = L_mat_tmp;
      }
    }

    //  Updating beta
    // beta_tmp =  std::fabs(beta_it + R::runif(-.0005, .0005));
    // ll_prop = arma::sum((arma::log(PP_lik_cur + beta_tmp)));
    // mh_log = ll_prop - beta_tmp * LM_slice + R::dgamma(beta_tmp, 0.01, 1.0 / 0.01, 1) -
    //   (ll_cur - beta_it * LM_slice + R::dgamma(beta_it, 0.01, 1.0 / 0.01, 1));
    // if(R::runif(0, 1) < exp(mh_log)){
    //   ////Rcpp::Rcout <<"MH accept" << std::endl;
    //   ll_cur = ll_prop;
    //   beta_it = beta_tmp;
    // }

    //  Updating mu_alpha

    mean_mu_post = ((mean_mu_alpha / var_mu_alpha) + (arma::sum(log_alpha_it) / var_log_alpha)) /
      ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    var_mu_post = 1.0 / ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    mu_alpha_it = R::rnorm(mean_mu_post, sqrt(var_mu_post));

    //Putting new values in the linked list
    update_node_values(head, centers_it, log_alpha_it, sigma_cube_it, L_mat);

    //  Updating death rates
    if(n_cent_it > 0){
      DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
    }

    // Storing MCMC samples
    centers_sample(i, 0) = centers_it;
    log_alpha_sample(i, 0) = log_alpha_it.t();
    sigma_sample(i, 0) = sigma_cube_it;
    // beta_sample(i) = beta_it;
    mu_alpha_sample(i) = mu_alpha_it;
    sample_n_cent(i) = n_cent_it;
    sample_cum_int(i) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  }
  return Rcpp::List::create(Rcpp::Named("log_alpha_sample") = log_alpha_sample,
                            Rcpp::Named("centers_sample") = centers_sample,
                            Rcpp::Named("sigmas_sample") = sigma_sample,
                            Rcpp::Named("mu_alpha_sample") = mu_alpha_sample,
                            Rcpp::Named("n_centers_sample") = sample_n_cent,
                            Rcpp::Named("cumulative_intensity_sample") = sample_cum_int,
                            Rcpp::Named("n_bd_events") = bd_event_vec,
                            Rcpp::Named("bd_vt") = vt_vec);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Doing a version with novel BD proposal surface
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//' Bayesian SNCP fit using BD-MCMC with informative proposal surface
//'
//' Run a BD-MCMC chain for a SNCP (mvn dispersion density) on lung ct data using the kernal-smoothed point pattern as the proposal surface in the BD process
//'
//'
//'
//'
//'
//' @param obs_points matrix of coordinates for observed points
//' @param mean_mu_alpha prior mean for mu_alpha
//' @param sd_log_alpha prior std. dev. for log_alphas
//' @param sd_prop_alpha std. dev. for random walk proposal for log_alphas
//' @param beta prior estimate of homogeneous Poisson intensity function for cluster centers
//' @param n_it number of iterations to run MCMC
//' @param window_hw half-width of random walk cube for proposing new cluster center locations
//' @param df_iw_prior degrees of freedom for inverse-wishart prior on dispersion matrices
//' @param df_iw_prop degrees of freedom for inverse-wishart proposal on dispersion matrices
//' @param sigma_prior prior for dispersion matrices sigma
//' @param lung_data matrix that contains coordinates of all points within 2d slice of lung CT that is lung tissue
//' @param var_mu_alpha prior variance on mu_alpha
//' @param pen_dist distance for strauss process repulsion
//' @param pen_val vaule between 0 and 1 that is the penalty for 2 cluster centers being within pen_dist of each other
//' @param n_cent_init number of clusters to start with
//' @param prior_n_cent prior value for number of clusters
//' @param max_bd_events max events (births + deaths) to allow at each iteration of BD-MCMC
//' @param max_bd_vt max ammount of virtual time to spend in BD process at each BD_MCMC iteration
//' @param sigma_smooth smoothing bandwidth for BD proposal surface
//'
//' @author Brian Vestal
//'
//' @return
//' Returns a named list containing posterior samples of model parameters with the following elements:
//' \itemize{
//' \item log_alpha_sample = list whose elements are vectors containing the estimated log-alpha estimates for the clusters present in that iteration
//' \item centers_sample = list whose elements matrices containing the estimated center locations for the clusters present in that iteration
//' \item sigmas_sample = list whose elements are lists containing the estimated covariance matrix estimates for the clusters present in that iteration
//' \item mu_alpha_sample = vector of mu_alpha samples
//' \item eta_sample = vector of eta samples,
//' \item n_centers_sample = vector that contains the number of estimated clusters at each iteration,
//' \item cumulative_intensity_sample = vector that contains the cumulative intensity estimate for the SNCP at each iteration
//' }
//'
//' @export
// [[Rcpp::export]]

Rcpp::List sncp_bdmcmc_smooth(arma::mat obs_points,
                              double mean_mu_alpha,
                              double sd_log_alpha,
                              double sd_prop_alpha,
                              double beta,
                              int n_it,
                              double window_hw,
                              int df_iw_prior,
                              int df_iw_prop,
                              arma::mat sigma_prior,
                              arma::mat lung_data,
                              double var_mu_alpha,
                              double pen_dist,
                              double pen_val,
                              int n_cent_init,
                              double prior_n_cent,
                              int max_bd_events,
                              double max_bd_vt,
                              double sigma_smooth = 100){
  int i, j, xdim = obs_points.n_cols, n_cent_it, idx_birth, n_points = obs_points.n_rows;
  double ll_cur, ll_prop, mh_log, log_alpha_prop, beta_it = beta, mu_alpha_it, log_pen_val = log(pen_val);
  arma::vec log_alpha_it, beta_sample(n_it), mu_alpha_sample(n_it), xwin(2), ywin(2);
  //arma::mat log_alpha_sample(n_it, n_cent);
  arma::umat P_mat;
  arma::uvec P_vec_tmp, bd_idx_inf;
  arma::rowvec center_prop;
  arma::mat centers_it(n_cent_init, xdim), L_mat, L_mat_tmp, centers_tmp, sigma_tmp, L_mat_tmp_BD;
  arma::field<arma::mat> centers_sample(n_it, 1), log_alpha_sample(n_it, 1);
  arma::cube sigma_cube_it(xdim, xdim, n_cent_init);
  arma::field<arma::cube> sigma_sample(n_it, 1);
  double beta_tmp, var_log_alpha = pow(sd_log_alpha, 2), mean_mu_post, var_mu_post, LM_slice = lung_data.n_rows;
  arma::vec lung_data_idx = arma::linspace(0, LM_slice - 1, LM_slice), BD_probs;
  arma::vec idx_centers_init, DR_vec(n_cent_init), log_alpha_tmp_BD, PP_lik_cur, PP_lik_prop, L_vec_prop, bd_prob_vec;
  double DR_tot, prior_n_cent_log = log(prior_n_cent / LM_slice);
  bool bd_flag;
  double bd_vt, BR_tot = 1;
  int n_bd_events, bd_kill_idx;
  struct Node *head = new Node;
  Node *kill_node;
  arma::vec sample_n_cent(n_it), sample_cum_int(n_it);
  //Rcpp::Rcout << "Initialization of variables OK" << std::endl;

  xwin(0) = lung_data.col(0).min();
  xwin(1) = lung_data.col(0).max();
  ywin(0) = lung_data.col(1).min();
  ywin(1) = lung_data.col(1).max();

  BD_probs = fastDMVNorm_diag_norm_sum(lung_data, obs_points, sigma_smooth);
  BD_probs = BD_probs / arma::sum(BD_probs);

  idx_centers_init = RcppArmadillo::sample(lung_data_idx, n_cent_init, false, BD_probs);
  //Rcpp::Rcout << "idx_centers_init = " << idx_centers_init << std::endl;
  //Rcpp::Rcout << "Sample for initial center_idx OK" << std::endl;
  mu_alpha_it = R::rnorm(mean_mu_alpha, sqrt(var_mu_alpha));
  initNode(head, lung_data.row(idx_centers_init(0)), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log(BD_probs(idx_centers_init(0))));
  for(i = 1; i < n_cent_init; i++){
    addNode(head, lung_data.row(idx_centers_init(i)), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log(BD_probs(idx_centers_init(i))));
  }
  update_Ps(head, n_cent_init, pen_dist, xdim);
  //log_alpha_it = arma::randn(n_cent_init) * sd_log_alpha + mean_mu_alpha;


  sigma_cube_it = get_sigmas(head, n_cent_init, xdim);
  log_alpha_it = get_log_alphas(head, n_cent_init);
  centers_it = get_centers(head, n_cent_init, xdim);

  //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
  //Rcpp::Rcout << "log_alphas_it = " << log_alpha_it << std::endl;
  //Rcpp::Rcout << "sigma_cube_it = " << sigma_cube_it << std::endl;

  sigma_sample(0, 0) = sigma_cube_it;
  log_alpha_sample(0, 0) = log_alpha_it;
  centers_sample(0, 0) = centers_it;
  mu_alpha_sample(0) = mu_alpha_it;
  beta_sample(0) = beta_it;
  sample_n_cent(0) = n_cent_init;
  sample_cum_int(0) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  //Rcpp::Rcout << "Storage of initials OK" << std::endl;
  //L_mat = dmvnrm_vec_arma_cube(obs_points, centers_it, sigma_cube_it);
  //Rcpp::Rcout << "L_mat creation OK" << std::endl;
  PP_lik_cur = calc_PP_lik_vec(head);
  ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
  //Rcpp::Rcout << "ll_cur OK" << std::endl;
  //P_mat = P_mat_gen1(centers_it, pen_dist);
  //Rcpp::Rcout << "P_mat OK" << std::endl;

  //  Initial BD death rates
  DR_vec = calc_DR_vec(head, n_cent_init, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);

  //Rcpp::Rcout << "Initial death rates OK" << std::endl;
  //Rcpp::Rcout << "DR_vec = " << DR_vec << std::endl;

  n_cent_it = n_cent_init;

  //   Birth-death process:
  for(i = 1; i < n_it; ++i){
    bd_flag = true;
    bd_vt = 0;
    n_bd_events = 0;
    while(bd_flag & (n_bd_events < max_bd_events)){
      //Rcpp::Rcout << "n_bd_events = " << n_bd_events << std::endl;
      //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
      //Rcpp::Rcout << "log_alpha_it = " << log_alpha_it << std::endl;
      //display(head);
      DR_tot = arma::sum(DR_vec);

      //Rcpp::Rcout << "DR_tot = " << DR_tot << std::endl;
      //Rcpp::Rcout << "DR_vec = " << DR_tot << std::endl;
      if(DR_vec.has_inf()){
        bd_idx_inf = arma::find_nonfinite(DR_vec);
        //Rcpp::Rcout << "Infinite DR, bd_idx_inf = " << bd_idx_inf << std::endl;
        if(bd_idx_inf.n_elem < 2){
          bd_kill_idx = bd_idx_inf(0);
        }
        else{
          bd_kill_idx = RcppArmadillo::sample(bd_idx_inf, 1, false)(0);
        }

        kill_node = searchNode(head, bd_kill_idx);
        if(n_cent_it < 2){
          PP_lik_cur.zeros();
        }
        else{
          PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
        }
        ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
        deleteNode(&head, kill_node);
        n_bd_events++;
        n_cent_it--;
        //Rcpp::Rcout << "Node killed succesfully (inf)" << std::endl;
      }
      else{
        bd_vt += R::rexp(1.0 / (BR_tot + DR_tot));
        //Rcpp::Rcout << "bd_vt = " << bd_vt << std::endl;
        if(bd_vt < max_bd_vt){
          if(R::runif(0, BR_tot + DR_tot) < BR_tot){
            idx_birth = RcppArmadillo::sample(lung_data_idx, 1, false)(0);
            insertFront(&head, lung_data.row(idx_birth), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log(BD_probs(idx_birth)));
            PP_lik_cur += (head -> alpha) * (head -> mvn_dens);
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            n_bd_events++;
            n_cent_it++;
            //Rcpp::Rcout << "Node birthed succesfully" << std::endl;
          }
          else{
            bd_prob_vec = DR_vec / DR_tot;
            bd_kill_idx = RcppArmadillo::sample(arma::linspace(0, n_cent_it - 1, n_cent_it), 1, false, bd_prob_vec)(0);
            kill_node = searchNode(head, bd_kill_idx);
            if(n_cent_it < 2){
              PP_lik_cur.zeros();
            }
            else{
              PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
            }
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            deleteNode(&head, kill_node);
            n_bd_events++;
            n_cent_it--;
            //Rcpp::Rcout << "Node killed succesfully" << std::endl;
          }
        }
        else{
          bd_flag = false;
        }
      }
      if(bd_flag){
        if(n_cent_it < 1){
          DR_vec.set_size(1);
          DR_vec(0) = 0;
        }
        else{
          update_Ps(head, n_cent_it, pen_dist, xdim);
          //Rcpp::Rcout << "Ps updated" << std::endl;
          DR_vec.set_size(n_cent_it);
          DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
          //Rcpp::Rcout << "DRs updated succesfully" << std::endl;
        }
      }
    }

    //Rcpp::Rcout << "Left BD Process" << std::endl;

    L_mat.set_size(n_points, n_cent_it);
    L_mat = get_L_mat(head, n_cent_it, n_points);

    log_alpha_it.set_size(n_cent_it);
    log_alpha_it = get_log_alphas(head, n_cent_it);

    centers_it.set_size(n_cent_it, xdim);
    centers_it = get_centers(head, n_cent_it, xdim);
    P_mat = P_mat_gen1(centers_it, pen_dist);

    sigma_cube_it.set_size(xdim, xdim, n_cent_it);
    sigma_cube_it = get_sigmas(head, n_cent_it, xdim);

    //PP_lik_cur = calc_PP_lik_vec(head);
    //ll_cur = arma::sum(arma::log(PP_lik_cur + beta_it));
    for(j = 0; j < n_cent_it; j++){
      //log_alpha_temp = log_alpha_it;
      log_alpha_prop = R::rnorm(log_alpha_it(j), sd_prop_alpha);
      //log_alpha_temp(j) = log_alpha_prop;
      PP_lik_prop = PP_lik_cur + (exp(log_alpha_prop) - exp(log_alpha_it(j))) * L_mat.col(j);
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop - ll_cur -
        exp(log_alpha_prop) + exp(log_alpha_it(j)) +
        (1.0/(2.0 * pow(sd_log_alpha, 2))) * (pow(log_alpha_it(j) - mean_mu_alpha, 2) - pow(log_alpha_prop - mean_mu_alpha, 2));
      if(R::runif(0, 1) < exp(mh_log)){
        //Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        log_alpha_it(j) = log_alpha_prop;
      }
    }

    // Doing center location MH steps
    for(j = 0; j < n_cent_it; j++){
      center_prop = centers_it.row(j);
      center_prop(0) += R::runif(-window_hw, window_hw);
      center_prop(1) += R::runif(-window_hw, window_hw);
      if(center_prop(0) < xwin(1) & center_prop(0) > xwin(0) &
         center_prop(1) < ywin(1) & center_prop(1) > ywin(0)){
        //L_mat_tmp = L_mat;
        //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        L_vec_prop = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        P_vec_tmp = P_mat_gen2(centers_it, center_prop, pen_dist);
        //P_vec_tmp.shed_row(j);
        PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
        ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
        mh_log = ll_prop - ll_cur + (arma::sum(P_vec_tmp) - P_vec_tmp(j) - (arma::sum(P_mat.col(j)) - 1)) * log_pen_val;
        if(R::runif(0, 1) < exp(mh_log)){
          PP_lik_cur = PP_lik_prop;
          ll_cur = ll_prop;
          centers_it.row(j) = center_prop;
          L_mat.col(j) = L_vec_prop;
          P_mat = P_mat_gen1(centers_it, pen_dist);
        }
      }
    }

    // Doing sigma matrix MH steps
    for(j = 0; j < n_cent_it; j++){
      sigma_tmp = riwish_arma(df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j));
      /*
       log(diwish(W = sigma_mat_it[[j]], v = df_inv_wish_prior,
       S = (df_inv_wish_prior - 3) * S_mat_prior)) + ## Prior
       log(diwish(W = sigma_prop_j, v = df_inv_wish_prop,
       S = (df_inv_wish_prop - 3) * sigma_mat_it[[j]]))
       */
      //L_mat_tmp = L_mat;
      //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      L_vec_prop = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop +
        log(diwish_arma(sigma_tmp, df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) +
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prop, (df_iw_prop - 3.0) * sigma_tmp)) -
        ll_cur -
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) -
        log(diwish_arma(sigma_tmp, df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j)));
      if(R::runif(0, 1) < exp(mh_log)){
        ////Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        L_mat.col(j) = L_vec_prop;
        sigma_cube_it.slice(j) = sigma_tmp;
        //L_mat = L_mat_tmp;
      }
    }

    //  Updating beta
    beta_tmp =  std::fabs(beta_it + R::runif(-.0005, .0005));
    ll_prop = arma::sum((arma::log(PP_lik_cur + beta_tmp)));
    mh_log = ll_prop - beta_tmp * LM_slice + R::dgamma(beta_tmp, 0.01, 1.0 / 0.01, 1) -
      (ll_cur - beta_it * LM_slice + R::dgamma(beta_it, 0.01, 1.0 / 0.01, 1));
    if(R::runif(0, 1) < exp(mh_log)){
      ////Rcpp::Rcout <<"MH accept" << std::endl;
      ll_cur = ll_prop;
      beta_it = beta_tmp;
    }

    //  Updating mu_alpha

    mean_mu_post = ((mean_mu_alpha / var_mu_alpha) + (arma::sum(log_alpha_it) / var_log_alpha)) /
      ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    var_mu_post = 1.0 / ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    mu_alpha_it = R::rnorm(mean_mu_post, sqrt(var_mu_post));

    //Putting new values in the linked list
    update_node_values(head, centers_it, log_alpha_it, sigma_cube_it, L_mat);

    //  Updating death rates
    if(n_cent_it > 0){
      DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
    }

    // Storing MCMC samples
    centers_sample(i, 0) = centers_it;
    log_alpha_sample(i, 0) = log_alpha_it.t();
    sigma_sample(i, 0) = sigma_cube_it;
    beta_sample(i) = beta_it;
    mu_alpha_sample(i) = mu_alpha_it;
    sample_n_cent(i) = n_cent_it;
    sample_cum_int(i) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  }
  return Rcpp::List::create(Rcpp::Named("log_alpha_sample") = log_alpha_sample,
                            Rcpp::Named("centers_sample") = centers_sample,
                            Rcpp::Named("sigmas_sample") = sigma_sample,
                            Rcpp::Named("mu_alpha_sample") = mu_alpha_sample,
                            Rcpp::Named("beta_sample") = beta_sample,
                            Rcpp::Named("n_centers_sample") = sample_n_cent,
                            Rcpp::Named("cumulative_intensity_sample") = sample_cum_int);
}


///////////////////////////////////////////////////////////////////////////////////////
/////     Doing a version with a uniform dispersion density instead      //////////////
///////////////////////////////////////////////////////////////////////////////////////

//' Bayesian SNCP fit using BD-MCMC
//'
//' Run a BD-MCMC chain for a SNCP (uniform dispersion density) on lung ct data with uniform proposal surface for BD process
//'
//'
//'
//'
//'
//' @param obs_points numeric matrix of coordinates for observed points
//' @param mean_mu_alpha prior mean for mu_alpha
//' @param sd_log_alpha prior std. dev. for log_alphas
//' @param sd_prop_alpha std. dev. for random walk proposal for log_alphas
//' @param beta prior estimate of homogeneous Poisson intensity function for cluster centers
//' @param n_it number of iterations to run MCMC
//' @param window_hw half-width of random walk cube for proposing new cluster center locations
//' @param df_iw_prior degrees of freedom for inverse-wishart prior on dispersion matrices
//' @param df_iw_prop degrees of freedom for inverse-wishart proposal on dispersion matrices
//' @param sigma_prior prior for dispersion matrices sigma
//' @param lung_data matrix that contains coordinates of all points within 2d slice of lung CT that is lung tissue
//' @param var_mu_alpha prior variance on mu_alpha
//' @param pen_dist distance for strauss process repulsion
//' @param pen_val vaule between 0 and 1 that is the penalty for 2 cluster centers being within pen_dist of each other
//' @param n_cent_init number of clusters to start with
//' @param prior_n_cent prior value for number of clusters
//' @param max_bd_events max events (births + deaths) to allow at each iteration of BD-MCMC
//' @param max_bd_vt max ammount of virtual time to spend in BD process at each BD_MCMC iteration
//'
//' @author Brian Vestal
//'
//' @return
//' Returns a named list containing posterior samples of model parameters with the following elements:
//' \itemize{
//' \item log_alpha_sample = list whose elements are vectors containing the estimated log-alpha estimates for the clusters present in that iteration
//' \item centers_sample = list whose elements matrices containing the estimated center locations for the clusters present in that iteration
//' \item sigmas_sample = list whose elements are lists containing the estimated covariance matrix estimates for the clusters present in that iteration
//' \item mu_alpha_sample = vector of mu_alpha samples
//' \item eta_sample = vector of eta samples,
//' \item n_centers_sample = vector that contains the number of estimated clusters at each iteration,
//' \item cumulative_intensity_sample = vector that contains the cumulative intensity estimate for the SNCP at each iteration
//' }
//'
//' @export
// [[Rcpp::export]]

Rcpp::List sncp_bdmcmc_unif(arma::mat obs_points,
                            double mean_mu_alpha,
                            double sd_log_alpha,
                            double sd_prop_alpha,
                            double beta,
                            int n_it,
                            double window_hw,
                            int df_iw_prior,
                            int df_iw_prop,
                            arma::mat sigma_prior,
                            arma::mat lung_data,
                            double var_mu_alpha,
                            double pen_dist,
                            double pen_val,
                            int n_cent_init,
                            double prior_n_cent,
                            int max_bd_events,
                            double max_bd_vt){
  int i, j, xdim = obs_points.n_cols, n_cent_it, idx_birth, n_points = obs_points.n_rows;
  double ll_cur, ll_prop, mh_log, log_alpha_prop, beta_it = beta, mu_alpha_it, log_pen_val = log(pen_val);
  arma::vec log_alpha_it, beta_sample(n_it), mu_alpha_sample(n_it), xwin(2), ywin(2);
  //arma::mat log_alpha_sample(n_it, n_cent);
  arma::umat P_mat;
  arma::uvec P_vec_tmp, bd_idx_inf;
  arma::rowvec center_prop;
  arma::mat centers_it(n_cent_init, xdim), L_mat, L_mat_tmp, centers_tmp, sigma_tmp, L_mat_tmp_BD;
  arma::field<arma::mat> centers_sample(n_it, 1), log_alpha_sample(n_it, 1);
  arma::cube sigma_cube_it(xdim, xdim, n_cent_init);
  arma::field<arma::cube> sigma_sample(n_it, 1);
  double beta_tmp, var_log_alpha = pow(sd_log_alpha, 2), mean_mu_post, var_mu_post, LM_slice = lung_data.n_rows;
  arma::vec lung_data_idx = arma::linspace(0, LM_slice - 1, LM_slice), BD_probs = arma::ones(LM_slice) / LM_slice;
  arma::vec idx_centers_init, DR_vec(n_cent_init), log_alpha_tmp_BD, PP_lik_cur, PP_lik_prop, L_vec_prop, bd_prob_vec;
  double log_bd_surf = log(BD_probs(0)), DR_tot, prior_n_cent_log = log(prior_n_cent / LM_slice);
  bool bd_flag;
  double bd_vt, BR_tot = 1;
  int n_bd_events, bd_kill_idx;
  struct Node *head = new Node;
  Node *kill_node;
  arma::vec sample_n_cent(n_it), sample_cum_int(n_it);
  //Rcpp::Rcout << "Initialization of variables OK" << std::endl;

  xwin(0) = lung_data.col(0).min();
  xwin(1) = lung_data.col(0).max();
  ywin(0) = lung_data.col(1).min();
  ywin(1) = lung_data.col(1).max();

  idx_centers_init = RcppArmadillo::sample(lung_data_idx, n_cent_init, false, BD_probs);
  //Rcpp::Rcout << "idx_centers_init = " << idx_centers_init << std::endl;
  //Rcpp::Rcout << "Sample for initial center_idx OK" << std::endl;
  mu_alpha_it = R::rnorm(mean_mu_alpha, sqrt(var_mu_alpha));
  initNode_u(head, lung_data.row(idx_centers_init(0)), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points);
  for(i = 1; i < n_cent_init; i++){
    addNode_u(head, lung_data.row(idx_centers_init(i)), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points);
  }
  update_Ps(head, n_cent_init, pen_dist, xdim);
  //log_alpha_it = arma::randn(n_cent_init) * sd_log_alpha + mean_mu_alpha;


  sigma_cube_it = get_sigmas(head, n_cent_init, xdim);
  log_alpha_it = get_log_alphas(head, n_cent_init);
  centers_it = get_centers(head, n_cent_init, xdim);

  //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
  //Rcpp::Rcout << "log_alphas_it = " << log_alpha_it << std::endl;
  //Rcpp::Rcout << "sigma_cube_it = " << sigma_cube_it << std::endl;

  sigma_sample(0, 0) = sigma_cube_it;
  log_alpha_sample(0, 0) = log_alpha_it;
  centers_sample(0, 0) = centers_it;
  mu_alpha_sample(0) = mu_alpha_it;
  beta_sample(0) = beta_it;
  sample_n_cent(0) = n_cent_init;
  sample_cum_int(0) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  //Rcpp::Rcout << "Storage of initials OK" << std::endl;
  //L_mat = dmvnrm_vec_arma_cube(obs_points, centers_it, sigma_cube_it);
  //Rcpp::Rcout << "L_mat creation OK" << std::endl;
  PP_lik_cur = calc_PP_lik_vec(head);
  ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
  //Rcpp::Rcout << "ll_cur OK" << std::endl;
  //P_mat = P_mat_gen1(centers_it, pen_dist);
  //Rcpp::Rcout << "P_mat OK" << std::endl;

  //  Initial BD death rates
  DR_vec = calc_DR_vec2(head, n_cent_init, PP_lik_cur, ll_cur, beta_it, n_points, log_bd_surf, log_pen_val, prior_n_cent_log);

  //Rcpp::Rcout << "Initial death rates OK" << std::endl;
  //Rcpp::Rcout << "DR_vec = " << DR_vec << std::endl;

  n_cent_it = n_cent_init;

  //   Birth-death process:
  for(i = 1; i < n_it; ++i){
    bd_flag = true;
    bd_vt = 0;
    n_bd_events = 0;
    while(bd_flag & (n_bd_events < max_bd_events)){
      //Rcpp::Rcout << "n_bd_events = " << n_bd_events << std::endl;
      //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
      //Rcpp::Rcout << "log_alpha_it = " << log_alpha_it << std::endl;
      //display(head);
      DR_tot = arma::sum(DR_vec);

      //Rcpp::Rcout << "DR_tot = " << DR_tot << std::endl;
      //Rcpp::Rcout << "DR_vec = " << DR_tot << std::endl;
      if(DR_vec.has_inf()){
        bd_idx_inf = arma::find_nonfinite(DR_vec);
        //Rcpp::Rcout << "Infinite DR, bd_idx_inf = " << bd_idx_inf << std::endl;
        if(bd_idx_inf.n_elem < 2){
          bd_kill_idx = bd_idx_inf(0);
        }
        else{
          bd_kill_idx = RcppArmadillo::sample(bd_idx_inf, 1, false)(0);
        }

        kill_node = searchNode(head, bd_kill_idx);
        if(n_cent_it < 2){
          PP_lik_cur.zeros();
        }
        else{
          PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
        }
        ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
        deleteNode(&head, kill_node);
        n_bd_events++;
        n_cent_it--;
        //Rcpp::Rcout << "Node killed succesfully (inf)" << std::endl;
      }
      else{
        bd_vt += R::rexp(1.0 / (BR_tot + DR_tot));
        //Rcpp::Rcout << "bd_vt = " << bd_vt << std::endl;
        if(bd_vt < max_bd_vt){
          if(R::runif(0, BR_tot + DR_tot) < BR_tot){
            idx_birth = RcppArmadillo::sample(lung_data_idx, 1, false)(0);
            insertFront_u(&head, lung_data.row(idx_birth), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points);
            PP_lik_cur += (head -> alpha) * (head -> mvn_dens);
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            n_bd_events++;
            n_cent_it++;
            //Rcpp::Rcout << "Node birthed succesfully" << std::endl;
          }
          else{
            bd_prob_vec = DR_vec / DR_tot;
            bd_kill_idx = RcppArmadillo::sample(arma::linspace(0, n_cent_it - 1, n_cent_it), 1, false, bd_prob_vec)(0);
            kill_node = searchNode(head, bd_kill_idx);
            if(n_cent_it < 2){
              PP_lik_cur.zeros();
            }
            else{
              PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
            }
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            deleteNode(&head, kill_node);
            n_bd_events++;
            n_cent_it--;
            //Rcpp::Rcout << "Node killed succesfully" << std::endl;
          }
        }
        else{
          bd_flag = false;
        }
      }
      if(bd_flag){
        if(n_cent_it < 1){
          DR_vec.set_size(1);
          DR_vec(0) = 0;
        }
        else{
          update_Ps(head, n_cent_it, pen_dist, xdim);
          //Rcpp::Rcout << "Ps updated" << std::endl;
          DR_vec.set_size(n_cent_it);
          DR_vec = calc_DR_vec2(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_bd_surf, log_pen_val, prior_n_cent_log);
          //Rcpp::Rcout << "DRs updated succesfully" << std::endl;
        }
      }
    }

    //Rcpp::Rcout << "Left BD Process" << std::endl;

    L_mat.set_size(n_points, n_cent_it);
    L_mat = get_L_mat(head, n_cent_it, n_points);

    log_alpha_it.set_size(n_cent_it);
    log_alpha_it = get_log_alphas(head, n_cent_it);

    centers_it.set_size(n_cent_it, xdim);
    centers_it = get_centers(head, n_cent_it, xdim);
    P_mat = P_mat_gen1(centers_it, pen_dist);

    sigma_cube_it.set_size(xdim, xdim, n_cent_it);
    sigma_cube_it = get_sigmas(head, n_cent_it, xdim);

    //PP_lik_cur = calc_PP_lik_vec(head);
    //ll_cur = arma::sum(arma::log(PP_lik_cur + beta_it));
    for(j = 0; j < n_cent_it; j++){
      //log_alpha_temp = log_alpha_it;
      log_alpha_prop = R::rnorm(log_alpha_it(j), sd_prop_alpha);
      //log_alpha_temp(j) = log_alpha_prop;
      PP_lik_prop = PP_lik_cur + (exp(log_alpha_prop) - exp(log_alpha_it(j))) * L_mat.col(j);
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop - ll_cur -
        exp(log_alpha_prop) + exp(log_alpha_it(j)) +
        (1.0/(2.0 * pow(sd_log_alpha, 2))) * (pow(log_alpha_it(j) - mean_mu_alpha, 2) - pow(log_alpha_prop - mean_mu_alpha, 2));
      if(R::runif(0, 1) < exp(mh_log)){
        //Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        log_alpha_it(j) = log_alpha_prop;
      }
    }

    // Doing center location MH steps
    for(j = 0; j < n_cent_it; j++){
      center_prop = centers_it.row(j);
      center_prop(0) += R::runif(-window_hw, window_hw);
      center_prop(1) += R::runif(-window_hw, window_hw);
      if(center_prop(0) < xwin(1) & center_prop(0) > xwin(0) &
         center_prop(1) < ywin(1) & center_prop(1) > ywin(0)){
        //L_mat_tmp = L_mat;
        //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        //L_vec_prop = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        L_vec_prop = mahaInt_bin(obs_points, center_prop.t(), sigma_cube_it.slice(j), 1, false);
        P_vec_tmp = P_mat_gen2(centers_it, center_prop, pen_dist);
        //P_vec_tmp.shed_row(j);
        PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
        ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
        mh_log = ll_prop - ll_cur + (arma::sum(P_vec_tmp) - P_vec_tmp(j) - (arma::sum(P_mat.col(j)) - 1)) * log_pen_val;
        if(R::runif(0, 1) < exp(mh_log)){
          PP_lik_cur = PP_lik_prop;
          ll_cur = ll_prop;
          centers_it.row(j) = center_prop;
          L_mat.col(j) = L_vec_prop;
          P_mat = P_mat_gen1(centers_it, pen_dist);
        }
      }
    }

    // Doing sigma matrix MH steps
    for(j = 0; j < n_cent_it; j++){
      sigma_tmp = riwish_arma(df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j));
      /*
       log(diwish(W = sigma_mat_it[[j]], v = df_inv_wish_prior,
       S = (df_inv_wish_prior - 3) * S_mat_prior)) + ## Prior
       log(diwish(W = sigma_prop_j, v = df_inv_wish_prop,
       S = (df_inv_wish_prop - 3) * sigma_mat_it[[j]]))
       */
      //L_mat_tmp = L_mat;
      //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      //L_vec_prop = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      L_vec_prop = mahaInt_bin(obs_points, centers_it.row(j).t(), sigma_tmp, 1, false);
      PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop +
        log(diwish_arma(sigma_tmp, df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) +
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prop, (df_iw_prop - 3.0) * sigma_tmp)) -
        ll_cur -
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) -
        log(diwish_arma(sigma_tmp, df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j)));
      if(R::runif(0, 1) < exp(mh_log)){
        ////Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        L_mat.col(j) = L_vec_prop;
        sigma_cube_it.slice(j) = sigma_tmp;
        //L_mat = L_mat_tmp;
      }
    }

    //  Updating beta
    beta_tmp =  std::fabs(beta_it + R::runif(-.0005, .0005));
    ll_prop = arma::sum((arma::log(PP_lik_cur + beta_tmp)));
    mh_log = ll_prop - beta_tmp * LM_slice + R::dgamma(beta_tmp, 0.01, 1.0 / 0.01, 1) -
      (ll_cur - beta_it * LM_slice + R::dgamma(beta_it, 0.01, 1.0 / 0.01, 1));
    if(R::runif(0, 1) < exp(mh_log)){
      ////Rcpp::Rcout <<"MH accept" << std::endl;
      ll_cur = ll_prop;
      beta_it = beta_tmp;
    }

    //  Updating mu_alpha

    mean_mu_post = ((mean_mu_alpha / var_mu_alpha) + (arma::sum(log_alpha_it) / var_log_alpha)) /
      ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    var_mu_post = 1.0 / ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    mu_alpha_it = R::rnorm(mean_mu_post, sqrt(var_mu_post));

    //Putting new values in the linked list
    update_node_values(head, centers_it, log_alpha_it, sigma_cube_it, L_mat);

    //  Updating death rates
    if(n_cent_it > 0){
      DR_vec = calc_DR_vec2(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_bd_surf, log_pen_val, prior_n_cent_log);
    }

    // Storing MCMC samples
    centers_sample(i, 0) = centers_it;
    log_alpha_sample(i, 0) = log_alpha_it.t();
    sigma_sample(i, 0) = sigma_cube_it;
    beta_sample(i) = beta_it;
    mu_alpha_sample(i) = mu_alpha_it;
    sample_n_cent(i) = n_cent_it;
    sample_cum_int(i) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  }
  return Rcpp::List::create(Rcpp::Named("log_alpha_sample") = log_alpha_sample,
                            Rcpp::Named("centers_sample") = centers_sample,
                            Rcpp::Named("sigmas_sample") = sigma_sample,
                            Rcpp::Named("mu_alpha_sample") = mu_alpha_sample,
                            Rcpp::Named("beta_sample") = beta_sample,
                            Rcpp::Named("n_centers_sample") = sample_n_cent,
                            Rcpp::Named("cumulative_intensity_sample") = sample_cum_int);
}





//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Doing a version for continuous points (i.e. not pixelated)
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

//' Bayesian SNCP fit using BD-MCMC
//'
//' Run a BD-MCMC chain for a SNCP (mvn dispersion density) on lung ct data using a uniform proposal surface for the BD process
//'
//'
//'
//'
//'
//' @param obs_points matrix of coordinates for observed points
//' @param mean_mu_alpha prior mean for mu_alpha
//' @param sd_log_alpha prior std. dev. for log_alphas
//' @param sd_prop_alpha std. dev. for random walk proposal for log_alphas
//' @param beta prior estimate of homogeneous Poisson intensity function for cluster centers
//' @param n_it number of iterations to run MCMC
//' @param window_hw half-width of random walk cube for proposing new cluster center locations
//' @param df_iw_prior degrees of freedom for inverse-wishart prior on dispersion matrices
//' @param df_iw_prop degrees of freedom for inverse-wishart proposal on dispersion matrices
//' @param sigma_prior prior for dispersion matrices sigma
//' @param xwin vector that has the min and max x values of the observation window
//' @param ywin vector that has the min and max y values of the observation window
//' @param var_mu_alpha prior variance on mu_alpha
//' @param pen_dist distance for strauss process repulsion
//' @param pen_val vaule between 0 and 1 that is the penalty for 2 cluster centers being within pen_dist of each other
//' @param n_cent_init number of clusters to start with
//' @param prior_n_cent prior value for number of clusters
//' @param max_bd_events max events (births + deaths) to allow at each iteration of BD-MCMC
//' @param max_bd_vt max ammount of virtual time to spend in BD process at each BD_MCMC iteration
//'
//' @author Brian Vestal
//'
//' @return
//' Returns a named list containing posterior samples of model parameters with the following elements:
//' \itemize{
//' \item log_alpha_sample = list whose elements are vectors containing the estimated log-alpha estimates for the clusters present in that iteration
//' \item centers_sample = list whose elements matrices containing the estimated center locations for the clusters present in that iteration
//' \item sigmas_sample = list whose elements are lists containing the estimated covariance matrix estimates for the clusters present in that iteration
//' \item mu_alpha_sample = vector of mu_alpha samples
//' \item eta_sample = vector of eta samples,
//' \item n_centers_sample = vector that contains the number of estimated clusters at each iteration,
//' \item cumulative_intensity_sample = vector that contains the cumulative intensity estimate for the SNCP at each iteration
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List sncp_mcmc_cont2(arma::mat obs_points,
                           arma::mat centers,
                           double mean_mu_alpha,
                           double sd_log_alpha,
                           double sd_prop_alpha,
                           double beta,
                           int n_it,
                           double window_hw,
                           int df_iw_prior,
                           int df_iw_prop,
                           arma::mat sigma_prior,
                           arma::vec xwin,
                           arma::vec ywin,
                           double var_mu_alpha,
                           double pen_dist,
                           double pen_val,
                           double prior_n_cent,
                           int max_bd_events,
                           double max_bd_vt){
  int i, j, xdim = obs_points.n_cols, n_cent_it, n_points = obs_points.n_rows, n_cent_init = centers.n_rows;
  double ll_cur, ll_prop, mh_log, log_alpha_prop, beta_it = beta, mu_alpha_it, log_pen_val = log(pen_val);
  arma::vec log_alpha_it, beta_sample(n_it), mu_alpha_sample(n_it);
  //arma::mat log_alpha_sample(n_it, n_cent);
  arma::umat P_mat;
  arma::uvec P_vec_tmp, bd_idx_inf;
  arma::rowvec center_prop;
  arma::mat centers_it(n_cent_init, xdim), L_mat, L_mat_tmp, centers_tmp, sigma_tmp, L_mat_tmp_BD;
  arma::field<arma::mat> centers_sample(n_it, 1);
  arma::mat log_alpha_sample(n_it, n_cent_init);
  arma::cube sigma_cube_it(xdim, xdim, n_cent_init);
  arma::field<arma::cube> sigma_sample(n_it, 1);
  double beta_tmp, var_log_alpha = pow(sd_log_alpha, 2), mean_mu_post, var_mu_post, LM_slice = (xwin(1) - xwin(0)) * (ywin(1) - ywin(0));
  arma::vec idx_centers_init, DR_vec(n_cent_init), log_alpha_tmp_BD, PP_lik_cur, PP_lik_prop, L_vec_prop, bd_prob_vec;
  double log_bd_surf = log(1.0 / LM_slice), DR_tot, prior_n_cent_log = log(prior_n_cent / LM_slice);
  bool bd_flag;
  double bd_vt, BR_tot = 1;
  int n_bd_events, bd_kill_idx;
  struct Node *head = new Node;
  Node *kill_node;
  arma::vec sample_n_cent(n_it), sample_cum_int(n_it);
  arma::rowvec cent_add(2);
  arma::vec bd_event_vec(n_it), vt_vec(n_it);
  //Rcpp::Rcout << "Initialization of variables OK" << std::endl;

  // xwin(0) = lung_data.col(0).min();
  // xwin(1) = lung_data.col(0).max();
  // ywin(0) = lung_data.col(1).min();
  // ywin(1) = lung_data.col(1).max();
  // double x_range = xwin(1) - xwin(0), y_range = ywin(1) - ywin(0);
  pen_val = 1;
  pen_dist = 0;
  mu_alpha_it = R::rnorm(mean_mu_alpha, sqrt(var_mu_alpha));
  initNode(head, centers.row(0), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  for(i = 1; i < n_cent_init; i++){
    addNode(head, centers.row(i), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  }
  //update_Ps(head, n_cent_init, pen_dist, xdim);
  //log_alpha_it = arma::randn(n_cent_init) * sd_log_alpha + mean_mu_alpha;


  sigma_cube_it = get_sigmas(head, n_cent_init, xdim);
  log_alpha_it = get_log_alphas(head, n_cent_init);
  centers_it = get_centers(head, n_cent_init, xdim);

  //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
  //Rcpp::Rcout << "log_alphas_it = " << log_alpha_it << std::endl;
  //Rcpp::Rcout << "sigma_cube_it = " << sigma_cube_it << std::endl;

  sigma_sample(0, 0) = sigma_cube_it;
  log_alpha_sample.row(0) = log_alpha_it.t();
  centers_sample(0, 0) = centers_it;
  mu_alpha_sample(0) = mu_alpha_it;
  beta_sample(0) = beta_it;
  sample_n_cent(0) = n_cent_init;
  sample_cum_int(0) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;


  //L_mat = dmvnrm_vec_arma_cube(obs_points, centers_it, sigma_cube_it);
  //Rcpp::Rcout << "L_mat creation OK" << std::endl;
  PP_lik_cur = calc_PP_lik_vec(head);
  ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
  //Rcpp::Rcout << "ll_cur OK" << std::endl;
  //P_mat = P_mat_gen1(centers_it, pen_dist);
  //Rcpp::Rcout << "P_mat OK" << std::endl;

  //  Initial BD death rates
  //DR_vec = calc_DR_vec(head, n_cent_init, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);

  //Rcpp::Rcout << "Initial death rates OK" << std::endl;
  //Rcpp::Rcout << "DR_vec = " << DR_vec << std::endl;

  n_cent_it = n_cent_init;
  L_mat.set_size(n_points, n_cent_it);
  L_mat = get_L_mat(head, n_cent_it, n_points);
  //Rcpp::Rcout << "Storage of initials OK" << std::endl;
  //   Birth-death process:
  for(i = 1; i < n_it; ++i){
    bd_flag = false;
    bd_vt = 0;
    n_bd_events = 0;
    while(bd_flag & (n_bd_events < max_bd_events)){
      //Rcpp::Rcout << "n_bd_events = " << n_bd_events << std::endl;
      //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
      //Rcpp::Rcout << "log_alpha_it = " << log_alpha_it << std::endl;
      //display(head);
      DR_tot = arma::sum(DR_vec);

      //Rcpp::Rcout << "DR_tot = " << DR_tot << std::endl;
      //Rcpp::Rcout << "DR_vec = " << DR_tot << std::endl;
      if(DR_vec.has_inf()){
        bd_idx_inf = arma::find_nonfinite(DR_vec);
        //Rcpp::Rcout << "Infinite DR, bd_idx_inf = " << bd_idx_inf << std::endl;
        if(bd_idx_inf.n_elem < 2){
          bd_kill_idx = bd_idx_inf(0);
        }
        else{
          bd_kill_idx = RcppArmadillo::sample(bd_idx_inf, 1, false)(0);
        }

        kill_node = searchNode(head, bd_kill_idx);
        if(n_cent_it < 2){
          PP_lik_cur.zeros();
        }
        else{
          PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
        }
        ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
        deleteNode(&head, kill_node);
        n_bd_events++;
        n_cent_it--;
        //Rcpp::Rcout << "Node killed succesfully (inf)" << std::endl;
      }
      else{
        bd_vt += R::rexp(1.0 / (BR_tot + DR_tot));
        //Rcpp::Rcout << "bd_vt = " << bd_vt << std::endl;
        if(bd_vt < max_bd_vt){
          if(R::runif(0, BR_tot + DR_tot) < BR_tot){
            // idx_birth = RcppArmadillo::sample(lung_data_idx, 1, false)(0);
            cent_add(0) = R::runif(xwin(0), xwin(1));
            cent_add(1) = R::runif(ywin(0), ywin(1));
            insertFront(&head, cent_add, mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
            PP_lik_cur += (head -> alpha) * (head -> mvn_dens);
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            n_bd_events++;
            n_cent_it++;
            //Rcpp::Rcout << "Node birthed succesfully" << std::endl;
          }
          else{
            bd_prob_vec = DR_vec / DR_tot;
            bd_kill_idx = RcppArmadillo::sample(arma::linspace(0, n_cent_it - 1, n_cent_it), 1, false, bd_prob_vec)(0);
            kill_node = searchNode(head, bd_kill_idx);
            if(n_cent_it < 2){
              PP_lik_cur.zeros();
            }
            else{
              PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
            }
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            deleteNode(&head, kill_node);
            n_bd_events++;
            n_cent_it--;
            //Rcpp::Rcout << "Node killed succesfully" << std::endl;
          }
        }
        else{
          bd_flag = false;
        }
      }
      if(bd_flag){
        if(n_cent_it < 1){
          DR_vec.set_size(1);
          DR_vec(0) = 0;
        }
        else{
          update_Ps(head, n_cent_it, pen_dist, xdim);
          //Rcpp::Rcout << "Ps updated" << std::endl;
          DR_vec.set_size(n_cent_it);
          DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
          //Rcpp::Rcout << "DRs updated succesfully" << std::endl;
        }
      }
    }

    //Rcpp::Rcout << "Left BD Process" << std::endl;
    bd_event_vec(i) = n_bd_events;
    vt_vec(i) = bd_vt;

    // L_mat = get_L_mat(head, n_cent_it, n_points);
    //
    // log_alpha_it.set_size(n_cent_it);
    // log_alpha_it = get_log_alphas(head, n_cent_it);
    //
    // centers_it.set_size(n_cent_it, xdim);
    // centers_it = get_centers(head, n_cent_it, xdim);
    // P_mat = P_mat_gen1(centers_it, pen_dist);
    //
    // sigma_cube_it.set_size(xdim, xdim, n_cent_it);
    // sigma_cube_it = get_sigmas(head, n_cent_it, xdim);

    //PP_lik_cur = calc_PP_lik_vec(head);
    //ll_cur = arma::sum(arma::log(PP_lik_cur + beta_it));
    for(j = 0; j < n_cent_it; j++){
      //log_alpha_temp = log_alpha_it;
      log_alpha_prop = R::rnorm(log_alpha_it(j), sd_prop_alpha);
      //log_alpha_temp(j) = log_alpha_prop;
      PP_lik_prop = PP_lik_cur + (exp(log_alpha_prop) - exp(log_alpha_it(j))) * L_mat.col(j);
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop - ll_cur -
        exp(log_alpha_prop) + exp(log_alpha_it(j)) +
        (1.0/(2.0 * pow(sd_log_alpha, 2))) * (pow(log_alpha_it(j) - mean_mu_alpha, 2) - pow(log_alpha_prop - mean_mu_alpha, 2));
      if(R::runif(0, 1) < exp(mh_log)){
        //Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        log_alpha_it(j) = log_alpha_prop;
      }
    }

    //Rcpp::Rcout << "doing centers" << std::endl;
    // Doing center location MH steps
    for(j = 0; j < n_cent_it; j++){
      center_prop = centers_it.row(j);
      center_prop(0) += R::runif(-window_hw, window_hw);
      center_prop(1) += R::runif(-window_hw, window_hw);
      if(center_prop(0) < xwin(1) & center_prop(0) > xwin(0) &
         center_prop(1) < ywin(1) & center_prop(1) > ywin(0)){
        //L_mat_tmp = L_mat;
        //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        L_vec_prop = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
        //P_vec_tmp = P_mat_gen2(centers_it, center_prop, pen_dist);
        //P_vec_tmp.shed_row(j);
        PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
        ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
        mh_log = ll_prop - ll_cur;
        if(R::runif(0, 1) < exp(mh_log)){
          PP_lik_cur = PP_lik_prop;
          ll_cur = ll_prop;
          centers_it.row(j) = center_prop;
          L_mat.col(j) = L_vec_prop;
          //P_mat = P_mat_gen1(centers_it, pen_dist);
        }
      }
    }

    //Rcpp::Rcout << "doing sigmas" << std::endl;
    // Doing sigma matrix MH steps
    for(j = 0; j < n_cent_it; j++){
      sigma_tmp = riwish_arma(df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j));
      /*
       log(diwish(W = sigma_mat_it[[j]], v = df_inv_wish_prior,
       S = (df_inv_wish_prior - 3) * S_mat_prior)) + ## Prior
       log(diwish(W = sigma_prop_j, v = df_inv_wish_prop,
       S = (df_inv_wish_prop - 3) * sigma_mat_it[[j]]))
       */
      //L_mat_tmp = L_mat;
      //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      L_vec_prop = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop +
        log(diwish_arma(sigma_tmp, df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) +
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prop, (df_iw_prop - 3.0) * sigma_tmp)) -
        ll_cur -
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) -
        log(diwish_arma(sigma_tmp, df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j)));
      if(R::runif(0, 1) < exp(mh_log)){
        ////Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        L_mat.col(j) = L_vec_prop;
        sigma_cube_it.slice(j) = sigma_tmp;
        //L_mat = L_mat_tmp;
      }
    }

    //  Updating beta
    //Rcpp::Rcout << "doing beta" << std::endl;
    beta_tmp =  std::fabs(beta_it + R::runif(-.0005, .0005));
    ll_prop = arma::sum((arma::log(PP_lik_cur + beta_tmp)));
    mh_log = ll_prop - beta_tmp * LM_slice + R::dgamma(beta_tmp, 0.01, 1.0 / 0.01, 1) -
      (ll_cur - beta_it * LM_slice + R::dgamma(beta_it, 0.01, 1.0 / 0.01, 1));
    if(R::runif(0, 1) < exp(mh_log)){
      ////Rcpp::Rcout <<"MH accept" << std::endl;
      ll_cur = ll_prop;
      beta_it = beta_tmp;
    }

    //  Updating mu_alpha

    //Rcpp::Rcout << "doing mu_alpha" << std::endl;
    mean_mu_post = ((mean_mu_alpha / var_mu_alpha) + (arma::sum(log_alpha_it) / var_log_alpha)) /
      ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    var_mu_post = 1.0 / ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    mu_alpha_it = R::rnorm(mean_mu_post, sqrt(var_mu_post));

    //Putting new values in the linked list
    // update_node_values(head, centers_it, log_alpha_it, sigma_cube_it, L_mat);
    //
    // //  Updating death rates
    // if(n_cent_it > 0){
    //   DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
    // }

    // Storing MCMC samples
    //Rcpp::Rcout << "storing results" << std::endl;
    centers_sample(i, 0) = centers_it;
    log_alpha_sample.row(i) = log_alpha_it.t();
    sigma_sample(i, 0) = sigma_cube_it;
    beta_sample(i) = beta_it;
    mu_alpha_sample(i) = mu_alpha_it;
    sample_n_cent(i) = n_cent_it;
    sample_cum_int(i) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;

  }
  return Rcpp::List::create(Rcpp::Named("log_alpha_sample") = log_alpha_sample,
                            Rcpp::Named("centers_sample") = centers_sample,
                            Rcpp::Named("sigmas_sample") = sigma_sample,
                            Rcpp::Named("mu_alpha_sample") = mu_alpha_sample,
                            Rcpp::Named("beta_sample") = beta_sample,
                            Rcpp::Named("n_centers_sample") = sample_n_cent,
                            Rcpp::Named("cumulative_intensity_sample") = sample_cum_int,
                            Rcpp::Named("n_bd_events") = bd_event_vec,
                            Rcpp::Named("bd_vt") = vt_vec);
}

//////////////////////////////////////////////////////////////////////////////

//' Bayesian SNCP fit using BD-MCMC
//'
//' Run a BD-MCMC chain for a SNCP (mvn dispersion density) on lung ct data using a uniform proposal surface for the BD process
//'
//'
//'
//'
//'
//' @param obs_points matrix of coordinates for observed points
//' @param mean_mu_alpha prior mean for mu_alpha
//' @param sd_log_alpha prior std. dev. for log_alphas
//' @param sd_prop_alpha std. dev. for random walk proposal for log_alphas
//' @param beta prior estimate of homogeneous Poisson intensity function for cluster centers
//' @param n_it number of iterations to run MCMC
//' @param window_hw half-width of random walk cube for proposing new cluster center locations
//' @param df_iw_prior degrees of freedom for inverse-wishart prior on dispersion matrices
//' @param df_iw_prop degrees of freedom for inverse-wishart proposal on dispersion matrices
//' @param sigma_prior prior for dispersion matrices sigma
//' @param xwin vector that has the min and max x values of the observation window
//' @param ywin vector that has the min and max y values of the observation window
//' @param var_mu_alpha prior variance on mu_alpha
//' @param pen_dist distance for strauss process repulsion
//' @param pen_val vaule between 0 and 1 that is the penalty for 2 cluster centers being within pen_dist of each other
//' @param n_cent_init number of clusters to start with
//' @param prior_n_cent prior value for number of clusters
//' @param max_bd_events max events (births + deaths) to allow at each iteration of BD-MCMC
//' @param max_bd_vt max ammount of virtual time to spend in BD process at each BD_MCMC iteration
//'
//' @author Brian Vestal
//'
//' @return
//' Returns a named list containing posterior samples of model parameters with the following elements:
//' \itemize{
//' \item log_alpha_sample = list whose elements are vectors containing the estimated log-alpha estimates for the clusters present in that iteration
//' \item centers_sample = list whose elements matrices containing the estimated center locations for the clusters present in that iteration
//' \item sigmas_sample = list whose elements are lists containing the estimated covariance matrix estimates for the clusters present in that iteration
//' \item mu_alpha_sample = vector of mu_alpha samples
//' \item eta_sample = vector of eta samples,
//' \item n_centers_sample = vector that contains the number of estimated clusters at each iteration,
//' \item cumulative_intensity_sample = vector that contains the cumulative intensity estimate for the SNCP at each iteration
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List sncp_mcmc_cont_fixed(arma::mat obs_points,
                           arma::mat centers,
                           double mean_mu_alpha,
                           double sd_log_alpha,
                           double sd_prop_alpha,
                           double beta,
                           int n_it,
                           double window_hw,
                           int df_iw_prior,
                           int df_iw_prop,
                           arma::mat sigma_prior,
                           arma::vec xwin,
                           arma::vec ywin,
                           double var_mu_alpha,
                           double pen_dist,
                           double pen_val,
                           double prior_n_cent,
                           int max_bd_events,
                           double max_bd_vt){
  int i, j, xdim = obs_points.n_cols, n_cent_it, n_points = obs_points.n_rows, n_cent_init = centers.n_rows;
  double ll_cur, ll_prop, mh_log, log_alpha_prop, beta_it = beta, mu_alpha_it, log_pen_val = log(pen_val);
  arma::vec log_alpha_it, beta_sample(n_it), mu_alpha_sample(n_it);
  //arma::mat log_alpha_sample(n_it, n_cent);
  arma::umat P_mat;
  arma::uvec P_vec_tmp, bd_idx_inf;
  arma::rowvec center_prop;
  arma::mat centers_it(n_cent_init, xdim), L_mat, L_mat_tmp, centers_tmp, sigma_tmp, L_mat_tmp_BD;
  arma::field<arma::mat> centers_sample(n_it, 1);
  arma::mat log_alpha_sample(n_it, n_cent_init);
  arma::cube sigma_cube_it(xdim, xdim, n_cent_init);
  arma::field<arma::cube> sigma_sample(n_it, 1);
  double beta_tmp, var_log_alpha = pow(sd_log_alpha, 2), mean_mu_post, var_mu_post, LM_slice = (xwin(1) - xwin(0)) * (ywin(1) - ywin(0));
  arma::vec idx_centers_init, DR_vec(n_cent_init), log_alpha_tmp_BD, PP_lik_cur, PP_lik_prop, L_vec_prop, bd_prob_vec;
  double log_bd_surf = log(1.0 / LM_slice), DR_tot, prior_n_cent_log = log(prior_n_cent / LM_slice);
  bool bd_flag;
  double bd_vt, BR_tot = 1;
  int n_bd_events, bd_kill_idx;
  struct Node *head = new Node;
  Node *kill_node;
  arma::vec sample_n_cent(n_it), sample_cum_int(n_it);
  arma::rowvec cent_add(2);
  arma::vec bd_event_vec(n_it), vt_vec(n_it);
  //Rcpp::Rcout << "Initialization of variables OK" << std::endl;

  // xwin(0) = lung_data.col(0).min();
  // xwin(1) = lung_data.col(0).max();
  // ywin(0) = lung_data.col(1).min();
  // ywin(1) = lung_data.col(1).max();
  // double x_range = xwin(1) - xwin(0), y_range = ywin(1) - ywin(0);
  pen_val = 1;
  pen_dist = 0;
  log_pen_val = 0;
  mu_alpha_it = R::rnorm(mean_mu_alpha, sqrt(var_mu_alpha));
  initNode(head, centers.row(0), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  for(i = 1; i < n_cent_init; i++){
    addNode(head, centers.row(i), mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
  }
  //update_Ps(head, n_cent_init, pen_dist, xdim);
  //log_alpha_it = arma::randn(n_cent_init) * sd_log_alpha + mean_mu_alpha;


  sigma_cube_it = get_sigmas(head, n_cent_init, xdim);
  log_alpha_it = get_log_alphas(head, n_cent_init);
  centers_it = get_centers(head, n_cent_init, xdim);

  //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
  //Rcpp::Rcout << "log_alphas_it = " << log_alpha_it << std::endl;
  //Rcpp::Rcout << "sigma_cube_it = " << sigma_cube_it << std::endl;

  sigma_sample(0, 0) = sigma_cube_it;
  log_alpha_sample.row(0) = log_alpha_it.t();
  centers_sample(0, 0) = centers_it;
  mu_alpha_sample(0) = mu_alpha_it;
  beta_sample(0) = beta_it;
  sample_n_cent(0) = n_cent_init;
  sample_cum_int(0) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;


  //L_mat = dmvnrm_vec_arma_cube(obs_points, centers_it, sigma_cube_it);
  //Rcpp::Rcout << "L_mat creation OK" << std::endl;
  PP_lik_cur = calc_PP_lik_vec(head);
  ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
  //Rcpp::Rcout << "ll_cur OK" << std::endl;
  //P_mat = P_mat_gen1(centers_it, pen_dist);
  //Rcpp::Rcout << "P_mat OK" << std::endl;

  //  Initial BD death rates
  DR_vec = calc_DR_vec(head, n_cent_init, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);

  arma::mat DR_sample(n_it, n_cent_init);
  DR_sample.row(0) = DR_vec.t();
  //Rcpp::Rcout << "Initial death rates OK" << std::endl;
  //Rcpp::Rcout << "DR_vec = " << DR_vec << std::endl;

  n_cent_it = n_cent_init;
  L_mat.set_size(n_points, n_cent_it);
  L_mat = get_L_mat(head, n_cent_it, n_points);
  //Rcpp::Rcout << "Storage of initials OK" << std::endl;
  //   Birth-death process:
  for(i = 1; i < n_it; ++i){
    bd_flag = false;
    bd_vt = 0;
    n_bd_events = 0;
    while(bd_flag & (n_bd_events < max_bd_events)){
      //Rcpp::Rcout << "n_bd_events = " << n_bd_events << std::endl;
      //Rcpp::Rcout << "centers_it = " << centers_it << std::endl;
      //Rcpp::Rcout << "log_alpha_it = " << log_alpha_it << std::endl;
      //display(head);
      DR_tot = arma::sum(DR_vec);

      //Rcpp::Rcout << "DR_tot = " << DR_tot << std::endl;
      //Rcpp::Rcout << "DR_vec = " << DR_tot << std::endl;
      if(DR_vec.has_inf()){
        bd_idx_inf = arma::find_nonfinite(DR_vec);
        //Rcpp::Rcout << "Infinite DR, bd_idx_inf = " << bd_idx_inf << std::endl;
        if(bd_idx_inf.n_elem < 2){
          bd_kill_idx = bd_idx_inf(0);
        }
        else{
          bd_kill_idx = RcppArmadillo::sample(bd_idx_inf, 1, false)(0);
        }

        kill_node = searchNode(head, bd_kill_idx);
        if(n_cent_it < 2){
          PP_lik_cur.zeros();
        }
        else{
          PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
        }
        ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
        deleteNode(&head, kill_node);
        n_bd_events++;
        n_cent_it--;
        //Rcpp::Rcout << "Node killed succesfully (inf)" << std::endl;
      }
      else{
        bd_vt += R::rexp(1.0 / (BR_tot + DR_tot));
        //Rcpp::Rcout << "bd_vt = " << bd_vt << std::endl;
        if(bd_vt < max_bd_vt){
          if(R::runif(0, BR_tot + DR_tot) < BR_tot){
            // idx_birth = RcppArmadillo::sample(lung_data_idx, 1, false)(0);
            cent_add(0) = R::runif(xwin(0), xwin(1));
            cent_add(1) = R::runif(ywin(0), ywin(1));
            insertFront(&head, cent_add, mu_alpha_it, sd_log_alpha, sigma_prior, df_iw_prior, obs_points, log_bd_surf);
            PP_lik_cur += (head -> alpha) * (head -> mvn_dens);
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            n_bd_events++;
            n_cent_it++;
            //Rcpp::Rcout << "Node birthed succesfully" << std::endl;
          }
          else{
            bd_prob_vec = DR_vec / DR_tot;
            bd_kill_idx = RcppArmadillo::sample(arma::linspace(0, n_cent_it - 1, n_cent_it), 1, false, bd_prob_vec)(0);
            kill_node = searchNode(head, bd_kill_idx);
            if(n_cent_it < 2){
              PP_lik_cur.zeros();
            }
            else{
              PP_lik_cur -= kill_node -> alpha * kill_node -> mvn_dens;
            }
            ll_cur = arma::sum((arma::log(PP_lik_cur + beta_it)));
            deleteNode(&head, kill_node);
            n_bd_events++;
            n_cent_it--;
            //Rcpp::Rcout << "Node killed succesfully" << std::endl;
          }
        }
        else{
          bd_flag = false;
        }
      }
      if(bd_flag){
        if(n_cent_it < 1){
          DR_vec.set_size(1);
          DR_vec(0) = 0;
        }
        else{
          update_Ps(head, n_cent_it, pen_dist, xdim);
          //Rcpp::Rcout << "Ps updated" << std::endl;
          DR_vec.set_size(n_cent_it);
          DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
          //Rcpp::Rcout << "DRs updated succesfully" << std::endl;
        }
      }
    }

    //Rcpp::Rcout << "Left BD Process" << std::endl;
    bd_event_vec(i) = n_bd_events;
    vt_vec(i) = bd_vt;

    // L_mat = get_L_mat(head, n_cent_it, n_points);
    //
    // log_alpha_it.set_size(n_cent_it);
    // log_alpha_it = get_log_alphas(head, n_cent_it);
    //
    // centers_it.set_size(n_cent_it, xdim);
    // centers_it = get_centers(head, n_cent_it, xdim);
    // P_mat = P_mat_gen1(centers_it, pen_dist);
    //
    // sigma_cube_it.set_size(xdim, xdim, n_cent_it);
    // sigma_cube_it = get_sigmas(head, n_cent_it, xdim);

    //PP_lik_cur = calc_PP_lik_vec(head);
    //ll_cur = arma::sum(arma::log(PP_lik_cur + beta_it));
    for(j = 0; j < n_cent_it; j++){
      //log_alpha_temp = log_alpha_it;
      log_alpha_prop = R::rnorm(log_alpha_it(j), sd_prop_alpha);
      //log_alpha_temp(j) = log_alpha_prop;
      PP_lik_prop = PP_lik_cur + (exp(log_alpha_prop) - exp(log_alpha_it(j))) * L_mat.col(j);
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop - ll_cur -
        exp(log_alpha_prop) + exp(log_alpha_it(j)) +
        (1.0/(2.0 * pow(sd_log_alpha, 2))) * (pow(log_alpha_it(j) - mean_mu_alpha, 2) - pow(log_alpha_prop - mean_mu_alpha, 2));
      if(R::runif(0, 1) < exp(mh_log)){
        //Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        log_alpha_it(j) = log_alpha_prop;
      }
    }

    //Rcpp::Rcout << "doing centers" << std::endl;
    // Doing center location MH steps
    // for(j = 0; j < n_cent_it; j++){
    //   center_prop = centers_it.row(j);
    //   center_prop(0) += R::runif(-window_hw, window_hw);
    //   center_prop(1) += R::runif(-window_hw, window_hw);
    //   if(center_prop(0) < xwin(1) & center_prop(0) > xwin(0) &
    //      center_prop(1) < ywin(1) & center_prop(1) > ywin(0)){
    //     //L_mat_tmp = L_mat;
    //     //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
    //     L_vec_prop = dmvnrm_vec_arma_1f(obs_points, center_prop.t(), sigma_cube_it.slice(j));
    //     //P_vec_tmp = P_mat_gen2(centers_it, center_prop, pen_dist);
    //     //P_vec_tmp.shed_row(j);
    //     PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
    //     ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
    //     mh_log = ll_prop - ll_cur;
    //     if(R::runif(0, 1) < exp(mh_log)){
    //       PP_lik_cur = PP_lik_prop;
    //       ll_cur = ll_prop;
    //       centers_it.row(j) = center_prop;
    //       L_mat.col(j) = L_vec_prop;
    //       //P_mat = P_mat_gen1(centers_it, pen_dist);
    //     }
    //   }
    // }

    //Rcpp::Rcout << "doing sigmas" << std::endl;
    // Doing sigma matrix MH steps
    for(j = 0; j < n_cent_it; j++){
      sigma_tmp = riwish_arma(df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j));
      /*
      log(diwish(W = sigma_mat_it[[j]], v = df_inv_wish_prior,
                 S = (df_inv_wish_prior - 3) * S_mat_prior)) + ## Prior
      log(diwish(W = sigma_prop_j, v = df_inv_wish_prop,
                 S = (df_inv_wish_prop - 3) * sigma_mat_it[[j]]))
      */
      //L_mat_tmp = L_mat;
      //L_mat_tmp.col(j) = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      L_vec_prop = dmvnrm_vec_arma_1f(obs_points, centers_it.row(j).t(), sigma_tmp);
      PP_lik_prop = PP_lik_cur + exp(log_alpha_it(j)) * (L_vec_prop - L_mat.col(j));
      ll_prop = arma::sum((arma::log(PP_lik_prop + beta_it)));
      mh_log = ll_prop +
        log(diwish_arma(sigma_tmp, df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) +
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prop, (df_iw_prop - 3.0) * sigma_tmp)) -
        ll_cur -
        log(diwish_arma(sigma_cube_it.slice(j), df_iw_prior, (df_iw_prior - 3.0) * sigma_prior)) -
        log(diwish_arma(sigma_tmp, df_iw_prop, (df_iw_prop - 3.0) * sigma_cube_it.slice(j)));
      if(R::runif(0, 1) < exp(mh_log)){
        ////Rcpp::Rcout <<"MH accept" << std::endl;
        PP_lik_cur = PP_lik_prop;
        ll_cur = ll_prop;
        L_mat.col(j) = L_vec_prop;
        sigma_cube_it.slice(j) = sigma_tmp;
        //L_mat = L_mat_tmp;
      }
    }

    //  Updating beta
    //Rcpp::Rcout << "doing beta" << std::endl;
    beta_tmp =  std::fabs(beta_it + R::runif(-.0005, .0005));
    ll_prop = arma::sum((arma::log(PP_lik_cur + beta_tmp)));
    mh_log = ll_prop - beta_tmp * LM_slice + R::dgamma(beta_tmp, 0.01, 1.0 / 0.01, 1) -
      (ll_cur - beta_it * LM_slice + R::dgamma(beta_it, 0.01, 1.0 / 0.01, 1));
    if(R::runif(0, 1) < exp(mh_log)){
      ////Rcpp::Rcout <<"MH accept" << std::endl;
      ll_cur = ll_prop;
      beta_it = beta_tmp;
    }

    //  Updating mu_alpha

    //Rcpp::Rcout << "doing mu_alpha" << std::endl;
    mean_mu_post = ((mean_mu_alpha / var_mu_alpha) + (arma::sum(log_alpha_it) / var_log_alpha)) /
      ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    var_mu_post = 1.0 / ((1.0 / var_mu_alpha) + (n_cent_it / var_log_alpha));
    mu_alpha_it = R::rnorm(mean_mu_post, sqrt(var_mu_post));

    //Putting new values in the linked list
    update_node_values(head, centers_it, log_alpha_it, sigma_cube_it, L_mat);
    //
    // //  Updating death rates
    if(n_cent_it > 0){
      DR_vec = calc_DR_vec(head, n_cent_it, PP_lik_cur, ll_cur, beta_it, n_points, log_pen_val, prior_n_cent_log);
    }

    // Storing MCMC samples
    //Rcpp::Rcout << "storing results" << std::endl;
    centers_sample(i, 0) = centers_it;
    log_alpha_sample.row(i) = log_alpha_it.t();
    sigma_sample(i, 0) = sigma_cube_it;
    beta_sample(i) = beta_it;
    mu_alpha_sample(i) = mu_alpha_it;
    sample_n_cent(i) = n_cent_it;
    sample_cum_int(i) = arma::sum(arma::exp(log_alpha_it)) + beta_it * LM_slice;
    DR_sample.row(i) = DR_vec.t();

  }
  return Rcpp::List::create(Rcpp::Named("log_alpha_sample") = log_alpha_sample,
                            Rcpp::Named("centers_sample") = centers_sample,
                            Rcpp::Named("sigmas_sample") = sigma_sample,
                            Rcpp::Named("mu_alpha_sample") = mu_alpha_sample,
                            Rcpp::Named("beta_sample") = beta_sample,
                            Rcpp::Named("n_centers_sample") = sample_n_cent,
                            Rcpp::Named("DR_sample") = DR_sample,
                            Rcpp::Named("cumulative_intensity_sample") = sample_cum_int,
                            Rcpp::Named("n_bd_events") = bd_event_vec,
                            Rcpp::Named("bd_vt") = vt_vec);
}

//' Fast p-dist2 function
//'
//' Run p-dist
//'
//' This is where you write details on the function...
//'
//' more details....
//'
//' @param A matrix of points
//' @param B matrix of points
//'
//' @author Brian Vestal
//'
//' @return
//' Returns a matrix with distances
//'
//' @export
// [[Rcpp::export]]
arma::mat fastPdist2(arma::mat A, arma::mat B) {

  arma::colvec An =  arma::sum(square(A),1);
  arma::colvec Bn =  arma::sum(square(B),1);

  arma::mat C = -2.0 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();

  return sqrt(abs(C));
}
