#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

//' @title Hidden Markov Model Parameter Estimation
 //' @description Estimates the parameters of a Hidden Markov Model (HMM) using the EM algorithm.
 //' @param observations A NumericMatrix where each row represents an observation at a time step, and columns are features.
 //' @param n_states An integer specifying the number of hidden states.
 //' @param dist_type A string specifying the observation distribution type. Supported values are "poisson" and "normal".
 //' @param max_iter An integer specifying the maximum number of iterations for the EM algorithm. Default is 100.
 //' @return A list containing:
 //'   - `A`: A NumericMatrix representing the state transition probability matrix.
 //'   - `pi`: A NumericVector representing the initial state distribution.
 //'   - `mu`: A NumericMatrix representing the mean of the observation distributions.
 //'   - `sigma`: A NumericMatrix representing the standard deviation of the observation distributions (only for normal distribution).
 //' @examples
 //' {
 //' observations <- matrix(runif(100), nrow = 10, ncol = 10)
 //' result <- hmm_em(observations, n_states = 3, dist_type = "normal", max_iter = 50)
 //' print(result)
 //' }
 //' @name hmm
 //' @export
 // [[Rcpp::export]]
 List hmm_em(NumericMatrix observations, int n_states, std::string dist_type = "normal", int max_iter = 100) {
   int n_obs = observations.nrow();
   int n_features = observations.ncol();
   
   // Initialize parameters
   NumericMatrix A(n_states, n_states); // Transition matrix
   NumericVector pi(n_states, 1.0 / n_states); // Initial state probabilities
   NumericMatrix mu(n_states, n_features); // Mean matrix
   NumericMatrix sigma(n_states, n_features); // Standard deviation matrix (for normal distribution)
   
   // Uniform initialization for A
   for (int i = 0; i < n_states; ++i) {
     for (int j = 0; j < n_states; ++j) {
       A(i, j) = 1.0 / n_states;
     }
   }
   
   // Random initialization for means and standard deviations
   for (int k = 0; k < n_states; ++k) {
     for (int d = 0; d < n_features; ++d) {
       mu(k, d) = R::runif(0, 10);
       sigma(k, d) = R::runif(1, 5);
     }
   }
   
   // Variables for the E-step
   NumericMatrix gamma(n_obs, n_states); // Posterior probability matrix
   
   // EM iteration
   for (int iter = 0; iter < max_iter; ++iter) {
     // --- E-step ---
     for (int t = 0; t < n_obs; ++t) {
       double row_sum = 0.0;
       for (int k = 0; k < n_states; ++k) {
         double likelihood = 1.0;
         for (int d = 0; d < n_features; ++d) {
           if (dist_type == "poisson") {
             likelihood *= R::dpois(observations(t, d), mu(k, d), false);
           } else if (dist_type == "normal") {
             likelihood *= R::dnorm(observations(t, d), mu(k, d), sigma(k, d), false);
           }
         }
         gamma(t, k) = pi[k] * likelihood;
         row_sum += gamma(t, k);
       }
       // Normalize
       for (int k = 0; k < n_states; ++k) {
         gamma(t, k) /= row_sum;
       }
     }
     
     // --- M-step ---
     // Update initial state probabilities
     for (int k = 0; k < n_states; ++k) {
       pi[k] = gamma(0, k);
     }
     
     // Update transition matrix
     for (int i = 0; i < n_states; ++i) {
       for (int j = 0; j < n_states; ++j) {
         double num = 0.0, denom = 0.0;
         for (int t = 0; t < n_obs - 1; ++t) {
           num += gamma(t, i) * gamma(t + 1, j);
           denom += gamma(t, i);
         }
         A(i, j) = num / denom;
       }
     }
     
     // Update means and standard deviations
     for (int k = 0; k < n_states; ++k) {
       for (int d = 0; d < n_features; ++d) {
         double num_mu = 0.0, denom_mu = 0.0, num_sigma = 0.0;
         for (int t = 0; t < n_obs; ++t) {
           num_mu += gamma(t, k) * observations(t, d);
           denom_mu += gamma(t, k);
         }
         mu(k, d) = num_mu / denom_mu;
         
         if (dist_type == "normal") {
           for (int t = 0; t < n_obs; ++t) {
             num_sigma += gamma(t, k) * std::pow(observations(t, d) - mu(k, d), 2);
           }
           sigma(k, d) = std::sqrt(num_sigma / denom_mu);
         }
       }
     }
   }
   
   return List::create(
     Named("A") = A,
     Named("pi") = pi,
     Named("mu") = mu,
     Named("sigma") = sigma
   );
 }
