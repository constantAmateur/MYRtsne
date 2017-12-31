#include <Rcpp.h>
#include "tsne.h"

using namespace Rcpp;

// Function that runs the Barnes-Hut implementation of t-SNE
// [[Rcpp::export]]
Rcpp::List Rtsne_cpp(NumericMatrix X, int no_dims_in, double perplexity_in, 
                     double theta_in, bool verbose, int max_iter, 
                     bool distance_precomputed, NumericMatrix Y_in, bool init, 
                     int stop_lying_iter_in, int mom_switch_iter_in,
                     double momentum_in, double final_momentum_in, 
                     double eta_in, double exaggeration_factor_in) {

  int origN, N, D, no_dims = no_dims_in;

	double  *data;
  TSNE* tsne = new TSNE();
  double perplexity = perplexity_in;
  double theta = theta_in;
  int stop_lying_iter = stop_lying_iter_in;
  int mom_switch_iter = mom_switch_iter_in;
  double momentum = momentum_in;
  double final_momentum = final_momentum_in;
  double eta = eta_in;
  double exaggeration_factor = exaggeration_factor_in;
  bool exact = (theta == .0) ? true : false;
  int K = (int) (3*perplexity);
  
  origN = X.nrow();
  D = X.ncol();
    
	data = (double*) calloc(D * origN, sizeof(double));
    if(data == NULL) { Rcpp::stop("Memory allocation failed!\n"); }
    for (int i = 0; i < origN; i++){
        for (int j = 0; j < D; j++){
            data[i*D+j] = X(i,j);
        }
    }
    //Extra things to store
    double* betas = (double*) calloc(origN, sizeof(double));
    double *Ps = (double*) calloc(origN*origN,sizeof(double));
    if(betas == NULL || Ps == NULL) {Rcpp::stop("Memory allocation failed!\n");}
    
    // Make dummy landmarks
    N = origN;
    if (verbose) Rprintf("Read the %i x %i data matrix successfully!\n", N, D);
    int* landmarks = (int*) malloc(N * sizeof(int));
    if(landmarks == NULL) { Rcpp::stop("Memory allocation failed!\n"); }
    for(int n = 0; n < N; n++) landmarks[n] = n;

		double* Y = (double*) malloc(N * no_dims * sizeof(double));
		double* costs = (double*) calloc(N, sizeof(double));
		double* itercosts = (double*) calloc((int)(ceil(max_iter/50.0)), sizeof(double));
    if(Y == NULL || costs == NULL ) { Rcpp::stop("Memory allocation failed!\n"); }
    
    // Initialize solution (randomly)
    if (init) {
      for (int i = 0; i < N; i++){
        for (int j = 0; j < no_dims; j++){
          Y[i*no_dims+j] = Y_in(i,j);
        }
      }
      if (verbose) Rprintf("Using user supplied starting positions\n");
    }
    
    // Run tsne
		tsne->run(data, N, D, Y, no_dims, perplexity, theta, verbose, max_iter, costs, distance_precomputed, 
            itercosts, init, stop_lying_iter, mom_switch_iter, momentum, final_momentum, eta, exaggeration_factor,betas,Ps);

  	// Save the results
    Rcpp::NumericMatrix Yr(N, no_dims);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < no_dims; j++){
            Yr(i,j) = Y[i*no_dims+j];
        }
    }
    
    Rcpp::NumericVector costsr(N);
    for (int i = 0; i < N; i++){
      costsr(i) = costs[i];
    }
    Rcpp::NumericVector itercostsr((int)(ceil(max_iter/50.0)));
    for (int i = 0; i < (int)(ceil(max_iter/50.0)); i++) {
      itercostsr(i) = itercosts[i];
    }
    //The extra return things
    Rcpp::NumericVector betasr(N);
    for(int i=0;i < N;i++) betasr(i) = betas[i];
    Rcpp::NumericMatrix Psr(N,N);
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        Psr(i,j) = Ps[i*N+j];
      }
    }
    //Calculate distances in resulting map
    double* DD = (double*) malloc(N*N*sizeof(double));
    if(DD==NULL){Rcpp::stop("Memory allocation failed!\n");}
    tsne->computeSquaredEuclideanDistance(Y,N,no_dims,DD);
    Rcpp::NumericMatrix Qr(N,N);
    double sum_Q = .0;
    for(int i=0; i<N ;i++){
      for(int j=0;j<N;j++){
        if(i!=j){
          Qr(i,j) = 1/(1+DD[i*N+j]);
          sum_Q += Qr(i,j);
        }
      }
    }
    ////Do normalisation
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        if(i!=j) Qr(i,j) = Qr(i,j)/sum_Q;
      }
    }
    free(DD); DD=NULL;
 
    free(Ps); Ps = NULL;
    free(data); data = NULL;
  	 free(Y); Y = NULL;
		free(costs); costs = NULL;
      free(betas); betas = NULL;
		free(landmarks); landmarks = NULL;
    delete(tsne);
    
    Rcpp::List output = Rcpp::List::create(Rcpp::_["theta"]=theta, 
                                           Rcpp::_["perplexity"]=perplexity, 
                                           Rcpp::_["N"]=N,
                                           Rcpp::_["origD"]=D,
                                           Rcpp::_["Y"]=Yr, 
                                           Rcpp::_["costs"]=costsr, 
                                           Rcpp::_["itercosts"]=itercostsr,
                                           Rcpp::_["stop_lying_iter"]=stop_lying_iter, 
                                           Rcpp::_["mom_switch_iter"]=mom_switch_iter, 
                                           Rcpp::_["momentum"]=momentum, 
                                           Rcpp::_["final_momentum"]=final_momentum, 
                                           Rcpp::_["eta"]=eta, 
                                           Rcpp::_["exaggeration_factor"]=exaggeration_factor,
                                           Rcpp::_["betas"]=betasr,
                                           Rcpp::_["P"]=Psr,
                                           Rcpp::_["Q"]=Qr);
    return output; 
}
