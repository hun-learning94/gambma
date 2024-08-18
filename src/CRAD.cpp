#include "CRAD.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::vec armapmax(const arma::vec &x,
                   const double &y){
  arma::vec res(x.n_elem);
  for(unsigned i{0}; i<x.n_elem; i++){
    if(x(i) > y){ res(i) = x(i); } else {res(i) = y;}
  }
  return res;
}

arma::mat CRADNS_1d_cpp(const arma::vec& x, const arma::vec &knot, 
                        bool knotalive, double bdmargin){
  unsigned N{x.n_elem};
  arma::mat B;
  if(knotalive){
    // augment knots with boundary knots
    arma::vec augknot(knot.n_elem+2);
    
    // augknot(0) = bdmargin;
    // augknot(augknot.n_elem-1) = 1.0 - bdmargin;
    augknot(0) = bdmargin - 0.01;
    augknot(augknot.n_elem-1) = 1.0 - bdmargin + 0.01;
    
    augknot.subvec(1, augknot.n_elem - 2) = knot;
    unsigned nk{augknot.n_elem};
    B.set_size(N, nk - 1);
    B.each_col() = 
      (-armapmax(arma::pow(arma::abs(x - augknot(nk-1)), 3), 0) + 
      armapmax(arma::pow(arma::abs(x - augknot(0)), 3), 0)) / 
      (augknot(nk-1) - augknot(0));
    for(unsigned i{0}; i < (nk-2); i++){
      // B.col(i + 1) -= 
      //   (armapmax(arma::pow(arma::abs(x - augknot(nk-1)), 3), 0) - 
      //   armapmax(arma::pow(arma::abs(x - augknot(i)), 3),0)) /
      //   (augknot(nk-1) - augknot(i));
      B.col(i + 1) -=
        (armapmax(arma::pow(arma::abs(x - augknot(nk-1)), 3), 0) -
        armapmax(arma::pow(arma::abs(x - augknot(i+1)), 3),0)) /
        (augknot(nk-1) - augknot(i+1));
    }
    B.col(0) = x;
  } else {
    B.set_size(N, 1);
    B.col(0) = x;
  }
  return B;
}


arma::mat CRAD_1d_cpp(const arma::vec& x, const arma::vec &knot, 
                      bool knotalive){
  unsigned N{x.n_elem};
  arma::mat B;
  if(knotalive){
    unsigned nk{knot.n_elem};
    B.set_size(N, nk + 3);
    for(unsigned i{0}; i < 3; i++){ B.col(i) = arma::pow(x, i+1); }
    for(unsigned i{0}; i < nk; i++){ B.col(i + 3) = armapmax(arma::pow(arma::abs(x - knot(i)), 3), 0); }
  } else {
    B.set_size(N, 3);
    for(unsigned i{0}; i < 3; i++){ B.col(i) = arma::pow(x, i+1); }
  }
  return B;
}


Rcpp::List CRAD_cpp(const arma::mat &X, 
                    const arma::mat &X_lin, 
                    const arma::vec &knots, 
                    const arma::uvec &knotsidx,
                    bool NS,
                    double bdmargin){
  unsigned p{X.n_cols};
  unsigned p_lin{X_lin.n_cols};
  unsigned N{X.n_rows};
  if(NS == 0) bdmargin = 0.0;
  
  // basis dimension:
  // ## 1. CRAD
  // 0 knot : u, u2, u3
  // 1 knot : u, u2, u3, abs(u - knot1)^3,
  // ...
  // ## 2. NS
  // 0 knot : u (1 bases)
  // 1 knot : augment 2 bd knots (e.g. 0.05, 0.95) and impose bd condition -> 3 knots - 1 (2 bases)
  // 2 knot : augment 2 bd knots (e.g. 0.05, 0.95) and impose bd condition -> 4 knots - 1 (3 bases)
  // k (> 1) knots : likewise, k + 2 (boundary) - 1 (natural cond) = k + 1 bases                                 
  // knots are placed excluding the boundary region, i.e. in (0.05, 0.95)
  // isknot is irrelevant here, as different knot sets require different basis expansion,
  // as opposed to CRAD where knot in/exclusion = column in/exclusion in design matrix
  
  arma::uvec numknots(p);
  arma::uvec numbasis(p); // number of columns associated with p's smooth covariate
  arma::uvec knotalive(p, arma::fill::ones);
  
  if(NS){
    for(unsigned i{0}; i < p; i++){
      numknots(i) = arma::accu(knotsidx == (i+1));
      if(numknots(i) == 0){
        numbasis(i) = 1; // u
        knotalive(i) = 0;
      } else {
        numknots(i) += 2; // add boundary knots
        numbasis(i) = numknots(i) - 1;
      }
    }
    
  } else {
    for(unsigned i{0}; i < p; i++){
      numknots(i) = arma::accu(knotsidx == (i+1));
      if(numknots(i) == 0){
        numbasis(i) = 3; // u, u^2, u^3
        knotalive(i) = 0;
      } else {
        numbasis(i) = numknots(i) + 3;
      }
    }
  }
  
  arma::uvec betaidx(arma::accu(numbasis) + p_lin, arma::fill::zeros);
  unsigned start{p_lin}, end{p_lin};
  for(unsigned i{0}; i < p; i++){
    end = start + numbasis(i) - 1;
    betaidx.subvec(start, end) += (i+1);
    start = end + 1;
  }
  
  // expand basis column by column
  arma::mat B(N, arma::accu(numbasis));
  if(NS){
    for(unsigned i{0}; i < p; i++){
      B.cols(arma::find(betaidx == (i+1)) - p_lin) = 
        CRADNS_1d_cpp(X.col(i), 
                      knots.elem(arma::find(knotsidx == (i+1))),
                      knotalive(i),
                      bdmargin);
    }
  } else {
    for(unsigned i{0}; i < p; i++){
      B.cols(arma::find(betaidx == (i+1)) - p_lin) = 
        CRAD_1d_cpp(X.col(i), 
                    knots.elem(arma::find(knotsidx == (i+1))),
                    knotalive(i));
    }
  }
  
  // get column-wise mean and center
  arma::vec Bcolmean = arma::mean(B, 0).as_col();
  B.each_row() -= Bcolmean.t();
  B = arma::join_horiz(X_lin, B);
  
  return Rcpp::List::create(
    Rcpp::Named("X") = B,
    Rcpp::Named("NS") = NS,
    Rcpp::Named("means") = Bcolmean,
    Rcpp::Named("knots") = knots,
    Rcpp::Named("knotsidx") = knotsidx,
    Rcpp::Named("betaidx") = betaidx,
    Rcpp::Named("knotalive") = knotalive,
    Rcpp::Named("bdmargin") = bdmargin
  );
}

arma::mat CRAD_test_cpp(const arma::mat &testX, 
                        const arma::mat &X_lin, 
                        const Rcpp::List &CRADlist){
  
  arma::vec knots = Rcpp::as<arma::vec>(CRADlist["knots"]);
  arma::uvec knotsidx = Rcpp::as<arma::uvec>(CRADlist["knotsidx"]);
  arma::uvec betaidx = Rcpp::as<arma::uvec>(CRADlist["betaidx"]);
  arma::vec Bcolmean = Rcpp::as<arma::vec>(CRADlist["means"]);
  arma::uvec knotalive = Rcpp::as<arma::uvec>(CRADlist["knotalive"]);
  bool NS = CRADlist["NS"];
  double bdmargin = (double) CRADlist["bdmargin"];
  
  unsigned p{testX.n_cols};
  unsigned N{testX.n_rows};
  unsigned p_lin{X_lin.n_cols};
  
  // expand basis column by column
  arma::mat B(N, betaidx.n_elem  - p_lin);
  if(NS){
    for(unsigned i{0}; i < p; i++){
      B.cols(arma::find(betaidx == (i+1)) - p_lin) = 
        CRADNS_1d_cpp(testX.col(i), 
                      knots.elem(arma::find(knotsidx == (i+1))),
                      knotalive(i),
                      bdmargin);
    }
  } else {
    for(unsigned i{0}; i < p; i++){
      B.cols(arma::find(betaidx == (i+1)) - p_lin) = 
        CRAD_1d_cpp(testX.col(i), 
                    knots.elem(arma::find(knotsidx == (i+1))),
                    knotalive(i));
    }
  }
  
  // get column-wise mean and center
  B.each_row() -= Bcolmean.t();
  B = arma::join_horiz(arma::mat(N, p_lin, arma::fill::zeros),
                       B);
  
  return B;
}