#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include "CRAD.h"

void BtrBpr_update(const int &BDRL,
                   const arma::vec &xtr,
                   const arma::vec &xpr,
                   const unsigned &p,
                   arma::mat &B_trP,
                   arma::mat &B_prP,
                   arma::vec &knotsP,
                   arma::uvec &knotsidxP,
                   arma::uvec &betaidxP,
                   const double &newknot,
                   const unsigned &oldknotidx){
  if((BDRL != 0) && (BDRL != 1) && (BDRL != 2)) Rcpp::stop("Invalied BDRL");
  double bdknot_l = .0;
  double bdknot_r = 1.0;
  // B.each_col() = (-arma::pow(arma::abs(x - augknot(nk-1)), 3) + arma::pow(arma::abs(x - augknot(0)), 3)) / 
  //   (augknot(nk-1) - augknot(0));
  // for(unsigned i{0}; i < (nk-2); i++){
  //   B.col(i + 1) -= (arma::pow(arma::abs(x - augknot(nk-1)), 3) - arma::pow(arma::abs(x - augknot(i)), 3)) / 
  //     (augknot(nk-1) - augknot(i));
  // }
  if(BDRL == 1){ // death -> kill existing column
    // Rcpp::Rcout << "oldknotidx " << oldknotidx << ", knotsP.n_elem " << knotsP.n_elem << 
    //   ", knotsidxP.n_elem " << knotsidxP.n_elem << 
    //   ", betaidxP.n_elem " << betaidxP.n_elem << "\n";
    
    knotsP.shed_row(oldknotidx);
    knotsidxP.shed_row(oldknotidx);
    betaidxP.shed_row(oldknotidx);
    B_trP.shed_col(oldknotidx);
    B_prP.shed_col(oldknotidx);
    
  } else { // birth, relocate -> add new column
    if(BDRL == 2){ // relocate -> kill existing column
      // Rcpp::Rcout << "oldknotidx " << oldknotidx << ", knotsP.n_elem " << knotsP.n_elem << 
      //   ", knotsidxP.n_elem " << knotsidxP.n_elem << 
      //   ", betaidxP.n_elem " << betaidxP.n_elem << "\n";
      knotsP.shed_row(oldknotidx);
      knotsidxP.shed_row(oldknotidx);
      betaidxP.shed_row(oldknotidx);
      B_trP.shed_col(oldknotidx);
      B_prP.shed_col(oldknotidx);
    }
    
    // if(knotsP(knotsP.n_elem-1) < 0.0001){
    //   Rcpp::Rcout << BDRL << "\n";
    //   Rcpp::Rcout << knotsP.t() << "\n";
    //   Rcpp::stop("shedding");
    // }
    
    // expand new column
    // Rcpp::Rcout << "Expanding new col \n";
    arma::vec Btr_newcol = 
      (-armapmax(arma::pow(arma::abs(xtr - bdknot_r), 3), 0) + armapmax(arma::pow(arma::abs(xtr - bdknot_l), 3), 0)) / (bdknot_r - bdknot_l) -
      (armapmax(arma::pow(arma::abs(xtr - bdknot_r), 3), 0) - armapmax(arma::pow(arma::abs(xtr - newknot), 3), 0)) / (bdknot_r - newknot);
    if(Btr_newcol.has_nan()) Rcpp::stop("nan B_trP");
    double Btr_newcol_mean = arma::as_scalar(arma::mean(Btr_newcol));
    Btr_newcol = Btr_newcol - Btr_newcol_mean;
    arma::vec Bpr_newcol = 
      (-armapmax(arma::pow(arma::abs(xpr - bdknot_r), 3), 0) + armapmax(arma::pow(arma::abs(xpr - bdknot_l), 3), 0)) / (bdknot_r - bdknot_l) -
      (armapmax(arma::pow(arma::abs(xpr - bdknot_r), 3), 0) - armapmax(arma::pow(arma::abs(xpr - newknot), 3), 0)) / (bdknot_r - newknot);
    if(Bpr_newcol.has_nan()) Rcpp::stop("nan B_prP");
    Bpr_newcol = Bpr_newcol - Btr_newcol_mean;
    
    // append to knotsP, knotsidxP, B_trP, B_prP
    // Rcpp::Rcout << "Appending new col \n";
    knotsP.insert_rows(knotsP.n_elem, arma::vec(1, arma::fill::value(newknot)));
    knotsidxP.insert_rows(knotsidxP.n_elem, arma::uvec(1, arma::fill::value(p+1)));
    betaidxP.insert_rows(betaidxP.n_elem, arma::uvec(1, arma::fill::value(p+1)));
    // Rcpp::Rcout << "B_trP.n_rows " << B_trP.n_rows<< ", Btr_newcol.n_rows " << Btr_newcol.n_rows  << "\n";
    // Rcpp::Rcout << "B_prP.n_rows " << B_prP.n_rows<< ", Bpr_newcol.n_rows " << Bpr_newcol.n_rows  << "\n";
    B_trP.insert_cols(B_trP.n_cols, Btr_newcol);
    B_prP.insert_cols(B_prP.n_cols, Bpr_newcol);
    
    // if(knotsP(knotsP.n_elem-1) < 0.0001){
    //   Rcpp::Rcout << BDRL << "\n";
    //   Rcpp::Rcout << knotsP.t() << "\n";
    //   Rcpp::stop("inserting");
    // }
    
    // sort columns 2) by increasing betaidxP
    // Rcpp::Rcout << "Sorting new col \n";
    // arma::uvec sortidx = sort_index(betaidxP);
    // knotsP = knotsP.elem(sortidx);
    // knotsidxP = knotsidxP.elem(sortidx);
    // betaidxP = betaidxP.elem(sortidx);
    // B_trP = B_trP.cols(sortidx);
    // B_prP = B_prP.cols(sortidx);
    // if(knotsP(knotsP.n_elem-1) < 0.0001){
    //   Rcpp::Rcout << BDRL << "\n";
    //   Rcpp::Rcout << knotsP.t() << "\n";
    //   Rcpp::stop("sorting");
    // }
  }
}


void RJMCMC(double &PtoC, 
            double &CtoP,
            const double &nu, 
            const double &bir_p, 
            const double &dea_p, 
            const arma::vec &xtr,
            const arma::vec &xpr,
            const unsigned &p,
            const unsigned &maxk,
            arma::mat &B_trP,
            arma::mat &B_prP,
            arma::vec &knotsP,
            arma::uvec &knotsidxP,
            arma::uvec &betaidxP)
{
  unsigned idx{0}, oldknotidx{0};
  double dice, idx_knot, new_knot, sum_bir;
  arma::uvec subknotsidx = arma::find((knotsidxP == p+1) && (knotsP > 0));
  arma::vec subknots = knotsP(subknotsidx);
  arma::vec subknotsP;
  int timer{0};
  unsigned BDRL{0};
  
  // throw dice
  dice = arma::randu<double>();
  // if(subknots.n_elem == 1 && (1.0-dice) < dea_p){ // k = 0 is not allowed
  //   while((1.0-dice) < dea_p) dice = arma::randu<double>();
  // }
  // if(subknots.n_elem >= maxk){ // birth not allowed
  //   while(dice < bir_p) dice = arma::randu<double>();
  // }
  // Rcpp::Rcout << dice << "\n";
  // Rcpp::Rcout << bir_p << "\n";
  if(dice < bir_p){
    // birth step ///////////////////////////////////////////////////////////////////////////////////////
    BDRL = 0; 
    // Rcpp::Rcout << "birth chosen \n";
    // 1) choose one knot and sample new knot (sample idx in 0, 1, ..., knots.n_elem - 1)
    try{idx = arma::randi<unsigned>(arma::distr_param(0, subknots.n_elem - 1));}
    catch(...){
      Rcpp::Rcout<< "birth subknots " << subknots.t()<<"\n" ;
      Rcpp::stop("rand i error");
    }
    idx_knot = subknots(idx);
    new_knot = R::rbeta(nu * idx_knot, nu * (1.0 - idx_knot));
    new_knot = std::round(new_knot * 50) / 50;
    if(std::abs(new_knot - 1.0)<0.001) Rcpp::stop("\n");
    if(std::abs(new_knot - 0.0)<0.001) Rcpp::stop("\n");
    
    // 2) get knotsP
    subknotsP = arma::vec(subknots.n_elem + 1, arma::fill::zeros);
    subknotsP.subvec(0, subknots.n_elem - 1) = subknots;
    subknotsP(subknots.n_elem) = new_knot;
    subknotsP = arma::sort(subknotsP);
    
    // Rcpp::Rcout << subknotsP.n_elem << "\n";
    // Rcpp::Rcout << maxk << "\n";
    
    if(subknotsP.n_elem > maxk) Rcpp::stop("\n");
    
    // 3) calculate CtoP, PtoC
    sum_bir = 0.0;
    for(unsigned i{0}; i<subknots.n_elem; i++){
      sum_bir += R::dbeta(new_knot, idx_knot*nu, (1.0 - idx_knot)*nu, false);
    }
    sum_bir = std::min(1.0, sum_bir);
    CtoP = bir_p * sum_bir / (double) subknots.n_elem;
    PtoC = dea_p / (double) (1 + subknots.n_elem);
    // Rcout << "sub done3 \n";
    
    if(std::isnan(CtoP) != 0 || std::isnan(PtoC) != 0 ||
       std::isfinite(log(CtoP)) != 1 || std::isfinite(log(PtoC)) != 1 ){
      // Rcpp::Rcout << "nan/inf occured in birth \n";
      // Rcpp::Rcout << "sum_bir " << sum_bir << '\n';
      // Rcpp::Rcout << "subknots.n_elem " << subknots.n_elem << '\n';
      // Rcpp::Rcout << "CtoP " << R::dbeta(new_knot, idx_knot*nu, (1.0 - idx_knot)*nu, false) << "\n";
      // Rcpp::Rcout << "PtoC " << R::dbeta(idx_knot, new_knot*nu, (1.0 - new_knot)*nu, false) << "\n";
      Rcpp::stop("\n");
    }
  } else if ((1.0 - dice) < dea_p){
    // death step ///////////////////////////////////////////////////////////////////////////////////////
    BDRL = 1; 
    // Rcpp::Rcout << "death chosen \n";
    if(subknots.n_elem == 1) Rcpp::stop("");
    // 1) choose one knot and kill
    try{
      idx = arma::randi<unsigned>(arma::distr_param(0, subknots.n_elem - 1));
      oldknotidx = subknotsidx(idx);
    } catch(...){
      Rcpp::Rcout<< "death subknots " << subknots.t()<<"\n" ;
      Rcpp::stop("rand i error");
    }
    new_knot = subknots(idx);
    // Rcout << "sub done1 \n";
    
    // 2) get knotsP, knotsidxP
    subknotsP = subknots;
    subknotsP.shed_row(idx);
    // Rcout << "sub done2 \n";
    
    // 3) calculate CtoP, PtoC
    sum_bir = 0.0;
    for(unsigned i{0}; i<subknotsP.n_elem; i++){
      sum_bir += R::dbeta(new_knot, subknotsP(i)*nu, (1.0 - subknotsP(i))*nu, false);
    }
    sum_bir = std::min(1.0, sum_bir);
    CtoP = dea_p / (double) subknots.n_elem;
    PtoC = bir_p * sum_bir / (double) (subknots.n_elem - 1);
    // Rcout << "sub done3 \n";
    
    if(std::isnan(CtoP) != 0 || std::isnan(PtoC) != 0 ||
       std::isfinite(log(CtoP)) != 1 || std::isfinite(log(PtoC)) != 1 ){
      // Rcpp::Rcout << "nan/inf occured in death \n";
      // Rcpp::Rcout << "sum_bir " << sum_bir << '\n';
      // Rcpp::Rcout << "subknots.n_elem " << subknots.n_elem << '\n';
      Rcpp::stop("\n");
    }
  } else {
    // relocate step /////////////////////////////////////////////////////////////////////////////////////
    BDRL = 2; // Rcpp::Rcout << "relocation chosen \n";
    // 1) choose one knot and relocate
    try{
      idx = arma::randi<unsigned>(arma::distr_param(0, subknots.n_elem - 1));
      oldknotidx = subknotsidx(idx);
    } catch(...) {
      Rcpp::Rcout<< "relocation subknots " << subknots.t()<<"\n" ;
      Rcpp::stop("rand i error");
    }
    idx_knot = subknots(idx);
    new_knot = R::rbeta(nu * idx_knot, nu * (1.0 - idx_knot));
    new_knot = std::round(new_knot * 50) / 50;
    if(std::abs(new_knot - 1.0)<0.001) Rcpp::stop("\n");
    if(std::abs(new_knot - 0.0)<0.010) Rcpp::stop("\n");
    
    // 2) get knotsP, knotsidxP
    subknotsP = subknots;
    subknotsP(idx) = new_knot; // relocate
    subknotsP = arma::sort(subknotsP);
    // Rcout << "sub done2 \n";
    
    CtoP = (1.0 - dea_p - bir_p) * 
      R::dbeta(new_knot, idx_knot*nu, (1.0 - idx_knot)*nu, false) / (double) subknots.n_elem;
    PtoC = (1.0 - dea_p - bir_p) * 
      R::dbeta(idx_knot, new_knot*nu, (1.0 - new_knot)*nu, false)/ (double) subknots.n_elem;
    
    if(std::isnan(CtoP) != 0 || std::isnan(PtoC) != 0 ||
       std::isfinite(log(CtoP)) != 1 || std::isfinite(log(PtoC)) != 1 ){
      // Rcpp::Rcout << "nan/inf occured in relocation \n";
      // Rcpp::Rcout << "subknots.n_elem " << subknots.n_elem << '\n';
      // Rcpp::Rcout << "new_knot " << new_knot << ", idx_knot " << idx_knot << "\n";
      // Rcpp::Rcout << R::dbeta(new_knot, idx_knot*nu, (1.0 - idx_knot)*nu, false) << "\n";
      // Rcpp::Rcout << R::dbeta(idx_knot, new_knot*nu, (1.0 - new_knot)*nu, false) << "\n";
      Rcpp::stop("\n");
    }
  }
  
  BtrBpr_update(BDRL,
                xtr,
                xpr,
                p,
                B_trP,
                B_prP,
                knotsP,
                knotsidxP,
                betaidxP,
                new_knot,
                oldknotidx);
  // Rcpp::Rcout << "BtrBpr update done \n";
}