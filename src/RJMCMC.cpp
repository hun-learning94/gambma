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
  
  arma::uvec sortidx;
  double bdknot_l = -.01;
  double bdknot_r = 1.01;
  // B.each_col() = (-arma::pow(arma::abs(x - augknot(nk-1)), 3) + arma::pow(arma::abs(x - augknot(0)), 3)) / 
  //   (augknot(nk-1) - augknot(0));
  // for(unsigned i{0}; i < (nk-2); i++){
  //   B.col(i + 1) -= (arma::pow(arma::abs(x - augknot(nk-1)), 3) - arma::pow(arma::abs(x - augknot(i)), 3)) / 
  //     (augknot(nk-1) - augknot(i));
  // }
  if(BDRL == 1){ // death -> kill existing column
    
    knotsP.shed_row(oldknotidx);
    knotsidxP.shed_row(oldknotidx);
    betaidxP.shed_row(oldknotidx);
    B_trP.shed_col(oldknotidx);
    B_prP.shed_col(oldknotidx);
    
  } else { // birth, relocate -> add new column
    if(BDRL == 2){ // relocate -> kill existing column
    
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
    unsigned insrtpos = oldknotidx;
    if(oldknotidx == 0) insrtpos = 1;
    knotsP.insert_rows(insrtpos, arma::vec(1, arma::fill::value(newknot)));
    knotsidxP.insert_rows(insrtpos, arma::uvec(1, arma::fill::value(p+1)));
    betaidxP.insert_rows(insrtpos, arma::uvec(1, arma::fill::value(p+1)));
    B_trP.insert_cols(insrtpos, Btr_newcol);
    B_prP.insert_cols(insrtpos, Bpr_newcol);
    
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
  arma::vec subknotsP, subknots;
  if(subknotsidx.is_empty()){
    subknots.empty();
  } else {
    subknots = knotsP(subknotsidx);
  }
  // int timer{0};
  unsigned BDRL{0};
  
  
  // throw dice
  // birth: BDRL = 0, death: BDRL = 1; relocate: BDRL = 2;
  // note that subknots always include 0 meaning linear basis term
  dice = arma::randu<double>();
  if(dice < bir_p){ BDRL = 0;
  } else if((1.0 - dice) < dea_p){ BDRL = 1;
  } else { BDRL = 2; }
  
  int fuck{0};
  double distance = 50.0;
  
  if(BDRL == 0){
    // birth step ///////////////////////////////////////////////////////////////////////////////////////
    if(subknots.n_elem >= (maxk+1)) Rcpp::stop("");
    if(subknots.n_elem > 0){
      try{idx = arma::randi<unsigned>(arma::distr_param(0, subknots.n_elem - 1));}
      catch(...){
        Rcpp::Rcout<< "birth subknots " << subknots.t()<<"\n" ;
        Rcpp::stop("rand i error");
      }
      idx_knot = subknots(idx);
      oldknotidx = subknotsidx(idx);
      new_knot = R::rbeta(nu * idx_knot, nu * (1.0 - idx_knot));
      new_knot = std::round(new_knot * distance) / distance;
      while(std::abs(new_knot - 1.0)<0.01 || std::abs(new_knot - 0.0)<0.01){
        new_knot = R::rbeta(nu * idx_knot, nu * (1.0 - idx_knot));
        new_knot = std::round(new_knot * distance) / distance;
        fuck++;
        if(fuck > 20) Rcpp::stop("propose failed");
      }
      subknotsP = arma::vec(subknots.n_elem + 1, arma::fill::zeros);
      subknotsP.subvec(0, subknots.n_elem - 1) = subknots;
      subknotsP(subknots.n_elem) = new_knot;
      // subknotsP = arma::sort(subknotsP);
      
      sum_bir = 0.0;
      for(unsigned i{0}; i<subknots.n_elem; i++){
        sum_bir += R::dbeta(new_knot, subknots(i)*nu, (1.0 - subknots(i))*nu, false);
      }
      sum_bir = std::min(1.0, sum_bir);
      CtoP = bir_p * sum_bir / (double) subknots.n_elem;
      PtoC = dea_p / (double) (1 + subknots.n_elem);  
      
    } else {
      // numknot 0 -> 1
      new_knot = 0.001 + arma::randu<double>() * (0.999-0.001);
      new_knot = std::round(new_knot * distance) / distance;
      while(std::abs(new_knot - 1.0)<0.01 || std::abs(new_knot - 0.0)<0.01){
        new_knot = 0.001 + arma::randu<double>() * (0.999-0.001);
        new_knot = std::round(new_knot * distance) / distance;
        fuck++;
        if(fuck > 20) Rcpp::stop("propose failed");
      }
      subknotsP.set_size(1);
      subknotsP(0) = new_knot;
      
      CtoP = bir_p;
      PtoC = dea_p;
    }
    
    if(std::isnan(CtoP) != 0 || std::isnan(PtoC) != 0 ||
       std::isfinite(log(CtoP)) != 1 || std::isfinite(log(PtoC)) != 1 ){
      // Rcpp::Rcout << "nan/inf occured in birth \n";
      // Rcpp::Rcout << "sum_bir " << sum_bir << '\n';
      // Rcpp::Rcout << "subknots.n_elem " << subknots.n_elem << '\n';
      // Rcpp::Rcout << "CtoP " << R::dbeta(new_knot, idx_knot*nu, (1.0 - idx_knot)*nu, false) << "\n";
      // Rcpp::Rcout << "PtoC " << R::dbeta(idx_knot, new_knot*nu, (1.0 - new_knot)*nu, false) << "\n";
      Rcpp::stop("\n");
    }
  } else if (BDRL == 1){
    // death step ///////////////////////////////////////////////////////////////////////////////////////
    // Rcpp::Rcout << "death chosen \n";
    if(subknots.n_elem == 0) Rcpp::stop("");
    // 1) choose one knot and kill
    try{
      idx = arma::randi<unsigned>(arma::distr_param(0, subknots.n_elem - 1));
    } catch(...){
      Rcpp::Rcout<< "death subknots " << subknots.t()<<"\n" ;
      Rcpp::stop("rand i error");
    }
    oldknotidx = subknotsidx(idx);
    new_knot = subknots(idx);
    subknotsP = subknots;
    subknotsP.shed_row(idx);
    
    if(subknots.n_elem > 1){
      sum_bir = 0.0;
      for(unsigned i{0}; i<subknotsP.n_elem; i++){
        sum_bir += R::dbeta(new_knot, subknotsP(i)*nu, (1.0 - subknotsP(i))*nu, false);
      }
      sum_bir = std::min(1.0, sum_bir);
      CtoP = dea_p / (double) subknots.n_elem;
      PtoC = bir_p * sum_bir / (double) (subknots.n_elem - 1);
      
    } else { // numknot 1 -> 0
      // Rcpp::Rcout << "propose killing the last knot \n";
      CtoP = dea_p;
      PtoC = bir_p;
    }
    
    if(std::isnan(CtoP) != 0 || std::isnan(PtoC) != 0 ||
       std::isfinite(log(CtoP)) != 1 || std::isfinite(log(PtoC)) != 1 ){
      // Rcpp::Rcout << "nan/inf occured in death \n";
      // Rcpp::Rcout << "sum_bir " << sum_bir << '\n';
      // Rcpp::Rcout << "subknots.n_elem " << subknots.n_elem << '\n';
      Rcpp::stop("\n");
    }
  } else {
    // relocate step /////////////////////////////////////////////////////////////////////////////////////
    // BDRL = 2; // Rcpp::Rcout << "relocation chosen \n";
    // 1) choose one knot and relocate
    if(subknots.n_elem == 0) Rcpp::stop("");
    try{
      idx = arma::randi<unsigned>(arma::distr_param(0, subknots.n_elem - 1));
    } catch(...) {
      // Rcpp::Rcout<< "relocation subknots " << subknots.t()<<"\n" ;
      Rcpp::stop("rand i error");
    }
    oldknotidx = subknotsidx(idx);
    idx_knot = subknots(idx);
    new_knot = R::rbeta(nu * idx_knot, nu * (1.0 - idx_knot));
    new_knot = std::round(new_knot * distance) / distance;
    if(std::abs(new_knot - 1.0)<0.01) Rcpp::stop("\n");
    if(std::abs(new_knot - 0.0)<0.01) Rcpp::stop("\n");
    
    // 2) get knotsP, knotsidxP
    subknotsP = subknots;
    subknotsP(idx) = new_knot; // relocate
    // subknotsP = arma::sort(subknotsP);
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
  // Rcpp::Rcout << "\n";
  // Rcpp::Rcout << "BDRL " << BDRL << " oldknotidx "<< oldknotidx <<  " new_knot " << new_knot << "\n";
  // Rcpp::Rcout << knotsP.t() << "\n";
  // Rcpp::Rcout << knotsidxP.t() << "\n";
  // Rcpp::Rcout << betaidxP.t() << "\n";
  // Rcpp::Rcout << B_trP.row(0) << "\n";
  // Rcpp::Rcout << B_prP.row(0) << "\n";
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
  // Rcpp::Rcout << knotsP.t() << "\n";
  // Rcpp::Rcout << knotsidxP.t() << "\n";
  // Rcpp::Rcout << betaidxP.t() << "\n";
  // Rcpp::Rcout << B_trP.row(0) << "\n";
  // Rcpp::Rcout << B_prP.row(0) << "\n";
}