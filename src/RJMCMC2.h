#ifndef __RJMCMC2__
#define __RJMCMC2__
#include <RcppArmadillo.h>

void BtrBpr_update2(const int &BDRL,
                   const unsigned &p,
                   arma::mat &B_trP,
                   arma::mat &B_prP,
                   arma::vec &knotsP,
                   arma::uvec &knotsidxP,
                   arma::uvec &betaidxP,
                   const double &newknot,
                   const unsigned &oldknotidx,
                   const arma::mat& X,
                   const arma::mat& X_pr,
                   const arma::mat& XLin);

void RJMCMC2(double &PtoC, 
            double &CtoP,
            const double &nu, 
            const double &bir_p, 
            const double &dea_p, 
            const unsigned &p,
            const unsigned &maxk,
            arma::mat &B_trP,
            arma::mat &B_prP,
            arma::vec &knotsP,
            arma::uvec &knotsidxP,
            arma::uvec &betaidxP,
            const arma::mat& X,
            const arma::mat& X_pr,
            const arma::mat& XLin);

#endif // __RJMCMC2__