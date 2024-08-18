#ifndef __RJMCMC__
#define __RJMCMC__
#include <RcppArmadillo.h>

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
                   const unsigned &oldknotidx);

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
            arma::uvec &betaidxP);

#endif // __RJMCMC__