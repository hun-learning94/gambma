#ifndef __gambmsSupportFunctions__
#define __gambmsSupportFunctions__
#include <RcppArmadillo.h>
#include "bayesFactors.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// functions for EVENKNOTS
////////////////////////////////////////////////////////////////////////////////
// repPermu_arma
// a recursion formula to list out all the possible knot combinations
void repPermu_void(const unsigned &maxk, 
                   const unsigned &P, 
                   arma::umat &OUT, 
                   unsigned depth, 
                   unsigned &rowidx);
arma::umat repPermu_arma(const arma::uvec &maxk, 
                         const unsigned &P);


////////////////////////////////////////////////////////////////////////////////
// knotnums_to_idx
// (3, 4, 1, 2) -> number in base maxk
int64_t knotnums_to_idx(const arma::ivec &knotnums, 
                        const unsigned &maxk);

////////////////////////////////////////////////////////////////////////////////
// idx_to_knotnums
// inverse of knotnums_to_idx
arma::ivec idx_to_knotnums(const int64_t &idx, 
                           const unsigned &knotnums_nelem, 
                           const unsigned &maxk);

////////////////////////////////////////////////////////////////////////////////
// knotsandidx
// knotnums to location and index of knots
void knotsandidx(arma::vec &knots,
                 arma::uvec &knotsidx,
                 const arma::uvec &knotnums,
                 const arma::mat &X);

////////////////////////////////////////////////////////////////////////////////
// knotprior 
// componentwise independent Poisson prior for the number of knots
void knotPrior(double &comp3,
               const arma::uvec &betaidx,
               const arma::uvec &maxk,
               const arma::vec &Lambda,
               const bool& isgrid);

////////////////////////////////////////////////////////////////////////////////
// MATX_TO_LPY
// from B_tr the design matrix, yield model evidence, mle and rootJBetaHat 
// input: y, B_tr, offset, glmWeight, 
//        knotnums,
//        Rglm(y, B_tr, offset, glmWeight, EtaHat)
//        Lambda, familyLink, gprior, aa, bb, ss, gg
// output: lpy, comp1, comp2, comp3, r2Qm, mle, rootJBetaHat
// comp1: related to model fit (r2 for Gaussian)
// comp2: related to penalty, differ by gprior
// comp3: related to prior on the number of knots
void MATX_TO_LPY(double &lpy,
                 double &comp1,
                 double &comp2,
                 double &comp3,
                 double &r2Qm,
                 arma::vec &mle,
                 arma::vec &EtaHat,
                 arma::mat &rootJBetaHat,
                 const arma::vec &y,
                 const double &glmWeight,
                 const arma::mat &B_tr,
                 const arma::uvec &betaidx,
                 const arma::vec &offset,
                 const arma::vec &Lambda,
                 const unsigned& familyLink,
                 const unsigned& gprior,
                 const double &aa, 
                 const double &bb, 
                 const double &ss, 
                 const double &gg,
                 const Rcpp::Function &Rglm,
                 const bool& isgrid,
                 const arma::uvec &maxk);

void KNOT_TO_LPY(const arma::uvec &knotnums,
                 arma::vec &knots,
                 arma::uvec &knotsidx,
                 double &lpy,
                 double &comp1,
                 double &comp2,
                 double &comp3,
                 double &r2Qm,
                 arma::vec &mle,
                 arma::vec &EtaHat,
                 arma::mat &rootJBetaHat,
                 const arma::vec &y,
                 const double &glmWeight,
                 const arma::mat &X,
                 const arma::mat &XLin,
                 const arma::vec &offset,
                 const arma::vec &Lambda,
                 const unsigned& familyLink,
                 const unsigned& gprior,
                 const double &aa, 
                 const double &bb, 
                 const double &ss, 
                 const double &gg,
                 const Rcpp::Function &Rglm,
                 const bool& isgrid,
                 const arma::uvec &maxk);


////////////////////////////////////////////////////////////////////////////////
// MATX_to_SAMPLE
// from (mle, EtaHat, rootJBetaHat, r2Qm), produce new posterior sample
// input: y, glmWeight, mle, EtaHat, rootJBetaHat, r2Qm, 
//        B_tr, B_pr, betaidx, 
//        familyLink, gprior, aa, bb, ss, gg
// output: PredSmooths, PredLinears, phi, g
// we do reuse stored mle rootJBetaHat but recreate CRAD basis as
// it is too costly memory-wise to store B_tr for every knot configuration
void MATX_TO_SAMPLE(arma::vec &FittedSmooths,
                    arma::vec &PredSmooths,
                    arma::vec &PredLinears,
                    double &phi,
                    double &g_,
                    double &v_,
                    double &t_,
                    double &q_,
                    double &l_,
                    double &m_,
                    const arma::vec &y,
                    const double &glmWeight,
                    const arma::vec &mle,
                    const arma::vec &EtaHat,
                    const arma::mat &rootJBetaHat,
                    const double &r2Qm,
                    const arma::mat &B_tr,
                    const arma::mat &B_pr,
                    const arma::uvec &betaidx,
                    const double &P,
                    const double &PLin,
                    const unsigned& familyLink,
                    const unsigned& gprior,
                    const double &aa, 
                    const double &bb, 
                    const double &ss, 
                    const double &gg,
                    const Rcpp::Function &nearPDres,
                    const Rcpp::Function &crossproduct,
                    const Rcpp::Function &Rbasechol,
                    const bool &storeFit);

void KNOT_TO_SAMPLE(const arma::uvec &knotnums,
                    arma::vec &knots,
                    arma::uvec &knotsidx,
                    arma::vec &FittedSmooths,
                    arma::vec &PredSmooths,
                    arma::vec &PredLinears,
                    double &phi,
                    double &g_,
                    double &v_,
                    double &t_,
                    double &q_,
                    double &l_,
                    double &m_,
                    const arma::vec &y,
                    const double &glmWeight,
                    const arma::mat &X,
                    const arma::mat &X_pr,
                    const arma::mat &XLin,
                    const arma::vec &mle,
                    const arma::vec &EtaHat,
                    const arma::mat &rootJBetaHat,
                    const double &r2Qm,
                    const unsigned& familyLink,
                    const unsigned& gprior,
                    const double &aa, 
                    const double &bb, 
                    const double &ss, 
                    const double &gg,
                    const Rcpp::Function &nearPDres,
                    const Rcpp::Function &crossproduct,
                    const Rcpp::Function &Rbasechol,
                    const bool& storeFit);

void MC_POSTERIOR(arma::mat &FITTEDSMOOTHS,
                  arma::mat &PREDSMOOTHS,
                  arma::mat &PREDLINEARS,
                  Rcpp::List &MCKNOTLOCS,
                  Rcpp::List &MCKNOTLOCSIDX,
                  arma::vec &PHI,
                  arma::vec &G,
                  const unsigned &numMCcandidate,
                  const unsigned &MCiter,
                  const arma::vec &WEIGHT,
                  const arma::vec &y,
                  const double &glmWeight,
                  const arma::mat &X,
                  const arma::mat &X_pr,
                  const arma::mat &XLin,
                  const arma::umat &KNOTNUMS_MCcandidate,
                  const Rcpp::List &MLE,
                  const arma::mat &ETAHAT,
                  const Rcpp::List &ROOTJBETAHAT,
                  const arma::vec &R2QM,
                  const unsigned& familyLink,
                  const unsigned& gprior,
                  const double &aa, 
                  const double &bb, 
                  const double &ss, 
                  const double &gg,
                  const Rcpp::Function &nearPDres,
                  const Rcpp::Function &crossproduct,
                  const Rcpp::Function &Rbasechol,
                  const bool& storeFit);

////////////////////////////////////////////////////////////////////////////////
// mark non-knot basis terms (linear etc)
arma::ivec markIsKont(const arma::uvec &betaidx, 
                      const unsigned &howmanybasisterms);

////////////////////////////////////////////////////////////////////////////////
// z to knots (VSKNOT, FREEKNOT)
void zToKnots(arma::vec &knotsZ, 
              arma::uvec &knotsidxZ,
              arma::uvec &betaidxZ,
              arma::mat &B_trZ,
              arma::mat &B_prZ,
              const arma::ivec &z,
              const arma::uvec &betaidx,
              const arma::mat &B_tr,
              const arma::mat &B_pr,
              const arma::ivec &isknot,
              const arma::vec &knots,
              const arma::uvec &knotsidx);

#endif // __gambmsSupportFunctions__