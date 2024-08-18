#include <RcppArmadillo.h>
#include <chrono>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include "gambmsSupportFunctions.h"
#include "gambmsVS.h"
#include "RJMCMC.h"

Rcpp::List gambmsFREE(const arma::vec &y,
                      const double &glmWeight,
                      const arma::mat &X,
                      const arma::mat &X_pr,
                      const arma::mat &XLin,
                      const arma::vec &offset,
                      const arma::uvec maxk,
                      const arma::vec &Lambda,
                      const unsigned& familyLink,
                      const unsigned& gprior,
                      const double &aa, 
                      const double &bb, 
                      const double &ss, 
                      const double &gg,
                      const int &initS,
                      const int &MCMCiter,
                      const int &thin,
                      double bir_p, double dea_p, double nu,
                      const Rcpp::Function &Rglm,
                      const Rcpp::Function &nearPDres,
                      const bool& storeFit,
                      const bool& forceLin,
                      double& linProb,
                      unsigned printiter){
  
  
  Rcpp::Environment base("package:base");
  Rcpp::Function crossproduct = base["crossprod"];
  Rcpp::Function Rbasechol = base["chol"];
  unsigned P{X.n_cols}; // num of predictors
  unsigned N{X.n_rows};
  unsigned PLin{XLin.n_cols};
  unsigned np{X_pr.n_rows};
  
  // Rcpp::Rcout << "done 1\n";
  
  // choose initial knot via VSknot
  arma::vec EtaHatCurr = etastart(glmWeight, y, familyLink);
  arma::vec EtaHatProp = EtaHatCurr;
  double lpyCurr, lpyProp;
  double comp1Curr, comp1Prop;
  double comp2Curr, comp2Prop;
  double comp3Curr, comp3Prop;
  double r2QmProp, r2QmCurr;
  arma::vec knotsCurr,knotsProp;
  arma::uvec knotsIdxCurr, knotsIdxProp;
  arma::uvec betaIdxCurr, betaIdxProp;
  arma::vec mleCurr, mleProp;
  arma::mat BtrCurr, BprCurr, BtrProp, BprProp;
  arma::mat rootJBetaHatCurr, rootJBetaHatProp;
  
  Rcpp::List the_initial_knot = gambmsVS(y, glmWeight, X, X_pr, XLin, offset, maxk, Lambda,
                                         familyLink, gprior, aa, bb, ss, gg,
                                         initS, Rglm, nearPDres, true, false, 
                                         forceLin, linProb, initS+100);
  
  knotsCurr = Rcpp::as<arma::vec>(the_initial_knot["knotsMax"]);
  knotsIdxCurr = Rcpp::as<arma::uvec>(the_initial_knot["knotsidxMax"]);
  betaIdxCurr = Rcpp::as<arma::uvec>(the_initial_knot["betaidxMax"]);
  BtrCurr = Rcpp::as<arma::mat>(the_initial_knot["B_trMax"]);
  BprCurr = Rcpp::as<arma::mat>(the_initial_knot["B_prMax"]);
  Rcpp::Rcout << "Initial sampling done, numknots is " << knotsCurr.n_elem - (PLin + P) << "\n";
  
  // for(unsigned p{0}; p < P; p++){
  //   if(arma::accu(knotsIdxCurr == (p+1)) == 0){
  //     BtrBpr_update(0,
  //                   X.col(p),
  //                   X_pr.col(p),
  //                   p,
  //                   BtrCurr,
  //                   BprCurr,
  //                   knotsCurr,
  //                   knotsIdxCurr,
  //                   betaIdxCurr,
  //                   0.5,
  //                   0);
  //   }
  // }
  
  MATX_TO_LPY(lpyCurr, 
              comp1Curr, 
              comp2Curr, 
              comp3Curr, 
              r2QmCurr,
              mleCurr, 
              EtaHatCurr, 
              rootJBetaHatCurr,
              y, 
              glmWeight, 
              BtrCurr, 
              betaIdxCurr,
              offset, 
              Lambda, 
              familyLink, 
              gprior, 
              aa, bb, ss, gg, 
              Rglm, false, forceLin, linProb, maxk);
  
  // placeholders for each MCMCiter
  double phi{1.0}, g_{static_cast<double>(N)};
  double v_{0.1},t_{0.1},q_{0.1},l_{0.1},m_{0.1};
  arma::vec PredSmooths, PredLinears, FittedSmooths;
  PredSmooths.set_size(P*np);FittedSmooths.set_size(P*N); PredLinears.set_size(PLin);
  Rcpp::List KNOTS(MCMCiter);
  Rcpp::List KNOTSIDX(MCMCiter);
  arma::mat FITTEDSMOOTHS(MCMCiter, P*N);
  arma::mat PREDSMOOTHS(MCMCiter, P*np);
  arma::mat PREDLINEARS(MCMCiter, PLin);
  arma::vec G(MCMCiter), PHI(MCMCiter);
  arma::vec R2QM(MCMCiter), COMP1(MCMCiter), COMP2(MCMCiter), COMP3(MCMCiter), LPY(MCMCiter);
  
  // miscellaneous for MH accept
  double CtoP, PtoC;
  double num_of_acpt{0.0};
  double o{0.0};
  double num_failed{.0};
  double totaliter{.0};
  std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point toc = tic;
  std::chrono::steady_clock::time_point tic1 = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point toc1 = tic1;
  std::chrono::steady_clock::time_point tic2 = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point toc2 = tic2;
  std::chrono::duration<double, std::milli> TOTAL = toc - tic;
  std::chrono::duration<double, std::milli> TOTAL1 = toc - tic;
  std::chrono::duration<double, std::milli> TOTAL2 = toc - tic;
  double TOTALiter{.0};
  double TOTAL1iter{.0};
  double TOTAL2iter{.0};
  // bool Acpted = false;
  
  arma::vec proposedLPY(MCMCiter*thin*P);
  arma::vec currentLPY(MCMCiter*thin*P);
  arma::vec PropToCurr(MCMCiter*thin*P);
  arma::vec CurrToProp(MCMCiter*thin*P);
  arma::uvec ACPTED(MCMCiter*thin*P);
  arma::vec BtrNcol(MCMCiter);
  
  for(int s{0}; s < MCMCiter; s++){
    tic = std::chrono::steady_clock::now(); ////// <- TIMER
    for(int ss{0}; ss < thin; ss++){
      for(unsigned p{0}; p < P; p++){
        // try{
        //   proposedLPY(static_cast<int>(totaliter)) = lpyProp;
        //   currentLPY(static_cast<int>(totaliter)) = lpyCurr;
        //   PropToCurr(static_cast<int>(totaliter)) = std::log(PtoC);
        //   CurrToProp(static_cast<int>(totaliter)) = std::log(CtoP);
        //   BtrNcol(static_cast<int>(totaliter)) = BtrProp.n_cols;
        //   // Rcpp::Rcout << "assign 1 done \n";
        // } catch(...){}
        totaliter++;
        
        BtrProp = BtrCurr; BprProp = BprCurr; 
        knotsProp = knotsCurr; knotsIdxProp = knotsIdxCurr; betaIdxProp = betaIdxCurr;
        
        try{
          tic1 = std::chrono::steady_clock::now(); ////// <- TIMER
          RJMCMC(PtoC, CtoP, 
                 nu, bir_p, dea_p, 
                 X.col(p), X_pr.col(p), p, maxk(p),
                 BtrProp, BprProp, 
                 knotsProp, knotsIdxProp, betaIdxProp);
          toc1 = std::chrono::steady_clock::now(); ////// <- TIMER
          TOTAL1 += toc1 - tic1;
          TOTAL1iter++;
        } catch(...) {
          toc1 = std::chrono::steady_clock::now(); ////// <- TIMER
          TOTAL1 += toc1 - tic1;
          TOTAL1iter++;
          num_failed++; 
          continue;
        }
        
        try{
          tic2 = std::chrono::steady_clock::now(); ////// <- TIMER
          // Rcpp::Rcout << BtrProp.row(0).t() << "\n";
          MATX_TO_LPY(lpyProp, comp1Prop, comp2Prop, comp3Prop, r2QmProp,
                      mleProp, EtaHatProp, rootJBetaHatProp,
                      y, glmWeight, 
                      BtrProp, betaIdxProp,
                      offset, Lambda, 
                      familyLink, gprior, aa, bb, ss, gg, 
                      Rglm, false, forceLin, linProb, maxk);
          toc2 = std::chrono::steady_clock::now(); ////// <- TIMER
          TOTAL2 += toc2 - tic2;
          TOTAL2iter++;
        } catch(...) {
          toc2 = std::chrono::steady_clock::now(); ////// <- TIMER
          TOTAL2 += toc2 - tic2;
          TOTAL2iter++;
          num_failed++; 
          continue;
        }
        
        if(std::isnan(r2QmProp)) Rcpp::stop("nan occured");
        
        // 2. coin toss to make a jump or not
        o = (lpyProp - lpyCurr) + (std::log(PtoC) - std::log(CtoP));
        if(std::isnan(o) != 0){
          // Rcpp::Rcout << "nan occured \n";
          // Rcpp::Rcout << "lpy_prop : " << lpyProp << '\n';
          // Rcpp::Rcout << "lpy_curr : " << lpyCurr << '\n';
          // Rcpp::Rcout << "log(PtoC) : " << std::log(PtoC) << '\n';
          // Rcpp::Rcout << "log(CtoP) : " << std::log(CtoP) << '\n';
          // num_failed++;
          continue;
        }
        
        // if(knotsProp.n_elem < knotsCurr.n_elem){
        //   Rcpp::Rcout << "death proposed \n";
        //   Rcpp::Rcout << "knotsProp.n_elem " << knotsProp.n_elem-1 << "\n";
        //   Rcpp::Rcout << "knotsCurr.n_elem " << knotsCurr.n_elem-1 << "\n";
        //   Rcpp::Rcout << "lpyProp " << lpyProp <<
        //     " std::log(PtoC) " << std::log(PtoC) << "\n";
        //   Rcpp::Rcout << "lpyCurr " << lpyCurr <<
        //     " std::log(CtoP) " << std::log(CtoP) << "\n";
        // }
        
        // Rcpp::Rcout << "o calculated \n";
        if(arma::randu<double>() < std::exp(o)){ // if accepted
          // if(knotsProp.n_elem == 0){
          //   Rcpp::Rcout << "no knot selected \n";
          //   Rcpp::Rcout << knotsCurr.n_elem << " " << knotsProp.n_elem << "\n";
          // }
          BtrCurr = BtrProp; 
          BprCurr = BprProp; 
          knotsCurr = knotsProp; 
          knotsIdxCurr = knotsIdxProp; 
          betaIdxCurr = betaIdxProp;
          lpyCurr = lpyProp;
          comp1Curr = comp1Prop; 
          comp2Curr = comp2Prop; 
          comp3Curr = comp3Prop;
          r2QmCurr = r2QmProp; 
          mleCurr = mleProp; 
          rootJBetaHatCurr = rootJBetaHatProp;
          ++num_of_acpt;
          // ACPTED(static_cast<int>(totaliter)) = true;
        } else { // if not accepted, do nothing
          // ACPTED(static_cast<int>(totaliter)) = false;
          continue;
        }
        
        
      }
    }
    
    // 3. sample from model_curr
    // Rcpp::Rcout << "knotsC " << knotsC.t() << "\n";
    // Rcpp::Rcout << "BtrCurr.n_cols " << BtrCurr.n_cols << "\n";
    // Rcpp::Rcout << "betaidxC " << betaidxC.t() << "\n";
    // Rcpp::Rcout << "trying MATX_TO_SAMPLE \n";
    MATX_TO_SAMPLE(FittedSmooths, PredSmooths, PredLinears,
                   phi, g_, v_, t_, q_, l_, m_,
                   y, glmWeight, mleCurr, EtaHatCurr, rootJBetaHatCurr, r2QmCurr, 
                   BtrCurr, BprCurr, betaIdxCurr, 
                   P, PLin,
                   familyLink, gprior, 
                   aa, bb, ss, gg,
                   nearPDres, crossproduct, Rbasechol, storeFit);
    KNOTSIDX(s) = knotsIdxCurr;
    KNOTS(s) = knotsCurr;
    LPY(s) = lpyCurr;
    R2QM(s) = r2QmCurr;
    COMP1(s) = comp1Curr;
    COMP2(s) = comp2Curr;
    COMP3(s) = comp3Curr;
    // BtrNcol(s) = BtrCurr.n_cols;
    
    // make prediction, store results
    FITTEDSMOOTHS.row(s) = FittedSmooths.t();
    PREDSMOOTHS.row(s) = PredSmooths.t();
    PREDLINEARS.row(s) = PredLinears.t();
    G(s) = g_; 
    PHI(s) = phi;
    
    // Rcpp::Rcout << "R2 " << r2QmProp << "\n";
    toc = std::chrono::steady_clock::now(); ////// <- TIMER
    TOTAL += toc - tic;
    TOTALiter++;
    
    if(s % printiter == 0 && s > 0){
      Rcpp::Rcout << "Finished " << s << 
        "th iter, the most recent numknots is " << knotsCurr.n_elem - (PLin + P) << 
          ", Accepted " << num_of_acpt  << " / " << totaliter  << ", MLE failed "<< num_failed << '\n'; 
      Rcpp::Rcout << " - Avg time per single update " << std::chrono::duration<double>(TOTAL).count() / (TOTALiter) << "\n";
      // Rcpp::Rcout << "  - Proposal (Rev. Jump) " << P* thin * std::chrono::duration<double>(TOTAL1).count() / TOTAL1iter << "\n";
      // Rcpp::Rcout << "  - Density evaluation (MLE) " << P* thin * std::chrono::duration<double>(TOTAL2).count() / TOTAL2iter << "\n";
      // Rcpp::Rcout << "  - The rest is squandered on the failed attempts to jump " << "\n";
      // Rcpp::Rcout << " - Total num of atomic iterations " << (int) totaliter << "\n";
    }
    
  }
  
  
  Rcpp::List OUT;
  OUT["KNOTS"] = KNOTS;
  OUT["KNOTSIDX"] = KNOTSIDX;
  OUT["lpy"] = LPY;
  OUT["r2Qm"] = R2QM;
  OUT["comp1"] = COMP1;
  OUT["comp2"] = COMP2;
  OUT["comp3"] = COMP3;
  OUT["G"] = G;
  OUT["PHI"] = PHI;
  OUT["FITTEDSMOOTHS"] = FITTEDSMOOTHS;
  OUT["PREDSMOOTHS"] = PREDSMOOTHS;
  OUT["PREDLINEARS"] = PREDLINEARS;
  OUT["AcceptProp"] = num_of_acpt / totaliter;
  // OUT["fuckyou"] = Rcpp::DataFrame::create(
  //   _["proposedLPY"] = proposedLPY,
  //   _["currentLPY"] = currentLPY,
  //   _["PropToCurr"] = PropToCurr,
  //   _["CurrToProp"] = CurrToProp,
  //   _["BtrNcol"] = BtrNcol
  // );
  OUT["maxk"] = maxk;
  return OUT;
}
























