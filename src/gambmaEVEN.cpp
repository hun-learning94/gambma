#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include "gambmaSupportFunctions.h"

//[[Rcpp::export]]
Rcpp::List gambmaEVEN(const arma::vec &y,
                      const double &glmWeight,
                      const arma::mat &X,
                      const arma::mat &X_pr,
                      const arma::mat &XLin,
                      const arma::vec &offset,
                      const arma::uvec &maxk,
                      const arma::vec &Lambda,
                      const unsigned& familyLink,
                      const unsigned& gprior,
                      const double &aa, 
                      const double &bb, 
                      const double &ss, 
                      const double &gg,
                      const bool& enumerate,
                      unsigned numMCcandidate,
                      const unsigned &MCiter,
                      const unsigned &MCMCiter,
                      const Rcpp::Function &Rglm,
                      const Rcpp::Function &nearPDres,
                      const bool& storeFit,
                      unsigned printiter){
  Rcpp::Environment base("package:base");
  Rcpp::Function crossproduct = base["crossprod"];
  Rcpp::Function Rbasechol = base["chol"];
  
  unsigned P{X.n_cols}; // num of predictors
  unsigned N{X.n_rows};
  unsigned PLin{XLin.n_cols};
  unsigned np{X_pr.n_rows};
  
  // calculate model space dimension 
  arma::uvec knotnums = arma::uvec(P, arma::fill::ones);
  unsigned maxk_max = arma::max(maxk); // base of knotsidx
  int64_t modelspace = arma::prod(maxk+1);
  // Rcpp::Rcout << "dim of model space " << modelspace << "\n";
  
  // whether to enumerate or not
  Rcpp::String whichmethod;
  if(enumerate){
    whichmethod = "Enumerate";
    if(modelspace < numMCcandidate) {numMCcandidate = modelspace;}
    if(modelspace > std::pow(30, 3)){
      Rcpp::Rcout << "too large to enumerate, switching to MH \n";
      whichmethod = "MH";
    }
  } else {whichmethod = "MH";}
  
  Rcpp::List OUT;
  OUT["whichmethod"] = whichmethod;
  
  // 1. Direct Enumeration
  if(whichmethod == "Enumerate"){
    OUT["enumerate"] = true;
    
    arma::umat KNOTNUMS = repPermu_arma(maxk, P); // matrix of all enumerations
    modelspace = KNOTNUMS.n_rows;
    arma::umat KNOTNUMS_MCcandidate(P, numMCcandidate);
    Rcpp::Rcout << "Model space " << modelspace << " direct enumeration \n";
    
    // set up placeholders to store 
    // MLE (Rcpp::List), EtaHat (arma::mat), rootJBetaHat(Rcpp::List), 
    // r2Qm, lpy, comp1, comp2, comp3 (all arma::vec)
    // store only numMCcandidate models with highest lpys
    Rcpp::List MLE(numMCcandidate), ROOTJBETAHAT(numMCcandidate);
    arma::mat ETAHAT(N, numMCcandidate);
    arma::vec R2QM(numMCcandidate), COMP1(numMCcandidate), COMP2(numMCcandidate), COMP3(numMCcandidate);
    double preposterous{-1000000.0};
    arma::vec LPY(numMCcandidate, arma::fill::value(preposterous));
    
    double lpy, minlpy{preposterous}, comp1, comp2, comp3, r2Qm;
    arma::vec EtaHat = etastart(glmWeight, y, familyLink);
    arma::vec mle;
    arma::mat rootJBetaHat;
    
    unsigned whichmin{0};
    arma::vec knots;
    arma::uvec knotsidx;
    arma::mat B_tr, B_pr;
    
    // 1) search all model space to find highest MODELS
    for(unsigned i{0}; i < modelspace; i++){
      
      if(i % printiter == 0 && i > 0){Rcpp::Rcout << "Enumerating ... " << i << "\n";} 
      try{
        KNOT_TO_LPY(KNOTNUMS.row(i).t(), knots, knotsidx, 
                    lpy, comp1, comp2, comp3, r2Qm, mle, EtaHat, rootJBetaHat,
                    y, glmWeight, X, XLin, offset, Lambda, 
                    familyLink, gprior, aa, bb, ss, gg,
                    Rglm, false, maxk);
        // Rcpp::Rcout << comp2 << "\n";
      } catch(...) {
        continue;
      }
      // Rcpp::Rcout << "Enumerating ... " << i << "\n";
      // filter the highest numMCcandidate ones
      if(i < numMCcandidate){
        KNOTNUMS_MCcandidate.col(i) = KNOTNUMS.row(i).t();
        MLE(i) = mle; ROOTJBETAHAT(i) = rootJBetaHat;
        ETAHAT.col(i) = EtaHat;
        R2QM(i) = r2Qm; COMP1(i) = comp1; COMP2(i) = comp2; COMP3(i) = comp3;
        LPY(i) = lpy; 
        // Rcpp::Rcout << comp2 << "\n";
      } else {
        if(i == numMCcandidate){ minlpy = arma::min(LPY); whichmin = arma::index_min(LPY); }
        if(lpy > minlpy){
          KNOTNUMS_MCcandidate.col(whichmin) = KNOTNUMS.row(i).t();
          MLE(whichmin) = mle; ROOTJBETAHAT(whichmin) = rootJBetaHat;
          ETAHAT.col(whichmin) = EtaHat;
          R2QM(whichmin) = r2Qm; COMP1(whichmin) = comp1; COMP2(whichmin) = comp2; COMP2(whichmin) = comp3;
          LPY(whichmin) = lpy; 
          minlpy = arma::min(LPY); whichmin = arma::index_min(LPY);
        }
      }
    }
    Rcpp::Rcout << "Enumerate done \n";
    
    OUT["KNOTNUMS"] = KNOTNUMS_MCcandidate.t();
    OUT["lpy"] = LPY;
    OUT["r2Qm"] = R2QM;
    OUT["comp1"] = COMP1;
    // Rcpp::Rcout << COMP2.t() << "\n";
    OUT["comp2"] = COMP2;
    OUT["comp3"] = COMP3;
    LPY = LPY - arma::max(LPY);
    arma::vec WEIGHT = arma::exp(LPY); // normalize weights
    WEIGHT = WEIGHT / accu(WEIGHT);
    
    // 2) MC sampling of models
    arma::mat FITTEDSMOOTHS(MCiter, P*N);
    arma::mat PREDSMOOTHS(MCiter, P*np);
    arma::mat PREDLINEARS(MCiter, PLin);
    Rcpp::List MCKNOTLOCS(MCiter);
    Rcpp::List MCKNOTLOCSIDX(MCiter);
    arma::vec G(MCiter), PHI(MCiter);
    
    Rcpp::Rcout << "getting MC_posterior \n";
    MC_POSTERIOR(FITTEDSMOOTHS, PREDSMOOTHS, PREDLINEARS, MCKNOTLOCS, MCKNOTLOCSIDX, PHI, G,
                 numMCcandidate, MCiter, WEIGHT, 
                 y, glmWeight, X, X_pr, XLin, KNOTNUMS_MCcandidate,
                 MLE, ETAHAT, ROOTJBETAHAT, R2QM, 
                 familyLink, gprior,
                 aa, bb, ss, gg, 
                 nearPDres, crossproduct, Rbasechol, storeFit);
    
    OUT["G"] = G;
    OUT["PHI"] = PHI;
    OUT["FITTEDSMOOTHS"] = FITTEDSMOOTHS;
    OUT["PREDSMOOTHS"] = PREDSMOOTHS;
    OUT["PREDLINEARS"] = PREDLINEARS;
    OUT["KNOTLOCS"] = MCKNOTLOCS;
    OUT["KNOTLOCSIDX"] = MCKNOTLOCSIDX;
  }
  
  // 2. Random Walk Metropolis Hastings 
  if(whichmethod == "MH"){
    OUT["enumerate"] = false;
    Rcpp::Rcout << "Too large modelspace, trying MH ...\n";
    
    bool store = true;
    if(P > 6){
      store = false;
      Rcpp::Rcout << "model space too large, do not store \n";
    }
    
    // placeholders for each MCMCiter
    double phi{1.0}, g_{static_cast<double>(N)};
    double v_{0.1},t_{0.1},q_{0.1},l_{0.1},m_{0.1};
    arma::vec FittedSmooths, PredSmooths, PredLinears;
    FittedSmooths.set_size(P*N); PredSmooths.set_size(P*np); PredLinears.set_size(PLin); 
    arma::mat FITTEDSMOOTHS(MCMCiter, P*N);
    arma::mat PREDSMOOTHS(MCMCiter, P*np);
    arma::mat PREDLINEARS(MCMCiter, PLin);
    arma::umat KNOTNUMS(P, MCMCiter);
    arma::vec G(MCMCiter), PHI(MCMCiter);
    arma::vec R2QM(MCMCiter), COMP1(MCMCiter), COMP2(MCMCiter), COMP3(MCMCiter);
    double preposterous{-1000000.0};
    arma::vec LPY(MCMCiter, arma::fill::value(preposterous));
    
    // placeholders for knots to be saved
    double totaliter{(double)MCMCiter * (double)P};
    double numaccepted{0.0};
    unsigned numrecycled{0}, numstored{0};
    std::vector<int64_t> str_knotnums_idx(totaliter);
    for(unsigned i{0}; i < (totaliter); i++){str_knotnums_idx[i] = -1;}
    arma::vec str_lpy(totaliter), str_r2Qm(totaliter);
    arma::mat str_EtaHat(N, totaliter);
    Rcpp::List str_mle(totaliter);
    Rcpp::List str_rootJBetaHat(totaliter);
    
    // miscellaneous variables for saving
    bool foundya{false};
    int64_t knotnums_idx{0};
    int str_idx{0};
    int propknotnum{0}; 
    int knotnum_lb{0};
    arma::vec knotnum_ub = arma::conv_to<arma::vec>::from(maxk);
    bool cointoss{true};
    int jumpsize{3};
    double maxlpy;
    arma::ivec maxknotnums;
    double o{0.0};
    
    // draw initial model
    arma::ivec knotnumsCurr(P, arma::fill::value(1));
    arma::ivec knotnumsProp = knotnumsCurr;
    arma::uvec knotsidx; arma::vec knots;
    double lpyCurr, lpyProp, comp1Curr, comp1Prop, comp2Curr, comp2Prop, comp3Curr, comp3Prop, r2QmProp, r2QmCurr;
    arma::vec knotsCurr,knotsProp, mleCurr,mleProp;
    arma::vec EtaHatCurr = etastart(glmWeight, y, familyLink);
    arma::vec EtaHatProp = EtaHatCurr;
    arma::mat rootJBetaHatCurr, rootJBetaHatProp;
    try{
      KNOT_TO_LPY(arma::conv_to<arma::uvec>::from(knotnumsProp), knots, knotsidx, 
                  lpyProp, comp1Prop, comp2Prop, comp3Prop, r2QmProp, mleProp, EtaHatProp, rootJBetaHatProp,
                  y, glmWeight, X, XLin, offset, Lambda, 
                  familyLink, gprior, aa, bb, ss, gg,
                  Rglm, false, maxk);
    } catch(...) {
      Rcpp::stop("initial model with a single knot has failed ... duh?");
    }
    if(true){ // initial model always accepted
      knotnumsCurr = knotnumsProp;
      lpyCurr = lpyProp;
      mleCurr = mleProp;
      rootJBetaHatCurr = rootJBetaHatProp;
      EtaHatCurr = EtaHatProp;
      r2QmCurr = r2QmProp;
      comp1Curr = comp1Prop;
      comp2Curr = comp2Prop;
      comp3Curr = comp3Prop;
      // numaccepted++;
    }
    maxlpy= lpyCurr;
    unsigned MLEfailed{0};
    unsigned realtotaliter{0};
    unsigned count{0};
    
    // iterate through P's for 1 iter
    for(unsigned s{0}; s < MCMCiter; s++){
      // iterate through maxk for each P
      // Rcpp::Rcout << "s " << s << "\n";
      for(unsigned p{0}; p < P; p++){ 
        // Rcpp::Rcout << "p " << p << "\n";
        realtotaliter++;
        // propose
        propknotnum = -1.0;
        while((propknotnum < knotnum_lb) || (propknotnum > knotnum_ub(p))){
          cointoss = (arma::randu<double>() > 0.5);
          if(jumpsize > 1){ propknotnum = knotnumsCurr(p) + arma::randi<int>(arma::distr_param(1, jumpsize)) * (2*cointoss - 1);
          } else { propknotnum = knotnumsCurr(p) + jumpsize * (2*cointoss - 1);}
        }
        
        // if within 0, 1, 2, ..., maxk, evaluate the proposed
        knotnumsProp = knotnumsCurr;
        knotnumsProp(p) = propknotnum;
        
        // Rcpp::Rcout << "jump proposed \n";
        
        // 1. get lpyProp
        if(store){
          knotnums_idx = knotnums_to_idx(knotnumsProp, maxk_max);
          // Rcpp::Rcout << knotnums_idx << "\n";
          // Rcpp::Rcout << knotnumsProp << "\n";
          for(unsigned i{0}; i < realtotaliter; i++){
            if(knotnums_idx == str_knotnums_idx[i]){foundya = true; str_idx = i;numrecycled++; break;} else {foundya = false;}
          }
        }
        // Rcpp::Rcout << "finding done \n";
        
        if((foundya == false) || (store == false)){ 
          // first time visited
          try{
            // Rcpp::Rcout << "getting prop model\n";
            KNOT_TO_LPY(arma::conv_to<arma::uvec>::from(knotnumsProp), knotsProp, knotsidx, 
                        lpyProp, comp1Prop, comp2Prop, comp3Prop, r2QmProp, mleProp, EtaHatProp, rootJBetaHatProp,
                        y, glmWeight, X, XLin, offset, Lambda, 
                        familyLink, gprior, aa, bb, ss, gg,
                        Rglm, false, maxk);
          }
          catch(...){
            MLEfailed++;
            continue;
          }
          // if new model, store it
          if(store){
            str_knotnums_idx[count] = knotnums_to_idx(knotnumsProp, maxk_max);
            str_lpy(count) = lpyProp;
            str_mle(count) = mleProp;
            str_rootJBetaHat(count) = rootJBetaHatProp;
            str_EtaHat.col(count) = EtaHatProp;
            str_r2Qm(count) = r2QmProp;
            numstored++;
          }
        } else { 
          // previously visited
          lpyProp = str_lpy(str_idx);
          mleProp = Rcpp::as<arma::vec>(str_mle(str_idx));
          rootJBetaHatProp = Rcpp::as<arma::mat>(str_rootJBetaHat(str_idx));
          EtaHatProp = str_EtaHat.col(str_idx);
          r2QmProp = str_r2Qm(str_idx);
        }
        
        // 2. coin toss to make a jump or not
        o = lpyProp - lpyCurr;
        if(std::isnan(o) != 0){
          Rcpp::Rcout << "nan occured \n";
          Rcpp::Rcout << "lpyProp : " << lpyProp << '\n';
          Rcpp::Rcout << "lpyCurr : " << lpyCurr << '\n';
          Rcpp::stop("\n");
        }
        // Rcpp::Rcout << "o calculated \n";
        if(arma::randu<double>() < std::exp(o)){
          // if accepted
          knotnumsCurr = knotnumsProp;
          lpyCurr = lpyProp;
          mleCurr = mleProp;
          rootJBetaHatCurr = rootJBetaHatProp;
          EtaHatCurr = EtaHatProp;
          r2QmCurr = r2QmProp;
          comp1Curr = comp1Prop;
          comp2Curr = comp2Prop;
          comp3Curr = comp3Prop;
          numaccepted++;
          if(lpyCurr > maxlpy){
            maxlpy = lpyCurr; maxknotnums = knotnumsCurr;
          }
        } else { // if not accepted, do nothing
          continue;
        }
        
        count++;
      }
      
      // Rcpp::Rcout << "now store \n";
      
      LPY(s) = lpyCurr;
      R2QM(s) = r2QmCurr;
      COMP1(s) = comp1Curr;
      COMP2(s) = comp2Curr;
      COMP3(s) = comp3Curr;
      KNOTNUMS.col(s) = arma::conv_to<arma::uvec>::from(knotnumsCurr);
      
      // given theMODEL is selected, sample fit from it
      // Rcpp::Rcout << "getting sample from model \n";
      KNOT_TO_SAMPLE(arma::conv_to<arma::uvec>::from(knotnumsCurr), knots, knotsidx, 
                     FittedSmooths, PredSmooths, PredLinears, phi, g_, v_, t_,q_,l_,m_, 
                     y, glmWeight, X, X_pr, XLin, 
                     mleCurr, EtaHatCurr, rootJBetaHatCurr, r2QmCurr,
                     familyLink, gprior, aa, bb, ss, gg,
                     nearPDres, crossproduct, Rbasechol, storeFit);
      FITTEDSMOOTHS.row(s) = FittedSmooths.t();
      PREDSMOOTHS.row(s) = PredSmooths.t();
      PREDLINEARS.row(s) = PredLinears.t();
      G(s) = g_; PHI(s) = phi;
      
      if((s+1) % printiter == 0) {
        Rcpp::Rcout <<  "MH iter " << (s+1) << " Accept rate " << numaccepted / (realtotaliter) <<
          " (recycled " << numrecycled << "/" << realtotaliter << " = "  << (double)numrecycled / (double)(realtotaliter) << "," << 
            " stored " << numstored <<", MLE failed " << MLEfailed << ") \n";
        // Rcpp::Rcout << "MAP knotnums " << maxknotnums << "\n";
        // Rcpp::Rcout << "most recent knotnums " << knotnumsCurr << "\n";
      }
    }
    OUT["KNOTNUMS"] = KNOTNUMS.t();
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
    OUT["AcceptProp"] = numaccepted / (totaliter);
  }
  
  return OUT;
  
}





