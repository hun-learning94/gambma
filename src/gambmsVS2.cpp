#include <RcppArmadillo.h>
#include <chrono>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include "CRAD2.h"
#include "gambmsSupportFunctions.h"

void ztoknots(const arma::ivec &z,
              const arma::ivec &isknot,
              const arma::vec &knots,
              const arma::uvec &knotsidx,
              arma::vec &knotsP,
              arma::uvec &knotsidxP){
  arma::vec knotsAug = arma::vec(z.n_elem, arma::fill::zeros);
  arma::uvec knotsidxAug = arma::uvec(z.n_elem, arma::fill::zeros);
  knotsAug.elem(arma::find(isknot > 0)) = knots;
  knotsidxAug.elem(arma::find(isknot > 0)) = knotsidx;
  knotsP = knotsAug.elem(arma::find((z>0)&&(isknot>0)));
  knotsidxP = knotsidxAug.elem(arma::find((z>0)&&(isknot>0)));
}

void constructnew(const arma::mat &X,
                  const arma::mat &X_pr,
                  const arma::ivec &isknot,
                  const arma::vec &knots,
                  const arma::uvec &knotsidx,
                  const arma::ivec &z,
                  const arma::mat &B_tr_z,
                  const arma::mat &B_pr_z,
                  const arma::uvec &betaidx_z,
                  const arma::ivec &zP,
                  arma::mat &B_tr_zP,
                  arma::mat &B_pr_zP,
                  arma::uvec &betaidx_zP,
                  int &numflags){

  // get knots_z, knotsidx_z, knots_zP, knotsidx_zP
  arma::vec knots_z, knots_zP;
  arma::uvec knotsidx_z, knotsidx_zP;
  ztoknots(z, isknot, knots, knotsidx, knots_z, knotsidx_z);
  ztoknots(zP, isknot, knots, knotsidx, knots_zP, knotsidx_zP);

  // get p (starts from zero)
  arma::uvec knotsidxAug = arma::uvec(z.n_elem, arma::fill::zeros);
  knotsidxAug.elem(arma::find(isknot > 0)) = knotsidx;
  arma::uvec whichp = knotsidxAug.elem(arma::find(z-zP!=0));
  unsigned p = whichp(0)-1;
  arma::vec xtr = X.col(p);
  arma::vec xpr = X_pr.col(p);

  arma::vec knots_z_sub = knots_z.elem(arma::find(knotsidx_z == (p+1)));
  arma::uvec knotsidx_z_sub = knotsidx_z.elem(arma::find(knotsidx_z == (p+1)));
  arma::vec knots_zP_sub = knots_zP.elem(arma::find(knotsidx_zP == (p+1)));
  arma::uvec knotsidx_zP_sub = knotsidx_zP.elem(arma::find(knotsidx_zP == (p+1)));

  double upperknot{.0}, maxknot{.0};
  bool flag{false}, knotalive{true};
  if(knots_z_sub.n_elem == 0) {
    upperknot = .0;
    flag = true; numflags++;
  } else if(knots_zP_sub.n_elem == 0) {
    maxknot = .0;
    knotalive = false;
    flag = true; numflags++;
  } else {
    upperknot = arma::max(knots_z_sub);
    maxknot = arma::max(knots_zP_sub);
    if(maxknot != upperknot){ flag = true; numflags++; }
  }


  // construct basis
  B_tr_zP = B_tr_z.cols(arma::find(betaidx_z != (p+1)));
  B_pr_zP = B_pr_z.cols(arma::find(betaidx_z != (p+1)));
  betaidx_zP = betaidx_z.rows(arma::find(betaidx_z != (p+1)));
  arma::mat B_tr_zP_sub, B_pr_zP_sub;
  arma::uvec betaidx_zP_sub;
  double newknot;
  arma::uvec newknotidx;

  // if(flag){
  if(flag){
    // Rcpp::Rcout << "construct new basis for pth covariate \n";
    B_tr_zP_sub = CRADNS_1d_cpp2(xtr, knots_zP_sub, knotalive, 0);
    arma::vec Bcolmean = arma::mean(B_tr_zP_sub, 0).as_col();
    B_tr_zP_sub.each_row() -= Bcolmean.t();
    B_pr_zP_sub = CRADNS_1d_cpp2(xpr, knots_zP_sub, knotalive,0);
    B_pr_zP_sub.each_row() -= Bcolmean.t();
    betaidx_zP_sub = arma::uvec(B_tr_zP_sub.n_cols, arma::fill::value(p+1));
    // Rcpp::Rcout << "construct new basis for pth covariate DONE \n";
  } else {
    B_tr_zP_sub = B_tr_z.cols(arma::find(betaidx_z == (p+1)));
    B_pr_zP_sub = B_pr_z.cols(arma::find(betaidx_z == (p+1)));
    betaidx_zP_sub = betaidx_z.rows(arma::find(betaidx_z == (p+1)));

    if(knots_z_sub.n_elem < knots_zP_sub.n_elem){

      // Rcpp::Rcout << "add one knot to pth covariate \n";
      // identify which is newknot among knots_zP_sub
      newknotidx = arma::find(arma::abs(knots_zP_sub - arma::join_vert(knots_z_sub, arma::vec(1, arma::fill::zeros))) > .0);
      newknot = knots_zP_sub(newknotidx(0));
      // double bdknot_l = -.01;
      double bdknot_r = 1.01;

      arma::vec Btr_newcol =
        (-armapmax2(arma::pow(arma::abs(xtr - bdknot_r), 3), 0) + armapmax2(arma::pow(arma::abs(xtr - newknot), 3), 0)) / (bdknot_r - newknot) -
        (armapmax2(arma::pow(arma::abs(xtr - maxknot), 3), 0) - armapmax2(arma::pow(arma::abs(xtr - bdknot_r), 3), 0)) / (bdknot_r - maxknot);
      if(Btr_newcol.has_nan()) Rcpp::stop("nan B_trP");
      double Btr_newcol_mean = arma::as_scalar(arma::mean(Btr_newcol));
      Btr_newcol = Btr_newcol - Btr_newcol_mean;
      arma::vec Bpr_newcol =
        (-armapmax2(arma::pow(arma::abs(xpr - bdknot_r), 3), 0) + armapmax2(arma::pow(arma::abs(xpr - newknot), 3), 0)) / (bdknot_r - newknot) -
        (armapmax2(arma::pow(arma::abs(xpr - maxknot), 3), 0) - armapmax2(arma::pow(arma::abs(xpr - bdknot_r), 3), 0)) / (bdknot_r - maxknot);
      if(Bpr_newcol.has_nan()) Rcpp::stop("nan B_prP");
      Bpr_newcol = Bpr_newcol - Btr_newcol_mean;

      B_tr_zP_sub.insert_cols(newknotidx(0)+1, Btr_newcol);
      B_pr_zP_sub.insert_cols(newknotidx(0)+1, Bpr_newcol);
      betaidx_zP_sub.insert_rows(newknotidx(0)+1, arma::uvec(1, arma::fill::value(p+1)));

      // Rcpp::Rcout << "add one knot to pth covariate DONE \n";
    } else {
      // Rcpp::Rcout << "delete one knot from pth covariate \n";
      // identify which knot to kill among knots_z_sub
      newknotidx = arma::find(arma::abs(knots_z_sub - arma::join_vert(knots_zP_sub, arma::vec(1, arma::fill::zeros))) > .0);
      // Rcpp::Rcout << newknotidx(0) << "\n";
      B_tr_zP_sub.shed_col(newknotidx(0)+1);
      B_pr_zP_sub.shed_col(newknotidx(0)+1);
      betaidx_zP_sub.shed_row(newknotidx(0)+1);
      // Rcpp::Rcout << "delete one knot from pth covariate DONE \n";

    }
  }
  B_tr_zP = arma::join_horiz(B_tr_zP, B_tr_zP_sub);
  B_pr_zP = arma::join_horiz(B_pr_zP, B_pr_zP_sub);
  betaidx_zP = arma::join_vert(betaidx_zP, betaidx_zP_sub);

  // Rcpp::Rcout << "Sorting new col \n";
  // arma::uvec sortidx;
  // sortidx = arma::sort_index(betaidx_zP);
  // betaidx_zP = betaidx_zP.elem(sortidx);
  // B_tr_zP = B_tr_zP.cols(sortidx);
  // B_pr_zP = B_pr_zP.cols(sortidx);
}



Rcpp::List gambmsVS2(const arma::vec &y,
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
                    const unsigned &MCMCiter,
                    const Rcpp::Function &Rglm,
                    const Rcpp::Function &nearPDres,
                    bool getmeMAP,
                    const bool& storeFit,
                    unsigned printiter){

  Rcpp::Environment base("package:base");
  Rcpp::Function crossproduct = base["crossprod"];
  Rcpp::Function Rbasechol = base["chol"];

  unsigned P{X.n_cols}; // num of predictors
  unsigned N{X.n_rows};
  unsigned PLin{XLin.n_cols};
  unsigned np{X_pr.n_rows};

  // do grid expansions for grid of knots
  arma::vec knots = arma::vec(arma::accu(maxk), arma::fill::zeros);
  arma::uvec knotsidx = arma::uvec(knots.n_elem, arma::fill::zeros); // knotsidx = (1,1,...,1, 2, ...2, ...), same dim as knots
  arma::vec qts; double bdmargin{.0};
  unsigned start{0}, end{0};
  for(unsigned i{0}; i<P; i++){
    end = start + maxk(i) - 1;
    // qts = arma::linspace(bdmargin, 1.0 - bdmargin, maxk(i) + 2); // all predictors rescaled to [0,1]
    // qts = qts.subvec(1, maxk(i)); // n_elem = maxk
    qts = arma::linspace(bdmargin, 1.0 - bdmargin, maxk(i)); // all predictors rescaled to [0,1]
    knots.subvec(start, end) = arma::quantile(arma::unique(X.col(i)), qts);
    knotsidx.subvec(start, end) += (i + 1);
    start = end + 1;
  }

  // do CRAD basis expansion, define isknot by hand
  bool NS{true};
  int howmanybasisterms = (NS) ? 1 : 3;
  Rcpp::List CRAD_tr = CRAD_cpp2(X, XLin, knots, knotsidx, NS, bdmargin); // full model
  arma::mat B_tr = Rcpp::as<arma::mat>(CRAD_tr["X"]); // expanded for train
  arma::mat B_pr = CRAD_test_cpp2(X_pr, XLin, CRAD_tr);
  arma::uvec betaidx = Rcpp::as<arma::uvec>(CRAD_tr["betaidx"]);
  arma::ivec isknot = markIsKont(betaidx, howmanybasisterms); // mark non-knot basis terms (linear etc)

  Rcpp::List CRAD_trP = CRAD_tr;
  arma::vec knotsP = knots;
  arma::uvec knotsidxP = knotsidx;
  arma::mat B_trP = B_tr;
  arma::mat B_prP = B_pr;
  arma::uvec betaidxP = betaidx;
  int numflags{0};

  // placeholders for each MCMCiter
  double phi{1.0}, g_{static_cast<double>(N)};
  double v_{0.1},t_{0.1},q_{0.1},l_{0.1},m_{0.1};
  arma::vec FittedSmooths, PredSmooths, PredLinears;
  FittedSmooths.set_size(P*N); PredSmooths.set_size(P*np); PredLinears.set_size(PLin);
  arma::mat FITTEDSMOOTHS(MCMCiter, P*N);
  arma::mat PREDSMOOTHS(MCMCiter, P*np);
  arma::mat PREDLINEARS(MCMCiter, PLin);
  arma::imat Z(isknot.n_elem, MCMCiter, arma::fill::zeros);
  arma::vec G(MCMCiter), PHI(MCMCiter);
  arma::vec R2QM(MCMCiter), COMP1(MCMCiter), COMP2(MCMCiter), COMP3(MCMCiter), LPY(MCMCiter);

  // choose initial set of knots
  arma::ivec z(isknot.n_elem, arma::fill::ones);
  arma::ivec zP = z;
  arma::vec EtaHatCurr = etastart(glmWeight, y, familyLink);
  arma::vec EtaHatProp = EtaHatCurr;
  double lpyCurr, lpyProp, comp1Curr, comp1Prop, comp2Curr, comp2Prop, comp3Curr, comp3Prop, r2QmProp, r2QmCurr;
  arma::vec knotsCurr,knotsProp, mleCurr, mleProp;
  arma::mat rootJBetaHatCurr, rootJBetaHatProp;
  bool initialdrawn = false;
  int initialtry = 0;

  arma::ivec isknot_sub = isknot.elem(arma::find(zP > 0));
  arma::ivec z_sub(isknot_sub.n_elem, arma::fill::ones);
  arma::ivec zP_sub(isknot_sub.n_elem, arma::fill::ones);


  // Rcpp::Rcout << "getting initial prop model\n";
  while(!initialdrawn){
    z.elem(arma::find(isknot > 0)) = Rcpp::as<arma::ivec>(Rcpp::rbinom(arma::accu(isknot), 1, 0.25));
    zP = z;
    try{
      ztoknots(zP, isknot, knots, knotsidx, knotsP, knotsidxP);
      CRAD_trP = CRAD_cpp2(X, XLin, knotsP, knotsidxP, NS, bdmargin); // full model
      B_trP = Rcpp::as<arma::mat>(CRAD_trP["X"]); // expanded for train
      B_prP = CRAD_test_cpp2(X_pr, XLin, CRAD_trP);
      betaidxP = Rcpp::as<arma::uvec>(CRAD_trP["betaidx"]);

      MATX_TO_LPY(lpyCurr, comp1Curr, comp2Curr, comp3Curr, r2QmCurr,
                  mleCurr, EtaHatCurr, rootJBetaHatCurr,
                  y, glmWeight,
                  B_trP,
                  betaidxP,
                  offset, Lambda,
                  familyLink, gprior, aa, bb, ss, gg,
                  Rglm, true, maxk);
      initialdrawn = true; break;
    } catch(...) {
      initialtry++;
      Rcpp::Rcout << "Initial failed for " << initialtry << " times \n";
      // Rcpp::stop("");
      if(initialtry > 1000) break;
    }
  }
  if(initialdrawn == false) Rcpp::stop("");
  // Rcpp::Rcout << "got initial prop model\n";
  B_tr = B_trP;
  B_pr = B_prP;
  betaidx = betaidxP;
  z = zP;

  // store MAP knot configuration
  double lpyMax = lpyCurr;
  arma::ivec zMax = z;

  // get yo clock ready
  std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point toc = tic;
  std::chrono::steady_clock::time_point tic1 = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point toc1 = tic1;
  std::chrono::steady_clock::time_point tic2 = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point toc2 = tic2;
  std::chrono::duration<double, std::milli> TOTAL = toc - tic;
  std::chrono::duration<double, std::milli> TOTAL1 = toc - tic;
  std::chrono::duration<double, std::milli> TOTAL2 = toc - tic;

  // iterate through MCMCiter
  // sub iterate through P
  arma::uvec sPermu; unsigned sItr{0};
  // sub sub iterate through knots in random
  arma::uvec ssPermu; int ssCnt{-1}; unsigned ssHlpr{0}, ssItr{0}, MLEfailed{0}; double o;

  for(unsigned i{0}; i<MCMCiter; i++){
    tic = std::chrono::steady_clock::now(); ////// <- TIMER
    sPermu = arma::randperm(P);
    for(unsigned j{0}; j<P; j++){
      sItr = sPermu(j); // current sPermu(j)th variable
      ssPermu = arma::randperm(maxk(sItr));
      ssHlpr = 0; for(unsigned l{0}; l < sItr; l++){ ssHlpr += (maxk(l) + 1); }
      start = ssHlpr; end = ssHlpr + maxk(sItr) - 1; // set of knots for current variable

      for(unsigned k{0}; k<maxk(sItr); k++){
        ssCnt++;
        ssItr = ssPermu(k); // current knot

        // propose and evaluated
        zP = z;
        zP(PLin + ssHlpr + ssItr + 1) = 1 - zP(PLin + ssHlpr + ssItr + 1);

        tic1 = std::chrono::steady_clock::now(); ////// <- TIMER
        try{
          constructnew(X, X_pr, isknot, knots, knotsidx,
                       z, B_tr, B_pr, betaidx,
                       zP, B_trP, B_prP, betaidxP,
                       numflags);
          MATX_TO_LPY(lpyProp, comp1Prop, comp2Prop, comp3Prop, r2QmProp,
                      mleProp, EtaHatProp, rootJBetaHatProp,
                      y, glmWeight,
                      B_trP,
                      betaidxP,
                      offset, Lambda,
                      familyLink, gprior, aa, bb, ss, gg,
                      Rglm, true, maxk);
        } catch(...) {
          MLEfailed++;
          // Rcpp::Rcout << "mle failed \n";
          // Rcpp::stop("\n");
          continue;
        }
        toc1 = std::chrono::steady_clock::now(); ////// <- TIMER
        TOTAL1 += toc1 - tic1;

        // coin toss to make a jump or not
        o = lpyProp - lpyCurr;
        if(std::isnan(o) != 0){
          Rcpp::Rcout << "nan occured \n";
          // Rcpp::Rcout << "betaidxP.t()" << betaidxP.t() <<" \n";
          Rcpp::Rcout << "lpyProp : " << lpyProp << '\n'; Rcpp::Rcout << "lpyCurr : " << lpyCurr << '\n';
          Rcpp::stop("\n");
          continue;
        }
        if(arma::randu<double>() < std::exp(o)){
          // if accepted
          z = zP;
          lpyCurr = lpyProp;
          comp1Curr = comp1Prop; comp2Curr = comp2Prop; comp3Curr = comp3Prop;
          r2QmCurr = r2QmProp; mleCurr = mleProp; rootJBetaHatCurr = rootJBetaHatProp;
          B_tr = B_trP;
          B_pr = B_prP;
          betaidx = betaidxP;

          if(getmeMAP){if(lpyCurr > lpyMax){ lpyMax = lpyCurr; zMax = z; } }
        } else { // if not accepted, do nothing
          continue;
        }

      }

    }

    LPY(i) = lpyCurr;
    R2QM(i) = r2QmCurr;
    COMP1(i) = comp1Curr;
    COMP2(i) = comp2Curr;
    COMP3(i) = comp3Curr;
    Z.col(i) = z;

    // sample
    tic2 = std::chrono::steady_clock::now(); ////// <- TIMER
    MATX_TO_SAMPLE(FittedSmooths, PredSmooths, PredLinears,
                   phi, g_, v_, t_, q_, l_, m_,
                   y, glmWeight, mleCurr, EtaHatCurr, rootJBetaHatCurr, r2QmCurr,
                   B_tr,
                   B_pr,
                   betaidx,
                   P, PLin,
                   familyLink, gprior,
                   aa, bb, ss, gg,
                   nearPDres, crossproduct, Rbasechol, storeFit);
    toc2 = std::chrono::steady_clock::now(); ////// <- TIMER
    TOTAL2 += toc2 - tic2;

    FITTEDSMOOTHS.row(i) = FittedSmooths.t();
    PREDSMOOTHS.row(i) = PredSmooths.t();
    PREDLINEARS.row(i) = PredLinears.t();
    G(i) = g_;
    PHI(i) = phi;

    toc = std::chrono::steady_clock::now(); ////// <- TIMER
    TOTAL += toc - tic;

    if((i+1) % printiter == 0) {
      Rcpp::Rcout <<  "Gibbs iter " << (i+1) << ", MLE failed " << MLEfailed <<  "\n";
      Rcpp::Rcout << "Avg time per single update " << std::chrono::duration<double>(TOTAL).count() / i << "\n";
      Rcpp::Rcout << " -  MLE and model evidence " << std::chrono::duration<double>(TOTAL1).count() / i << "\n";
      Rcpp::Rcout << " -  Posterior sampling " << std::chrono::duration<double>(TOTAL2).count() / i << "\n";
      Rcpp::Rcout << "new basis "<< numflags << " / " << ssCnt << "\n";
      // Rcpp::Rcout << "Total num of atomic iterations " << (int) ssItr << "\n";
      // Rcpp::Rcout << zP.t() << "\n";
      // Rcpp::Rcout << isknot.t() << "\n";
      // Rcpp::Rcout << betaidx.t() << "\n";
    }
  }

  Rcpp::List OUT;
  OUT["KNOTS"] = Z.t();
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

  if(getmeMAP){
    arma::vec knotsMax;
    arma::uvec knotsidxMax, betaidxMax;
    arma::mat B_trMax, B_prMax;
    zToKnots(knotsMax, knotsidxMax, betaidxMax, B_trMax, B_prMax, zMax,
             betaidx, B_tr, B_pr, isknot, knots, knotsidx);
    OUT["knotsMax"] = knotsMax;
    OUT["knotsidxMax"] = knotsidxMax;
    OUT["B_trMax"] = B_trMax;
    OUT["B_prMax"] = B_prMax;
    OUT["betaidxMax"] = betaidxMax;
  }

  return OUT;
}













































