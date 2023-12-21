#include <RcppArmadillo.h>
#include <chrono>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include "gambmsFREE.h"
#include "gambmsFREE2.h"
#include "gambmsVS.h"
#include "gambmsVS2.h"
#include "gambmsEVEN.h"

////////////////////////////////////////////////////////////////////////////////////
// EVEN OUTPUT /////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// OUT["whichmethod"] = whichmethod;
// 1. Direct Enumeration
// OUT["enumerate"] = true;
// OUT["KNOTNUMS"] = KNOTNUMS_MCcandidate.t();
// OUT["lpy"] = LPY;
// OUT["r2Qm"] = R2QM;
// OUT["comp1"] = COMP1;
// OUT["comp2"] = COMP2;
// OUT["comp3"] = COMP3;
// OUT["G"] = G;
// OUT["PHI"] = PHI;
// OUT["FITTEDSMOOTHS"] = FITTEDSMOOTHS;
// OUT["PREDSMOOTHS"] = PREDSMOOTHS;
// OUT["PREDLINEARS"] = PREDLINEARS;
// OUT["KNOTLOCS"] = MCKNOTLOCS;
// OUT["KNOTLOCSIDX"] = MCKNOTLOCSIDX;

// 2. Random Walk Metropolis Hastings 
// OUT["enumerate"] = false;
// OUT["KNOTNUMS"] = KNOTNUMS.t();
// OUT["lpy"] = LPY;
// OUT["r2Qm"] = R2QM;
// OUT["comp1"] = COMP1;
// OUT["comp2"] = COMP2;
// OUT["comp3"] = COMP3;
// OUT["G"] = G;
// OUT["PHI"] = PHI;
// OUT["FITTEDSMOOTHS"] = FITTEDSMOOTHS;
// OUT["PREDSMOOTHS"] = PREDSMOOTHS;
// OUT["PREDLINEARS"] = PREDLINEARS;

////////////////////////////////////////////////////////////////////////////////////
// VS OUTPUT ///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// OUT["KNOTS"] = Z.t();
// OUT["lpy"] = LPY;
// OUT["r2Qm"] = R2QM;
// OUT["comp1"] = COMP1;
// OUT["comp2"] = COMP2;
// OUT["comp3"] = COMP3;
// OUT["G"] = G;
// OUT["PHI"] = PHI;
// OUT["PREDSMOOTHS"] = PREDSMOOTHS;
// OUT["PREDLINEARS"] = PREDLINEARS;

////////////////////////////////////////////////////////////////////////////////////
// FREE OUTPUT /////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// OUT["KNOTS"] = KNOTS;
// OUT["KNOTSIDX"] = KNOTSIDX;
// OUT["lpy"] = LPY;
// OUT["r2Qm"] = R2QM;
// OUT["comp1"] = COMP1;
// OUT["comp2"] = COMP2;
// OUT["comp3"] = COMP3;
// OUT["G"] = G;
// OUT["PHI"] = PHI;
// OUT["FITTEDSMOOTHS"] = FITTEDSMOOTHS;
// OUT["PREDSMOOTHS"] = PREDSMOOTHS;
// OUT["PREDLINEARS"] = PREDLINEARS;
