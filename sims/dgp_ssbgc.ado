
//****************************************************************************/
*! $Id$
*! Data generating processes A-G from Setoguchi et al. (2008) and other articles
//    Setoguchi S, Schneeweiss S, Brookhart MA, Glynn RJ, Cook EF. Evaluating uses of data mining techniques
//       in propensity score estimation: a simulation study. Pharmacoepidemiology and Drug Safety. 2008;17(6):546-555.
//       doi:10.1002/pds.1555
//    Leacy FP, Stuart EA. On the joint use of propensity and prognostic scores in estimation of the average
//       treatment effect on the treated: a simulation study. Statistics in Medicine. 2014;33(20):3488-3508.
//       doi:10.1002/sim.6030
//    Lee BK, Lessler J, Stuart EA. Weight Trimming and Propensity Score Weighting. PLOS ONE. 2011;6(3):e18174.
//       doi:10.1371/journal.pone.0018174
//    Lee BK, Lessler J, Stuart EA. Improving propensity score weighting using machine learning. Statistics in
//        Medicine. 2009;29(3):337-346. doi:10.1002/sim.3782
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.

program define dgp_ssbgc

  syntax anything(name=scenario id="scenario") [, n(int 200)]

  // ---------
  // Draw W_i
  // ---------

  matrix m = ( 0,  0,  0,  0,  0,  0,  0,  0,  0, 0)
  matrix C = ( 1, ///
               0,  1, ///
               0,  0,  1, ///
               0,  0,  0,  1, ///
              .2,  0,  0,  0,  1, ///
               0, .9,  0,  0,  0,  1, ///
               0,  0,  0,  0,  0,  0,  1, ///
               0,  0, .2,  0,  0,  0,  0,  1, ///
               0,  0,  0, .9,  0,  0,  0,  0,  1, ///
               0,  0,  0,  0,  0,  0,  0,  0,  0,  1)
  drawnorm    w1  w2  w3  w4  w5  w6  w7  w8  w9 w10 , n(`n') corr(C) cstorage(lower) clear

  // dichotomize 6 of the variables (will attenuate correlations above)
  // (the other four are continuous)
  foreach v in w1 w3 w5 w6 w8 w9 {
    assert !mi(`v')
    summ `v', mean
    qui replace `v'=(`v'>r(mean))
    recast byte `v'
  }

  // add labels
  label var w1   "w1: Confounder (binary)"
  label var w2   "w2: Confounder (continuous)"
  label var w3   "w3: Confounder (binary)"
  label var w4   "w4: Confounder (continuous)"
  label var w5   "w5: Predictor of treatment assignment only (binary)"
  label var w6   "w6: Predictor of treatment assignment only (binary)"
  label var w7   "w7: Predictor of treatment assignment only (continuous)"
  label var w8   "w8: Predictor of outcome only (binary)"
  label var w9   "w9: Predictor of outcome only (binary)"
  label var w10  "w10: Predictor of outcome only (continuous)"

  // -------------
  // Coefficients
  // -------------

  sca _b0 =  0
  sca _b1 =  0.8
  sca _b2 = -0.25
  sca _b3 =  0.6
  sca _b4 = -0.4
  sca _b5 = -0.8
  sca _b6 = -0.5
  sca _b7 =  0.7
  sca _a0 = -3.85
  sca _a1 =  0.3
  sca _a2 = -0.36
  sca _a3 = -0.73
  sca _a4 = -0.2
  sca _a5 =  0.71
  sca _a6 = -0.19
  sca _a7 =  0.26
  sca _g1 = -0.4 // true effect

  // ---------------------------------------------
  // True propensity score models & outcome model
  // ---------------------------------------------

  // Scenario A (a model with additivity and linearity (Setoguchi et al. 2008)
  if      ("`scenario'"=="A")  gen ps = (1 + exp(-(_b0 + _b1 * w1 + _b2 * w2 + _b3 * w3 + _b4 * w4 + _b5 * w5 + _b6 * w6 + _b7 * w7)))^-1
  // Scenario B (a model with mild non-linearity (Setoguchi et al. 2008)
  else if ("`scenario'"=="B")  gen ps = (1 + exp(-(_b0 + _b1 * w1 + _b2 * w2 + _b3 * w3 + _b4 * w4 + _b5 * w5 + _b6 * w6 + _b7 * w7 + _b2 * w2 * w2)))^-1
  // Scenario C (a model with moderate non-linearity (Setoguchi et al. 2008)
  else if ("`scenario'"=="C")  gen ps = (1 + exp(-(_b0 + _b1 * w1 + _b2 * w2 + _b3 * w3 + _b4 * w4 + _b5 * w5 + _b6 * w6 + _b7 * w7 + _b2 * w2 * w2 + _b4 * w4 * w4 + _b7 * w7 * w7)))^-1
  // Scenario D (a model with mild non-additivity (Setoguchi et al. 2008)
  else if ("`scenario'"=="D")  gen ps = (1 + exp(-(_b0 + _b1 * w1 + _b2 * w2 + _b3 * w3 + _b4 * w4 + _b5 * w5 + _b6 * w6 + _b7 * w7 + _b1 * 0.5 * w1 * w3 + _b2 * 0.7 * w2 * w4 + _b4 * 0.5 * w4 * w5 + _b5 * 0.5 * w5 * w6)))^-1
  // Scenario E (a model with mild non-additivity and non-linearity (Setoguchi et al. 2008)
  else if ("`scenario'"=="E")  gen ps = (1 + exp(-(_b0 + _b1 * w1 + _b2 * w2 + _b3 * w3 + _b4 * w4 + _b5 * w5 + _b6 * w6 + _b7 * w7 + _b2 * w2 * w2 + _b1 * 0.5 * w1 * w3 + _b2 * 0.7 * w2 * w4 + _b4 * 0.5 * w4 * w5 + _b5 * 0.5 * w5 * w6)))^-1
  // Scenario F (a model with moderate non-additivity (Setoguchi et al. 2008)
  else if ("`scenario'"=="F")  gen ps = (1 + exp(-(_b0 + _b1 * w1 + _b2 * w2 + _b3 * w3 + _b4 * w4 + _b5 * w5 + _b6 * w6 + _b7 * w7 + _b1 * 0.5 * w1 * w3 + _b2 * 0.7 * w2 * w4 + _b3 * 0.5 * w3 * w5 + _b4 * 0.7 * w4 * w6 + _b5 * 0.5 * w5 * w7 + _b1 * 0.5 * w1 * w6 + _b2 * 0.7 * w2 * w3 + _b4 * 0.5 * w4 * w5 + _b5 * 0.5 * w5 * w6)))^-1
  // Scenario F (a model with moderate non-additivity (Lee et all. (2011, appendix)
  else if ("`scenario'"=="F2") gen ps = (1 + exp(-(_b0 + _b1 * w1 + _b2 * w2 + _b3 * w3 + _b4 * w4 + _b5 * w5 + _b6 * w6 + _b7 * w7 + _b1 * 0.5 * w1 * w3 + _b2 * 0.7 * w2 * w4 + _b3 * 0.5 * w3 * w5 + _b4 * 0.7 * w4 * w6 + _b5 * 0.5 * w5 * w7 + _b1 * 0.5 * w1 * w6 + _b2 * 0.7 * w2 * w3 + _b3 * 0.5 * w3 * w4 + _b4 * 0.5 * w4 * w5 + _b5 * 0.5 * w5 * w6)))^-1
  // Scenario G (a model with moderate non-additivity and non-linearity (Setoguchi et al. 2008)
  else if ("`scenario'"=="G")  gen ps = (1 + exp(-(_b0 + _b1 * w1 + _b2 * w2 + _b3 * w3 + _b4 * w4 + _b5 * w5 + _b6 * w6 + _b7 * w7 + _b2 * w2 * w2 + _b4 * w4 * w4 + _b7 * w7 * w7 + _b1 * 0.5 * w1 * w3 + _b2 * 0.7 * w2 * w4 + _b3 * 0.5 * w3 * w5 + _b4 * 0.7 * w4 * w6 + _b5 * 0.5 * w5 * w7 + _b1 * 0.5 * w1 * w6 + _b2 * 0.7 * w2 * w3 + _b3 * 0.5 * w3 * w4 + _b4 * 0.5 * w4 * w5 + _b5 * 0.5 * w5 * w6)))^-1

  else {
    di as err "Invalid scenario"
    error 198
  }
  label var ps "ps: True propensity score"

  gen byte a = runiform() < ps
  label var a "a: Treatment assignment"

  gen y = _a0 + _a1 * w1 + _a2 * w2 + _a3 * w3 + _a4 * w4 + _a5 * w8 + _a6 * w9 + _a7 * w10 + _g1 * a
  label var y "y: Observed outcome"

  format %6.3fc w1  w2  w3  w4  w5  w6  w7  w8  w9 w10 ps y
  format %12.0f a

  if      ("`scenario'"=="A")  label data "DGP A: Additivity and linearity (Setoguchi et al. 2008)"
  else if ("`scenario'"=="B")  label data "DGP B: Mild non-linearity (Setoguchi et al. 2008)"
  else if ("`scenario'"=="C")  label data "DGP C: Moderate non-linearity (Setoguchi et al. 2008)"
  else if ("`scenario'"=="D")  label data "DGP D: Mild non-additivity (Setoguchi et al. 2008)"
  else if ("`scenario'"=="E")  label data "DGP E: Mild non-additivity and non-linearity (Setoguchi et al. 2008)"
  else if ("`scenario'"=="F")  label data "DGP F: Moderate non-additivity (Setoguchi et al. 2008)"
  else if ("`scenario'"=="F2") label data "DGP F: Moderate non-additivity (Lee et all. 2011)"
  else if ("`scenario'"=="G")  label data "DGP G: Moderate non-additivity and non-linearity (Setoguchi et al. 2008)"

end
