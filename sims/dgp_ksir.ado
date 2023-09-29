*! $Id$
*! Data generating processe (DGP) senarios from:
//    Kang, J. D. and Schafer, J. L. (2007) Demystifying double robustness: a comparison of
//            alternative strategies for estimating a population mean from incomplete data
//            (with discussions). Statist. Sci., 22, 523–539.
//    Imai K, Ratkovic M. (2014) Covariate balancing propensity score. Journal of the Royal
//            Statistical Society: Series B (Statistical Methodology, 76, 1, 243–263.
//
// Note: The dataset in memomry will be replaced
//
*! By Keith Kranker
// Last updated $Date$

program define dgp_ksir
  version 15.1
  syntax, n(int)         /// Speficy sample size (T+C)
    [ HISTogram          /// Make a histogram of the propensity scores (requires psmatch2)
      irversion          /// use Imai & Ratkovic's formula for X4 instead of of Kang & Schafer's formula
    ]

  // true covariates
  matrix m = ( 0,  0,  0,  0)
  matrix C = ( 1, ///
               0,  1, ///
               0,  0,  1, ///
               0,  0,  0,  1)
  drawnorm z1 z2 z3 z4, n(`n') corr(C) cstorage(lower) clear

  // observed covariates
  gen x1 = exp(z1 / 2),
  gen x2 = z2 / (1 + exp(z1)) + 10
  gen x3 = (z1 * z3 / 25 + 0.6)^3
  if ("`irversion'"=="irversion") {
    gen x4 = (z1 + z4 + 20)^2   // I&R article has  x4 = (z1 + z4 + 20)^2
  }
  else {
    gen x4 = (z2 + z4 + 20)^2   // K&S article has  x4 = (z2 + z4 + 20)^2
  }

  // add labels
  label var z1   "z1: Confounder (not seen by analyst)"
  label var z2   "z2: Confounder (not seen by analyst)"
  label var z3   "z3: Confounder (not seen by analyst)"
  label var z4   "z4: Confounder (not seen by analyst)"
  label var x1   "x1: Confounder (seen by analyst)"
  label var x2   "x2: Confounder (seen by analyst)"
  label var x3   "x3: Confounder (seen by analyst)"
  label var x4   "x4: Confounder (seen by analyst)"

  // True propensity score models & outcome model
  gen e  = rnormal(0, 1)
  gen y  = 210 + 27.4 * z1 + 13.7 * z2 + 13.7 * z3 + 13.7 * z4 + e
  gen ps = invlogit(-z1 + 0.5 * z2 - 0.25 * z3 - 0.1 * z4)
  gen byte a = runiform() < ps
  if ("`irversion'"=="irversion") {
    summ y if a
    replace y = y - r(mean) + 210
  }
  label var e  "e: Error term (not seen by analyst), distributed normal(0,1)"
  label var ps "ps: True propensity score"
  label var a "a: Treatment assignment"
  label var y "y: Observed outcome"
  note   : In this scenario, the true impact of treatment on the outcome is zero. 
  note   : Comparing the (unweighted) difference in mean outcomes between the treatment and comparison groups will yield an estimated impact of about -20, given the selection bias.

  format %5.3fc z1-z4 x1-x4 ps y e
  cap lab drop dgp_ksir_tc
  label define dgp_ksir_tc 0 "Comparison" 1 "Treatment"
  label val a dgp_ksir_tc

  // optionally, make a histogram
  if ("`histogram'"!="") {
    psgraph, bin(30) treated(a) pscore(ps) name(ps, replace)
  }
end

exit

// to test this program, try:
set seed 1
dgp_ksir, n(5000000)
summ, sep(4)
collapse (mean) z1-z4 x1-x4 ps y, by(a) fast
list

dgp_ksir, n(5000) hist
summ, sep(4)
tabstat y, by(a)

logit a x1-x4
predict naive_ps, pr
gen logit_ps       = logit(ps)
gen logit_naive_ps = logit(naive_ps)
graph box logit_*_ps, over(a) nooutsides name(ps_box) ylab(-6(2)2)
