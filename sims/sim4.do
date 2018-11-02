if ("`: environment computername'"=="M116" | "`: environment computername'"=="NJ1STATA1") cd "U:\Stata\Ado\Devel\gmatch"
else cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"

clear all
cls

local sim = 4
cap log close sim`sim'
set linesize 180
cap mkdir  sims/sim`sim'
log using "sims/sim`sim'/logfile.log", replace name(sim`sim')

*! $Id$
*! Simulation setup #4
// Refining the Simulation -- larger CV targets, add augmenteed models
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.
mac list _sim

// DGP is saved in a separate .ado file
adopath ++ "./sims"
which dgp_ssbgc
which onerep
which sim_reshape

// options
parallel setclusters 6
local commonopts atet iter(80) cformat(%9.3fc) pformat(%5.3f) sformat(%7.3f) pooledvariance

// control simulations
set seed `sim'  //  1 for simulation 1, 2 for simulation 2, etc.
local reps 1000


// ------------------------------------------------------------------------
// 4-A Variouis estmators with impact = -.09, N=2000
// ------------------------------------------------------------------------

onerep A E G, impact(-0.09) n(2000) ///
  estimators(cbps) cvtargets(99 85 50) `commonopts'

parallel sim, expr(_b) reps(`reps') processors(4): ///
  onerep A E G, impact(-0.09) n(2000) ///
    estimators(cbps) cvtargets(99 85 50) `commonopts'

sim_reshape
save sims/sim`sim'/sim`sim'a.dta, replace

foreach v of var impact_est-wgt_max impact_est_var {
  tabstat `v', by(result) s(N mean p50 sd) nototal
  graph bar (mean) `v', over(estimator) asyvar over(dgp) title(`v') ytitle("") xsize(12) name(`v', replace)
  graph export "sims/sim`sim'/sim`sim'a_`v'.png", replace
}


// ------------------------------------------------------------------------
// 4-B Variouis estmators with impact = -.09, N=2000, with additional noise
// ------------------------------------------------------------------------

onerep A E G, impact(-0.09) n(2000) ///
  estimators(cbps) cvtargets(99 85 50) `commonopts' noise(.5)

parallel sim, expr(_b) reps(`reps') processors(4): ///
  onerep A E G, impact(-0.09) n(2000) ///
    estimators(cbps) cvtargets(99 85 50) `commonopts' noise(.5)

sim_reshape
save sims/sim`sim'/sim`sim'b.dta, replace

foreach v of var impact_est-wgt_max impact_est_var {
  tabstat `v', by(result) s(N mean p50 sd) nototal
  graph bar (mean) `v', over(estimator) asyvar over(dgp) title(`v') ytitle("") xsize(12) name(`v', replace)
  graph export "sims/sim`sim'/sim`sim'b_`v'.png", replace
}


// ------------------------------------------------------------------------------------
// 4-C Variouis estmators with impact = -.09, N=2000, with augmented regression models
// ------------------------------------------------------------------------------------

onerep A E G, impact(-0.09) n(2000) ///
  estimators(cbps) cvtargets(99 85 50) `commonopts' augmented

parallel sim, expr(_b) reps(`reps') processors(4): ///
  onerep A E G, impact(-0.09) n(2000) ///
    estimators(cbps) cvtargets(99 85 50) `commonopts' augmented

sim_reshape
save sims/sim`sim'/sim`sim'c.dta, replace

foreach v of var impact_est-wgt_max impact_est_var {
  tabstat `v', by(result) s(N mean p50 sd) nototal
  graph bar (mean) `v', over(estimator) asyvar over(dgp) title(`v') ytitle("") xsize(12) name(`v', replace)
  graph export "sims/sim`sim'/sim`sim'c_`v'.png", replace
}


// ------------------------------------------------------------------------------------------------------------
// 4-D Variouis estmators with impact = -.09, N=2000, with additional noise and augmented regression models
// ------------------------------------------------------------------------------------------------------------

onerep A E G, impact(-0.09) n(2000) ///
  estimators(cbps) cvtargets(99 85 50) `commonopts' noise(.5) augmented

parallel sim, expr(_b) reps(`reps') processors(4): ///
  onerep A E G, impact(-0.09) n(2000) ///
    estimators(cbps) cvtargets(99 85 50) `commonopts' noise(.5) augmented

sim_reshape
save sims/sim`sim'/sim`sim'd.dta, replace

foreach v of var impact_est-wgt_max impact_est_var {
  tabstat `v', by(result) s(N mean p50 sd) nototal
  graph bar (mean) `v', over(estimator) asyvar over(dgp) title(`v') ytitle("") xsize(12) name(`v', replace)
  graph export "sims/sim`sim'/sim`sim'd_`v'.png", replace
}

log close sim`sim'
