cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"
clear all
cls

local sim = 2
cap log close sim`sim'
set linesize 180
cap mkdir  sims/sim`sim'
log using "sims/sim`sim'/logfile.log", replace name(sim`sim')

*! $Id$
*! Simulation setup #1
// Here I'm trying to find combinations of N/impact that are "marginally" powered
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
parallel setclusters `c(processors_mach)'
local commonopts atet iter(20) cformat(%9.3fc) pformat(%5.3f) sformat(%7.3f) pooledvariance

// control simulations
set seed `sim'  //  1 for simulation 1, 2 for simulation 2, etc.
local reps 1000

// ------------------------------------------------------------------------
// 2-A Variouis estmators with impact = -.075, N=2000
// ------------------------------------------------------------------------

onerep A, impact(-.075) n(2000) ///
     estimators(raw ipw_true_ps ipw ipw_te cbps stdprogdiff) ///
     cvtargets(99 97 95(-5)70) `commonopts'

parallel sim, expr(_b) reps(`reps') processors(1): ///
  onerep A, impact(-.075) n(2000) ///
     estimators(raw ipw_true_ps ipw ipw_te cbps stdprogdiff) ///
     cvtargets(99 97 95(-5)70) `commonopts'

sim_reshape
save sims/sim`sim'/sim`sim'a.dta, replace

foreach v of var impact_est-wgt_max {
  tabstat `v', by(estimator) s(N mean p50 sd)
  graph bar (mean) `v', over(estimator) name(`v', replace)
  graph export "sims/sim`sim'/sim`sim'a_`v'.png", replace
}

log close sim`sim'
