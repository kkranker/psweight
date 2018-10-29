cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"
clear all
cls

local sim = 3
cap log close sim`sim'
set linesize 180
cap mkdir  sims/sim`sim'
log using "sims/sim`sim'/logfile.log", replace name(sim`sim')

*! $Id$
*! Simulation setup #1
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
parallel setclusters `c(processors_mach)'
local commonopts atet iter(20) cformat(%9.3fc) pformat(%5.3f) sformat(%7.3f) pooledvariance

// control simulations
set seed `sim'  //  1 for simulation 1, 2 for simulation 2, etc.
local reps 200

local reps 4

// ------------------------------------------------------------------------
// 3-A Variouis estmators with impact = -.09, N=2000
// ------------------------------------------------------------------------

onerep A, impact(-0.09) n(2000) ///
    estimators(raw ipw_true_ps ipw cbps stdprogdiff) ///
    cvtargets(99.9999 99.999 99.99 99.9(-.1)99 97 95(-10)85) `commonopts'

parallel sim, expr(_b) reps(`reps') processors(1): ///
  onerep A, impact(-0.09) n(2000) ///
     estimators(raw ipw_true_ps ipw cbps stdprogdiff) ///
     cvtargets(99.9999 99.999 99.99 99.9(-.1)99 97 95(-10)85) `commonopts'

sim_reshape
save sims/sim`sim'/sim`sim'a.dta, replace

foreach v of var impact_est-wgt_max {
  tabstat `v', by(estimator) s(N mean p50 sd)
  graph bar (mean) `v', over(estimator, label(alternate labsize(tiny))) name(`v', replace)
  graph export "sims/sim`sim'/sim`sim'a_`v'.png", replace
}

log close sim`sim'
