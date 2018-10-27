cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"
clear all
cls

local sim = 1
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
local commonopts atet iter(50) cformat(%9.3fc) pformat(%5.3f) sformat(%7.3f) pooledvariance

// control simulations
set seed `sim'  //  1 for simulation 1, 2 for simulation 2, etc.
local reps 1000

// ------------------------------------------------------------------------
// 1-A Variouis sample sizes with impact = -.10
// ------------------------------------------------------------------------

parallel sim, expr(_b) reps(`reps') processors(1): ///
  onerep A, impact(-.075) estimators(cbps) ///
  n(8000 500 4000 2000 1750 1500 1250 1000 800 600 400 200 100) `commonopts'

sim_reshape
save sims/sim`sim'/sim`sim'a.dta, replace

foreach v of var impact_est-wgt_max {
  tabstat `v', by(N) s(N mean p50 sd)
}

graph bar (mean) reject_0, over(N) name(g1a)
graph export "sims/sim`sim'/sim`sim'a_power.png", replace

// ------------------------------------------------------------------------
// 1-B Variouis imapct sizes with N=2000
// ------------------------------------------------------------------------

parallel sim, expr(_b) reps(`reps') processors(1): ///
  onerep A, impact(-.05(-.01)-.15 -.07(-.005)-.06) estimators(cbps) n(2000) `commonopts'

sim_reshape
save sims/sim`sim'/sim`sim'b.dta, replace

foreach v of var impact_est-wgt_max {
  tabstat `v', by(true) s(N mean p50 sd)
}

graph bar (mean) reject_0, over(true) name(g1b)
graph export "sims/sim`sim'/sim`sim'b_power.png", replace

// ------------------------------------------------------------------------
// 1-C All DPGs (A to G) with n=2000, impact = -.075
// ------------------------------------------------------------------------

parallel sim, expr(_b) reps(`reps') processors(1): ///
  onerep A B C D E F G, impact(-.075) estimators(cbps) n(2000) `commonopts'

sim_reshape
save sims/sim`sim'/sim`sim'c.dta, replace

foreach v of var impact_est-wgt_max {
  tabstat `v', by(dgp) s(N mean p50 sd)
}

graph bar (mean) reject_0, over(dgp) name(g1c)
graph export "sims/sim`sim'/sim`sim'c_power.png", replace

log close sim`sim'
