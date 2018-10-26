cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"
clear all
cls

cap log close sim1
//set linesize 180
cap mkdir  sims/sim1
log using "sims/sim1/logfile.log", replace name(sim1)

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

// DGP is saved in a separate .ado file
  adopath ++ "./sims"
  which dgp_ssbgc
  which onerep
  which sim_reshape

// options
  parallel setclusters `c(processors_mach)'
  set cformat %9.3fc
  set pformat %5.3f
  set sformat %7.3f
  set type float
  set tracedepth 1
  set maxiter 50

// control simulations
set seed 1  //  1 for simulation 1, 2 for simulation 2, etc.
local reps 500

// ------------------------------------------------------------------------
// 1-A Variouis sample sizes with impact = .10
// ------------------------------------------------------------------------

parallel sim, expr(_b) reps(`reps') processors(1): ///
  onerep A, impact(-.10) estimators(ipw) n(8000 500 4000 2000 1750 1500 1250 1000 800 600 400 200 100)

sim_reshape
save sims/sim1/sim1a.dta, replace
graph bar (mean) reject_0, over(N) over(estimator) name(g1a)
graph export sims/sim1/sim1a_power.png, replace


// ------------------------------------------------------------------------
// 1-B Variouis sample sizes with impact = .10
// ------------------------------------------------------------------------

parallel sim, expr(_b) reps(`reps') processors(1): ///
  onerep A, impact(-.05(-.01)-.15) estimators(ipw) n(2000)

sim_reshape
save sims/sim1/sim1b.dta, replace
graph bar (mean) reject_0, over(true) name(g1b)
graph export sims/sim1/sim1b_power.png, replace


// ------------------------------------------------------------------------
// 1-C Variouis sample sizes with impact = .10
// ------------------------------------------------------------------------

parallel sim, expr(_b) reps(`reps') processors(1): ///
  onerep A B C D E F G, impact(-.07) estimators(ipw) n(2000)

sim_reshape
save sims/sim1/sim1c.dta, replace
graph bar (mean) reject_0, over(dgp) name(g1c)
graph export sims/sim1/sim1c_power.png, replace

log close sim1
