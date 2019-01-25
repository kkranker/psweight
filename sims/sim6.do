if ("`: environment computername'"=="M116" | "`: environment computername'"=="NJ1STATA1") cd "U:\Stata\Ado\Devel\gmatch"
else cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"

clear all
cls

local sim = 6
cap log close sim`sim'
set linesize 180
cap mkdir  sims/sim`sim'
log using "sims/sim`sim'/logfile.log", replace name(sim`sim')

*! $Id$
*! BPO CBPS Simulation
// Refining the Simulation -- add nois in one of the confounders
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
parallel setclusters `=min(6, c(processors_max)-1)'
local commonopts atet vce(robust) iter(80) cformat(%9.3fc) pformat(%5.3f) sformat(%7.3f) pooledvariance

// control simulations
set seed `sim'  //  1 for simulation 1, 2 for simulation 2, etc.
local reps 400
set matsize 5000


// ------------------------------------------------------------------------
// A Variouis estmators with impact = -.09, N=2000, wnois(.10)
// ------------------------------------------------------------------------

local subsection A

onerep G K1, impact(-0.09) n(2000) ///
  estimators(cbps) cvtargets(99 85 50) `commonopts' histogram

parallel sim, expr(_b) reps(`reps') processors(4): ///
  onerep G K1, impact(-0.09) n(2000) ///
    estimators(cbps) cvtargets(99 85 50) `commonopts'

sim_reshape
save sims/sim`sim'/Data_`subsection'.dta, replace

foreach v of var impact_est-wgt_max impact_est_var {
  tabstat `v', by(result) s(N mean p50 sd) nototal
  graph bar (mean) `v', over(estimator) asyvar over(dgp) title(`v') ytitle("") xsize(12) name(`v', replace)
  graph export "sims/sim`sim'/Figure_`subsection'_`v'.png", replace
}

log close sim`sim'
