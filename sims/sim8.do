if ("`: environment computername'"=="M116" | "`: environment computername'"=="NJ1STATA1") cd "U:\Stata\Ado\Devel\gmatch"
else cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"

clear all
cls

local sim = 8
cap log close sim`sim'
set linesize 180
cap mkdir  sims/sim`sim'
log using "sims/sim`sim'/logfile.log", replace name(sim`sim')

*! $Id$
*! BPO CBPS Simulation
// Try using the K-S and I-R simulation setup
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

// control simulations
set seed `sim'  //  1 for simulation 1, 2 for simulation 2, etc.
local reps 500
set matsize 5000
set maxiter 80

// options
parallel setclusters `=min(6, c(processors_max)-1)'
local commonopts ate vce(robust) iter(`c(maxiter)') cformat(%9.3fc) pformat(%5.3f) sformat(%7.3f) pooledvariance


// ------------------------------------------------------------------------
// Various estmators K-S and I-R setup
// ------------------------------------------------------------------------

local subsection A

onerep_ksir, n(200 1000) ///
  estimators(ipw ipwcbps cbps) cvtargets(99 90 80) `commonopts'

parallel sim, expr(_b) reps(`reps') processors(4): `e(cmdline)'

sim_reshape, dropaugsuffix
gen subsection = "Means only"
save sims/sim`sim'/Data_`subsection'.dta, replace

// ------------------------------------------------------------------------
// Various estmators K-S and I-R setup
// augmented
// ------------------------------------------------------------------------

local subsection B

onerep_ksir, n(200 1000) ///
  estimators(ipw ipwcbps cbps) cvtargets(99 90 80) `commonopts' augmented

parallel sim, expr(_b) reps(`reps') processors(4): `e(cmdline)'

sim_reshape, dropaugsuffix
gen subsection = "augmented"
save sims/sim`sim'/Data_`subsection'.dta, replace

// ------------------------------------------------------------------------
// Various estmators K-S and I-R setup
// augmented truepscore
// ------------------------------------------------------------------------

local subsection C

onerep_ksir, n(200 1000) ///
  estimators(ipw ipwcbps cbps) cvtargets(99 90 80) `commonopts' augmented truepscore

parallel sim, expr(_b) reps(`reps') processors(4): `e(cmdline)'

sim_reshape, dropaugsuffix
gen subsection = "augmented truepscore"
save sims/sim`sim'/Data_`subsection'.dta, replace

// ------------------------------------------------------------------------
// Various estmators K-S and I-R setup
// trueoutcome
// ------------------------------------------------------------------------

local subsection D

onerep_ksir, n(200 1000) ///
  estimators(ipw ipwcbps cbps) cvtargets(99 90 80) `commonopts' trueoutcome

parallel sim, expr(_b) reps(`reps') processors(4): `e(cmdline)'

sim_reshape, dropaugsuffix
gen subsection = "trueoutcome"
save sims/sim`sim'/Data_`subsection'.dta, replace

// ------------------------------------------------------------------------
// Various estmators K-S and I-R setup
// trueoutcome truepscore
// ------------------------------------------------------------------------

local subsection E

onerep_ksir, n(200 1000) ///
  estimators(ipw ipwcbps cbps) cvtargets(99 90 80) `commonopts' trueoutcome truepscore

parallel sim, expr(_b) reps(`reps') processors(4): `e(cmdline)'

sim_reshape, dropaugsuffix
gen subsection = "trueoutcome truepscore"
save sims/sim`sim'/Data_`subsection'.dta, replace

// ------------------------------------------------------------------------
// Combine results and summarize
// ------------------------------------------------------------------------

clear
append using ///
  sims/sim`sim'/Data_A.dta ///
  sims/sim`sim'/Data_B.dta ///
  sims/sim`sim'/Data_C.dta ///
  sims/sim`sim'/Data_D.dta ///
  sims/sim`sim'/Data_E.dta ///
  , gen(append)
labmask append, values(subsection)

foreach v of var impact_est-wgt_max impact_est_var {
  bys append: tabstat `v', by(result) s(N mean p50 sd) nototal
  graph bar (mean) `v', asyvar over(append) by(estimator, row(1) noyrescale title(`v')) ytitle("") xsize(12) name(`v', replace)
  graph export "sims/sim`sim'/Figure_`subsection'_`v'.png", replace
}

log close sim`sim'
