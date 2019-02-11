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
local reps  1000
set matsize 5000
set maxiter 50
set scheme mpr

// options
parallel setclusters `=min(6, c(processors_max)-1)'
local simopts    expr(_b) reps(\`reps') processors(`=c(processors_max)')
local commonopts n(50 200 1000) ///
                 estimators(ipw ipwcbps cbps) ///
         cvtargets(99 98 97 95 93 90 85 80 75 50) ///
         ate ///
         vce(robust) iter(\`c(maxiter)') cformat(%9.3fc) pformat(%5.3f) sformat(%7.3f) pooledvariance ///
         quietly

// ------------------------------------------------------------------------
// difference in means
// ------------------------------------------------------------------------

local subsection A

onerep_ksir, `commonopts'

parallel sim, `simopts': onerep_ksir, `commonopts'

sim_reshape, dropaugsuffix
gen subsection = "Means only"
save sims/sim`sim'/Data_`subsection'.dta, replace


// ------------------------------------------------------------------------
// augmented
// ------------------------------------------------------------------------

local subsection B

onerep_ksir, `commonopts' augmented

parallel sim, `simopts': onerep_ksir, `commonopts' augmented

sim_reshape, dropaugsuffix
gen subsection = "augmented"
save sims/sim`sim'/Data_`subsection'.dta, replace

// ------------------------------------------------------------------------
// augmented truepscore
// ------------------------------------------------------------------------

local subsection C

onerep_ksir, `commonopts' augmented truepscore

parallel sim, `simopts': onerep_ksir, `commonopts' augmented truepscore

sim_reshape, dropaugsuffix
gen subsection = "truepscore"
save sims/sim`sim'/Data_`subsection'.dta, replace

// ------------------------------------------------------------------------
// trueoutcome
// ------------------------------------------------------------------------

local subsection D

onerep_ksir, `commonopts' trueoutcome

parallel sim, `simopts': onerep_ksir, `commonopts' trueoutcome

sim_reshape, dropaugsuffix
gen subsection = "trueoutcome"
save sims/sim`sim'/Data_`subsection'.dta, replace

// ------------------------------------------------------------------------
// trueoutcome truepscore
// ------------------------------------------------------------------------

local subsection E

onerep_ksir, `commonopts' trueoutcome truepscore

parallel sim, `simopts': `e(cmdline)'

sim_reshape, dropaugsuffix
gen subsection = "trueoutcome truepscore"
save sims/sim`sim'/Data_`subsection'.dta, replace

// ------------------------------------------------------------------------
// Combine results and summarize
// ------------------------------------------------------------------------

/* I mannually added est_string in each dataset */

clear
append using ///
  sims/sim`sim'/Data_A.dta ///
  sims/sim`sim'/Data_B.dta ///
  sims/sim`sim'/Data_C.dta ///
  sims/sim`sim'/Data_D.dta ///
  sims/sim`sim'/Data_E.dta ///
  , gen(sctn)

labmask sctn, values(subsection)

drop estimator
label define est_consol ///
 1 IPW ///
 2 CBPS ///
 3 IPWCBPS ///
 4 CBPS99 ///
 5 CBPS98 ///
 6 CBPS97 ///
 7 CBPS95 ///
 8 CBPS93 ///
 9 CBPS90 ///
 10 CBPS85 ///
 11 CBPS80 ///
 12 CBPS75 ///
 13 CBPS50
encode est_string, gen(estimator) label(est_consol)
drop est_string subsection

table estimator sctn, by(N) c(count error_sqr)

foreach v of var impact_est-wgt_max impact_est_var {
  // bys sctn: tabstat `v', by(result) s(N mean p50 sd) nototal
  table estimator sctn, by(N) c(mean error_sqr)
  graph bar (mean) `v', asyvar over(estimator) over(N)  legend(col(3)) by(sctn, row(1) yrescale title(`v')) ytitle("") xsize(12) name(`v', replace)
  graph export "sims/sim`sim'/Figure_`subsection'_`v'.png", replace
}

log close sim`sim'
