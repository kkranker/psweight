if ("`: environment computername'"=="M116" | "`: environment computername'"=="NJ1STATA1") cd "U:\Stata\Ado\Devel\gmatch"
else cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"

clear all
cls

local sim = 9
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
which elasticregress
which randomforest

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
                 estimators(raw ipw_true_ps ipw ipwcbps cbps elastic rf) ///
         cvtargets(99 98 97 95 93 90 85 80 75 50) ///
         ate ///
         vce(robust) iter(\`c(maxiter)') cformat(%9.3fc) pformat(%5.3f) sformat(%7.3f) pooledvariance ///
         alpha(.5) numfolds(10) ///
         iterations(100) depth(12) lsize(5) numvars(4) ///
         quietly

// ------------------------------------------------------------------------
// difference in means
// and augmented
// ------------------------------------------------------------------------

onerep_ksir, `commonopts' augmented

parallel sim, `simopts': onerep_ksir, `commonopts' augmented

sim_reshape
table estimator augmented N, c(count error_sqr)

foreach v of var impact_est-wgt_max impact_est_var {
  // bys augmented: tabstat `v', by(result) s(N mean p50 sd) nototal
  table estimator augmented N, c(mean error_sqr)
  graph bar (mean) `v', asyvar over(estimator) over(N) legend(col(3)) by(augmented, row(1) yrescale title(`v')) ytitle("") xsize(12) name(`v', replace)
  graph export "sims/sim`sim'/Figure_`v'.png", replace
}

log close sim`sim'
