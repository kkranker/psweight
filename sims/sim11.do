if ("`: environment computername'"=="M116" | "`: environment computername'"=="NJ1STATA1") cd "U:\Stata\Ado\Devel\gmatch"
else cd "C:\Users\kkranker\Documents\Stata\psweight\code-psweight\"

clear all
cls

local sim = 11
cap log close sim`sim'
set linesize 180
cap mkdir  sims/sim`sim'
log using "sims/sim`sim'/logfile.log", append name(sim`sim')

*! $Id$
*! BPO CBPS Simulation
// Try using the K-S and I-R simulation setup
// Compare to IPW, CBPS (both versions), RF (with RF paremters picked through X-validation)
// (Sim 11 is what I was trying to run in Sim 10. I've just fixed the issue with matsize.)
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
adopath ++ "../gridsearchcv"
which dgp_ksir
which onerep_ksir
which sim_reshape
which randomforest

// control simulations
local reps 5000
local Nrange 25 50 100 200 1000
local cvreps 100
local ntrees 250
local grid (depth: 1 5 10 20) (lsize: 5 10 20 40) (numvars: 1 2 3 4)

// other settings
set seed `sim'  //  1 for simulation 1, 2 for simulation 2, etc.
set matsize `c(max_matsize)'
set maxiter 50
set scheme mpr


// ------------------------------------------------------------------------
// cross validation for RF model parameters
// ------------------------------------------------------------------------

// SEE sim10/logfile.log for cross validation results.
// (Sim 11 is what I was trying to run in Sim 10. I've just fixed the issue with matsize.)

local bestparams = "depth(5) lsize(20) numvars(4)"


// ------------------------------------------------------------------------
// setup simulations
// ------------------------------------------------------------------------

// options
parallel setclusters `=min(8, c(processors_max)-1)'
local simopts    expr(_b) reps(\`reps') processors(`=c(processors_max)')
local commonopts n(`Nrange') ///
                 estimators(ebalance trim discard) ///
         ate ///
         vce(robust) iter(\`c(maxiter)') cformat(%9.3fc) pformat(%5.3f) sformat(%7.3f) pooledvariance ///
         rfopts(iterations(200) `bestparams') ///
         // quietly

// ------------------------------------------------------------------------
// multiple matching approaches
// difference in means and augmented
// ------------------------------------------------------------------------

onerep_ksir, `commonopts' augmented
drop _all

parallel sim, `simopts': onerep_ksir, `commonopts' augmented

append using sims/sim`sim'/Data_Unprocessed.dta, gen(old_sims)

save sims/sim`sim'/Data_Unprocessed_append.dta, replace

// ------------------------------------------------------------------------
// summarize results
// ------------------------------------------------------------------------

sim_reshape
save sims/sim`sim'/Data_append.dta, replace

set linesize 220

table estimator augmented N, c(count error_sqr)

foreach v of var impact_est-wgt_max impact_est_var {
  di _n(2) as res `"`v'  `:var lab `v''"'
  table estimator augmented N, c(mean `v')
  graph bar (mean) `v', asyvar over(estimator) over(N) legend(col(3)) by(augmented, row(1) yrescale title(`v')) legend(col(5)) ytitle("") xsize(12) name(`v', replace)
  graph export "sims/sim`sim'/Figure_`v'.png", replace
}

log close sim`sim'
