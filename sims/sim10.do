if ("`: environment computername'"=="M116" | "`: environment computername'"=="NJ1STATA1") cd "U:\Stata\Ado\Devel\gmatch"
else cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"

clear all
cls

local sim = 10
cap log close sim`sim'
set linesize 180
cap mkdir  sims/sim`sim'
log using "sims/sim`sim'/logfile.log", replace name(sim`sim')

*! $Id:$
*! BPO CBPS Simulation
// Try using the K-S and I-R simulation setup
// Compare to IPW, CBPS (both versions), RF
// RF settings picked through X-validation
//
*! By Keith Kranker
// Last updated $Date:$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.

mac list _sim

// DGP is saved in a separate .ado file
adopath ++ "./sims"
adopath ++ "../gridsearchcv"
which dgp_ssbgc
which onerep
which sim_reshape
which randomforest
which gridsearchcv

// control simulations
local reps 2000
local Nrange 25 50 100 200 1000
local cvreps 100
local ntrees 250
local grid (depth: 1 5 10 20) (lsize: 5 10 20 40) (numvars: 1 2 3 4)

// other settings
set seed `sim'  //  1 for simulation 1, 2 for simulation 2, etc.
set matsize 5000
set maxiter 50
set scheme mpr


// ------------------------------------------------------------------------
// cross validation for RF model parameters
// ------------------------------------------------------------------------

// I'm being lazy and just using a global to capture the results
mata: resulttable = J(0,3,"")

// run once
rf_xval_ksir `grid', n(`Nrange') rfopts(iterations(`ntrees'))
mata: resulttable

simulate, reps(`cvreps'): `e(cmd)'
mata: resulttable

// move results from Mata global back into Stata
drop _all
mata: stata("set obs " + strofreal(rows(resulttable)))
gen strL N = ""
gen strL score = ""
gen strL params = ""
mata: st_sstore(.,.,resulttable)
save "cv_results.dta", replace

// find the one that was picked the most
destring N score, replace
bysort N (params): gen denom  = _N
by     N  params : gen n_best = _N
gsort  N -n_best -score
by N: gen tag_best_n = _n == 1
// best param for each option
list N params n_best denom if tag_best_n, sep(0) noobs

drop tag_best_n denom n_best
bys params : gen n_best = _n
gsort -n_best -score
gen denom = _N
keep in 1
// best overall
list params n_best denom, sep(0) noobs
local bestparams = params[1]


// ------------------------------------------------------------------------
// setup simulations
// ------------------------------------------------------------------------

// options
parallel setclusters `=min(6, c(processors_max)-1)'
local simopts    expr(_b) reps(\`reps') processors(`=c(processors_max)')
local commonopts n(`Nrange') ///
                 estimators(raw ipw_true_ps ipw ipwcbps cbps rf) ///
         cvtargets(99.9 99 98 97 95 93 90 85 80 75 50) ///
         ate ///
         vce(robust) iter(\`c(maxiter)') cformat(%9.3fc) pformat(%5.3f) sformat(%7.3f) pooledvariance ///
         rfopts(iterations(200) `bestparams') ///
         // quietly

// ------------------------------------------------------------------------
// multiple matching approaches
// difference in means and augmented
// ------------------------------------------------------------------------

onerep_ksir, `commonopts' augmented

parallel sim, `simopts': onerep_ksir, `commonopts' augmented


// ------------------------------------------------------------------------
// summarize results
// ------------------------------------------------------------------------

set linesize 220

sim_reshape
save sims/sim`sim'/Data.dta, replace

table estimator augmented N, c(count error_sqr)

foreach v of var impact_est-wgt_max impact_est_var {
  di _n(2) as res `"`v'  `:var lab `v''"'
  table estimator augmented N, c(mean `v')
  graph bar (mean) `v', asyvar over(estimator) over(N) legend(col(3)) by(augmented, row(1) yrescale title(`v')) ytitle("") xsize(12) name(`v', replace)
  graph export "sims/sim`sim'/Figure_`v'.png", replace
}

log close sim`sim'
