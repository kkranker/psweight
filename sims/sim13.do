if ("`: environment computername'"=="M116" | "`: environment computername'"=="NJ1STATA1") cd "U:\Stata\Ado\Devel\gmatch"
else cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"

clear all
cls

local sim = 13
cap log close sim`sim'
set linesize 180
cap mkdir  sims/sim`sim'
log using "sims/sim`sim'/logfile.log", replace name(sim`sim')

*! $Id$
*! BPO CBPS Simulation
// See if the typo on page 215 of the I-R article explains differences between our and their results.
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
which dgp_ksir
which onerep_ksir_means
which sim_reshape

// control simulations
local reps 2500

// other settings
set seed `sim'  //  1 for simulation 1, 2 for simulation 2, etc.
set matsize `c(max_matsize)'
set maxiter 50
set scheme mpr

// ------------------------------------------------------------------------
// setup simulations
// ------------------------------------------------------------------------

// options
parallel setclusters `=min(6, c(processors_max)-1)'
local simopts    expr(_b) reps(\`reps') processors(`=c(processors_max)')
local commonopts n(200 1000) ///
                 estimators(ipw_true_ps ipw ipwcbps cbps) ///
                 iter(`c(maxiter)') cformat(%9.3fc) pformat(%5.3f) sformat(%7.3f) ///
                 quietly

// ------------------------------------------------------------------------
// multiple matching approaches
// difference in means and augmented
// with the I-R definition of X4
// ------------------------------------------------------------------------

onerep_ksir_means, `commonopts' irversion
drop _all

parallel sim, `simopts': onerep_ksir_means, `commonopts' irversion

qui compress
save sims/sim`sim'/Data_Unprocessed_1.dta, replace

sim_reshape
gen irversion = 1

tempfile f1
save "`f1'", replace
drop _all

// ------------------------------------------------------------------------
// multiple matching approaches
// difference in means and augmented
// with the K-S definition of X4
// ------------------------------------------------------------------------

onerep_ksir_means, `commonopts'
drop _all

parallel sim, `simopts': onerep_ksir_means, `commonopts'

qui compress
save sims/sim`sim'/Data_Unprocessed_0.dta, replace

sim_reshape
gen irversion = 0


// ------------------------------------------------------------------------
// summarize results
// ------------------------------------------------------------------------

label define irversion 0 "Kang-Shafler" 1 "Ima-Ratkovic"
label val irversion irversion
append using "`f1'"

qui compress
save sims/sim`sim'/Data.dta, replace

assert augmented[1]==augmented
assert augmented[1]==0
label define aug 0 "IPW-POP", replace

table irversion estimator, by(N) c(count bias)

foreach v of var impact_est rmse impact_est_var bias rmse  pct_bias mae {
  di _n(2) as res `"`v'  `:var lab `v''"'
  table irversion estimator, by(N) c(mean `v')
}

log close sim`sim'
