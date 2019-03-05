if ("`: environment computername'"=="M116" | "`: environment computername'"=="NJ1STATA1") cd "U:\Stata\Ado\Devel\gmatch"
else cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"

clear all
cls

// simulation folder with data
local sim = 11

cap log close sim`sim'_tables
set linesize 180
cap mkdir  sims/sim`sim'
log using "sims/sim`sim'/summary_tables.log", replace name(sim`sim'_tables)

*! $Id$
*! BPO CBPS Simulation
// Repackage the results from simulation #11 into ready-to-use tables.
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.

mac list _sim
qui adopath ++ ./sims/

// there was a mistake in the old version on sim_reshape.ado; I'm re-running here with the fixed version
if 0 {
  use using sims/sim`sim'/Data_Unprocessed.dta, clear
  sim_reshape
  save sims/sim`sim'/Data.dta, replace
}
else {
  use sims/sim`sim'/Data.dta, clear
}

// sort variables into groups
unab stats : impact_est-wgt_kurtosis Nt
unab sumstats : rmse impact_est_var
unab scenario: result result dataset dgp N estimator augmented true
unab  allleft : _all
local allleft : list allleft - stats
local allleft : list allleft - sumstats
local allleft : list allleft - scenario
mac list _allleft
assert "`allleft'"=="rep"

// update labels
label list aug est
/*     aug:
           0 Unadjusted
           1 Reg. adjusted
       est:
           1 RAW
           2 IPW_TRUE_PS
           3 STDPROGDIFF
           4 IPW_TE
           5 ELASTIC
           6 RF
           7 IPWCBPS
           8 CBPS
           9 CBPS99_9
          10 CBPS99
          11 CBPS98
          12 CBPS97
          13 CBPS95
          14 CBPS93
          15 CBPS90
          16 CBPS85
          17 CBPS80
          18 CBPS75
          19 CBPS50
          20 IPW     */

label define aug 0 "IPW" 1 "WLS", modify

recode estimator ///
           (20 =   1 "Logit") ///
           ( 6 =   2 "RF") ///
           ( 5 =   3 "ELASTIC") ///
           ( 3 =   4 "PROG") ///
           ( 7 =   5 "CBPS-O") ///
           ( 8 =   6 "CBPS-J") ///
           ( 9 =   7 "CBPS-99.9%") ///
           (10 =   8 "CBPS-99%") ///
           (11 =   9 "CBPS-98%") ///
           (12 =  10 "CBPS-97%") ///
           (13 =  11 "CBPS-95%") ///
           (14 =  12 "CBPS-93%") ///
           (15 =  13 "CBPS-90%") ///
           (16 =  14 "CBPS-85%") ///
           (17 =  15 "CBPS-80%") ///
           (18 =  16 "CBPS-75%") ///
           (19 =  17 "CBPS-50%") ///
           ( 2 =  18 "TRUE") ///
           ( 4 =  19 "IPW_TE") ///
           ( 1 =  20 "RAW") ///
           (else = 999 "something is missing here") ///
           , gen(estimator_recode)
recode estimator_ ///
             (  6 = 100  ) ///
             (  7 =  99.9) ///
             (  8 =  99  ) ///
             (  9 =  98  ) ///
             ( 10 =  97  ) ///
             ( 11 =  95  ) ///
             ( 12 =  93  ) ///
             ( 13 =  90  ) ///
             ( 14 =  85  ) ///
             ( 15 =  80  ) ///
             ( 16 =  75  ) ///
             ( 17 =  50  ) ///
             (else = 0), ///
             gen(cbps_pct)
tablist estimator_ cbps_pct estimator, sep(0) sort(v)
assert estimator_recode != 99
order  estimator_ , after(estimator)
drop estimator
rename estimator_ estimator
local scenario `macval(scenario)' cbps_pct

// drop the PROG estimator
// the PROG estimator was a disaster -- converged rarely and not at all with 3 of the 5 Ns
table  N aug estimator if estimator==4, row col c(count rep count impact_est)
drop   if estimator==4

// box plots of impact estimates
levelsof N, local(Narray)
set scheme s2color
foreach i of local Narray {
  **graph hbox impact_est if N==`i', over(estimator) by(augmented, compact) nooutsides name(n_`i')
}

// collapse down to 1 row per scenario/estimator
isid rep `scenario'
collapse (count) n_reps=impact_est n_reps_attempt=rep ///
         (mean) `stats' ///
         (firstnm) `sumstats' ///
         (max) rep ///
         , by(result `scenario')

// format numbers so tables look good
replace reject_0 = reject_0*100
replace covered  = covered*100
format %7.0fc n_reps*
format %7.3fc power_zstat_0 p_0 bal_max_asd bal_mean_asd wgt_sd wgt_cv wgt_skewness wgt_kurtosis
format %7.2fc impact_est sd_error bias error_sqr rmse impact_est_var
format %7.1fc reject_0 covered wgt_max

local stats rmse // bias  impact_est_var


// wide tables - columns are N's and reg-adjusted
foreach v of local stats {
  di _n(2) as res `"`v'  `:var lab `v''"'
  tabdisp estimator N augmented, c(`v')
}

// tall tables - columns are N's and reg-adjusted
foreach v of local stats {
  di _n(2) as res `"`v'  `:var lab `v''"'
  tabdisp augmented N, by(estimator) c(`v')
}

// separate tables for aug=0/1
foreach v of local stats {
  di _n(2) as res `"`v'  `:var lab `v''"'
  bys augmented: tabdisp estimator N, c(`v')
}



log close sim`sim'_tables
