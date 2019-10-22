if ("`: environment computername'"=="M116" | "`: environment computername'"=="NJ1STATA1") cd "U:\Stata\Ado\Devel\gmatch"
else cd "C:\Users\kkranker\Documents\Stata\psweight\code-psweight"

clear all
cls

// simulation folder with data
local sim = 11

cap log close sim`sim'_tables
set linesize 180
cap mkdir  sims/sim`sim'
log using "sims/sim`sim'/summary_tables.log", replace name(sim`sim'_tables)

*! BPO CBPS Simulation
// Repackage the results from simulation #11 into ready-to-use tables.
//
*! By Keith Kranker
//
// Note - This code was "hacked." I ran most of the simulations back for the first paper draft.
// Then, after we got the R&R from TAS, I ran ***ONLY*** three additional estimators.
// Then in this file I patched the old and new data files together and summarized the results.
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.

mac list _sim
qui adopath ++ ./sims/

// there was a mistake in the old version on sim_reshape.ado; I'm re-running here with the fixed version
use sims/sim`sim'/Data_append.dta, clear


// sort variables into groups
unab stats : impact_est-wgt_kurtosis Nt mae pct_bias
unab sumstats : rmse impact_est_var
unab scenario: result dataset dgp N estimator augmented true
unab  allleft : _all
local allleft : list allleft - stats
local allleft : list allleft - sumstats
local allleft : list allleft - scenario
mac list _allleft
assert "`allleft'"=="rep"

// update labels
label list aug est
qui {
  /*
  aug:
           0 Unadjusted
           1 Reg. adjusted
  est:
           1 RAW
           2 IPW_TRUE_PS
           3 STDPROGDIFF
           4 IPW_TE
           5 IPW
           6 TRIM
           7 EBALANCE
           8 RF
           9 IPWCBPS
          10 CBPS
          11 CBPS99_9
          12 CBPS99
          13 CBPS98
          14 CBPS97
          15 CBPS95
          16 CBPS93
          17 CBPS90
          18 CBPS85
          19 CBPS80
          20 CBPS75
          21 CBPS50
          22 DISCARD
  */
}
label define aug 0 "IPW" 1 "WLS", modify

recode estimator ///
           ( 5 =   1 "Logit regression") ///
           ( 6 =   2 "Logit regression, with trimming") ///
           ( 8 =   4 "Random forest classifier") ///
           ( 9 =   5 "CBPS, over identified") ///
           (10 =   6 "CBPS, just identified") ///
           (12 =   7 "PCBPS 99%") ///
           (14 =   8 "PCBPS 97%") ///
           (15 =   9 "PCBPS 95%") ///
           (17 =  10 "PCBPS 90%") ///
           (20 =  11 "PCBPS 75%") ///
           ( 2 =  12 "True propensity score") ///
           ( 1 =  13 "None (unweighted data)") ///
           (else = 999 "<discarded>") ///
           , gen(estimator_recode)
recode estimator_ ///
             (  6 = 100  ) ///
             (  7 =  99  ) ///
             (  8 =  97  ) ///
             (  9 =  95  ) ///
             ( 10 =  90  ) ///
             ( 11 =  75  ) ///
             (else = 0), ///
             gen(cbps_pct)
tablist estimator_ cbps_pct estimator, sep(0) sort(v)
order  estimator_ , after(estimator)
rename estimator old_estimator
local scenario `macval(scenario)' cbps_pct

// select cells to show in the final tables
// drop "prog";
// drop several different % for PCBPS -- all these different targets were overkill
// the 99.9 doesn't clearly do anything -- virtually the same as CBPS in many settings
// also drop N = 25 / 100 because our editor wanted a simpler table
// in some tables/figures, only show 2 sample sizes
local ifstmnt !inlist(estimator, 999) & inlist(N, 50, 200, 1000)
local 2ns     inlist(N, 50)

// box plots of impact estimates
levelsof N, local(Narray)
set scheme s2manual_KAK
graph hbox impact_est if `ifstmnt' & `2ns', ///
  over(estimator) nooutsides ///
  ytitle("Distribution of impact estimates" "(true impact = 0)", size(vsmall)) ysize(6.5) ylab(0) xsize(6.5) ///
  by(N augmented, norescale colfirst row(2) iscale(*.6) iylabel title("") note("Outside values not shown", size(vsmall))) note("")
graph save   sims/sim`sim'/hbox_impact_est.gph, replace
graph export sims/sim`sim'/hbox_impact_est.emf, replace
graph export sims/sim`sim'/hbox_impact_est.png, replace

// collapse down to 1 row per scenario/estimator
isid rep `scenario'
collapse (count) n_reps=impact_est n_reps_attempt=rep ///
         (mean) `stats' ///
         (firstnm) `sumstats' ///
         (max) n_reps_attempt_alt = rep ///
         , by(`scenario' old_estimator)
assert n_reps_attempt_alt == n_reps_attempt
drop   n_reps_attempt_alt

// format numbers so tables look good
replace reject_0 = reject_0*100
replace covered  = covered*100
format %7.1fc reject_0 covered wgt_max
format %7.0fc n_reps*
format %7.3fc power_zstat_0 p_0 bal_max_asd bal_mean_asd wgt_sd wgt_cv wgt_skewness wgt_kurtosis
format %7.2fc impact_est sd_error bias error_sqr rmse impact_est_var

// save dataset
save "sims/sim`sim'/1-row summary.dta", replace

// wide tables with ALL the results - columns are N's and reg-adjusted
foreach v of var n_reps* `stats' `sumstats' {
  di _n(2) as res `"`v'  `:var lab `v''"'
  tabdisp old_estimator N augmented, c(`v') format(`:format `v'') stubwidth(13)
}


di _n(10) "Tables with selected results" _n(10)

// tall table with columns as N's
foreach v of var rmse bias impact_est_var {
  di _n(2) as res `"`v'  `:var lab `v''"'
  preserve
  keep if `ifstmnt'
  keep    estimator N augmented `v'
  tabdisp estimator N augmented, c(`v') format(`:format `v'') cellwidth(6) csep(2) stubwidth(11)

  // rebuild table and export to LaTeX snippet
  isid estimator augmented N
  rename `v' `v'_
  reshape wide `v'_ , i(estimator augmented) j(N)
  rename (`v'_*) (`v'_*_)
  reshape wide `v'_*, i(estimator) j(augmented)
  unab vvv: `v'_*
  mata: st_matrix("dta", st_data(., "`vvv'"))
  mata: st_matrixcolstripe("dta", (J(`:list sizeof vvv', 1, ""), tokens("`vvv'")'))
  mata: st_matrixrowstripe("dta", (J(`c(N)', 1, ""), st_vlmap("estimator_recode", st_data(., "estimator"))))
  frmttable using "sims/sim`sim'/LaTeX summary.tex", s(dta) tex `=cond("`v'"=="rmse", "replace", "addtable")' fragment plain sdec(2) coljust(r{r}r) posttext("columns refer to N and 0/1 for regression adjustment.")
  restore
}

// table with columns for each statistic
preserve
  keep if `ifstmnt' & `2ns'
  keep estimator N augmented wgt_cv bal_mean_asd bal_max_asd
  keep if aug==0 // these stats are all identical for aug and non-aug estimators
  local c=0
  foreach v of var           wgt_cv bal_mean_asd bal_max_asd {
    rename `v' stat`++c'
    local def `macval(def)' `c' "`v'"
  }
  qui reshape long stat, i(estimator N augmented)
  label define statname `def'
  label val _j statname
  tabdisp estimator  N _j, by(augmented) c(stat) cellwidth(6) csep(2) stubwidth(13)

// rebuild and export to LaTeX snippet
restore, preserve
  keep if `ifstmnt' & `2ns'
  keep estimator N augmented wgt_cv bal_mean_asd bal_max_asd
  keep if aug==0 // these stats are all identical for aug and non-aug estimators
  reshape wide wgt_cv bal_mean_asd bal_max_asd , i(estimator) j(N)
  unab vvv: wgt_cv* bal_mean_asd* bal_max_asd*
  mata: st_matrix("dta", st_data(., "`vvv'"))
  mata: st_matrixcolstripe("dta", (J(`:list sizeof vvv', 1, ""), tokens("`vvv'")'))
  mata: st_matrixrowstripe("dta", (J(`c(N)', 1, ""), st_vlmap("estimator_recode", st_data(., "estimator"))))
  frmttable using "sims/sim`sim'/LaTeX summary.tex", s(dta) tex addtable fragment sdec(3) coljust(r{r}r)

restore

log close sim`sim'_tables
