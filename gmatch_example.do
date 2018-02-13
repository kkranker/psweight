mac drop _all
clear all
cls
set varabbrev off
set scheme mpr_blue
set linesize 160
set maxiter 100
cap log close gmatch_example
cap log close gmatch_example_ado
local makegraphs = 01
cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\"

log using gmatch_example.log, name(gmatch_example) replace

// Multiple-equation models: An introduction and potential applications to our work at Mathematica
// Design and methods “brown bag” workshop
// April 8, 2015
// Keith Kranker

// Stata_code_2_IPW.do This program includes all the examples in the powerpoint slides, plus more.
// This program includes examples of how to code inverse propensity weighting (IPW) estimators using Stata's GMM command

// include C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\gmatchclass.mata
do C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\cr_lgmatch.do

version 15.1
set type double
di as txt "Current user: `c(username)'" _n "Environment: `c(os)' `c(machine_type)' `: environment computername'" _n "Stata: `c(stata_version)'" cond(c(stata_version)==c(version),""," (set to version `c(version)')") _n "Date: " c(current_date) " " c(current_time)


************************************************************************************
* Describe/summarize the example datasets
************************************************************************************

*** Input data file (simple_cattaneo_data) comes from the program named Make_example_datasets.do (in C:\Users\kkranker\Documents\Stata\Multiple-Equation-Models)

use "C:\Users\kkranker\Documents\Stata\Multiple-Equation-Models\simple_cattaneo_data.dta"
desc, short
notes _dta
summ, sep(0)
tab treat treat_cat, mi
corr treat y1 y1_binary

local if if _n<=500
set seed 1
gen wgt = max(.1,rnormal(2,.4))
gen fwgt = round(rnormal(2,.4))
// forvalues i = 20/200 {
forvalues i = 90/95 {
  gen x`i' = rnormal()
}
// expand 5e4 if touse
// expand 1e5 if touse

// // give the sample poor overlap
// tab2    treat x1
// replace x1 = 1 if  treat & runiform()<.85
// replace x1 = 0 if !treat & runiform()<.85
// tab2    treat x1


local depvars = "y1 y1_binary"
local treatvar = "treat"
local varlist = "x1 i.x2 i.x3 x4 x5 x6 x7 x9*"
// local varlist = "x*"
// local varlist = "x1 ib0.x2"
//local wgtvar = "wgt"
local wgtvar = "wgt"
local tousevar = "touse"
local estimate = "atet"


// some automatic parsing based on options above
if "`wgtvar'"!="" local wgtexp "[iw=`wgtvar']"
mark    `tousevar' `if' `in' `wgtexp'
markout `tousevar' `depvars' `treatvar' `varlist'
_rmdcoll `treatvar' `varlist' if `tousevar' `wgtexp', expand
// _rmcoll `treatvar' `varlist' if `tousevar' `wgtexp', expand logit touse(`tousevar')
// fvexpand `varlist' if `tousevar'
local varlist `r(varlist)'
forvalues j=1/`: list sizeof varlist' {
  local v : word `j' of `varlist'
  _ms_parse_parts `v'
  if !r(omit) local varlist1 `"`varlist1' `v'"'
}
local varlist : copy local varlist1

fvrevar `varlist'
export delimited `treatvar' `r(varlist)' `wgtvar' using testfile.csv if `tousevar', replace nolabel

mata:

depvars  = st_local("depvars" )
treatvar = st_local("treatvar")
wgtvar   = st_local("wgtvar"  )
varlist  = st_local("varlist" )
tousevar = st_local("tousevar")
estimate = st_local("estimate")


// ****************************
// * UNWEIGHTED DATA EXAMPLES *
// ****************************

D = gmatch()
D.set(treatvar, varlist, tousevar)
if (depvars!="") D.set_Y(depvars,tousevar)

// Misc balance measures
  D.diff()
  D.stddiff()
  D.mean_asd()
  D.stddiff(1)
  D.stddiff(0)
  D.varratio()
  D.prognosticdiff()

  "Balance table before matching"
  table = D.balancetable(1)

// Replicate CBPS

  "--- ATE (not overidentified) ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , ate      logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbps("ate", "cbps", 2, 0)

  "--- ATE overidentified ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , ate over logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbps("ate", "cbps", 2, 1)

  "--- ATET (not overidentified) ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , att      logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbps("atet", "cbps", 2, 0)

  "--- ATET overidentified ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , att over logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbps("atet", "cbps", 2, 1)

// Other objective functions
  D.cbps("atet","mean_sd_sq",1)
  D.cbps("atet","sd_sq",1)

// IPW
  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations")
  stata("di _b[POmean:0.treat]")
  stata("tebalance summarize")
  stata("tebalance summarize, baseline")

  D.ipw("atet")
  D.pomean()
  table = D.balancetable(3)
  D.reweight()
  table = D.balancetable(3)

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , ate aequations")
  ipw = D.ipw("ate")

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations tlevel(0) control(1)")
  ipw = D.ipw("ateu")

// tradeoff between CBPS-like balance and variance in weights
  cbpsweight = D.cbps("atet","cbps", 1, 0)
  cbpsweight = D.cbps("atet","cbps", 1, 0, (1,.75,6))
  cbpsweight = D.cbps("atet","cbps", 1, 0, (1,.50,6))

  cbpsweight = D.cbps("atet","mean_sd_sq", 1, 0, (1,.75,6))
  cbpsweight = D.cbps("atet","mean_sd_sq", 1, 0, (1,.50,6))

mata drop D


// **************************
// * WEIGHTED DATA EXAMPLES *
// **************************

DW = gmatch()
DW.set(treatvar, varlist, tousevar, wgtvar)
if (depvars!="") DW.set_Y(depvars,tousevar)

// Misc balance measures
  DW.diff()
  DW.stddiff()
  DW.mean_asd()
  DW.stddiff(1)
  DW.stddiff(0)
  DW.varratio()
  DW.prognosticdiff()

  "Balance table before matching"
  temp = DW.balancetable(1)

// Replicate CBPS

  "--- ATE (not overidentified) ---"; ""; ""
    DW.cbps("ate", "cbps", 2, 0)

  "--- ATE overidentified ---"; ""; ""
    DW.cbps("ate", "cbps", 2, 1)

  "--- ATET (not overidentified) ---"; ""; ""
    DW.cbps("atet", "cbps", 2, 0)

  "--- ATET overidentified ---"; ""; ""
    DW.cbps("atet", "cbps", 2, 1)

// Other objective functions
  DW.cbps("atet","mean_sd_sq",1)
  DW.cbps("atet","sd_sq",1)

  DW.prognosticdiff()
  DW.reweight()
  DW.prognosticdiff()

// IPW

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations")
  stata("di _b[POmean:0.treat]")
  stata("tebalance summarize")  
  stata("tebalance summarize, baseline")  // I noticed the sum of weights in tebalance are weird
  
  DW.ipw("atet")
  table = DW.balancetable(3)
  DW.pomean()
  DW.reweight()
  table = DW.balancetable(3)

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], ate aequations")
  DW.ipw("ate")

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations tlevel(0) control(1)")
  DW.ipw("ateu")

// tradeoff between CBPS-like balance and variance in weights

  DW.cbps("atet","cbps", 1, 0)
  DW.cbps("atet","cbps", 1, 0, (1,.75,6))
  DW.cbps("atet","cbps", 1, 0, (1,.50,6))

  DW.cbps("atet","mean_sd_sq", 1, 0, (1,.75,6))
  DW.cbps("atet","mean_sd_sq", 1, 0, (1,.50,6))

end  // end of Mata block

set tracedepth 2
//set trace on
gmatch `treatvar' `varlist' if `tousevar'              , cbps
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], cbps
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], sd
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], sd cvopt(1 .6667 4)

// tradeoff between CBPS-like balance and variance in weights
gmatch `treatvar' `varlist' if `tousevar' , atet cbps treatvariance
gmatch `treatvar' `varlist' if `tousevar' , atet cbps treatvariance cvopt(1 .75 6)
gmatch `treatvar' `varlist' if `tousevar' , atet cbps treatvariance cvopt(1 .50 6)
gmatch `treatvar' `varlist' if `tousevar' , atet mean_sd_sq treatvariance cvopt(1 .75 6)
gmatch `treatvar' `varlist' if `tousevar' , atet mean_sd_sq treatvariance cvopt(1 .50 6)


// **************************
// * WEIGHTED DATA EXAMPLES *
// **************************

// Replicate CBPS
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate cbps pooledvariance
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate cbps ipw pooledvariance
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps pooledvariance
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps ipw pooledvariance

// Other objective functions
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet sd_sq treatvariance

// IPW
teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations
di _b[POmean:0.treat]
tebalance summarize
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet ipw treatvariance
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet ipw averagevariance
teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], ate aequations
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate ipw treatvariance
teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations tlevel(0) control(1)
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ateu ipw treatvariance

// tradeoff between CBPS-like balance and variance in weights
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps treatvariance
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps treatvariance cvopt(1 .75 6)
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps treatvariance cvopt(1 .50 6)
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance cvopt(1 .75 6)
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance cvopt(1 .50 6)

log close gmatch_example_ado

