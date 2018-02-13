mac drop _all
clear all
cls
set varabbrev off
set scheme mpr_blue
set linesize 160
set maxiter 100
cap log close gmatch_example
local makegraphs = 01
cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\"

log using gmatch_example.log, name(gmatch_example) replace

// Multiple-equation models: An introduction and potential applications to our work at Mathematica
// Design and methods �brown bag� workshop
// April 8, 2015
// Keith Kranker

// Stata_code_2_IPW.do This program includes all the examples in the powerpoint slides, plus more.
// This program includes examples of how to code inverse propensity weighting (IPW) estimators using Stata's GMM command

include C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\gmatchclass.mata

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
D.set( st_local("treatvar"),st_local("varlist"), st_local("tousevar"))
if (depvars!="") D.set_Y(st_local("depvars"),st_local("tousevar"))

M = gmatch()
M.clone(D)

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
    cbpsweight = M.cbps("ate", "cbps", 2, 0)

  "--- ATE overidentified ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , ate over logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    cbpsweight = M.cbps("ate", "cbps", 2, 1)

  "--- ATET (not overidentified) ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , att      logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    cbpsweight = M.cbps("atet", "cbps", 2, 0)

  "--- ATET overidentified ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , att over logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    cbpsweight = M.cbps("atet", "cbps", 2, 1)

// Other objective functions
  cbpsweight = M.cbps("atet","mean_sd_sq",1)
  cbpsweight = M.cbps("atet","sd_sq",1)
  cbpsweight = M.cbps("atet","mean_asd") // I took these out; should error
  cbpsweight = M.cbps("atet","max_asd")  // I took these out; should error  
//  cbpsweight = M.cbps("atet","mean_sd_sq_cv",1, (1,1,6)

// IPW
  M = gmatch()
  M.clone(D)
  stata("qui teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations")
  iwpweight = D.ipw("atet")

  M.reweight(iwpweight)

  stata("di _b[POmean:0.treat]")
  M.pomean()

  stata("tebalance summarize")
  table = M.balancetable(3)

  stata("tebalance summarize, baseline")
  M.reweight()
  table = M.balancetable(3)

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations")
  stata("predict pscore1, tlevel(1) ")
  stata("list `treatvar' pscore1 in 1/20, nolab ")
  ipw = D.ipw("atet")

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , ate aequations")
  ipw = D.ipw("ate")

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations tlevel(0) control(1)")
  ipw = D.ipw("ateu")

// tradeoff between CBPS-like balance and variance in weights
  cbpsweight = M.cbps("atet","cbps", 1, 0)
  cbpsweight = M.cbps("atet","cbps", 1, 0, (1,.75,6))
  cbpsweight = M.cbps("atet","cbps", 1, 0, (1,.50,6))
  
  cbpsweight = M.cbps("atet","mean_sd_sq", 1, 0, (1,.75,6))
  cbpsweight = M.cbps("atet","mean_sd_sq", 1, 0, (1,.50,6))
   
mata drop D M


// **************************
// * WEIGHTED DATA EXAMPLES *
// **************************

DW = gmatch()
DW.set(st_local("treatvar"),st_local("varlist"), st_local("tousevar"), st_local("wgtvar"))
if (depvars!="") DW.set_Y(st_local("depvars"),st_local("tousevar"))

MW = gmatch()
MW.clone(DW)

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
    cbpsweight = MW.cbps("ate", "cbps", 2, 0)

  "--- ATE overidentified ---"; ""; ""
    cbpsweight = MW.cbps("ate", "cbps", 2, 1)

  "--- ATET (not overidentified) ---"; ""; ""
    cbpsweight = MW.cbps("atet", "cbps", 2, 0)

  "--- ATET overidentified ---"; ""; ""
    cbpsweight = MW.cbps("atet", "cbps", 2, 1)
  
// Other objective functions
  cbpsweight = MW.cbps("atet","mean_sd_sq",1)
  cbpsweight = MW.cbps("atet","sd_sq",1)

  MW.prognosticdiff()
  MW.reweight(cbpsweight)
  MW.prognosticdiff()
  
// IPW

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations")
  stata("tebalance summarize, baseline")  // I noticed the sum of weights in tebalance are weird
  table = DW.balancetable(3)

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations")
  ipw = DW.ipw("atet")
  stata("tebalance summarize")  // I noticed the sum of weights in tebalance are weird
  DW.reweight(ipw)
  table = DW.balancetable(3)

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], ate aequations")
  DW.reweight()
  ipw = DW.ipw("ate")

// tradeoff between CBPS-like balance and variance in weights

  mata drop MW
  MW = gmatch()
  MW.clone(DW)
  
  cbpsweight = MW.cbps("atet","cbps", 1, 0)
  cbpsweight = MW.cbps("atet","cbps", 1, 0, (1,.75,6))
  cbpsweight = MW.cbps("atet","cbps", 1, 0, (1,.50,6))
  
  cbpsweight = MW.cbps("atet","mean_sd_sq", 1, 0, (1,.75,6))
  cbpsweight = MW.cbps("atet","mean_sd_sq", 1, 0, (1,.50,6))
   
  

end  // end of Mata block

log close gmatch_example
