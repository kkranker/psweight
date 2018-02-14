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

//****************************************************************************/
*! $Id$
*! Generalization of IPW and CBPS estimators
*! Example files
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.
//*****************************************************************************/

version 15.1
set type double
di as txt "Current user: `c(username)'" _n "Environment: `c(os)' `c(machine_type)' `: environment computername'" _n "Stata: `c(stata_version)'" cond(c(stata_version)==c(version),""," (set to version `c(version)')") _n "Date: " c(current_date) " " c(current_time)

// do C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\gmatch_one_time_setup.do

which gmatch
which gmatchcall
mata: mata describe using lgmatch

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

local depvars = "y1 y1_binary"
local treatvar = "treat"
local varlist = "x1 i.x2 i.x3 x4 x5 x6 x7 x9*"
local wgtvar = "wgt"
local tousevar = "touse"
local estimate = "atet"


// some automatic parsing based on options above, since Mata doesn't have this stuff
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

// file for testing in R
if 0 {
  fvrevar `varlist'
  export delimited `treatvar' `r(varlist)' `wgtvar' using testfile.csv if `tousevar', replace nolabel
}


// *******************************
// * RUN MODELS DIRECTLY IN MATA *
// *******************************

mata:

depvars  = st_local("depvars" )
treatvar = st_local("treatvar")
wgtvar   = st_local("wgtvar"  )
varlist  = st_local("varlist" )
tousevar = st_local("tousevar")
estimate = st_local("estimate")

// * UNWEIGHTED DATA EXAMPLES *

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
    D.balanceresults("ate",2)

  "--- ATE overidentified ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , ate over logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbps("ate", "cbps", 2, 1)
    D.balanceresults("ate",2)

  "--- ATET (not overidentified) ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , att      logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbps("atet", "cbps", 2, 0)
    D.balanceresults("atet",2)

  "--- ATET overidentified ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , att over logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbps("atet", "cbps", 2, 1)
    D.balanceresults("atet",2)

// Other objective functions
  D.cbps("atet","mean_sd_sq",1)
  D.balanceresults("atet",1)

  D.cbps("atet","sd_sq",1)
  D.balanceresults("atet",1)

// IPW
  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations")
  stata("di _b[POmean:0.treat]")
  stata("tebalance summarize")
  stata("tebalance summarize, baseline")

  D.ipw("atet")
  D.balanceresults("atet")
  table = D.balancetable(3)
  D.reweight()
  table = D.balancetable(3)

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , ate aequations")
  D.ipw("ate")
  D.balanceresults("ate")

  stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations tlevel(0) control(1)")
  D.ipw("ateu")
  D.balanceresults("ateu")

// tradeoff between CBPS-like balance and variance in weights
  D.cbps("atet","cbps", 1, 0)
  D.balanceresults("atet",1)

  D.cbps("atet","cbps", 1, 0, (1,.75,6))
  D.balanceresults("atet",1)

  D.cbps("atet","cbps", 1, 0, (1,.50,6))
  D.balanceresults("atet",1)

  D.cbps("atet","mean_sd_sq", 1, 0, (1,.75,6))
  D.balanceresults("atet",1)

  D.cbps("atet","mean_sd_sq", 1, 0, (1,.50,6))
  D.balanceresults("atet",1)

mata drop D


// * WEIGHTED DATA EXAMPLES *

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
    DW.balanceresults("ate",2)

  "--- ATE overidentified ---"; ""; ""
    DW.cbps("ate", "cbps", 2, 1)
    DW.balanceresults("ate",2)

  "--- ATET (not overidentified) ---"; ""; ""
    DW.cbps("atet", "cbps", 2, 0)
    DW.balanceresults("atet",2)

  "--- ATET overidentified ---"; ""; ""
    DW.cbps("atet", "cbps", 2, 1)
    DW.balanceresults("atet",2)

// Other objective functions
  DW.cbps("atet","mean_sd_sq",1)
  DW.balanceresults("atet",1)

  DW.cbps("atet","sd_sq",1)
  DW.balanceresults("atet",1)

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
  DW.balanceresults("atet",1)

  DW.cbps("atet","cbps", 1, 0, (1,.75,6))
  DW.balanceresults("atet",1)

  DW.cbps("atet","cbps", 1, 0, (1,.50,6))
  DW.balanceresults("atet",1)

  DW.cbps("atet","mean_sd_sq", 1, 0, (1,.75,6))
  DW.balanceresults("atet",1)

  DW.cbps("atet","mean_sd_sq", 1, 0, (1,.50,6))
  DW.balanceresults("atet",1)


end  // end of Mata block

log close gmatch_example


// *************************************
// * NOW RE-RUN WITH .ADO FILE VERSION *
// *************************************

log using gmatch_example_ado.log, name(gmatch_example_ado) replace

//****************************************************************************/
*! $Id$
*! Generalization of IPW and CBPS estimators
*! Example files
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.
//*****************************************************************************/

di as txt "Current user: `c(username)'" _n "Environment: `c(os)' `c(machine_type)' `: environment computername'" _n "Stata: `c(stata_version)'" cond(c(stata_version)==c(version),""," (set to version `c(version)')") _n "Date: " c(current_date) " " c(current_time)
desc, short
summ `treatvar' `varlist' `tousevar' `wgtvar' `depvars'


// * UNWEIGHTED DATA EXAMPLES *

// Replicate CBPS
cbps `treatvar' `varlist' if `tousevar' , ate logit optimization_technique("nr") evaluator_type("gf1")
gmatch `treatvar' `varlist' if `tousevar' , ate cbps pooledvariance
gmatchcall balanceresults()

cbps `treatvar' `varlist' if `tousevar' , ate over logit optimization_technique("nr") evaluator_type("gf1")
gmatch `treatvar' `varlist' if `tousevar' , ate cbps ipw pooledvariance
gmatchcall balanceresults()

cbps `treatvar' `varlist' if `tousevar' , att logit optimization_technique("nr") evaluator_type("gf1")
gmatch `treatvar' `varlist' if `tousevar' , atet cbps pooledvariance
gmatchcall balanceresults()

cbps `treatvar' `varlist' if `tousevar' , att over logit optimization_technique("nr") evaluator_type("gf1")
gmatch `treatvar' `varlist' if `tousevar' , atet cbps ipw pooledvariance
gmatchcall balanceresults()

// After calling gmatch, the data is stored in a class instance named gmatch_ado_most_recent
// You can print any of the public functions or variables to the screen with gmatchcall. For example:
gmatchcall diff()
gmatchcall balancetable()

// Other objective functions
gmatch `treatvar' `varlist' if `tousevar' , atet mean_sd_sq treatvariance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet sd_sq treatvariance difficult
gmatchcall balanceresults()

// IPW
teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations
di _b[POmean:0.treat]
tebalance summarize
tebalance summarize, baseline
gmatch `treatvar' `varlist' if `tousevar' , atet ipw treatvariance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet ipw averagevariance
gmatchcall balanceresults()

teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , ate aequations
gmatch `treatvar' `varlist' if `tousevar' , ate ipw treatvariance
gmatchcall balanceresults()

teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations tlevel(0) control(1)
gmatch `treatvar' `varlist' if `tousevar' , ateu ipw treatvariance
gmatchcall balanceresults()

// tradeoff between CBPS-like balance and variance in weights
gmatch `treatvar' `varlist' if `tousevar' , atet cbps treatvariance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet cbps treatvariance cvopt(1 .75 6)
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet cbps treatvariance cvopt(1 .50 6)
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet mean_sd_sq treatvariance cvopt(1 .75 6) difficult
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet mean_sd_sq treatvariance cvopt(1 .50 6) difficult
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet cbps ipw treatvariance cvopt(1 .75 6)
gmatchcall balanceresults()


// * WEIGHTED DATA EXAMPLES *

// Replicate CBPS
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate cbps pooledvariance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate cbps ipw pooledvariance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps pooledvariance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps ipw pooledvariance
gmatchcall balanceresults()

// Other objective functions
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance difficult
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet sd_sq treatvariance difficult
gmatchcall balanceresults()

// IPW
teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations
di _b[POmean:0.treat]
tebalance summarize
tebalance summarize, baseline

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet ipw treatvariance
gmatchcall balanceresults()
gmatchcall reweight()
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet ipw averagevariance
gmatchcall balanceresults()

teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], ate aequations
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate ipw treatvariance
gmatchcall balanceresults()

teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations tlevel(0) control(1)
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ateu ipw treatvariance

// tradeoff between CBPS-like balance and variance in weights
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps treatvariance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps treatvariance cvopt(1 .75 6)
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps treatvariance cvopt(1 .50 6)
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance cvopt(1 .75 6)
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance cvopt(1 .50 6)
gmatchcall balanceresults()

log close gmatch_example_ado

cap nois beep
