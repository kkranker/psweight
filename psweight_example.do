clear all
cls
set varabbrev off
set scheme mpr_blue
set linesize 160
set maxiter 100
cap log close psweight_example
cap log close psweight_example_ado
cap log close psweight_example_R
local makegraphs = 01
cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\"

do C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\psweight_one_time_setup.do

log using psweight_example.log, name(psweight_example) replace

// ***************************************************************************
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
// ***************************************************************************

version 15.1
set type double
di as txt "Current user: `c(username)'" _n "Environment: `c(os)' `c(machine_type)' `: environment computername'" _n "Stata: `c(stata_version)'" cond(c(stata_version)==c(version),""," (set to version `c(version)')") _n "Date: " c(current_date) " " c(current_time)

which psweight
mata: mata describe using lpsweight

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
local varlist_orig : copy local varlist

_rmcoll `treatvar' `varlist' if `tousevar' `wgtexp', expand logit touse(`tousevar')
local varlist `r(varlist)'
gettoken trash varlist : varlist


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

D = psweight()
D.set(treatvar, varlist, tousevar)
if (depvars!="") D.set_Y(depvars,tousevar)

// Misc balance measures
  D.diff()
  D.stddiff()
  D.mean_asd()
  D.stddiff(1)
  D.stddiff(0)
  D.varratio()
  D.progdiff()

  "Balance table before matching"
  table = D.balancetable(1)

// Replicate CBPS

  "--- ATE (not overidentified) ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , ate      logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbps("ate", 2)
    D.balanceresults("ate",2)

  "--- ATE overidentified ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , ate over logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbpsoid("ate", 2)
    D.balanceresults("ate",2)

  "--- ATET (not overidentified) ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , att      logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbps("atet", 2)
    D.balanceresults("atet",2)

  "--- ATET overidentified ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , att over logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbpsoid("atet", 2)
    D.balanceresults("atet",2)

// Other objective functions
  st_local("mlopts", "difficult nonrtolerance")
  D.psweight("atet","mean_sd_sq",1)
  D.balanceresults("atet",1)

  D.psweight("atet","sd_sq",1)
  D.balanceresults("atet",1)

  if (depvars!="") {
    D.psweight("atet","stdprogdiff",1)
    D.balanceresults("atet",1)
  }
  st_local("mlopts", "")

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
  D.cbps("atet", 1)
  D.balanceresults("atet",1)

  D.psweight("atet","cbps", 1, (1,.75,6))
  D.balanceresults("atet",1)

  D.psweight("atet","cbps", 1, (1,.50,6))
  D.balanceresults("atet",1)

  st_local("mlopts", "difficult nonrtolerance")
  D.psweight("atet","mean_sd_sq", 1, (1,.75,6))
  D.balanceresults("atet",1)

  D.psweight("atet","mean_sd_sq", 1, (1,.50,6))
  D.balanceresults("atet",1)

// furthermore, you can target skiwness, kurtosis, or max weight
  D.psweight("atet","cbps", 1, (1,.50,6,1,0.7,2))
  D.balanceresults("atet",1)
  st_local("mlopts", "")

mata drop D


// * WEIGHTED DATA EXAMPLES *

DW = psweight()
DW.set(treatvar, varlist, tousevar, wgtvar)
if (depvars!="") DW.set_Y(depvars,tousevar)

// Misc balance measures
  DW.diff()
  DW.stddiff()
  DW.mean_asd()
  DW.stddiff(1)
  DW.stddiff(0)
  DW.varratio()
  DW.progdiff()

  "Balance table before matching"
  temp = DW.balancetable(1)

// Replicate CBPS

  "--- ATE (not overidentified) ---"; ""; ""
    DW.cbps("ate", 2)
    DW.balanceresults("ate",2)

  "--- ATE overidentified ---"; ""; ""
    DW.cbpsoid("ate", 2)
    DW.balanceresults("ate",2)

  "--- ATET (not overidentified) ---"; ""; ""
    DW.cbps("atet", 2)
    DW.balanceresults("atet",2)

  "--- ATET overidentified ---"; ""; ""
    DW.cbpsoid("atet", 2)
    DW.balanceresults("atet",2)

// Other objective functions
  st_local("mlopts", "difficult nonrtolerance")
  DW.psweight("atet","mean_sd_sq",1)
  DW.balanceresults("atet",1)

  DW.psweight("atet","sd_sq",1)
  DW.balanceresults("atet",1)

  if (depvars!="") {
    DW.psweight("atet","stdprogdiff",1)
    DW.balanceresults("atet",1)
  }
  st_local("mlopts", "")

  DW.progdiff()
  DW.reweight()
  DW.progdiff()

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

  DW.cbps("atet", 1)
  DW.balanceresults("atet",1)

  DW.psweight("atet","cbps", 1, (1,.75,6))
  DW.balanceresults("atet",1)

  DW.psweight("atet","cbps", 1, (1,.50,6))
  DW.balanceresults("atet",1)

  st_local("mlopts", " nonrtolerance")
  DW.psweight("atet","mean_sd_sq", 1, (1,.75,6))
  DW.balanceresults("atet",1)

  DW.psweight("atet","mean_sd_sq", 1, (1,.50,6))
  DW.balanceresults("atet",1)

  DW.psweight("atet","mean_sd_sq", 1)
  DW.balanceresults("atet",1)
  st_local("mlopts", "")

  // you can also target skiwness, kurtosis, or max weight, but it's finicky

end  // end of Mata block

log close psweight_example


// *************************************
// * NOW RE-RUN WITH .ADO FILE VERSION *
// *************************************

log using psweight_example_ado.log, name(psweight_example_ado) replace

di as txt "Current user: `c(username)'" _n "Environment: `c(os)' `c(machine_type)' `: environment computername'" _n "Stata: `c(stata_version)'" cond(c(stata_version)==c(version),""," (set to version `c(version)')") _n "Date: " c(current_date) " " c(current_time)
desc, short
local varlist : copy local varlist_orig
summ `treatvar' `varlist' `tousevar' `wgtvar' `depvars'


// balance before matching
psweight `treatvar' `varlist' if `tousevar' , balanceonly

// * UNWEIGHTED DATA EXAMPLES *

// Replicate CBPS
cbps `treatvar' `varlist' if `tousevar' , ate logit optimization_technique("nr") evaluator_type("gf1")
psweight `treatvar' `varlist' if `tousevar' , ate cbps pooledvariance
psweight call balanceresults()
return list

cbps `treatvar' `varlist' if `tousevar' , ate over logit optimization_technique("nr") evaluator_type("gf1")
psweight `treatvar' `varlist' if `tousevar' , ate cbps ipw pooledvariance
psweight call balanceresults()

cbps `treatvar' `varlist' if `tousevar' , att logit optimization_technique("nr") evaluator_type("gf1")
psweight `treatvar' `varlist' if `tousevar' , atet cbps pooledvariance
psweight call balanceresults()

cbps `treatvar' `varlist' if `tousevar' , att over logit optimization_technique("nr") evaluator_type("gf1")
psweight `treatvar' `varlist' if `tousevar' , atet cbps ipw pooledvariance
psweight call balanceresults()

// After calling psweight, the data is stored in a class instance named psweight_ado_most_recent
// You can print any of the public functions or variables to the screen with psweight call. For example:
psweight call diff()
psweight call balancetable()

// Other objective functions
psweight `treatvar' `varlist' if `tousevar' , atet mean_sd_sq treatvariance difficult nonrtolerance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' , atet sd_sq treatvariance difficult nonrtolerance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' , atet stdprogdiff depvars(`depvars') treatvariance difficult nonrtolerance
psweight call balanceresults()


// IPW
teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations
di _b[POmean:0.treat]
tebalance summarize
tebalance summarize, baseline
psweight `treatvar' `varlist' if `tousevar' , atet ipw treatvariance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' , atet ipw averagevariance
psweight call balanceresults()

teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , ate aequations
psweight `treatvar' `varlist' if `tousevar' , ate ipw treatvariance
psweight call balanceresults()

teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations tlevel(0) control(1)
psweight `treatvar' `varlist' if `tousevar' , ateu ipw treatvariance
psweight call balanceresults()

// tradeoff between CBPS-like balance and variance in weights
psweight `treatvar' `varlist' if `tousevar' , atet cbps treatvariance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' , atet cbps treatvariance cvtarget(1 .75 6)
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' , atet cbps treatvariance cvtarget(1 .50 6)
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' , atet mean_sd_sq treatvariance cvtarget(1 .75 6) difficult nonrtolerance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' , atet mean_sd_sq treatvariance cvtarget(1 .50 6) difficult nonrtolerance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' , atet cbps ipw treatvariance cvtarget(1 .75 6)
psweight call balanceresults()

// furthermore, you can target skiwness, kurtosis, or max weight
psweight `treatvar' `varlist' if `tousevar', atet cbps treatvariance cvtarget(1 .5 6) skewtarget(1 0.07 2)
psweight call balanceresults()

capture nois {
psweight `treatvar' `varlist' if `tousevar', atet cbps treatvariance maxtarget(1 10 2) difficult nonrtolerance
psweight call balanceresults()
}

// * WEIGHTED DATA EXAMPLES *

// Replicate CBPS
psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate cbps pooledvariance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate cbps ipw pooledvariance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps pooledvariance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps ipw pooledvariance
psweight call balanceresults()

// Other objective functions
psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance difficult nonrtolerance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet sd_sq treatvariance difficult nonrtolerance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet stdprogdiff depvars(`depvars') treatvariance difficult nonrtolerance
psweight call balanceresults()

// IPW
teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations
di _b[POmean:0.treat]
tebalance summarize
tebalance summarize, baseline

psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet ipw treatvariance
psweight call balanceresults()
psweight call reweight()
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet ipw averagevariance
psweight call balanceresults()

teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], ate aequations
psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate ipw treatvariance
psweight call balanceresults()

teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations tlevel(0) control(1)
psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ateu ipw treatvariance

// tradeoff between CBPS-like balance and variance in weights
psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps treatvariance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps treatvariance cvtarget(1 .75 6)
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps treatvariance cvtarget(1 .50 6)
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance cvtarget(1 .75 6) difficult nonrtolerance
psweight call balanceresults()

psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance cvtarget(1 .50 6) difficult nonrtolerance
psweight call balanceresults()

// you can also target skiwness, kurtosis, or max weight, but it's finicky

capture nois {
psweight `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance maxtarget(1 10 2) difficult nonrtolerance
psweight call balanceresults()
}

log close psweight_example_ado


// *******************************
// * BENCHMARK A FEW MODELS IN R *
// *******************************

log using psweight_example_R.log, name(psweight_example_R) replace

fvrevar `varlist'
export delimited `treatvar' `r(varlist)' `wgtvar' using "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\testfile.csv" if `tousevar', replace nolabel

rsource, terminator(END_OF_R)
  mydata <- read.csv("C:\\Users\\kkranker\\Documents\\Stata\\Ado\\Devel\\psweight\\testfile.csv", stringsAsFactors = F);
  library(CBPS);
  fit_ATE       <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 0, method='exact', standardize=stdall);
  summary(fit_ATE);
  print(  fit_ATE$weights[1:10]);
  balance(fit_ATE);
  fit_ATE_over  <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 0, method="over", standardize=stdall);
  summary(fit_ATE_over);
  print(  fit_ATE_over$weights[1:10]);
  balance(fit_ATE_over);
  fit_ATET      <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 1, method="exact", standardize=stdall);
  summary(fit_ATET);
  print(  fit_ATET$weights[1:10]);
  balance(fit_ATET);
  fit_ATET_over <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 1, method="over", standardize=stdall);
  summary(fit_ATET_over);
  print(  fit_ATET_over$weights[1:10]);
  balance(fit_ATET_over);
  W_fit_ATE       <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 0, method='exact', standardize=stdall, sample.weights=mydata$wgt);
  summary(W_fit_ATE);
  print(  W_fit_ATE$weights[1:10]);
  balance(W_fit_ATE);
  W_fit_ATE_over  <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 0, method="over", standardize=stdall, sample.weights=mydata$wgt);
  summary(W_fit_ATE_over);
  print(  W_fit_ATE_over$weights[1:10]);
  balance(W_fit_ATE_over);
  W_fit_ATET      <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 1, method="exact", standardize=stdall, sample.weights=mydata$wgt);
  summary(W_fit_ATET);
  print(  W_fit_ATET$weights[1:10]);
  balance(W_fit_ATET);
  W_fit_ATET_over <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 1, method="over", standardize=stdall, sample.weights=mydata$wgt);
  summary(W_fit_ATET_over);
  print(  W_fit_ATET_over$weights[1:10]);
  balance(W_fit_ATET_over);
  q();
END_OF_R
erase "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\testfile.csv"


log close psweight_example_R

cap nois beep
