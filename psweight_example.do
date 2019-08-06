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
cd "C:\Users\kkranker\Documents\Stata\psweight\code-psweight\"

do C:\Users\kkranker\Documents\Stata\psweight\code-psweight\_build.do

log using psweight_example.log, name(psweight_example) replace

// ***************************************************************************
*! psweight_example.do
*! IPW- and CBPS-type propensity score reweighting, with extentions
*! Examples
//
*! By Keith Kranker
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

if 1 {

************************************************************************************
* Simple examples in the help file (psweight.sthlp)
************************************************************************************

//  Setup
webuse cattaneo2

//  Balance before reweighting
psweight balanceonly mbsmoke mmarried mage fbaby medu, ntab

//  Estimate the average treatment effect of smoking on birthweight, using a logit model to predict treatment status
psweight ipw mbsmoke mmarried mage fbaby medu
psweight call balanceresults()

psweight call mystddiff = stddiff()
mata: mystddiff

//  Estimate the average treatment effect on the treated with CBPS
psweight cbps mbsmoke mmarried mage fbaby medu, atet
psweight call balanceresults()

//  Estimate the average treatment effect on the treated with Penalized CBPS
psweight pcbps mbsmoke mmarried mage fbaby medu, atet cvtarget(1 .5 6)
psweight call balanceresults()


************************************************************************************
* Examples in the help file (psweight_class.sthlp)
************************************************************************************


//  Setup
webuse cattaneo2, clear
gen byte touse=1

mata:

// Create an instance of the class, tell it where the data are
P = psweight()
P.st_set("mbsmoke", "mmarried mage fbaby medu", "touse")

//  Balance before reweighting
P.balancetable(2)

//  Estimate the average treatment effect of smoking on birthweight, using a logit model to predict treatment status
P.ipw()
P.balanceresults("ate", 1)

//  Estimate the average treatment effect on the treated with CBPS
P.cbps("atet")
P.balanceresults("atet", 1)

//  Estimate the average treatment effect on the treated with Penalized CBPS
P.solve("atet", "cbps", 2, (1, .5, 6))
P.balanceresults("atet", 1)

end // end of Mata

} // end of simple examples


************************************************************************************
* Describe/summarize the example datasets
************************************************************************************

*** Input data file (simple_cattaneo_data) comes from the program named Make_example_datasets.do (in C:\Users\kkranker\Documents\Stata\Multiple-Equation-Models)

use "C:\Users\kkranker\Documents\Stata\Multiple-Equation-Models\simple_cattaneo_data.dta", clear
desc, short
notes _dta
summ, sep(0)
tab treat treat_cat, mi
corr treat y1 y1_binary

local if if _n<=500
set seed 1
gen wgt = max(.1, rnormal(2,.4))
gen fwgt = round(rnormal(2,.4))
// forvalues i = 20/200 {
forvalues i = 90/95 {
  gen x`i' = rnormal()
}

local depvarlist = "y1 y1_binary"
local treatvar = "treat"
local varlist = "x1 i.x2 i.x3 x4 x5 x6 x7 x9*"
local wgtvar = "wgt"
local tousevar = "touse"
local estimate = "atet"


// some automatic parsing based on options above, since Mata doesn't have this stuff
if "`wgtvar'"!="" local wgtexp "[iw=`wgtvar']"
mark    `tousevar' `if' `in' `wgtexp'
markout `tousevar' `depvarlist' `treatvar' `varlist'
local varlist_orig : copy local varlist

_rmcoll `treatvar' `varlist' if `tousevar' `wgtexp', expand logit touse(`tousevar')
local varlist `r(varlist)'
gettoken trash varlist : varlist


// *******************************
// * RUN MODELS DIRECTLY IN MATA *
// *******************************

if 1 {

mata:

depvarlist = st_local("depvarlist")
treatvar   = st_local("treatvar")
wgtvar     = st_local("wgtvar")
varlist    = st_local("varlist")
tousevar   = st_local("tousevar")
estimate   = st_local("estimateerror")
// * UNWEIGHTED DATA EXAMPLES *

D = psweight()
D.st_set(treatvar, varlist, tousevar)
if (depvarlist!="") D.st_set_depvars(depvarlist, tousevar)

// Misc balance measures
  D.diff()
  D.stddiff()
  D.stddiff(1)
  D.mean_asd(1)
  D.stddiff(1)
  D.stddiff(0)
  D.varratio()
  D.progdiff(1)
  D.get_N(1)

  "Balance table before matching"
  table = D.balancetable(1)

// Replicate CBPS

  "--- ATE (not overidentified) ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , ate      logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbps("ate", 2)
    D.balanceresults("ate", 2)
    D.get_N(1)

  "--- ATE overidentified ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , ate over logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbpsoid("ate", 2)
    D.balanceresults("ate", 2)

  "--- ATET (not overidentified) ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , att      logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbps("atet", 2)
    D.balanceresults("atet", 2)

  "--- ATET overidentified ---"; ""; ""
    stata(`"cbps `treatvar' `varlist' if `tousevar' , att over logit optimization_technique("nr") evaluator_type("gf1")"')
    stata(`"cbps_imbalance"')
    D.cbpsoid("atet", 2)
    D.balanceresults("atet", 2)

// Other objective functions
  st_local("mlopts", "difficult nonrtolerance")
  D.solve("atet","mean_sd_sq", 1)
  D.balanceresults("atet", 1)

  D.solve("atet","sd_sq", 1)
  D.balanceresults("atet", 1)

  if (depvarlist!="") {
    D.solve("atet","stdprogdiff", 1)
    D.balanceresults("atet", 1)
  }
  st_local("mlopts", "")

// IPW
  stata("teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' , atet aequations")
  stata("di _b[POmean:0.treat]")
  stata("tebalance summarize")
  stata("tebalance summarize, baseline")

  D.ipw("atet")
  D.balanceresults("atet", 1)
  table = D.balancetable(3)
  D.reweight()
  table = D.balancetable(3)

  stata("teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' , ate aequations")
  D.ipw("ate")
  D.balanceresults("ate", 1)

  stata("teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' , atet aequations tlevel(0) control(1)")
  D.ipw("ateu")
  D.balanceresults("ateu", 1)

// tradeoff between CBPS-like balance and variance in weights
  D.cbps("atet", 1)
  D.balanceresults("atet", 1)

  D.solve("atet","cbps", 1, (1,.75, 6))
  D.balanceresults("atet", 1)

  D.solve("atet","cbps", 1, (1,.50, 6))
  D.balanceresults("atet", 1)

  st_local("mlopts", "difficult nonrtolerance")
  D.solve("atet","mean_sd_sq", 1, (1,.75, 6))
  D.balanceresults("atet", 1)

  D.solve("atet","mean_sd_sq", 1, (1,.50, 6))
  D.balanceresults("atet", 1)

// furthermore, you can target skiwness, kurtosis, or max weight
  D.solve("atet","cbps", 1, (1,.50, 6, 1, 0.7, 2))
  D.balanceresults("atet", 1)
  st_local("mlopts", "")

Dcpy = psweight()
Dcpy.clone(D)
Dcpy.balancetable(2)
mata drop Dcpy
mata drop D


// * WEIGHTED DATA EXAMPLES *

DW = psweight()
DW.st_set(treatvar, varlist, tousevar, wgtvar)
if (depvarlist!="") DW.st_set_depvars(depvarlist, tousevar)

// Misc balance measures
  DW.diff()
  DW.stddiff()
  DW.stddiff(1)
  DW.mean_asd(1)
  DW.stddiff(1)
  DW.stddiff(0)
  DW.varratio()
  DW.progdiff(1)
  DW.get_N(1)

  "Balance table before matching"
  temp = DW.balancetable(1)

// Replicate CBPS

  "--- ATE (not overidentified) ---"; ""; ""
    DW.cbps("ate", 2)
    DW.balanceresults("ate", 2)
    DW.get_N(1)

  "--- ATE overidentified ---"; ""; ""
    DW.cbpsoid("ate", 2)
    DW.balanceresults("ate", 2)

  "--- ATET (not overidentified) ---"; ""; ""
    DW.cbps("atet", 2)
    DW.balanceresults("atet", 2)

  "--- ATET overidentified ---"; ""; ""
    DW.cbpsoid("atet", 2)
    DW.balanceresults("atet", 2)

// Other objective functions
  st_local("mlopts", "difficult nonrtolerance")
  DW.solve("atet","mean_sd_sq", 1)
  DW.balanceresults("atet", 1)

  DW.solve("atet","sd_sq", 1)
  DW.balanceresults("atet", 1)

  if (depvarlist!="") {
    DW.solve("atet","stdprogdiff", 1)
    DW.balanceresults("atet", 1)
  }
  st_local("mlopts", "")

  DW.progdiff(1)
  DW.reweight()
  DW.progdiff(1)

// IPW

  stata("teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations")
  stata("di _b[POmean:0.treat]")
  stata("tebalance summarize")
  stata("tebalance summarize, baseline")  // I noticed the sum of weights in tebalance are weird

  DW.ipw("atet")
  table = DW.balancetable(3)
  DW.pomean()
  DW.reweight()
  table = DW.balancetable(3)

  stata("teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], ate aequations")
  DW.ipw("ate")

  stata("teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations tlevel(0) control(1)")
  DW.ipw("ateu")

// tradeoff between CBPS-like balance and variance in weights

  DW.cbps("atet", 1)
  DW.balanceresults("atet", 1)

  DW.solve("atet","cbps", 1, (1,.75, 6))
  DW.balanceresults("atet", 1)

  DW.solve("atet","cbps", 1, (1,.50, 6))
  DW.balanceresults("atet", 1)

  st_local("mlopts", " nonrtolerance")
  DW.solve("atet","mean_sd_sq", 1, (1,.75, 6))
  DW.balanceresults("atet", 1)

  DW.solve("atet","mean_sd_sq", 1, (1,.50, 6))
  DW.balanceresults("atet", 1)

  DW.solve("atet","mean_sd_sq", 1)
  DW.balanceresults("atet", 1)
  st_local("mlopts", "")

  // you can also target skiwness, kurtosis, or max weight, but it's finicky

end  // end of Mata block

} // end of Mata examples

log off psweight_example


// *************************************
// * NOW RE-RUN WITH .ADO FILE VERSION *
// *************************************

if 1 {

log using psweight_example_ado.log, name(psweight_example_ado) replace

di as txt "Current user: `c(username)'" _n "Environment: `c(os)' `c(machine_type)' `: environment computername'" _n "Stata: `c(stata_version)'" cond(c(stata_version)==c(version),""," (set to version `c(version)')") _n "Date: " c(current_date) " " c(current_time)
desc, short
local varlist : copy local varlist_orig
summ `treatvar' `varlist' `tousevar' `wgtvar' `depvarlist'
if (trim("`depvarlist'") != "") local depvaropt depvarlist(`depvarlist')
ereturn clear
 return clear

// balance before matching
psweight balanceonly `treatvar' `varlist' if `tousevar', ntable
ereturn list
 return list

psweight balanceonly `treatvar' `varlist' if `tousevar' [iw=`wgtvar']
psweight balanceonly `treatvar' `varlist' if `tousevar', mweight(`wgtvar') ate

// * UNWEIGHTED DATA EXAMPLES *

// Replicate CBPS
cbps `treatvar' `varlist' if `tousevar' , ate logit optimization_technique("nr") evaluator_type("gf1")
ereturn clear
 return clear

psweight cbps `treatvar' `varlist' if `tousevar' , ate pooledvariance `depvaropt' ntable
ereturn list
 return list

psweight call balanceresults()
ereturn list
return list
psweight // test replay

cbps `treatvar' `varlist' if `tousevar' , ate over logit optimization_technique("nr") evaluator_type("gf1")
psweight cbpsoid `treatvar' `varlist' if `tousevar' , ate pooledvariance  `depvaropt'
psweight call balanceresults()

cbps `treatvar' `varlist' if `tousevar' , att logit optimization_technique("nr") evaluator_type("gf1")
psweight cbps `treatvar' `varlist' if `tousevar' , atet pooledvariance  `depvaropt'
psweight call balanceresults()

cbps `treatvar' `varlist' if `tousevar' , att over logit optimization_technique("nr") evaluator_type("gf1")
psweight cbpsoid `treatvar' `varlist' if `tousevar' , atet pooledvariance  `depvaropt'
psweight call balanceresults()

// After calling psweight, the data is stored in a class instance named psweight_ado_most_recent
// You can print any of the public functions or variables to the screen with psweight call. For example:
psweight call diff()
psweight call balancetable()

// Other objective functions
psweight mean_sd_sq `treatvar' `varlist' if `tousevar' , atet treatvariance difficult nonrtolerance `depvaropt'
psweight call balanceresults()

psweight sd_sq `treatvar' `varlist' if `tousevar' , atet treatvariance difficult nonrtolerance `depvaropt'
psweight call balanceresults()

psweight stdprogdiff `treatvar' `varlist' if `tousevar' , atet treatvariance difficult nonrtolerance `depvaropt'
psweight call balanceresults()


// IPW
teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' , atet aequations
di _b[POmean:0.treat]
tebalance summarize
tebalance summarize, baseline
psweight ipw `treatvar' `varlist' if `tousevar' , atet averagevariance `depvaropt'
psweight call balanceresults()

psweight ipw `treatvar' `varlist' if `tousevar' , atet averagevariance `depvaropt'
psweight call balanceresults()

teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' , ate aequations
psweight ipw `treatvar' `varlist' if `tousevar' , ate pooledvariance
psweight call balanceresults()

teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' , atet aequations tlevel(0) control(1)
psweight ipw `treatvar' `varlist' if `tousevar' , ateu pooledvariance
psweight call balanceresults()

// tradeoff between CBPS-like balance and variance in weights
psweight cbps `treatvar' `varlist' if `tousevar' , atet treatvariance `depvaropt'
psweight call balanceresults()

psweight pcbps `treatvar' `varlist' if `tousevar' , atet treatvariance cvtarget(1 .75 6) `depvaropt'
psweight call balanceresults()

psweight cbps `treatvar' `varlist' if `tousevar' , atet treatvariance cvtarget(1 .50 6) `depvaropt' // could say pcbps or cbps
psweight call balanceresults()

psweight mean_sd_sq `treatvar' `varlist' if `tousevar' , atet treatvariance cvtarget(1 .75 6) difficult nonrtolerance `depvaropt'
psweight call balanceresults()

psweight mean_sd_sq `treatvar' `varlist' if `tousevar' , atet treatvariance cvtarget(1 .50 6) difficult nonrtolerance `depvaropt'
psweight call balanceresults()

psweight cbpsoid `treatvar' `varlist' if `tousevar' , atet treatvariance cvtarget(1 .75 6) `depvaropt'
psweight call balanceresults()

// furthermore, you can target skiwness, kurtosis, or max weight
psweight pcbps `treatvar' `varlist' if `tousevar', atet treatvariance cvtarget(1 .5 6) skewtarget(1 0.07 2) `depvaropt'
psweight call balanceresults()

capture nois {
psweight pcbps `treatvar' `varlist' if `tousevar', atet treatvariance maxtarget(1 10 2) difficult nonrtolerance `depvaropt'
psweight call balanceresults()
}

// * WEIGHTED DATA EXAMPLES *

// Replicate CBPS
psweight cbps `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate pooledvariance `depvaropt'
psweight call balanceresults()

psweight cbpsoid `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate pooledvariance `depvaropt'
psweight call balanceresults()

psweight cbps `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet pooledvariance `depvaropt'
psweight call balanceresults()

psweight cbpsoid `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet pooledvariance `depvaropt'
psweight call balanceresults()

// Other objective functions
psweight mean_sd_sq `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet treatvariance difficult nonrtolerance `depvaropt'
psweight call balanceresults()

psweight sd_sq `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet treatvariance difficult nonrtolerance `depvaropt'
psweight call balanceresults()

psweight stdprogdiff `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet treatvariance difficult nonrtolerance `depvaropt'
psweight call balanceresults()

// IPW
teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations
di _b[POmean:0.treat]
tebalance summarize
tebalance summarize, baseline

psweight ipw `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet treatvariance `depvaropt' ntable
psweight call balanceresults()
psweight call reweight()
psweight call balanceresults()

psweight ipw `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet averagevariance `depvaropt'
psweight call balanceresults()

teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], ate aequations
psweight ipw `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ate treatvariance
psweight call balanceresults()

teffects ipw (`:word 1 of `depvarlist'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations tlevel(0) control(1)
psweight ipw `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], ateu treatvariance `depvaropt'

// tradeoff between CBPS-like balance and variance in weights
psweight cbps `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet  treatvariance `depvaropt'
psweight call balanceresults()

psweight pcbps `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet treatvariance cvtarget(1 .75 6) `depvaropt'
psweight call balanceresults()

psweight pcbps `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet treatvariance cvtarget(1 .50 6) `depvaropt'
psweight call balanceresults()

psweight mean_sd_sq `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet treatvariance cvtarget(1 .75 6) difficult nonrtolerance `depvaropt'
psweight call balanceresults()

psweight mean_sd_sq `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet treatvariance cvtarget(1 .50 6) difficult nonrtolerance `depvaropt'
psweight call balanceresults()

// you can also target skiwness, kurtosis, or max weight, but it's finicky

capture nois {
psweight mean_sd_sq `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet treatvariance maxtarget(1 10 2) difficult nonrtolerance `depvaropt'
psweight call balanceresults()
}

log close psweight_example_ado

} // end of ado examples

// *******************************
// * BENCHMARK A FEW MODELS IN R *
// *******************************

if 0 {

log using psweight_example_R.log, name(psweight_example_R) replace

fvrevar `varlist'
tempfile csvout
export delimited `treatvar' `r(varlist)' `wgtvar' using "C:\Users\kkranker\Documents\Stata\psweight\code-psweight\testfile.csv" if `tousevar', replace nolabel

rsource, terminator(END_OF_R) lsource
  mydata <- read.csv("C:\\Users\\kkranker\\Documents\\Stata\\psweight\\code-psweight\\testfile.csv", stringsAsFactors = F);
  library(CBPS);
  summary(mydata);
  fit_ATE         <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 0, method='exact', standardize=TRUE);
  summary(fit_ATE);
  print(  fit_ATE$weights[1:10]);
  balance(fit_ATE);
  fit_ATE_over    <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 0, method="over", standardize=TRUE);
  summary(fit_ATE_over);
  print(  fit_ATE_over$weights[1:10]);
  balance(fit_ATE_over);
  fit_ATET        <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 1, method="exact", standardize=TRUE);
  summary(fit_ATET);
  print(  fit_ATET$weights[1:10]);
  balance(fit_ATET);
  fit_ATET_over   <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 1, method="over", standardize=TRUE);
  summary(fit_ATET_over);
  print(  fit_ATET_over$weights[1:10]);
  balance(fit_ATET_over);
  W_fit_ATE       <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 0, method='exact', standardize=TRUE, sample.weights=mydata$wgt);
  summary(W_fit_ATE);
  print(  W_fit_ATE$weights[1:10]);
  balance(W_fit_ATE);
  W_fit_ATE_over  <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 0, method="over", standardize=TRUE, sample.weights=mydata$wgt);
  summary(W_fit_ATE_over);
  print(  W_fit_ATE_over$weights[1:10]);
  balance(W_fit_ATE_over);
  W_fit_ATET      <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 1, method="exact", standardize=TRUE, sample.weights=mydata$wgt);
  summary(W_fit_ATET);
  print(  W_fit_ATET$weights[1:10]);
  balance(W_fit_ATET);
  W_fit_ATET_over <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 1, method="over", standardize=TRUE, sample.weights=mydata$wgt);
  summary(W_fit_ATET_over);
  print(  W_fit_ATET_over$weights[1:10]);
  balance(W_fit_ATET_over);
  q();
END_OF_R
erase "C:\Users\kkranker\Documents\Stata\psweight\code-psweight\testfile.csv"

log close psweight_example_R

} // end of R benchmarks

log on psweight_example
log close psweight_example
