clear all
cls
set varabbrev off
set scheme mpr_blue
set linesize 160
set maxiter 100
cap log close gmatch_example
cap log close gmatch_example_ado
cap log close gmatch_example_R
local makegraphs = 01

cd "C:/Users/kkranker/Downloads/Ado-clone/Devel/gmatch"

do gmatch_one_time_setup.do

log using gmatch_example.log, name(gmatch_example) replace

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

D = gmatch()
D.set(treatvar, varlist, tousevar)
if (depvars!="") D.set_Y(depvars,tousevar)

// Misc balance measures
  D.means0()
  D.means1()
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
  D.gmatch("atet","mean_sd_sq",1)
  D.balanceresults("atet",1)

  D.gmatch("atet","sd_sq",1)
  D.balanceresults("atet",1)

  if (depvars!="") {
    D.gmatch("atet","stdprogdiff",1)
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
  D.means0()
  D.means1()
  D.diff()
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

  D.gmatch("atet","cbps", 1, (1,.75,6))
  D.balanceresults("atet",1)

  D.gmatch("atet","cbps", 1, (1,.50,6))
  D.balanceresults("atet",1)

  st_local("mlopts", "difficult nonrtolerance")
  D.gmatch("atet","mean_sd_sq", 1, (1,.75,6))
  D.balanceresults("atet",1)

  D.gmatch("atet","mean_sd_sq", 1, (1,.50,6))
  D.balanceresults("atet",1)

// furthermore, you can target skiwness, kurtosis, or max weight
  D.gmatch("atet","cbps", 1, (1,.50,6,1,0.7,2))
  D.balanceresults("atet",1)
  st_local("mlopts", "")

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
  DW.gmatch("atet","mean_sd_sq",1)
  DW.balanceresults("atet",1)

  DW.gmatch("atet","sd_sq",1)
  DW.balanceresults("atet",1)

  if (depvars!="") {
    DW.gmatch("atet","stdprogdiff",1)
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

  DW.gmatch("atet","cbps", 1, (1,.75,6))
  DW.balanceresults("atet",1)

  DW.gmatch("atet","cbps", 1, (1,.50,6))
  DW.balanceresults("atet",1)

  st_local("mlopts", " nonrtolerance")
  DW.gmatch("atet","mean_sd_sq", 1, (1,.75,6))
  DW.balanceresults("atet",1)

  DW.gmatch("atet","mean_sd_sq", 1, (1,.50,6))
  DW.balanceresults("atet",1)

  DW.gmatch("atet","mean_sd_sq", 1)
  DW.balanceresults("atet",1)
  st_local("mlopts", "")

  // you can also target skiwness, kurtosis, or max weight, but it's finicky

end  // end of Mata block

log close gmatch_example


// *************************************
// * NOW RE-RUN WITH .ADO FILE VERSION *
// *************************************

log using gmatch_example_ado.log, name(gmatch_example_ado) replace

di as txt "Current user: `c(username)'" _n "Environment: `c(os)' `c(machine_type)' `: environment computername'" _n "Stata: `c(stata_version)'" cond(c(stata_version)==c(version),""," (set to version `c(version)')") _n "Date: " c(current_date) " " c(current_time)
desc, short
local varlist : copy local varlist_orig
summ `treatvar' `varlist' `tousevar' `wgtvar' `depvars'


// * UNWEIGHTED DATA EXAMPLES *

// Replicate CBPS
cbps `treatvar' `varlist' if `tousevar' , ate logit optimization_technique("nr") evaluator_type("gf1")
gmatch `treatvar' `varlist' if `tousevar' , ate cbps pooledvariance
gmatchcall balanceresults()
gmatchcall means0()
return list

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
gmatch `treatvar' `varlist' if `tousevar' , atet mean_sd_sq treatvariance difficult nonrtolerance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet sd_sq treatvariance difficult nonrtolerance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet stdprogdiff depvars(`depvars') treatvariance difficult nonrtolerance
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

gmatch `treatvar' `varlist' if `tousevar' , atet cbps treatvariance cvtarget(1 .75 6)
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet cbps treatvariance cvtarget(1 .50 6)
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet mean_sd_sq treatvariance cvtarget(1 .75 6) difficult nonrtolerance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet mean_sd_sq treatvariance cvtarget(1 .50 6) difficult nonrtolerance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' , atet cbps ipw treatvariance cvtarget(1 .75 6)
gmatchcall balanceresults()

// furthermore, you can target skiwness, kurtosis, or max weight
gmatch `treatvar' `varlist' if `tousevar', atet cbps treatvariance cvtarget(1 .5 6) skewtarget(1 0.07 2)
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
gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance difficult nonrtolerance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet sd_sq treatvariance difficult nonrtolerance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet stdprogdiff depvars(`depvars') treatvariance difficult nonrtolerance
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

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps treatvariance cvtarget(1 .75 6)
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet cbps treatvariance cvtarget(1 .50 6)
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance cvtarget(1 .75 6) difficult nonrtolerance
gmatchcall balanceresults()

gmatch `treatvar' `varlist' if `tousevar' [iw=`wgtvar'], atet mean_sd_sq treatvariance cvtarget(1 .50 6) difficult nonrtolerance
gmatchcall balanceresults()

log close gmatch_example_ado


// *******************************
// * BENCHMARK A FEW MODELS IN R *
// *******************************

log using gmatch_example_R.log, name(gmatch_example_R) replace

fvrevar `varlist'
local rvarlist = r(varlist)
local rvarlist: list clean rvarlist
export delimited `treatvar' `rvarlist' `wgtvar' using "testfile.csv" if `tousevar', replace nolabel

local rvarlist: subinstr local rvarlist "__" "X__", all
local rvarlist: subinstr local rvarlist " " " + ", all
mac list _rvarlist

rsource, terminator(END_OF_R) lsource  roptions(`" --vanilla --args "`treatvar' ~ `rvarlist'" "')

  trailargs<-commandArgs(trailingOnly=TRUE);
  trailargs[1];
  getwd()
  .libPaths('c:/Users/kkranker/documents/R/myRlib');
  .libPaths();
  library(CBPS);

  mydata <- read.csv('testfile.csv', stringsAsFactors = F);
  print(summary(mydata));

  fit_ATE       <- CBPS(as.formula(trailargs[1]), data = mydata, ATT = 0, method='exact', standardize=FALSE);
  summary(fit_ATE);
  print(  fit_ATE$weights[1:10]);
  balance(fit_ATE);

  fit_ATE_over  <- CBPS(as.formula(trailargs[1]), data = mydata, ATT = 0, method='over', standardize=FALSE);
  summary(fit_ATE_over);
  print(  fit_ATE_over$weights[1:10]);
  balance(fit_ATE_over);

  fit_ATET      <- CBPS(as.formula(trailargs[1]), data = mydata, ATT = 1, method='exact', standardize=FALSE);
  summary(fit_ATET);
  print(  fit_ATET$weights[1:10]);
  balance(fit_ATET);

  fit_ATET_over <- CBPS(as.formula(trailargs[1]), data = mydata, ATT = 1, method='over', standardize=FALSE);
  summary(fit_ATET_over);
  print(  fit_ATET_over$weights[1:10]);
  balance(fit_ATET_over);

  W_fit_ATE       <- CBPS(as.formula(trailargs[1]), data = mydata, ATT = 0, method='exact', standardize=FALSE, sample.weights=mydata$wgt);
  summary(W_fit_ATE);
  print(  W_fit_ATE$weights[1:10]);
  balance(W_fit_ATE);

  W_fit_ATE_over  <- CBPS(as.formula(trailargs[1]), data = mydata, ATT = 0, method='over', standardize=FALSE, sample.weights=mydata$wgt);
  summary(W_fit_ATE_over);
  print(  W_fit_ATE_over$weights[1:10]);
  balance(W_fit_ATE_over);

  W_fit_ATET      <- CBPS(as.formula(trailargs[1]), data = mydata, ATT = 1, method='exact', standardize=FALSE, sample.weights=mydata$wgt);
  summary(W_fit_ATET);
  print(  W_fit_ATET$weights[1:10]);
  balance(W_fit_ATET);

  W_fit_ATET_over <- CBPS(as.formula(trailargs[1]), data = mydata, ATT = 1, method='over', standardize=FALSE, sample.weights=mydata$wgt);
  summary(W_fit_ATET_over);
  print(  W_fit_ATET_over$weights[1:10]);
  balance(W_fit_ATET_over);

  q();
END_OF_R


log close gmatch_example_R

cap nois beep


