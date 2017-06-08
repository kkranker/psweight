clear all
mac drop _all
cls
set varabbrev off
set scheme mpr_blue
set linesize 160
cap log close _all
local makegraphs = 01
cd "C:\Users\kkranker\Documents\Stata\Multiple-Equation-Models\"


// Multiple-equation models: An introduction and potential applications to our work at Mathematica
// Design and methods “brown bag” workshop
// April 8, 2015
// Keith Kranker

// Stata_code_2_IPW.do This program includes all the examples in the powerpoint slides, plus more.
// This program includes examples of how to code inverse propensity weighting (IPW) estimators using Stata's GMM command

version 14.1
set type double
di as txt "Current user: `c(username)'" _n "Environment: `c(os)' `c(machine_type)' `: environment computername'" _n "Stata: `c(stata_version)'" cond(c(stata_version)==c(version),""," (set to version `c(version)')") _n "Date: " c(current_date) " " c(current_time)



************************************************************************************
* Describe/summarize the example datasets
************************************************************************************

*** Input data file (simple_cattaneo_data) comes from the program named Make_example_datasets.do (in C:\Users\kkranker\Documents\Stata\Multiple-Equation-Models)

use simple_cattaneo_data
desc, short
notes _dta
summ, sep(0)
tab treat treat_cat, mi
corr treat y1 y1_binary

/*

************************************************************************************
* 5. IPW w/ binary treatment
************************************************************************************




* moment evaluator program for IPW - link the two equations
program gmm_ipw
	version 13.1
	syntax varlist if [fw aw iw pw], at(name) depvar(varname numeric) treat(varname numeric)
	gettoken resid1 resid2: varlist
	tempvar xb1 pr wgt xb2
	quietly {
		// Logit model  (propensity score model)
		matrix score double `xb1' = `at' `if', eq(logit)
		gen double `pr' = max(exp(`xb') / (1 + exp(`xb')), c(epsdouble)) `if'  // accounts for propensity scores near zero or one, but otherwise identical to the following: .gen `pr' = invlogit(`xb1') `if'
		replace `resid1' = `treat' - `pr' `if'

		// Construct weights
		gen `wgt' = cond(`treat', 1, `pr'/(1-`pr')) `if'
		summ `wgt' `if' [`weight'`exp'], meanonly
		replace `wgt' = `wgt' / r(mean)

		// Weighted OLS model
		matrix score `xb2' = `at' `if', eq(ols)
		replace `resid2' = `wgt' * (`depvar' - `xb2') `if'
	}
end
gmm gmm_ipw, equations(logit ols)  ///
	 instruments(logit: x1 x2 x3 x4 x5 x6 x7) instruments(ols: treat) ///
	 parameters(logit:x1 logit:x2 logit:x3 logit:x4 logit:x5 logit:x6 logit:x7 logit:_cons ols:treat ols:_cons) ///
	 depvar(y1) treat(treat) winitial(unadjusted, indep) onestep nolog

* teffects IPW (ATET)
teffects ipw (y1) (treat x1 x2 x3 x4 x5 x6 x7), aequations atet nolog

*/

local if if _n<=20
set seed 1
gen wgt = max(.1,rnormal(2,.4))
gen fwgt = round(rnormal(2,.4))
// forvalues i = 20/200 {
forvalues i = 20/25 {
gen x`i' = rnormal()
}
// expand 5e4 if touse
// expand 1e5 if touse

local depvar = "y1"
local treatvar = "treat"
local varlist = "x1 i.x2 i.x3 x4 x5 x6 x7 x2*"
// local varlist = "x*"
// local varlist = "x1 ib0.x2"
//local wgtvar = "wgt"
//local wgtvar = "fwgt"
local tousevar = "touse"
local estimate = "atet"

// some automatic parsing based on options above
if "`wgtvar'"!="" local wgtexp "[iw=`wgtvar']"
mark    `tousevar' `if' `in' `wgtexp'
markout `tousevar' `depvar' `treatvar' `varlist'
_rmdcoll `treatvar' `varlist' if `tousevar' `wgtexp', noconstant expand
// fvexpand `varlist' if `tousevar'
local varlist `r(varlist)'
forvalues j=1/`: list sizeof varlist' {
  local v : word `j' of `varlist'
  _ms_parse_parts `v'
  if !r(omit) local varlist1 `"`varlist1' `v'"'
}
local varlist : copy local varlist1


mata:

depvar   = st_local("depvar"  )
treatvar = st_local("treatvar")
wgtvar   = st_local("wgtvar"  )
varlist  = tokens(st_local("varlist"))
tousevar = st_local("tousevar")
estimate = st_local("estimate")

covars   = 0 // set =1 to calculate covariancess

X=T=Y=X0=X1=Y0=.
st_view(X, ., varlist, tousevar)
st_view(T, ., treatvar, tousevar)
st_select(X0, X, !T)
st_select(X1, X,  T)
if (depvar!="") {
  st_view(Y, ., depvar, tousevar)
  st_select(Y0, Y, !T)
}

if (wgtvar=="") W=W0=W1=1
else {
  W=W0=W1=.
  st_view(W, ., wgtvar, tousevar)
  st_select(W0, W, !T)
  st_select(W1, W,  T)
}
// to normalize weights
// w_norm = w :/ (rows(w) / quadcolsum(w))

// sample sizes
if (W==1) {
  N0 = N0_raw = rows(X0)
  N1 = N1_raw = rows(X1)
}
else {
  N0 = quadcolsum(W0)
  N1 = quadcolsum(W1)
  N0_raw = rows(X0)
  N1_raw = rows(X1)
}

// means
means0 = mean(X0,W0)
means1 = mean(X1,W1)

// Define function colvariance(x,w) == diagonal(quadvariance(x,w))'
// This function can be a lot faster than quadvariance, especially when you have lots of columns.
// Optionally, you can provide weights and/or provide a rowvector with the column means.
// For testing, mreldif(colvariance(X, w), diagonal(quadvariance(X, w))') should be small
real rowvector colvariance(real matrix X,| real colvector w, real rowvector Xbar)
{
  real rowvector v
  if (args()==1) w = 1
  if (args()<3) Xbar = mean(X,w)

  if (w==1) v = quadcolsum( (X:-Xbar):^2)     / (rows(X)-1)
  else      v = quadcolsum(((X:-Xbar):^2):*w) / (quadcolsum(w)-1)
  return(v)
}

if (covars) {
  covariance0 = quadvariance(X0,W0)
  covariance1 = quadvariance(X1,W1)
  variance0   = diagonal(covariance0)'
  variance1   = diagonal(covariance1)'
}
else {
  variance0 = colvariance(X0,W0,means0)
  variance1 = colvariance(X1,W1,means1)
}

//  For testing, this should be small:
mreldif(colvariance(X0, W0), diagonal(quadvariance(X0, W0))')

// variable-by-variable measures of imbalance
diff     = means1 :- means0
std_diff = diff :/ sqrt((variance1:+variance0)/2)
ratio    = variance1 :/ variance0
if (covars) covratio = covariance1 :/ covariance0

// scalar summaries
mean_asd = mean(abs(std_diff'))
min_asd  =  min(abs(std_diff ))
max_asd  =  max(abs(std_diff ))

// quick display
( N1_raw+N0_raw, N1+N0 \ N1_raw, N1 \ N0_raw, N0)
(varlist', strofreal(round((means0 \ means1 \ diff \ std_diff \ variance0 \ variance1 \ ratio)' , .0001)))
(mean_asd , min_asd , max_asd)

// if (covars) {
//   covariance0
//   covariance1
//   covarratio
// }


end

exit
qui teffects ipw (`depvar') (`treatvar' `varlist') if `tousevar' `wgtexp', `estimate' aequations
set tracedepth 3
//set trace on
tebalance summarize, baseline
mat list r(table)
tebalance summarize
mat list r(table)
mat list r(size)
