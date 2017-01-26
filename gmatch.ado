clear all
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

gen touse = _n<=20
set seed 1
gen wgt = max(.1,rnormal(2,.4))
forvalues i = 20/200 {
gen x`i' = rnormal()
}
//expand 1e5 if touse
expand 5e4 if touse

mata:

depvar = "y1"
treatvar = "treat"
wgtvar = ""
wgtvar = "wgt"
// varlist = "x1 i.x2 i.x3 x4 x5 x6 x7 x20-x100"
varlist = "x*"
// varlist = "x1 i.x2"
tousevar = "touse"
estimate = "atet"

st_view(X=., ., varlist, tousevar)
st_view(T=., ., treatvar, tousevar)
st_view(Y=., ., depvar, tousevar)
st_select(X0=., X, !T)
st_select(X1=., X,  T)
st_select(Y0=., Y, !T)

if (wgtvar=="") W=W0=W1=1
else {
st_view(W=., ., wgtvar, tousevar)
st_select(W0=., W, !T)
st_select(W1=., W,  T)
}

// (T,Y,X)
// X1
// (Y0,X0)

X0bar = mean(X0,W0)
X1bar = mean(X1,W1)
// diff  = X0bar :- X1bar
// diff

// to normalize weights 
// w_norm = w :/ (rows(w) / quadcolsum(w))

// colvariance(x,w) equivalent to diagonal(quadvariance(x,w))'
// but it is a lot faster when there are lots of columns
real rowvector colvariance(real matrix X,| real colvector w, real rowvector means) {
  real rowvector v
  if (args()<2) w = 1
  if (args()<3) means = mean(X,w)
  
  // cand this be done with cross() or quadcross()? they are supposed to be quicker w/ views
  // the problem is that the () mean that I create this giant matrix before taking its sum
  if (w==1) v = quadcolsum( (X:-means):^2)     / (rows(X)-1)
  else      v = quadcolsum(((X:-means):^2):*w) / (quadcolsum(w)-1)
  return(v)
}


timer_on(1)
variances = colvariance(X,W) // , X0bar)
// colvariance(X0,W0) // , X0bar)
// colvariance(X1,W1) // , X1bar)
timer_off(1)

/*
timer_on(2)
diagonal(quadvariance(X, W))'
// diagonal(quadvariance(X0, W0))'
// diagonal(quadvariance(X1, W1))'

// : var   = meanvariance(X)
// : means = var[1,.]
// : var   = var[|2,1 \ .,.|]

timer_off(2)
*/
mreldif( variances, diagonal(quadvariance(X, W))' )

timer()

end

qui sum x1 if touse [iw=wgt]
di r(Var)
tabstat x* if touse [iw=wgt], s(mean var n) col(v)
