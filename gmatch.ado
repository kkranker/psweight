mac drop _all
// clear all
// cls
set varabbrev off
set scheme mpr_blue
set linesize 160
set maxiter 100
cap log close gmatch_example
local makegraphs = 01
cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\"

log using gmatch.log, name(gmatch_example) replace

// Multiple-equation models: An introduction and potential applications to our work at Mathematica
// Design and methods “brown bag” workshop
// April 8, 2015
// Keith Kranker

// Stata_code_2_IPW.do This program includes all the examples in the powerpoint slides, plus more.
// This program includes examples of how to code inverse propensity weighting (IPW) estimators using Stata's GMM command

// include C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\gmatchclass.mata

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

/*

************************************************************************************
* 5. IPW w/ binary treatment
************************************************************************************




* moment evaluator program for IPW - link the two equations
program gmm_ipw
	version 13.1
	syntax varlist if [fw aw iw pw], at(name) depvars(varname numeric) treat(varname numeric)
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
		replace `resid2' = `wgt' * (`depvars' - `xb2') `if'
	}
end
gmm gmm_ipw, equations(logit ols)  ///
	 instruments(logit: x1 x2 x3 x4 x5 x6 x7) instruments(ols: treat) ///
	 parameters(logit:x1 logit:x2 logit:x3 logit:x4 logit:x5 logit:x6 logit:x7 logit:_cons ols:treat ols:_cons) ///
	 depvars(y1) treat(treat) winitial(unadjusted, indep) onestep nolog

* teffects IPW (ATET)
teffects ipw (y1) (treat x1 x2 x3 x4 x5 x6 x7), aequations atet nolog

*/

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


// give the sample poor overlap
tab2    treat x1
replace x1 = 1 if  treat & runiform()<.85
replace x1 = 0 if !treat & runiform()<.85
tab2    treat x1 


local depvars = "y1 y1_binary"
local treatvar = "treat"
local varlist = "x1 i.x2 i.x3 x4 x5 x6 x7 x9*"
// local varlist = "x*"
// local varlist = "x1 ib0.x2"
//local wgtvar = "wgt"
local wgtvar = "fwgt"
local tousevar = "touse"
local estimate = "atet"

// some automatic parsing based on options above
if "`wgtvar'"!="" local wgtexp "[iw=`wgtvar']"
mark    `tousevar' `if' `in' `wgtexp'
markout `tousevar' `depvars' `treatvar' `varlist'
_rmdcoll `treatvar' `varlist' if `tousevar' `wgtexp', noconstant expand
// _rmcoll `treatvar' `varlist' if `tousevar' `wgtexp', noconstant expand logit touse(`tousevar')
// fvexpand `varlist' if `tousevar'
local varlist `r(varlist)'
forvalues j=1/`: list sizeof varlist' {
  local v : word `j' of `varlist'
  _ms_parse_parts `v'
  if !r(omit) local varlist1 `"`varlist1' `v'"'
}
local varlist : copy local varlist1


mata:

depvars  = st_local("depvars"  )
treatvar = st_local("treatvar")
wgtvar   = st_local("wgtvar"  )
varlist  = st_local("varlist" )
tousevar = st_local("tousevar")
estimate = st_local("estimate")

// ****************************
// * UNWEIGHTED DATA EXAMPLES *
// ****************************

D = gmatch()
D.set(st_local("treatvar"),st_local("varlist") ,st_local("tousevar"))
if (depvars!="") D.set_Y(st_local("depvars"),st_local("tousevar"))

// Misc balance measures
D.diff()
D.stddiff()
D.mean_asd()
D.stddiff(1)
D.stddiff(0)
D.varratio()
D.prognosticdiff()
table = D.balancetable(1)


// CBPS
stata(`"cbps `treatvar' `varlist' if `tousevar' , att logit optimization_technique("nr") evaluator_type("gf1")"')
stata(`"cbps_imbalance"')

M = gmatch()
M.clone(D)
// cbpsweight = M.cbps("atet","mean_sd_sq",1)
// cbpsweight = D.cbps("atet","mean_asd")
// cbpsweight = D.cbps("atet","max_asd")
// cbpsweight = D.cbps("atet","cbpslossST")
// cbpsweight = D.cbps("atet","cbpsloss")
cbpsweight = D.cbps("atet","cbpseval_port",3)
cbpsweight = D.cbps("atet","sd_sq",1)
cbpsweight = D.cbps("atet","sd_sq_ent",1)
cbpsweight = D.cbps("atet","sd_sq_cv",1)

// cbpsweight = D.cbps("atet","cbpslossOID")
// cbpsweight = D.cbps("atet","cbpsloss")
// cbpsweight = D.cbps("atet","cbpsloss2")


_error("Stop here")

// IPW
M = gmatch()
M.clone(D)
stata("qui teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations")
iwpweight = D.ipw("atet")
M.multweight(iwpweight)

stata("di _b[POmean:0.treat]")
M.pomean()

stata("tebalance summarize")
table = D.balancetable(3)
table = M.balancetable(3)

_error("stop")

stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations")
stata("predict pscore1, tlevel(1) ")
stata("list `treatvar' pscore1 in 1/20, nolab ")
ipw = D.ipw("atet")

stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , ate aequations")
ipw = D.ipw("ate")

stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' , atet aequations tlevel(0) control(1)")
ipw = D.ipw("ateu")


// **************************
// * WEIGHTED DATA EXAMPLES *
// **************************

D = gmatch()
D.set(st_local("treatvar"),st_local("varlist") ,st_local("tousevar"),st_local("wgtvar"))
if (depvars!="") D.set_Y(st_local("depvars"),st_local("tousevar"))

D.stddiff(1)
D.mean_asd
D.min_asd
D.max_asd
D.prognosticdiff()

stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations")
stata("tebalance summarize, baseline")
table = D.balancetable()
_error("weights from tebalance are weird")

stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], atet aequations")
ipw = D.ipw("atet")

stata("teffects ipw (`:word 1 of `depvars'') (`treatvar' `varlist') if `tousevar' [iw=`wgtvar'], ate aequations")
ipw = D.ipw("ate")


/*
// this is wrong
diff_invD     = diff * invsym(quadvariance(X,W))
diff_invDdiag = diff * invsym(diag(diagvariance(X,W)))

*/

end

qui teffects ipw (`depvars') (`treatvar' `varlist') if `tousevar' `wgtexp', `estimate' aequations
tebalance summarize, baseline
mat list r(table)
tebalance summarize
mat list r(table)
mat list r(size)


log close gmatch_example

