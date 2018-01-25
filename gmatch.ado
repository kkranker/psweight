clear all
mac drop _all
cls
set varabbrev off
set scheme mpr_blue
set linesize 160
set maxiter 25
cap log close _all
local makegraphs = 01
cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\"


// Multiple-equation models: An introduction and potential applications to our work at Mathematica
// Design and methods “brown bag” workshop
// April 8, 2015
// Keith Kranker

// Stata_code_2_IPW.do This program includes all the examples in the powerpoint slides, plus more.
// This program includes examples of how to code inverse propensity weighting (IPW) estimators using Stata's GMM command

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

local depvar = "y1"
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
markout `tousevar' `depvar' `treatvar' `varlist'
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

include C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\gmatchclass.mata

mata:

depvar   = st_local("depvar"  )
treatvar = st_local("treatvar")
wgtvar   = st_local("wgtvar"  )
varlist  = st_local("varlist" )
tousevar = st_local("tousevar")
estimate = st_local("estimate")

covars   = 0 // set =1 to calculate covariancess


D = gmatch()
D.set_X(st_local("varlist") ,st_local("tousevar"))
D.set_T(st_local("treatvar"),st_local("tousevar"))

D.diff()
D.stddiff()
D.stddiff(1)
_error("stop")

if (depvar!="") D.set_Y(st_local("depvar"),st_local("tousevar"))
if (wgtvar!="") D.set_W(st_local("wgtvar"),st_local("tousevar"),0)
D.diff()
D.stddiff()


// to normalize weights


/* testing */ mreldif(diagvariance(X0, W0), diagonal(quadvariance(X0, W0))')

if (covars) {
  covariance0 = quadvariance(X0,W0)
  covariance1 = quadvariance(X1,W1)
  variance0   = diagonal(covariance0)'
  variance1   = diagonal(covariance1)'
}
else {
  variance0 = diagvariance(X0,W0,means0)
  variance1 = diagvariance(X1,W1,means1)
}


// variable-by-variable measures of imbalance
std_diff = diff :/ sqrt((variance1:+variance0)/2)
ratio    = variance1 :/ variance0
if (covars) covratio = covariance1 :/ covariance0

// scalar summaries
mean_asd = mean(abs(std_diff'))
min_asd  =  min(abs(std_diff ))
max_asd  =  max(abs(std_diff ))

// Define function to run OLS regression model (coefficients only)
// A contant term is included in the regression, but its coefficient
// is dropped from output (so number of coefficients matches number of columns)
real matrix olsbeta(real colvector y, real matrix X, | real colvector w)
{
  real colvector beta
  real scalar C
  real matrix XX, Xy
  if (args()<3) w=1
  C = cols(X)
  XX = quadcross(X, 1, w, X, 1)
  Xy = quadcross(X, 1, w, y, 0)
  beta  = invsym(XX,++C)*Xy
  return(beta[1..--C]')
}

// Define a function to get logit model coefficients
// Source: https://www.stata.com/statalist/archive/2010-10/msg01188.html
//    From   jpitblado@stata.com (Jeff Pitblado, StataCorp LP)
//    To   statalist@hsphsun2.harvard.edu
//    Subject   Re: st: pointing to a class member function (likelihood) with optimize() in mata
//    Date   Thu, 28 Oct 2010 10:19:51 -0500
class logit_model {
        real colvector  y
        real matrix     X
        void eval()
}
void logit_model::eval( real    scalar          todo,
                        real    rowvector       beta,
                        real    scalar          lnf,
                        real    rowvector       g,
                        real    matrix          H)
{
        real colvector  pm
        real colvector  xb
        real colvector  lj
        real colvector  dllj
        real colvector  d2llj

        pm      = 2*(this.y :!= 0) :- 1
        xb      = this.X*beta'

        lj      = invlogit(pm:*xb)
        if (any(lj :== 0)) {
                lnf = .
                return
        }
        lnf = quadcolsum(ln(lj))
        if (todo == 0) return

        dllj    = pm :* invlogit(-pm:*xb)
        if (missing(dllj)) {
                lnf = .
                return
        }
        g       = quadcross(dllj, X)
        if (todo == 1) return

        d2llj   = abs(dllj) :* lj
        if (missing(d2llj)) {
                lnf = .
                return
        }
        H       = - quadcross(X, d2llj, X)
}

void logit_eval(        real    scalar          todo,
                        real    rowvector       beta,
                        class   logit_model     M,
                        real    scalar          lnf,
                        real    rowvector       g,
                        real    matrix          H)
{
        M.eval(todo,beta,lnf,g,H)
}

M = logit_model()
M.y=T
M.X=(X,J(rows(X),1,1))

S = optimize_init()
optimize_init_evaluator(S, &logit_eval())
optimize_init_evaluatortype(S, "d2")
optimize_init_argument(S, 1, M)
optimize_init_params(S, J(1,cols(X)+1,0))
stata("logit `treatvar' `varlist' if `tousevar'")
optimize(S)


// Differences times the coefficients of OLS on Y using the treated observations
invsym(quadvariance(X,W))
diff_beta0    = diff :* olsbeta(Y0, X0, W0)

// stata("regress "+depvar+" "+invtokens(varlist)+" if 0=="+treatvar+" & 1=="+tousevar); olsbeta(Y0, X0, W0)
diff_alpha    = diff :* olsbeta(T, X, W)

// this is wrong
diff_invD     = diff * invsym(quadvariance(X,W))
diff_invDdiag = diff * invsym(diag(diagvariance(X,W)))



// quick display
( N1_raw+N0_raw, N1+N0 \ N1_raw, N1 \ N0_raw, N0)
(varlist', strofreal(round((means0 \ means1 \ diff \ std_diff \ diff_beta0 \ diff_alpha \ diff_invD \ diff_invDdiag \ variance0 \ variance1 \ ratio)' , .0001)))
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
