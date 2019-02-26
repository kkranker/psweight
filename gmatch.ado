//****************************************************************************/
*! $Id$
*! Generalization of IPW and CBPS estimators
*! Stata command to estimate models
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.
//*****************************************************************************/

program define gmatch, eclass byable(onecall)
  version 15.1
  if replay() {
    if ("`e(cmd2)'" != "gmatch") error 301
    if _by() error 190
    if `"`0'"'=="" local 0 ","
    ereturn display
    exit
  }
  if _by() {
    local BY `"by `_byvars'`_byrc0':"'
  }

  `BY' Estimate `0'
end

program Estimate, eclass sortpreserve byable(recall)
  version 15.1

  // standard syntax parsing
  local cmdline : copy local 0
  syntax varlist(min=2 numeric fv) [if] [in] [fw iw/], ///
          [ DEPvars(varlist numeric) /// outcome variables (if any)
            ate atet ateu /// to fill in est
            ipw cbps mean_sd sd mean_sd_sq sd_sq STDProgdiff /// to fill in fctn and oid
            TREatvariance CONtrolvariance POOledvariance Averagevariance /// to fill in denominator
            cvtarget(numlist min=3 max=3) skewtarget(numlist min=3 max=3) kurttarget(numlist min=3 max=3) maxtarget(numlist min=3 max=3) ///
            from(name) /// starting values for maximization
            BALanceonly MWeight(varname numeric) /// just checks balance (skips reweighting)
            * ] //  display and ml options are allowed

  marksample tousevar
  _get_diopts diopts options, `options'
  get_matrix_table_options  , `options' `diopts'
  local matrix_table_options = s(opts)
  mlopts      mlopts        , `options'  // rest is not specified, so any other options will cause error
  if ("`weight'"!="") {
    tempvar wgtvar
    qui gen double `wgtvar'=`exp'
    local wgtexp [`weight'=`exp']
  }

  // check treatment variable
  gettoken treatvar varlist: varlist
  _fv_check_depvar `treatvar'
  cap assert inlist(`treatvar',0,1) if `tousevar'
  if _rc {
    di as err `"The treatment variable (`treatvar') must be a dummy variable."'
    error 125
  }
  sum `treatvar' if `tousevar', mean
  cap assert 0<r(mean) & r(mean) <1
  if _rc {
    di as err `"The treatment variable (`treatvar') must be a dummy variable with >1 treatment obs and >1 control obs."'
    error 125
  }
  if ("`balanceonly'"=="" & "`mweight'"!="") {
    di as err `"The mweight(`treatvar') option is only applicable with balanceonly."'
    error 198
  }
  else if ("`mweight'"!="") {
    markout `tousevar' `mweight'
  }

  // mark collinear variables
  if ("`balanceonly'"=="balanceonly") _rmcoll `treatvar' `varlist' if `tousevar' `wgtexp', expand
  else  _rmcoll `treatvar' `varlist' if `tousevar' `wgtexp', expand logit touse(`tousevar')
  local varlist `r(varlist)'
  gettoken trash varlist : varlist

  // check type of dependent variables (if any)
  foreach v of local depvars {
    markout `tousevar' `v'
    _fv_check_depvar `v'
  }

  // parse the "est" options
  local est "`ate'`atet'`ateu'"
  if ("`est'"=="") local est ate
  else if (!inlist("`est'", "ate", "atet", "ateu")) {
    di as err `"Specify one of the following: ate, atet, or ateu"'
    error 198
  }

  // parse the "fctn" and "oid" options
  if ("`mean_sd'"=="mean_sd") local mean_sd_sq mean_sd_sq
  if ("`sd'"=="sd")           local sd_sq sd_sq
  local fctn "`ipw'`cbps'`mean_sd_sq'`sd_sq'`stdprogdiff'"
  if ("`balanceonly'"=="balanceonly" & "`fctn'"!="") di "`fctn' ignored"
  else if ("`fctn'"=="") local fctn cbps
  else if (!inlist("`fctn'", "ipw", "cbps", "ipwcbps", "mean_sd_sq", "sd_sq","stdprogdiff")) {
    di as err `"Specify a valid combination of options: ipw, cbps, mean_sd, sd, or stdprogdiff."'
    error 198
  }

  // parse the "cvopt" option
  if (!mi("`maxtarget'")  & mi("`kurttarget'")) local kurttarget "0 0 2"
  if (!mi("`kurttarget'") & mi("`skewtarget'")) local skewtarget "0 0 2"
  if (!mi("`skewtarget'") & mi("`cvtarget'"))   local cvtarget   "0 0 2"
  local cvopt "`cvtarget' `skewtarget' `kurttarget' `maxtarget'"
  local cvopt : list clean cvopt
  if ("`balanceonly'"=="balanceonly" & "`cvopt'"!="") di "`cvopt' ignored"
  else if (!inlist(`: list sizeof cvopt',0,3,6,9,12)) {
    di as error `"cvopt() requires 3, 6, 9, 12 elements"'
    error 198
  }

  // parse the "denominator" options
  local denominator "`treatvariance'`controlvariance'`pooledvariance'`averagevariance'"
  if ("`denominator'"=="")                     local denominator = 1
  else if ("`denominator'"=="controlvariance") local denominator = 0
  else if ("`denominator'"=="treatvariance")   local denominator = 1
  else if ("`denominator'"=="pooledvariance")  local denominator = 2
  else if ("`denominator'"=="averagevariance") local denominator = 3
  else {
    di as err `"Specify one of the following: controlvariance, treatvariance, pooledvariance, or averagevariance"'
    error 198
  }

  // clear existing results  (varnames match those from psmatch2)
  foreach v in _weight _weight_mtch _pscore _treated {
    cap drop `v'
    if ("`balanceonly'"=="balanceonly") continue
    qui gen double `v' = .
    format %7.3g `v'
  }
  if ("`balanceonly'"!="balanceonly") ereturn clear
  return  clear

  // balanceonly option just prints balance and then end the program
  if ("`balanceonly'"=="balanceonly") {
    mata: BalanceOnly()
    exit
  }

  // switch over to Mata, helper function runs the main function
  mata: Estimate()

  // print results to screen
  di as txt _n "Propensity score model coefficients" _c
  di as txt _col(52) "Number of obs" _col(67) "=" _col(69) as res %10.0fc `gmatch_N_out'
  di as txt "Generalization of IPW/CPBS-type reweigting"
  if      ("`fctn'"=="ipw"        ) di as txt "Loss = IPW" _c
  else if ("`fctn'"=="cbps"       ) di as txt "Loss = CBPS" _c
  else if ("`fctn'"=="ipwcbps"    ) di as txt "Loss = CBPS + IPW (overidentified)" _c
  else if ("`fctn'"=="mean_sd_sq" ) di as txt "Loss = mean(stddiff())^2" _c
  else if ("`fctn'"=="sd_sq"      ) di as txt "Loss = sum(stddiff()^2)" _c
  else if ("`fctn'"=="stdprogdiff") di as txt "Loss = sum(stdprogdiff()^2)" _c
  tokenize `cvopt'
  if ("`1'"!="")  di as txt   " + `1'*abs(wgt_cv()-`2')^`3')" _c
  if ("`4'"!="")  di as txt   " + `4'*abs(wgt_skewness()-`5')^`6')" _c
  if ("`7'"!="")  di as txt   " + `7'*abs(wgt_kurtosis()-`8')^`9')" _c
  if ("`10'"!="") di as txt   " + `10'*abs(wgt_max()-`11')^`12')" _c
  di ""
  ereturn post `gmatch_beta_out' `wgtexp', obs(`gmatch_N_out') buildfvinfo esample(`tousevar')
  ereturn local est                       = "`est'"
  ereturn local fctn                      = "`fctn'"
  ereturn local depvar                    = "`treatvar'"
  ereturn local varlist                   = "`varlist'"
  ereturn local cmd                       = "gmatch"
  ereturn local cmdline                   = "gmatch `cmdline'"
  if ("`weight'"!="") ereturn local wtype = "`weight'"
  if ("`wexp'"!="")   ereturn local wexp  = "`wexp'"
  ereturn scalar denominator              = `denominator'
  if ("`cvopt'"!="")  ereturn local cvopt = "`cvopt'"
  _coef_table, `diopts'

  // print distribution of weights to screen
  di as txt _n "New variables (unweighted summary statistics)"
  tabstat _weight _weight_mtch _pscore if e(sample), by(_treated) c(s) s(N mean sd min p1 p10 p25 p50 p75 p90 p99 max) format
  return clear

end

program define get_matrix_table_options, sclass
  syntax [, formats(passthru) NOOMITted vsquish NOEMPTYcells BASElevels ALLBASElevels NOFVLABel fvwrap(passthru) fvwrapon(passthru) nolstretch *]
  sreturn local opts = strrtrim(stritrim(`"`formats' `noomitted' `vsquish' `noemptycells' `baselevels' `passthru' `allbaselevels' `nofvlabel' `fvwrap' `fvwrapon' `nolstretch'"'))
end

// DEFINE MATA FUNCTIONS
version 15.1
mata:
mata set matastrict on
mata set matafavor speed

// helper function to move Stata locals into Mata and call the main function
void Estimate()
{
  external class gmatch scalar gmatch_ado_most_recent
  string scalar    treatvar, varlist, tousevar, wgtvar, depvars
  string scalar    est, fctn
  real   scalar    denominator, oid
  real   rowvector cvopt
  transmorphic temp

  treatvar    = st_local("treatvar")
  varlist     = st_local("varlist")
  tousevar    = st_local("tousevar")
  wgtvar      = st_local("wgtvar")
  depvars     = st_local("depvars")

  est         = st_local("est")
  fctn        = st_local("fctn")
  denominator = strtoreal(st_local("denominator"))
  if  (st_local("cvopt")!="") {
    cvopt       = strtoreal(tokens(st_local("cvopt")))
  }
  else cvopt = J(1,0,.)

  gmatch_ado_most_recent = gmatch()
  if  (wgtvar!="") gmatch_ado_most_recent.set(treatvar, varlist, tousevar, wgtvar)
  else             gmatch_ado_most_recent.set(treatvar, varlist, tousevar)
  if (depvars!="") gmatch_ado_most_recent.set_Y(depvars,tousevar)

  temp = gmatch_ado_most_recent.gmatch(est, fctn, denominator, cvopt)
  gmatch_ado_most_recent.get_scores("_weight _weight_mtch _pscore _treated", tousevar)

}

// helper function to move Stata locals into Mata and call the main function
void BalanceOnly()
{
  external class gmatch scalar gmatch_ado_most_recent
  string scalar    treatvar, varlist, tousevar, wgtvar, depvars
  real   scalar    denominator
  transmorphic temp

  treatvar    = st_local("treatvar")
  varlist     = st_local("varlist")
  tousevar    = st_local("tousevar")
  wgtvar      = st_local("wgtvar")
  depvars     = st_local("depvars")
  mweightvar  = st_local("mweight")
  est         = st_local("est")
  denominator = strtoreal(st_local("denominator"))

  gmatch_ado_most_recent = gmatch()
  if  (wgtvar!="") gmatch_ado_most_recent.set(treatvar, varlist, tousevar, wgtvar)
  else             gmatch_ado_most_recent.set(treatvar, varlist, tousevar)
  if (depvars!="") gmatch_ado_most_recent.set_Y(depvars,tousevar)
  if (mweightvar!="") {
    gmatch_ado_most_recent.userweight(mweightvar, tousevar)
  }
  else if (wgtvar!="") {
    gmatch_ado_most_recent.userweight()
  }
  temp = gmatch_ado_most_recent.balanceresults(est, denominator)
}

end
