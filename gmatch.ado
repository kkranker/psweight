//****************************************************************************/
*! $Id: gmatch.ado,v 56130091156f 2018/03/02 17:22:24 kkranker $
*! Generalization of IPW and CBPS estimators
*! Stata command to estimate models
//
*! By Keith Kranker
// Last updated $Date: 2018/03/02 17:22:24 $
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

program Estimate, eclass sortpreserve
  version 15.1

  // standard syntax parsing
  local cmdline : copy local 0
  syntax varlist(min=2 numeric fv) [if] [in] [fw iw/], ///
          [ depvars(varlist numeric) /// outcome variables (if any)
            ate atet ateu /// to fill in est
            ipw cbps mean_sd sd mean_sd_sq sd_sq STDProgdiff /// to fill in fctn
            TREatvariance CONtrolvariance POOledvariance Averagevariance /// to fill in denominator
            cvtarget(numlist min=3 max=3) skewtarget(numlist min=3 max=3) kurttarget(numlist min=3 max=3) ///
            * ] //  display and ml options are allowed

  marksample tousevar
  _get_diopts diopts options, `options'
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

  // mark collinear variables
  _rmcoll `treatvar' `varlist' if `tousevar' `wgtexp', expand logit touse(`tousevar')
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
  if ("`fctn'"=="") local fctn cbps
  else if (!inlist("`fctn'", "ipw", "cbps", "ipwcbps", "mean_sd_sq", "sd_sq","stdprogdiff")) {
    di as err `"Specify a valid combination of options: ipw, cbps, mean_sd, sd, or stdprogdiff."'
    error 198
  }

  // parse the "cvopt" option
  // Instead of having lots of options, I just pass the mata functions one big vector of up to 9 elements.
  // It has 0 elements if none of these options are provided,
  //        3 elements if you're just using cvtarget(),
  //        6 elements if you're using skewtarget(), and
  //        9 elements if you're using kurttarget().
  //        If you are using skewtarget() without cvtarget(), then defaults keep cvtarget() from having any effect.
  //        If you are using kurttarget() without skewtarget() or cvtarget(), then defaults keep cvtarget() adn skewtarget() from having any effect.
  if (!mi("`kurttarget'") & mi("`skewtarget'")) local skewtarget "0 0 2"
  if (!mi("`skewtarget'") & mi("`cvtarget'"))   local cvtarget   "0 0 2"
  local cvopt "`cvtarget' `skewtarget' `kurttarget'"
  local cvopt : list clean cvopt
  if (!inlist(`: list sizeof cvopt',0,3,6,9)) {
    di as error `"cvopt() requires 3, 6, or 9 elements"'
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
    qui gen double `v' = .
    format %7.3g `v'
  }
  ereturn clear
  return  clear

  // switch over to Mata, run helper function with runs the main function
  mata: Estimate()
  local converged = r(converged)

  // print results to screen
  di as txt _n "Propensity score model coefficients" _c
  di as txt _col(52) "Number of obs" _col(67) "=" _col(69) as res %10.0fc `gmatch_N_out'
  di as txt "Generalization of IPW/CPBS-type reweighting"
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
  di ""
  ereturn post `gmatch_beta_out' `wgtexp', obs(`gmatch_N_out') buildfvinfo esample(`tousevar')
  _coef_table, `diopts'
  ereturn local est                       = "`est'"
  ereturn local fctn                      = "`: list clean fctn'"
  ereturn local depvar                    = "`: list clean treatvar'"
  ereturn local varlist                   = "`: list clean varlist'"
  ereturn local cmd                       = "gmatch"
  ereturn local cmdline                   = "gmatch `cmdline'"
  if ("`weight'"!="") ereturn local wtype = "`weight'"
  if ("`wexp'"!="")   ereturn local wexp  = "`wexp'"
  ereturn scalar denominator              = `denominator'
  ereturn local cvopt                     = "`cvopt'"
  ereturn scalar converged                = `converged'

  // print distribution of weights to screen
  di as txt _n "New variables (unweighted summary statistics)"
  tabstat _weight _weight_mtch _pscore if e(sample), by(_treated) c(s) s(N mean sd min p1 p10 p25 p50 p75 p90 p99 max) format
  return clear

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
  real   scalar    denominator
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

  temp = gmatch_ado_most_recent.gmatch(est, fctn, denominator, cvopt) // temp is not used; it just keeps the betas from being printed
  gmatch_ado_most_recent.get_scores("_weight _weight_mtch _pscore _treated", tousevar)
}
end
