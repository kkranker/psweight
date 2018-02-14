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

program Estimate, eclass sortpreserve
  version 15.1

  // standard syntax parsing
  syntax varlist(min=2 numeric fv ts) [if] [in] [fw iw/], ///
          [ depvars(varlist numeric) /// outcome variables (if any)
            ate atet ateu /// to fill in est
            ipw cbps mean_sd sd mean_sd_sq sd_sq /// to fill in fctn and oid
            TREatvariance CONtrolvariance POOledvariance Averagevariance /// to fill in denominator
            cvtarget(numlist min=3 max=3) skewtarget(numlist min=3 max=3) kurttarget(numlist min=3 max=3) maxtarget(numlist min=3 max=3) ///
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

  // remove collinear variables   
  _rmcoll `treatvar' `varlist' if `tousevar' `wgtexp', expand logit touse(`tousevar')
  local varlist `r(varlist)'
  
  // Just in case there's some oddball situation where Mata would read in the blank columns for omitted variables
  fvexpand `varlist' if `tousevar'
  local varlist `r(varlist)'
  gettoken trash varlist : varlist
  forvalues j=1/`: list sizeof varlist' {
    local v : word `j' of `varlist'
    _ms_parse_parts `v'
    if !r(omit) local varlist1 `"`varlist1' `v'"'
  }
  local varlist : copy local varlist1

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
  if ("`fctn'"=="mean_sd") local mean_sd_sq mean_sd_sq
  else if ("`fctn'"=="sd") local sd_sq sd_sq
  local fctn "`ipw'`cbps'`mean_sd_sq'`sd_sq'"
  if ("`fctn'"=="") local fctn sd_sq
  else if (!inlist("`fctn'", "ipw", "cbps", "ipwcbps", "mean_sd_sq", "sd_sq")) {
    di as err `"Specify a valid combination of options: ipw, cbps, mean_sd, sd."'
    error 198
  }
  if ("`fctn'"=="ipwcbps") {
    local fctn cbps
    local oid = 1
  }
  else {
    local oid = 0
  }

  // parse the "cvopt" option
  if (!mi("`maxtarget'")  & mi("`kurttarget'")) local kurttarget "0 0 2"
  if (!mi("`kurttarget'") & mi("`skewtarget'")) local skewtarget "0 0 2"
  if (!mi("`skewtarget'") & mi("`cvtarget'"))   local cvtarget   "0 0 2"
  local cvopt "`cvtarget' `skewtarget' `kurttarget' `maxtarget'"
  local cvopt : list clean cvopt
  if (!inlist(`: list sizeof cvopt',0,3,6,9,12)) {
    di as error `"cvopt() requires 3, 6, 9, 12 elements"'
    error 198    
  }
  if (trim("`cvopt'")!="") mac list _cvopt

  // parse the "denominator" options
  local denominator "`treatvariance'`controlvariance'`pooledvariance'`averagevariance'"
  if ("`denominator'"=="")                local denominator = 1
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

  // print results to screen
  di as txt _n "Propensity score model coefficients" _c
  di as txt _col(52) "Number of obs" _col(67) "=" _col(69) as res %10.0fc `gmatch_N_out'
  di as txt "Generalization of IPW/CPBS-type reweigting"
  if      ("`fctn'"=="ipw"       ) di as txt "Loss = IPW" _c
  else if ("`fctn'"=="cbps"      ) di as txt "Loss = CBPS" _c
  else if ("`fctn'"=="ipwcbps"   ) di as txt "Loss = CBPS + IPW (overidentified)" _c
  else if ("`fctn'"=="mean_sd_sq") di as txt "Loss = mean(stddiff())^2" _c
  else if ("`fctn'"=="sd_sq"     ) di as txt "Loss = stddiff():^2" _c
  if ("`cvopt'"!="") {
     gettoken a b : cvopt
     gettoken b c : b
     gettoken c d : c
     di as txt   " + `a'*abs(CV-`b')^`c')"
  }
  else di ""
  ereturn post `gmatch_beta_out' `wgtexp', obs(`gmatch_N_out') buildfvinfo esample(`tousevar')
  ereturn local est          = "`est'"
  ereturn local fctn         = "`fctn'"
  ereturn scalar denominator = `denominator'
  ereturn scalar oid         = `oid'
  ereturn local cvopt        = "`cvopt'"
  _coef_table, `diopts'

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
  oid         = strtoreal(st_local("oid"))
  if  (st_local("cvopt")!="") {  
    cvopt       = strtoreal(tokens(st_local("cvopt")))
  }
  else cvopt = J(1,0,.)

  gmatch_ado_most_recent = gmatch()
  if  (wgtvar!="") gmatch_ado_most_recent.set(treatvar, varlist, tousevar, wgtvar)
  else             gmatch_ado_most_recent.set(treatvar, varlist, tousevar)
  if (depvars!="") gmatch_ado_most_recent.set_Y(depvars,tousevar)

  temp = gmatch_ado_most_recent.cbps(est, fctn, denominator, oid, cvopt)
  gmatch_ado_most_recent.get_scores("_weight _weight_mtch _pscore _treated", tousevar)

}
end
