//****************************************************************************/
*! $Id$
*! Generalization of IPW and CBPS estimators
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
          [ depvars(varlist numeric)        /// outcome variables (if any)
            ate atet ateu                   /// to fill in est
            ipw cbps mean_sd sd mean_sd_sq sd_sq            /// to fill in fctn and oid
            TREatvariance CONtrolvariance POOledvariance Averagevariance /// to fill in denominator
            cvopt(numlist min=3 max=3)      ///
            /// vce(passthru) optimize?
            * ] //  display options are allowed

  marksample tousevar
  _get_diopts diopts /*options*/, `options'
  if ("`weight'"!="") {
    tempvar wgtvar
    gen double `wgtvar'=`exp'
/* */ //    local wgtexp [`weight'=`exp']
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

  // check type of dependent variables (if any)
  foreach v of local depvars {
    markout `tousevar' `v'
    _fv_check_depvar `v'
  }

  // parse the "est" options
  local est "`ate'`atet'`ateu'"
  if ("`est'"=="") local est ate
  else if (!inlist("`est'", "ate", "atet", "ateu") {
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
  if ("`cvopt'"=="") local cvopt "0 0 0"

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


  ereturn clear
  return  clear
  mata: Estimate()

/*
  // add to ereturn
  local retlist `retlist' N_g g_min g_max g_avg Tbar0 Tbar1 Tcon med_n depvar treat mvarlist panelvar ivar predict cmd2 properties marginsnotok wgtexp gmatch_eq_n gmatch_eq_names sigmas
  foreach ret of local retlist {
    if inlist("`ret'","predict","properties","marginsnotok") continue // I know I'm overwriting these locals
    if (`"`e(`ret')'"'=="." | `"`e(`ret')'"'=="")            continue // otherwise, return an error if the local/scalar is already defined
    di as error `"Programming error: e(`ret') is already defined:  e(`ret')=`e(`ret')'"'
    error 152
  }
  ereturn scalar N_g             = `N_g'
  ereturn scalar g_min           = `g_min'
  ereturn scalar g_max           = `g_max'
  ereturn scalar g_avg           = `g_avg'
  ereturn scalar Tbar0           = `tbar0'
  ereturn scalar Tbar1           = `tbar1'
  ereturn scalar Tcon            = (`g_min'==`g_max')
  ereturn scalar med_n           = `med_n'
  ereturn local  depvar          "`depvar'"
  ereturn local  treat           "`treat'"
  ereturn local  mvarlist        "`mvarlist'"
  ereturn local  panelvar        "`panelvar'"
  ereturn local  ivar            "`panelvar'"
  ereturn local  predict         "gmatch_predict"
  ereturn local  cmd2            "gmatch"
  ereturn local  properties      "`e(properties)' mi"
  ereturn local  marginsnotok    "E U UE SCore STDP XBU FE XBFE Ratio"
  ereturn local  wgtexp          "`wgtexp'"
  ereturn local  gmatch_eq_n     = `gmatch_eq_n'
  ereturn local  gmatch_eq_names = "`gmatch_eq_names'"


  // display sigma_u and sigma_e and rho
  if ("`sigmas'"=="sigmas") {
    local p=0
    forvalues i=1/`med_n' {
      forvalues t=0/1 {
        ereturn scalar `m`i''_`t'_sigma_u  = exp([`m`i''_`t'_lnsig_u]_b[_cons])
        ereturn scalar `m`i''_`t'_sigma_e  = exp([`m`i''_`t'_lnsig_e]_b[_cons])
        ereturn scalar `m`i''_`t'_rho      = exp([`m`i''_`t'_lnsig_u]_b[_cons])^2 / (exp([`m`i''_`t'_lnsig_u]_b[_cons])^2 + exp([`m`i''_`t'_lnsig_e]_b[_cons])^2)
        ereturn hidden local diparm`++p' `m`i''_`t'_lnsig_u, exp label("`m`i''_`t'_sigma_u") ci(log) noprob
        ereturn hidden local diparm`++p' `m`i''_`t'_lnsig_e, exp label("`m`i''_`t'_sigma_e") ci(log) noprob
        ereturn hidden local diparm`++p' `m`i''_`t'_lnsig_u `m`i''_`t'_lnsig_e, label(`m`i''_`t'_rho) func(exp(@1)^2/(exp(@1)^2+exp(@2)^2)) der( 2*exp(@1)*(exp(@2)/(exp(@1)^2+exp(@2)^2))^2 -2*exp(@2)*(exp(@1)/(exp(@1)^2+exp(@2)^2))^2 ) ci(probit) noprob
      }
    }
    di "Ancillary parameters:"
    _coef_table, neq(0) notest `diopts'
    ereturn scalar sigmas = 1
  }
  else {
    ereturn scalar sigmas = 0
  }

  // post-estimation commands
  if (`gmatch_eq_n'>0) di as txt _n "Post-estimation commands. Results will be stored in e(EQ_estimate), e(EQ_se), and e(EQ_cross_var), where EQ is the name of the equation.)"
  foreach eq of local gmatch_eq_names {

    // pull out matrix for equation `eq'
    matrix `submat_b' = e(b)
    matrix `submat_b' = `submat_b'[1,"`eq':"]
    matrix `submat_V' = e(V)
    matrix `submat_V' = `submat_V'["`eq':","`eq':"]

    // for gamma_D, we only want the site#treat interaction terms -- the second half of the matrix
    if  (`med_n'==2) & "`eq'"=="gamma_D" {
      local C = colsof(`submat_b')/2+1
      matrix `submat_b' = `submat_b'[1     ,`C'...]
      matrix `submat_V' = `submat_V'[`C'...,`C'...]
    }

    mata: st_local("this_lincom_exp", gmatch_build_lincom_exp("`submat_b'"))

    di _n as txt `"Estimator for the ``eq'_title' (from equation `eq')"' _n ///
          as inp `">> . lincom `this_lincom_exp'"'
    lincom `this_lincom_exp'

    ereturn scalar `eq'_estimate          = r(estimate)
    ereturn scalar `eq'_se                = r(se)
    if !mi(r(df)) ereturn scalar `eq'_df  = r(df)
    ereturn local `eq'_title              `"``eq'_title'"'
    ereturn local `eq'_cmd                `"lincom `this_lincom_exp'"'

    di _n as txt `"Estimate of the cross-site variance for the ``eq'_title' (from equation `eq')"'
    mata: st_local("thisvariance", strofreal( gmatch_cross_site_variance("`submat_b'","`submat_V'")))
    di as res `thisvariance'

    ereturn scalar `eq'_cross_var = `thisvariance'
  }

// this code can eventually be deleted
  // check that I actually returned everthing I thought I wanted to
  foreach ret of local retlist {
    if inlist("`ret'","wgtexp","gmatch_eq_names") continue // these are optional
    if (`"`e(`ret')'"'!="." | `"`e(`ret')'"'!="") continue // otherwise, throw error if the local/scalar is not defined
    di as error `"Programming error: e(`ret') is not defined:  e(`ret')=`e(`ret')'"'
    error 152
  }

  // checks difference initial values and final estimates; throws error if there are are large differences
  if (`med_n'==2) di _n "Initial values for the gamma equations are zero. We expect to see several iterations before achieving convergence."
  else            di _n "Initial values should be very close to final; we expect the model to converge quickly."
  mata: gmatch_compare_initval("`from'", `inittolerance')
*/
end

// LOAD THE CLASS
mata: mata mlib index

// DEFINE MATA FUNCTIONS
version 15.1
mata:
mata set matastrict on
mata set matafavor speed


// helper function to move Stata locals into Mata and call the main function
void Estimate()
{
  string scalar    treatvar, varlist, tousevar, wgtvar, depvars
  string scalar    est, fctn
  real   scalar    denominator, oid
  real   rowvector cvopt
  
  treatvar    = st_local("treatvar")
  varlist     = st_local("varlist")
  tousevar    = st_local("tousevar")
  wgtvar      = st_local("wgtvar")
  depvars     = st_local("depvars")
  
  est         = st_local("est")
  fctn        = st_local("fctn")
  denominator = strtoreal(st_local("denominator"))
  oid         = strtoreal(st_local("oid"))
  cvopt       = strtoreal(tokens(st_local("cvopt")))

  "(Loading data into mata.)"
  class gmatch scalar D
  D = gmatch()
  if  (wgtvar!="") D.set(treatvar, varlist, tousevar, wgtvar)
  else             D.set(treatvar, varlist, tousevar)
  if (depvars!="") D.set_Y(depvars,tousevar)
  
  "(Estimating model.)"
  D.cbps(est, fctn, denominator, oid, cvopt)
}
end

