//****************************************************************************/
*! $Id$
*! IPW- and CBPS-type propensity score reweighting, with extentions
*! Stata command to estimate models
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.
//*****************************************************************************/

program define psweight, eclass byable(onecall)
  version 15.1
  if _by() {
    local BY `"by `_byvars'`_byrc0':"'
  }
  local cmdline `"psweight `0'"'
	local cmdline : list retokenize cmdline
  if replay() {
    if ("`e(cmd)'" != "psweight") error 301
    if _by() error 190
    if `"`0'"'=="" local 0 ","
    ereturn display
    exit
  }
  gettoken fctn 0 : 0, parse(" ,")
  if ("`fctn'"=="call") {
    if _by() error 190
    psweightcall `0'
  }
  else {
    `BY' Estimate `fctn' `0'
  }
  ereturn local cmd "psweight"
  ereturn local cmdline `"`cmdline'"'
end

program Estimate, eclass sortpreserve byable(recall)
  version 15.1

  // first word is the name of the "subroutine"
  gettoken fctn 0 : 0, parse(" ,")
  if      ("`fctn'"=="mean_sd") local fctn mean_sd_sq
  else if ("`fctn'"=="sd")      local fctn sd_sq
  else if !inlist(`"`fctn'"', "balanceonly", "pcbps", "ipw", "cbps", "cbpsoid", "mean_sd_sq", "sd_sq", "stdprogdiff") {
    di as error `""`fctn'" subcommand invalid"'
    error 198
  }

  // standard syntax parsing
  syntax varlist(min=2 numeric fv) [if] [in] [fw iw/], ///
          [ DEPvars(varlist numeric) /// outcome variables (if any)
            ate atet ateu /// to fill in est
            TREatvariance CONtrolvariance POOledvariance Averagevariance /// to fill in denominator
            cvtarget(numlist min=3 max=3) skewtarget(numlist min=3 max=3) kurttarget(numlist min=3 max=3) maxtarget(numlist min=3 max=3) ///
            from(name) /// starting values for maximization
            MWeight(varname numeric) /// 'matching weights' for balanceonly option
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
  if ("`fctn'"!="balanceonly" & "`mweight'"!="") {
    di as err `"The mweight(`treatvar') option is only applicable with balanceonly."'
    error 198
  }
  else if ("`mweight'"!="") {
    markout `tousevar' `mweight'
  }

  // mark collinear variables
  if ("`fctn'"=="balanceonly") _rmcoll `treatvar' `varlist' if `tousevar' `wgtexp', expand
  else                         _rmcoll `treatvar' `varlist' if `tousevar' `wgtexp', expand logit touse(`tousevar')
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

  // parse the "cvopt" option
  if (!mi("`maxtarget'")  & mi("`kurttarget'")) local kurttarget "0 0 2"
  if (!mi("`kurttarget'") & mi("`skewtarget'")) local skewtarget "0 0 2"
  if (!mi("`skewtarget'") & mi("`cvtarget'"))   local cvtarget   "0 0 2"
  local cvopt "`cvtarget' `skewtarget' `kurttarget' `maxtarget'"
  local cvopt : list clean cvopt
  if ("`fctn'"=="balanceonly" & "`cvopt'"!="") {
    di as error "`cvopt' not allowed with `fctn' subcommand"
    error 198
  }
  else if ("`fctn'"=="pcbps" & "`cvopt'"=="") {
    di as error `"cvtarget(), skewtarget(), or kurttarget() required with pcbps subcommand"'
    error 198
  }
  else if (!inlist(`: list sizeof cvopt',0,3,6,9,12)) {
    di as error `"cvopt() requires 3, 6, 9, 12 elements"'
    error 198
  }
  if ("`fctn'"=="pcbps") local fctn cbps // pcbps is a synonym of cbps with cvopt()

  // parse the "denominator" options
  local denominator "`treatvariance'`controlvariance'`pooledvariance'`averagevariance'"
  if ("`denominator'"=="")                     local denominator = 2
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
  if ("`fctn'"!="balanceonly") {
    ereturn clear
    cap mata: mata drop psweight_ado_most_recent
  }
  return clear

  // balanceonly option just prints balance and then end the program
  if ("`fctn'"=="balanceonly") {
    mata: Estimate(0)
    ereturn local est                       = "`est'"
    ereturn local depvar                    = "`treatvar'"
    ereturn local varlist                   = "`varlist'"
    ereturn scalar balanceonly              = 1
    if ("`weight'"!="") ereturn local wtype = "`weight'"
    if ("`wexp'"!="")   ereturn local wexp  = "`wexp'"
    ereturn local mataobj                   = "psweight_ado_most_recent"
    exit
  }

  // switch over to Mata, helper function runs the main function
  mata: Estimate(1)

  // print results to screen
  di as txt _n "Propensity score model coefficients" _c
  di as txt _col(52) "Number of obs" _col(67) "=" _col(69) as res %10.0fc `psweight_N_out'
  di as txt "Propensity score reweigting"
  if      ("`fctn'"=="ipw"        ) di as txt "Loss = IPW" _c
  else if ("`fctn'"=="cbps"       ) di as txt "Loss = CBPS (just identified)" _c
  else if ("`fctn'"=="cbpsoid"    ) di as txt "Loss = CBPS (over identified)" _c
  else if ("`fctn'"=="mean_sd_sq" ) di as txt "Loss = mean(stddiff())^2" _c
  else if ("`fctn'"=="sd_sq"      ) di as txt "Loss = sum(stddiff()^2)" _c
  else if ("`fctn'"=="stdprogdiff") di as txt "Loss = sum(stdprogdiff()^2)" _c
  tokenize `cvopt'
  if ("`1'"!="")  di as txt   " + `1'*abs(wgt_cv()-`2')^`3')" _c
  if ("`4'"!="")  di as txt   " + `4'*abs(wgt_skewness()-`5')^`6')" _c
  if ("`7'"!="")  di as txt   " + `7'*abs(wgt_kurtosis()-`8')^`9')" _c
  if ("`10'"!="") di as txt   " + `10'*abs(wgt_max()-`11')^`12')" _c
  di ""
  ereturn post `psweight_beta_out' `wgtexp', obs(`psweight_N_out') buildfvinfo esample(`tousevar')
  ereturn local est                       = "`est'"
  ereturn local fctn                      = "`fctn'"
  ereturn local depvar                    = "`treatvar'"
  ereturn local varlist                   = "`varlist'"
  ereturn scalar balanceonly              = 0
  if ("`weight'"!="") ereturn local wtype = "`weight'"
  if ("`wexp'"!="")   ereturn local wexp  = "`wexp'"
  ereturn local mataobj                   = "psweight_ado_most_recent"
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

program define psweightcall, rclass
  mata: `e(mataobj)'.`0'
  return add
end


// DEFINE MATA FUNCTIONS
version 15.1
mata:
mata set matastrict on
mata set matafavor speed

class psweightado extends psweight
{
  protected:
    string scalar    est
    string scalar    fctn
    real   scalar    denominator
    real   rowvector cvopt
  public:
    void   set_opts()
    void   userweight()
}

void psweightado::set_opts(string scalar    est_in,
                           string scalar    fctn_in,
                           real   scalar    denominator_in,
                           real   rowvector cvopt_in)
{
  this.est = est_in
  this.fctn = fctn_in
  this.denominator = denominator_in
  this.cvopt = cvopt_in
}

// sets this.W appropriately for the balanceonly option in the .ado file
void psweightado::userweight(| string scalar wgtvar, string scalar tousevar)
{
  if (args()==0 | wgtvar=="") this.reweight()
  else if (args()==2) {
    real colvector userweight
    userweight=.
    st_view(userweight, ., wgtvar, tousevar)
    if (length(userweight)!=length(this.T)) _error("Unexpected dimension for " + wgtvar)
    this.reweight(userweight)
  }
  else _error("userweight() requires 0 or 2 arguments")
}

// these functiosn are just wrappers
void           psweightado::balanceresults() return(this.super.balanceresults(this.est, this.denominator))
real rowvector psweightado::psweight()       return(this.super.psweight(this.est, this.fctn, this.denominator, this.cvopt))
real rowvector psweightado::ipw()            return(this.super.ipw(this.est))
real rowvector psweightado::cbps()           return(this.super.cbps(this.est, this.denominator))
real rowvector psweightado::cbpsoid()        return(this.super.cbpsoid(this.est, this.denominator))
real rowvector psweightado::stddiff()        return(this.super.stddiff(this.denominator))
real rowvector psweightado::varratio()       return(this.super.varratio(this.denominator))
real rowvector psweightado::progdiff()       return(this.super.progdiff(this.denominator))
real rowvector psweightado::stdprogdiff()    return(this.super.stdprogdiff(this.denominator))
real scalar    psweightado::mean_asd()       return(this.super.mean_asd(this.denominator))
real scalar    psweightado::max_asd()        return(this.super.max_asd(this.denominator))
real scalar    psweightado::wgt_cv()         return(this.super.wgt_cv(this.est))
real scalar    psweightado::wgt_sd()         return(this.super.wgt_sd(this.est))
real scalar    psweightado::wgt_skewness()   return(this.super.wgt_skewness(this.est))
real scalar    psweightado::wgt_kurtosis()   return(this.super.wgt_kurtosis(this.est))
real scalar    psweightado::wgt_max()        return(this.super.wgt_max(this.est))
real matrix    psweightado::balancetable()   return(this.super.balancetable(this.denominator))

// helper function to move Stata locals into Mata and call the main function
// reweight == 1: calcualte inverse propensity weights
//          == 0: just calcuate balance
void Estimate(real scalar reweight)
{
  external class   psweightado scalar psweight_ado_most_recent
  string scalar    treatvar, varlist, tousevar, wgtvar, depvars
  string scalar    est, fctn, mweightvar
  real   scalar    denominator
  real   rowvector cvopt
  transmorphic temp

  // access key parameters from Stata locals
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

  // initialize class and read in data and parameters
  psweight_ado_most_recent = psweightado()
  if  (wgtvar!="") psweight_ado_most_recent.set(treatvar, varlist, tousevar, wgtvar)
  else             psweight_ado_most_recent.set(treatvar, varlist, tousevar)
  if (depvars!="") psweight_ado_most_recent.set_Y(depvars,tousevar)
  psweight_ado_most_recent.set_opts(est, fctn, denominator, cvopt)

  // compute invere probabily weights into child class
  if (reweight) {
    temp = psweight_ado_most_recent.psweight(est, fctn, denominator, cvopt)
  }

  // just compute balance
  else {
    mweightvar  = st_local("mweight")
    psweight_ado_most_recent.userweight(mweightvar, tousevar)
    temp = psweight_ado_most_recent.balanceresults(est, denominator)
  }

  // stick obs-specific weigths and such into Stata vaiables
  psweight_ado_most_recent.get_scores("_weight _weight_mtch _pscore _treated", tousevar)
}

end
