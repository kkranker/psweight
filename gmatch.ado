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

/*  QA COMMENTS:

From: Liz Potamites
Sent: Wednesday, March 28, 2018 1:01 PM
To: Keith Kranker <KKranker@mathematica-mpr.com>
Subject: gmatch comments

Hey Keith, I’m sorry I wasn’t able to do a more thorough job on this. When you do a file compare,
you’ll see that my comments are very superficial. Perhaps I can take another deep-dive with fresh
eyes next year and/or if/when you submit this to the common code repo (or elsewhere). I only looked
through the code and log files. I did not play around with testing the code myself.

Here’s a summary of what you’ll see when you file compare:

•	gmatch_one_time_setup.do – added comments with regard to the problem I had using this
•	gmatch.ado – some clarifying questions, none of them crucial
•	gmatch_summarize_results.do – only one line (I fixed the order of two columns that were reversed – Lauren is already using a modified version)
•	gmatchclass.mata
  -	mostly just typo-fixes in your comments
  -	some clarifying questions about the comments
  -	a few clarifying questions about the code
  -	only substantive thing is the non-crucial thing I asked you about when we met about the possible
    (slight) differences between your formulas and formulas in cbps.ado but maybe I’m missing something
    there (either way doesn’t actually matter since as we discussed the N’s are just scaling factor, right?)

Two more general questions:

A.	I might have asked you this when we met. Do you want any special error messages or flag when
convergence is not achieved?

B.	I didn’t have time to really think this through more carefully but while I was reviewing the mata
class I couldn’t help wondering if perhaps there could be less functions in the class and if that would
make it more readable. I know more functions allows them each to be a smaller building block which I
think programmers generally greatly value but as I was trying to parse things, I felt like I was having
to skip around the code and do a lot of crtl-F, chasing down each building block. So maybe some things
could be combined? (I wrote this note to myself while I was deep into the code a few weeks ago. It’s too
bad I didn’t have time to write down an example illustrating what I was thinking because to be honest I
no longer know if this comment is at all valid, it’s just a general feeling.)

Three project specific notes:

I.	You probably want to ask Lauren to clean up the project folders (lots of gmatch_run_2.do,
gmatch_run_3.do…different versions of variable_lists.do). I’m sure she knows what it all is and probably
she already has this on her to-do list but obviously you don’t want to leave things looking like this.

II.	Even after things are cleaned up, a readme would be helpful (and probably better if she does it
rather than you, so she can write down the stuff that a user actually needs to know about running the
programs and the overall process – where changes might be made in future rounds and stuff like that).

III.	I know you guys are ok with the “black box” –ness of the loss function and you are happy with
the set of weights that the process came up with but a few notes about the grids you chose to search
over and the values you chose to use (and not search over) would be helpful for someone looking at this
next year. I’m not saying there is anything wrong with what you did but you did not choose these numbers
with a random number generator (I don’t think) so just a few notes might be helpful both for this project
and possible future users thinking about what parameters they want to use.

*/


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
            ipw cbps mean_sd sd mean_sd_sq sd_sq STDProgdiff /// to fill in fctn and oid
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
  if (!mi("`maxtarget'")  & mi("`kurttarget'")) local kurttarget "0 0 2"
  if (!mi("`kurttarget'") & mi("`skewtarget'")) local skewtarget "0 0 2"
  if (!mi("`skewtarget'") & mi("`cvtarget'"))   local cvtarget   "0 0 2"
  local cvopt "`cvtarget' `skewtarget' `kurttarget' `maxtarget'"
  local cvopt : list clean cvopt
  if (!inlist(`: list sizeof cvopt',0,3,6,9,12)) {
    di as error `"cvopt() requires 3, 6, 9, or 12 elements"' // LP: given above won't cvopt always have 9 or 12 items?
    error 198
  }

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
  di as txt "Generalization of IPW/CPBS-type reweighting"
  if      ("`fctn'"=="ipw"        ) di as txt "Loss = IPW" _c
  else if ("`fctn'"=="cbps"       ) di as txt "Loss = CBPS" _c
  else if ("`fctn'"=="ipwcbps"    ) di as txt "Loss = CBPS + IPW (overidentified)" _c
  else if ("`fctn'"=="mean_sd_sq" ) di as txt "Loss = mean(stddiff())^2" _c
  else if ("`fctn'"=="sd_sq"      ) di as txt "Loss = sum(stddiff()^2)" _c
  else if ("`fctn'"=="stdprogdiff") di as txt "Loss = sum(stdprogdiff()^2)" _c
  tokenize `cvopt'
  if ("`1'"!="")  di as txt   " + `1'*abs(wgt_cv()-`2')^`3')" _c // LP: how is 1 ever missing given line parsing of cvopt above, seems like it would be zero
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

  temp = gmatch_ado_most_recent.gmatch(est, fctn, denominator, cvopt) // LP: where is temp used?
  gmatch_ado_most_recent.get_scores("_weight _weight_mtch _pscore _treated", tousevar)
  // LP: the variable names could change but the order of the concepts cannot, right?
}
end
