*! $Id$
*! Reshapes the outpute from simmulate: from wide to tall
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.

program define sim_reshape
  version 15.1
  di _n as txt "Before reshape:" _n
  desc, short
  qui compress
  unab varlist : _all
  unab prefixes : *_b_N
  local prefixes: subinstr local prefixes "_b_N" "", all
  di as txt "`: list sizeof prefixes' prefixes:"
  mac list _prefixes
  foreach pref1 of local prefixes {
    unab  ests1: `pref1'_*_b_impact_est
    local ests1: subinstr local ests1 "`pref1'_" "", all
    local ests1: subinstr local ests1 "_b_impact_est" "", all
    local ests : list ests | ests1
  }
  local ests : subinstr local ests "aug" "", all
  local ests : list uniq ests
  di as txt "`: list sizeof ests' estimators:"
  mac list _ests
  local ee = 0
  local estdefine `++ee' "RAW" `++ee' "IPW_TRUE_PS" `++ee' "STDPROGDIFF" `++ee' "IPW_TE" `++ee' "IPW" `ee' "ELASTIC" `++ee' "RF" `++ee' "IPWCBPS" `++ee' "CBPS"
  foreach est of local ests {
    if !regexm(strupper("`est'"), "CBPS[0-9].*") continue
    local estdefine `estdefine' `++ee' `=strupper("`est'")'
  }
  format %7.4f _all
  format %7.1f *_wgt_max*
  format %7.0g *_reject* *covered* *_N* *_Nt*

  local s = 0
  tempfile stack
  gen rep = _n
  preserve
  foreach prefix of local prefixes {
    foreach est of local ests {
      foreach aug in "" "AUG" {
        restore, preserve
        // di as res "`prefix'/`est':" _c
        cap nois keep rep `prefix'_b_* `prefix'_`est'_b_*
        if !_rc {
          rename `prefix'_* *
          rename `est'_* *
          rename  b_* *
          gen dgp_txt = regexs(1) if regexm("`prefix'", "^([A-Z]*)_([0-9]*)$")
          gen dataset = regexs(2) if regexm("`prefix'", "^([A-Z]*)_([0-9]*)$")
          gen augmented = ("`aug'"=="AUG")
          gen est_txt = "`est'"
          if (`++s'>1) qui append using "`stack'"
          qui save "`stack'", replace
        }
      }
    }
  }
  restore, not
  use "`stack'", clear

  // strings to numeric
  label define dgp 1 "A" 2 "B" 3 "C" 4 "D" 5 "E" 6 "F" 7 "G" 8 "I" 9 "K"
  encode dgp_txt, gen(dgp) label(dgp)
  label define est `estdefine'
  qui replace est_txt = upper(est_txt)
  encode est_txt, gen(estimator) label(est)
  qui destring dataset, replace
  drop dgp_txt est_txt
  label define aug 0 "Unadjusted" 1 "Reg. adjusted"
  label val augmented aug

  // "result" is a unique ID combiining "dataset" and "estimator"
  // dataset is one loop of "one rep" -- one per scenario per impact per
  // sanmple size.  However, each daaset can have multiple estimators
  cap nois isid rep dgp true N estimator augmented
  cap nois isid rep dataset    estimator augmented
  egen result = group(dataset dgp true N estimator augmented), label missing
  isid rep result

  // claculate variance of impact_est
  // store on the last row
  tempvar mean_bias MSE count
  bys result (rep): egen double `mean_bias' = mean(bias)
  by  result (rep): egen double `MSE' = mean(error_sqr)
  by  result (rep): egen double `count' = count(impact_est)
  by  result (rep): gen  double impact_est_var = (`MSE' - `mean_bias'^2) * `count' / (`count' - 1) if _n == _N
  format `:format impact_est' impact_est_var
  drop `mean_bias' `MSE' `count'
  gen double rmse = sqrt(error_sqr)
  order rmse, after(error_sqr)

  // checks, cleanup
  order result dataset dgp true N estimator augmented rep, first
  sort  result dataset dgp true N estimator augmented rep
  foreach v of var _all {
    label var `v'
  }
  qui compress

  // save/summarize
  di _n as txt "After reshape:" _n
  desc
  summ, sep(0)
  table estimator N true, by(augmented dgp) concise
end
