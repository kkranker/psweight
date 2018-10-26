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
  syntax [using/] [, clear replace]
  if (`"`using'"'!="") use "`using'", `clear'
  desc, short
  qui compress
  unab varlist : _all
  unab prefixes : *_b_N
  local prefixes: subinstr local prefixes "_b_N" "", all
  di as txt "There are `: list sizeof prefixes' prefixes"
  mac list _prefixes

  local pref1 : word 1 of `prefixes'
  unab ests: `pref1'_*_b_impact_est
  local ests: subinstr local ests "`pref1'_" "", all
  local ests: subinstr local ests "_b_impact_est" "", all
  di as txt "There are `: list sizeof ests' ests"
  mac list _ests

  // preserve
  local s = 0
  tempfile stack
  gen rep = _n
  preserve
  foreach prefix of local prefixes {
    foreach est of local ests {
      restore, preserve
      di "`prefix'/`est':"
      keep rep `prefix'_b_* `prefix'_`est'_b_*
      qui desc, varlist
      di as txt "   from: " as res =r(varlist)
      rename `prefix'_* *
      rename `est'_* *
      rename  b_* *
      gen dgp_txt = regexs(1) if regexm("`prefix'", "^([A-G]*)_([0-9]*)$")
      gen c       = regexs(2) if regexm("`prefix'", "^([A-G]*)_([0-9]*)$")
      gen est_txt = "`est'"
      qui desc, varlist
      di as txt "     to: " as res =r(varlist)
      if `++s'>1 append using "`stack'"
      qui save "`stack'", replace
    }
  }
  restore, not
  use "`stack'", clear

  foreach v of var _all {
    label var `v'
  }
  label define dgp 1 "A" 2 "B" 3 "C" 4 "D" 5 "E" 6 "F" 7 "G"
  label define est 1 "OLS" 2 "IPW" 3 "CBPS"
  encode dgp_txt, gen(dgp) label(dgp) noextend
  replace est_txt = upper(est_txt)
  encode est_txt, gen(est) label(est)
  qui destring c, replace

  drop dgp_txt est_txt
  order rep c dgp est, first
  isid rep c est
  sort rep c est

  if ("`replace'"!="" & `"`using'"'!="") save "`using'", replace
  desc
  codebook rep c dgp est, compact
  table true N est, by(dgp) concise

end
