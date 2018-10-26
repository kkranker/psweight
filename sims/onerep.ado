*! $Id$
*! One replication of the Simulations
// User specifies
//   - one or more simulation DGP (A to G)
//   - one or more sample sizes of dataset(s)
//   - one or more sizes of impacts
//   - models to run (on each dataset)
// The command then loops over DGPs, Ns, and Impacts.
// In each loop it
//   - generates the data (calls dgp_ssbgc.ado)
//   - estimates one or more models
//   - computes various statistics
// All the results are returned in one wide matrix (e(b))
//
*! By Keith Kranker
// Last updated $Date$
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.

program define onerep, eclass
  version 15.1

  local cmd `"onerep `0'"'
  syntax namelist(name=scenariolist id="scenario(s)"), ///
         n(numlist min=1 >0 integer sort) ///
         Impacts(numlist sort) ///
         ESTimators(namelist)

  tempname _b_ add
  foreach scenario of local scenariolist {
    local c = 0
    foreach impact of local impacts {
      local L : word count `n'
      dgp_ssbgc `scenario', n(`: word `L' of `n'') impact(`impact')

      local l = `L'
      while `l' > 0 {
        local thisN : word `l' of `n'
        local prefix = strtoname(trim("`scenario'_`++c'"))

        keep in 1/`thisN'
        di _n(3) as txt _dup(20) "-" _n ///
                 as txt "`prefix'" _n ///
                 as txt _dup(20) "-" _n ///
                 as txt "Scenario: " as res "`scenario'" _n ///
                 as txt "Impact:   " as res "`impact'" _n ///
                 as txt "N:        " as res "`thisN'" _n ///
                 as txt _dup(20) "-"  _n(2)

        qui count if a
        matrix          `add' = (`thisN', r(N), _g1)
        matrix colnames `add' = N Nt true
        matrix coleq    `add' = `prefix'
        matrix          `_b_' = (nullmat(`_b_'), `add')

        // IPW model
        local e "raw"
        if `: list e in estimators' {
          regress y i.a, noheader
          addstats `_b_' 1.a `prefix'_ols
        }

        // IPW model
        local e "ipw"
        if `: list e in estimators' {
          gmatch a w1-w10, atet ipw pooledvariance
          regress y i.a [aw=_weight], noheader
          addstats `_b_' 1.a `prefix'_ipw
        }

        local --l
      }
    }
  }

  ereturn clear
  mat list `_b_'
  ereturn post `_b_'
  ereturn local cmd "onerep"
  ereturn local cmdline "`cmd'"
  ereturn local cmdline "Simulation results"
end


// helper program to add impact esimate (named "coef" to matrix
program define addstats
  version 15.1
  args matname coef eqname
  tempname add cell

  mat `add' = (_b[`coef'], /// beta
               _se[`coef'], /// standard error
               (_b[`coef'] - _g1), /// impact estimate - true
               (_b[`coef'] - _g1)^2) // (impact estimate - true)^2
  mat colnames `add' = impact_est sd_error bias error_sqr

  cap nois {
    return clear
    test _b[`coef']=0
    mat `cell' = (r(p), (r(p) <= (1-c(clevel)/100)))
    mat colnames `cell' = p_0 reject_0
    mat `add' = (`add', `cell')
  }

  cap nois {
    return clear
    test _b[`coef']=_g1
    mat `cell' = (r(p), (r(p) <= (1-c(clevel)/100)))
    mat colnames `cell' = p_g1 reject_g1
    mat `add' = (`add', `cell')
  }

  cap nois {
    gmatchcall balanceresults()
    mat `cell'  = (r(max_asd), /// Maximum absolute standardized diff.
                   r(mean_asd), /// Mean absolute standardized diff.
                   r(wgt_sd), /// S.D. of matching weights:
                   r(wgt_cv), /// C.V. of matching weights:
                   r(wgt_skewness), /// Skewness of matching weights:
                   r(wgt_kurtosis), /// Kurtosis of matching weights:
                   r(wgt_max)) // Maximum matching weight:
    mat colnames `cell' = bal_max_asd bal_mean_asd wgt_sd wgt_cv wgt_skewness wgt_kurtosis wgt_max //
    mat `add' = (`add', `cell')
}


  matrix coleq `add' = `eqname'
  matrix `matname' = (nullmat(`matname'), `add')
  mat list `add'

end
