*! $Id$
*! One replication of the Simulations with DGP_KSIR.ADO
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

program define rf_xval_ksir, eclass
  version 15.1

  local cmd `"rf_xval_ksir `0'"'

  syntax anything(name=cvparams),  ///
         n(numlist min=1 >0 integer sort) /// loop over one or more sample sizes
       [ TRUEpscore /// use Z, rather than X in p-score model
         QUIETly /// supress a lot of the output
         rfopts(string) /// passed to randomforest
         *] // display options passed through
  _get_diopts diopts , `options'

  // process options
  // throw error if invalid options
  local estimators rf
  if ("`truepscore'"=="truepscore") {
    local pscorevarlist "z1-z4"
  }
  else {
    local pscorevarlist "x1-x4"
  }
  if ("`quietly'"=="") local quietly noisily
  set matsize 2000

  tempname _b_ add from
  local c = 0
  local scenario I
    local impact = 0
      local L : word count `n'

      dgp_ksir, n(`: word `L' of `n'') `histogram'

      scalar _g1 = `impact'

      local l = `L'
      while (`l' > 0) {
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

        // X-validate super-- parameters for the Randomforest model used to estimate P-scores

          qui `quietly' gridsearchcv `cvparams', scoring("neg_mean_squared_error"): ///
             randomforest a `pscorevarlist', type(class) `rfopts' `diopts'

          mata: resulttable = resulttable \  ("`thisN'", "`r(best_score)'", "`r(best_params)'")

        local --l
      }

  ereturn local cmd `"`cmd'"'

end
