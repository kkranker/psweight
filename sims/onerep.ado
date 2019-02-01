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
  syntax namelist(name=scenariolist id="scenario(s)"), /// DGP process
         n(numlist min=1 >0 integer sort) /// loop over one or more sample sizes
         Impacts(numlist sort) /// loop over one or more impact estimates
         ESTimators(namelist) /// pick wich estimators to include (IPW, CBPS, etc.) see below
       [ NOIse(passthru) /// passed to DGP
         WNOIse(passthru) /// passed to DGP
         wcoef(passthru) /// passed to DGP
         HISTogram /// passed to DGP
         CVTargets(numlist >=10 <=100) /// passed to gmatch
         ate atet ateu /// atet is the default; passed to reweighting commands (teffects, gmatch, etc)
         AUGmented /// run OLS models to estimate impacts
         vce(passthru) /// e.g., add robust standard errors in outcome models
         *] // display options passed everywhere; remaining options passed to gmatch.ado (if applicable)
  _get_diopts diopts options, `options'

  // throw error if invalid options
  local valid_est raw ipw_true_ps ipw ipw_te stdprogdiff cbps
  if !`:list estimators in valid_est' {
    di as error "estimators(`: list estimators-valid_est') invalid"
    error 198
  }
  if ("`ate'`atet'`ateu'"=="") {
    local atet atet
  }
  if ("`augmented'"=="augmented") {
    local omvarlist "w1-w10"
    local aug "aug"
  }
  set matsize 2000

  tempname _b_ add from
  local c = 0
  foreach scenario of local scenariolist {
    foreach impact of local impacts {
      local L : word count `n'
      dgp_ssbgc `scenario', n(`: word `L' of `n'') impact(`impact') `noise' `wnoise' `wcoef' `histogram'

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

        // Difference in means ("raw")
        local e "raw"
        if `: list e in estimators' {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `augmented'" _n(2)
          regress y i.a `omvarlist', `vce' noheader `diopts'
          addstats `_b_' 1.a `prefix'_`e'`aug'
        }

        // Difference in means, weighted using true propensity scores ("true")
        local e "ipw_true_ps"
        if `: list e in estimators' cap nois {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `augmented'" _n(2)
          tempvar trueW
          if "`ate'`atet'`ateu'"=="ate" {
            gen `trueW' = cond(a, 1/ps, 1/(1-ps))
            summ `trueW' if a
            qui replace `trueW' = `trueW' / r(mean) if a
            summ `trueW' if !a
            qui replace `trueW' = `trueW' / r(mean) if !a
          }
          else if "`ate'`atet'`ateu'"=="atet" {
            gen `trueW' = cond(a, 1, ps/(1-ps))
            summ `trueW' if !a
            qui replace `trueW' = `trueW' / r(mean) if !a
          }
          else if "`ate'`atet'`ateu'"=="atet" {
            gen `trueW' = cond(a, (1-ps)/ps, 1)
            summ `trueW' if a
            qui replace `trueW' = `trueW' / r(mean) if a
          }
          else {
            di as error "`ate'`atet'`ateu' invalid"
            error 198
          }
          regress y i.a `omvarlist' [aw=`trueW'], `vce' noheader `diopts'
          addstats `_b_' 1.a `prefix'_`e'`aug'
        }

        // IPW model (with teffects)
        local e "ipw_te"
        if `: list e in estimators' cap nois {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `augmented'" _n(2)
          if ("`augmented'"=="augmented") teffects aipw (y `omvarlist') (a w1-w10), `ate'`atet'`ateu' aeq `diopts'
          else teffects ipw (y) (a w1-w10), `ate'`atet'`ateu' aeq `diopts'
          matrix `from' = e(b)
          matrix `from' = `from'[1, "TME1:"]
          local fromopt from(`from')
          addstats `_b_' ATET:r1vs0.a `prefix'_`e'`aug'
        }

        // IPW model (with gmatch.ado)
        local e "ipw"
        if `: list e in estimators' cap nois  {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `augmented'" _n(2)
          gmatch a w1-w10, ipw `ate'`atet'`ateu' `fromopt' `options' `diopts'
          matrix `from' = e(b)
          local fromopt from(`from')
          regress y i.a `omvarlist' [aw=_weight], `vce' noheader `diopts'
          addstats `_b_' 1.a `prefix'_`e'`aug'
        }

        // Minimize difference in prognostic scores model (with gmatch.ado)
        local e "stdprogdiff"
        if `: list e in estimators' cap nois {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `augmented'" _n(2)
          gmatch a w1-w10, stdprogdiff depvar(y) `ate'`atet'`ateu' `fromopt' `options' `diopts'
          matrix `from' = e(b)
          local fromopt from(`from')
          regress y i.a `omvarlist' [aw=_weight], `vce' noheader `diopts'
          addstats `_b_' 1.a `prefix'_`e'`aug'
          local cvcbps = r(wgt_cv)
          mac list _cvcbps
        }

        // CBPS model (with gmatch.ado)
        local e "cbps"
        if `: list e in estimators' {
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `augmented'" _n(2)
          gmatch a w1-w10, cbps `ate'`atet'`ateu' `fromopt' `options' `diopts'
          matrix `from' = e(b)
          local fromopt from(`from')
          regress y i.a `omvarlist' [aw=_weight], `vce' noheader `diopts'
          addstats `_b_' 1.a `prefix'_`e'`aug'
          local cvcbps = r(wgt_cv)
          mac list _cvcbps
        }

        // CBPS model (with gmatch.ado), with CV at X%
        foreach cut of local cvtargets {
          cap nois {
          local e = strtoname("cbps`cut'")
          local confirm cbps
          if (!`: list confirm in estimators') continue
          if (`cut'==100) continue
          local cvtarget = (`cvcbps'*`cut'/100)
          di _n(2) as txt "`prefix' with estimator: " as res "`e'" as txt " `augmented'" _n ///
                   as txt "CV target: cvtarget(20 " as res %7.4f `cvtarget' as txt " 6)" _n(2)
          gmatch a w1-w10, cbps cvtarget(20 `cvtarget' 6) `ate'`atet'`ateu' `fromopt' `options' `diopts'
          regress y i.a `omvarlist' [aw=_weight], `vce' noheader `diopts'
          addstats `_b_' 1.a `prefix'_`e'`aug'
          }
        }

        local --l
      }
    }
  }

  ereturn clear
  // mat list `_b_'
  ereturn post `_b_'
  ereturn local cmd "onerep"
  ereturn local cmdline "`cmd'"
  ereturn local cmdline "Simulation results"
  ereturn display, `diopts'
end


// helper program to add impact esimate (named "coef" to matrix
program define addstats
  version 15.1
  args matname coef eqname
  tempname add cell

  // impact estimate, standard error, bias, MSE
  mat `add' = (_b[`coef'], /// beta
               _se[`coef'], /// standard error
               (_b[`coef'] - _g1), /// impact estimate - true
               ((_b[`coef'] - _g1)^2)) // (impact estimate - true)^2
  mat colnames `add' = impact_est sd_error bias error_sqr

  // power to detect true effect given SE
  cap nois {
    mata: power_zstat(st_matrix("`add'")[2], st_numscalar("_g1"), .05, 2)
    mat `cell' = (r(beta))
    mat colnames `cell' = power_zstat_0
    mat `add' = (`add', `cell')
  }

  // Ho: null of effect = 0
  // lincom is just a conventient way to get p-value and CIs.
  cap nois {
    return clear
    lincom _b[`coef']
    mat `cell' = (r(p), /// grab p-value
                  (r(p) <= (1 - c(clevel) / 100)), /// reject null?
                  ((r(lb) <= _g1) & (_g1 <= r(ub)))) // coverage
    mat colnames `cell' = p_0 reject_0 covered
    mat `add' = (`add', `cell')
  }

  cap nois {
    gmatchcall balanceresults()
    mat `cell'  = (r(max_asd), /// Maximum absolute standardized diff.
                   r(mean_asd), /// Mean absolute standardized diff.
                   r(wgt_sd), /// S.D. of matching weights
                   r(wgt_cv), /// C.V. of matching weights
                   r(wgt_skewness), /// Skewness of matching weights
                   r(wgt_kurtosis), /// Kurtosis of matching weights
                   r(wgt_max)) // Maximum matching weight
    mat colnames `cell' = bal_max_asd bal_mean_asd wgt_sd wgt_cv wgt_skewness wgt_kurtosis wgt_max //
    mat `add' = (`add', `cell')
  }

  matrix coleq `add' = `eqname'
  matrix `matname' = (nullmat(`matname'), `add')
  // mat list `add'

end

mata:
  mata set matastrict on
  void power_zstat(real scalar se, real scalar effect, real scalar alpha , real scalar sides) {
    if (alpha<=0 | alpha>=1)  _error("The third argument (alpha) must be between 0 and 1")
    if (sides!=2 & sides!=1)  _error("The fourth argument (sides) must be 1 or 2")
    beta = normal(abs(effect/se) + invnormal(alpha/sides))
    st_numscalar("r(beta)", beta)
    strofreal(round(beta#100,.01)) + " percent power to detect an effect of " + strofreal(effect) + ", assuming a standard error of " + strofreal(se) + " and a " + strofreal(sides) + "-sided test and alpha = " + strofreal(alpha)
  }
end
